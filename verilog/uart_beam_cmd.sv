//=============================================================================
// UART beam command receiver
//
// Frame format, UART 8N1:
//   0xAA 0x55 cmd seq len_l len_h payload... crc_l crc_h
//
// CRC:
//   CRC-16/CCITT-FALSE, init 0xFFFF, polynomial 0x1021.
//   Covered bytes: cmd, seq, len_l, len_h, payload.
//
// Commands:
//   0x10 WRITE_PARAMS: payload = 25 * {phase_l, phase_h[1:0], amplitude}
//   0x11 START:        len = 0, enable ultrasound output
//   0x12 STOP:         len = 0, disable ultrasound output immediately
//   0x13 SOFT_RESET:   len = 0, pulse ultrasound reset and reload active params
//=============================================================================

`timescale 1ns / 1ps

module uart_rx_8n1 #(
    parameter integer CLK_FREQ = 50_000_000,
    parameter integer BAUD     = 921_600
) (
    input  wire       clk,
    input  wire       rst_n,
    input  wire       rx,
    output reg  [7:0] data,
    output reg        data_valid,
    output reg        frame_error
);

    localparam integer CLKS_PER_BIT = CLK_FREQ / BAUD;
    localparam integer HALF_BIT     = CLKS_PER_BIT / 2;

    localparam [2:0] RX_IDLE  = 3'd0;
    localparam [2:0] RX_START = 3'd1;
    localparam [2:0] RX_DATA  = 3'd2;
    localparam [2:0] RX_STOP  = 3'd3;

    reg [2:0] state;
    reg [15:0] clk_cnt;
    reg [2:0] bit_idx;
    reg [7:0] shift_reg;
    reg rx_meta;
    reg rx_sync;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            rx_meta <= 1'b1;
            rx_sync <= 1'b1;
        end else begin
            rx_meta <= rx;
            rx_sync <= rx_meta;
        end
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state       <= RX_IDLE;
            clk_cnt     <= 16'd0;
            bit_idx     <= 3'd0;
            shift_reg   <= 8'd0;
            data        <= 8'd0;
            data_valid  <= 1'b0;
            frame_error <= 1'b0;
        end else begin
            data_valid  <= 1'b0;
            frame_error <= 1'b0;

            case (state)
                RX_IDLE: begin
                    clk_cnt <= 16'd0;
                    bit_idx <= 3'd0;
                    if (!rx_sync)
                        state <= RX_START;
                end

                RX_START: begin
                    if (clk_cnt == HALF_BIT) begin
                        if (!rx_sync) begin
                            clk_cnt <= 16'd0;
                            state   <= RX_DATA;
                        end else begin
                            state <= RX_IDLE;
                        end
                    end else begin
                        clk_cnt <= clk_cnt + 1'b1;
                    end
                end

                RX_DATA: begin
                    if (clk_cnt == CLKS_PER_BIT - 1) begin
                        clk_cnt <= 16'd0;
                        shift_reg[bit_idx] <= rx_sync;
                        if (bit_idx == 3'd7) begin
                            bit_idx <= 3'd0;
                            state   <= RX_STOP;
                        end else begin
                            bit_idx <= bit_idx + 1'b1;
                        end
                    end else begin
                        clk_cnt <= clk_cnt + 1'b1;
                    end
                end

                RX_STOP: begin
                    if (clk_cnt == CLKS_PER_BIT - 1) begin
                        clk_cnt <= 16'd0;
                        state   <= RX_IDLE;
                        if (rx_sync) begin
                            data       <= shift_reg;
                            data_valid <= 1'b1;
                        end else begin
                            frame_error <= 1'b1;
                        end
                    end else begin
                        clk_cnt <= clk_cnt + 1'b1;
                    end
                end

                default: state <= RX_IDLE;
            endcase
        end
    end

endmodule

module uart_beam_cmd #(
    parameter integer CLK_FREQ = 50_000_000,
    parameter integer BAUD     = 921_600,
    parameter integer CH_NUM   = 25
) (
    input  wire       clk,
    input  wire       rst_n,
    input  wire       uart_rx,
    input  wire       carrier_tick,

    output reg [7:0]  amplitude [0:CH_NUM-1],
    output reg [9:0]  phase     [0:CH_NUM-1],
    output reg        param_valid,

    output reg        ultrasound_en,
    output reg        ultrasound_soft_rst,

    output reg        frame_ok,
    output reg        frame_error,
    output reg [7:0]  last_cmd,
    output reg [7:0]  last_seq,
    output wire [7:0] status
);

    localparam [7:0] SOF0 = 8'hAA;
    localparam [7:0] SOF1 = 8'h55;

    localparam [7:0] CMD_WRITE_PARAMS = 8'h10;
    localparam [7:0] CMD_START        = 8'h11;
    localparam [7:0] CMD_STOP         = 8'h12;
    localparam [7:0] CMD_SOFT_RESET   = 8'h13;

    localparam integer PARAM_BYTES = CH_NUM * 3;
    localparam [15:0] PARAM_LEN    = CH_NUM * 3;
    localparam integer MAX_PAYLOAD = PARAM_BYTES;

    localparam [3:0] ST_SOF0    = 4'd0;
    localparam [3:0] ST_SOF1    = 4'd1;
    localparam [3:0] ST_CMD     = 4'd2;
    localparam [3:0] ST_SEQ     = 4'd3;
    localparam [3:0] ST_LEN_L   = 4'd4;
    localparam [3:0] ST_LEN_H   = 4'd5;
    localparam [3:0] ST_PAYLOAD = 4'd6;
    localparam [3:0] ST_CRC_L   = 4'd7;
    localparam [3:0] ST_CRC_H   = 4'd8;

    wire [7:0] rx_byte;
    wire       rx_valid;
    wire       rx_frame_error;

    uart_rx_8n1 #(
        .CLK_FREQ(CLK_FREQ),
        .BAUD    (BAUD)
    ) u_uart_rx (
        .clk         (clk),
        .rst_n       (rst_n),
        .rx          (uart_rx),
        .data        (rx_byte),
        .data_valid  (rx_valid),
        .frame_error (rx_frame_error)
    );

    reg [3:0]  state;
    reg [7:0]  cmd_reg;
    reg [7:0]  seq_reg;
    reg [15:0] len_reg;
    reg [15:0] payload_idx;
    reg [15:0] crc_reg;
    reg [7:0]  crc_l_reg;
    reg [7:0]  payload_buf [0:MAX_PAYLOAD-1];

    reg [7:0]  shadow_amp   [0:CH_NUM-1];
    reg [9:0]  shadow_phase [0:CH_NUM-1];
    reg        shadow_valid;
    reg        reload_pending;
    reg        active_loaded;

    integer i;

    reg [3:0] error_code;

    assign status = {
        error_code,
        shadow_valid,
        active_loaded,
        reload_pending,
        ultrasound_en
    };

    function automatic [15:0] crc16_next_byte;
        input [15:0] crc_in;
        input [7:0]  data_in;
        reg [15:0] crc;
        integer bit_i;
        begin
            crc = crc_in ^ {data_in, 8'h00};
            for (bit_i = 0; bit_i < 8; bit_i = bit_i + 1) begin
                if (crc[15])
                    crc = (crc << 1) ^ 16'h1021;
                else
                    crc = (crc << 1);
            end
            crc16_next_byte = crc;
        end
    endfunction

    function automatic is_supported_cmd;
        input [7:0] cmd;
        begin
            is_supported_cmd = (cmd == CMD_WRITE_PARAMS) ||
                               (cmd == CMD_START) ||
                               (cmd == CMD_STOP) ||
                               (cmd == CMD_SOFT_RESET);
        end
    endfunction

    function automatic is_valid_len;
        input [7:0]  cmd;
        input [15:0] len;
        begin
            case (cmd)
                CMD_WRITE_PARAMS: is_valid_len = (len == PARAM_LEN);
                CMD_START,
                CMD_STOP,
                CMD_SOFT_RESET:   is_valid_len = (len == 16'd0);
                default:          is_valid_len = 1'b0;
            endcase
        end
    endfunction

    task automatic accept_frame;
        integer ch;
        reg [9:0] next_phase;
        begin
            frame_ok <= 1'b1;
            last_cmd <= cmd_reg;
            last_seq <= seq_reg;
            error_code <= 4'd0;

            case (cmd_reg)
                CMD_WRITE_PARAMS: begin
                    for (ch = 0; ch < CH_NUM; ch = ch + 1) begin
                        next_phase = {payload_buf[ch*3 + 1][1:0], payload_buf[ch*3 + 0]};
                        shadow_phase[ch] <= next_phase;
                        shadow_amp[ch]   <= payload_buf[ch*3 + 2];
                    end
                    shadow_valid <= 1'b1;
                end

                CMD_START: begin
                    ultrasound_en <= 1'b1;
                    reload_pending <= active_loaded;
                end

                CMD_STOP: begin
                    ultrasound_en <= 1'b0;
                end

                CMD_SOFT_RESET: begin
                    ultrasound_en         <= 1'b0;
                    ultrasound_soft_rst   <= 1'b1;
                    reload_pending        <= active_loaded;
                end

                default: begin
                    frame_error <= 1'b1;
                    error_code  <= 4'd1;
                end
            endcase
        end
    endtask

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state                <= ST_SOF0;
            cmd_reg              <= 8'd0;
            seq_reg              <= 8'd0;
            len_reg              <= 16'd0;
            payload_idx          <= 16'd0;
            crc_reg              <= 16'hFFFF;
            crc_l_reg            <= 8'd0;
            param_valid          <= 1'b0;
            ultrasound_en        <= 1'b0;
            ultrasound_soft_rst  <= 1'b0;
            frame_ok             <= 1'b0;
            frame_error          <= 1'b0;
            last_cmd             <= 8'd0;
            last_seq             <= 8'd0;
            shadow_valid         <= 1'b0;
            reload_pending       <= 1'b0;
            active_loaded        <= 1'b0;
            error_code           <= 4'd0;

            for (i = 0; i < CH_NUM; i = i + 1) begin
                amplitude[i]     <= 8'd255;
                phase[i]         <= 10'd0;
                shadow_amp[i]    <= 8'd255;
                shadow_phase[i]  <= 10'd0;
            end
            for (i = 0; i < MAX_PAYLOAD; i = i + 1) begin
                payload_buf[i] <= 8'd0;
            end
        end else begin
            param_valid         <= 1'b0;
            ultrasound_soft_rst <= 1'b0;
            frame_ok            <= 1'b0;
            frame_error         <= rx_frame_error;

            if (rx_frame_error) begin
                state      <= ST_SOF0;
                error_code <= 4'd4;
            end

            if (carrier_tick) begin
                if (shadow_valid) begin
                    for (i = 0; i < CH_NUM; i = i + 1) begin
                        amplitude[i] <= shadow_amp[i];
                        phase[i]     <= shadow_phase[i];
                    end
                    shadow_valid   <= 1'b0;
                    active_loaded  <= 1'b1;
                    reload_pending <= 1'b0;
                    param_valid    <= 1'b1;
                end else if (reload_pending && active_loaded) begin
                    reload_pending <= 1'b0;
                    param_valid    <= 1'b1;
                end
            end

            if (rx_valid && !rx_frame_error) begin
                case (state)
                    ST_SOF0: begin
                        if (rx_byte == SOF0)
                            state <= ST_SOF1;
                    end

                    ST_SOF1: begin
                        if (rx_byte == SOF1)
                            state <= ST_CMD;
                        else if (rx_byte == SOF0)
                            state <= ST_SOF1;
                        else
                            state <= ST_SOF0;
                    end

                    ST_CMD: begin
                        cmd_reg <= rx_byte;
                        crc_reg <= crc16_next_byte(16'hFFFF, rx_byte);
                        if (is_supported_cmd(rx_byte))
                            state <= ST_SEQ;
                        else begin
                            frame_error <= 1'b1;
                            error_code  <= 4'd1;
                            state <= ST_SOF0;
                        end
                    end

                    ST_SEQ: begin
                        seq_reg <= rx_byte;
                        crc_reg <= crc16_next_byte(crc_reg, rx_byte);
                        state <= ST_LEN_L;
                    end

                    ST_LEN_L: begin
                        len_reg[7:0] <= rx_byte;
                        crc_reg <= crc16_next_byte(crc_reg, rx_byte);
                        state <= ST_LEN_H;
                    end

                    ST_LEN_H: begin
                        len_reg[15:8] <= rx_byte;
                        crc_reg <= crc16_next_byte(crc_reg, rx_byte);
                        payload_idx <= 16'd0;
                        if (!is_valid_len(cmd_reg, {rx_byte, len_reg[7:0]})) begin
                            frame_error <= 1'b1;
                            error_code  <= 4'd2;
                            state <= ST_SOF0;
                        end else if ({rx_byte, len_reg[7:0]} == 16'd0) begin
                            state <= ST_CRC_L;
                        end else begin
                            state <= ST_PAYLOAD;
                        end
                    end

                    ST_PAYLOAD: begin
                        payload_buf[payload_idx] <= rx_byte;
                        crc_reg <= crc16_next_byte(crc_reg, rx_byte);
                        if (payload_idx == len_reg - 1)
                            state <= ST_CRC_L;
                        payload_idx <= payload_idx + 1'b1;
                    end

                    ST_CRC_L: begin
                        crc_l_reg <= rx_byte;
                        state <= ST_CRC_H;
                    end

                    ST_CRC_H: begin
                        if ({rx_byte, crc_l_reg} == crc_reg) begin
                            accept_frame();
                        end else begin
                            frame_error <= 1'b1;
                            error_code  <= 4'd3;
                        end
                        state <= ST_SOF0;
                    end

                    default: state <= ST_SOF0;
                endcase
            end
        end
    end

endmodule
