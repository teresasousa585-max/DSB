`timescale 1ns / 1ps

module uart_rx #(
    parameter CLK_FRE    = 100,     // System clock in MHz.
    parameter BAUD_RATE  = 115200,
    parameter DATA_WIDTH = 8
) (
    input clk,
    input rst_n,
    input i_uart_rx,

    output reg [DATA_WIDTH-1:0] o_uart_data,
    output reg                  o_rx_done
);

    localparam integer CLK_COUNT  = ((CLK_FRE * 1000000) + (BAUD_RATE / 2)) / BAUD_RATE;
    localparam integer HALF_COUNT = CLK_COUNT / 2;

    localparam [1:0] S_IDLE  = 2'd0;
    localparam [1:0] S_START = 2'd1;
    localparam [1:0] S_DATA  = 2'd2;
    localparam [1:0] S_STOP  = 2'd3;

    (* ASYNC_REG = "TRUE" *) reg rx_meta;
    (* ASYNC_REG = "TRUE" *) reg rx_sync;
    reg rx_sync_d;

    reg [1:0] state;
    reg [31:0] baud_cnt;
    reg [3:0] bit_cnt;
    reg [DATA_WIDTH-1:0] rx_shift;

    wire rx_falling = rx_sync_d & ~rx_sync;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            rx_meta   <= 1'b1;
            rx_sync   <= 1'b1;
            rx_sync_d <= 1'b1;
        end else begin
            rx_meta   <= i_uart_rx;
            rx_sync   <= rx_meta;
            rx_sync_d <= rx_sync;
        end
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state       <= S_IDLE;
            baud_cnt    <= 32'd0;
            bit_cnt     <= 4'd0;
            rx_shift    <= {DATA_WIDTH{1'b0}};
            o_uart_data <= {DATA_WIDTH{1'b0}};
            o_rx_done   <= 1'b0;
        end else begin
            o_rx_done <= 1'b0;

            case (state)
                S_IDLE: begin
                    baud_cnt <= 32'd0;
                    bit_cnt  <= 4'd0;
                    if (rx_falling) begin
                        state <= S_START;
                    end
                end

                S_START: begin
                    if (baud_cnt == HALF_COUNT - 1) begin
                        baud_cnt <= 32'd0;
                        if (!rx_sync) begin
                            state <= S_DATA;
                        end else begin
                            state <= S_IDLE;
                        end
                    end else begin
                        baud_cnt <= baud_cnt + 1'b1;
                    end
                end

                S_DATA: begin
                    if (baud_cnt == CLK_COUNT - 1) begin
                        baud_cnt <= 32'd0;
                        rx_shift[bit_cnt] <= rx_sync;

                        if (bit_cnt == DATA_WIDTH - 1) begin
                            bit_cnt <= 4'd0;
                            state   <= S_STOP;
                        end else begin
                            bit_cnt <= bit_cnt + 1'b1;
                        end
                    end else begin
                        baud_cnt <= baud_cnt + 1'b1;
                    end
                end

                S_STOP: begin
                    if (baud_cnt == CLK_COUNT - 1) begin
                        baud_cnt <= 32'd0;
                        state    <= S_IDLE;

                        if (rx_sync) begin
                            o_uart_data <= rx_shift;
                            o_rx_done   <= 1'b1;
                        end
                    end else begin
                        baud_cnt <= baud_cnt + 1'b1;
                    end
                end

                default: begin
                    state <= S_IDLE;
                end
            endcase
        end
    end

endmodule
