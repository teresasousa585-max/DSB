`timescale 1ns / 1ps

module pwm32_generator (
    input  logic               clk,
    input  logic               rst_n,
    input  logic signed [15:0] audio_in,
    input  logic [7:0]         amplitude [0:31],
    input  logic [11:0]        phase_del [0:31],
    output logic [31:0]        pwm_out
);

    // 100 MHz / 2500 = 40 kHz carrier.
    logic [11:0] carrier_cnt;
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n)
            carrier_cnt <= 12'd0;
        else if (carrier_cnt == 12'd2499)
            carrier_cnt <= 12'd0;
        else
            carrier_cnt <= carrier_cnt + 1'b1;
    end

    // Signed 16-bit audio is offset into the sqrt-mapping ROM address range.
    logic [15:0] rom_addr;
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n)
            rom_addr <= 16'd0;
        else
            rom_addr <= audio_in + 16'd32768;
    end

    logic [11:0] duty_base_d2;

    sqrt_mapping_rom u_math_rom (
        .clka  (clk),
        .addra (rom_addr),
        .ena   (1'b1),
        .douta (duty_base_d2)
    );

    genvar i;
    generate
        for (i = 0; i < 32; i++) begin : gen_ch
            logic [19:0] duty_ch_mult;
            assign duty_ch_mult = duty_base_d2 * {12'h0, amplitude[i]};

            logic [11:0] duty_ch_d3;
            logic [11:0] center_d3;
            always_ff @(posedge clk or negedge rst_n) begin
                if (!rst_n) begin
                    duty_ch_d3   <= 12'd0;
                    center_d3    <= 12'd0;
                end else begin
                    duty_ch_d3   <= duty_ch_mult[19:8];
                    center_d3    <= (phase_del[i] >= 12'd2500) ? (phase_del[i] - 12'd2500) : phase_del[i];
                end
            end

            logic [11:0] half_duty_d3;
            logic [12:0] start_centered_raw;
            logic [11:0] start_centered;
            assign half_duty_d3 = duty_ch_d3 >> 1;
            assign start_centered_raw = ({1'b0, center_d3} >= {1'b0, half_duty_d3}) ?
                                        ({1'b0, center_d3} - {1'b0, half_duty_d3}) :
                                        ({1'b0, center_d3} + 13'd2500 - {1'b0, half_duty_d3});
            assign start_centered = start_centered_raw[11:0];

            logic [11:0] duty_ch_d4;
            logic [11:0] start_d4;
            logic [12:0] end_raw_d4;
            always_ff @(posedge clk or negedge rst_n) begin
                if (!rst_n) begin
                    duty_ch_d4 <= 12'd0;
                    start_d4   <= 12'd0;
                    end_raw_d4 <= 13'd0;
                end else begin
                    duty_ch_d4 <= duty_ch_d3;
                    start_d4   <= start_centered;
                    end_raw_d4 <= {1'b0, start_centered} + duty_ch_d3;
                end
            end

            // phase_del is the pulse center; commit the next channel window only on a carrier boundary.
            logic [11:0] duty_locked;
            logic [11:0] start_locked;
            logic [11:0] end_locked;
            logic        wrap_locked;
            always_ff @(posedge clk or negedge rst_n) begin
                if (!rst_n) begin
                    duty_locked  <= 12'd0;
                    start_locked <= 12'd0;
                    end_locked   <= 12'd0;
                    wrap_locked  <= 1'b0;
                end else if (carrier_cnt == 12'd2499) begin
                    duty_locked  <= duty_ch_d4;
                    start_locked <= start_d4;
                    end_locked   <= (end_raw_d4 >= 13'd2500) ? (end_raw_d4 - 13'd2500) : end_raw_d4[11:0];
                    wrap_locked  <= (end_raw_d4 >= 13'd2500);
                end
            end

            assign pwm_out[i] = (duty_locked == 12'd0) ? 1'b0 :
                                (!wrap_locked) ? (carrier_cnt >= start_locked && carrier_cnt < end_locked) :
                                                 (carrier_cnt >= start_locked || carrier_cnt < end_locked);
        end
    endgenerate

endmodule
