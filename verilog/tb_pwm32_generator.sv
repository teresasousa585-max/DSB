`timescale 1ns / 1ps

module tb_pwm32_generator;

    logic clk;
    logic rst_n;
    logic signed [15:0] audio_in;
    logic [7:0] amplitude [0:31];
    logic [11:0] phase_del [0:31];
    logic [31:0] pwm_out;

    integer pass_cnt;
    integer error_cnt;
    integer i;
    integer seen_high;
    integer seen_low;

    pwm32_generator dut (
        .clk       (clk),
        .rst_n     (rst_n),
        .audio_in  (audio_in),
        .amplitude (amplitude),
        .phase_del (phase_del),
        .pwm_out   (pwm_out)
    );

    initial begin
        clk = 1'b0;
        forever #5 clk = ~clk;
    end

    function automatic has_unknown;
        input [31:0] value;
        begin
            has_unknown = (^value === 1'bx);
        end
    endfunction

    task automatic check;
        input condition;
        input [1023:0] message;
        begin
            if (condition) begin
                pass_cnt = pass_cnt + 1;
                $display("[PASS] %0s", message);
            end else begin
                error_cnt = error_cnt + 1;
                $display("[FAIL] %0s", message);
            end
        end
    endtask

    initial begin
        $dumpfile("tb_pwm32_generator.vcd");
        $dumpvars(0, tb_pwm32_generator);

        pass_cnt = 0;
        error_cnt = 0;
        audio_in = 16'sd0;
        for (i = 0; i < 32; i = i + 1) begin
            amplitude[i] = 8'd0;
            phase_del[i] = 12'd0;
        end

        rst_n = 1'b0;
        repeat (8) @(posedge clk);
        check(pwm_out == 32'd0, "PWM outputs are low during reset");
        rst_n = 1'b1;

        repeat (3000) @(posedge clk);
        check(!has_unknown(pwm_out), "PWM outputs are not X after reset release");
        check(pwm_out == 32'd0, "zero amplitudes keep PWM outputs low");

        amplitude[0] = 8'd255;
        phase_del[0] = 12'd100;
        repeat (3000) @(posedge clk);

        seen_high = 0;
        seen_low = 0;
        repeat (3000) begin
            @(posedge clk);
            if (pwm_out[0])
                seen_high = 1;
            else
                seen_low = 1;
        end
        check(seen_high && seen_low, "channel 0 toggles with nonzero amplitude");
        check(pwm_out[31:1] == 31'd0, "other zero-amplitude channels stay low");

        amplitude[0] = 8'd0;
        repeat (3000) @(posedge clk);
        check(pwm_out[0] == 1'b0, "channel 0 returns low after amplitude clears");

        $display("Simulation summary: passes=%0d errors=%0d", pass_cnt, error_cnt);
        if (error_cnt == 0)
            $display("RESULT: ALL TESTS PASSED");
        else
            $display("RESULT: %0d CHECK(S) FAILED", error_cnt);
        $finish;
    end

    initial begin
        #2000000;
        $display("[FAIL] simulation timeout");
        $finish;
    end

endmodule

module sqrt_mapping_rom (
    input  wire        clka,
    input  wire [15:0] addra,
    input  wire        ena,
    output logic [11:0] douta
);
    always_ff @(posedge clka) begin
        if (ena)
            douta <= 12'd800;
    end

    // Keep address visible to lint/simulation without changing behavior.
    wire _unused_addr = ^addra;
endmodule
