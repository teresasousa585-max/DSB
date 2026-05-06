`timescale 1ns / 1ps

module tb_top;

    localparam integer CLK_PERIOD_NS = 20;
    localparam integer UART_BIT_NS   = 1080; // 50 MHz / 921600 -> 54 clocks

    reg clk;
    reg rst_n;
    reg uart_rx;
    wire [7:0] cmd_status;
    wire [7:0] last_cmd;
    wire [7:0] last_seq;
    wire frame_ok;
    wire frame_error;
    wire ultrasound_en;
    wire adc_cs;
    wire adc_clk;
    wire adc_data;
    wire [24:0] pwm_out;

    integer pass_cnt;
    integer error_cnt;
    integer ok_count;
    integer err_count;
    integer param_count;
    reg [7:0] seq;
    integer i;
    integer seen_pwm;
    integer start_ok;
    integer start_err;
    integer start_param;
    reg [15:0] tx_crc;
    reg [9:0] expected_phase [0:24];
    reg [7:0] expected_amp   [0:24];

    assign adc_data = 1'b1;

    top_ultrasound_array dut (
        .clk            (clk),
        .rst_n          (rst_n),
        .uart_rx        (uart_rx),
        .cmd_status     (cmd_status),
        .last_cmd       (last_cmd),
        .last_seq       (last_seq),
        .frame_ok       (frame_ok),
        .frame_error    (frame_error),
        .ultrasound_en  (ultrasound_en),
        .adc_cs         (adc_cs),
        .adc_clk        (adc_clk),
        .adc_data       (adc_data),
        .pwm_out        (pwm_out)
    );

    initial begin
        clk = 1'b0;
        forever #(CLK_PERIOD_NS / 2) clk = ~clk;
    end

    initial begin
        rst_n   = 1'b0;
        uart_rx = 1'b1;
        repeat (20) @(posedge clk);
        rst_n = 1'b1;
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            ok_count    <= 0;
            err_count   <= 0;
            param_count <= 0;
        end else begin
            if (frame_ok)
                ok_count <= ok_count + 1;
            if (frame_error)
                err_count <= err_count + 1;
            if (dut.u_uart_beam_cmd.param_valid)
                param_count <= param_count + 1;
        end
    end

    function automatic [15:0] crc16_next_byte;
        input [15:0] crc_in;
        input [7:0] data_in;
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

    task automatic wait_ok_after;
        input integer old_count;
        input [1023:0] message;
        integer timeout;
        begin
            timeout = 0;
            while (ok_count <= old_count && timeout < 200000) begin
                @(posedge clk);
                timeout = timeout + 1;
            end
            check(ok_count > old_count, message);
        end
    endtask

    task automatic wait_err_after;
        input integer old_count;
        input [1023:0] message;
        integer timeout;
        begin
            timeout = 0;
            while (err_count <= old_count && timeout < 200000) begin
                @(posedge clk);
                timeout = timeout + 1;
            end
            check(err_count > old_count, message);
        end
    endtask

    task automatic wait_param_after;
        input integer old_count;
        input [1023:0] message;
        integer timeout;
        begin
            timeout = 0;
            while (param_count <= old_count && timeout < 200000) begin
                @(posedge clk);
                timeout = timeout + 1;
            end
            check(param_count > old_count, message);
        end
    endtask

    task automatic uart_send_byte;
        input [7:0] value;
        integer bit_i;
        begin
            uart_rx = 1'b0;
            #(UART_BIT_NS);
            for (bit_i = 0; bit_i < 8; bit_i = bit_i + 1) begin
                uart_rx = value[bit_i];
                #(UART_BIT_NS);
            end
            uart_rx = 1'b1;
            #(UART_BIT_NS);
        end
    endtask

    task automatic send_crc_byte;
        input [7:0] value;
        begin
            tx_crc = crc16_next_byte(tx_crc, value);
            uart_send_byte(value);
        end
    endtask

    task automatic send_cmd_no_payload;
        input [7:0] cmd;
        input [7:0] seq_in;
        begin
            tx_crc = 16'hFFFF;
            uart_send_byte(8'hAA);
            uart_send_byte(8'h55);
            send_crc_byte(cmd);
            send_crc_byte(seq_in);
            send_crc_byte(8'd0);
            send_crc_byte(8'd0);
            uart_send_byte(tx_crc[7:0]);
            uart_send_byte(tx_crc[15:8]);
        end
    endtask

    task automatic send_bad_crc_cmd;
        input [7:0] cmd;
        input [7:0] seq_in;
        begin
            tx_crc = 16'hFFFF;
            uart_send_byte(8'hAA);
            uart_send_byte(8'h55);
            send_crc_byte(cmd);
            send_crc_byte(seq_in);
            send_crc_byte(8'd0);
            send_crc_byte(8'd0);
            uart_send_byte(~tx_crc[7:0]);
            uart_send_byte(tx_crc[15:8]);
        end
    endtask

    task automatic send_write_params;
        input [7:0] seq_in;
        input [9:0] phase_seed;
        input [7:0] amp_seed;
        integer ch;
        reg [9:0] p;
        reg [7:0] a;
        begin
            tx_crc = 16'hFFFF;
            uart_send_byte(8'hAA);
            uart_send_byte(8'h55);
            send_crc_byte(8'h10);
            send_crc_byte(seq_in);
            send_crc_byte(8'd75);
            send_crc_byte(8'd0);
            for (ch = 0; ch < 25; ch = ch + 1) begin
                p = phase_seed + (ch * 10'd37);
                a = amp_seed + ch;
                expected_phase[ch] = p;
                expected_amp[ch]   = a;
                send_crc_byte(p[7:0]);
                send_crc_byte({6'd0, p[9:8]});
                send_crc_byte(a);
            end
            uart_send_byte(tx_crc[7:0]);
            uart_send_byte(tx_crc[15:8]);
        end
    endtask

    initial begin
        $dumpfile("tb_top.vcd");
        $dumpvars(0, tb_top);

        pass_cnt = 0;
        error_cnt = 0;
        seq = 1;

        $display("============================================================");
        $display("UART controlled ultrasonic array testbench");
        $display("============================================================");

        @(posedge rst_n);
        repeat (10) @(posedge clk);

        check(ultrasound_en == 1'b0, "ultrasound output disabled after reset");
        check(pwm_out == 25'd0, "PWM outputs are low after reset");
        check(cmd_status[7:4] == 4'd0, "command status has no latched error after reset");

        $display("\n[Test] write full 25-channel phase/amplitude frame");
        start_ok = ok_count;
        start_param = param_count;
        send_write_params(seq, 10'd25, 8'd160);
        wait_ok_after(start_ok, "WRITE_PARAMS frame accepted");
        wait_param_after(start_param, "WRITE_PARAMS committed on carrier boundary");
        check(last_cmd == 8'h10, "last successful command is WRITE_PARAMS");
        check(last_seq == seq, "last sequence matches WRITE_PARAMS");
        for (i = 0; i < 25; i = i + 1) begin
            check(dut.u_uart_beam_cmd.phase[i] == expected_phase[i], "active phase matches UART payload");
            check(dut.u_uart_beam_cmd.amplitude[i] == expected_amp[i], "active amplitude matches UART payload");
        end
        seq = seq + 1;

        $display("\n[Test] start ultrasound output");
        start_ok = ok_count;
        send_cmd_no_payload(8'h11, seq);
        wait_ok_after(start_ok, "START frame accepted");
        repeat (4) @(posedge clk);
        check(ultrasound_en == 1'b1, "ultrasound_en asserted by START");
        seen_pwm = 0;
        repeat (3000) begin
            @(posedge clk);
            if (pwm_out != 25'd0)
                seen_pwm = 1;
        end
        check(seen_pwm == 1, "PWM becomes active after START and ADC envelope");
        seq = seq + 1;

        $display("\n[Test] stop ultrasound output");
        start_ok = ok_count;
        send_cmd_no_payload(8'h12, seq);
        wait_ok_after(start_ok, "STOP frame accepted");
        repeat (8) @(posedge clk);
        check(ultrasound_en == 1'b0, "ultrasound_en cleared by STOP");
        check(pwm_out == 25'd0, "PWM outputs forced low after STOP");
        seq = seq + 1;

        $display("\n[Test] update params while stopped, then restart");
        start_ok = ok_count;
        start_param = param_count;
        send_write_params(seq, 10'd300, 8'd120);
        wait_ok_after(start_ok, "second WRITE_PARAMS frame accepted");
        wait_param_after(start_param, "second WRITE_PARAMS committed while stopped");
        check(pwm_out == 25'd0, "PWM stays low while stopped after parameter update");
        for (i = 0; i < 25; i = i + 1) begin
            check(dut.u_uart_beam_cmd.phase[i] == expected_phase[i], "second active phase matches UART payload");
            check(dut.u_uart_beam_cmd.amplitude[i] == expected_amp[i], "second active amplitude matches UART payload");
        end
        seq = seq + 1;

        start_ok = ok_count;
        send_cmd_no_payload(8'h11, seq);
        wait_ok_after(start_ok, "START after stopped update accepted");
        repeat (4) @(posedge clk);
        check(ultrasound_en == 1'b1, "ultrasound_en asserted after restart");
        seq = seq + 1;

        $display("\n[Test] soft reset ultrasound chain");
        start_ok = ok_count;
        start_param = param_count;
        send_cmd_no_payload(8'h13, seq);
        wait_ok_after(start_ok, "SOFT_RESET frame accepted");
        wait_param_after(start_param, "SOFT_RESET reloads active parameters");
        repeat (12) @(posedge clk);
        check(ultrasound_en == 1'b0, "SOFT_RESET leaves ultrasound disabled");
        check(pwm_out == 25'd0, "PWM outputs low after SOFT_RESET");
        seq = seq + 1;

        $display("\n[Test] bad CRC is rejected");
        start_err = err_count;
        start_ok = ok_count;
        send_bad_crc_cmd(8'h11, seq);
        wait_err_after(start_err, "bad CRC frame rejected");
        check(ok_count == start_ok, "bad CRC frame does not increment ok count");
        check(cmd_status[7:4] == 4'd3, "status reports CRC error code");
        check(ultrasound_en == 1'b0, "bad START frame does not enable ultrasound");

        repeat (200) @(posedge clk);

        $display("\n============================================================");
        $display("Simulation summary: passes=%0d errors=%0d", pass_cnt, error_cnt);
        if (error_cnt == 0)
            $display("RESULT: ALL TESTS PASSED");
        else
            $display("RESULT: %0d CHECK(S) FAILED", error_cnt);
        $display("============================================================");

        $finish;
    end

    initial begin
        #5_000_000;
        $display("[FAIL] simulation timeout");
        $finish;
    end

endmodule
