`timescale 1ns / 1ps

module tb_uart_protocol_parser;

    logic clk;
    logic rst_n;
    logic [7:0] rx_data;
    logic rx_done;
    logic [7:0] amplitude [0:31];
    logic [11:0] phase [0:31];
    logic update_pulse;

    integer pass_cnt;
    integer error_cnt;
    integer update_count;
    integer ch;

    uart_protocol_parser dut (
        .clk            (clk),
        .rst_n          (rst_n),
        .rx_data        (rx_data),
        .rx_done        (rx_done),
        .o_amplitude    (amplitude),
        .o_phase        (phase),
        .o_update_pulse (update_pulse)
    );

    initial begin
        clk = 1'b0;
        forever #5 clk = ~clk;
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n)
            update_count <= 0;
        else if (update_pulse)
            update_count <= update_count + 1;
    end

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

    task automatic send_byte;
        input [7:0] value;
        begin
            @(negedge clk);
            rx_data = value;
            rx_done = 1'b1;
            @(negedge clk);
            rx_done = 1'b0;
            repeat (2) @(negedge clk);
        end
    endtask

    task automatic send_n32_frame;
        input bad_checksum;
        input bad_tail;
        integer i;
        reg [7:0] sum;
        reg [11:0] ph;
        reg [7:0] amp;
        begin
            sum = 8'h01;
            send_byte(8'hAA);
            send_byte(8'hBB);
            send_byte(8'h01);

            for (i = 0; i < 32; i = i + 1) begin
                amp = i + 8'd1;
                sum = sum + amp;
                send_byte(amp);
            end

            for (i = 0; i < 32; i = i + 1) begin
                ph = (i * 77) & 12'hFFF;
                sum = sum + ph[7:0];
                send_byte(ph[7:0]);
                sum = sum + {4'd0, ph[11:8]};
                send_byte({4'd0, ph[11:8]});
            end

            send_byte(bad_checksum ? (sum ^ 8'h55) : sum);
            send_byte(bad_tail ? 8'h00 : 8'h0D);
            send_byte(8'h0A);
        end
    endtask

    task automatic send_partial_frame;
        integer i;
        begin
            send_byte(8'hAA);
            send_byte(8'hBB);
            send_byte(8'h01);
            for (i = 0; i < 5; i = i + 1)
                send_byte(i[7:0]);
        end
    endtask

    initial begin
        $dumpfile("tb_uart_protocol_parser.vcd");
        $dumpvars(0, tb_uart_protocol_parser);

        pass_cnt = 0;
        error_cnt = 0;
        rx_data = 8'd0;
        rx_done = 1'b0;
        rst_n = 1'b0;
        repeat (5) @(posedge clk);
        rst_n = 1'b1;
        repeat (3) @(posedge clk);

        check(dut.o_amplitude[0] == 8'd255, "reset default amplitude is full scale");
        check(dut.o_phase[0] == 12'd0, "reset default phase is zero");

        send_n32_frame(1'b0, 1'b0);
        repeat (5) @(posedge clk);
        check(update_count == 1, "valid N32 frame produces one update pulse");
        for (ch = 0; ch < 32; ch = ch + 1) begin
            check(dut.o_amplitude[ch] == ch + 1, "amplitude matches payload");
            check(dut.o_phase[ch] == ((ch * 77) & 12'hFFF), "phase matches payload");
        end

        send_n32_frame(1'b1, 1'b0);
        repeat (5) @(posedge clk);
        check(update_count == 1, "bad checksum does not update");

        send_n32_frame(1'b0, 1'b1);
        repeat (5) @(posedge clk);
        check(update_count == 1, "bad tail does not update");

        rst_n = 1'b0;
        repeat (3) @(posedge clk);
        rst_n = 1'b1;
        repeat (3) @(posedge clk);
        send_partial_frame();
        repeat (20) @(posedge clk);
        check(update_count == 0, "short packet does not update");

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
