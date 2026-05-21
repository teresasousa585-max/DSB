`timescale 1ns / 1ps

module tb_uart_rx_protocol;

    localparam integer CLK_FREQ_HZ = 100_000_000;
    localparam integer BAUD_RATE   = 115200;
    localparam integer BIT_CLKS    = CLK_FREQ_HZ / BAUD_RATE;

    logic clk;
    logic rst_n;
    logic uart_line;
    wire [7:0] uart_byte;
    wire uart_done;
    logic [7:0] amplitude [0:31];
    logic [11:0] phase [0:31];
    logic update_pulse;

    integer pass_cnt;
    integer error_cnt;
    integer uart_done_count;
    integer update_count;
    integer ch;
    reg [7:0] captured_data [0:127];

    uart_rx #(
        .CLK_FRE   (100),
        .BAUD_RATE (BAUD_RATE)
    ) u_uart_rx (
        .clk         (clk),
        .rst_n       (rst_n),
        .i_uart_rx   (uart_line),
        .o_uart_data (uart_byte),
        .o_rx_done   (uart_done)
    );

    uart_protocol_parser u_parser (
        .clk            (clk),
        .rst_n          (rst_n),
        .rx_data        (uart_byte),
        .rx_done        (uart_done),
        .o_amplitude    (amplitude),
        .o_phase        (phase),
        .o_update_pulse (update_pulse)
    );

    initial begin
        clk = 1'b0;
        forever #5 clk = ~clk;
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            uart_done_count <= 0;
            update_count <= 0;
        end else begin
            if (uart_done) begin
                captured_data[uart_done_count] <= uart_byte;
                uart_done_count <= uart_done_count + 1;
            end
            if (update_pulse)
                update_count <= update_count + 1;
        end
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

    task automatic uart_send_byte;
        input [7:0] value;
        integer bit_idx;
        begin
            uart_line = 1'b0;
            repeat (BIT_CLKS) @(posedge clk);

            for (bit_idx = 0; bit_idx < 8; bit_idx = bit_idx + 1) begin
                uart_line = value[bit_idx];
                repeat (BIT_CLKS) @(posedge clk);
            end

            uart_line = 1'b1;
            repeat (BIT_CLKS) @(posedge clk);
        end
    endtask

    task automatic send_n32_frame;
        integer i;
        reg [7:0] sum;
        reg [7:0] amp;
        reg [11:0] ph;
        begin
            sum = 8'h01;
            uart_send_byte(8'hAA);
            uart_send_byte(8'hBB);
            uart_send_byte(8'h01);

            for (i = 0; i < 32; i = i + 1) begin
                amp = i + 8'd1;
                sum = sum + amp;
                uart_send_byte(amp);
            end

            for (i = 0; i < 32; i = i + 1) begin
                ph = (i * 77) & 12'hFFF;
                sum = sum + ph[7:0];
                uart_send_byte(ph[7:0]);
                sum = sum + {4'd0, ph[11:8]};
                uart_send_byte({4'd0, ph[11:8]});
            end

            uart_send_byte(sum);
            uart_send_byte(8'h0D);
            uart_send_byte(8'h0A);
        end
    endtask

    initial begin
        $dumpfile("tb_uart_rx_protocol.vcd");
        $dumpvars(0, tb_uart_rx_protocol);

        pass_cnt = 0;
        error_cnt = 0;
        uart_line = 1'b1;
        rst_n = 1'b0;

        repeat (20) @(posedge clk);
        rst_n = 1'b1;
        repeat (20) @(posedge clk);

        send_n32_frame();
        repeat (BIT_CLKS * 3) @(posedge clk);

        check(uart_done_count == 102, "uart_rx emits one done pulse for each frame byte");
        check(update_count == 1, "parser emits one update pulse after serial N32 frame");

        for (ch = 0; ch < 32; ch = ch + 1) begin
            check(u_parser.o_amplitude[ch] == ch + 1, "amplitude matches serial frame payload");
            check(u_parser.o_phase[ch] == ((ch * 77) & 12'hFFF), "phase matches serial frame payload");
        end

        $display("Simulation summary: passes=%0d errors=%0d uart_done_count=%0d update_count=%0d",
                 pass_cnt, error_cnt, uart_done_count, update_count);
        if (error_cnt == 0)
            $display("RESULT: ALL TESTS PASSED");
        else begin
            $display("RESULT: %0d CHECK(S) FAILED", error_cnt);
            for (ch = 0; ch < uart_done_count && ch < 32; ch = ch + 1)
                $display("captured[%0d]=0x%02h", ch, captured_data[ch]);
        end
        $finish;
    end

    initial begin
        #12000000;
        $display("[FAIL] simulation timeout");
        $finish;
    end

endmodule
