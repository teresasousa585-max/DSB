`timescale 1ns / 1ps

module tb_mcp3201_driver;
    localparam [11:0] ADC_CODE = 12'hA5C;

    reg clk = 1'b0;
    reg rst_n = 1'b0;

    wire adc_cs;
    wire adc_clk;
    tri  adc_data;
    reg  adc_data_drv = 1'bz;

    wire [11:0] sample_data;
    wire        sample_valid;
    wire        sample_busy;

    integer rise_count = 0;
    integer last_rise_count = 0;
    integer fall_count = 0;

    assign adc_data = adc_data_drv;

    always #10 clk = ~clk;

    mcp3201_driver u_dut (
        .clk          (clk),
        .rst_n        (rst_n),
        .adc_cs       (adc_cs),
        .adc_clk      (adc_clk),
        .adc_data     (adc_data),
        .sample_data  (sample_data),
        .sample_valid (sample_valid),
        .sample_busy  (sample_busy)
    );

    always @(negedge adc_cs) begin
        rise_count  = 0;
        fall_count  = 0;
        adc_data_drv = 1'bz;
    end

    always @(posedge adc_clk) begin
        if (!adc_cs)
            rise_count = rise_count + 1;
    end

    always @(posedge adc_cs) begin
        last_rise_count = rise_count;
        adc_data_drv = 1'bz;
    end

    always @(negedge adc_clk) begin
        if (!adc_cs) begin
            fall_count = fall_count + 1;
            case (fall_count)
                1:  adc_data_drv = 1'bz;
                2:  adc_data_drv = 1'b0;
                3:  adc_data_drv = ADC_CODE[11];
                4:  adc_data_drv = ADC_CODE[10];
                5:  adc_data_drv = ADC_CODE[9];
                6:  adc_data_drv = ADC_CODE[8];
                7:  adc_data_drv = ADC_CODE[7];
                8:  adc_data_drv = ADC_CODE[6];
                9:  adc_data_drv = ADC_CODE[5];
                10: adc_data_drv = ADC_CODE[4];
                11: adc_data_drv = ADC_CODE[3];
                12: adc_data_drv = ADC_CODE[2];
                13: adc_data_drv = ADC_CODE[1];
                14: adc_data_drv = ADC_CODE[0];
                15: adc_data_drv = ADC_CODE[1];
                default: adc_data_drv = 1'b0;
            endcase
        end
    end

    always @(posedge clk) begin
        if (sample_valid) begin
            #1;
            if (sample_data !== ADC_CODE) begin
                $display("FAIL: sample_data=%03h expected=%03h", sample_data, ADC_CODE);
                $finish;
            end
            if (last_rise_count != 16) begin
                $display("FAIL: adc_clk rising edges=%0d expected=16", last_rise_count);
                $finish;
            end
            $display("PASS: sample_data=%03h adc_clk rising edges=%0d", sample_data, last_rise_count);
            $finish;
        end
    end

    initial begin
        $dumpfile("tb_mcp3201_driver.vcd");
        $dumpvars(0, tb_mcp3201_driver);

        repeat (5) @(posedge clk);
        rst_n = 1'b1;
    end

    initial begin
        #100000;
        $display("FAIL: timeout");
        $finish;
    end
endmodule
