`timescale 1ns / 1ps
// 语言: SystemVerilog
// 功能: I2S 主机接收器 (16-bit, 48kHz)

module i2s_rx_master (
    input  logic        rst_n,
    input  logic        mclk,        // 输入主时钟 12.288MHz

    // 连接到 WM8978 的 I2S 物理引脚
    output logic        bclk,        // 位时钟 (1.536MHz = mclk / 8)
    output logic        lrck,        // 左右声道时钟 (48kHz = bclk / 32)
    input  logic        adcdat,      // ADC 串行数据输入

    // FPGA 内部输出总线
    output logic signed [15:0] left_data,  
    output logic signed [15:0] right_data, 
    output logic               data_valid  // 脉冲: 一帧(左右)接收完毕
);

    // 1. MCLK 分频产生 BCLK (8分频)
    logic [2:0] div_cnt;
    always_ff @(posedge mclk or negedge rst_n) begin
        if (!rst_n) div_cnt <= 3'd0;
        else        div_cnt <= div_cnt + 1'b1;
    end

    // 生成纯净的 BCLK 物理时钟 (占空比 50%)
    always_ff @(posedge mclk or negedge rst_n) begin
        if (!rst_n) bclk <= 1'b0;
        else if (div_cnt == 3'd3) bclk <= 1'b0; // 下降沿
        else if (div_cnt == 3'd7) bclk <= 1'b1; // 上升沿
    end

    // 提取 BCLK 的边沿触发信号，用于内部逻辑同步
    wire bclk_fall = (div_cnt == 3'd3);
    wire bclk_rise = (div_cnt == 3'd7);

    // 2. 位计数器与 LRCK 生成
    logic [4:0] bit_cnt; // 0~31 计数，共 32 个 bit

    always_ff @(posedge mclk or negedge rst_n) begin
        if (!rst_n) begin
            bit_cnt <= 5'd0;
            lrck    <= 1'b0;
        end else if (bclk_fall) begin
            // 在 BCLK 下降沿切换状态，确保建立时间
            bit_cnt <= bit_cnt + 1'b1;
            
            // I2S 规范: 0 为左声道，1 为右声道
            // 在 bit_cnt == 0 时 LRCK 变 0; 在 bit_cnt == 16 时 LRCK 变 1
            if (bit_cnt == 5'd31)      lrck <= 1'b0; 
            else if (bit_cnt == 5'd15) lrck <= 1'b1; 
        end
    end

    // 3. 移位寄存器抓取数据 (处理 1-Bit 延迟)
    logic [15:0] shift_reg;

    always_ff @(posedge mclk or negedge rst_n) begin
        if (!rst_n) begin
            shift_reg  <= 16'd0;
            left_data  <= 16'd0;
            right_data <= 16'd0;
            data_valid <= 1'b0;
        end else begin
            data_valid <= 1'b0; // 默认拉低脉冲
            
            if (bclk_rise) begin
                // 在 BCLK 上升沿采样串行数据
                shift_reg <= {shift_reg[14:0], adcdat};
                // 完美对齐的采样点
                // 因为有 1-bit 延迟，左声道的最后一个 bit(LSB) 会在 bit_cnt == 16 时采样完毕
                if (bit_cnt == 5'd16) begin
                    left_data <= {shift_reg[14:0], adcdat};
                end
                // 右声道的最后一个 bit(LSB) 会在 bit_cnt == 0 时采样完毕
                else if (bit_cnt == 5'd0) begin
                    right_data <= {shift_reg[14:0], adcdat};
                    data_valid <= 1'b1; // 此时一帧数据全部完整，触发有效脉冲
                end
            end
        end
    end

endmodule