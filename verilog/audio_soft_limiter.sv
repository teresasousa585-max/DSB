`timescale 1ns / 1ps

module audio_soft_limiter (
    input  logic               clk,
    input  logic               rst_n,
    input  logic signed [21:0] audio_in,   // 放大后的输入信号 (*8 之后)
    output logic signed [15:0] audio_out   // 软限幅后的输出
);

    // ==========================================
    // Stage 1: 取绝对值与符号位提取
    // ==========================================
    logic        sign_d1;
    logic [21:0] abs_val_d1;

    always_ff @(posedge clk) begin
        sign_d1 <= audio_in[21]; 
        if (audio_in[21]) 
            abs_val_d1 <= -audio_in; 
        else 
            abs_val_d1 <= audio_in;
    end

    // ==========================================
    // Stage 2: 钳位寻址与符号打拍
    // ==========================================
    logic [14:0] rom_addr_d2;
    logic        sign_d2;

    always_ff @(posedge clk) begin
        sign_d2 <= sign_d1;
        // 防呆保护：极限超量程时，死锁在最高地址
        if (abs_val_d1 > 22'd32767)
            rom_addr_d2 <= 15'd32767;
        else
            rom_addr_d2 <= abs_val_d1[14:0];
    end

    // ==========================================
    // Stage 3: 例化 Vivado ROM IP 核
    // ==========================================
    logic [15:0] rom_data_d3;
    
    soft_knee_rom_ip u_soft_limiter_rom (
        .clka  (clk),           // 输入时钟
        .addra (rom_addr_d2),   // 15位查表地址
        .douta (rom_data_d3)    // 16位平滑后的数据输出 (固定延迟 1 拍或 2 拍，取决于你的 IP 配置)
    );

    // 维持符号位流水线，与 ROM 输出对齐
    // 注意：如果你的 ROM IP 勾选了 Output Register，这里可能需要打两拍 (sign_d3_reg)
    logic sign_d3;
    always_ff @(posedge clk) begin
        sign_d3 <= sign_d2;
    end

    // ==========================================
    // Stage 4: 恢复符号位，输出最终音频
    // ==========================================
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            audio_out <= 16'sd0;
        end else begin
            if (sign_d3)
                audio_out <= -rom_data_d3; // 负数
            else
                audio_out <= rom_data_d3;  // 正数
        end
    end

endmodule