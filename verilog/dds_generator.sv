`timescale 1ns / 1ps
// 语言: Verilog-2001

module dds_generator (
    input  wire        clk,      // 100MHz 系统时钟
    input  wire        rst_n,
    output reg  signed [15:0] sin_1k  // 输出总线拓宽为 16-bit
);

    // 1kHz FTW = (1000 * 2^32) / 100_000_000 = 42950
    localparam [31:0] FTW_1K = 32'd42950;
    
    reg [31:0] phase_acc;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) phase_acc <= 32'd0;
        else        phase_acc <= phase_acc + FTW_1K;
    end

    // 高 8 位作为全局相位 (0~255)
    wire [7:0] phase_idx = phase_acc[31:24];
    
    // 严格取低 6 位作为 0~63 的查表物理地址，杜绝越界
    wire [5:0] lut_idx = phase_idx[5:0]; 
    
    //寄存器宽度升为 15 位，存储真正的 16 位正弦波第一象限绝对值
    // 公式: round(32767 * sin(pi/2 * i/64))
    reg [14:0] quarter_lut [0:63];
    initial begin
        quarter_lut[0]=0;     quarter_lut[1]=804;   quarter_lut[2]=1608;  quarter_lut[3]=2410;
        quarter_lut[4]=3211;  quarter_lut[5]=4011;  quarter_lut[6]=4807;  quarter_lut[7]=5601;
        quarter_lut[8]=6392;  quarter_lut[9]=7179;  quarter_lut[10]=7961; quarter_lut[11]=8739;
        quarter_lut[12]=9511; quarter_lut[13]=10278;quarter_lut[14]=11038;quarter_lut[15]=11792;
        quarter_lut[16]=12539;quarter_lut[17]=13278;quarter_lut[18]=14009;quarter_lut[19]=14732;
        quarter_lut[20]=15446;quarter_lut[21]=16150;quarter_lut[22]=16845;quarter_lut[23]=17530;
        quarter_lut[24]=18204;quarter_lut[25]=18867;quarter_lut[26]=19519;quarter_lut[27]=20159;
        quarter_lut[28]=20787;quarter_lut[29]=21402;quarter_lut[30]=22004;quarter_lut[31]=22594;
        quarter_lut[32]=23169;quarter_lut[33]=23731;quarter_lut[34]=24278;quarter_lut[35]=24811;
        quarter_lut[36]=25329;quarter_lut[37]=25831;quarter_lut[38]=26318;quarter_lut[39]=26789;
        quarter_lut[40]=27244;quarter_lut[41]=27683;quarter_lut[42]=28105;quarter_lut[43]=28510;
        quarter_lut[44]=28897;quarter_lut[45]=29268;quarter_lut[46]=29621;quarter_lut[47]=29955;
        quarter_lut[48]=30272;quarter_lut[49]=30571;quarter_lut[50]=30851;quarter_lut[51]=31113;
        quarter_lut[52]=31356;quarter_lut[53]=31580;quarter_lut[54]=31785;quarter_lut[55]=31970;
        quarter_lut[56]=32137;quarter_lut[57]=32284;quarter_lut[58]=32412;quarter_lut[59]=32520;
        quarter_lut[60]=32609;quarter_lut[61]=32678;quarter_lut[62]=32727;quarter_lut[63]=32757;
    end

    // phase_idx[6] 控制象限翻转，实现镜像平滑
    // 取出的数据线拓宽到 15 位
    wire [14:0] lut_val = (phase_idx[6]) ? quarter_lut[63 - lut_idx] : quarter_lut[lut_idx];
    
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) sin_1k <= 16'd0;
        else begin
            // phase_idx[7] 决定正负半周
            // {1'b0, lut_val} 会自动拼成 16 位的正数，前面加负号直接完成 16位补码转换
            if (phase_idx[7]) sin_1k <= -{1'b0, lut_val};
            else              sin_1k <=  {1'b0, lut_val};
        end
    end
endmodule