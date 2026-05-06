//////////////////////////////////////////////////////////////////////////////
// 模块名称  : beam_controller
// 功能描述  : 5x5超声阵列波束控制器
//             计算25路换能器的相位和幅度参数，支持波束偏转、聚焦和空域加窗
// 作者      : Auto-generated
// 版本      : V1.0
// 日期      : 2024
//////////////////////////////////////////////////////////////////////////////
//
// 阵列几何参数：
//   - 5x5方阵，阵元间距 d = 16mm
//   - 阵列中心在 (2,2)
//   - 超声波长 lambda = 8.575mm (f=40kHz, c=343m/s)
//   - d/lambda = 1.866
//
// 相位量化：
//   - 10位表示 0~1023 对应 0~360度
//   - 相位分辨率 = 360/1024 ~ 0.35度
//
// 幅度量化：
//   - 8位表示 0~255
//   - 255 对应满幅度 1.0
//
//////////////////////////////////////////////////////////////////////////////

module beam_controller (
    input        clk,       // 50MHz系统时钟
    input        rst_n,     // 低电平复位（异步复位同步释放）

    // 控制模式
    input [2:0]  mode,      // 模式选择
                            // 000: 同相驱动（全0相位，全幅度）
                            // 001: 波束偏转（单角度偏转）
                            // 010: 波束聚焦（单距离聚焦）
                            // 011: 偏转+聚焦
                            // 100: 同相+加窗（抑制栅瓣）

    // 控制参数
    input signed [7:0]  steer_angle,    // 偏转角度（度，-30~+30，有符号）
    input [12:0]        focus_dist_mm,  // 聚焦距离（mm，100~5000）
    input [2:0]         window_type,    // 窗函数类型
                                        // 0: 矩形窗（全1）
                                        // 1: 汉宁窗
                                        // 2: 汉明窗
                                        // 3: 布莱克曼窗

    // 参数更新触发
    input               param_update,   // 高电平脉冲更新参数

    // 输出（25路参数）
    output reg [7:0]  amplitude [0:24], // 25路幅度 0~255
    output reg [9:0]  phase     [0:24], // 25路相位 0~1023
    output reg        param_valid       // 参数有效标志
);

    //=========================================================================
    // 参数定义
    //=========================================================================
    parameter NUM_ELEMENTS       = 25;      // 阵元总数 5x5
    parameter PHASE_BITS         = 10;      // 相位位数
    parameter AMP_BITS           = 8;       // 幅度位数
    parameter PHASE_MAX          = 1023;    // 相位最大值 (2^10 - 1)
    parameter AMP_MAX            = 255;     // 幅度最大值 (2^8 - 1)

    // d/lambda = 16/8.575 = 1.866
    // phase_increment = 1024 * d * sin(theta) / lambda
    //                 = 1024 * 1.866 * sin(theta)

    // 聚焦LUT参数: 距离100~5000mm，步进50mm，共99个点
    parameter FOCUS_LUT_STEP     = 50;
    parameter FOCUS_LUT_POINTS   = 99;

    //=========================================================================
    // 状态机定义
    //=========================================================================
    localparam IDLE     = 3'd0;
    localparam LOAD     = 3'd1;
    localparam CALC     = 3'd2;
    localparam DONE_S   = 3'd3;

    //=========================================================================
    // 内部寄存器
    //=========================================================================
    reg [2:0]  state;               // 状态机
    reg [2:0]  mode_reg;            // 锁存模式
    reg signed [7:0] steer_reg;     // 锁存偏转角度
    reg [12:0] focus_dist_reg;      // 锁存聚焦距离
    reg [2:0]  window_reg;          // 锁存窗函数类型
    reg [4:0]  calc_idx;            // 计算索引 0~24
    reg signed [9:0] steer_phase_incr;  // 偏转相位增量

    //=========================================================================
    // 窗函数LUT（5x5，8位定点，255=1.0）
    // 索引：idx = m*5 + n，m,n in {0,1,2,3,4}
    //=========================================================================

    // 矩形窗（全1）
    wire [7:0] window_rect [0:24];
    genvar gi;
    generate
        for (gi = 0; gi < 25; gi = gi + 1) begin : gen_rect
            assign window_rect[gi] = 8'd255;
        end
    endgenerate

    // 汉宁窗 LUT
    wire [7:0] window_hann [0:24];
    assign window_hann[0]  = 8'd0;    assign window_hann[1]  = 8'd16;
    assign window_hann[2]  = 8'd32;   assign window_hann[3]  = 8'd16;
    assign window_hann[4]  = 8'd0;    assign window_hann[5]  = 8'd16;
    assign window_hann[6]  = 8'd64;   assign window_hann[7]  = 8'd95;
    assign window_hann[8]  = 8'd64;   assign window_hann[9]  = 8'd16;
    assign window_hann[10] = 8'd32;   assign window_hann[11] = 8'd95;
    assign window_hann[12] = 8'd143;  assign window_hann[13] = 8'd95;
    assign window_hann[14] = 8'd32;   assign window_hann[15] = 8'd16;
    assign window_hann[16] = 8'd64;   assign window_hann[17] = 8'd95;
    assign window_hann[18] = 8'd64;   assign window_hann[19] = 8'd16;
    assign window_hann[20] = 8'd0;    assign window_hann[21] = 8'd16;
    assign window_hann[22] = 8'd32;   assign window_hann[23] = 8'd16;
    assign window_hann[24] = 8'd0;

    // 汉明窗 LUT (alpha=0.54)
    wire [7:0] window_hamm [0:24];
    assign window_hamm[0]  = 8'd6;    assign window_hamm[1]  = 8'd23;
    assign window_hamm[2]  = 8'd35;   assign window_hamm[3]  = 8'd23;
    assign window_hamm[4]  = 8'd6;    assign window_hamm[5]  = 8'd23;
    assign window_hamm[6]  = 8'd83;   assign window_hamm[7]  = 8'd129;
    assign window_hamm[8]  = 8'd83;   assign window_hamm[9]  = 8'd23;
    assign window_hamm[10] = 8'd35;   assign window_hamm[11] = 8'd129;
    assign window_hamm[12] = 8'd255;  assign window_hamm[13] = 8'd129;
    assign window_hamm[14] = 8'd35;   assign window_hamm[15] = 8'd23;
    assign window_hamm[16] = 8'd83;   assign window_hamm[17] = 8'd129;
    assign window_hamm[18] = 8'd83;   assign window_hamm[19] = 8'd23;
    assign window_hamm[20] = 8'd6;    assign window_hamm[21] = 8'd23;
    assign window_hamm[22] = 8'd35;   assign window_hamm[23] = 8'd23;
    assign window_hamm[24] = 8'd6;

    // 布莱克曼窗 LUT
    wire [7:0] window_blk [0:24];
    assign window_blk[0]  = 8'd0;     assign window_blk[1]  = 8'd4;
    assign window_blk[2]  = 8'd7;     assign window_blk[3]  = 8'd4;
    assign window_blk[4]  = 8'd0;     assign window_blk[5]  = 8'd4;
    assign window_blk[6]  = 8'd32;    assign window_blk[7]  = 8'd57;
    assign window_blk[8]  = 8'd32;    assign window_blk[9]  = 8'd4;
    assign window_blk[10] = 8'd7;     assign window_blk[11] = 8'd57;
    assign window_blk[12] = 8'd255;   assign window_blk[13] = 8'd57;
    assign window_blk[14] = 8'd7;     assign window_blk[15] = 8'd4;
    assign window_blk[16] = 8'd32;    assign window_blk[17] = 8'd57;
    assign window_blk[18] = 8'd32;    assign window_blk[19] = 8'd4;
    assign window_blk[20] = 8'd0;     assign window_blk[21] = 8'd4;
    assign window_blk[22] = 8'd7;     assign window_blk[23] = 8'd4;
    assign window_blk[24] = 8'd0;

    //=========================================================================
    // 阵元行位置LUT（相对中心，单位：16mm步进）
    // m=0:-2, m=1:-1, m=2:0, m=3:+1, m=4:+2
    //=========================================================================
    wire signed [2:0] element_row [0:24];
    assign element_row[0]  = -3'd2;  assign element_row[1]  = -3'd2;
    assign element_row[2]  = -3'd2;  assign element_row[3]  = -3'd2;
    assign element_row[4]  = -3'd2;  assign element_row[5]  = -3'd1;
    assign element_row[6]  = -3'd1;  assign element_row[7]  = -3'd1;
    assign element_row[8]  = -3'd1;  assign element_row[9]  = -3'd1;
    assign element_row[10] =  3'd0;  assign element_row[11] =  3'd0;
    assign element_row[12] =  3'd0;  assign element_row[13] =  3'd0;
    assign element_row[14] =  3'd0;  assign element_row[15] =  3'd1;
    assign element_row[16] =  3'd1;  assign element_row[17] =  3'd1;
    assign element_row[18] =  3'd1;  assign element_row[19] =  3'd1;
    assign element_row[20] =  3'd2;  assign element_row[21] =  3'd2;
    assign element_row[22] =  3'd2;  assign element_row[23] =  3'd2;
    assign element_row[24] =  3'd2;

    //=========================================================================
    // 阵元列位置LUT（相对中心，单位：16mm步进）
    // n=0:-2, n=1:-1, n=2:0, n=3:+1, n=4:+2
    //=========================================================================
    wire signed [2:0] element_col [0:24];
    assign element_col[0]  = -3'd2;  assign element_col[1]  = -3'd1;
    assign element_col[2]  =  3'd0;  assign element_col[3]  =  3'd1;
    assign element_col[4]  =  3'd2;  assign element_col[5]  = -3'd2;
    assign element_col[6]  = -3'd1;  assign element_col[7]  =  3'd0;
    assign element_col[8]  =  3'd1;  assign element_col[9]  =  3'd2;
    assign element_col[10] = -3'd2;  assign element_col[11] = -3'd1;
    assign element_col[12] =  3'd0;  assign element_col[13] =  3'd1;
    assign element_col[14] =  3'd2;  assign element_col[15] = -3'd2;
    assign element_col[16] = -3'd1;  assign element_col[17] =  3'd0;
    assign element_col[18] =  3'd1;  assign element_col[19] =  3'd2;
    assign element_col[20] = -3'd2;  assign element_col[21] = -3'd1;
    assign element_col[22] =  3'd0;  assign element_col[23] =  3'd1;
    assign element_col[24] =  3'd2;

    //=========================================================================
    // 偏转相位增量LUT
    // 相邻阵元(16mm间距)的相位差 = 1024 * d * sin(theta) / lambda
    // 角度范围 -30~+30度，步进1度，共61个点
    // LUT索引 = steer_angle + 30 (映射-30~+30到0~60)
    //=========================================================================
    reg signed [9:0] steer_phase_lut [0:60];
    initial begin
        steer_phase_lut[ 0] = -10'sd955;  // theta=-30
        steer_phase_lut[ 1] = -10'sd926;  // theta=-29
        steer_phase_lut[ 2] = -10'sd897;  // theta=-28
        steer_phase_lut[ 3] = -10'sd867;  // theta=-27
        steer_phase_lut[ 4] = -10'sd838;  // theta=-26
        steer_phase_lut[ 5] = -10'sd807;  // theta=-25
        steer_phase_lut[ 6] = -10'sd777;  // theta=-24
        steer_phase_lut[ 7] = -10'sd747;  // theta=-23
        steer_phase_lut[ 8] = -10'sd716;  // theta=-22
        steer_phase_lut[ 9] = -10'sd685;  // theta=-21
        steer_phase_lut[10] = -10'sd653;  // theta=-20
        steer_phase_lut[11] = -10'sd622;  // theta=-19
        steer_phase_lut[12] = -10'sd590;  // theta=-18
        steer_phase_lut[13] = -10'sd559;  // theta=-17
        steer_phase_lut[14] = -10'sd527;  // theta=-16
        steer_phase_lut[15] = -10'sd495;  // theta=-15
        steer_phase_lut[16] = -10'sd462;  // theta=-14
        steer_phase_lut[17] = -10'sd430;  // theta=-13
        steer_phase_lut[18] = -10'sd397;  // theta=-12
        steer_phase_lut[19] = -10'sd365;  // theta=-11
        steer_phase_lut[20] = -10'sd332;  // theta=-10
        steer_phase_lut[21] = -10'sd299;  // theta= -9
        steer_phase_lut[22] = -10'sd266;  // theta= -8
        steer_phase_lut[23] = -10'sd233;  // theta= -7
        steer_phase_lut[24] = -10'sd200;  // theta= -6
        steer_phase_lut[25] = -10'sd167;  // theta= -5
        steer_phase_lut[26] = -10'sd133;  // theta= -4
        steer_phase_lut[27] = -10'sd100;  // theta= -3
        steer_phase_lut[28] = -10'sd67;  // theta= -2
        steer_phase_lut[29] = -10'sd33;  // theta= -1
        steer_phase_lut[30] =  10'sd0;   // theta= +0
        steer_phase_lut[31] =  10'sd33;   // theta= +1
        steer_phase_lut[32] =  10'sd67;   // theta= +2
        steer_phase_lut[33] =  10'sd100;   // theta= +3
        steer_phase_lut[34] =  10'sd133;   // theta= +4
        steer_phase_lut[35] =  10'sd167;   // theta= +5
        steer_phase_lut[36] =  10'sd200;   // theta= +6
        steer_phase_lut[37] =  10'sd233;   // theta= +7
        steer_phase_lut[38] =  10'sd266;   // theta= +8
        steer_phase_lut[39] =  10'sd299;   // theta= +9
        steer_phase_lut[40] =  10'sd332;   // theta=+10
        steer_phase_lut[41] =  10'sd365;   // theta=+11
        steer_phase_lut[42] =  10'sd397;   // theta=+12
        steer_phase_lut[43] =  10'sd430;   // theta=+13
        steer_phase_lut[44] =  10'sd462;   // theta=+14
        steer_phase_lut[45] =  10'sd495;   // theta=+15
        steer_phase_lut[46] =  10'sd527;   // theta=+16
        steer_phase_lut[47] =  10'sd559;   // theta=+17
        steer_phase_lut[48] =  10'sd590;   // theta=+18
        steer_phase_lut[49] =  10'sd622;   // theta=+19
        steer_phase_lut[50] =  10'sd653;   // theta=+20
        steer_phase_lut[51] =  10'sd685;   // theta=+21
        steer_phase_lut[52] =  10'sd716;   // theta=+22
        steer_phase_lut[53] =  10'sd747;   // theta=+23
        steer_phase_lut[54] =  10'sd777;   // theta=+24
        steer_phase_lut[55] =  10'sd807;   // theta=+25
        steer_phase_lut[56] =  10'sd838;   // theta=+26
        steer_phase_lut[57] =  10'sd867;   // theta=+27
        steer_phase_lut[58] =  10'sd897;   // theta=+28
        steer_phase_lut[59] =  10'sd926;   // theta=+29
        steer_phase_lut[60] =  10'sd955;   // theta=+30
    end

    //=========================================================================
    // 聚焦相位LUT
    // 聚焦距离100~5000mm，步进50mm，共99个点
    // LUT索引 = (focus_dist_mm - 100) / 50
    //=========================================================================
    reg [9:0] focus_phase_lut [0:98] [0:24];
    initial begin
        // focus_dist=100mm (idx=0)
        focus_phase_lut[ 0][ 0] = 10'd142;
        focus_phase_lut[ 0][ 1] = 10'd741;
        focus_phase_lut[ 0][ 2] = 10'd597;
        focus_phase_lut[ 0][ 3] = 10'd741;
        focus_phase_lut[ 0][ 4] = 10'd142;
        focus_phase_lut[ 0][ 5] = 10'd741;
        focus_phase_lut[ 0][ 6] = 10'd302;
        focus_phase_lut[ 0][ 7] = 10'd152;
        focus_phase_lut[ 0][ 8] = 10'd302;
        focus_phase_lut[ 0][ 9] = 10'd741;
        focus_phase_lut[ 0][10] = 10'd597;
        focus_phase_lut[ 0][11] = 10'd152;
        focus_phase_lut[ 0][12] = 10'd0;
        focus_phase_lut[ 0][13] = 10'd152;
        focus_phase_lut[ 0][14] = 10'd597;
        focus_phase_lut[ 0][15] = 10'd741;
        focus_phase_lut[ 0][16] = 10'd302;
        focus_phase_lut[ 0][17] = 10'd152;
        focus_phase_lut[ 0][18] = 10'd302;
        focus_phase_lut[ 0][19] = 10'd741;
        focus_phase_lut[ 0][20] = 10'd142;
        focus_phase_lut[ 0][21] = 10'd741;
        focus_phase_lut[ 0][22] = 10'd597;
        focus_phase_lut[ 0][23] = 10'd741;
        focus_phase_lut[ 0][24] = 10'd142;
        // focus_dist=150mm (idx=1)
        focus_phase_lut[ 1][ 0] = 10'd797;
        focus_phase_lut[ 1][ 1] = 10'd502;
        focus_phase_lut[ 1][ 2] = 10'd403;
        focus_phase_lut[ 1][ 3] = 10'd502;
        focus_phase_lut[ 1][ 4] = 10'd797;
        focus_phase_lut[ 1][ 5] = 10'd502;
        focus_phase_lut[ 1][ 6] = 10'd203;
        focus_phase_lut[ 1][ 7] = 10'd102;
        focus_phase_lut[ 1][ 8] = 10'd203;
        focus_phase_lut[ 1][ 9] = 10'd502;
        focus_phase_lut[ 1][10] = 10'd403;
        focus_phase_lut[ 1][11] = 10'd102;
        focus_phase_lut[ 1][12] = 10'd0;
        focus_phase_lut[ 1][13] = 10'd102;
        focus_phase_lut[ 1][14] = 10'd403;
        focus_phase_lut[ 1][15] = 10'd502;
        focus_phase_lut[ 1][16] = 10'd203;
        focus_phase_lut[ 1][17] = 10'd102;
        focus_phase_lut[ 1][18] = 10'd203;
        focus_phase_lut[ 1][19] = 10'd502;
        focus_phase_lut[ 1][20] = 10'd797;
        focus_phase_lut[ 1][21] = 10'd502;
        focus_phase_lut[ 1][22] = 10'd403;
        focus_phase_lut[ 1][23] = 10'd502;
        focus_phase_lut[ 1][24] = 10'd797;
        // focus_dist=200mm (idx=2)
        focus_phase_lut[ 2][ 0] = 10'd604;
        focus_phase_lut[ 2][ 1] = 10'd379;
        focus_phase_lut[ 2][ 2] = 10'd304;
        focus_phase_lut[ 2][ 3] = 10'd379;
        focus_phase_lut[ 2][ 4] = 10'd604;
        focus_phase_lut[ 2][ 5] = 10'd379;
        focus_phase_lut[ 2][ 6] = 10'd152;
        focus_phase_lut[ 2][ 7] = 10'd76;
        focus_phase_lut[ 2][ 8] = 10'd152;
        focus_phase_lut[ 2][ 9] = 10'd379;
        focus_phase_lut[ 2][10] = 10'd304;
        focus_phase_lut[ 2][11] = 10'd76;
        focus_phase_lut[ 2][12] = 10'd0;
        focus_phase_lut[ 2][13] = 10'd76;
        focus_phase_lut[ 2][14] = 10'd304;
        focus_phase_lut[ 2][15] = 10'd379;
        focus_phase_lut[ 2][16] = 10'd152;
        focus_phase_lut[ 2][17] = 10'd76;
        focus_phase_lut[ 2][18] = 10'd152;
        focus_phase_lut[ 2][19] = 10'd379;
        focus_phase_lut[ 2][20] = 10'd604;
        focus_phase_lut[ 2][21] = 10'd379;
        focus_phase_lut[ 2][22] = 10'd304;
        focus_phase_lut[ 2][23] = 10'd379;
        focus_phase_lut[ 2][24] = 10'd604;
        // focus_dist=250mm (idx=3)
        focus_phase_lut[ 3][ 0] = 10'd485;
        focus_phase_lut[ 3][ 1] = 10'd304;
        focus_phase_lut[ 3][ 2] = 10'd244;
        focus_phase_lut[ 3][ 3] = 10'd304;
        focus_phase_lut[ 3][ 4] = 10'd485;
        focus_phase_lut[ 3][ 5] = 10'd304;
        focus_phase_lut[ 3][ 6] = 10'd122;
        focus_phase_lut[ 3][ 7] = 10'd61;
        focus_phase_lut[ 3][ 8] = 10'd122;
        focus_phase_lut[ 3][ 9] = 10'd304;
        focus_phase_lut[ 3][10] = 10'd244;
        focus_phase_lut[ 3][11] = 10'd61;
        focus_phase_lut[ 3][12] = 10'd0;
        focus_phase_lut[ 3][13] = 10'd61;
        focus_phase_lut[ 3][14] = 10'd244;
        focus_phase_lut[ 3][15] = 10'd304;
        focus_phase_lut[ 3][16] = 10'd122;
        focus_phase_lut[ 3][17] = 10'd61;
        focus_phase_lut[ 3][18] = 10'd122;
        focus_phase_lut[ 3][19] = 10'd304;
        focus_phase_lut[ 3][20] = 10'd485;
        focus_phase_lut[ 3][21] = 10'd304;
        focus_phase_lut[ 3][22] = 10'd244;
        focus_phase_lut[ 3][23] = 10'd304;
        focus_phase_lut[ 3][24] = 10'd485;
        // focus_dist=300mm (idx=4)
        focus_phase_lut[ 4][ 0] = 10'd405;
        focus_phase_lut[ 4][ 1] = 10'd254;
        focus_phase_lut[ 4][ 2] = 10'd203;
        focus_phase_lut[ 4][ 3] = 10'd254;
        focus_phase_lut[ 4][ 4] = 10'd405;
        focus_phase_lut[ 4][ 5] = 10'd254;
        focus_phase_lut[ 4][ 6] = 10'd102;
        focus_phase_lut[ 4][ 7] = 10'd51;
        focus_phase_lut[ 4][ 8] = 10'd102;
        focus_phase_lut[ 4][ 9] = 10'd254;
        focus_phase_lut[ 4][10] = 10'd203;
        focus_phase_lut[ 4][11] = 10'd51;
        focus_phase_lut[ 4][12] = 10'd0;
        focus_phase_lut[ 4][13] = 10'd51;
        focus_phase_lut[ 4][14] = 10'd203;
        focus_phase_lut[ 4][15] = 10'd254;
        focus_phase_lut[ 4][16] = 10'd102;
        focus_phase_lut[ 4][17] = 10'd51;
        focus_phase_lut[ 4][18] = 10'd102;
        focus_phase_lut[ 4][19] = 10'd254;
        focus_phase_lut[ 4][20] = 10'd405;
        focus_phase_lut[ 4][21] = 10'd254;
        focus_phase_lut[ 4][22] = 10'd203;
        focus_phase_lut[ 4][23] = 10'd254;
        focus_phase_lut[ 4][24] = 10'd405;
        // focus_dist=350mm (idx=5)
        focus_phase_lut[ 5][ 0] = 10'd348;
        focus_phase_lut[ 5][ 1] = 10'd218;
        focus_phase_lut[ 5][ 2] = 10'd174;
        focus_phase_lut[ 5][ 3] = 10'd218;
        focus_phase_lut[ 5][ 4] = 10'd348;
        focus_phase_lut[ 5][ 5] = 10'd218;
        focus_phase_lut[ 5][ 6] = 10'd87;
        focus_phase_lut[ 5][ 7] = 10'd44;
        focus_phase_lut[ 5][ 8] = 10'd87;
        focus_phase_lut[ 5][ 9] = 10'd218;
        focus_phase_lut[ 5][10] = 10'd174;
        focus_phase_lut[ 5][11] = 10'd44;
        focus_phase_lut[ 5][12] = 10'd0;
        focus_phase_lut[ 5][13] = 10'd44;
        focus_phase_lut[ 5][14] = 10'd174;
        focus_phase_lut[ 5][15] = 10'd218;
        focus_phase_lut[ 5][16] = 10'd87;
        focus_phase_lut[ 5][17] = 10'd44;
        focus_phase_lut[ 5][18] = 10'd87;
        focus_phase_lut[ 5][19] = 10'd218;
        focus_phase_lut[ 5][20] = 10'd348;
        focus_phase_lut[ 5][21] = 10'd218;
        focus_phase_lut[ 5][22] = 10'd174;
        focus_phase_lut[ 5][23] = 10'd218;
        focus_phase_lut[ 5][24] = 10'd348;
        // focus_dist=400mm (idx=6)
        focus_phase_lut[ 6][ 0] = 10'd305;
        focus_phase_lut[ 6][ 1] = 10'd191;
        focus_phase_lut[ 6][ 2] = 10'd153;
        focus_phase_lut[ 6][ 3] = 10'd191;
        focus_phase_lut[ 6][ 4] = 10'd305;
        focus_phase_lut[ 6][ 5] = 10'd191;
        focus_phase_lut[ 6][ 6] = 10'd76;
        focus_phase_lut[ 6][ 7] = 10'd38;
        focus_phase_lut[ 6][ 8] = 10'd76;
        focus_phase_lut[ 6][ 9] = 10'd191;
        focus_phase_lut[ 6][10] = 10'd153;
        focus_phase_lut[ 6][11] = 10'd38;
        focus_phase_lut[ 6][12] = 10'd0;
        focus_phase_lut[ 6][13] = 10'd38;
        focus_phase_lut[ 6][14] = 10'd153;
        focus_phase_lut[ 6][15] = 10'd191;
        focus_phase_lut[ 6][16] = 10'd76;
        focus_phase_lut[ 6][17] = 10'd38;
        focus_phase_lut[ 6][18] = 10'd76;
        focus_phase_lut[ 6][19] = 10'd191;
        focus_phase_lut[ 6][20] = 10'd305;
        focus_phase_lut[ 6][21] = 10'd191;
        focus_phase_lut[ 6][22] = 10'd153;
        focus_phase_lut[ 6][23] = 10'd191;
        focus_phase_lut[ 6][24] = 10'd305;
        // focus_dist=450mm (idx=7)
        focus_phase_lut[ 7][ 0] = 10'd271;
        focus_phase_lut[ 7][ 1] = 10'd170;
        focus_phase_lut[ 7][ 2] = 10'd136;
        focus_phase_lut[ 7][ 3] = 10'd170;
        focus_phase_lut[ 7][ 4] = 10'd271;
        focus_phase_lut[ 7][ 5] = 10'd170;
        focus_phase_lut[ 7][ 6] = 10'd68;
        focus_phase_lut[ 7][ 7] = 10'd34;
        focus_phase_lut[ 7][ 8] = 10'd68;
        focus_phase_lut[ 7][ 9] = 10'd170;
        focus_phase_lut[ 7][10] = 10'd136;
        focus_phase_lut[ 7][11] = 10'd34;
        focus_phase_lut[ 7][12] = 10'd0;
        focus_phase_lut[ 7][13] = 10'd34;
        focus_phase_lut[ 7][14] = 10'd136;
        focus_phase_lut[ 7][15] = 10'd170;
        focus_phase_lut[ 7][16] = 10'd68;
        focus_phase_lut[ 7][17] = 10'd34;
        focus_phase_lut[ 7][18] = 10'd68;
        focus_phase_lut[ 7][19] = 10'd170;
        focus_phase_lut[ 7][20] = 10'd271;
        focus_phase_lut[ 7][21] = 10'd170;
        focus_phase_lut[ 7][22] = 10'd136;
        focus_phase_lut[ 7][23] = 10'd170;
        focus_phase_lut[ 7][24] = 10'd271;
        // focus_dist=500mm (idx=8)
        focus_phase_lut[ 8][ 0] = 10'd244;
        focus_phase_lut[ 8][ 1] = 10'd153;
        focus_phase_lut[ 8][ 2] = 10'd122;
        focus_phase_lut[ 8][ 3] = 10'd153;
        focus_phase_lut[ 8][ 4] = 10'd244;
        focus_phase_lut[ 8][ 5] = 10'd153;
        focus_phase_lut[ 8][ 6] = 10'd61;
        focus_phase_lut[ 8][ 7] = 10'd31;
        focus_phase_lut[ 8][ 8] = 10'd61;
        focus_phase_lut[ 8][ 9] = 10'd153;
        focus_phase_lut[ 8][10] = 10'd122;
        focus_phase_lut[ 8][11] = 10'd31;
        focus_phase_lut[ 8][12] = 10'd0;
        focus_phase_lut[ 8][13] = 10'd31;
        focus_phase_lut[ 8][14] = 10'd122;
        focus_phase_lut[ 8][15] = 10'd153;
        focus_phase_lut[ 8][16] = 10'd61;
        focus_phase_lut[ 8][17] = 10'd31;
        focus_phase_lut[ 8][18] = 10'd61;
        focus_phase_lut[ 8][19] = 10'd153;
        focus_phase_lut[ 8][20] = 10'd244;
        focus_phase_lut[ 8][21] = 10'd153;
        focus_phase_lut[ 8][22] = 10'd122;
        focus_phase_lut[ 8][23] = 10'd153;
        focus_phase_lut[ 8][24] = 10'd244;
        // focus_dist=550mm (idx=9)
        focus_phase_lut[ 9][ 0] = 10'd222;
        focus_phase_lut[ 9][ 1] = 10'd139;
        focus_phase_lut[ 9][ 2] = 10'd111;
        focus_phase_lut[ 9][ 3] = 10'd139;
        focus_phase_lut[ 9][ 4] = 10'd222;
        focus_phase_lut[ 9][ 5] = 10'd139;
        focus_phase_lut[ 9][ 6] = 10'd56;
        focus_phase_lut[ 9][ 7] = 10'd28;
        focus_phase_lut[ 9][ 8] = 10'd56;
        focus_phase_lut[ 9][ 9] = 10'd139;
        focus_phase_lut[ 9][10] = 10'd111;
        focus_phase_lut[ 9][11] = 10'd28;
        focus_phase_lut[ 9][12] = 10'd0;
        focus_phase_lut[ 9][13] = 10'd28;
        focus_phase_lut[ 9][14] = 10'd111;
        focus_phase_lut[ 9][15] = 10'd139;
        focus_phase_lut[ 9][16] = 10'd56;
        focus_phase_lut[ 9][17] = 10'd28;
        focus_phase_lut[ 9][18] = 10'd56;
        focus_phase_lut[ 9][19] = 10'd139;
        focus_phase_lut[ 9][20] = 10'd222;
        focus_phase_lut[ 9][21] = 10'd139;
        focus_phase_lut[ 9][22] = 10'd111;
        focus_phase_lut[ 9][23] = 10'd139;
        focus_phase_lut[ 9][24] = 10'd222;
        // focus_dist=600mm (idx=10)
        focus_phase_lut[10][ 0] = 10'd204;
        focus_phase_lut[10][ 1] = 10'd127;
        focus_phase_lut[10][ 2] = 10'd102;
        focus_phase_lut[10][ 3] = 10'd127;
        focus_phase_lut[10][ 4] = 10'd204;
        focus_phase_lut[10][ 5] = 10'd127;
        focus_phase_lut[10][ 6] = 10'd51;
        focus_phase_lut[10][ 7] = 10'd25;
        focus_phase_lut[10][ 8] = 10'd51;
        focus_phase_lut[10][ 9] = 10'd127;
        focus_phase_lut[10][10] = 10'd102;
        focus_phase_lut[10][11] = 10'd25;
        focus_phase_lut[10][12] = 10'd0;
        focus_phase_lut[10][13] = 10'd25;
        focus_phase_lut[10][14] = 10'd102;
        focus_phase_lut[10][15] = 10'd127;
        focus_phase_lut[10][16] = 10'd51;
        focus_phase_lut[10][17] = 10'd25;
        focus_phase_lut[10][18] = 10'd51;
        focus_phase_lut[10][19] = 10'd127;
        focus_phase_lut[10][20] = 10'd204;
        focus_phase_lut[10][21] = 10'd127;
        focus_phase_lut[10][22] = 10'd102;
        focus_phase_lut[10][23] = 10'd127;
        focus_phase_lut[10][24] = 10'd204;
        // focus_dist=650mm (idx=11)
        focus_phase_lut[11][ 0] = 10'd188;
        focus_phase_lut[11][ 1] = 10'd117;
        focus_phase_lut[11][ 2] = 10'd94;
        focus_phase_lut[11][ 3] = 10'd117;
        focus_phase_lut[11][ 4] = 10'd188;
        focus_phase_lut[11][ 5] = 10'd117;
        focus_phase_lut[11][ 6] = 10'd47;
        focus_phase_lut[11][ 7] = 10'd24;
        focus_phase_lut[11][ 8] = 10'd47;
        focus_phase_lut[11][ 9] = 10'd117;
        focus_phase_lut[11][10] = 10'd94;
        focus_phase_lut[11][11] = 10'd24;
        focus_phase_lut[11][12] = 10'd0;
        focus_phase_lut[11][13] = 10'd24;
        focus_phase_lut[11][14] = 10'd94;
        focus_phase_lut[11][15] = 10'd117;
        focus_phase_lut[11][16] = 10'd47;
        focus_phase_lut[11][17] = 10'd24;
        focus_phase_lut[11][18] = 10'd47;
        focus_phase_lut[11][19] = 10'd117;
        focus_phase_lut[11][20] = 10'd188;
        focus_phase_lut[11][21] = 10'd117;
        focus_phase_lut[11][22] = 10'd94;
        focus_phase_lut[11][23] = 10'd117;
        focus_phase_lut[11][24] = 10'd188;
        // focus_dist=700mm (idx=12)
        focus_phase_lut[12][ 0] = 10'd175;
        focus_phase_lut[12][ 1] = 10'd109;
        focus_phase_lut[12][ 2] = 10'd87;
        focus_phase_lut[12][ 3] = 10'd109;
        focus_phase_lut[12][ 4] = 10'd175;
        focus_phase_lut[12][ 5] = 10'd109;
        focus_phase_lut[12][ 6] = 10'd44;
        focus_phase_lut[12][ 7] = 10'd22;
        focus_phase_lut[12][ 8] = 10'd44;
        focus_phase_lut[12][ 9] = 10'd109;
        focus_phase_lut[12][10] = 10'd87;
        focus_phase_lut[12][11] = 10'd22;
        focus_phase_lut[12][12] = 10'd0;
        focus_phase_lut[12][13] = 10'd22;
        focus_phase_lut[12][14] = 10'd87;
        focus_phase_lut[12][15] = 10'd109;
        focus_phase_lut[12][16] = 10'd44;
        focus_phase_lut[12][17] = 10'd22;
        focus_phase_lut[12][18] = 10'd44;
        focus_phase_lut[12][19] = 10'd109;
        focus_phase_lut[12][20] = 10'd175;
        focus_phase_lut[12][21] = 10'd109;
        focus_phase_lut[12][22] = 10'd87;
        focus_phase_lut[12][23] = 10'd109;
        focus_phase_lut[12][24] = 10'd175;
        // focus_dist=750mm (idx=13)
        focus_phase_lut[13][ 0] = 10'd163;
        focus_phase_lut[13][ 1] = 10'd102;
        focus_phase_lut[13][ 2] = 10'd81;
        focus_phase_lut[13][ 3] = 10'd102;
        focus_phase_lut[13][ 4] = 10'd163;
        focus_phase_lut[13][ 5] = 10'd102;
        focus_phase_lut[13][ 6] = 10'd41;
        focus_phase_lut[13][ 7] = 10'd20;
        focus_phase_lut[13][ 8] = 10'd41;
        focus_phase_lut[13][ 9] = 10'd102;
        focus_phase_lut[13][10] = 10'd81;
        focus_phase_lut[13][11] = 10'd20;
        focus_phase_lut[13][12] = 10'd0;
        focus_phase_lut[13][13] = 10'd20;
        focus_phase_lut[13][14] = 10'd81;
        focus_phase_lut[13][15] = 10'd102;
        focus_phase_lut[13][16] = 10'd41;
        focus_phase_lut[13][17] = 10'd20;
        focus_phase_lut[13][18] = 10'd41;
        focus_phase_lut[13][19] = 10'd102;
        focus_phase_lut[13][20] = 10'd163;
        focus_phase_lut[13][21] = 10'd102;
        focus_phase_lut[13][22] = 10'd81;
        focus_phase_lut[13][23] = 10'd102;
        focus_phase_lut[13][24] = 10'd163;
        // focus_dist=800mm (idx=14)
        focus_phase_lut[14][ 0] = 10'd153;
        focus_phase_lut[14][ 1] = 10'd95;
        focus_phase_lut[14][ 2] = 10'd76;
        focus_phase_lut[14][ 3] = 10'd95;
        focus_phase_lut[14][ 4] = 10'd153;
        focus_phase_lut[14][ 5] = 10'd95;
        focus_phase_lut[14][ 6] = 10'd38;
        focus_phase_lut[14][ 7] = 10'd19;
        focus_phase_lut[14][ 8] = 10'd38;
        focus_phase_lut[14][ 9] = 10'd95;
        focus_phase_lut[14][10] = 10'd76;
        focus_phase_lut[14][11] = 10'd19;
        focus_phase_lut[14][12] = 10'd0;
        focus_phase_lut[14][13] = 10'd19;
        focus_phase_lut[14][14] = 10'd76;
        focus_phase_lut[14][15] = 10'd95;
        focus_phase_lut[14][16] = 10'd38;
        focus_phase_lut[14][17] = 10'd19;
        focus_phase_lut[14][18] = 10'd38;
        focus_phase_lut[14][19] = 10'd95;
        focus_phase_lut[14][20] = 10'd153;
        focus_phase_lut[14][21] = 10'd95;
        focus_phase_lut[14][22] = 10'd76;
        focus_phase_lut[14][23] = 10'd95;
        focus_phase_lut[14][24] = 10'd153;
        // focus_dist=850mm (idx=15)
        focus_phase_lut[15][ 0] = 10'd144;
        focus_phase_lut[15][ 1] = 10'd90;
        focus_phase_lut[15][ 2] = 10'd72;
        focus_phase_lut[15][ 3] = 10'd90;
        focus_phase_lut[15][ 4] = 10'd144;
        focus_phase_lut[15][ 5] = 10'd90;
        focus_phase_lut[15][ 6] = 10'd36;
        focus_phase_lut[15][ 7] = 10'd18;
        focus_phase_lut[15][ 8] = 10'd36;
        focus_phase_lut[15][ 9] = 10'd90;
        focus_phase_lut[15][10] = 10'd72;
        focus_phase_lut[15][11] = 10'd18;
        focus_phase_lut[15][12] = 10'd0;
        focus_phase_lut[15][13] = 10'd18;
        focus_phase_lut[15][14] = 10'd72;
        focus_phase_lut[15][15] = 10'd90;
        focus_phase_lut[15][16] = 10'd36;
        focus_phase_lut[15][17] = 10'd18;
        focus_phase_lut[15][18] = 10'd36;
        focus_phase_lut[15][19] = 10'd90;
        focus_phase_lut[15][20] = 10'd144;
        focus_phase_lut[15][21] = 10'd90;
        focus_phase_lut[15][22] = 10'd72;
        focus_phase_lut[15][23] = 10'd90;
        focus_phase_lut[15][24] = 10'd144;
        // focus_dist=900mm (idx=16)
        focus_phase_lut[16][ 0] = 10'd136;
        focus_phase_lut[16][ 1] = 10'd85;
        focus_phase_lut[16][ 2] = 10'd68;
        focus_phase_lut[16][ 3] = 10'd85;
        focus_phase_lut[16][ 4] = 10'd136;
        focus_phase_lut[16][ 5] = 10'd85;
        focus_phase_lut[16][ 6] = 10'd34;
        focus_phase_lut[16][ 7] = 10'd17;
        focus_phase_lut[16][ 8] = 10'd34;
        focus_phase_lut[16][ 9] = 10'd85;
        focus_phase_lut[16][10] = 10'd68;
        focus_phase_lut[16][11] = 10'd17;
        focus_phase_lut[16][12] = 10'd0;
        focus_phase_lut[16][13] = 10'd17;
        focus_phase_lut[16][14] = 10'd68;
        focus_phase_lut[16][15] = 10'd85;
        focus_phase_lut[16][16] = 10'd34;
        focus_phase_lut[16][17] = 10'd17;
        focus_phase_lut[16][18] = 10'd34;
        focus_phase_lut[16][19] = 10'd85;
        focus_phase_lut[16][20] = 10'd136;
        focus_phase_lut[16][21] = 10'd85;
        focus_phase_lut[16][22] = 10'd68;
        focus_phase_lut[16][23] = 10'd85;
        focus_phase_lut[16][24] = 10'd136;
        // focus_dist=950mm (idx=17)
        focus_phase_lut[17][ 0] = 10'd129;
        focus_phase_lut[17][ 1] = 10'd80;
        focus_phase_lut[17][ 2] = 10'd64;
        focus_phase_lut[17][ 3] = 10'd80;
        focus_phase_lut[17][ 4] = 10'd129;
        focus_phase_lut[17][ 5] = 10'd80;
        focus_phase_lut[17][ 6] = 10'd32;
        focus_phase_lut[17][ 7] = 10'd16;
        focus_phase_lut[17][ 8] = 10'd32;
        focus_phase_lut[17][ 9] = 10'd80;
        focus_phase_lut[17][10] = 10'd64;
        focus_phase_lut[17][11] = 10'd16;
        focus_phase_lut[17][12] = 10'd0;
        focus_phase_lut[17][13] = 10'd16;
        focus_phase_lut[17][14] = 10'd64;
        focus_phase_lut[17][15] = 10'd80;
        focus_phase_lut[17][16] = 10'd32;
        focus_phase_lut[17][17] = 10'd16;
        focus_phase_lut[17][18] = 10'd32;
        focus_phase_lut[17][19] = 10'd80;
        focus_phase_lut[17][20] = 10'd129;
        focus_phase_lut[17][21] = 10'd80;
        focus_phase_lut[17][22] = 10'd64;
        focus_phase_lut[17][23] = 10'd80;
        focus_phase_lut[17][24] = 10'd129;
        // focus_dist=1000mm (idx=18)
        focus_phase_lut[18][ 0] = 10'd122;
        focus_phase_lut[18][ 1] = 10'd76;
        focus_phase_lut[18][ 2] = 10'd61;
        focus_phase_lut[18][ 3] = 10'd76;
        focus_phase_lut[18][ 4] = 10'd122;
        focus_phase_lut[18][ 5] = 10'd76;
        focus_phase_lut[18][ 6] = 10'd31;
        focus_phase_lut[18][ 7] = 10'd15;
        focus_phase_lut[18][ 8] = 10'd31;
        focus_phase_lut[18][ 9] = 10'd76;
        focus_phase_lut[18][10] = 10'd61;
        focus_phase_lut[18][11] = 10'd15;
        focus_phase_lut[18][12] = 10'd0;
        focus_phase_lut[18][13] = 10'd15;
        focus_phase_lut[18][14] = 10'd61;
        focus_phase_lut[18][15] = 10'd76;
        focus_phase_lut[18][16] = 10'd31;
        focus_phase_lut[18][17] = 10'd15;
        focus_phase_lut[18][18] = 10'd31;
        focus_phase_lut[18][19] = 10'd76;
        focus_phase_lut[18][20] = 10'd122;
        focus_phase_lut[18][21] = 10'd76;
        focus_phase_lut[18][22] = 10'd61;
        focus_phase_lut[18][23] = 10'd76;
        focus_phase_lut[18][24] = 10'd122;
        // focus_dist=1050mm (idx=19)
        focus_phase_lut[19][ 0] = 10'd116;
        focus_phase_lut[19][ 1] = 10'd73;
        focus_phase_lut[19][ 2] = 10'd58;
        focus_phase_lut[19][ 3] = 10'd73;
        focus_phase_lut[19][ 4] = 10'd116;
        focus_phase_lut[19][ 5] = 10'd73;
        focus_phase_lut[19][ 6] = 10'd29;
        focus_phase_lut[19][ 7] = 10'd15;
        focus_phase_lut[19][ 8] = 10'd29;
        focus_phase_lut[19][ 9] = 10'd73;
        focus_phase_lut[19][10] = 10'd58;
        focus_phase_lut[19][11] = 10'd15;
        focus_phase_lut[19][12] = 10'd0;
        focus_phase_lut[19][13] = 10'd15;
        focus_phase_lut[19][14] = 10'd58;
        focus_phase_lut[19][15] = 10'd73;
        focus_phase_lut[19][16] = 10'd29;
        focus_phase_lut[19][17] = 10'd15;
        focus_phase_lut[19][18] = 10'd29;
        focus_phase_lut[19][19] = 10'd73;
        focus_phase_lut[19][20] = 10'd116;
        focus_phase_lut[19][21] = 10'd73;
        focus_phase_lut[19][22] = 10'd58;
        focus_phase_lut[19][23] = 10'd73;
        focus_phase_lut[19][24] = 10'd116;
        // focus_dist=1100mm (idx=20)
        focus_phase_lut[20][ 0] = 10'd111;
        focus_phase_lut[20][ 1] = 10'd69;
        focus_phase_lut[20][ 2] = 10'd56;
        focus_phase_lut[20][ 3] = 10'd69;
        focus_phase_lut[20][ 4] = 10'd111;
        focus_phase_lut[20][ 5] = 10'd69;
        focus_phase_lut[20][ 6] = 10'd28;
        focus_phase_lut[20][ 7] = 10'd14;
        focus_phase_lut[20][ 8] = 10'd28;
        focus_phase_lut[20][ 9] = 10'd69;
        focus_phase_lut[20][10] = 10'd56;
        focus_phase_lut[20][11] = 10'd14;
        focus_phase_lut[20][12] = 10'd0;
        focus_phase_lut[20][13] = 10'd14;
        focus_phase_lut[20][14] = 10'd56;
        focus_phase_lut[20][15] = 10'd69;
        focus_phase_lut[20][16] = 10'd28;
        focus_phase_lut[20][17] = 10'd14;
        focus_phase_lut[20][18] = 10'd28;
        focus_phase_lut[20][19] = 10'd69;
        focus_phase_lut[20][20] = 10'd111;
        focus_phase_lut[20][21] = 10'd69;
        focus_phase_lut[20][22] = 10'd56;
        focus_phase_lut[20][23] = 10'd69;
        focus_phase_lut[20][24] = 10'd111;
        // focus_dist=1150mm (idx=21)
        focus_phase_lut[21][ 0] = 10'd106;
        focus_phase_lut[21][ 1] = 10'd66;
        focus_phase_lut[21][ 2] = 10'd53;
        focus_phase_lut[21][ 3] = 10'd66;
        focus_phase_lut[21][ 4] = 10'd106;
        focus_phase_lut[21][ 5] = 10'd66;
        focus_phase_lut[21][ 6] = 10'd27;
        focus_phase_lut[21][ 7] = 10'd13;
        focus_phase_lut[21][ 8] = 10'd27;
        focus_phase_lut[21][ 9] = 10'd66;
        focus_phase_lut[21][10] = 10'd53;
        focus_phase_lut[21][11] = 10'd13;
        focus_phase_lut[21][12] = 10'd0;
        focus_phase_lut[21][13] = 10'd13;
        focus_phase_lut[21][14] = 10'd53;
        focus_phase_lut[21][15] = 10'd66;
        focus_phase_lut[21][16] = 10'd27;
        focus_phase_lut[21][17] = 10'd13;
        focus_phase_lut[21][18] = 10'd27;
        focus_phase_lut[21][19] = 10'd66;
        focus_phase_lut[21][20] = 10'd106;
        focus_phase_lut[21][21] = 10'd66;
        focus_phase_lut[21][22] = 10'd53;
        focus_phase_lut[21][23] = 10'd66;
        focus_phase_lut[21][24] = 10'd106;
        // focus_dist=1200mm (idx=22)
        focus_phase_lut[22][ 0] = 10'd102;
        focus_phase_lut[22][ 1] = 10'd64;
        focus_phase_lut[22][ 2] = 10'd51;
        focus_phase_lut[22][ 3] = 10'd64;
        focus_phase_lut[22][ 4] = 10'd102;
        focus_phase_lut[22][ 5] = 10'd64;
        focus_phase_lut[22][ 6] = 10'd25;
        focus_phase_lut[22][ 7] = 10'd13;
        focus_phase_lut[22][ 8] = 10'd25;
        focus_phase_lut[22][ 9] = 10'd64;
        focus_phase_lut[22][10] = 10'd51;
        focus_phase_lut[22][11] = 10'd13;
        focus_phase_lut[22][12] = 10'd0;
        focus_phase_lut[22][13] = 10'd13;
        focus_phase_lut[22][14] = 10'd51;
        focus_phase_lut[22][15] = 10'd64;
        focus_phase_lut[22][16] = 10'd25;
        focus_phase_lut[22][17] = 10'd13;
        focus_phase_lut[22][18] = 10'd25;
        focus_phase_lut[22][19] = 10'd64;
        focus_phase_lut[22][20] = 10'd102;
        focus_phase_lut[22][21] = 10'd64;
        focus_phase_lut[22][22] = 10'd51;
        focus_phase_lut[22][23] = 10'd64;
        focus_phase_lut[22][24] = 10'd102;
        // focus_dist=1250mm (idx=23)
        focus_phase_lut[23][ 0] = 10'd98;
        focus_phase_lut[23][ 1] = 10'd61;
        focus_phase_lut[23][ 2] = 10'd49;
        focus_phase_lut[23][ 3] = 10'd61;
        focus_phase_lut[23][ 4] = 10'd98;
        focus_phase_lut[23][ 5] = 10'd61;
        focus_phase_lut[23][ 6] = 10'd24;
        focus_phase_lut[23][ 7] = 10'd12;
        focus_phase_lut[23][ 8] = 10'd24;
        focus_phase_lut[23][ 9] = 10'd61;
        focus_phase_lut[23][10] = 10'd49;
        focus_phase_lut[23][11] = 10'd12;
        focus_phase_lut[23][12] = 10'd0;
        focus_phase_lut[23][13] = 10'd12;
        focus_phase_lut[23][14] = 10'd49;
        focus_phase_lut[23][15] = 10'd61;
        focus_phase_lut[23][16] = 10'd24;
        focus_phase_lut[23][17] = 10'd12;
        focus_phase_lut[23][18] = 10'd24;
        focus_phase_lut[23][19] = 10'd61;
        focus_phase_lut[23][20] = 10'd98;
        focus_phase_lut[23][21] = 10'd61;
        focus_phase_lut[23][22] = 10'd49;
        focus_phase_lut[23][23] = 10'd61;
        focus_phase_lut[23][24] = 10'd98;
        // focus_dist=1300mm (idx=24)
        focus_phase_lut[24][ 0] = 10'd94;
        focus_phase_lut[24][ 1] = 10'd59;
        focus_phase_lut[24][ 2] = 10'd47;
        focus_phase_lut[24][ 3] = 10'd59;
        focus_phase_lut[24][ 4] = 10'd94;
        focus_phase_lut[24][ 5] = 10'd59;
        focus_phase_lut[24][ 6] = 10'd24;
        focus_phase_lut[24][ 7] = 10'd12;
        focus_phase_lut[24][ 8] = 10'd24;
        focus_phase_lut[24][ 9] = 10'd59;
        focus_phase_lut[24][10] = 10'd47;
        focus_phase_lut[24][11] = 10'd12;
        focus_phase_lut[24][12] = 10'd0;
        focus_phase_lut[24][13] = 10'd12;
        focus_phase_lut[24][14] = 10'd47;
        focus_phase_lut[24][15] = 10'd59;
        focus_phase_lut[24][16] = 10'd24;
        focus_phase_lut[24][17] = 10'd12;
        focus_phase_lut[24][18] = 10'd24;
        focus_phase_lut[24][19] = 10'd59;
        focus_phase_lut[24][20] = 10'd94;
        focus_phase_lut[24][21] = 10'd59;
        focus_phase_lut[24][22] = 10'd47;
        focus_phase_lut[24][23] = 10'd59;
        focus_phase_lut[24][24] = 10'd94;
        // focus_dist=1350mm (idx=25)
        focus_phase_lut[25][ 0] = 10'd91;
        focus_phase_lut[25][ 1] = 10'd57;
        focus_phase_lut[25][ 2] = 10'd45;
        focus_phase_lut[25][ 3] = 10'd57;
        focus_phase_lut[25][ 4] = 10'd91;
        focus_phase_lut[25][ 5] = 10'd57;
        focus_phase_lut[25][ 6] = 10'd23;
        focus_phase_lut[25][ 7] = 10'd11;
        focus_phase_lut[25][ 8] = 10'd23;
        focus_phase_lut[25][ 9] = 10'd57;
        focus_phase_lut[25][10] = 10'd45;
        focus_phase_lut[25][11] = 10'd11;
        focus_phase_lut[25][12] = 10'd0;
        focus_phase_lut[25][13] = 10'd11;
        focus_phase_lut[25][14] = 10'd45;
        focus_phase_lut[25][15] = 10'd57;
        focus_phase_lut[25][16] = 10'd23;
        focus_phase_lut[25][17] = 10'd11;
        focus_phase_lut[25][18] = 10'd23;
        focus_phase_lut[25][19] = 10'd57;
        focus_phase_lut[25][20] = 10'd91;
        focus_phase_lut[25][21] = 10'd57;
        focus_phase_lut[25][22] = 10'd45;
        focus_phase_lut[25][23] = 10'd57;
        focus_phase_lut[25][24] = 10'd91;
        // focus_dist=1400mm (idx=26)
        focus_phase_lut[26][ 0] = 10'd87;
        focus_phase_lut[26][ 1] = 10'd55;
        focus_phase_lut[26][ 2] = 10'd44;
        focus_phase_lut[26][ 3] = 10'd55;
        focus_phase_lut[26][ 4] = 10'd87;
        focus_phase_lut[26][ 5] = 10'd55;
        focus_phase_lut[26][ 6] = 10'd22;
        focus_phase_lut[26][ 7] = 10'd11;
        focus_phase_lut[26][ 8] = 10'd22;
        focus_phase_lut[26][ 9] = 10'd55;
        focus_phase_lut[26][10] = 10'd44;
        focus_phase_lut[26][11] = 10'd11;
        focus_phase_lut[26][12] = 10'd0;
        focus_phase_lut[26][13] = 10'd11;
        focus_phase_lut[26][14] = 10'd44;
        focus_phase_lut[26][15] = 10'd55;
        focus_phase_lut[26][16] = 10'd22;
        focus_phase_lut[26][17] = 10'd11;
        focus_phase_lut[26][18] = 10'd22;
        focus_phase_lut[26][19] = 10'd55;
        focus_phase_lut[26][20] = 10'd87;
        focus_phase_lut[26][21] = 10'd55;
        focus_phase_lut[26][22] = 10'd44;
        focus_phase_lut[26][23] = 10'd55;
        focus_phase_lut[26][24] = 10'd87;
        // focus_dist=1450mm (idx=27)
        focus_phase_lut[27][ 0] = 10'd84;
        focus_phase_lut[27][ 1] = 10'd53;
        focus_phase_lut[27][ 2] = 10'd42;
        focus_phase_lut[27][ 3] = 10'd53;
        focus_phase_lut[27][ 4] = 10'd84;
        focus_phase_lut[27][ 5] = 10'd53;
        focus_phase_lut[27][ 6] = 10'd21;
        focus_phase_lut[27][ 7] = 10'd11;
        focus_phase_lut[27][ 8] = 10'd21;
        focus_phase_lut[27][ 9] = 10'd53;
        focus_phase_lut[27][10] = 10'd42;
        focus_phase_lut[27][11] = 10'd11;
        focus_phase_lut[27][12] = 10'd0;
        focus_phase_lut[27][13] = 10'd11;
        focus_phase_lut[27][14] = 10'd42;
        focus_phase_lut[27][15] = 10'd53;
        focus_phase_lut[27][16] = 10'd21;
        focus_phase_lut[27][17] = 10'd11;
        focus_phase_lut[27][18] = 10'd21;
        focus_phase_lut[27][19] = 10'd53;
        focus_phase_lut[27][20] = 10'd84;
        focus_phase_lut[27][21] = 10'd53;
        focus_phase_lut[27][22] = 10'd42;
        focus_phase_lut[27][23] = 10'd53;
        focus_phase_lut[27][24] = 10'd84;
        // focus_dist=1500mm (idx=28)
        focus_phase_lut[28][ 0] = 10'd82;
        focus_phase_lut[28][ 1] = 10'd51;
        focus_phase_lut[28][ 2] = 10'd41;
        focus_phase_lut[28][ 3] = 10'd51;
        focus_phase_lut[28][ 4] = 10'd82;
        focus_phase_lut[28][ 5] = 10'd51;
        focus_phase_lut[28][ 6] = 10'd20;
        focus_phase_lut[28][ 7] = 10'd10;
        focus_phase_lut[28][ 8] = 10'd20;
        focus_phase_lut[28][ 9] = 10'd51;
        focus_phase_lut[28][10] = 10'd41;
        focus_phase_lut[28][11] = 10'd10;
        focus_phase_lut[28][12] = 10'd0;
        focus_phase_lut[28][13] = 10'd10;
        focus_phase_lut[28][14] = 10'd41;
        focus_phase_lut[28][15] = 10'd51;
        focus_phase_lut[28][16] = 10'd20;
        focus_phase_lut[28][17] = 10'd10;
        focus_phase_lut[28][18] = 10'd20;
        focus_phase_lut[28][19] = 10'd51;
        focus_phase_lut[28][20] = 10'd82;
        focus_phase_lut[28][21] = 10'd51;
        focus_phase_lut[28][22] = 10'd41;
        focus_phase_lut[28][23] = 10'd51;
        focus_phase_lut[28][24] = 10'd82;
        // focus_dist=1550mm (idx=29)
        focus_phase_lut[29][ 0] = 10'd79;
        focus_phase_lut[29][ 1] = 10'd49;
        focus_phase_lut[29][ 2] = 10'd39;
        focus_phase_lut[29][ 3] = 10'd49;
        focus_phase_lut[29][ 4] = 10'd79;
        focus_phase_lut[29][ 5] = 10'd49;
        focus_phase_lut[29][ 6] = 10'd20;
        focus_phase_lut[29][ 7] = 10'd10;
        focus_phase_lut[29][ 8] = 10'd20;
        focus_phase_lut[29][ 9] = 10'd49;
        focus_phase_lut[29][10] = 10'd39;
        focus_phase_lut[29][11] = 10'd10;
        focus_phase_lut[29][12] = 10'd0;
        focus_phase_lut[29][13] = 10'd10;
        focus_phase_lut[29][14] = 10'd39;
        focus_phase_lut[29][15] = 10'd49;
        focus_phase_lut[29][16] = 10'd20;
        focus_phase_lut[29][17] = 10'd10;
        focus_phase_lut[29][18] = 10'd20;
        focus_phase_lut[29][19] = 10'd49;
        focus_phase_lut[29][20] = 10'd79;
        focus_phase_lut[29][21] = 10'd49;
        focus_phase_lut[29][22] = 10'd39;
        focus_phase_lut[29][23] = 10'd49;
        focus_phase_lut[29][24] = 10'd79;
        // focus_dist=1600mm (idx=30)
        focus_phase_lut[30][ 0] = 10'd76;
        focus_phase_lut[30][ 1] = 10'd48;
        focus_phase_lut[30][ 2] = 10'd38;
        focus_phase_lut[30][ 3] = 10'd48;
        focus_phase_lut[30][ 4] = 10'd76;
        focus_phase_lut[30][ 5] = 10'd48;
        focus_phase_lut[30][ 6] = 10'd19;
        focus_phase_lut[30][ 7] = 10'd10;
        focus_phase_lut[30][ 8] = 10'd19;
        focus_phase_lut[30][ 9] = 10'd48;
        focus_phase_lut[30][10] = 10'd38;
        focus_phase_lut[30][11] = 10'd10;
        focus_phase_lut[30][12] = 10'd0;
        focus_phase_lut[30][13] = 10'd10;
        focus_phase_lut[30][14] = 10'd38;
        focus_phase_lut[30][15] = 10'd48;
        focus_phase_lut[30][16] = 10'd19;
        focus_phase_lut[30][17] = 10'd10;
        focus_phase_lut[30][18] = 10'd19;
        focus_phase_lut[30][19] = 10'd48;
        focus_phase_lut[30][20] = 10'd76;
        focus_phase_lut[30][21] = 10'd48;
        focus_phase_lut[30][22] = 10'd38;
        focus_phase_lut[30][23] = 10'd48;
        focus_phase_lut[30][24] = 10'd76;
        // focus_dist=1650mm (idx=31)
        focus_phase_lut[31][ 0] = 10'd74;
        focus_phase_lut[31][ 1] = 10'd46;
        focus_phase_lut[31][ 2] = 10'd37;
        focus_phase_lut[31][ 3] = 10'd46;
        focus_phase_lut[31][ 4] = 10'd74;
        focus_phase_lut[31][ 5] = 10'd46;
        focus_phase_lut[31][ 6] = 10'd19;
        focus_phase_lut[31][ 7] = 10'd9;
        focus_phase_lut[31][ 8] = 10'd19;
        focus_phase_lut[31][ 9] = 10'd46;
        focus_phase_lut[31][10] = 10'd37;
        focus_phase_lut[31][11] = 10'd9;
        focus_phase_lut[31][12] = 10'd0;
        focus_phase_lut[31][13] = 10'd9;
        focus_phase_lut[31][14] = 10'd37;
        focus_phase_lut[31][15] = 10'd46;
        focus_phase_lut[31][16] = 10'd19;
        focus_phase_lut[31][17] = 10'd9;
        focus_phase_lut[31][18] = 10'd19;
        focus_phase_lut[31][19] = 10'd46;
        focus_phase_lut[31][20] = 10'd74;
        focus_phase_lut[31][21] = 10'd46;
        focus_phase_lut[31][22] = 10'd37;
        focus_phase_lut[31][23] = 10'd46;
        focus_phase_lut[31][24] = 10'd74;
        // focus_dist=1700mm (idx=32)
        focus_phase_lut[32][ 0] = 10'd72;
        focus_phase_lut[32][ 1] = 10'd45;
        focus_phase_lut[32][ 2] = 10'd36;
        focus_phase_lut[32][ 3] = 10'd45;
        focus_phase_lut[32][ 4] = 10'd72;
        focus_phase_lut[32][ 5] = 10'd45;
        focus_phase_lut[32][ 6] = 10'd18;
        focus_phase_lut[32][ 7] = 10'd9;
        focus_phase_lut[32][ 8] = 10'd18;
        focus_phase_lut[32][ 9] = 10'd45;
        focus_phase_lut[32][10] = 10'd36;
        focus_phase_lut[32][11] = 10'd9;
        focus_phase_lut[32][12] = 10'd0;
        focus_phase_lut[32][13] = 10'd9;
        focus_phase_lut[32][14] = 10'd36;
        focus_phase_lut[32][15] = 10'd45;
        focus_phase_lut[32][16] = 10'd18;
        focus_phase_lut[32][17] = 10'd9;
        focus_phase_lut[32][18] = 10'd18;
        focus_phase_lut[32][19] = 10'd45;
        focus_phase_lut[32][20] = 10'd72;
        focus_phase_lut[32][21] = 10'd45;
        focus_phase_lut[32][22] = 10'd36;
        focus_phase_lut[32][23] = 10'd45;
        focus_phase_lut[32][24] = 10'd72;
        // focus_dist=1750mm (idx=33)
        focus_phase_lut[33][ 0] = 10'd70;
        focus_phase_lut[33][ 1] = 10'd44;
        focus_phase_lut[33][ 2] = 10'd35;
        focus_phase_lut[33][ 3] = 10'd44;
        focus_phase_lut[33][ 4] = 10'd70;
        focus_phase_lut[33][ 5] = 10'd44;
        focus_phase_lut[33][ 6] = 10'd17;
        focus_phase_lut[33][ 7] = 10'd9;
        focus_phase_lut[33][ 8] = 10'd17;
        focus_phase_lut[33][ 9] = 10'd44;
        focus_phase_lut[33][10] = 10'd35;
        focus_phase_lut[33][11] = 10'd9;
        focus_phase_lut[33][12] = 10'd0;
        focus_phase_lut[33][13] = 10'd9;
        focus_phase_lut[33][14] = 10'd35;
        focus_phase_lut[33][15] = 10'd44;
        focus_phase_lut[33][16] = 10'd17;
        focus_phase_lut[33][17] = 10'd9;
        focus_phase_lut[33][18] = 10'd17;
        focus_phase_lut[33][19] = 10'd44;
        focus_phase_lut[33][20] = 10'd70;
        focus_phase_lut[33][21] = 10'd44;
        focus_phase_lut[33][22] = 10'd35;
        focus_phase_lut[33][23] = 10'd44;
        focus_phase_lut[33][24] = 10'd70;
        // focus_dist=1800mm (idx=34)
        focus_phase_lut[34][ 0] = 10'd68;
        focus_phase_lut[34][ 1] = 10'd42;
        focus_phase_lut[34][ 2] = 10'd34;
        focus_phase_lut[34][ 3] = 10'd42;
        focus_phase_lut[34][ 4] = 10'd68;
        focus_phase_lut[34][ 5] = 10'd42;
        focus_phase_lut[34][ 6] = 10'd17;
        focus_phase_lut[34][ 7] = 10'd8;
        focus_phase_lut[34][ 8] = 10'd17;
        focus_phase_lut[34][ 9] = 10'd42;
        focus_phase_lut[34][10] = 10'd34;
        focus_phase_lut[34][11] = 10'd8;
        focus_phase_lut[34][12] = 10'd0;
        focus_phase_lut[34][13] = 10'd8;
        focus_phase_lut[34][14] = 10'd34;
        focus_phase_lut[34][15] = 10'd42;
        focus_phase_lut[34][16] = 10'd17;
        focus_phase_lut[34][17] = 10'd8;
        focus_phase_lut[34][18] = 10'd17;
        focus_phase_lut[34][19] = 10'd42;
        focus_phase_lut[34][20] = 10'd68;
        focus_phase_lut[34][21] = 10'd42;
        focus_phase_lut[34][22] = 10'd34;
        focus_phase_lut[34][23] = 10'd42;
        focus_phase_lut[34][24] = 10'd68;
        // focus_dist=1850mm (idx=35)
        focus_phase_lut[35][ 0] = 10'd66;
        focus_phase_lut[35][ 1] = 10'd41;
        focus_phase_lut[35][ 2] = 10'd33;
        focus_phase_lut[35][ 3] = 10'd41;
        focus_phase_lut[35][ 4] = 10'd66;
        focus_phase_lut[35][ 5] = 10'd41;
        focus_phase_lut[35][ 6] = 10'd17;
        focus_phase_lut[35][ 7] = 10'd8;
        focus_phase_lut[35][ 8] = 10'd17;
        focus_phase_lut[35][ 9] = 10'd41;
        focus_phase_lut[35][10] = 10'd33;
        focus_phase_lut[35][11] = 10'd8;
        focus_phase_lut[35][12] = 10'd0;
        focus_phase_lut[35][13] = 10'd8;
        focus_phase_lut[35][14] = 10'd33;
        focus_phase_lut[35][15] = 10'd41;
        focus_phase_lut[35][16] = 10'd17;
        focus_phase_lut[35][17] = 10'd8;
        focus_phase_lut[35][18] = 10'd17;
        focus_phase_lut[35][19] = 10'd41;
        focus_phase_lut[35][20] = 10'd66;
        focus_phase_lut[35][21] = 10'd41;
        focus_phase_lut[35][22] = 10'd33;
        focus_phase_lut[35][23] = 10'd41;
        focus_phase_lut[35][24] = 10'd66;
        // focus_dist=1900mm (idx=36)
        focus_phase_lut[36][ 0] = 10'd64;
        focus_phase_lut[36][ 1] = 10'd40;
        focus_phase_lut[36][ 2] = 10'd32;
        focus_phase_lut[36][ 3] = 10'd40;
        focus_phase_lut[36][ 4] = 10'd64;
        focus_phase_lut[36][ 5] = 10'd40;
        focus_phase_lut[36][ 6] = 10'd16;
        focus_phase_lut[36][ 7] = 10'd8;
        focus_phase_lut[36][ 8] = 10'd16;
        focus_phase_lut[36][ 9] = 10'd40;
        focus_phase_lut[36][10] = 10'd32;
        focus_phase_lut[36][11] = 10'd8;
        focus_phase_lut[36][12] = 10'd0;
        focus_phase_lut[36][13] = 10'd8;
        focus_phase_lut[36][14] = 10'd32;
        focus_phase_lut[36][15] = 10'd40;
        focus_phase_lut[36][16] = 10'd16;
        focus_phase_lut[36][17] = 10'd8;
        focus_phase_lut[36][18] = 10'd16;
        focus_phase_lut[36][19] = 10'd40;
        focus_phase_lut[36][20] = 10'd64;
        focus_phase_lut[36][21] = 10'd40;
        focus_phase_lut[36][22] = 10'd32;
        focus_phase_lut[36][23] = 10'd40;
        focus_phase_lut[36][24] = 10'd64;
        // focus_dist=1950mm (idx=37)
        focus_phase_lut[37][ 0] = 10'd63;
        focus_phase_lut[37][ 1] = 10'd39;
        focus_phase_lut[37][ 2] = 10'd31;
        focus_phase_lut[37][ 3] = 10'd39;
        focus_phase_lut[37][ 4] = 10'd63;
        focus_phase_lut[37][ 5] = 10'd39;
        focus_phase_lut[37][ 6] = 10'd16;
        focus_phase_lut[37][ 7] = 10'd8;
        focus_phase_lut[37][ 8] = 10'd16;
        focus_phase_lut[37][ 9] = 10'd39;
        focus_phase_lut[37][10] = 10'd31;
        focus_phase_lut[37][11] = 10'd8;
        focus_phase_lut[37][12] = 10'd0;
        focus_phase_lut[37][13] = 10'd8;
        focus_phase_lut[37][14] = 10'd31;
        focus_phase_lut[37][15] = 10'd39;
        focus_phase_lut[37][16] = 10'd16;
        focus_phase_lut[37][17] = 10'd8;
        focus_phase_lut[37][18] = 10'd16;
        focus_phase_lut[37][19] = 10'd39;
        focus_phase_lut[37][20] = 10'd63;
        focus_phase_lut[37][21] = 10'd39;
        focus_phase_lut[37][22] = 10'd31;
        focus_phase_lut[37][23] = 10'd39;
        focus_phase_lut[37][24] = 10'd63;
        // focus_dist=2000mm (idx=38)
        focus_phase_lut[38][ 0] = 10'd61;
        focus_phase_lut[38][ 1] = 10'd38;
        focus_phase_lut[38][ 2] = 10'd31;
        focus_phase_lut[38][ 3] = 10'd38;
        focus_phase_lut[38][ 4] = 10'd61;
        focus_phase_lut[38][ 5] = 10'd38;
        focus_phase_lut[38][ 6] = 10'd15;
        focus_phase_lut[38][ 7] = 10'd8;
        focus_phase_lut[38][ 8] = 10'd15;
        focus_phase_lut[38][ 9] = 10'd38;
        focus_phase_lut[38][10] = 10'd31;
        focus_phase_lut[38][11] = 10'd8;
        focus_phase_lut[38][12] = 10'd0;
        focus_phase_lut[38][13] = 10'd8;
        focus_phase_lut[38][14] = 10'd31;
        focus_phase_lut[38][15] = 10'd38;
        focus_phase_lut[38][16] = 10'd15;
        focus_phase_lut[38][17] = 10'd8;
        focus_phase_lut[38][18] = 10'd15;
        focus_phase_lut[38][19] = 10'd38;
        focus_phase_lut[38][20] = 10'd61;
        focus_phase_lut[38][21] = 10'd38;
        focus_phase_lut[38][22] = 10'd31;
        focus_phase_lut[38][23] = 10'd38;
        focus_phase_lut[38][24] = 10'd61;
        // focus_dist=2050mm (idx=39)
        focus_phase_lut[39][ 0] = 10'd60;
        focus_phase_lut[39][ 1] = 10'd37;
        focus_phase_lut[39][ 2] = 10'd30;
        focus_phase_lut[39][ 3] = 10'd37;
        focus_phase_lut[39][ 4] = 10'd60;
        focus_phase_lut[39][ 5] = 10'd37;
        focus_phase_lut[39][ 6] = 10'd15;
        focus_phase_lut[39][ 7] = 10'd7;
        focus_phase_lut[39][ 8] = 10'd15;
        focus_phase_lut[39][ 9] = 10'd37;
        focus_phase_lut[39][10] = 10'd30;
        focus_phase_lut[39][11] = 10'd7;
        focus_phase_lut[39][12] = 10'd0;
        focus_phase_lut[39][13] = 10'd7;
        focus_phase_lut[39][14] = 10'd30;
        focus_phase_lut[39][15] = 10'd37;
        focus_phase_lut[39][16] = 10'd15;
        focus_phase_lut[39][17] = 10'd7;
        focus_phase_lut[39][18] = 10'd15;
        focus_phase_lut[39][19] = 10'd37;
        focus_phase_lut[39][20] = 10'd60;
        focus_phase_lut[39][21] = 10'd37;
        focus_phase_lut[39][22] = 10'd30;
        focus_phase_lut[39][23] = 10'd37;
        focus_phase_lut[39][24] = 10'd60;
        // focus_dist=2100mm (idx=40)
        focus_phase_lut[40][ 0] = 10'd58;
        focus_phase_lut[40][ 1] = 10'd36;
        focus_phase_lut[40][ 2] = 10'd29;
        focus_phase_lut[40][ 3] = 10'd36;
        focus_phase_lut[40][ 4] = 10'd58;
        focus_phase_lut[40][ 5] = 10'd36;
        focus_phase_lut[40][ 6] = 10'd15;
        focus_phase_lut[40][ 7] = 10'd7;
        focus_phase_lut[40][ 8] = 10'd15;
        focus_phase_lut[40][ 9] = 10'd36;
        focus_phase_lut[40][10] = 10'd29;
        focus_phase_lut[40][11] = 10'd7;
        focus_phase_lut[40][12] = 10'd0;
        focus_phase_lut[40][13] = 10'd7;
        focus_phase_lut[40][14] = 10'd29;
        focus_phase_lut[40][15] = 10'd36;
        focus_phase_lut[40][16] = 10'd15;
        focus_phase_lut[40][17] = 10'd7;
        focus_phase_lut[40][18] = 10'd15;
        focus_phase_lut[40][19] = 10'd36;
        focus_phase_lut[40][20] = 10'd58;
        focus_phase_lut[40][21] = 10'd36;
        focus_phase_lut[40][22] = 10'd29;
        focus_phase_lut[40][23] = 10'd36;
        focus_phase_lut[40][24] = 10'd58;
        // focus_dist=2150mm (idx=41)
        focus_phase_lut[41][ 0] = 10'd57;
        focus_phase_lut[41][ 1] = 10'd36;
        focus_phase_lut[41][ 2] = 10'd28;
        focus_phase_lut[41][ 3] = 10'd36;
        focus_phase_lut[41][ 4] = 10'd57;
        focus_phase_lut[41][ 5] = 10'd36;
        focus_phase_lut[41][ 6] = 10'd14;
        focus_phase_lut[41][ 7] = 10'd7;
        focus_phase_lut[41][ 8] = 10'd14;
        focus_phase_lut[41][ 9] = 10'd36;
        focus_phase_lut[41][10] = 10'd28;
        focus_phase_lut[41][11] = 10'd7;
        focus_phase_lut[41][12] = 10'd0;
        focus_phase_lut[41][13] = 10'd7;
        focus_phase_lut[41][14] = 10'd28;
        focus_phase_lut[41][15] = 10'd36;
        focus_phase_lut[41][16] = 10'd14;
        focus_phase_lut[41][17] = 10'd7;
        focus_phase_lut[41][18] = 10'd14;
        focus_phase_lut[41][19] = 10'd36;
        focus_phase_lut[41][20] = 10'd57;
        focus_phase_lut[41][21] = 10'd36;
        focus_phase_lut[41][22] = 10'd28;
        focus_phase_lut[41][23] = 10'd36;
        focus_phase_lut[41][24] = 10'd57;
        // focus_dist=2200mm (idx=42)
        focus_phase_lut[42][ 0] = 10'd56;
        focus_phase_lut[42][ 1] = 10'd35;
        focus_phase_lut[42][ 2] = 10'd28;
        focus_phase_lut[42][ 3] = 10'd35;
        focus_phase_lut[42][ 4] = 10'd56;
        focus_phase_lut[42][ 5] = 10'd35;
        focus_phase_lut[42][ 6] = 10'd14;
        focus_phase_lut[42][ 7] = 10'd7;
        focus_phase_lut[42][ 8] = 10'd14;
        focus_phase_lut[42][ 9] = 10'd35;
        focus_phase_lut[42][10] = 10'd28;
        focus_phase_lut[42][11] = 10'd7;
        focus_phase_lut[42][12] = 10'd0;
        focus_phase_lut[42][13] = 10'd7;
        focus_phase_lut[42][14] = 10'd28;
        focus_phase_lut[42][15] = 10'd35;
        focus_phase_lut[42][16] = 10'd14;
        focus_phase_lut[42][17] = 10'd7;
        focus_phase_lut[42][18] = 10'd14;
        focus_phase_lut[42][19] = 10'd35;
        focus_phase_lut[42][20] = 10'd56;
        focus_phase_lut[42][21] = 10'd35;
        focus_phase_lut[42][22] = 10'd28;
        focus_phase_lut[42][23] = 10'd35;
        focus_phase_lut[42][24] = 10'd56;
        // focus_dist=2250mm (idx=43)
        focus_phase_lut[43][ 0] = 10'd54;
        focus_phase_lut[43][ 1] = 10'd34;
        focus_phase_lut[43][ 2] = 10'd27;
        focus_phase_lut[43][ 3] = 10'd34;
        focus_phase_lut[43][ 4] = 10'd54;
        focus_phase_lut[43][ 5] = 10'd34;
        focus_phase_lut[43][ 6] = 10'd14;
        focus_phase_lut[43][ 7] = 10'd7;
        focus_phase_lut[43][ 8] = 10'd14;
        focus_phase_lut[43][ 9] = 10'd34;
        focus_phase_lut[43][10] = 10'd27;
        focus_phase_lut[43][11] = 10'd7;
        focus_phase_lut[43][12] = 10'd0;
        focus_phase_lut[43][13] = 10'd7;
        focus_phase_lut[43][14] = 10'd27;
        focus_phase_lut[43][15] = 10'd34;
        focus_phase_lut[43][16] = 10'd14;
        focus_phase_lut[43][17] = 10'd7;
        focus_phase_lut[43][18] = 10'd14;
        focus_phase_lut[43][19] = 10'd34;
        focus_phase_lut[43][20] = 10'd54;
        focus_phase_lut[43][21] = 10'd34;
        focus_phase_lut[43][22] = 10'd27;
        focus_phase_lut[43][23] = 10'd34;
        focus_phase_lut[43][24] = 10'd54;
        // focus_dist=2300mm (idx=44)
        focus_phase_lut[44][ 0] = 10'd53;
        focus_phase_lut[44][ 1] = 10'd33;
        focus_phase_lut[44][ 2] = 10'd27;
        focus_phase_lut[44][ 3] = 10'd33;
        focus_phase_lut[44][ 4] = 10'd53;
        focus_phase_lut[44][ 5] = 10'd33;
        focus_phase_lut[44][ 6] = 10'd13;
        focus_phase_lut[44][ 7] = 10'd7;
        focus_phase_lut[44][ 8] = 10'd13;
        focus_phase_lut[44][ 9] = 10'd33;
        focus_phase_lut[44][10] = 10'd27;
        focus_phase_lut[44][11] = 10'd7;
        focus_phase_lut[44][12] = 10'd0;
        focus_phase_lut[44][13] = 10'd7;
        focus_phase_lut[44][14] = 10'd27;
        focus_phase_lut[44][15] = 10'd33;
        focus_phase_lut[44][16] = 10'd13;
        focus_phase_lut[44][17] = 10'd7;
        focus_phase_lut[44][18] = 10'd13;
        focus_phase_lut[44][19] = 10'd33;
        focus_phase_lut[44][20] = 10'd53;
        focus_phase_lut[44][21] = 10'd33;
        focus_phase_lut[44][22] = 10'd27;
        focus_phase_lut[44][23] = 10'd33;
        focus_phase_lut[44][24] = 10'd53;
        // focus_dist=2350mm (idx=45)
        focus_phase_lut[45][ 0] = 10'd52;
        focus_phase_lut[45][ 1] = 10'd33;
        focus_phase_lut[45][ 2] = 10'd26;
        focus_phase_lut[45][ 3] = 10'd33;
        focus_phase_lut[45][ 4] = 10'd52;
        focus_phase_lut[45][ 5] = 10'd33;
        focus_phase_lut[45][ 6] = 10'd13;
        focus_phase_lut[45][ 7] = 10'd7;
        focus_phase_lut[45][ 8] = 10'd13;
        focus_phase_lut[45][ 9] = 10'd33;
        focus_phase_lut[45][10] = 10'd26;
        focus_phase_lut[45][11] = 10'd7;
        focus_phase_lut[45][12] = 10'd0;
        focus_phase_lut[45][13] = 10'd7;
        focus_phase_lut[45][14] = 10'd26;
        focus_phase_lut[45][15] = 10'd33;
        focus_phase_lut[45][16] = 10'd13;
        focus_phase_lut[45][17] = 10'd7;
        focus_phase_lut[45][18] = 10'd13;
        focus_phase_lut[45][19] = 10'd33;
        focus_phase_lut[45][20] = 10'd52;
        focus_phase_lut[45][21] = 10'd33;
        focus_phase_lut[45][22] = 10'd26;
        focus_phase_lut[45][23] = 10'd33;
        focus_phase_lut[45][24] = 10'd52;
        // focus_dist=2400mm (idx=46)
        focus_phase_lut[46][ 0] = 10'd51;
        focus_phase_lut[46][ 1] = 10'd32;
        focus_phase_lut[46][ 2] = 10'd25;
        focus_phase_lut[46][ 3] = 10'd32;
        focus_phase_lut[46][ 4] = 10'd51;
        focus_phase_lut[46][ 5] = 10'd32;
        focus_phase_lut[46][ 6] = 10'd13;
        focus_phase_lut[46][ 7] = 10'd6;
        focus_phase_lut[46][ 8] = 10'd13;
        focus_phase_lut[46][ 9] = 10'd32;
        focus_phase_lut[46][10] = 10'd25;
        focus_phase_lut[46][11] = 10'd6;
        focus_phase_lut[46][12] = 10'd0;
        focus_phase_lut[46][13] = 10'd6;
        focus_phase_lut[46][14] = 10'd25;
        focus_phase_lut[46][15] = 10'd32;
        focus_phase_lut[46][16] = 10'd13;
        focus_phase_lut[46][17] = 10'd6;
        focus_phase_lut[46][18] = 10'd13;
        focus_phase_lut[46][19] = 10'd32;
        focus_phase_lut[46][20] = 10'd51;
        focus_phase_lut[46][21] = 10'd32;
        focus_phase_lut[46][22] = 10'd25;
        focus_phase_lut[46][23] = 10'd32;
        focus_phase_lut[46][24] = 10'd51;
        // focus_dist=2450mm (idx=47)
        focus_phase_lut[47][ 0] = 10'd50;
        focus_phase_lut[47][ 1] = 10'd31;
        focus_phase_lut[47][ 2] = 10'd25;
        focus_phase_lut[47][ 3] = 10'd31;
        focus_phase_lut[47][ 4] = 10'd50;
        focus_phase_lut[47][ 5] = 10'd31;
        focus_phase_lut[47][ 6] = 10'd12;
        focus_phase_lut[47][ 7] = 10'd6;
        focus_phase_lut[47][ 8] = 10'd12;
        focus_phase_lut[47][ 9] = 10'd31;
        focus_phase_lut[47][10] = 10'd25;
        focus_phase_lut[47][11] = 10'd6;
        focus_phase_lut[47][12] = 10'd0;
        focus_phase_lut[47][13] = 10'd6;
        focus_phase_lut[47][14] = 10'd25;
        focus_phase_lut[47][15] = 10'd31;
        focus_phase_lut[47][16] = 10'd12;
        focus_phase_lut[47][17] = 10'd6;
        focus_phase_lut[47][18] = 10'd12;
        focus_phase_lut[47][19] = 10'd31;
        focus_phase_lut[47][20] = 10'd50;
        focus_phase_lut[47][21] = 10'd31;
        focus_phase_lut[47][22] = 10'd25;
        focus_phase_lut[47][23] = 10'd31;
        focus_phase_lut[47][24] = 10'd50;
        // focus_dist=2500mm (idx=48)
        focus_phase_lut[48][ 0] = 10'd49;
        focus_phase_lut[48][ 1] = 10'd31;
        focus_phase_lut[48][ 2] = 10'd24;
        focus_phase_lut[48][ 3] = 10'd31;
        focus_phase_lut[48][ 4] = 10'd49;
        focus_phase_lut[48][ 5] = 10'd31;
        focus_phase_lut[48][ 6] = 10'd12;
        focus_phase_lut[48][ 7] = 10'd6;
        focus_phase_lut[48][ 8] = 10'd12;
        focus_phase_lut[48][ 9] = 10'd31;
        focus_phase_lut[48][10] = 10'd24;
        focus_phase_lut[48][11] = 10'd6;
        focus_phase_lut[48][12] = 10'd0;
        focus_phase_lut[48][13] = 10'd6;
        focus_phase_lut[48][14] = 10'd24;
        focus_phase_lut[48][15] = 10'd31;
        focus_phase_lut[48][16] = 10'd12;
        focus_phase_lut[48][17] = 10'd6;
        focus_phase_lut[48][18] = 10'd12;
        focus_phase_lut[48][19] = 10'd31;
        focus_phase_lut[48][20] = 10'd49;
        focus_phase_lut[48][21] = 10'd31;
        focus_phase_lut[48][22] = 10'd24;
        focus_phase_lut[48][23] = 10'd31;
        focus_phase_lut[48][24] = 10'd49;
        // focus_dist=2550mm (idx=49)
        focus_phase_lut[49][ 0] = 10'd48;
        focus_phase_lut[49][ 1] = 10'd30;
        focus_phase_lut[49][ 2] = 10'd24;
        focus_phase_lut[49][ 3] = 10'd30;
        focus_phase_lut[49][ 4] = 10'd48;
        focus_phase_lut[49][ 5] = 10'd30;
        focus_phase_lut[49][ 6] = 10'd12;
        focus_phase_lut[49][ 7] = 10'd6;
        focus_phase_lut[49][ 8] = 10'd12;
        focus_phase_lut[49][ 9] = 10'd30;
        focus_phase_lut[49][10] = 10'd24;
        focus_phase_lut[49][11] = 10'd6;
        focus_phase_lut[49][12] = 10'd0;
        focus_phase_lut[49][13] = 10'd6;
        focus_phase_lut[49][14] = 10'd24;
        focus_phase_lut[49][15] = 10'd30;
        focus_phase_lut[49][16] = 10'd12;
        focus_phase_lut[49][17] = 10'd6;
        focus_phase_lut[49][18] = 10'd12;
        focus_phase_lut[49][19] = 10'd30;
        focus_phase_lut[49][20] = 10'd48;
        focus_phase_lut[49][21] = 10'd30;
        focus_phase_lut[49][22] = 10'd24;
        focus_phase_lut[49][23] = 10'd30;
        focus_phase_lut[49][24] = 10'd48;
        // focus_dist=2600mm (idx=50)
        focus_phase_lut[50][ 0] = 10'd47;
        focus_phase_lut[50][ 1] = 10'd29;
        focus_phase_lut[50][ 2] = 10'd24;
        focus_phase_lut[50][ 3] = 10'd29;
        focus_phase_lut[50][ 4] = 10'd47;
        focus_phase_lut[50][ 5] = 10'd29;
        focus_phase_lut[50][ 6] = 10'd12;
        focus_phase_lut[50][ 7] = 10'd6;
        focus_phase_lut[50][ 8] = 10'd12;
        focus_phase_lut[50][ 9] = 10'd29;
        focus_phase_lut[50][10] = 10'd24;
        focus_phase_lut[50][11] = 10'd6;
        focus_phase_lut[50][12] = 10'd0;
        focus_phase_lut[50][13] = 10'd6;
        focus_phase_lut[50][14] = 10'd24;
        focus_phase_lut[50][15] = 10'd29;
        focus_phase_lut[50][16] = 10'd12;
        focus_phase_lut[50][17] = 10'd6;
        focus_phase_lut[50][18] = 10'd12;
        focus_phase_lut[50][19] = 10'd29;
        focus_phase_lut[50][20] = 10'd47;
        focus_phase_lut[50][21] = 10'd29;
        focus_phase_lut[50][22] = 10'd24;
        focus_phase_lut[50][23] = 10'd29;
        focus_phase_lut[50][24] = 10'd47;
        // focus_dist=2650mm (idx=51)
        focus_phase_lut[51][ 0] = 10'd46;
        focus_phase_lut[51][ 1] = 10'd29;
        focus_phase_lut[51][ 2] = 10'd23;
        focus_phase_lut[51][ 3] = 10'd29;
        focus_phase_lut[51][ 4] = 10'd46;
        focus_phase_lut[51][ 5] = 10'd29;
        focus_phase_lut[51][ 6] = 10'd12;
        focus_phase_lut[51][ 7] = 10'd6;
        focus_phase_lut[51][ 8] = 10'd12;
        focus_phase_lut[51][ 9] = 10'd29;
        focus_phase_lut[51][10] = 10'd23;
        focus_phase_lut[51][11] = 10'd6;
        focus_phase_lut[51][12] = 10'd0;
        focus_phase_lut[51][13] = 10'd6;
        focus_phase_lut[51][14] = 10'd23;
        focus_phase_lut[51][15] = 10'd29;
        focus_phase_lut[51][16] = 10'd12;
        focus_phase_lut[51][17] = 10'd6;
        focus_phase_lut[51][18] = 10'd12;
        focus_phase_lut[51][19] = 10'd29;
        focus_phase_lut[51][20] = 10'd46;
        focus_phase_lut[51][21] = 10'd29;
        focus_phase_lut[51][22] = 10'd23;
        focus_phase_lut[51][23] = 10'd29;
        focus_phase_lut[51][24] = 10'd46;
        // focus_dist=2700mm (idx=52)
        focus_phase_lut[52][ 0] = 10'd45;
        focus_phase_lut[52][ 1] = 10'd28;
        focus_phase_lut[52][ 2] = 10'd23;
        focus_phase_lut[52][ 3] = 10'd28;
        focus_phase_lut[52][ 4] = 10'd45;
        focus_phase_lut[52][ 5] = 10'd28;
        focus_phase_lut[52][ 6] = 10'd11;
        focus_phase_lut[52][ 7] = 10'd6;
        focus_phase_lut[52][ 8] = 10'd11;
        focus_phase_lut[52][ 9] = 10'd28;
        focus_phase_lut[52][10] = 10'd23;
        focus_phase_lut[52][11] = 10'd6;
        focus_phase_lut[52][12] = 10'd0;
        focus_phase_lut[52][13] = 10'd6;
        focus_phase_lut[52][14] = 10'd23;
        focus_phase_lut[52][15] = 10'd28;
        focus_phase_lut[52][16] = 10'd11;
        focus_phase_lut[52][17] = 10'd6;
        focus_phase_lut[52][18] = 10'd11;
        focus_phase_lut[52][19] = 10'd28;
        focus_phase_lut[52][20] = 10'd45;
        focus_phase_lut[52][21] = 10'd28;
        focus_phase_lut[52][22] = 10'd23;
        focus_phase_lut[52][23] = 10'd28;
        focus_phase_lut[52][24] = 10'd45;
        // focus_dist=2750mm (idx=53)
        focus_phase_lut[53][ 0] = 10'd44;
        focus_phase_lut[53][ 1] = 10'd28;
        focus_phase_lut[53][ 2] = 10'd22;
        focus_phase_lut[53][ 3] = 10'd28;
        focus_phase_lut[53][ 4] = 10'd44;
        focus_phase_lut[53][ 5] = 10'd28;
        focus_phase_lut[53][ 6] = 10'd11;
        focus_phase_lut[53][ 7] = 10'd6;
        focus_phase_lut[53][ 8] = 10'd11;
        focus_phase_lut[53][ 9] = 10'd28;
        focus_phase_lut[53][10] = 10'd22;
        focus_phase_lut[53][11] = 10'd6;
        focus_phase_lut[53][12] = 10'd0;
        focus_phase_lut[53][13] = 10'd6;
        focus_phase_lut[53][14] = 10'd22;
        focus_phase_lut[53][15] = 10'd28;
        focus_phase_lut[53][16] = 10'd11;
        focus_phase_lut[53][17] = 10'd6;
        focus_phase_lut[53][18] = 10'd11;
        focus_phase_lut[53][19] = 10'd28;
        focus_phase_lut[53][20] = 10'd44;
        focus_phase_lut[53][21] = 10'd28;
        focus_phase_lut[53][22] = 10'd22;
        focus_phase_lut[53][23] = 10'd28;
        focus_phase_lut[53][24] = 10'd44;
        // focus_dist=2800mm (idx=54)
        focus_phase_lut[54][ 0] = 10'd44;
        focus_phase_lut[54][ 1] = 10'd27;
        focus_phase_lut[54][ 2] = 10'd22;
        focus_phase_lut[54][ 3] = 10'd27;
        focus_phase_lut[54][ 4] = 10'd44;
        focus_phase_lut[54][ 5] = 10'd27;
        focus_phase_lut[54][ 6] = 10'd11;
        focus_phase_lut[54][ 7] = 10'd5;
        focus_phase_lut[54][ 8] = 10'd11;
        focus_phase_lut[54][ 9] = 10'd27;
        focus_phase_lut[54][10] = 10'd22;
        focus_phase_lut[54][11] = 10'd5;
        focus_phase_lut[54][12] = 10'd0;
        focus_phase_lut[54][13] = 10'd5;
        focus_phase_lut[54][14] = 10'd22;
        focus_phase_lut[54][15] = 10'd27;
        focus_phase_lut[54][16] = 10'd11;
        focus_phase_lut[54][17] = 10'd5;
        focus_phase_lut[54][18] = 10'd11;
        focus_phase_lut[54][19] = 10'd27;
        focus_phase_lut[54][20] = 10'd44;
        focus_phase_lut[54][21] = 10'd27;
        focus_phase_lut[54][22] = 10'd22;
        focus_phase_lut[54][23] = 10'd27;
        focus_phase_lut[54][24] = 10'd44;
        // focus_dist=2850mm (idx=55)
        focus_phase_lut[55][ 0] = 10'd43;
        focus_phase_lut[55][ 1] = 10'd27;
        focus_phase_lut[55][ 2] = 10'd21;
        focus_phase_lut[55][ 3] = 10'd27;
        focus_phase_lut[55][ 4] = 10'd43;
        focus_phase_lut[55][ 5] = 10'd27;
        focus_phase_lut[55][ 6] = 10'd11;
        focus_phase_lut[55][ 7] = 10'd5;
        focus_phase_lut[55][ 8] = 10'd11;
        focus_phase_lut[55][ 9] = 10'd27;
        focus_phase_lut[55][10] = 10'd21;
        focus_phase_lut[55][11] = 10'd5;
        focus_phase_lut[55][12] = 10'd0;
        focus_phase_lut[55][13] = 10'd5;
        focus_phase_lut[55][14] = 10'd21;
        focus_phase_lut[55][15] = 10'd27;
        focus_phase_lut[55][16] = 10'd11;
        focus_phase_lut[55][17] = 10'd5;
        focus_phase_lut[55][18] = 10'd11;
        focus_phase_lut[55][19] = 10'd27;
        focus_phase_lut[55][20] = 10'd43;
        focus_phase_lut[55][21] = 10'd27;
        focus_phase_lut[55][22] = 10'd21;
        focus_phase_lut[55][23] = 10'd27;
        focus_phase_lut[55][24] = 10'd43;
        // focus_dist=2900mm (idx=56)
        focus_phase_lut[56][ 0] = 10'd42;
        focus_phase_lut[56][ 1] = 10'd26;
        focus_phase_lut[56][ 2] = 10'd21;
        focus_phase_lut[56][ 3] = 10'd26;
        focus_phase_lut[56][ 4] = 10'd42;
        focus_phase_lut[56][ 5] = 10'd26;
        focus_phase_lut[56][ 6] = 10'd11;
        focus_phase_lut[56][ 7] = 10'd5;
        focus_phase_lut[56][ 8] = 10'd11;
        focus_phase_lut[56][ 9] = 10'd26;
        focus_phase_lut[56][10] = 10'd21;
        focus_phase_lut[56][11] = 10'd5;
        focus_phase_lut[56][12] = 10'd0;
        focus_phase_lut[56][13] = 10'd5;
        focus_phase_lut[56][14] = 10'd21;
        focus_phase_lut[56][15] = 10'd26;
        focus_phase_lut[56][16] = 10'd11;
        focus_phase_lut[56][17] = 10'd5;
        focus_phase_lut[56][18] = 10'd11;
        focus_phase_lut[56][19] = 10'd26;
        focus_phase_lut[56][20] = 10'd42;
        focus_phase_lut[56][21] = 10'd26;
        focus_phase_lut[56][22] = 10'd21;
        focus_phase_lut[56][23] = 10'd26;
        focus_phase_lut[56][24] = 10'd42;
        // focus_dist=2950mm (idx=57)
        focus_phase_lut[57][ 0] = 10'd41;
        focus_phase_lut[57][ 1] = 10'd26;
        focus_phase_lut[57][ 2] = 10'd21;
        focus_phase_lut[57][ 3] = 10'd26;
        focus_phase_lut[57][ 4] = 10'd41;
        focus_phase_lut[57][ 5] = 10'd26;
        focus_phase_lut[57][ 6] = 10'd10;
        focus_phase_lut[57][ 7] = 10'd5;
        focus_phase_lut[57][ 8] = 10'd10;
        focus_phase_lut[57][ 9] = 10'd26;
        focus_phase_lut[57][10] = 10'd21;
        focus_phase_lut[57][11] = 10'd5;
        focus_phase_lut[57][12] = 10'd0;
        focus_phase_lut[57][13] = 10'd5;
        focus_phase_lut[57][14] = 10'd21;
        focus_phase_lut[57][15] = 10'd26;
        focus_phase_lut[57][16] = 10'd10;
        focus_phase_lut[57][17] = 10'd5;
        focus_phase_lut[57][18] = 10'd10;
        focus_phase_lut[57][19] = 10'd26;
        focus_phase_lut[57][20] = 10'd41;
        focus_phase_lut[57][21] = 10'd26;
        focus_phase_lut[57][22] = 10'd21;
        focus_phase_lut[57][23] = 10'd26;
        focus_phase_lut[57][24] = 10'd41;
        // focus_dist=3000mm (idx=58)
        focus_phase_lut[58][ 0] = 10'd41;
        focus_phase_lut[58][ 1] = 10'd25;
        focus_phase_lut[58][ 2] = 10'd20;
        focus_phase_lut[58][ 3] = 10'd25;
        focus_phase_lut[58][ 4] = 10'd41;
        focus_phase_lut[58][ 5] = 10'd25;
        focus_phase_lut[58][ 6] = 10'd10;
        focus_phase_lut[58][ 7] = 10'd5;
        focus_phase_lut[58][ 8] = 10'd10;
        focus_phase_lut[58][ 9] = 10'd25;
        focus_phase_lut[58][10] = 10'd20;
        focus_phase_lut[58][11] = 10'd5;
        focus_phase_lut[58][12] = 10'd0;
        focus_phase_lut[58][13] = 10'd5;
        focus_phase_lut[58][14] = 10'd20;
        focus_phase_lut[58][15] = 10'd25;
        focus_phase_lut[58][16] = 10'd10;
        focus_phase_lut[58][17] = 10'd5;
        focus_phase_lut[58][18] = 10'd10;
        focus_phase_lut[58][19] = 10'd25;
        focus_phase_lut[58][20] = 10'd41;
        focus_phase_lut[58][21] = 10'd25;
        focus_phase_lut[58][22] = 10'd20;
        focus_phase_lut[58][23] = 10'd25;
        focus_phase_lut[58][24] = 10'd41;
        // focus_dist=3050mm (idx=59)
        focus_phase_lut[59][ 0] = 10'd40;
        focus_phase_lut[59][ 1] = 10'd25;
        focus_phase_lut[59][ 2] = 10'd20;
        focus_phase_lut[59][ 3] = 10'd25;
        focus_phase_lut[59][ 4] = 10'd40;
        focus_phase_lut[59][ 5] = 10'd25;
        focus_phase_lut[59][ 6] = 10'd10;
        focus_phase_lut[59][ 7] = 10'd5;
        focus_phase_lut[59][ 8] = 10'd10;
        focus_phase_lut[59][ 9] = 10'd25;
        focus_phase_lut[59][10] = 10'd20;
        focus_phase_lut[59][11] = 10'd5;
        focus_phase_lut[59][12] = 10'd0;
        focus_phase_lut[59][13] = 10'd5;
        focus_phase_lut[59][14] = 10'd20;
        focus_phase_lut[59][15] = 10'd25;
        focus_phase_lut[59][16] = 10'd10;
        focus_phase_lut[59][17] = 10'd5;
        focus_phase_lut[59][18] = 10'd10;
        focus_phase_lut[59][19] = 10'd25;
        focus_phase_lut[59][20] = 10'd40;
        focus_phase_lut[59][21] = 10'd25;
        focus_phase_lut[59][22] = 10'd20;
        focus_phase_lut[59][23] = 10'd25;
        focus_phase_lut[59][24] = 10'd40;
        // focus_dist=3100mm (idx=60)
        focus_phase_lut[60][ 0] = 10'd39;
        focus_phase_lut[60][ 1] = 10'd25;
        focus_phase_lut[60][ 2] = 10'd20;
        focus_phase_lut[60][ 3] = 10'd25;
        focus_phase_lut[60][ 4] = 10'd39;
        focus_phase_lut[60][ 5] = 10'd25;
        focus_phase_lut[60][ 6] = 10'd10;
        focus_phase_lut[60][ 7] = 10'd5;
        focus_phase_lut[60][ 8] = 10'd10;
        focus_phase_lut[60][ 9] = 10'd25;
        focus_phase_lut[60][10] = 10'd20;
        focus_phase_lut[60][11] = 10'd5;
        focus_phase_lut[60][12] = 10'd0;
        focus_phase_lut[60][13] = 10'd5;
        focus_phase_lut[60][14] = 10'd20;
        focus_phase_lut[60][15] = 10'd25;
        focus_phase_lut[60][16] = 10'd10;
        focus_phase_lut[60][17] = 10'd5;
        focus_phase_lut[60][18] = 10'd10;
        focus_phase_lut[60][19] = 10'd25;
        focus_phase_lut[60][20] = 10'd39;
        focus_phase_lut[60][21] = 10'd25;
        focus_phase_lut[60][22] = 10'd20;
        focus_phase_lut[60][23] = 10'd25;
        focus_phase_lut[60][24] = 10'd39;
        // focus_dist=3150mm (idx=61)
        focus_phase_lut[61][ 0] = 10'd39;
        focus_phase_lut[61][ 1] = 10'd24;
        focus_phase_lut[61][ 2] = 10'd19;
        focus_phase_lut[61][ 3] = 10'd24;
        focus_phase_lut[61][ 4] = 10'd39;
        focus_phase_lut[61][ 5] = 10'd24;
        focus_phase_lut[61][ 6] = 10'd10;
        focus_phase_lut[61][ 7] = 10'd5;
        focus_phase_lut[61][ 8] = 10'd10;
        focus_phase_lut[61][ 9] = 10'd24;
        focus_phase_lut[61][10] = 10'd19;
        focus_phase_lut[61][11] = 10'd5;
        focus_phase_lut[61][12] = 10'd0;
        focus_phase_lut[61][13] = 10'd5;
        focus_phase_lut[61][14] = 10'd19;
        focus_phase_lut[61][15] = 10'd24;
        focus_phase_lut[61][16] = 10'd10;
        focus_phase_lut[61][17] = 10'd5;
        focus_phase_lut[61][18] = 10'd10;
        focus_phase_lut[61][19] = 10'd24;
        focus_phase_lut[61][20] = 10'd39;
        focus_phase_lut[61][21] = 10'd24;
        focus_phase_lut[61][22] = 10'd19;
        focus_phase_lut[61][23] = 10'd24;
        focus_phase_lut[61][24] = 10'd39;
        // focus_dist=3200mm (idx=62)
        focus_phase_lut[62][ 0] = 10'd38;
        focus_phase_lut[62][ 1] = 10'd24;
        focus_phase_lut[62][ 2] = 10'd19;
        focus_phase_lut[62][ 3] = 10'd24;
        focus_phase_lut[62][ 4] = 10'd38;
        focus_phase_lut[62][ 5] = 10'd24;
        focus_phase_lut[62][ 6] = 10'd10;
        focus_phase_lut[62][ 7] = 10'd5;
        focus_phase_lut[62][ 8] = 10'd10;
        focus_phase_lut[62][ 9] = 10'd24;
        focus_phase_lut[62][10] = 10'd19;
        focus_phase_lut[62][11] = 10'd5;
        focus_phase_lut[62][12] = 10'd0;
        focus_phase_lut[62][13] = 10'd5;
        focus_phase_lut[62][14] = 10'd19;
        focus_phase_lut[62][15] = 10'd24;
        focus_phase_lut[62][16] = 10'd10;
        focus_phase_lut[62][17] = 10'd5;
        focus_phase_lut[62][18] = 10'd10;
        focus_phase_lut[62][19] = 10'd24;
        focus_phase_lut[62][20] = 10'd38;
        focus_phase_lut[62][21] = 10'd24;
        focus_phase_lut[62][22] = 10'd19;
        focus_phase_lut[62][23] = 10'd24;
        focus_phase_lut[62][24] = 10'd38;
        // focus_dist=3250mm (idx=63)
        focus_phase_lut[63][ 0] = 10'd38;
        focus_phase_lut[63][ 1] = 10'd24;
        focus_phase_lut[63][ 2] = 10'd19;
        focus_phase_lut[63][ 3] = 10'd24;
        focus_phase_lut[63][ 4] = 10'd38;
        focus_phase_lut[63][ 5] = 10'd24;
        focus_phase_lut[63][ 6] = 10'd9;
        focus_phase_lut[63][ 7] = 10'd5;
        focus_phase_lut[63][ 8] = 10'd9;
        focus_phase_lut[63][ 9] = 10'd24;
        focus_phase_lut[63][10] = 10'd19;
        focus_phase_lut[63][11] = 10'd5;
        focus_phase_lut[63][12] = 10'd0;
        focus_phase_lut[63][13] = 10'd5;
        focus_phase_lut[63][14] = 10'd19;
        focus_phase_lut[63][15] = 10'd24;
        focus_phase_lut[63][16] = 10'd9;
        focus_phase_lut[63][17] = 10'd5;
        focus_phase_lut[63][18] = 10'd9;
        focus_phase_lut[63][19] = 10'd24;
        focus_phase_lut[63][20] = 10'd38;
        focus_phase_lut[63][21] = 10'd24;
        focus_phase_lut[63][22] = 10'd19;
        focus_phase_lut[63][23] = 10'd24;
        focus_phase_lut[63][24] = 10'd38;
        // focus_dist=3300mm (idx=64)
        focus_phase_lut[64][ 0] = 10'd37;
        focus_phase_lut[64][ 1] = 10'd23;
        focus_phase_lut[64][ 2] = 10'd19;
        focus_phase_lut[64][ 3] = 10'd23;
        focus_phase_lut[64][ 4] = 10'd37;
        focus_phase_lut[64][ 5] = 10'd23;
        focus_phase_lut[64][ 6] = 10'd9;
        focus_phase_lut[64][ 7] = 10'd5;
        focus_phase_lut[64][ 8] = 10'd9;
        focus_phase_lut[64][ 9] = 10'd23;
        focus_phase_lut[64][10] = 10'd19;
        focus_phase_lut[64][11] = 10'd5;
        focus_phase_lut[64][12] = 10'd0;
        focus_phase_lut[64][13] = 10'd5;
        focus_phase_lut[64][14] = 10'd19;
        focus_phase_lut[64][15] = 10'd23;
        focus_phase_lut[64][16] = 10'd9;
        focus_phase_lut[64][17] = 10'd5;
        focus_phase_lut[64][18] = 10'd9;
        focus_phase_lut[64][19] = 10'd23;
        focus_phase_lut[64][20] = 10'd37;
        focus_phase_lut[64][21] = 10'd23;
        focus_phase_lut[64][22] = 10'd19;
        focus_phase_lut[64][23] = 10'd23;
        focus_phase_lut[64][24] = 10'd37;
        // focus_dist=3350mm (idx=65)
        focus_phase_lut[65][ 0] = 10'd37;
        focus_phase_lut[65][ 1] = 10'd23;
        focus_phase_lut[65][ 2] = 10'd18;
        focus_phase_lut[65][ 3] = 10'd23;
        focus_phase_lut[65][ 4] = 10'd37;
        focus_phase_lut[65][ 5] = 10'd23;
        focus_phase_lut[65][ 6] = 10'd9;
        focus_phase_lut[65][ 7] = 10'd5;
        focus_phase_lut[65][ 8] = 10'd9;
        focus_phase_lut[65][ 9] = 10'd23;
        focus_phase_lut[65][10] = 10'd18;
        focus_phase_lut[65][11] = 10'd5;
        focus_phase_lut[65][12] = 10'd0;
        focus_phase_lut[65][13] = 10'd5;
        focus_phase_lut[65][14] = 10'd18;
        focus_phase_lut[65][15] = 10'd23;
        focus_phase_lut[65][16] = 10'd9;
        focus_phase_lut[65][17] = 10'd5;
        focus_phase_lut[65][18] = 10'd9;
        focus_phase_lut[65][19] = 10'd23;
        focus_phase_lut[65][20] = 10'd37;
        focus_phase_lut[65][21] = 10'd23;
        focus_phase_lut[65][22] = 10'd18;
        focus_phase_lut[65][23] = 10'd23;
        focus_phase_lut[65][24] = 10'd37;
        // focus_dist=3400mm (idx=66)
        focus_phase_lut[66][ 0] = 10'd36;
        focus_phase_lut[66][ 1] = 10'd22;
        focus_phase_lut[66][ 2] = 10'd18;
        focus_phase_lut[66][ 3] = 10'd22;
        focus_phase_lut[66][ 4] = 10'd36;
        focus_phase_lut[66][ 5] = 10'd22;
        focus_phase_lut[66][ 6] = 10'd9;
        focus_phase_lut[66][ 7] = 10'd4;
        focus_phase_lut[66][ 8] = 10'd9;
        focus_phase_lut[66][ 9] = 10'd22;
        focus_phase_lut[66][10] = 10'd18;
        focus_phase_lut[66][11] = 10'd4;
        focus_phase_lut[66][12] = 10'd0;
        focus_phase_lut[66][13] = 10'd4;
        focus_phase_lut[66][14] = 10'd18;
        focus_phase_lut[66][15] = 10'd22;
        focus_phase_lut[66][16] = 10'd9;
        focus_phase_lut[66][17] = 10'd4;
        focus_phase_lut[66][18] = 10'd9;
        focus_phase_lut[66][19] = 10'd22;
        focus_phase_lut[66][20] = 10'd36;
        focus_phase_lut[66][21] = 10'd22;
        focus_phase_lut[66][22] = 10'd18;
        focus_phase_lut[66][23] = 10'd22;
        focus_phase_lut[66][24] = 10'd36;
        // focus_dist=3450mm (idx=67)
        focus_phase_lut[67][ 0] = 10'd35;
        focus_phase_lut[67][ 1] = 10'd22;
        focus_phase_lut[67][ 2] = 10'd18;
        focus_phase_lut[67][ 3] = 10'd22;
        focus_phase_lut[67][ 4] = 10'd35;
        focus_phase_lut[67][ 5] = 10'd22;
        focus_phase_lut[67][ 6] = 10'd9;
        focus_phase_lut[67][ 7] = 10'd4;
        focus_phase_lut[67][ 8] = 10'd9;
        focus_phase_lut[67][ 9] = 10'd22;
        focus_phase_lut[67][10] = 10'd18;
        focus_phase_lut[67][11] = 10'd4;
        focus_phase_lut[67][12] = 10'd0;
        focus_phase_lut[67][13] = 10'd4;
        focus_phase_lut[67][14] = 10'd18;
        focus_phase_lut[67][15] = 10'd22;
        focus_phase_lut[67][16] = 10'd9;
        focus_phase_lut[67][17] = 10'd4;
        focus_phase_lut[67][18] = 10'd9;
        focus_phase_lut[67][19] = 10'd22;
        focus_phase_lut[67][20] = 10'd35;
        focus_phase_lut[67][21] = 10'd22;
        focus_phase_lut[67][22] = 10'd18;
        focus_phase_lut[67][23] = 10'd22;
        focus_phase_lut[67][24] = 10'd35;
        // focus_dist=3500mm (idx=68)
        focus_phase_lut[68][ 0] = 10'd35;
        focus_phase_lut[68][ 1] = 10'd22;
        focus_phase_lut[68][ 2] = 10'd17;
        focus_phase_lut[68][ 3] = 10'd22;
        focus_phase_lut[68][ 4] = 10'd35;
        focus_phase_lut[68][ 5] = 10'd22;
        focus_phase_lut[68][ 6] = 10'd9;
        focus_phase_lut[68][ 7] = 10'd4;
        focus_phase_lut[68][ 8] = 10'd9;
        focus_phase_lut[68][ 9] = 10'd22;
        focus_phase_lut[68][10] = 10'd17;
        focus_phase_lut[68][11] = 10'd4;
        focus_phase_lut[68][12] = 10'd0;
        focus_phase_lut[68][13] = 10'd4;
        focus_phase_lut[68][14] = 10'd17;
        focus_phase_lut[68][15] = 10'd22;
        focus_phase_lut[68][16] = 10'd9;
        focus_phase_lut[68][17] = 10'd4;
        focus_phase_lut[68][18] = 10'd9;
        focus_phase_lut[68][19] = 10'd22;
        focus_phase_lut[68][20] = 10'd35;
        focus_phase_lut[68][21] = 10'd22;
        focus_phase_lut[68][22] = 10'd17;
        focus_phase_lut[68][23] = 10'd22;
        focus_phase_lut[68][24] = 10'd35;
        // focus_dist=3550mm (idx=69)
        focus_phase_lut[69][ 0] = 10'd34;
        focus_phase_lut[69][ 1] = 10'd22;
        focus_phase_lut[69][ 2] = 10'd17;
        focus_phase_lut[69][ 3] = 10'd22;
        focus_phase_lut[69][ 4] = 10'd34;
        focus_phase_lut[69][ 5] = 10'd22;
        focus_phase_lut[69][ 6] = 10'd9;
        focus_phase_lut[69][ 7] = 10'd4;
        focus_phase_lut[69][ 8] = 10'd9;
        focus_phase_lut[69][ 9] = 10'd22;
        focus_phase_lut[69][10] = 10'd17;
        focus_phase_lut[69][11] = 10'd4;
        focus_phase_lut[69][12] = 10'd0;
        focus_phase_lut[69][13] = 10'd4;
        focus_phase_lut[69][14] = 10'd17;
        focus_phase_lut[69][15] = 10'd22;
        focus_phase_lut[69][16] = 10'd9;
        focus_phase_lut[69][17] = 10'd4;
        focus_phase_lut[69][18] = 10'd9;
        focus_phase_lut[69][19] = 10'd22;
        focus_phase_lut[69][20] = 10'd34;
        focus_phase_lut[69][21] = 10'd22;
        focus_phase_lut[69][22] = 10'd17;
        focus_phase_lut[69][23] = 10'd22;
        focus_phase_lut[69][24] = 10'd34;
        // focus_dist=3600mm (idx=70)
        focus_phase_lut[70][ 0] = 10'd34;
        focus_phase_lut[70][ 1] = 10'd21;
        focus_phase_lut[70][ 2] = 10'd17;
        focus_phase_lut[70][ 3] = 10'd21;
        focus_phase_lut[70][ 4] = 10'd34;
        focus_phase_lut[70][ 5] = 10'd21;
        focus_phase_lut[70][ 6] = 10'd8;
        focus_phase_lut[70][ 7] = 10'd4;
        focus_phase_lut[70][ 8] = 10'd8;
        focus_phase_lut[70][ 9] = 10'd21;
        focus_phase_lut[70][10] = 10'd17;
        focus_phase_lut[70][11] = 10'd4;
        focus_phase_lut[70][12] = 10'd0;
        focus_phase_lut[70][13] = 10'd4;
        focus_phase_lut[70][14] = 10'd17;
        focus_phase_lut[70][15] = 10'd21;
        focus_phase_lut[70][16] = 10'd8;
        focus_phase_lut[70][17] = 10'd4;
        focus_phase_lut[70][18] = 10'd8;
        focus_phase_lut[70][19] = 10'd21;
        focus_phase_lut[70][20] = 10'd34;
        focus_phase_lut[70][21] = 10'd21;
        focus_phase_lut[70][22] = 10'd17;
        focus_phase_lut[70][23] = 10'd21;
        focus_phase_lut[70][24] = 10'd34;
        // focus_dist=3650mm (idx=71)
        focus_phase_lut[71][ 0] = 10'd34;
        focus_phase_lut[71][ 1] = 10'd21;
        focus_phase_lut[71][ 2] = 10'd17;
        focus_phase_lut[71][ 3] = 10'd21;
        focus_phase_lut[71][ 4] = 10'd34;
        focus_phase_lut[71][ 5] = 10'd21;
        focus_phase_lut[71][ 6] = 10'd8;
        focus_phase_lut[71][ 7] = 10'd4;
        focus_phase_lut[71][ 8] = 10'd8;
        focus_phase_lut[71][ 9] = 10'd21;
        focus_phase_lut[71][10] = 10'd17;
        focus_phase_lut[71][11] = 10'd4;
        focus_phase_lut[71][12] = 10'd0;
        focus_phase_lut[71][13] = 10'd4;
        focus_phase_lut[71][14] = 10'd17;
        focus_phase_lut[71][15] = 10'd21;
        focus_phase_lut[71][16] = 10'd8;
        focus_phase_lut[71][17] = 10'd4;
        focus_phase_lut[71][18] = 10'd8;
        focus_phase_lut[71][19] = 10'd21;
        focus_phase_lut[71][20] = 10'd34;
        focus_phase_lut[71][21] = 10'd21;
        focus_phase_lut[71][22] = 10'd17;
        focus_phase_lut[71][23] = 10'd21;
        focus_phase_lut[71][24] = 10'd34;
        // focus_dist=3700mm (idx=72)
        focus_phase_lut[72][ 0] = 10'd33;
        focus_phase_lut[72][ 1] = 10'd21;
        focus_phase_lut[72][ 2] = 10'd17;
        focus_phase_lut[72][ 3] = 10'd21;
        focus_phase_lut[72][ 4] = 10'd33;
        focus_phase_lut[72][ 5] = 10'd21;
        focus_phase_lut[72][ 6] = 10'd8;
        focus_phase_lut[72][ 7] = 10'd4;
        focus_phase_lut[72][ 8] = 10'd8;
        focus_phase_lut[72][ 9] = 10'd21;
        focus_phase_lut[72][10] = 10'd17;
        focus_phase_lut[72][11] = 10'd4;
        focus_phase_lut[72][12] = 10'd0;
        focus_phase_lut[72][13] = 10'd4;
        focus_phase_lut[72][14] = 10'd17;
        focus_phase_lut[72][15] = 10'd21;
        focus_phase_lut[72][16] = 10'd8;
        focus_phase_lut[72][17] = 10'd4;
        focus_phase_lut[72][18] = 10'd8;
        focus_phase_lut[72][19] = 10'd21;
        focus_phase_lut[72][20] = 10'd33;
        focus_phase_lut[72][21] = 10'd21;
        focus_phase_lut[72][22] = 10'd17;
        focus_phase_lut[72][23] = 10'd21;
        focus_phase_lut[72][24] = 10'd33;
        // focus_dist=3750mm (idx=73)
        focus_phase_lut[73][ 0] = 10'd33;
        focus_phase_lut[73][ 1] = 10'd20;
        focus_phase_lut[73][ 2] = 10'd16;
        focus_phase_lut[73][ 3] = 10'd20;
        focus_phase_lut[73][ 4] = 10'd33;
        focus_phase_lut[73][ 5] = 10'd20;
        focus_phase_lut[73][ 6] = 10'd8;
        focus_phase_lut[73][ 7] = 10'd4;
        focus_phase_lut[73][ 8] = 10'd8;
        focus_phase_lut[73][ 9] = 10'd20;
        focus_phase_lut[73][10] = 10'd16;
        focus_phase_lut[73][11] = 10'd4;
        focus_phase_lut[73][12] = 10'd0;
        focus_phase_lut[73][13] = 10'd4;
        focus_phase_lut[73][14] = 10'd16;
        focus_phase_lut[73][15] = 10'd20;
        focus_phase_lut[73][16] = 10'd8;
        focus_phase_lut[73][17] = 10'd4;
        focus_phase_lut[73][18] = 10'd8;
        focus_phase_lut[73][19] = 10'd20;
        focus_phase_lut[73][20] = 10'd33;
        focus_phase_lut[73][21] = 10'd20;
        focus_phase_lut[73][22] = 10'd16;
        focus_phase_lut[73][23] = 10'd20;
        focus_phase_lut[73][24] = 10'd33;
        // focus_dist=3800mm (idx=74)
        focus_phase_lut[74][ 0] = 10'd32;
        focus_phase_lut[74][ 1] = 10'd20;
        focus_phase_lut[74][ 2] = 10'd16;
        focus_phase_lut[74][ 3] = 10'd20;
        focus_phase_lut[74][ 4] = 10'd32;
        focus_phase_lut[74][ 5] = 10'd20;
        focus_phase_lut[74][ 6] = 10'd8;
        focus_phase_lut[74][ 7] = 10'd4;
        focus_phase_lut[74][ 8] = 10'd8;
        focus_phase_lut[74][ 9] = 10'd20;
        focus_phase_lut[74][10] = 10'd16;
        focus_phase_lut[74][11] = 10'd4;
        focus_phase_lut[74][12] = 10'd0;
        focus_phase_lut[74][13] = 10'd4;
        focus_phase_lut[74][14] = 10'd16;
        focus_phase_lut[74][15] = 10'd20;
        focus_phase_lut[74][16] = 10'd8;
        focus_phase_lut[74][17] = 10'd4;
        focus_phase_lut[74][18] = 10'd8;
        focus_phase_lut[74][19] = 10'd20;
        focus_phase_lut[74][20] = 10'd32;
        focus_phase_lut[74][21] = 10'd20;
        focus_phase_lut[74][22] = 10'd16;
        focus_phase_lut[74][23] = 10'd20;
        focus_phase_lut[74][24] = 10'd32;
        // focus_dist=3850mm (idx=75)
        focus_phase_lut[75][ 0] = 10'd32;
        focus_phase_lut[75][ 1] = 10'd20;
        focus_phase_lut[75][ 2] = 10'd16;
        focus_phase_lut[75][ 3] = 10'd20;
        focus_phase_lut[75][ 4] = 10'd32;
        focus_phase_lut[75][ 5] = 10'd20;
        focus_phase_lut[75][ 6] = 10'd8;
        focus_phase_lut[75][ 7] = 10'd4;
        focus_phase_lut[75][ 8] = 10'd8;
        focus_phase_lut[75][ 9] = 10'd20;
        focus_phase_lut[75][10] = 10'd16;
        focus_phase_lut[75][11] = 10'd4;
        focus_phase_lut[75][12] = 10'd0;
        focus_phase_lut[75][13] = 10'd4;
        focus_phase_lut[75][14] = 10'd16;
        focus_phase_lut[75][15] = 10'd20;
        focus_phase_lut[75][16] = 10'd8;
        focus_phase_lut[75][17] = 10'd4;
        focus_phase_lut[75][18] = 10'd8;
        focus_phase_lut[75][19] = 10'd20;
        focus_phase_lut[75][20] = 10'd32;
        focus_phase_lut[75][21] = 10'd20;
        focus_phase_lut[75][22] = 10'd16;
        focus_phase_lut[75][23] = 10'd20;
        focus_phase_lut[75][24] = 10'd32;
        // focus_dist=3900mm (idx=76)
        focus_phase_lut[76][ 0] = 10'd31;
        focus_phase_lut[76][ 1] = 10'd20;
        focus_phase_lut[76][ 2] = 10'd16;
        focus_phase_lut[76][ 3] = 10'd20;
        focus_phase_lut[76][ 4] = 10'd31;
        focus_phase_lut[76][ 5] = 10'd20;
        focus_phase_lut[76][ 6] = 10'd8;
        focus_phase_lut[76][ 7] = 10'd4;
        focus_phase_lut[76][ 8] = 10'd8;
        focus_phase_lut[76][ 9] = 10'd20;
        focus_phase_lut[76][10] = 10'd16;
        focus_phase_lut[76][11] = 10'd4;
        focus_phase_lut[76][12] = 10'd0;
        focus_phase_lut[76][13] = 10'd4;
        focus_phase_lut[76][14] = 10'd16;
        focus_phase_lut[76][15] = 10'd20;
        focus_phase_lut[76][16] = 10'd8;
        focus_phase_lut[76][17] = 10'd4;
        focus_phase_lut[76][18] = 10'd8;
        focus_phase_lut[76][19] = 10'd20;
        focus_phase_lut[76][20] = 10'd31;
        focus_phase_lut[76][21] = 10'd20;
        focus_phase_lut[76][22] = 10'd16;
        focus_phase_lut[76][23] = 10'd20;
        focus_phase_lut[76][24] = 10'd31;
        // focus_dist=3950mm (idx=77)
        focus_phase_lut[77][ 0] = 10'd31;
        focus_phase_lut[77][ 1] = 10'd19;
        focus_phase_lut[77][ 2] = 10'd15;
        focus_phase_lut[77][ 3] = 10'd19;
        focus_phase_lut[77][ 4] = 10'd31;
        focus_phase_lut[77][ 5] = 10'd19;
        focus_phase_lut[77][ 6] = 10'd8;
        focus_phase_lut[77][ 7] = 10'd4;
        focus_phase_lut[77][ 8] = 10'd8;
        focus_phase_lut[77][ 9] = 10'd19;
        focus_phase_lut[77][10] = 10'd15;
        focus_phase_lut[77][11] = 10'd4;
        focus_phase_lut[77][12] = 10'd0;
        focus_phase_lut[77][13] = 10'd4;
        focus_phase_lut[77][14] = 10'd15;
        focus_phase_lut[77][15] = 10'd19;
        focus_phase_lut[77][16] = 10'd8;
        focus_phase_lut[77][17] = 10'd4;
        focus_phase_lut[77][18] = 10'd8;
        focus_phase_lut[77][19] = 10'd19;
        focus_phase_lut[77][20] = 10'd31;
        focus_phase_lut[77][21] = 10'd19;
        focus_phase_lut[77][22] = 10'd15;
        focus_phase_lut[77][23] = 10'd19;
        focus_phase_lut[77][24] = 10'd31;
        // focus_dist=4000mm (idx=78)
        focus_phase_lut[78][ 0] = 10'd31;
        focus_phase_lut[78][ 1] = 10'd19;
        focus_phase_lut[78][ 2] = 10'd15;
        focus_phase_lut[78][ 3] = 10'd19;
        focus_phase_lut[78][ 4] = 10'd31;
        focus_phase_lut[78][ 5] = 10'd19;
        focus_phase_lut[78][ 6] = 10'd8;
        focus_phase_lut[78][ 7] = 10'd4;
        focus_phase_lut[78][ 8] = 10'd8;
        focus_phase_lut[78][ 9] = 10'd19;
        focus_phase_lut[78][10] = 10'd15;
        focus_phase_lut[78][11] = 10'd4;
        focus_phase_lut[78][12] = 10'd0;
        focus_phase_lut[78][13] = 10'd4;
        focus_phase_lut[78][14] = 10'd15;
        focus_phase_lut[78][15] = 10'd19;
        focus_phase_lut[78][16] = 10'd8;
        focus_phase_lut[78][17] = 10'd4;
        focus_phase_lut[78][18] = 10'd8;
        focus_phase_lut[78][19] = 10'd19;
        focus_phase_lut[78][20] = 10'd31;
        focus_phase_lut[78][21] = 10'd19;
        focus_phase_lut[78][22] = 10'd15;
        focus_phase_lut[78][23] = 10'd19;
        focus_phase_lut[78][24] = 10'd31;
        // focus_dist=4050mm (idx=79)
        focus_phase_lut[79][ 0] = 10'd30;
        focus_phase_lut[79][ 1] = 10'd19;
        focus_phase_lut[79][ 2] = 10'd15;
        focus_phase_lut[79][ 3] = 10'd19;
        focus_phase_lut[79][ 4] = 10'd30;
        focus_phase_lut[79][ 5] = 10'd19;
        focus_phase_lut[79][ 6] = 10'd8;
        focus_phase_lut[79][ 7] = 10'd4;
        focus_phase_lut[79][ 8] = 10'd8;
        focus_phase_lut[79][ 9] = 10'd19;
        focus_phase_lut[79][10] = 10'd15;
        focus_phase_lut[79][11] = 10'd4;
        focus_phase_lut[79][12] = 10'd0;
        focus_phase_lut[79][13] = 10'd4;
        focus_phase_lut[79][14] = 10'd15;
        focus_phase_lut[79][15] = 10'd19;
        focus_phase_lut[79][16] = 10'd8;
        focus_phase_lut[79][17] = 10'd4;
        focus_phase_lut[79][18] = 10'd8;
        focus_phase_lut[79][19] = 10'd19;
        focus_phase_lut[79][20] = 10'd30;
        focus_phase_lut[79][21] = 10'd19;
        focus_phase_lut[79][22] = 10'd15;
        focus_phase_lut[79][23] = 10'd19;
        focus_phase_lut[79][24] = 10'd30;
        // focus_dist=4100mm (idx=80)
        focus_phase_lut[80][ 0] = 10'd30;
        focus_phase_lut[80][ 1] = 10'd19;
        focus_phase_lut[80][ 2] = 10'd15;
        focus_phase_lut[80][ 3] = 10'd19;
        focus_phase_lut[80][ 4] = 10'd30;
        focus_phase_lut[80][ 5] = 10'd19;
        focus_phase_lut[80][ 6] = 10'd7;
        focus_phase_lut[80][ 7] = 10'd4;
        focus_phase_lut[80][ 8] = 10'd7;
        focus_phase_lut[80][ 9] = 10'd19;
        focus_phase_lut[80][10] = 10'd15;
        focus_phase_lut[80][11] = 10'd4;
        focus_phase_lut[80][12] = 10'd0;
        focus_phase_lut[80][13] = 10'd4;
        focus_phase_lut[80][14] = 10'd15;
        focus_phase_lut[80][15] = 10'd19;
        focus_phase_lut[80][16] = 10'd7;
        focus_phase_lut[80][17] = 10'd4;
        focus_phase_lut[80][18] = 10'd7;
        focus_phase_lut[80][19] = 10'd19;
        focus_phase_lut[80][20] = 10'd30;
        focus_phase_lut[80][21] = 10'd19;
        focus_phase_lut[80][22] = 10'd15;
        focus_phase_lut[80][23] = 10'd19;
        focus_phase_lut[80][24] = 10'd30;
        // focus_dist=4150mm (idx=81)
        focus_phase_lut[81][ 0] = 10'd29;
        focus_phase_lut[81][ 1] = 10'd18;
        focus_phase_lut[81][ 2] = 10'd15;
        focus_phase_lut[81][ 3] = 10'd18;
        focus_phase_lut[81][ 4] = 10'd29;
        focus_phase_lut[81][ 5] = 10'd18;
        focus_phase_lut[81][ 6] = 10'd7;
        focus_phase_lut[81][ 7] = 10'd4;
        focus_phase_lut[81][ 8] = 10'd7;
        focus_phase_lut[81][ 9] = 10'd18;
        focus_phase_lut[81][10] = 10'd15;
        focus_phase_lut[81][11] = 10'd4;
        focus_phase_lut[81][12] = 10'd0;
        focus_phase_lut[81][13] = 10'd4;
        focus_phase_lut[81][14] = 10'd15;
        focus_phase_lut[81][15] = 10'd18;
        focus_phase_lut[81][16] = 10'd7;
        focus_phase_lut[81][17] = 10'd4;
        focus_phase_lut[81][18] = 10'd7;
        focus_phase_lut[81][19] = 10'd18;
        focus_phase_lut[81][20] = 10'd29;
        focus_phase_lut[81][21] = 10'd18;
        focus_phase_lut[81][22] = 10'd15;
        focus_phase_lut[81][23] = 10'd18;
        focus_phase_lut[81][24] = 10'd29;
        // focus_dist=4200mm (idx=82)
        focus_phase_lut[82][ 0] = 10'd29;
        focus_phase_lut[82][ 1] = 10'd18;
        focus_phase_lut[82][ 2] = 10'd15;
        focus_phase_lut[82][ 3] = 10'd18;
        focus_phase_lut[82][ 4] = 10'd29;
        focus_phase_lut[82][ 5] = 10'd18;
        focus_phase_lut[82][ 6] = 10'd7;
        focus_phase_lut[82][ 7] = 10'd4;
        focus_phase_lut[82][ 8] = 10'd7;
        focus_phase_lut[82][ 9] = 10'd18;
        focus_phase_lut[82][10] = 10'd15;
        focus_phase_lut[82][11] = 10'd4;
        focus_phase_lut[82][12] = 10'd0;
        focus_phase_lut[82][13] = 10'd4;
        focus_phase_lut[82][14] = 10'd15;
        focus_phase_lut[82][15] = 10'd18;
        focus_phase_lut[82][16] = 10'd7;
        focus_phase_lut[82][17] = 10'd4;
        focus_phase_lut[82][18] = 10'd7;
        focus_phase_lut[82][19] = 10'd18;
        focus_phase_lut[82][20] = 10'd29;
        focus_phase_lut[82][21] = 10'd18;
        focus_phase_lut[82][22] = 10'd15;
        focus_phase_lut[82][23] = 10'd18;
        focus_phase_lut[82][24] = 10'd29;
        // focus_dist=4250mm (idx=83)
        focus_phase_lut[83][ 0] = 10'd29;
        focus_phase_lut[83][ 1] = 10'd18;
        focus_phase_lut[83][ 2] = 10'd14;
        focus_phase_lut[83][ 3] = 10'd18;
        focus_phase_lut[83][ 4] = 10'd29;
        focus_phase_lut[83][ 5] = 10'd18;
        focus_phase_lut[83][ 6] = 10'd7;
        focus_phase_lut[83][ 7] = 10'd4;
        focus_phase_lut[83][ 8] = 10'd7;
        focus_phase_lut[83][ 9] = 10'd18;
        focus_phase_lut[83][10] = 10'd14;
        focus_phase_lut[83][11] = 10'd4;
        focus_phase_lut[83][12] = 10'd0;
        focus_phase_lut[83][13] = 10'd4;
        focus_phase_lut[83][14] = 10'd14;
        focus_phase_lut[83][15] = 10'd18;
        focus_phase_lut[83][16] = 10'd7;
        focus_phase_lut[83][17] = 10'd4;
        focus_phase_lut[83][18] = 10'd7;
        focus_phase_lut[83][19] = 10'd18;
        focus_phase_lut[83][20] = 10'd29;
        focus_phase_lut[83][21] = 10'd18;
        focus_phase_lut[83][22] = 10'd14;
        focus_phase_lut[83][23] = 10'd18;
        focus_phase_lut[83][24] = 10'd29;
        // focus_dist=4300mm (idx=84)
        focus_phase_lut[84][ 0] = 10'd28;
        focus_phase_lut[84][ 1] = 10'd18;
        focus_phase_lut[84][ 2] = 10'd14;
        focus_phase_lut[84][ 3] = 10'd18;
        focus_phase_lut[84][ 4] = 10'd28;
        focus_phase_lut[84][ 5] = 10'd18;
        focus_phase_lut[84][ 6] = 10'd7;
        focus_phase_lut[84][ 7] = 10'd4;
        focus_phase_lut[84][ 8] = 10'd7;
        focus_phase_lut[84][ 9] = 10'd18;
        focus_phase_lut[84][10] = 10'd14;
        focus_phase_lut[84][11] = 10'd4;
        focus_phase_lut[84][12] = 10'd0;
        focus_phase_lut[84][13] = 10'd4;
        focus_phase_lut[84][14] = 10'd14;
        focus_phase_lut[84][15] = 10'd18;
        focus_phase_lut[84][16] = 10'd7;
        focus_phase_lut[84][17] = 10'd4;
        focus_phase_lut[84][18] = 10'd7;
        focus_phase_lut[84][19] = 10'd18;
        focus_phase_lut[84][20] = 10'd28;
        focus_phase_lut[84][21] = 10'd18;
        focus_phase_lut[84][22] = 10'd14;
        focus_phase_lut[84][23] = 10'd18;
        focus_phase_lut[84][24] = 10'd28;
        // focus_dist=4350mm (idx=85)
        focus_phase_lut[85][ 0] = 10'd28;
        focus_phase_lut[85][ 1] = 10'd18;
        focus_phase_lut[85][ 2] = 10'd14;
        focus_phase_lut[85][ 3] = 10'd18;
        focus_phase_lut[85][ 4] = 10'd28;
        focus_phase_lut[85][ 5] = 10'd18;
        focus_phase_lut[85][ 6] = 10'd7;
        focus_phase_lut[85][ 7] = 10'd4;
        focus_phase_lut[85][ 8] = 10'd7;
        focus_phase_lut[85][ 9] = 10'd18;
        focus_phase_lut[85][10] = 10'd14;
        focus_phase_lut[85][11] = 10'd4;
        focus_phase_lut[85][12] = 10'd0;
        focus_phase_lut[85][13] = 10'd4;
        focus_phase_lut[85][14] = 10'd14;
        focus_phase_lut[85][15] = 10'd18;
        focus_phase_lut[85][16] = 10'd7;
        focus_phase_lut[85][17] = 10'd4;
        focus_phase_lut[85][18] = 10'd7;
        focus_phase_lut[85][19] = 10'd18;
        focus_phase_lut[85][20] = 10'd28;
        focus_phase_lut[85][21] = 10'd18;
        focus_phase_lut[85][22] = 10'd14;
        focus_phase_lut[85][23] = 10'd18;
        focus_phase_lut[85][24] = 10'd28;
        // focus_dist=4400mm (idx=86)
        focus_phase_lut[86][ 0] = 10'd28;
        focus_phase_lut[86][ 1] = 10'd17;
        focus_phase_lut[86][ 2] = 10'd14;
        focus_phase_lut[86][ 3] = 10'd17;
        focus_phase_lut[86][ 4] = 10'd28;
        focus_phase_lut[86][ 5] = 10'd17;
        focus_phase_lut[86][ 6] = 10'd7;
        focus_phase_lut[86][ 7] = 10'd3;
        focus_phase_lut[86][ 8] = 10'd7;
        focus_phase_lut[86][ 9] = 10'd17;
        focus_phase_lut[86][10] = 10'd14;
        focus_phase_lut[86][11] = 10'd3;
        focus_phase_lut[86][12] = 10'd0;
        focus_phase_lut[86][13] = 10'd3;
        focus_phase_lut[86][14] = 10'd14;
        focus_phase_lut[86][15] = 10'd17;
        focus_phase_lut[86][16] = 10'd7;
        focus_phase_lut[86][17] = 10'd3;
        focus_phase_lut[86][18] = 10'd7;
        focus_phase_lut[86][19] = 10'd17;
        focus_phase_lut[86][20] = 10'd28;
        focus_phase_lut[86][21] = 10'd17;
        focus_phase_lut[86][22] = 10'd14;
        focus_phase_lut[86][23] = 10'd17;
        focus_phase_lut[86][24] = 10'd28;
        // focus_dist=4450mm (idx=87)
        focus_phase_lut[87][ 0] = 10'd27;
        focus_phase_lut[87][ 1] = 10'd17;
        focus_phase_lut[87][ 2] = 10'd14;
        focus_phase_lut[87][ 3] = 10'd17;
        focus_phase_lut[87][ 4] = 10'd27;
        focus_phase_lut[87][ 5] = 10'd17;
        focus_phase_lut[87][ 6] = 10'd7;
        focus_phase_lut[87][ 7] = 10'd3;
        focus_phase_lut[87][ 8] = 10'd7;
        focus_phase_lut[87][ 9] = 10'd17;
        focus_phase_lut[87][10] = 10'd14;
        focus_phase_lut[87][11] = 10'd3;
        focus_phase_lut[87][12] = 10'd0;
        focus_phase_lut[87][13] = 10'd3;
        focus_phase_lut[87][14] = 10'd14;
        focus_phase_lut[87][15] = 10'd17;
        focus_phase_lut[87][16] = 10'd7;
        focus_phase_lut[87][17] = 10'd3;
        focus_phase_lut[87][18] = 10'd7;
        focus_phase_lut[87][19] = 10'd17;
        focus_phase_lut[87][20] = 10'd27;
        focus_phase_lut[87][21] = 10'd17;
        focus_phase_lut[87][22] = 10'd14;
        focus_phase_lut[87][23] = 10'd17;
        focus_phase_lut[87][24] = 10'd27;
        // focus_dist=4500mm (idx=88)
        focus_phase_lut[88][ 0] = 10'd27;
        focus_phase_lut[88][ 1] = 10'd17;
        focus_phase_lut[88][ 2] = 10'd14;
        focus_phase_lut[88][ 3] = 10'd17;
        focus_phase_lut[88][ 4] = 10'd27;
        focus_phase_lut[88][ 5] = 10'd17;
        focus_phase_lut[88][ 6] = 10'd7;
        focus_phase_lut[88][ 7] = 10'd3;
        focus_phase_lut[88][ 8] = 10'd7;
        focus_phase_lut[88][ 9] = 10'd17;
        focus_phase_lut[88][10] = 10'd14;
        focus_phase_lut[88][11] = 10'd3;
        focus_phase_lut[88][12] = 10'd0;
        focus_phase_lut[88][13] = 10'd3;
        focus_phase_lut[88][14] = 10'd14;
        focus_phase_lut[88][15] = 10'd17;
        focus_phase_lut[88][16] = 10'd7;
        focus_phase_lut[88][17] = 10'd3;
        focus_phase_lut[88][18] = 10'd7;
        focus_phase_lut[88][19] = 10'd17;
        focus_phase_lut[88][20] = 10'd27;
        focus_phase_lut[88][21] = 10'd17;
        focus_phase_lut[88][22] = 10'd14;
        focus_phase_lut[88][23] = 10'd17;
        focus_phase_lut[88][24] = 10'd27;
        // focus_dist=4550mm (idx=89)
        focus_phase_lut[89][ 0] = 10'd27;
        focus_phase_lut[89][ 1] = 10'd17;
        focus_phase_lut[89][ 2] = 10'd13;
        focus_phase_lut[89][ 3] = 10'd17;
        focus_phase_lut[89][ 4] = 10'd27;
        focus_phase_lut[89][ 5] = 10'd17;
        focus_phase_lut[89][ 6] = 10'd7;
        focus_phase_lut[89][ 7] = 10'd3;
        focus_phase_lut[89][ 8] = 10'd7;
        focus_phase_lut[89][ 9] = 10'd17;
        focus_phase_lut[89][10] = 10'd13;
        focus_phase_lut[89][11] = 10'd3;
        focus_phase_lut[89][12] = 10'd0;
        focus_phase_lut[89][13] = 10'd3;
        focus_phase_lut[89][14] = 10'd13;
        focus_phase_lut[89][15] = 10'd17;
        focus_phase_lut[89][16] = 10'd7;
        focus_phase_lut[89][17] = 10'd3;
        focus_phase_lut[89][18] = 10'd7;
        focus_phase_lut[89][19] = 10'd17;
        focus_phase_lut[89][20] = 10'd27;
        focus_phase_lut[89][21] = 10'd17;
        focus_phase_lut[89][22] = 10'd13;
        focus_phase_lut[89][23] = 10'd17;
        focus_phase_lut[89][24] = 10'd27;
        // focus_dist=4600mm (idx=90)
        focus_phase_lut[90][ 0] = 10'd27;
        focus_phase_lut[90][ 1] = 10'd17;
        focus_phase_lut[90][ 2] = 10'd13;
        focus_phase_lut[90][ 3] = 10'd17;
        focus_phase_lut[90][ 4] = 10'd27;
        focus_phase_lut[90][ 5] = 10'd17;
        focus_phase_lut[90][ 6] = 10'd7;
        focus_phase_lut[90][ 7] = 10'd3;
        focus_phase_lut[90][ 8] = 10'd7;
        focus_phase_lut[90][ 9] = 10'd17;
        focus_phase_lut[90][10] = 10'd13;
        focus_phase_lut[90][11] = 10'd3;
        focus_phase_lut[90][12] = 10'd0;
        focus_phase_lut[90][13] = 10'd3;
        focus_phase_lut[90][14] = 10'd13;
        focus_phase_lut[90][15] = 10'd17;
        focus_phase_lut[90][16] = 10'd7;
        focus_phase_lut[90][17] = 10'd3;
        focus_phase_lut[90][18] = 10'd7;
        focus_phase_lut[90][19] = 10'd17;
        focus_phase_lut[90][20] = 10'd27;
        focus_phase_lut[90][21] = 10'd17;
        focus_phase_lut[90][22] = 10'd13;
        focus_phase_lut[90][23] = 10'd17;
        focus_phase_lut[90][24] = 10'd27;
        // focus_dist=4650mm (idx=91)
        focus_phase_lut[91][ 0] = 10'd26;
        focus_phase_lut[91][ 1] = 10'd16;
        focus_phase_lut[91][ 2] = 10'd13;
        focus_phase_lut[91][ 3] = 10'd16;
        focus_phase_lut[91][ 4] = 10'd26;
        focus_phase_lut[91][ 5] = 10'd16;
        focus_phase_lut[91][ 6] = 10'd7;
        focus_phase_lut[91][ 7] = 10'd3;
        focus_phase_lut[91][ 8] = 10'd7;
        focus_phase_lut[91][ 9] = 10'd16;
        focus_phase_lut[91][10] = 10'd13;
        focus_phase_lut[91][11] = 10'd3;
        focus_phase_lut[91][12] = 10'd0;
        focus_phase_lut[91][13] = 10'd3;
        focus_phase_lut[91][14] = 10'd13;
        focus_phase_lut[91][15] = 10'd16;
        focus_phase_lut[91][16] = 10'd7;
        focus_phase_lut[91][17] = 10'd3;
        focus_phase_lut[91][18] = 10'd7;
        focus_phase_lut[91][19] = 10'd16;
        focus_phase_lut[91][20] = 10'd26;
        focus_phase_lut[91][21] = 10'd16;
        focus_phase_lut[91][22] = 10'd13;
        focus_phase_lut[91][23] = 10'd16;
        focus_phase_lut[91][24] = 10'd26;
        // focus_dist=4700mm (idx=92)
        focus_phase_lut[92][ 0] = 10'd26;
        focus_phase_lut[92][ 1] = 10'd16;
        focus_phase_lut[92][ 2] = 10'd13;
        focus_phase_lut[92][ 3] = 10'd16;
        focus_phase_lut[92][ 4] = 10'd26;
        focus_phase_lut[92][ 5] = 10'd16;
        focus_phase_lut[92][ 6] = 10'd7;
        focus_phase_lut[92][ 7] = 10'd3;
        focus_phase_lut[92][ 8] = 10'd7;
        focus_phase_lut[92][ 9] = 10'd16;
        focus_phase_lut[92][10] = 10'd13;
        focus_phase_lut[92][11] = 10'd3;
        focus_phase_lut[92][12] = 10'd0;
        focus_phase_lut[92][13] = 10'd3;
        focus_phase_lut[92][14] = 10'd13;
        focus_phase_lut[92][15] = 10'd16;
        focus_phase_lut[92][16] = 10'd7;
        focus_phase_lut[92][17] = 10'd3;
        focus_phase_lut[92][18] = 10'd7;
        focus_phase_lut[92][19] = 10'd16;
        focus_phase_lut[92][20] = 10'd26;
        focus_phase_lut[92][21] = 10'd16;
        focus_phase_lut[92][22] = 10'd13;
        focus_phase_lut[92][23] = 10'd16;
        focus_phase_lut[92][24] = 10'd26;
        // focus_dist=4750mm (idx=93)
        focus_phase_lut[93][ 0] = 10'd26;
        focus_phase_lut[93][ 1] = 10'd16;
        focus_phase_lut[93][ 2] = 10'd13;
        focus_phase_lut[93][ 3] = 10'd16;
        focus_phase_lut[93][ 4] = 10'd26;
        focus_phase_lut[93][ 5] = 10'd16;
        focus_phase_lut[93][ 6] = 10'd6;
        focus_phase_lut[93][ 7] = 10'd3;
        focus_phase_lut[93][ 8] = 10'd6;
        focus_phase_lut[93][ 9] = 10'd16;
        focus_phase_lut[93][10] = 10'd13;
        focus_phase_lut[93][11] = 10'd3;
        focus_phase_lut[93][12] = 10'd0;
        focus_phase_lut[93][13] = 10'd3;
        focus_phase_lut[93][14] = 10'd13;
        focus_phase_lut[93][15] = 10'd16;
        focus_phase_lut[93][16] = 10'd6;
        focus_phase_lut[93][17] = 10'd3;
        focus_phase_lut[93][18] = 10'd6;
        focus_phase_lut[93][19] = 10'd16;
        focus_phase_lut[93][20] = 10'd26;
        focus_phase_lut[93][21] = 10'd16;
        focus_phase_lut[93][22] = 10'd13;
        focus_phase_lut[93][23] = 10'd16;
        focus_phase_lut[93][24] = 10'd26;
        // focus_dist=4800mm (idx=94)
        focus_phase_lut[94][ 0] = 10'd25;
        focus_phase_lut[94][ 1] = 10'd16;
        focus_phase_lut[94][ 2] = 10'd13;
        focus_phase_lut[94][ 3] = 10'd16;
        focus_phase_lut[94][ 4] = 10'd25;
        focus_phase_lut[94][ 5] = 10'd16;
        focus_phase_lut[94][ 6] = 10'd6;
        focus_phase_lut[94][ 7] = 10'd3;
        focus_phase_lut[94][ 8] = 10'd6;
        focus_phase_lut[94][ 9] = 10'd16;
        focus_phase_lut[94][10] = 10'd13;
        focus_phase_lut[94][11] = 10'd3;
        focus_phase_lut[94][12] = 10'd0;
        focus_phase_lut[94][13] = 10'd3;
        focus_phase_lut[94][14] = 10'd13;
        focus_phase_lut[94][15] = 10'd16;
        focus_phase_lut[94][16] = 10'd6;
        focus_phase_lut[94][17] = 10'd3;
        focus_phase_lut[94][18] = 10'd6;
        focus_phase_lut[94][19] = 10'd16;
        focus_phase_lut[94][20] = 10'd25;
        focus_phase_lut[94][21] = 10'd16;
        focus_phase_lut[94][22] = 10'd13;
        focus_phase_lut[94][23] = 10'd16;
        focus_phase_lut[94][24] = 10'd25;
        // focus_dist=4850mm (idx=95)
        focus_phase_lut[95][ 0] = 10'd25;
        focus_phase_lut[95][ 1] = 10'd16;
        focus_phase_lut[95][ 2] = 10'd13;
        focus_phase_lut[95][ 3] = 10'd16;
        focus_phase_lut[95][ 4] = 10'd25;
        focus_phase_lut[95][ 5] = 10'd16;
        focus_phase_lut[95][ 6] = 10'd6;
        focus_phase_lut[95][ 7] = 10'd3;
        focus_phase_lut[95][ 8] = 10'd6;
        focus_phase_lut[95][ 9] = 10'd16;
        focus_phase_lut[95][10] = 10'd13;
        focus_phase_lut[95][11] = 10'd3;
        focus_phase_lut[95][12] = 10'd0;
        focus_phase_lut[95][13] = 10'd3;
        focus_phase_lut[95][14] = 10'd13;
        focus_phase_lut[95][15] = 10'd16;
        focus_phase_lut[95][16] = 10'd6;
        focus_phase_lut[95][17] = 10'd3;
        focus_phase_lut[95][18] = 10'd6;
        focus_phase_lut[95][19] = 10'd16;
        focus_phase_lut[95][20] = 10'd25;
        focus_phase_lut[95][21] = 10'd16;
        focus_phase_lut[95][22] = 10'd13;
        focus_phase_lut[95][23] = 10'd16;
        focus_phase_lut[95][24] = 10'd25;
        // focus_dist=4900mm (idx=96)
        focus_phase_lut[96][ 0] = 10'd25;
        focus_phase_lut[96][ 1] = 10'd16;
        focus_phase_lut[96][ 2] = 10'd12;
        focus_phase_lut[96][ 3] = 10'd16;
        focus_phase_lut[96][ 4] = 10'd25;
        focus_phase_lut[96][ 5] = 10'd16;
        focus_phase_lut[96][ 6] = 10'd6;
        focus_phase_lut[96][ 7] = 10'd3;
        focus_phase_lut[96][ 8] = 10'd6;
        focus_phase_lut[96][ 9] = 10'd16;
        focus_phase_lut[96][10] = 10'd12;
        focus_phase_lut[96][11] = 10'd3;
        focus_phase_lut[96][12] = 10'd0;
        focus_phase_lut[96][13] = 10'd3;
        focus_phase_lut[96][14] = 10'd12;
        focus_phase_lut[96][15] = 10'd16;
        focus_phase_lut[96][16] = 10'd6;
        focus_phase_lut[96][17] = 10'd3;
        focus_phase_lut[96][18] = 10'd6;
        focus_phase_lut[96][19] = 10'd16;
        focus_phase_lut[96][20] = 10'd25;
        focus_phase_lut[96][21] = 10'd16;
        focus_phase_lut[96][22] = 10'd12;
        focus_phase_lut[96][23] = 10'd16;
        focus_phase_lut[96][24] = 10'd25;
        // focus_dist=4950mm (idx=97)
        focus_phase_lut[97][ 0] = 10'd25;
        focus_phase_lut[97][ 1] = 10'd15;
        focus_phase_lut[97][ 2] = 10'd12;
        focus_phase_lut[97][ 3] = 10'd15;
        focus_phase_lut[97][ 4] = 10'd25;
        focus_phase_lut[97][ 5] = 10'd15;
        focus_phase_lut[97][ 6] = 10'd6;
        focus_phase_lut[97][ 7] = 10'd3;
        focus_phase_lut[97][ 8] = 10'd6;
        focus_phase_lut[97][ 9] = 10'd15;
        focus_phase_lut[97][10] = 10'd12;
        focus_phase_lut[97][11] = 10'd3;
        focus_phase_lut[97][12] = 10'd0;
        focus_phase_lut[97][13] = 10'd3;
        focus_phase_lut[97][14] = 10'd12;
        focus_phase_lut[97][15] = 10'd15;
        focus_phase_lut[97][16] = 10'd6;
        focus_phase_lut[97][17] = 10'd3;
        focus_phase_lut[97][18] = 10'd6;
        focus_phase_lut[97][19] = 10'd15;
        focus_phase_lut[97][20] = 10'd25;
        focus_phase_lut[97][21] = 10'd15;
        focus_phase_lut[97][22] = 10'd12;
        focus_phase_lut[97][23] = 10'd15;
        focus_phase_lut[97][24] = 10'd25;
        // focus_dist=5000mm (idx=98)
        focus_phase_lut[98][ 0] = 10'd24;
        focus_phase_lut[98][ 1] = 10'd15;
        focus_phase_lut[98][ 2] = 10'd12;
        focus_phase_lut[98][ 3] = 10'd15;
        focus_phase_lut[98][ 4] = 10'd24;
        focus_phase_lut[98][ 5] = 10'd15;
        focus_phase_lut[98][ 6] = 10'd6;
        focus_phase_lut[98][ 7] = 10'd3;
        focus_phase_lut[98][ 8] = 10'd6;
        focus_phase_lut[98][ 9] = 10'd15;
        focus_phase_lut[98][10] = 10'd12;
        focus_phase_lut[98][11] = 10'd3;
        focus_phase_lut[98][12] = 10'd0;
        focus_phase_lut[98][13] = 10'd3;
        focus_phase_lut[98][14] = 10'd12;
        focus_phase_lut[98][15] = 10'd15;
        focus_phase_lut[98][16] = 10'd6;
        focus_phase_lut[98][17] = 10'd3;
        focus_phase_lut[98][18] = 10'd6;
        focus_phase_lut[98][19] = 10'd15;
        focus_phase_lut[98][20] = 10'd24;
        focus_phase_lut[98][21] = 10'd15;
        focus_phase_lut[98][22] = 10'd12;
        focus_phase_lut[98][23] = 10'd15;
        focus_phase_lut[98][24] = 10'd24;
    end

    //=========================================================================
    // 组合逻辑：LUT索引计算和窗函数选择
    //=========================================================================

    // 偏转LUT索引： steer_angle(-30~+30) + 30 -> 0~60
    // 使用有符号运算确保负数角度正确
    wire signed [7:0] steer_idx_signed;
    wire [6:0] steer_lut_idx;
    assign steer_idx_signed = steer_reg + 8'sd30;
    assign steer_lut_idx = {1'b0, steer_idx_signed[6:0]};
    wire signed [9:0] steer_lut_val;
    assign steer_lut_val = steer_phase_lut[steer_lut_idx[5:0]];

    // 聚焦LUT索引：(focus_dist - 100) / 50 -> 0~98
    wire [12:0] focus_dist_minus100;
    wire [6:0] focus_lut_idx;
    assign focus_dist_minus100 = focus_dist_reg - 13'd100;
    assign focus_lut_idx = focus_dist_minus100 / 7'd50;
    wire [9:0] focus_lut_val;
    assign focus_lut_val = focus_phase_lut[focus_lut_idx[6:0]][calc_idx[4:0]];

    // 窗函数选择
    wire [7:0] window_val;
    assign window_val = (window_reg == 3'd1) ? window_hann[calc_idx[4:0]] :
                        (window_reg == 3'd2) ? window_hamm[calc_idx[4:0]] :
                        (window_reg == 3'd3) ? window_blk[calc_idx[4:0]] :
                                               window_rect[calc_idx[4:0]];

    // 下一计算索引
    wire [4:0] calc_idx_plus1;
    assign calc_idx_plus1 = calc_idx + 5'd1;

    //=========================================================================
    // 主状态机 (时序逻辑)
    //=========================================================================
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state            <= IDLE;
            mode_reg         <= 3'd0;
            steer_reg        <= 8'sd0;
            focus_dist_reg   <= 13'd500;
            window_reg       <= 3'd0;
            calc_idx         <= 5'd0;
            steer_phase_incr <= 10'sd0;
            param_valid      <= 1'b0;
        end else begin
            case (state)
                //---------------------------------------------------------------
                // IDLE: 等待参数更新
                //---------------------------------------------------------------
                IDLE: begin
                    param_valid <= 1'b0;
                    if (param_update) begin
                        mode_reg       <= mode;
                        steer_reg      <= steer_angle;
                        focus_dist_reg <= focus_dist_mm;
                        window_reg     <= window_type;
                        calc_idx       <= 5'd0;
                        state          <= LOAD;
                    end
                end

                //---------------------------------------------------------------
                // LOAD: 加载LUT值，准备计算
                //---------------------------------------------------------------
                LOAD: begin
                    case (mode_reg)
                        3'b000: begin  // 同相
                            steer_phase_incr <= 10'sd0;
                            state <= CALC;
                        end
                        3'b001: begin  // 偏转
                            if (steer_reg >= -8'sd30 && steer_reg <= 8'sd30)
                                steer_phase_incr <= steer_lut_val;
                            else
                                steer_phase_incr <= 10'sd0;
                            state <= CALC;
                        end
                        3'b010: begin  // 聚焦
                            steer_phase_incr <= 10'sd0;
                            state <= CALC;
                        end
                        3'b011: begin  // 偏转+聚焦
                            if (steer_reg >= -8'sd30 && steer_reg <= 8'sd30)
                                steer_phase_incr <= steer_lut_val;
                            else
                                steer_phase_incr <= 10'sd0;
                            state <= CALC;
                        end
                        3'b100: begin  // 同相+加窗
                            steer_phase_incr <= 10'sd0;
                            state <= CALC;
                        end
                        default: begin
                            steer_phase_incr <= 10'sd0;
                            state <= CALC;
                        end
                    endcase
                end

                //---------------------------------------------------------------
                // CALC: 逐路计算25个阵元 (每个时钟周期计算1路)
                //---------------------------------------------------------------
                CALC: begin
                    case (mode_reg)
                        // 同相模式：全0相位，全幅度
                        3'b000: begin
                            phase[calc_idx]     <= 10'd0;
                            amplitude[calc_idx] <= 8'd255;
                        end

                        // 波束偏转：行方向相位梯度
                        3'b001: begin
                            begin : steer_calc_block
                                reg signed [11:0] phase_signed;
                                reg signed [11:0] phase_mod;
                                phase_signed = element_col[calc_idx] * steer_phase_incr;
                                // 对1024取模，处理正负
                                if (phase_signed >= 12'sd0)
                                    phase_mod = phase_signed % 12'sd1024;
                                else begin
                                    phase_mod = 12'sd1024 + (phase_signed % 12'sd1024);
                                    if (phase_mod == 12'sd1024)
                                        phase_mod = 12'sd0;
                                end
                                phase[calc_idx]     <= phase_mod[9:0];
                                amplitude[calc_idx] <= 8'd255;
                            end
                        end

                        // 波束聚焦：从LUT读取相位
                        3'b010: begin
                            phase[calc_idx]     <= focus_lut_val;
                            amplitude[calc_idx] <= 8'd255;
                        end

                        // 偏转+聚焦：相位叠加
                        3'b011: begin
                            begin : sf_calc_block
                                reg signed [11:0] steer_ph;
                                reg signed [11:0] phase_total;
                                reg signed [11:0] phase_mod;
                                steer_ph   = element_col[calc_idx] * steer_phase_incr;
                                phase_total = steer_ph + $signed({1'b0, focus_lut_val});
                                if (phase_total >= 12'sd0)
                                    phase_mod = phase_total % 12'sd1024;
                                else begin
                                    phase_mod = 12'sd1024 + (phase_total % 12'sd1024);
                                    if (phase_mod == 12'sd1024)
                                        phase_mod = 12'sd0;
                                end
                                phase[calc_idx]     <= phase_mod[9:0];
                                amplitude[calc_idx] <= 8'd255;
                            end
                        end

                        // 同相+加窗
                        3'b100: begin
                            phase[calc_idx]     <= 10'd0;
                            amplitude[calc_idx] <= window_val;
                        end

                        // 默认
                        default: begin
                            phase[calc_idx]     <= 10'd0;
                            amplitude[calc_idx] <= 8'd255;
                        end
                    endcase

                    // 索引递增或完成
                    if (calc_idx == 5'd24)
                        state <= DONE_S;
                    else
                        calc_idx <= calc_idx_plus1;
                end

                //---------------------------------------------------------------
                // DONE: 输出有效
                //---------------------------------------------------------------
                DONE_S: begin
                    param_valid <= 1'b1;
                    state <= IDLE;
                end

                default: state <= IDLE;
            endcase
        end
    end

endmodule
