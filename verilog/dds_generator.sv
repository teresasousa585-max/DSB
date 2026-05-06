//////////////////////////////////////////////////////////////////////////////
// Module: dds_generator
// Description: DDS (Direct Digital Synthesis) Dual-Channel Sine Wave Generator
// 
// Features:
//   - Dual independent DDS channels using 32-bit phase accumulators
//   - Channel 1: 1kHz modulation signal (for DSB-AM testing)
//   - Channel 2: 40kHz carrier phase reference (for PWM synchronization)
//   - 256-point sine LUT with 12-bit signed output (-2048 ~ 2047)
//   - Carrier cycle tick signal for synchronizing 25-channel PWM generators
// 
// Parameters:
//   - SYS_CLK_FREQ: 50MHz system clock
//   - PHASE_WIDTH: 32-bit phase accumulator
//   - LUT_ADDR_WIDTH: 8-bit LUT address (256 points)
//   - DATA_WIDTH: 12-bit signed output
// 
// Frequency Tuning Word (FTW) Calculation:
//   FTW = Fout * 2^N / Fs
//   FTW_1k  = 1000  * 2^32 / 50e6 = 85899   (0x00014F8B)
//   FTW_40k = 40000 * 2^32 / 50e6 = 3435974 (0x00346DC6)
//
// Author: Auto-generated
// Date: 2024
//////////////////////////////////////////////////////////////////////////////

module dds_generator (
    input        clk,       // 50MHz system clock
    input        rst_n,     // Active-low asynchronous reset
    input        en,        // Global enable signal

    // 1kHz modulation signal output (for DSB-AM modulation)
    output reg signed [11:0] sin_1k,         // 1kHz sine wave: modulation signal
    output reg               sin_1k_valid,   // 1kHz output data valid flag

    // 40kHz carrier phase reference (for PWM carrier synchronization)
    output reg signed [11:0] sin_40k,        // 40kHz sine wave: carrier reference
    output reg               sin_40k_valid,  // 40kHz output data valid flag

    // Carrier cycle marker (used for synchronizing 25-channel PWM generator)
    // Outputs a single-clock pulse at the start of each 40kHz carrier cycle
    output reg               carrier_cycle_tick
);

    //-------------------------------------------------------------------------
    // Parameters
    //-------------------------------------------------------------------------
    parameter SYS_CLK_FREQ = 50_000_000;  // System clock frequency: 50MHz
    parameter PHASE_WIDTH  = 32;          // Phase accumulator bit width
    parameter LUT_DEPTH    = 256;         // LUT depth: 256 points per cycle
    parameter LUT_ADDR_W   = 8;           // LUT address width: log2(256) = 8
    parameter DATA_WIDTH   = 12;          // Output data bit width (signed)

    // Frequency Tuning Words (FTW) = Fout * 2^PHASE_WIDTH / SYS_CLK_FREQ
    // FTW_1k  = 1000  * 2^32 / 50MHz = 85899
    localparam [PHASE_WIDTH-1:0] FTW_1K  = 32'd85899;    // 0x00014F8B
    // FTW_40k = 40000 * 2^32 / 50MHz = 3435974
    localparam [PHASE_WIDTH-1:0] FTW_40K = 32'd3435974;  // 0x00346DC6

    // Phase accumulator overflow threshold for carrier cycle detection
    // When phase_acc > THRESHOLD, the next addition will cause overflow
    localparam [PHASE_WIDTH-1:0] CARRIER_OVERFLOW_THRESH = 32'hFFFFFFFF - FTW_40K + 1;

    //-------------------------------------------------------------------------
    // Internal Signals
    //-------------------------------------------------------------------------
    reg [PHASE_WIDTH-1:0] phase_acc_1k;   // 1kHz channel phase accumulator
    reg [PHASE_WIDTH-1:0] phase_acc_40k;  // 40kHz channel phase accumulator

    wire [LUT_ADDR_W-1:0] lut_addr_1k;    // 1kHz LUT address (upper 8 bits of phase)
    wire [LUT_ADDR_W-1:0] lut_addr_40k;   // 40kHz LUT address (upper 8 bits of phase)

    wire signed [DATA_WIDTH-1:0] lut_data_1k;   // LUT output for 1kHz channel
    wire signed [DATA_WIDTH-1:0] lut_data_40k;  // LUT output for 40kHz channel

    //-------------------------------------------------------------------------
    // 256-Point Sine Lookup Table (12-bit signed)
    // Formula: round(2047 * sin(2*pi*n/256)), n = 0,1,2,...,255
    // Range: -2048 ~ +2047
    // The LUT covers a complete 0~2*pi cycle:
    //   [0:63]   : 0 to pi/2     (1st quadrant, 0 to +2047)
    //   [64:127] : pi/2 to pi    (2nd quadrant, +2047 to 0)
    //   [128:191]: pi to 3pi/2   (3rd quadrant, 0 to -2047)
    //   [192:255]: 3pi/2 to 2pi  (4th quadrant, -2047 to 0)
    //-------------------------------------------------------------------------
    reg signed [DATA_WIDTH-1:0] sin_lut [0:LUT_DEPTH-1];

    initial begin
        // Initialize sine LUT with pre-calculated values
        // n=0~127: positive half cycle (0 to +2047 to 0)
        sin_lut[0]   = 12'd0;      sin_lut[1]   = 12'd50;     sin_lut[2]   = 12'd100;    sin_lut[3]   = 12'd151;
        sin_lut[4]   = 12'd201;    sin_lut[5]   = 12'd251;    sin_lut[6]   = 12'd300;    sin_lut[7]   = 12'd350;
        sin_lut[8]   = 12'd399;    sin_lut[9]   = 12'd449;    sin_lut[10]  = 12'd497;    sin_lut[11]  = 12'd546;
        sin_lut[12]  = 12'd594;    sin_lut[13]  = 12'd642;    sin_lut[14]  = 12'd690;    sin_lut[15]  = 12'd737;
        sin_lut[16]  = 12'd783;    sin_lut[17]  = 12'd830;    sin_lut[18]  = 12'd875;    sin_lut[19]  = 12'd920;
        sin_lut[20]  = 12'd965;    sin_lut[21]  = 12'd1009;   sin_lut[22]  = 12'd1052;   sin_lut[23]  = 12'd1095;
        sin_lut[24]  = 12'd1137;   sin_lut[25]  = 12'd1179;   sin_lut[26]  = 12'd1219;   sin_lut[27]  = 12'd1259;
        sin_lut[28]  = 12'd1299;   sin_lut[29]  = 12'd1337;   sin_lut[30]  = 12'd1375;   sin_lut[31]  = 12'd1411;
        sin_lut[32]  = 12'd1447;   sin_lut[33]  = 12'd1483;   sin_lut[34]  = 12'd1517;   sin_lut[35]  = 12'd1550;
        sin_lut[36]  = 12'd1582;   sin_lut[37]  = 12'd1614;   sin_lut[38]  = 12'd1644;   sin_lut[39]  = 12'd1674;
        sin_lut[40]  = 12'd1702;   sin_lut[41]  = 12'd1729;   sin_lut[42]  = 12'd1756;   sin_lut[43]  = 12'd1781;
        sin_lut[44]  = 12'd1805;   sin_lut[45]  = 12'd1828;   sin_lut[46]  = 12'd1850;   sin_lut[47]  = 12'd1871;
        sin_lut[48]  = 12'd1891;   sin_lut[49]  = 12'd1910;   sin_lut[50]  = 12'd1927;   sin_lut[51]  = 12'd1944;
        sin_lut[52]  = 12'd1959;   sin_lut[53]  = 12'd1973;   sin_lut[54]  = 12'd1986;   sin_lut[55]  = 12'd1997;
        sin_lut[56]  = 12'd2008;   sin_lut[57]  = 12'd2017;   sin_lut[58]  = 12'd2025;   sin_lut[59]  = 12'd2032;
        sin_lut[60]  = 12'd2037;   sin_lut[61]  = 12'd2041;   sin_lut[62]  = 12'd2045;   sin_lut[63]  = 12'd2046;
        sin_lut[64]  = 12'd2047;   sin_lut[65]  = 12'd2046;   sin_lut[66]  = 12'd2045;   sin_lut[67]  = 12'd2041;
        sin_lut[68]  = 12'd2037;   sin_lut[69]  = 12'd2032;   sin_lut[70]  = 12'd2025;   sin_lut[71]  = 12'd2017;
        sin_lut[72]  = 12'd2008;   sin_lut[73]  = 12'd1997;   sin_lut[74]  = 12'd1986;   sin_lut[75]  = 12'd1973;
        sin_lut[76]  = 12'd1959;   sin_lut[77]  = 12'd1944;   sin_lut[78]  = 12'd1927;   sin_lut[79]  = 12'd1910;
        sin_lut[80]  = 12'd1891;   sin_lut[81]  = 12'd1871;   sin_lut[82]  = 12'd1850;   sin_lut[83]  = 12'd1828;
        sin_lut[84]  = 12'd1805;   sin_lut[85]  = 12'd1781;   sin_lut[86]  = 12'd1756;   sin_lut[87]  = 12'd1729;
        sin_lut[88]  = 12'd1702;   sin_lut[89]  = 12'd1674;   sin_lut[90]  = 12'd1644;   sin_lut[91]  = 12'd1614;
        sin_lut[92]  = 12'd1582;   sin_lut[93]  = 12'd1550;   sin_lut[94]  = 12'd1517;   sin_lut[95]  = 12'd1483;
        sin_lut[96]  = 12'd1447;   sin_lut[97]  = 12'd1411;   sin_lut[98]  = 12'd1375;   sin_lut[99]  = 12'd1337;
        sin_lut[100] = 12'd1299;   sin_lut[101] = 12'd1259;   sin_lut[102] = 12'd1219;   sin_lut[103] = 12'd1179;
        sin_lut[104] = 12'd1137;   sin_lut[105] = 12'd1095;   sin_lut[106] = 12'd1052;   sin_lut[107] = 12'd1009;
        sin_lut[108] = 12'd965;    sin_lut[109] = 12'd920;    sin_lut[110] = 12'd875;    sin_lut[111] = 12'd830;
        sin_lut[112] = 12'd783;    sin_lut[113] = 12'd737;    sin_lut[114] = 12'd690;    sin_lut[115] = 12'd642;
        sin_lut[116] = 12'd594;    sin_lut[117] = 12'd546;    sin_lut[118] = 12'd497;    sin_lut[119] = 12'd449;
        sin_lut[120] = 12'd399;    sin_lut[121] = 12'd350;    sin_lut[122] = 12'd300;    sin_lut[123] = 12'd251;
        sin_lut[124] = 12'd201;    sin_lut[125] = 12'd151;    sin_lut[126] = 12'd100;    sin_lut[127] = 12'd50;

        // n=128~255: negative half cycle (0 to -2047 to 0)
        sin_lut[128] = 12'sd0;      sin_lut[129] = -12'sd50;    sin_lut[130] = -12'sd100;
        sin_lut[131] = -12'sd151;   sin_lut[132] = -12'sd201;   sin_lut[133] = -12'sd251;
        sin_lut[134] = -12'sd300;   sin_lut[135] = -12'sd350;   sin_lut[136] = -12'sd399;
        sin_lut[137] = -12'sd449;   sin_lut[138] = -12'sd497;   sin_lut[139] = -12'sd546;
        sin_lut[140] = -12'sd594;   sin_lut[141] = -12'sd642;   sin_lut[142] = -12'sd690;
        sin_lut[143] = -12'sd737;   sin_lut[144] = -12'sd783;   sin_lut[145] = -12'sd830;
        sin_lut[146] = -12'sd875;   sin_lut[147] = -12'sd920;   sin_lut[148] = -12'sd965;
        sin_lut[149] = -12'sd1009;  sin_lut[150] = -12'sd1052;  sin_lut[151] = -12'sd1095;
        sin_lut[152] = -12'sd1137;  sin_lut[153] = -12'sd1179;  sin_lut[154] = -12'sd1219;
        sin_lut[155] = -12'sd1259;  sin_lut[156] = -12'sd1299;  sin_lut[157] = -12'sd1337;
        sin_lut[158] = -12'sd1375;  sin_lut[159] = -12'sd1411;  sin_lut[160] = -12'sd1447;
        sin_lut[161] = -12'sd1483;  sin_lut[162] = -12'sd1517;  sin_lut[163] = -12'sd1550;
        sin_lut[164] = -12'sd1582;  sin_lut[165] = -12'sd1614;  sin_lut[166] = -12'sd1644;
        sin_lut[167] = -12'sd1674;  sin_lut[168] = -12'sd1702;  sin_lut[169] = -12'sd1729;
        sin_lut[170] = -12'sd1756;  sin_lut[171] = -12'sd1781;  sin_lut[172] = -12'sd1805;
        sin_lut[173] = -12'sd1828;  sin_lut[174] = -12'sd1850;  sin_lut[175] = -12'sd1871;
        sin_lut[176] = -12'sd1891;  sin_lut[177] = -12'sd1910;  sin_lut[178] = -12'sd1927;
        sin_lut[179] = -12'sd1944;  sin_lut[180] = -12'sd1959;  sin_lut[181] = -12'sd1973;
        sin_lut[182] = -12'sd1986;  sin_lut[183] = -12'sd1997;  sin_lut[184] = -12'sd2008;
        sin_lut[185] = -12'sd2017;  sin_lut[186] = -12'sd2025;  sin_lut[187] = -12'sd2032;
        sin_lut[188] = -12'sd2037;  sin_lut[189] = -12'sd2041;  sin_lut[190] = -12'sd2045;
        sin_lut[191] = -12'sd2046;  sin_lut[192] = -12'sd2047;  sin_lut[193] = -12'sd2046;
        sin_lut[194] = -12'sd2045;  sin_lut[195] = -12'sd2041;  sin_lut[196] = -12'sd2037;
        sin_lut[197] = -12'sd2032;  sin_lut[198] = -12'sd2025;  sin_lut[199] = -12'sd2017;
        sin_lut[200] = -12'sd2008;  sin_lut[201] = -12'sd1997;  sin_lut[202] = -12'sd1986;
        sin_lut[203] = -12'sd1973;  sin_lut[204] = -12'sd1959;  sin_lut[205] = -12'sd1944;
        sin_lut[206] = -12'sd1927;  sin_lut[207] = -12'sd1910;  sin_lut[208] = -12'sd1891;
        sin_lut[209] = -12'sd1871;  sin_lut[210] = -12'sd1850;  sin_lut[211] = -12'sd1828;
        sin_lut[212] = -12'sd1805;  sin_lut[213] = -12'sd1781;  sin_lut[214] = -12'sd1756;
        sin_lut[215] = -12'sd1729;  sin_lut[216] = -12'sd1702;  sin_lut[217] = -12'sd1674;
        sin_lut[218] = -12'sd1644;  sin_lut[219] = -12'sd1614;  sin_lut[220] = -12'sd1582;
        sin_lut[221] = -12'sd1550;  sin_lut[222] = -12'sd1517;  sin_lut[223] = -12'sd1483;
        sin_lut[224] = -12'sd1447;  sin_lut[225] = -12'sd1411;  sin_lut[226] = -12'sd1375;
        sin_lut[227] = -12'sd1337;  sin_lut[228] = -12'sd1299;  sin_lut[229] = -12'sd1259;
        sin_lut[230] = -12'sd1219;  sin_lut[231] = -12'sd1179;  sin_lut[232] = -12'sd1137;
        sin_lut[233] = -12'sd1095;  sin_lut[234] = -12'sd1052;  sin_lut[235] = -12'sd1009;
        sin_lut[236] = -12'sd965;   sin_lut[237] = -12'sd920;   sin_lut[238] = -12'sd875;
        sin_lut[239] = -12'sd830;   sin_lut[240] = -12'sd783;   sin_lut[241] = -12'sd737;
        sin_lut[242] = -12'sd690;   sin_lut[243] = -12'sd642;   sin_lut[244] = -12'sd594;
        sin_lut[245] = -12'sd546;   sin_lut[246] = -12'sd497;   sin_lut[247] = -12'sd449;
        sin_lut[248] = -12'sd399;   sin_lut[249] = -12'sd350;   sin_lut[250] = -12'sd300;
        sin_lut[251] = -12'sd251;   sin_lut[252] = -12'sd201;   sin_lut[253] = -12'sd151;
        sin_lut[254] = -12'sd100;   sin_lut[255] = -12'sd50;
    end

    //-------------------------------------------------------------------------
    // Phase Accumulator: 1kHz Modulation Channel
    // Accumulates FTW_1K (85899) every clock cycle at 50MHz
    // Phase wraps around naturally every 2^32 counts (50MHz/85899 ~ 1.000004kHz)
    //-------------------------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            phase_acc_1k <= {PHASE_WIDTH{1'b0}};  // Clear phase accumulator on reset
        end else if (en) begin
            // Phase accumulation: adds frequency tuning word each clock cycle
            // The upper bits determine the sine phase, lower bits provide fractional precision
            phase_acc_1k <= phase_acc_1k + FTW_1K;
        end
    end

    //-------------------------------------------------------------------------
    // Phase Accumulator: 40kHz Carrier Channel
    // Accumulates FTW_40K (3435974) every clock cycle at 50MHz
    // Phase wraps around naturally every 2^32 counts (50MHz/3435974 ~ 40.000002kHz)
    //-------------------------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            phase_acc_40k <= {PHASE_WIDTH{1'b0}};  // Clear phase accumulator on reset
        end else if (en) begin
            phase_acc_40k <= phase_acc_40k + FTW_40K;
        end
    end

    //-------------------------------------------------------------------------
    // Carrier Cycle Tick Detection
    // Detects when the 40kHz phase accumulator overflows (wraps around).
    // This happens when: phase_acc_40k > (2^32 - 1 - FTW_40K)
    // At this point, the next addition will cause the accumulator to wrap,
    // marking the beginning of a new carrier cycle.
    // The tick is used to synchronize 25-channel PWM generators.
    //-------------------------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            carrier_cycle_tick <= 1'b0;
        end else if (en) begin
            // Detect imminent overflow: if current phase + FTW would exceed 2^32-1
            if (phase_acc_40k > (32'hFFFFFFFF - FTW_40K))
                carrier_cycle_tick <= 1'b1;  // Pulse at start of new carrier cycle
            else
                carrier_cycle_tick <= 1'b0;
        end else begin
            carrier_cycle_tick <= 1'b0;
        end
    end

    //-------------------------------------------------------------------------
    // LUT Address Generation
    // Use upper LUT_ADDR_W (8) bits of phase accumulator as LUT index.
    // 32-bit phase -> upper 8 bits give 256 LUT addresses per full cycle.
    // Phase resolution: 2^(32-8) = 2^24 fractional bits for smooth frequency tuning.
    //-------------------------------------------------------------------------
    // LUT address is the upper 8 bits of the 32-bit phase accumulator
    assign lut_addr_1k  = phase_acc_1k[PHASE_WIDTH-1:PHASE_WIDTH-LUT_ADDR_W];   // bits [31:24]
    assign lut_addr_40k = phase_acc_40k[PHASE_WIDTH-1:PHASE_WIDTH-LUT_ADDR_W];  // bits [31:24]

    //-------------------------------------------------------------------------
    // Sine LUT Lookup (asynchronous read)
    // LUT is implemented as distributed RAM (reg array), read is combinational.
    // Register output in the next always block for synchronous output.
    //-------------------------------------------------------------------------
    assign lut_data_1k  = sin_lut[lut_addr_1k];
    assign lut_data_40k = sin_lut[lut_addr_40k];

    //-------------------------------------------------------------------------
    // Output Register Stage
    // Register the LUT output for clean, synchronous outputs.
    // valid signals are held high when enabled.
    //-------------------------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            sin_1k       <= {DATA_WIDTH{1'b0}};
            sin_1k_valid <= 1'b0;
            sin_40k      <= {DATA_WIDTH{1'b0}};
            sin_40k_valid<= 1'b0;
        end else if (en) begin
            // 1kHz channel output
            sin_1k       <= lut_data_1k;    // Sine value from LUT
            sin_1k_valid <= 1'b1;           // Data is always valid when enabled

            // 40kHz channel output
            sin_40k       <= lut_data_40k;   // Sine value from LUT
            sin_40k_valid <= 1'b1;          // Data is always valid when enabled
        end else begin
            sin_1k_valid  <= 1'b0;
            sin_40k_valid <= 1'b0;
        end
    end

endmodule
// End of module dds_generator
