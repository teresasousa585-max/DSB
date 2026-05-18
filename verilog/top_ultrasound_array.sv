`timescale 1ns / 1ps

module top_ultrasound_array (
    input  logic clk_50m,     
    input  logic async_rst_n, 

    input  logic uart_rx,
    input  logic key_in,

    // 32 路 PWM 直驱输出
    output logic [31:0] transducer_io,

    output logic led,

    // WM8978 I2C接口
    output logic i2c_scl,       // I2C 时钟
    inout  wire  i2c_sda,       // I2C 数据
    output logic config_led,    //I2C 配置成功指示灯

    // WM8978 I2S 接口
    output logic clk_mclk,
    output logic wm8978_bclk,
    output logic wm8978_lrck,
    input  logic wm8978_adcdat
);

    // 内部倍频到 100MHz 
    logic clk_100m;
    logic locked;
    mmcm_main u_mmcm (
        .clk_in1  (clk_50m),
        .clk_out1 (clk_100m),
        .clk_out2 (clk_mclk),//12.288MHZ
        .locked   (locked)
    );

    logic rst_n;
    assign rst_n = async_rst_n & locked;


    // 2. WM8978 I2C 初始化 (配置为从机, 16bit, 48k, 开启ALC)
    logic i2c_init_done, i2c_error;
    wm8978_i2c_init u_conf (
        .clk       (clk_100m),
        .rst_n     (rst_n),
        .i2c_scl   (i2c_scl),
        .i2c_sda   (i2c_sda),
        .init_done (i2c_init_done),
        .ack_error (i2c_error)
    );
// 逻辑：!(完成 && 无错) -> 只有完成且无错时输出0(亮)
assign config_led = !(i2c_init_done && !i2c_error);

    // 1. 按键消抖与模式切换
    logic key_flag;
    key_filter u_key (
        .clk       (clk_100m),
        .rst       (~rst_n),    
        .key_in    (key_in),
        .key_flag  (key_flag),  
        .key_state ()
    );

    logic audio_mode; // 0: 播放内部 DDS (1kHz), 1: 播放外部音频
    always_ff @(posedge clk_100m or negedge rst_n) begin
        if (!rst_n) begin
            audio_mode <= 1'b0;
        end else if (key_flag) begin
            audio_mode <= ~audio_mode;
        end
    end

    assign led = audio_mode;

    // 3. I2S 接收器 (12.288MHz 时钟域)
    logic signed [15:0] i2s_left_12m, i2s_right_12m;
    logic i2s_valid_12m;
    
    i2s_rx_master u_i2s (
        .rst_n      (rst_n),
        .mclk       (clk_mclk),
        .bclk       (wm8978_bclk),
        .lrck       (wm8978_lrck),
        .adcdat     (wm8978_adcdat),
        .left_data  (i2s_left_12m),
        .right_data (i2s_right_12m),
        .data_valid (i2s_valid_12m)  // 这是 12.288MHz 域的脉冲
    );

    //跨时钟域处理 (CDC: 12.288MHz -> 100MHz) 
    logic [2:0] valid_sync_100m;
    logic signed [15:0] i2s_left_100m;

    always_ff @(posedge clk_100m or negedge rst_n) begin
        if (!rst_n) begin
            valid_sync_100m <= 3'b0;
            i2s_left_100m   <= 16'd0;
        end else begin
            // 经典打拍：将 12M 域的 valid 信号同步到 100M 域
            valid_sync_100m <= {valid_sync_100m[1:0], i2s_valid_12m};
            
            // 边沿检测：当在 100M 域检测到 valid 的上升沿时
            // 此时 i2s_left_12m 的数据早就稳定了几十纳秒，读取绝对安全！
            if (valid_sync_100m[2:1] == 2'b01) begin
                i2s_left_100m <= i2s_left_12m;
            end
        end
    end

    // 3. DDS 1kHz 调制信号 ( 16位输出)
    logic signed [15:0] sin_1k_16b;
    dds_generator u_dds (
        .clk    (clk_100m),
        .rst_n  (rst_n),
        .sin_1k (sin_1k_16b)
    );

    // // 4. 音频流安全放大与多路选择 (MUX)
    // logic signed [15:0] final_audio_stream;
    // logic signed [21:0] scaled_audio; // 拓宽位宽防止计算溢出
    
    // // 给麦克风信号乘以 32 倍
    // assign scaled_audio = i2s_left_100m *4; 

    // // 饱和截断逻辑 (防爆音限幅器)
    // always_comb begin
    //     if (audio_mode == 1'b0) begin
    //         final_audio_stream = sin_1k_16b; // 1kHz 模式
    //     end else begin
    //         // 麦克风模式：判断放大后的数据是否超出了 16-bit 有符号数的极限
    //         if (scaled_audio > 22'sd32767) 
    //             final_audio_stream = 16'sd32767;     // 正向削峰
    //         else if (scaled_audio < -22'sd32768) 
    //             final_audio_stream = -16'sd32768;    // 负向削峰
    //         else 
    //             final_audio_stream = scaled_audio[15:0]; // 安全范围内，正常截取
    //     end
    // end

// =========================================================
    // 6. 音频流处理核心 (噪声门 -> 放大 -> 饱和截断)
    // =========================================================
    
    // Stage 1: 数字噪声门
    logic signed [15:0] gated_audio;
    localparam signed [15:0] NOISE_THRESHOLD = 16'sd200; // 根据实际底噪调整
    
    always_comb begin
        // 提取绝对值进行判断，滤除微小底噪
        if (i2s_left_100m > NOISE_THRESHOLD || i2s_left_100m < -NOISE_THRESHOLD)
            gated_audio = i2s_left_100m; 
        else
            gated_audio = 16'sd0; // 纯净待机
    end

    // Stage 2: 拓宽位宽并进行增益放大
    logic signed [21:0] scaled_audio;
    // 这里的放大倍数可以根据需要调整 (例如 4, 8, 16)
    // 既然 *1 太小，*8 足够响，建议保留 *8 或 *4

    assign scaled_audio = gated_audio; 
ila_0 u_ila_0 (
	.clk(clk_100m), // input wire clk


	.probe0(scaled_audio) // input wire [20:0] probe0
);
    // Stage 3: 多路选择与饱和限幅器 (防爆音)
    logic signed [15:0] final_audio_stream;
    
    always_comb begin
        if (audio_mode == 1'b0) begin
            final_audio_stream = sin_1k_16b; // 1kHz DDS 测试模式
        end else begin
            // 外部音频模式：执行硬限幅，绝对禁止数据溢出回绕
            if (scaled_audio > 22'sd32767) 
                final_audio_stream = 16'sd32767;      // 正向削峰
            else if (scaled_audio < -22'sd32768) 
                final_audio_stream = -16'sd32768;     // 负向削峰
            else 
                final_audio_stream = scaled_audio[15:0]; // 安全区间原样输出
        end
    end

    // 5. 串口协议解析 (幅度与相位控制)
    logic [7:0]  beam_amplitude [0:31];
    logic [11:0] beam_phase     [0:31];
    logic [7:0]  uart_byte;
    logic        uart_done;
    
    uart_rx #(.CLK_FRE(100),
     .BAUD_RATE(115200)) u_uart_rx (
        .clk(clk_100m), 
        .rst_n(rst_n),
         .i_uart_rx(uart_rx),
        .o_uart_data(uart_byte),
         .o_rx_done(uart_done)
    );

    uart_protocol_parser u_parser (
        .clk(clk_100m), 
        .rst_n(rst_n), 
        .rx_data(uart_byte),
         .rx_done(uart_done),
        .o_amplitude(beam_amplitude), 
        .o_phase(beam_phase)
    );

    // 6. 核心调制与驱动 (接收 16位音频流)
    pwm32_generator u_pwm32 (
        .clk       (clk_100m),
        .rst_n     (rst_n),
        .audio_in  (final_audio_stream), // 16-bit 接口
        .amplitude (beam_amplitude),
        .phase_del (beam_phase),
        .pwm_out   (transducer_io)
    );


endmodule