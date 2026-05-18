`timescale 1ns / 1ps
// 功能: WM8978G I2C 初始化配置驱动
// 架构: 4步过采样状态机, 100kHz SCL, 带有闭环 ACK 监测

module wm8978_i2c_init (
    input wire clk,   // 100MHz 系统时钟
    input wire rst_n, // 异步复位 (低有效)

    output reg  i2c_scl,  // I2C 时钟线
    inout  wire i2c_sda,  // I2C 数据线 (双向)

    output reg init_done,  // 初始化完成标志
    output reg ack_error   //  (高电平表示发生过NACK)
);

  // 100MHz / (250 * 4) = 100kHz I2C 标准频率
  parameter CLK_DIVIDE = 10'd250;
  parameter DEVICE_ID = 8'h34;  // WM8978 物理地址

  localparam IDLE         = 4'd0,
               START        = 4'd1,
               SEND_ID      = 4'd2,
               ACK_1        = 4'd3,
               SEND_BYTE1   = 4'd4, 
               ACK_2        = 4'd5,
               SEND_BYTE2   = 4'd6, 
               ACK_3        = 4'd7,
               STOP         = 4'd8,
               WAIT         = 4'd9,
               DONE         = 4'd10;

  reg [ 3:0] state;
  reg [ 9:0] clk_cnt;
  reg [ 1:0] i2c_clk_cnt;  // 实现 0~3 自然循环溢出
  reg [ 4:0] bit_cnt;
  reg [15:0] lut_data;
  reg [ 5:0] lut_index;
  reg [ 8:0] wait_cnt;

  reg        i2c_sda_reg;
  reg        i2c_sda_en;

  assign i2c_sda = i2c_sda_en ? i2c_sda_reg : 1'bz;

  // // WM8978 麦克风纯净采集版 (带硬件噪声门)
  // always @(*) begin
  //   case (lut_index)
  //     4'd0: lut_data = {7'h00, 9'h000};

  //     // --- 电源与模拟前端 ---
  //     4'd1: lut_data = {7'h01, 9'h01D};
  //     4'd2: lut_data = {7'h02, 9'h03F};

  //     // --- 接口与时钟 ---
  //     4'd3: lut_data = {7'h04, 9'h010};
  //     4'd4: lut_data = {7'h06, 9'h000};

  //     // --- 核心信号处理 ---
  //     4'd5: lut_data = {7'h0E, 9'h100};

  //     // --- ALC 与 噪声门 (除噪核心) ---
  //     // 降低了 ALC 的最大增益，防止底噪被过度放大
  //     4'd6: lut_data = {7'h20, 9'h128};  // R32: ALC 开启，最大增益 +23dB
  //     4'd7: lut_data = {7'h21, 9'h032};  // R33: ALC 目标电平 -12dBFS
  //     4'd8: lut_data = {7'h22, 9'h132};  // R34: ALC Limiter 模式
  //     // 【新增】当声音低于 -72dB 时直接静音，彻底消除环境沙沙声
  //     4'd9: lut_data = {7'h23, 9'h00B};  // R35: 开启 Noise Gate

  //     // --- 模拟输入路由 ---
  //     4'd10: lut_data = {7'h2C, 9'h033};
  //     4'd11: lut_data = {7'h2D, 9'h110};  // R45: 初始 0dB
  //     4'd12: lut_data = {7'h2E, 9'h110};

  //     // --- 适度 BOOST 增益 (+20dB) ---
  //     4'd13: lut_data = {7'h2F, 9'h170};
  //     4'd14: lut_data = {7'h30, 9'h170};

  //     default: lut_data = 16'hFFFF;
  //   endcase
  // end

  // WM8978 纯净音乐线路采集版 (3.5mm Line-In 专用)
  always @(*) begin
    case (lut_index)
      4'd0: lut_data = {7'h00, 9'h000};  // R0: 软复位

      // --- 电源管理 (关闭麦克风偏置，只留核心电压) ---
      // VMID=75k, BUFIOEN=1, BIASEN=1 -> 0x00F
      4'd1: lut_data = {7'h01, 9'h00F};
      // 只开启左右 ADC，关闭 PGA 放大器 (Line-In不需要放大) -> 0x003
      4'd2: lut_data = {7'h02, 9'h033};

      // --- 接口与时钟 (维持原样，FPGA 主机 16-bit I2S) ---
      4'd3: lut_data = {7'h04, 9'h010};
      4'd4: lut_data = {7'h06, 9'h000};

      // --- 核心信号处理 ---
      4'd5: lut_data = {7'h0E, 9'h100};  // R14: 开启 ADC 高通滤波器，滤除直流偏置

      // --- 动态控制 (听音乐必须关闭 ALC，保留音乐原本的强弱动态) ---
      4'd6: lut_data = {7'h20, 9'h000};  // R32: ALC 关闭

      // --- 模拟输入路由 (彻底切换到 3.5mm Line-In 通道) ---
      4'd7: lut_data = {7'h2C, 9'h000};  // R44: 断开所有 MIC 物理连接

      // R47/R48: 将 3.5mm 接口 (L2/R2) 直接连入 ADC，增益设为 0dB (Bit 6:4 = 101)
      4'd8: lut_data = {7'h2F, 9'h050};  // R47: 左声道 L2 直通 ADC (0dB)
      4'd9: lut_data = {7'h30, 9'h050};  // R48: 右声道 R2 直通 ADC (0dB)

      default: lut_data = 16'hFFFF;
    endcase
  end

  always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      state <= IDLE;
      clk_cnt <= 0;
      i2c_clk_cnt <= 0;
      bit_cnt <= 0;
      lut_index <= 0;
      wait_cnt <= 0;
      i2c_scl <= 1;
      i2c_sda_reg <= 1;
      i2c_sda_en <= 1;
      init_done <= 0;
      ack_error <= 0;  // 复位错误标志
    end else begin
      if (clk_cnt < CLK_DIVIDE - 1) begin
        clk_cnt <= clk_cnt + 1;
      end else begin
        clk_cnt <= 0;

        case (state)
          IDLE: begin
            if (lut_data != 16'hFFFF) state <= START;
            else begin
              state <= DONE;
              init_done <= 1;
            end
          end

          START: begin
            case (i2c_clk_cnt)
              0: begin
                i2c_sda_reg <= 1;
                i2c_sda_en  <= 1;
              end
              1: begin
                i2c_scl <= 1;
              end
              2: begin
                i2c_sda_reg <= 0;
              end
              3: begin
                i2c_scl <= 0;
                state   <= SEND_ID;
                bit_cnt <= 7;
              end
            endcase
            i2c_clk_cnt <= i2c_clk_cnt + 1'b1;
          end

          SEND_ID: begin
            case (i2c_clk_cnt)
              0: begin
                i2c_sda_reg <= DEVICE_ID[bit_cnt];
              end
              1: begin
                i2c_scl <= 1;
              end
              2: begin  /* 维持总线状态，满足数据保持时间要求 */
              end
              3: begin
                i2c_scl <= 0;
                if (bit_cnt > 0) bit_cnt <= bit_cnt - 1;
                else state <= ACK_1;
              end
            endcase
            i2c_clk_cnt <= i2c_clk_cnt + 1'b1;
          end

          ACK_1: begin
            case (i2c_clk_cnt)
              0: begin
                i2c_sda_en <= 0;
              end
              1: begin
                i2c_scl <= 1;
              end
              2: begin
                // 总线释放且SCL为高时采样，高电平代表 NACK
                if (i2c_sda == 1'b1) ack_error <= 1'b1;
              end
              3: begin
                i2c_scl <= 0;
                state   <= SEND_BYTE1;
                bit_cnt <= 7;
              end
            endcase
            i2c_clk_cnt <= i2c_clk_cnt + 1'b1;
          end

          SEND_BYTE1: begin
            case (i2c_clk_cnt)
              0: begin
                i2c_sda_en  <= 1;
                i2c_sda_reg <= (bit_cnt > 0) ? lut_data[bit_cnt+8] : lut_data[8];
              end
              1: begin
                i2c_scl <= 1;
              end
              2: begin  /* 维持总线状态，满足数据保持时间要求 */
              end
              3: begin
                i2c_scl <= 0;
                if (bit_cnt > 0) bit_cnt <= bit_cnt - 1;
                else state <= ACK_2;
              end
            endcase
            i2c_clk_cnt <= i2c_clk_cnt + 1'b1;
          end

          ACK_2: begin
            case (i2c_clk_cnt)
              0: begin
                i2c_sda_en <= 0;
              end
              1: begin
                i2c_scl <= 1;
              end
              2: begin
                // 采样应答
                if (i2c_sda == 1'b1) ack_error <= 1'b1;
              end
              3: begin
                i2c_scl <= 0;
                state   <= SEND_BYTE2;
                bit_cnt <= 7;
              end
            endcase
            i2c_clk_cnt <= i2c_clk_cnt + 1'b1;
          end

          SEND_BYTE2: begin
            case (i2c_clk_cnt)
              0: begin
                i2c_sda_en  <= 1;
                i2c_sda_reg <= lut_data[bit_cnt];
              end
              1: begin
                i2c_scl <= 1;
              end
              2: begin  /* 维持总线状态，满足数据保持时间要求 */
              end
              3: begin
                i2c_scl <= 0;
                if (bit_cnt > 0) bit_cnt <= bit_cnt - 1;
                else state <= ACK_3;
              end
            endcase
            i2c_clk_cnt <= i2c_clk_cnt + 1'b1;
          end

          ACK_3: begin
            case (i2c_clk_cnt)
              0: begin
                i2c_sda_en <= 0;
              end
              1: begin
                i2c_scl <= 1;
              end
              2: begin
                // 采样应答
                if (i2c_sda == 1'b1) ack_error <= 1'b1;
              end
              3: begin
                i2c_scl <= 0;
                state   <= STOP;
              end
            endcase
            i2c_clk_cnt <= i2c_clk_cnt + 1'b1;
          end

          STOP: begin
            case (i2c_clk_cnt)
              0: begin
                i2c_sda_en  <= 1;
                i2c_sda_reg <= 0;
              end
              1: begin
                i2c_scl <= 1;
              end
              2: begin
                i2c_sda_reg <= 1;
              end
              3: begin
                state <= WAIT;
                wait_cnt <= 0;
              end
            endcase
            i2c_clk_cnt <= i2c_clk_cnt + 1'b1;
          end

          WAIT: begin
            if (wait_cnt < 9'd400) begin
              wait_cnt <= wait_cnt + 1;
            end else begin
              wait_cnt <= 0;
              lut_index <= lut_index + 1;
              state <= IDLE;
            end
          end

          DONE: begin
            init_done <= 1;
          end
        endcase
      end
    end
  end
endmodule
