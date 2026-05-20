`timescale 1ns / 1ps
// 功能: WM8978G I2C 初始化配置驱动
// 架构: 4步过采样状态机, 100kHz SCL, 带有闭环 ACK 监测与自动熔断重试机制

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

  // WM8978 纯净音乐线路采集版配置
  always @(*) begin
    case (lut_index)
      4'd0: lut_data = {7'h00, 9'h000};  // R0: 软复位
      4'd1: lut_data = {7'h01, 9'h033};  // VMID=75k, BUFIOEN=1, BIASEN=1
      4'd2: lut_data = {7'h02, 9'h00F};  // 只开启左右 ADC，关闭 PGA
      4'd3: lut_data = {7'h04, 9'h010};  // 标准 I2S, 16-bit
      4'd4: lut_data = {7'h06, 9'h000};  // 时钟设置
      4'd5: lut_data = {7'h0E, 9'h100};  // 开启 ADC 高通滤波器，滤除直流
      4'd6: lut_data = {7'h20, 9'h000};  // ALC 关闭
      4'd7: lut_data = {7'h2C, 9'h000};  // 断开 MIC 物理连接
      4'd8: lut_data = {7'h2F, 9'h070};  // 左声道 L2 直通 ADC (6dB)
      4'd9: lut_data = {7'h30, 9'h070};  // 右声道 R2 直通 ADC (6dB)
      // R15/R16: ADC 数字音量控制 (ADCVOL)。
      // 9'h1FF 包含 ADCVU=1 (同步更新) 以及 VOL=FF (+25.875dB 满血拉满)
      4'd10: lut_data = {7'h0F, 9'h1FF};  // R15: 左声道 ADC 数字音量拉满
      4'd11: lut_data = {7'h10, 9'h1FF};  // R16: 右声道 ADC 数字音量拉满
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
      ack_error <= 0;
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
              2: begin  /* 维持数据 */
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
                if (i2c_sda == 1'b1) ack_error <= 1'b1;
              end  // 采样 NACK
              3: begin
                i2c_scl <= 0;
                // [核心修改点] 熔断拦截与自动重试机制
                if (ack_error) begin
                  state <= IDLE;
                  lut_index <= 0;  // 重头再来
                  ack_error <= 0;  // 清除错误标志
                end else begin
                  state   <= SEND_BYTE1;
                  bit_cnt <= 7;
                end
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
              2: begin  /* 维持数据 */
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
                if (i2c_sda == 1'b1) ack_error <= 1'b1;
              end  // 采样 NACK
              3: begin
                i2c_scl <= 0;
                // [核心修改点] 熔断拦截与自动重试机制
                if (ack_error) begin
                  state <= IDLE;
                  lut_index <= 0;
                  ack_error <= 0;
                end else begin
                  state   <= SEND_BYTE2;
                  bit_cnt <= 7;
                end
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
              2: begin  /* 维持数据 */
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
                if (i2c_sda == 1'b1) ack_error <= 1'b1;
              end  // 采样 NACK
              3: begin
                i2c_scl <= 0;
                // [核心修改点] 熔断拦截与自动重试机制
                if (ack_error) begin
                  state <= IDLE;
                  lut_index <= 0;
                  ack_error <= 0;
                end else begin
                  state <= STOP;
                end
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
