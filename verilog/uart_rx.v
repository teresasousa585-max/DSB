`timescale 1ns / 1ps
// 语言: Verilog-2001

module uart_rx #(
    parameter CLK_FRE    = 100,     // 主频 100MHz
    parameter BAUD_RATE  = 115200,  // 默认目标波特率改为 115200
    parameter DATA_WIDTH = 8        // 数据位宽
) (
    input clk,       // 系统时钟 (原 i_clk_sys)
    input rst_n,     // 低电平同步复位 (原 i_rst)
    input i_uart_rx, // UART接收引脚

    output reg [DATA_WIDTH-1 : 0] o_uart_data,  // 接收到的 1 字节数据
    output reg                    o_rx_done     // 接收完成单脉冲
);

  // 动态计算周期和采样中心点
  localparam CLK_COUNT = (CLK_FRE * 1000000) / BAUD_RATE;
  localparam MID_PT = CLK_COUNT / 2;

  // 三级打拍同步，消除跨时钟域亚稳态
  reg [2:0] rx_reg;
  always @(posedge clk or negedge rst_n) begin
    if (!rst_n) rx_reg <= 3'b111;
    else rx_reg <= {rx_reg[1:0], i_uart_rx};
  end

  wire rx_falling = (rx_reg[2] & !rx_reg[1]);  // 捕捉起始位下降沿
  wire rx_sync = rx_reg[2];  // 同步后的接收信号

  // 状态机定义
  localparam S_IDLE = 3'd0;
  localparam S_START = 3'd1;
  localparam S_DATA = 3'd2;
  localparam S_STOP = 3'd3;

  reg [2:0] state;
  reg [15:0] cycle_cnt;
  reg [3:0] bit_cnt;
  reg [DATA_WIDTH-1:0] rx_temp;

  // 多数表决采样寄存器
  reg smp1, smp2, smp3;
  // 表决逻辑：3次采样中至少有2次为1，结果才为1；否则为0
  wire bit_val = (smp1 & smp2) | (smp1 & smp3) | (smp2 & smp3);

  always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      state       <= S_IDLE;
      o_rx_done   <= 1'b0;
      o_uart_data <= 0;
      cycle_cnt   <= 0;
      bit_cnt     <= 0;
      smp1        <= 1'b1;
      smp2        <= 1'b1;
      smp3        <= 1'b1;
    end else begin
      o_rx_done <= 1'b0;  // 默认拉低

      case (state)
        S_IDLE: begin
          cycle_cnt <= 0;
          bit_cnt   <= 0;
          if (rx_falling) state <= S_START;
        end

        S_START: begin
          cycle_cnt <= cycle_cnt + 1'b1;
          // 自适应中心点采样
          if (cycle_cnt == MID_PT - 5) smp1 <= rx_sync;
          if (cycle_cnt == MID_PT) smp2 <= rx_sync;
          if (cycle_cnt == MID_PT + 5) smp3 <= rx_sync;

          if (cycle_cnt == CLK_COUNT - 1) begin
            cycle_cnt <= 0;
            if (!bit_val) state <= S_DATA;  // 起始位必须为0
            else state <= S_IDLE;
          end
        end

        S_DATA: begin
          cycle_cnt <= cycle_cnt + 1'b1;
          // 自适应中心点采样
          if (cycle_cnt == MID_PT - 5) smp1 <= rx_sync;
          if (cycle_cnt == MID_PT) smp2 <= rx_sync;
          if (cycle_cnt == MID_PT + 5) smp3 <= rx_sync;

          if (cycle_cnt == CLK_COUNT - 1) begin
            cycle_cnt <= 0;
            rx_temp[bit_cnt] <= bit_val;

            if (bit_cnt == DATA_WIDTH - 1) begin
              bit_cnt <= 0;
              state   <= S_STOP;
            end else begin
              bit_cnt <= bit_cnt + 1'b1;
            end
          end
        end

        S_STOP: begin
          cycle_cnt <= cycle_cnt + 1'b1;
          // 自适应中心点采样
          if (cycle_cnt == MID_PT - 5) smp1 <= rx_sync;
          if (cycle_cnt == MID_PT) smp2 <= rx_sync;
          if (cycle_cnt == MID_PT + 5) smp3 <= rx_sync;

          if (cycle_cnt == CLK_COUNT - 1) begin
            cycle_cnt <= 0;
            if (bit_val == 1'b1) begin
              o_uart_data <= rx_temp;
              o_rx_done   <= 1'b1;  // 接收成功
            end
            state <= S_IDLE;
          end
        end

        default: state <= S_IDLE;
      endcase
    end
  end
endmodule
