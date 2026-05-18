`timescale 1ns / 1ps

module key_filter #(
    parameter CLK_FREQ = 100_000_000,  // 系统时钟，默认50 MHz
    parameter DEBOUNCE_TIME_MS = 20  //消抖时间，默认 20 ms
) (
    input clk,
    input rst,  //连接经过复位桥处理的复位信号
    input key_in,
    output reg key_flag,  //按键有效标志（按下瞬间产生一个时钟周期的高脉冲）
    output reg key_state  //消除抖动后的稳定按键状态
);

  //计算计数值
  localparam CNT_MAX = (CLK_FREQ / 1000) * DEBOUNCE_TIME_MS;
  //计算位宽
  localparam CNT_WIDTH = $clog2(CNT_MAX);

  // 状态机定义（独热码）
  localparam IDLE = 4'b0001, FILTER0 = 4'b0010, DOWN = 4'b0100, FILTER1 = 4'b1000;

  //参数定义
  reg [CNT_WIDTH-1:0] cnt;  //计数器
  reg [3:0] state;  //状态机当前状态
  reg key_d0, key_d1;  //按键输入的双寄存器同步

  //1.双寄存器同步，防止亚稳态
  always @(posedge clk) begin
    if (rst) begin
      key_d0 <= 1'b1;
      key_d1 <= 1'b1;
    end else begin
      key_d0 <= key_in;
      key_d1 <= key_d0;
    end
  end

  //2.边缘检测逻辑
  wire nedge = (key_d1 == 1'b1) && (key_d0 == 1'b0);  //下降沿检测
  wire pedge = (key_d1 == 1'b0) && (key_d0 == 1'b1);  //上升沿检测

  //3.状态机
  always @(posedge clk) begin
    if (rst) begin
      state <= IDLE;
      cnt <= {CNT_WIDTH{1'b0}};
      key_flag <= 1'b0;
      key_state <= 1'b1;  //按键未按下时为高电平
    end else begin
      case (state)
        IDLE: begin
          key_flag <= 1'b0;  //清除按键标志
          if (nedge) begin  //检测到按键按下边沿
            state <= FILTER0;
            cnt   <= {CNT_WIDTH{1'b0}};  //计数器清零
          end
        end
        FILTER0: begin
          if (cnt < CNT_MAX - 1) begin
            cnt <= cnt + 1'b1;  //计数器递增
            if (pedge) begin
              state <= IDLE;  //抖动，返回空闲状态
              cnt   <= {CNT_WIDTH{1'b0}};  //计数器清零
            end
          end else begin
            state <= DOWN;  //消抖完成，进入按下状态
            cnt <= {CNT_WIDTH{1'b0}};  //计数器清零
            key_state <= 1'b0;  //更新按键状态为按下
            //  key_flag <= 1'b1;  //产生按键有效标志脉冲位置1
          end
        end
        DOWN: begin
          if (pedge) begin  //检测到按键释放边沿
            state <= FILTER1;
            cnt   <= {CNT_WIDTH{1'b0}};  //计数器清零
          end
        end
        FILTER1: begin
          if (cnt < CNT_MAX - 1) begin
            cnt <= cnt + 1'b1;  //计数器递增
            if (nedge) begin
              state <= DOWN;  //抖动，返回按下状态
              cnt   <= {CNT_WIDTH{1'b0}};  //计数器清零
            end
          end else begin
            state <= IDLE;  //消抖完成，返回空闲状态
            key_state <= 1'b1;  //更新按键状态为未按下
            key_flag <= 1'b1;  //产生按键有效标志脉冲位置2
          end
        end

        default: state <= IDLE;  //默认回到空闲状态
      endcase
    end
  end
endmodule
