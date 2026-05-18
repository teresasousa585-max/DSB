`timescale 1ns / 1ps
// 语言: SystemVerilog

module uart_protocol_parser (
    input  logic        clk,
    input  logic        rst_n,
    
    // 来自 UART RX 底层模块的信号
    input  logic [7:0]  rx_data,
    input  logic        rx_done,
    
    // 输出给 32 路 PWM 发生器的参数
    output logic [7:0]  o_amplitude [0:31],
    output logic [11:0] o_phase     [0:31],
    output logic        o_update_pulse // 当一组数据成功更新时，拉高一个时钟周期
);

    // 状态机定义
    typedef enum logic [3:0] {
        S_IDLE,       // 0: 等待包头1 (0xAA)
        S_HEAD2,      // 1: 等待包头2 (0xBB)
        S_CMD,        // 2: 接收指令码 (0x01)
        S_AMP,        // 3: 接收32个幅度字节
        S_PHASE_L,    // 4: 接收相位低字节 (Little-Endian)
        S_PHASE_H,    // 5: 接收相位高字节
        S_CHECKSUM,   // 6: 接收并比对校验和
        S_TAIL1,      // 7: 等待包尾1 (0x0D)
        S_TAIL2       // 8: 等待包尾2 (0x0A)
    } state_t;

    state_t state;

    // 内部寄存器
    logic [5:0] byte_cnt;           // 计数器 (0~31)
    logic [7:0] calc_sum;           // 本地计算的校验和
    logic [7:0] temp_phase_L;       // 暂存相位的低字节
    
    // 影子寄存器 (Shadow Registers) - 只有校验通过才会覆盖输出
    logic [7:0]  shadow_amplitude [0:31];
    logic [11:0] shadow_phase     [0:31];

    // 状态机逻辑
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= S_IDLE;
            byte_cnt <= '0;
            calc_sum <= '0;
            o_update_pulse <= 1'b0;
            
            // 复位时默认参数：幅度满载，相位为0 (指向正前方)
            for (int i = 0; i < 32; i++) begin
                o_amplitude[i] <= 8'd255;
                o_phase[i]     <= 12'd0;
            end
        end 
        else begin
            o_update_pulse <= 1'b0; // 默认拉低

            if (rx_done) begin
                case (state)
                    // ----------------------------------------
                    S_IDLE: begin
                        if (rx_data == 8'hAA) state <= S_HEAD2;
                    end
                    
                    // ----------------------------------------
                    S_HEAD2: begin
                        if (rx_data == 8'hBB) state <= S_CMD;
                        else state <= S_IDLE; // 假包头，退回
                    end
                    
                    // ----------------------------------------
                    S_CMD: begin
                        if (rx_data == 8'h01) begin
                            calc_sum <= rx_data; // 校验和从指令码开始累加
                            byte_cnt <= '0;
                            state    <= S_AMP;
                        end else begin
                            state <= S_IDLE; // 不认识的指令
                        end
                    end
                    
                    // ----------------------------------------
                    S_AMP: begin
                        calc_sum <= calc_sum + rx_data; // 累加校验和
                        shadow_amplitude[byte_cnt] <= rx_data;
                        
                        if (byte_cnt == 6'd31) begin
                            byte_cnt <= '0;
                            state    <= S_PHASE_L;
                        end else begin
                            byte_cnt <= byte_cnt + 1'b1;
                        end
                    end
                    
                    // ----------------------------------------
                    S_PHASE_L: begin
                        calc_sum <= calc_sum + rx_data;
                        temp_phase_L <= rx_data; // 暂存低字节
                        state <= S_PHASE_H;
                    end
                    
                    // ----------------------------------------
                    S_PHASE_H: begin
                        calc_sum <= calc_sum + rx_data;
                        // 拼合 12 位相位: {高位低4位, 低字节8位}
                        shadow_phase[byte_cnt] <= {rx_data[3:0], temp_phase_L};
                        
                        if (byte_cnt == 6'd31) begin
                            state <= S_CHECKSUM;
                        end else begin
                            byte_cnt <= byte_cnt + 1'b1;
                            state    <= S_PHASE_L;
                        end
                    end
                    
                    // ----------------------------------------
                    S_CHECKSUM: begin
                        if (rx_data == calc_sum) begin
                            state <= S_TAIL1; // 校验正确
                        end else begin
                            state <= S_IDLE;  // 校验失败，丢弃这包数据
                        end
                    end
                    
                    // ----------------------------------------
                    S_TAIL1: begin
                        if (rx_data == 8'h0D) state <= S_TAIL2;
                        else state <= S_IDLE;
                    end
                    
                    // ----------------------------------------
                    S_TAIL2: begin
                        if (rx_data == 8'h0A) begin
                            // 校验和包尾全部正确，将影子寄存器安全地倒入输出寄存器
                            for (int i = 0; i < 32; i++) begin
                                o_amplitude[i] <= shadow_amplitude[i];
                                o_phase[i]     <= shadow_phase[i];
                            end
                            o_update_pulse <= 1'b1; // 触发更新脉冲
                        end
                        state <= S_IDLE; // 无论成败，回到空闲状态准备下一帧
                    end
                    
                    // ----------------------------------------
                    default: state <= S_IDLE;
                endcase
            end
        end
    end
endmodule