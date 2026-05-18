`timescale 1ns / 1ps
// 语言: SystemVerilog

module pwm32_generator (
    input  logic               clk,
    input  logic               rst_n,
    input  logic signed [15:0] audio_in,         
    input  logic [7:0]         amplitude [0:31], // 扩展为 32
    input  logic [11:0]        phase_del [0:31], // 扩展为 32
    output logic [31:0]        pwm_out           // 扩展为 32
);

    // ================== 1. 载波生成器 ==================
    // 40kHz 载波计数器 (100MHz / 2500 = 40kHz)，计数范围 0~2499
    logic [11:0] carrier_cnt;
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) carrier_cnt <= 12'd0;
        else if (carrier_cnt == 12'd2499) carrier_cnt <= 12'd0;
        else carrier_cnt <= carrier_cnt + 1'b1;
    end

    // ================== 2. AM-SQRT 预畸变查表 (16-bit 深度) ==================
    // 将有符号音频 (-32768 ~ +32767) 抬高，转为无符号 ROM 地址 (0 ~ 65535)
    logic [15:0] rom_addr;
    always_ff @(posedge clk) begin
        rom_addr <= audio_in + 16'd32768; 
    end

    logic [11:0] duty_base_d2; // 从 ROM 读出的开完根号的基础占空比 (0~1250)
    
    // Vivado中需将此 IP 的深度设置为 65536，宽度为 12
    sqrt_mapping_rom u_math_rom (
        .clka(clk),
        .addra(rom_addr),
        .ena(1'b1),
        .douta(duty_base_d2)
    );

    // ================== 3. DSP 分发流水线 (32路并发) ==================
    genvar i;
    generate
        for (i = 0; i < 32; i++) begin : gen_ch
            
            // Stage 3: 通道加权 (波束赋形)
            logic [19:0] duty_ch_mult;
            assign duty_ch_mult = duty_base_d2 * {12'h0, amplitude[i]};
            
            logic [11:0] duty_ch_d3;
            logic [12:0] start_raw_d3; // 13位防溢出保护
            always_ff @(posedge clk) begin
                duty_ch_d3   <= duty_ch_mult[19:8]; // 右移 8 位 (等效于除以 256 归一化)
                
                // AM 模式下直接使用相位延迟
                start_raw_d3 <= {1'b0, phase_del[i]};
            end

            // Stage 4: 计算高电平落点并处理环形回绕
            logic [11:0] duty_ch_d4;
            logic [11:0] start_d4;
            logic [12:0] end_raw_d4;
            always_ff @(posedge clk) begin
                duty_ch_d4 <= duty_ch_d3;
                // 取模 2500，确保起点合法
                start_d4   <= (start_raw_d3 >= 13'd2500) ? (start_raw_d3 - 13'd2500) : start_raw_d3[11:0];
                // 终点 = 起点 + 占空比宽度
                end_raw_d4 <= ((start_raw_d3 >= 13'd2500) ? (start_raw_d3 - 13'd2500) : start_raw_d3[11:0]) + duty_ch_d3;
            end

            // Stage 5: 载波边界快门锁
            // 只有当载波计数器达到周期末尾(2499)时，才更新下一周期的占空比和起止点
            logic [11:0] duty_locked;
            logic [11:0] start_locked;
            logic [11:0] end_locked;
            logic        wrap_locked;
            always_ff @(posedge clk) begin
                if (carrier_cnt == 12'd2499) begin
                    duty_locked  <= duty_ch_d4;
                    start_locked <= start_d4;
                    // 处理跨周期回绕 (例如起点 2400，宽度 200，则终点为 100)
                    end_locked   <= (end_raw_d4 >= 13'd2500) ? (end_raw_d4 - 13'd2500) : end_raw_d4[11:0];
                    wrap_locked  <= (end_raw_d4 >= 13'd2500); 
                end
            end

            // 终极零延时组合逻辑输出
            assign pwm_out[i] = (duty_locked == 0) ? 1'b0 :
                                (!wrap_locked)     ? (carrier_cnt >= start_locked && carrier_cnt < end_locked) :
                                                     (carrier_cnt >= start_locked || carrier_cnt < end_locked);
        end
    endgenerate
endmodule