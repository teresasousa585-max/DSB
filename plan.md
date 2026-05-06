# FPGA超声阵列定向音频系统 — Verilog代码框架生成计划

## Stage 1: 生成所有Verilog模块（并行）
- mcp3201_driver.v — MCP3201 SPI ADC驱动
- pwm25_generator.v — 25路独立相位/幅度PWM（并行输出直接驱动换能器）
- dds_generator.v — DDS正弦波（1kHz调制 + 40kHz载波参考）
- dsb_modulator.v — DSB调幅器
- beam_controller.v — 波束相位/幅度控制器

## Stage 2: 生成顶层模块 + Testbench
- top_ultrasound_array.v — 顶层整合
- tb_top.v — 完整测试平台

## Stage 3: 验证与交付
- 代码review、时序检查
- 输出文件到 /mnt/agents/output/verilog/
