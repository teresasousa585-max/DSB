import math

# ==========================================
# 1. 核心参数配置 (16-bit 架构)
# ==========================================
ADC_MIN = 0
ADC_MAX = 65535
ADC_CENTER = 32768      # 16-bit 的静音中点
MODULATION_INDEX = 0.7 # 调制深度 70%
PWM_MAX = 1250          # 载波 2500，半周期最大占空比 1250 (对应 50% 占空比)

coe_filename = "am_arcsin_sqrt_mapping_64k.coe"
pwm_data = []

# 理想平方根包络的最大值 (对应 audio_norm = 1.0)
max_envelope = math.sqrt(1.0 + MODULATION_INDEX)

print("正在计算 16-bit ArcSin-SQRT 完美预畸变映射矩阵...")

for adc_val in range(ADC_MIN, ADC_MAX + 1):
    # 归一化：将 0~65535 映射为 -1.0 到 +1.0
    audio_norm = (adc_val - ADC_CENTER) / float(ADC_CENTER)
    
    if audio_norm < -1.0: audio_norm = -1.0
    if audio_norm > 1.0:  audio_norm = 1.0
    
    # ----------------------------------------------------
    # 核心算法修正：ArcSin + Sqrt 联合预畸变
    # ----------------------------------------------------
    # 1. 计算理想的平方根目标包络
    target_envelope = math.sqrt(1.0 + MODULATION_INDEX * audio_norm)
    
    # 2. 将包络归一化到 0.0 ~ 1.0 之间 (这就是我们期望的声压幅值)
    norm_amplitude = target_envelope / max_envelope
    
    # 3. 反三角函数预畸变 (核心救命代码！)
    # 因为 PWM 物理幅值 A = sin(pi * D)，我们需要占空比比例 D_ratio = arcsin(A) / pi
    # 由于 PWM_MAX=1250 对应的是半周期 (D_ratio = 0.5)，最大幅值对应 arcsin(1.0) = pi/2
    # 所以： pwm_val = (arcsin(norm_amplitude) / (pi / 2)) * PWM_MAX
    pwm_val = int(round((math.asin(norm_amplitude) / (math.pi / 2)) * PWM_MAX))
    
    # 限幅保护
    if pwm_val > PWM_MAX: pwm_val = PWM_MAX
    if pwm_val < 0:       pwm_val = 0
        
    pwm_data.append(pwm_val)

print(f"正在生成 Vivado 配置文件: {coe_filename} ...")
with open(coe_filename, 'w') as f:
    f.write('memory_initialization_radix=10;\n')
    f.write('memory_initialization_vector=\n')
    for i, val in enumerate(pwm_data):
        if i == len(pwm_data) - 1:
            f.write(f'{val};')
        else:
            f.write(f'{val},\n')

print("✅ 生成成功！请在 Vivado 中将 ROM 数据替换为此文件。你的超声波音响即将解锁真正的 HIFI 级音质！")