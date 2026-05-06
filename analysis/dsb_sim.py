# -*- coding: utf-8 -*-
"""
DSB信号与频谱仿真器
===================
FPGA超声阵列DSB定向音频系统的时域/频域/PWM/包络检波仿真。

输出:
    analysis_outputs/dsb_waveform.png
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import signal

from analysis.common import FC, FM, C, LAMBDA

# ============================================================================
# 仿真参数
# ============================================================================
FS = 10 * FC              # 采样率 400 kHz
T_DURATION = 2e-3         # 仿真时长 2 ms (2个1kHz周期)
M = 0.8                   # 调制指数

# PWM 参数
PWM_FS = 2_000_000        # 2 MHz，用于PWM细节仿真
PWM_FC = FC               # 40 kHz PWM载频

# 包络检波器参数
RC = 1.0 / (2 * np.pi * FM * 5)  # 低通RC时间常数，截止频率约5*fm


def ensure_output_dir(path: str) -> str:
    """确保输出目录存在，返回完整路径。"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def generate_dsb(fs: float, duration: float, m: float = M):
    """
    生成DSB调制信号。

    返回
    ----
    t : ndarray
        时间轴
    envelope : ndarray
        包络 [1 + m*cos(2*pi*fm*t)]
    carrier : ndarray
        载波 cos(2*pi*fc*t)
    s : ndarray
        DSB调制信号 envelope * carrier
    """
    t = np.arange(0, duration, 1 / fs)
    envelope = 1.0 + m * np.cos(2 * np.pi * FM * t)
    carrier = np.cos(2 * np.pi * FC * t)
    s = envelope * carrier
    return t, envelope, carrier, s


def compute_fft(s: np.ndarray, fs: float):
    """
    对信号做FFT，加Hann窗，返回单边谱频率与dB幅度。

    返回
    ----
    freqs : ndarray
        0 ~ fs/2 的频率轴
    mag_db : ndarray
        dB幅度 (以信号最大值为0 dB参考)
    """
    N = len(s)
    window = np.hanning(N)
    s_win = s * window
    S = np.fft.rfft(s_win)
    freqs = np.fft.rfftfreq(N, 1 / fs)
    # 幅度校正 (Hann窗能量补偿 ~1.63)
    mag = np.abs(S) * (2 / N) * 1.63
    # 避免log(0)
    mag_db = 20 * np.log10(mag + 1e-12)
    # 以信号最大值为0 dB参考
    mag_db -= np.max(mag_db)
    return freqs, mag_db


def compute_sideband_suppression(mag_db: np.ndarray, freqs: np.ndarray):
    """
    计算边带抑制比（载波与边带的dB差）。

    返回
    ----
    carrier_db : float
        载波(40kHz)幅度dB
    sideband_db : float
        边带(39kHz/41kHz)幅度dB
    suppression_db : float
        抑制比 (carrier_db - sideband_db)，应为负值表示边带低于载波
    """
    # 找40kHz附近峰值
    idx_c = np.argmin(np.abs(freqs - FC))
    carrier_db = mag_db[idx_c]

    # 找39kHz和41kHz附近峰值，取较高者
    idx_l = np.argmin(np.abs(freqs - (FC - FM)))
    idx_r = np.argmin(np.abs(freqs - (FC + FM)))
    sideband_db = max(mag_db[idx_l], mag_db[idx_r])

    suppression_db = carrier_db - sideband_db
    return carrier_db, sideband_db, suppression_db


def simulate_pwm(envelope: np.ndarray, t: np.ndarray, pwm_fs: float, pwm_fc: float):
    """
    模拟PWM占空比调制。

    占空比正比于瞬时包络（归一化到0~1）。
    返回PWM波形（0/1）以及占空比序列。
    """
    # 归一化包络到0~1作为占空比
    env_min = np.min(envelope)
    env_max = np.max(envelope)
    duty = (envelope - env_min) / (env_max - env_min + 1e-12)

    # 在更细的采样率下生成PWM
    # 对每个原始采样点，在PWM周期内生成高电平占空比*duty的方波
    samples_per_period = int(pwm_fs / pwm_fc)
    pwm_t = np.arange(0, t[-1] + 1 / pwm_fs, 1 / pwm_fs)
    pwm_wave = np.zeros_like(pwm_t)

    # 线性插值duty到pwm_t
    duty_interp = np.interp(pwm_t, t, duty)

    for i, d in enumerate(duty_interp):
        phase_in_period = (i % samples_per_period) / samples_per_period
        pwm_wave[i] = 1.0 if phase_in_period < d else 0.0

    return pwm_t, pwm_wave, duty


def envelope_detector(s: np.ndarray, fs: float, rc: float):
    """
    简单二极管包络检波器：全波整流 + 一阶RC低通滤波。

    返回
    ----
    rectified : ndarray
        全波整流信号
    recovered : ndarray
        恢复包络
    """
    # 全波整流
    rectified = np.abs(s)

    # 一阶IIR低通滤波器：y[n] = alpha * x[n] + (1-alpha) * y[n-1]
    # alpha = dt / (RC + dt)
    dt = 1.0 / fs
    alpha = dt / (rc + dt)
    recovered = np.zeros_like(rectified)
    recovered[0] = rectified[0]
    for n in range(1, len(rectified)):
        recovered[n] = alpha * rectified[n] + (1 - alpha) * recovered[n - 1]

    return rectified, recovered


def compute_thd(reference: np.ndarray, recovered: np.ndarray):
    """
    计算恢复包络相对于原始包络的失真度（归一化均方误差）。

    返回
    ----
    thd_db : float
        失真度 dB
    thd_percent : float
        失真度百分比
    """
    # 去除直流偏移后比较交流成分
    ref_ac = reference - np.mean(reference)
    rec_ac = recovered - np.mean(recovered)

    # 对齐增益
    gain = np.sum(ref_ac * rec_ac) / (np.sum(rec_ac ** 2) + 1e-12)
    rec_aligned = rec_ac * gain + np.mean(reference)

    error = reference - rec_aligned
    thd = np.sqrt(np.mean(error ** 2)) / (np.sqrt(np.mean(ref_ac ** 2)) + 1e-12)
    thd_percent = thd * 100.0
    thd_db = 20 * np.log10(thd + 1e-12)
    return thd_db, thd_percent


def run():
    # ========================================================================
    # 1. 时域波形仿真
    # ========================================================================
    t, envelope, carrier, s = generate_dsb(FS, T_DURATION, M)

    # ========================================================================
    # 2. FFT频谱
    # ========================================================================
    freqs, mag_db = compute_fft(s, FS)
    carrier_db, sideband_db, suppression_db = compute_sideband_suppression(mag_db, freqs)

    # ========================================================================
    # 3. PWM占空比调制
    # ========================================================================
    pwm_t, pwm_wave, duty = simulate_pwm(envelope, t, PWM_FS, PWM_FC)

    # ========================================================================
    # 4. 包络检波
    # ========================================================================
    rectified, recovered = envelope_detector(s, FS, RC)
    thd_db, thd_percent = compute_thd(envelope, recovered)

    # ========================================================================
    # 绘图
    # ========================================================================
    fig = plt.figure(figsize=(14, 12))

    # ------------------------------------------------------------------
    # 子图1: 包络
    # ------------------------------------------------------------------
    ax1 = fig.add_subplot(4, 1, 1)
    ax1.plot(t * 1e3, envelope, color="C0", linewidth=1.5)
    ax1.set_title("包络 envelope(t) = 1 + m·cos(2π·fm·t)")
    ax1.set_xlabel("时间 (ms)")
    ax1.set_ylabel("幅度")
    ax1.set_xlim(0, T_DURATION * 1e3)
    ax1.grid(True, linestyle="--", alpha=0.6)
    ax1.axhline(1.0, color="gray", linestyle=":", linewidth=0.8)

    # ------------------------------------------------------------------
    # 子图2: 载波局部放大 (0 ~ 0.1 ms)
    # ------------------------------------------------------------------
    ax2 = fig.add_subplot(4, 1, 2)
    zoom_end = int(0.1e-3 * FS)  # 0.1 ms 对应的采样点数
    ax2.plot(t[:zoom_end] * 1e3, carrier[:zoom_end], color="C1", linewidth=1.0)
    ax2.set_title(f"载波 carrier(t) = cos(2π·fc·t)  (局部放大 0~0.1 ms, fc={FC/1e3:.0f} kHz)")
    ax2.set_xlabel("时间 (ms)")
    ax2.set_ylabel("幅度")
    ax2.set_xlim(0, 0.1)
    ax2.grid(True, linestyle="--", alpha=0.6)

    # ------------------------------------------------------------------
    # 子图3: DSB调制输出
    # ------------------------------------------------------------------
    ax3 = fig.add_subplot(4, 1, 3)
    ax3.plot(t * 1e3, s, color="C2", linewidth=0.8, label="DSB信号 s(t)")
    ax3.plot(t * 1e3, envelope, color="C0", linewidth=1.0, linestyle="--", label="上包络")
    ax3.plot(t * 1e3, -envelope, color="C0", linewidth=1.0, linestyle="--", label="下包络")
    ax3.set_title("DSB调制输出 s(t) = envelope(t) · carrier(t)")
    ax3.set_xlabel("时间 (ms)")
    ax3.set_ylabel("幅度")
    ax3.set_xlim(0, T_DURATION * 1e3)
    ax3.grid(True, linestyle="--", alpha=0.6)
    ax3.legend(loc="upper right", fontsize=8)

    # ------------------------------------------------------------------
    # 子图4: FFT频谱
    # ------------------------------------------------------------------
    ax4 = fig.add_subplot(4, 1, 4)
    # 只显示 0 ~ 100 kHz
    f_max = 100_000
    idx_fmax = np.argmin(np.abs(freqs - f_max))
    ax4.plot(freqs[:idx_fmax] / 1e3, mag_db[:idx_fmax], color="C3", linewidth=1.0)

    # 标注峰值
    peak_freqs = [FC - FM, FC, FC + FM]
    peak_labels = [f"{ (FC - FM)/1e3:.0f}kHz", f"{FC/1e3:.0f}kHz", f"{(FC + FM)/1e3:.0f}kHz"]
    for pf, pl in zip(peak_freqs, peak_labels):
        idx = np.argmin(np.abs(freqs - pf))
        ax4.axvline(freqs[idx] / 1e3, color="gray", linestyle=":", linewidth=0.8)
        ax4.annotate(
            f"{pl}\n{mag_db[idx]:.1f} dB",
            xy=(freqs[idx] / 1e3, mag_db[idx]),
            xytext=(5, 10),
            textcoords="offset points",
            fontsize=7,
            arrowprops=dict(arrowstyle="->", color="black", lw=0.5),
        )

    ax4.set_title(
        f"FFT频谱 (Hann窗) | 边带抑制比: {suppression_db:.1f} dB | "
        f"THD(包络恢复): {thd_percent:.2f}% ({thd_db:.1f} dB)"
    )
    ax4.set_xlabel("频率 (kHz)")
    ax4.set_ylabel("幅度 (dB)")
    ax4.set_xlim(0, f_max / 1e3)
    ax4.set_ylim(-80, 5)
    ax4.grid(True, linestyle="--", alpha=0.6)

    plt.tight_layout()
    out_path = ensure_output_dir("analysis_outputs/dsb_waveform.png")
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    # ========================================================================
    # PWM占空比曲线 (单独子图，附加到同一张图或另存；这里另存一张详细图)
    # ========================================================================
    fig2, axes2 = plt.subplots(2, 1, figsize=(12, 6))

    # 一个1kHz周期内的占空比变化
    ax_pwm1 = axes2[0]
    one_period_samples = int(FS / FM)
    t_one = t[:one_period_samples] * 1e3
    duty_one = duty[:one_period_samples]
    ax_pwm1.plot(t_one, duty_one, color="C4", linewidth=1.5)
    ax_pwm1.set_title("PWM占空比调制 (一个1kHz调制周期内)")
    ax_pwm1.set_xlabel("时间 (ms)")
    ax_pwm1.set_ylabel("占空比")
    ax_pwm1.set_xlim(0, 1.0)
    ax_pwm1.grid(True, linestyle="--", alpha=0.6)

    # PWM频谱（边带围绕40kHz整数倍）
    ax_pwm2 = axes2[1]
    pwm_freqs, pwm_mag_db = compute_fft(pwm_wave, PWM_FS)
    # 显示 0 ~ 200 kHz
    f_max_pwm = 200_000
    idx_fmax_pwm = np.argmin(np.abs(pwm_freqs - f_max_pwm))
    ax_pwm2.plot(pwm_freqs[:idx_fmax_pwm] / 1e3, pwm_mag_db[:idx_fmax_pwm], color="C5", linewidth=0.8)
    # 标注40kHz及其倍频
    for n in range(1, 6):
        f_n = n * PWM_FC
        if f_n > f_max_pwm:
            break
        idx_n = np.argmin(np.abs(pwm_freqs - f_n))
        ax_pwm2.axvline(pwm_freqs[idx_n] / 1e3, color="gray", linestyle=":", linewidth=0.6)
        ax_pwm2.annotate(
            f"{f_n/1e3:.0f}kHz",
            xy=(pwm_freqs[idx_n] / 1e3, pwm_mag_db[idx_n]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=6,
        )
    ax_pwm2.set_title("PWM频谱 (边带围绕40kHz整数倍)")
    ax_pwm2.set_xlabel("频率 (kHz)")
    ax_pwm2.set_ylabel("幅度 (dB)")
    ax_pwm2.set_xlim(0, f_max_pwm / 1e3)
    ax_pwm2.set_ylim(-100, 5)
    ax_pwm2.grid(True, linestyle="--", alpha=0.6)

    plt.tight_layout()
    out_path2 = ensure_output_dir("analysis_outputs/dsb_pwm_detail.png")
    fig2.savefig(out_path2, dpi=200, bbox_inches="tight")
    plt.close(fig2)

    # ========================================================================
    # 包络检波对比图
    # ========================================================================
    fig3, ax3 = plt.subplots(figsize=(12, 4))
    ax3.plot(t * 1e3, envelope, color="C0", linewidth=1.5, label="原始包络")
    ax3.plot(t * 1e3, recovered, color="C6", linewidth=1.5, linestyle="--", label="恢复包络(检波后)")
    ax3.set_title(f"包络检波对比 | THD = {thd_percent:.2f}% ({thd_db:.1f} dB)")
    ax3.set_xlabel("时间 (ms)")
    ax3.set_ylabel("幅度")
    ax3.set_xlim(0, T_DURATION * 1e3)
    ax3.grid(True, linestyle="--", alpha=0.6)
    ax3.legend(loc="upper right")
    plt.tight_layout()
    out_path3 = ensure_output_dir("analysis_outputs/dsb_envelope_detection.png")
    fig3.savefig(out_path3, dpi=200, bbox_inches="tight")
    plt.close(fig3)

    # ========================================================================
    # 控制台输出
    # ========================================================================
    print("=" * 60)
    print("DSB 信号与频谱仿真结果")
    print("=" * 60)
    print(f"调制指数 m          : {M}")
    print(f"载波频率 fc         : {FC} Hz")
    print(f"调制频率 fm         : {FM} Hz")
    print(f"声速 C              : {C} m/s")
    print(f"波长 λ              : {LAMBDA*1e3:.3f} mm")
    print(f"采样率 fs           : {FS/1e6:.2f} MHz")
    print(f"PWM采样率           : {PWM_FS/1e6:.2f} MHz")
    print("-" * 60)
    print(f"载波幅度(40kHz)     : {carrier_db:.2f} dB")
    print(f"边带幅度(39/41kHz)  : {sideband_db:.2f} dB")
    print(f"边带抑制比          : {suppression_db:.1f} dB")
    print(f"PWM占空比范围       : {duty.min()*100:.1f}% ~ {duty.max()*100:.1f}%")
    print(f"包络恢复 THD        : {thd_percent:.2f}% ({thd_db:.1f} dB)")
    print("-" * 60)
    print(f"输出图像:")
    print(f"  1) {out_path}")
    print(f"  2) {out_path2}")
    print(f"  3) {out_path3}")
    print("=" * 60)


if __name__ == "__main__":
    run()
