"""
波阵面与阵列方向图仿真

功能：
1. 2D 声压场分布（近场/远场）
2. 极坐标方向图（主瓣、旁瓣、栅瓣分析）
3. 加窗效果对比
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from analysis.common import (
    C, D, LAMBDA, K, D_OVER_LAMBDA, PHASE_STEPS, AMP_MAX,
    N_ELEMENTS, ELEMENT_X, ELEMENT_Y, ELEMENT_ROWS, ELEMENT_COLS,
    STEER_MIN, STEER_MAX, STEER_POINTS,
    FOCUS_MIN_MM, FOCUS_MAX_MM, FOCUS_STEP_MM, FOCUS_POINTS,
)

# ============================================================================
# 物理计算核心
# ============================================================================

def compute_pressure_field(
    x_grid: np.ndarray,
    z_grid: np.ndarray,
    steer_angle_deg: float = 0.0,
    focus_dist_mm: float = None,
    window_amp: np.ndarray = None,
) -> np.ndarray:
    """
    计算声压场 P(x, z)。

    参数：
        x_grid, z_grid: 2D meshgrid 坐标（单位：m）
        steer_angle_deg: 偏转角度（度）
        focus_dist_mm: 聚焦距离（mm），None 表示不聚焦
        window_amp: 25 路幅度增益 (0..1)，None 表示全 1

    返回：
        与 x_grid 同形的复数声压数组
    """
    if window_amp is None:
        window_amp = np.ones(N_ELEMENTS, dtype=np.float64)

    # 偏转相位（线性相位梯度）
    steer_rad = np.deg2rad(steer_angle_deg)
    steer_phase = PHASE_STEPS * D * (ELEMENT_COLS * np.sin(steer_rad)) / LAMBDA
    # 也可以考虑俯仰角，但这里只在 x-z 平面仿真，所以只考虑 x 方向相位梯度

    # 聚焦相位
    if focus_dist_mm is not None:
        r0 = focus_dist_mm / 1000.0
        ri = np.sqrt(ELEMENT_X**2 + ELEMENT_Y**2 + r0**2)
        focus_phase = (PHASE_STEPS * (ri - r0) / LAMBDA) % PHASE_STEPS
    else:
        focus_phase = np.zeros(N_ELEMENTS)

    total_phase = steer_phase + focus_phase
    total_phase_rad = 2 * np.pi * total_phase / PHASE_STEPS

    # 计算每个阵元到每个场点的距离
    # x_grid.shape = (Nz, Nx)
    P = np.zeros_like(x_grid, dtype=np.complex128)

    for i in range(N_ELEMENTS):
        dx = x_grid - ELEMENT_X[i]
        dz = z_grid - 0.0  # 阵元在 z=0 平面
        # 注意：假设阵元在 y=0 平面，我们只仿真 x-z 平面（y=0）
        # 严格来说距离应该是 sqrt(dx^2 + dy^2 + dz^2)，但 dy=0
        r = np.sqrt(dx**2 + dz**2)
        # 避免除零
        r = np.where(r < 1e-6, 1e-6, r)

        k_r = K * r
        phase_i = total_phase_rad[i]
        amp_i = window_amp[i]

        P += amp_i * np.exp(1j * (phase_i + k_r)) / r

    return P


def compute_beam_pattern(
    angles_deg: np.ndarray,
    r_m: float = 5.0,
    steer_angle_deg: float = 0.0,
    focus_dist_mm: float = None,
    window_amp: np.ndarray = None,
) -> np.ndarray:
    """
    计算远场方向图（标准阵列因子公式）。

    使用远场近似，避免点源精确模型在 d/λ>0.5 时引入的栅瓣失真。

    参数：
        angles_deg: 角度数组（度），在 x-z 平面内扫描
        r_m: 距离（m），仅用于兼容性，不影响远场方向图形状
        steer_angle_deg: 偏转角度
        focus_dist_mm: 聚焦距离（远场下聚焦效应可忽略）
        window_amp: 25 路幅度增益

    返回：
        归一化声压幅度（线性，非 dB）
    """
    angles_rad = np.deg2rad(angles_deg)

    if window_amp is None:
        window_amp = np.ones(N_ELEMENTS)

    # 偏转相位（弧度）
    steer_rad = np.deg2rad(steer_angle_deg)
    steer_phase = K * ELEMENT_X * np.sin(steer_rad)

    # 聚焦相位（远场可忽略，保留接口兼容性）
    focus_phase = np.zeros(N_ELEMENTS)

    total_phase = steer_phase + focus_phase

    # 远场阵列因子：
    # AF(theta) = sum_i A_i * exp(j * (phase_i - k * x_i * sin(theta)))
    # 其中 -k*x_i*sin(theta) 是传播到角度 theta 的波程差相位
    P = np.zeros(len(angles_deg), dtype=np.complex128)
    for i in range(N_ELEMENTS):
        phase_i = total_phase[i]
        amp_i = window_amp[i]
        x_i = ELEMENT_X[i]
        prop_phase = -K * x_i * np.sin(angles_rad)
        P += amp_i * np.exp(1j * (phase_i + prop_phase))

    # 归一化
    mag = np.abs(P)
    mag = mag / mag.max()
    return mag


# ============================================================================
# 可视化
# ============================================================================

def plot_2d_pressure(
    out_dir: Path,
    steer_angle_deg: float = 0.0,
    focus_dist_mm: float = None,
    window_amp: np.ndarray = None,
    tag: str = "",
):
    """绘制 2D 声压场分布。"""
    # 网格：x = -500~+500 mm, z = 10~1000 mm
    x_mm = np.linspace(-500, 500, 401)
    z_mm = np.linspace(10, 1000, 400)
    x_m, z_m = np.meshgrid(x_mm / 1000.0, z_mm / 1000.0)

    P = compute_pressure_field(x_m, z_m, steer_angle_deg, focus_dist_mm, window_amp)
    intensity = np.abs(P)**2
    intensity_db = 20 * np.log10(np.abs(P) + 1e-12)
    intensity_db = np.clip(intensity_db, -40, 0)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # 线性强度
    ax = axes[0]
    im = ax.imshow(
        intensity,
        extent=[x_mm.min(), x_mm.max(), z_mm.min(), z_mm.max()],
        origin="lower",
        aspect="auto",
        cmap="hot",
    )
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title("声压强度 |P|² (线性)")
    plt.colorbar(im, ax=ax)

    # dB 强度
    ax = axes[1]
    im = ax.imshow(
        intensity_db,
        extent=[x_mm.min(), x_mm.max(), z_mm.min(), z_mm.max()],
        origin="lower",
        aspect="auto",
        cmap="viridis",
        vmin=-40,
        vmax=0,
    )
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title("声压强度 (dB)")
    plt.colorbar(im, ax=ax, label="dB")

    title_parts = []
    if steer_angle_deg != 0:
        title_parts.append(f"偏转 {steer_angle_deg}°")
    if focus_dist_mm is not None:
        title_parts.append(f"聚焦 {focus_dist_mm}mm")
    if tag:
        title_parts.append(tag)
    fig.suptitle(", ".join(title_parts) if title_parts else "同相驱动")

    fig.tight_layout()
    fname = out_dir / f"array_2d_{tag or 'default'}.png"
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"  已保存: {fname}")


def plot_polar_pattern(
    out_dir: Path,
    steer_angle_deg: float = 0.0,
    focus_dist_mm: float = None,
    tag: str = "",
):
    """绘制极坐标方向图，对比不同窗函数。"""
    angles = np.linspace(-90, 90, 1801)
    r_m = 5.0  # 远场 5m

    # 无窗（矩形窗）
    mag_rect = compute_beam_pattern(angles, r_m, steer_angle_deg, focus_dist_mm, window_amp=np.ones(N_ELEMENTS))

    # Hann
    from scipy.signal.windows import hann
    w_hann = np.outer(hann(5), hann(5)).flatten()
    w_hann = w_hann / w_hann.max()
    mag_hann = compute_beam_pattern(angles, r_m, steer_angle_deg, focus_dist_mm, window_amp=w_hann)

    # Hamming
    from scipy.signal.windows import hamming
    w_hamm = np.outer(hamming(5), hamming(5)).flatten()
    w_hamm = w_hamm / w_hamm.max()
    mag_hamm = compute_beam_pattern(angles, r_m, steer_angle_deg, focus_dist_mm, window_amp=w_hamm)

    # Blackman
    from scipy.signal.windows import blackman
    w_blk = np.outer(blackman(5), blackman(5)).flatten()
    w_blk = w_blk / w_blk.max()
    mag_blk = compute_beam_pattern(angles, r_m, steer_angle_deg, focus_dist_mm, window_amp=w_blk)

    # 转换为 dB
    mag_rect_db = 20 * np.log10(mag_rect + 1e-12)
    mag_hann_db = 20 * np.log10(mag_hann + 1e-12)
    mag_hamm_db = 20 * np.log10(mag_hamm + 1e-12)
    mag_blk_db = 20 * np.log10(mag_blk + 1e-12)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), subplot_kw={"projection": "polar"})

    for ax, data, title in [
        (axes[0], [(mag_rect_db, "矩形", "C0"), (mag_hann_db, "Hann", "C1"), (mag_hamm_db, "Hamming", "C2"), (mag_blk_db, "Blackman", "C3")], "方向图对比（dB）"),
    ]:
        for mag_db, label, color in data:
            ax.plot(np.deg2rad(angles), mag_db, label=label, color=color, linewidth=0.8)
        ax.set_ylim(-40, 0)
        ax.set_yticks([-30, -20, -10, 0])
        ax.set_title(title)
        ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1))
        ax.grid(True)

    # 第二个极坐标图：放大主瓣区域
    ax = axes[1]
    zoom_mask = (angles >= steer_angle_deg - 20) & (angles <= steer_angle_deg + 20)
    angles_zoom = angles[zoom_mask]
    for mag_db, label, color in [
        (mag_rect_db[zoom_mask], "矩形", "C0"),
        (mag_hann_db[zoom_mask], "Hann", "C1"),
        (mag_hamm_db[zoom_mask], "Hamming", "C2"),
        (mag_blk_db[zoom_mask], "Blackman", "C3"),
    ]:
        ax.plot(np.deg2rad(angles_zoom), mag_db, label=label, color=color, linewidth=1.2)
    ax.set_ylim(-40, 0)
    ax.set_title("主瓣细节（±20°）")
    ax.grid(True)

    # 标记理论栅瓣位置
    theta_g = np.rad2deg(np.arcsin(LAMBDA / D))  # ~32.5°
    for ax in axes:
        ax.axvline(np.deg2rad(theta_g), color="red", linestyle="--", alpha=0.5, linewidth=0.8)
        ax.axvline(np.deg2rad(-theta_g), color="red", linestyle="--", alpha=0.5, linewidth=0.8)

    fig.suptitle(f"5×5 阵列方向图 (d/λ={D_OVER_LAMBDA:.2f}, 偏转={steer_angle_deg}°)")
    fig.tight_layout()
    fname = out_dir / f"array_polar_{tag or 'default'}.png"
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"  已保存: {fname}")

    # 打印指标
    print("\n  方向图指标（矩形窗）：")
    peak_idx = np.argmax(mag_rect)
    peak_angle = angles[peak_idx]
    print(f"    峰值角度: {peak_angle:.2f}°")

    # -3 dB 波束宽度：在主瓣附近搜索，避免跨越到栅瓣区域
    # 方法：从主峰向两侧搜索，找到 mag_rect_db 首次低于 -3 dB 的角度
    left_idx = peak_idx
    while left_idx > 0 and mag_rect_db[left_idx] >= -3:
        left_idx -= 1
    right_idx = peak_idx
    while right_idx < len(angles) - 1 and mag_rect_db[right_idx] >= -3:
        right_idx += 1
    bw = angles[right_idx] - angles[left_idx]
    print(f"    -3 dB 波束宽度: {bw:.2f}°")

    # 最大旁瓣电平（排除主瓣附近 ±5°）
    mainlobe_mask = np.abs(angles - peak_angle) < 5
    sll_db = mag_rect_db[~mainlobe_mask].max()
    print(f"    最大旁瓣电平 (SLL): {sll_db:.2f} dB")
    print(f"    理论第一栅瓣位置: ±{theta_g:.2f}°")


def plot_window_comparison(out_dir: Path):
    """绘制加窗 vs 不加窗的 2D 场对比。"""
    focus_mm = 300
    tag_rect = "focus300_rect"
    tag_hann = "focus300_hann"

    from scipy.signal.windows import hann
    w_hann = np.outer(hann(5), hann(5)).flatten()
    w_hann = w_hann / w_hann.max()

    plot_2d_pressure(out_dir, steer_angle_deg=0, focus_dist_mm=focus_mm, window_amp=np.ones(N_ELEMENTS), tag=tag_rect)
    plot_2d_pressure(out_dir, steer_angle_deg=0, focus_dist_mm=focus_mm, window_amp=w_hann, tag=tag_hann)


# ============================================================================
# 主入口
# ============================================================================

def run():
    print("=" * 70)
    print("波阵面与阵列方向图仿真")
    print("=" * 70)

    out_dir = Path(__file__).resolve().parent.parent / "analysis_outputs"
    out_dir.mkdir(exist_ok=True)

    # 1. 同相驱动 2D 场
    print("\n[1/5] 同相驱动声压场 ...")
    plot_2d_pressure(out_dir, steer_angle_deg=0, focus_dist_mm=None, tag="broadside")

    # 2. 偏转 15° 的 2D 场
    print("\n[2/5] 偏转 15° 声压场 ...")
    plot_2d_pressure(out_dir, steer_angle_deg=15, focus_dist_mm=None, tag="steer15")

    # 3. 聚焦 300mm 的 2D 场
    print("\n[3/5] 聚焦 300mm 声压场 ...")
    plot_2d_pressure(out_dir, steer_angle_deg=0, focus_dist_mm=300, tag="focus300")

    # 4. 方向图对比
    print("\n[4/5] 方向图对比 ...")
    plot_polar_pattern(out_dir, steer_angle_deg=0, focus_dist_mm=None, tag="broadside")
    plot_polar_pattern(out_dir, steer_angle_deg=15, focus_dist_mm=None, tag="steer15")

    # 5. 加窗对比
    print("\n[5/5] 加窗效果对比 ...")
    plot_window_comparison(out_dir)

    print("\n" + "=" * 70)
    print("所有阵列仿真图已保存到 analysis_outputs/")
    print("=" * 70)


if __name__ == "__main__":
    run()
