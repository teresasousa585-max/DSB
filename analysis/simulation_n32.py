"""
N=32 同心圆环超声波阵列仿真代码
频率: 40kHz, 空气中
阵列结构: 中心1 + 三环(7/19.5mm + 9/37mm + 15/58mm) = 32个阵元
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j1

# ==================== 参数 ====================
c = 343.0          # 声速 (m/s)
f = 40000          # 频率 (Hz)
lam = c / f        # 波长 (m)
k = 2 * np.pi / lam  # 波数 (rad/m)
D = 16e-3          # 换能器直径 (m)
a = D / 2          # 换能器半径 (m)
MIN_DIST = 16.5e-3 # 最小中心距 (m)

# N=32 阵列参数
center = True
ring_radii = [19.5e-3, 37.0e-3, 58.0e-3]  # 环半径 (m)
ring_counts = [7, 9, 15]                     # 每环阵元数

# ==================== 函数 ====================

def generate_ring_array(center, radii, counts):
    """生成同心圆环阵元位置"""
    positions = []
    if center:
        positions.append([0.0, 0.0])
    for r, n in zip(radii, counts):
        for i in range(n):
            angle = 2 * np.pi * i / n
            positions.append([r * np.cos(angle), r * np.sin(angle)])
    return np.array(positions)

def element_pattern(theta, a, k_val):
    """圆形活塞换能器方向图"""
    t = np.where(np.abs(theta) < 1e-10, 1e-10, theta)
    x = k_val * a * np.sin(t)
    x = np.where(np.abs(x) < 1e-10, 1e-10, x)
    return np.abs(2 * j1(x) / x)

def array_pattern(positions, weights, theta, phi, k_val):
    """计算阵列方向图 (远场)"""
    AF = np.zeros_like(theta, dtype=complex)
    kx = k_val * np.sin(theta) * np.cos(phi)
    ky = k_val * np.sin(theta) * np.sin(phi)
    for i, (x, y) in enumerate(positions):
        AF += weights[i] * np.exp(1j * (kx * x + ky * y))
    return np.abs(AF)

def normalize_db(pattern):
    """归一化到dB"""
    m = np.max(np.abs(pattern))
    return 20 * np.log10(np.abs(pattern) / m + 1e-10) if m > 0 else np.zeros_like(pattern) - 100

def beamsteer_weights(positions, theta_steer, phi_steer, k_val):
    """计算波束偏转权重 (相位延迟)"""
    kx_s = k_val * np.sin(theta_steer) * np.cos(phi_steer)
    ky_s = k_val * np.sin(theta_steer) * np.sin(phi_steer)
    weights = np.ones(len(positions), dtype=complex)
    for i, (x, y) in enumerate(positions):
        phase = kx_s * x + ky_s * y
        weights[i] = np.exp(-1j * phase)
    return weights / np.sum(np.abs(weights))

def focus_weights(positions, x_focus, y_focus, z_focus, k_val):
    """计算聚焦权重"""
    weights = np.ones(len(positions), dtype=complex)
    r_ref = np.sqrt(x_focus**2 + y_focus**2 + z_focus**2)
    for i, (x, y) in enumerate(positions):
        r_focus = np.sqrt((x_focus - x)**2 + (y_focus - y)**2 + z_focus**2)
        phase = k_val * (r_focus - r_ref)
        weights[i] = np.exp(-1j * phase)
    return weights / np.sum(np.abs(weights))

def compute_pressure_field(positions, weights, x_grid, y_grid, z_dist, k_val, a_val):
    """计算近场声压分布"""
    P = np.zeros_like(x_grid, dtype=complex)
    for i, (x0, y0) in enumerate(positions):
        r = np.sqrt((x_grid - x0)**2 + (y_grid - y0)**2 + z_dist**2)
        theta_e = np.arctan2(np.sqrt((x_grid-x0)**2 + (y_grid-y0)**2), z_dist)
        ka = k_val * a_val
        ts = np.where(np.abs(theta_e) < 1e-10, 1e-10, theta_e)
        xe = ka * np.sin(ts)
        xe = np.where(np.abs(xe) < 1e-10, 1e-10, xe)
        elem = np.abs(2 * j1(xe) / xe)
        P += weights[i] * elem * np.exp(1j * k_val * r) / (r + 1e-10)
    return np.abs(P)

# ==================== 主程序 ====================
if __name__ == "__main__":
    # 生成阵列
    pos = generate_ring_array(True, ring_radii, ring_counts)
    N = len(pos)
    w_uniform = np.ones(N) / N

    print(f"N={N} 同心圆环阵列已生成")
    print(f"孔径: {2*max(ring_radii)*1000:.0f}mm")

    # 角度网格
    theta = np.linspace(0.01, np.pi/2, 901)
    theta_deg = theta * 180 / np.pi
    EP = element_pattern(theta, a, k)

    # 1. 远场方向图
    phi_angles = [0, np.pi/8, np.pi/4, np.pi/2]
    phi_labels = ["0°", "22.5°", "45°", "90°"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    for idx, (phi, label) in enumerate(zip(phi_angles, phi_labels)):
        AF = array_pattern(pos, w_uniform, theta, phi, k)
        db = normalize_db(AF * EP)
        psl = np.max(db[theta_deg > 5])
        axes[idx].plot(theta_deg, db, "b-", linewidth=2)
        axes[idx].set_title(f"φ = {label}, PSL = {psl:.1f} dB")
        axes[idx].set_xlim(0, 90)
        axes[idx].set_ylim(-50, 5)
        axes[idx].grid(True, alpha=0.3)
    plt.suptitle("N=32 Ring Array - Far-field Pattern")
    plt.tight_layout()
    plt.savefig("farfield_pattern.png", dpi=150)

    # 2. 波束偏转
    steer_angles = [0, 10, 20, 30]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    for idx, sa in enumerate(steer_angles):
        w_s = beamsteer_weights(pos, sa*np.pi/180, 0, k)
        AF = array_pattern(pos, w_s, theta, 0, k)
        db = normalize_db(AF * EP)
        axes[idx].plot(theta_deg, db, "b-", linewidth=2)
        axes[idx].axvline(x=sa, color="r", linestyle="--", alpha=0.5)
        axes[idx].set_title(f"Steered to {sa}°")
        axes[idx].set_xlim(0, 60)
        axes[idx].set_ylim(-35, 5)
        axes[idx].grid(True, alpha=0.3)
    plt.suptitle("Beam Steering Capability")
    plt.tight_layout()
    plt.savefig("beam_steering.png", dpi=150)

    # 3. 近场聚焦
    x_range, y_range = 100e-3, 100e-3
    x_grid = np.linspace(-x_range, x_range, 201)
    y_grid = np.linspace(-y_range, y_range, 201)
    X, Y = np.meshgrid(x_grid, y_grid)

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    for idx, z_f in enumerate([100e-3, 200e-3, 500e-3]):
        w_f = focus_weights(pos, 0, 0, z_f, k)
        P = compute_pressure_field(pos, w_f, X, Y, z_f, k, a)
        P_db = 20 * np.log10(P / np.max(P) + 1e-10)
        im = axes[0][idx].contourf(X*1000, Y*1000, P_db, levels=41, cmap="jet", extend="min")
        axes[0][idx].set_title(f"Focused @ z={z_f*1000:.0f}mm (XY)")
        axes[0][idx].set_aspect("equal")
        plt.colorbar(im, ax=axes[0][idx], shrink=0.6)
    for idx, z_f in enumerate([100e-3, 200e-3, 500e-3]):
        w_f = focus_weights(pos, 0, 0, z_f, k)
        x_xz = np.linspace(-80e-3, 80e-3, 161)
        z_xz = np.linspace(20e-3, 800e-3, 391)
        X_xz, Z_xz = np.meshgrid(x_xz, z_xz)
        Y_xz = np.zeros_like(X_xz)
        P_xz = np.zeros_like(X_xz, dtype=complex)
        for i, (x0, y0) in enumerate(pos):
            r_xz = np.sqrt((X_xz-x0)**2 + (Y_xz-y0)**2 + Z_xz**2)
            theta_xz = np.arctan2(np.sqrt((X_xz-x0)**2+(Y_xz-y0)**2), Z_xz)
            ka = k * a
            ts = np.where(np.abs(theta_xz)<1e-10, 1e-10, theta_xz)
            xe = ka*np.sin(ts)
            xe = np.where(np.abs(xe)<1e-10, 1e-10, xe)
            elem = np.abs(2*j1(xe)/xe)
            P_xz += w_f[i]*elem*np.exp(1j*k*r_xz)/(r_xz+1e-10)
        P_xz_db = 20*np.log10(np.abs(P_xz)/np.max(np.abs(P_xz))+1e-10)
        im = axes[1][idx].contourf(X_xz*1000, Z_xz*1000, P_xz_db, levels=41, cmap="jet", extend="min")
        axes[1][idx].set_title(f"Focused @ z={z_f*1000:.0f}mm (XZ)")
        axes[1][idx].axhline(y=z_f*1000, color="white", linestyle="--")
        plt.colorbar(im, ax=axes[1][idx], shrink=0.6)
    plt.suptitle("Focusing Capability")
    plt.tight_layout()
    plt.savefig("focusing.png", dpi=150)

    plt.show()
    print("仿真完成!")
