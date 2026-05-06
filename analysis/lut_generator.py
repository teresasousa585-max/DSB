"""
LUT 计算与对拍工具

1. 从 beam_controller.sv 中提取硬编码 LUT
2. 根据物理公式重新计算
3. 逐元素对比，生成差异报告
4. 输出可直接粘贴回 Verilog 的 generated_lut.sv
"""

import os
import re
from pathlib import Path

import numpy as np
from scipy.signal.windows import hann, hamming, blackman

from analysis.common import (
    C, D, LAMBDA, PHASE_STEPS,
    N_ELEMENTS, ELEMENT_ROWS, ELEMENT_COLS,
    STEER_MIN, STEER_MAX, STEER_STEP, STEER_POINTS,
    FOCUS_MIN_MM, FOCUS_MAX_MM, FOCUS_STEP_MM, FOCUS_POINTS,
    AMP_MAX,
)

# ============================================================================
# Verilog 解析
# ============================================================================

def _read_beam_controller():
    """读取 beam_controller.sv 全文。"""
    sv_path = Path(__file__).resolve().parent.parent / "verilog" / "beam_controller.sv"
    with open(sv_path, "r", encoding="utf-8") as f:
        return f.read()


def parse_steer_lut(text: str) -> np.ndarray:
    """提取 steer_phase_lut[0:60]（有符号十进制）。"""
    pattern = re.compile(r"steer_phase_lut\[\s*(\d+)\]\s*=\s*(-?10'sd\d+);")
    matches = pattern.findall(text)
    lut = np.zeros(STEER_POINTS, dtype=np.int16)
    for idx_str, val_str in matches:
        idx = int(idx_str)
        val = int(val_str.replace("-10'sd", "-").replace("10'sd", ""))
        lut[idx] = val
    return lut


def parse_focus_lut(text: str) -> np.ndarray:
    """提取 focus_phase_lut[0:98][0:24]（无符号十进制）。"""
    pattern = re.compile(r"focus_phase_lut\[\s*(\d+)\]\[\s*(\d+)\]\s*=\s*10'd(\d+);")
    matches = pattern.findall(text)
    lut = np.zeros((FOCUS_POINTS, N_ELEMENTS), dtype=np.uint16)
    for dist_idx_str, elem_idx_str, val_str in matches:
        dist_idx = int(dist_idx_str)
        elem_idx = int(elem_idx_str)
        val = int(val_str)
        lut[dist_idx, elem_idx] = val
    return lut


def parse_window_lut(text: str, name: str) -> np.ndarray:
    """提取 window_xxx[0:24]（无符号十进制）。"""
    pattern = re.compile(rf"assign\s+window_{name}\[(\d+)\]\s*=\s*8'd(\d+);")
    matches = pattern.findall(text)
    lut = np.zeros(N_ELEMENTS, dtype=np.uint8)
    for idx_str, val_str in matches:
        idx = int(idx_str)
        val = int(val_str)
        lut[idx] = val
    return lut


# ============================================================================
# LUT 重新计算
# ============================================================================

def compute_steer_lut() -> np.ndarray:
    """
    偏转相位增量 LUT。
    相邻阵元(16mm间距)的相位差 = 1024 * d * sin(theta) / lambda
    """
    thetas = np.arange(STEER_MIN, STEER_MAX + 1, STEER_STEP, dtype=np.float64)
    lut = np.round(PHASE_STEPS * D * np.sin(np.deg2rad(thetas)) / LAMBDA).astype(np.int16)
    return lut


def compute_focus_lut() -> np.ndarray:
    """
    聚焦相位 LUT。
    对于每个聚焦距离 r0 和每个阵元 i:
        ri = sqrt(xi^2 + yi^2 + r0^2)
        phase = round(1024 * (ri - r0) / lambda) mod 1024
    """
    focus_dists_mm = np.arange(FOCUS_MIN_MM, FOCUS_MAX_MM + 1, FOCUS_STEP_MM, dtype=np.float64)
    focus_dists_m = focus_dists_mm / 1000.0

    lut = np.zeros((FOCUS_POINTS, N_ELEMENTS), dtype=np.uint16)

    xs = ELEMENT_COLS * D
    ys = ELEMENT_ROWS * D

    for d_idx, r0 in enumerate(focus_dists_m):
        ri = np.sqrt(xs**2 + ys**2 + r0**2)
        phases = np.round(PHASE_STEPS * (ri - r0) / LAMBDA).astype(np.int64)
        phases = phases % PHASE_STEPS
        lut[d_idx, :] = phases.astype(np.uint16)

    return lut


def compute_window_lut(window_name: str) -> np.ndarray:
    """
    计算标准 2D 可分离窗函数。

    使用 scipy.signal.windows 的 1D 窗，做外积得到 2D 窗，
    然后归一化到 0..255 并四舍五入。

    注意：Verilog 中的硬编码窗值可能与标准公式有差异（近似或自定义），
    本函数仅提供标准值作为参考。
    """
    M = 5
    if window_name == "hann":
        w1d = hann(M)
    elif window_name == "hamm":
        w1d = hamming(M)
    elif window_name == "blk":
        w1d = blackman(M)
    else:
        raise ValueError(f"Unknown window: {window_name}")

    # 归一化到 0..1（scipy 默认最大值不一定是 1，需显式归一化）
    w1d = w1d / w1d.max()

    # 2D 可分离窗
    w2d = np.outer(w1d, w1d)

    # 转为一维 row-major，映射到 0..255
    lut = np.round(w2d.flatten() * AMP_MAX).astype(np.uint8)
    return lut


# ============================================================================
# 对比与报告
# ============================================================================

def compare_luts(name: str, expected: np.ndarray, actual: np.ndarray, tol: int = 0) -> bool:
    """对比两个 LUT，打印差异报告。"""
    diff = np.abs(expected.astype(np.int64) - actual.astype(np.int64))
    mismatches = np.where(diff > tol)
    mismatch_count = len(mismatches[0])

    total = expected.size
    if mismatch_count == 0:
        print(f"  [OK] {name}: 全部 {total} 个值匹配（容差 ±{tol}）")
        return True
    else:
        print(f"  [WARN] {name}: {mismatch_count}/{total} 个值不匹配（容差 ±{tol}）")
        # 最多打印前 10 条
        for i in range(min(mismatch_count, 10)):
            idx = tuple(m[i] for m in mismatches)
            print(f"    索引 {idx}: 期望={expected[idx]}, 实际={actual[idx]}, 差={diff[idx]}")
        if mismatch_count > 10:
            print(f"    ... 还有 {mismatch_count - 10} 处差异未显示")
        return False


def generate_sv_lut(steer_lut, focus_lut, window_luts, out_path):
    """生成可粘贴回 beam_controller.sv 的 LUT 代码。"""
    lines = []
    lines.append("// ==========================================================================")
    lines.append("// 自动生成的 LUT — 由 analysis/lut_generator.py 生成")
    lines.append("// ==========================================================================")
    lines.append("")

    # Steering LUT
    lines.append("    // 偏转相位增量LUT")
    lines.append("    initial begin")
    for i, val in enumerate(steer_lut):
        theta = STEER_MIN + i
        lines.append(f"        steer_phase_lut[{i:2d}] = {'-' if val < 0 else ''}10'sd{abs(val)};  // theta={theta:+d}")
    lines.append("    end")
    lines.append("")

    # Focus LUT
    lines.append("    // 聚焦相位LUT")
    lines.append("    initial begin")
    for d_idx in range(FOCUS_POINTS):
        dist_mm = FOCUS_MIN_MM + d_idx * FOCUS_STEP_MM
        lines.append(f"        // focus_dist={dist_mm}mm (idx={d_idx})")
        for e_idx in range(N_ELEMENTS):
            val = int(focus_lut[d_idx, e_idx])
            lines.append(f"        focus_phase_lut[{d_idx:2d}][{e_idx:2d}] = 10'd{val};")
    lines.append("    end")
    lines.append("")

    # Window LUTs
    for name, lut in window_luts.items():
        lines.append(f"    // {name} 窗 LUT")
        for i in range(N_ELEMENTS):
            lines.append(f"    assign window_{name}[{i:2d}] = 8'd{lut[i]};")
        lines.append("")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"  已生成 SV LUT 文件: {out_path}")


# ============================================================================
# 主入口
# ============================================================================

def run():
    print("=" * 70)
    print("LUT 计算与对拍")
    print("=" * 70)

    # 读取 Verilog
    print("\n[1/4] 解析 beam_controller.sv ...")
    sv_text = _read_beam_controller()
    sv_steer = parse_steer_lut(sv_text)
    sv_focus = parse_focus_lut(sv_text)
    sv_windows = {
        "hann": parse_window_lut(sv_text, "hann"),
        "hamm": parse_window_lut(sv_text, "hamm"),
        "blk": parse_window_lut(sv_text, "blk"),
    }
    print(f"  提取到 steering LUT: {sv_steer.shape}")
    print(f"  提取到 focus LUT: {sv_focus.shape}")
    print(f"  提取到窗口 LUT: hann={sv_windows['hann'].shape}, hamm={sv_windows['hamm'].shape}, blk={sv_windows['blk'].shape}")

    # 重新计算
    print("\n[2/4] 重新计算 LUT ...")
    calc_steer = compute_steer_lut()
    calc_focus = compute_focus_lut()
    calc_windows = {
        "hann": compute_window_lut("hann"),
        "hamm": compute_window_lut("hamm"),
        "blk": compute_window_lut("blk"),
    }

    # 对比
    print("\n[3/4] 对比结果 ...")
    all_ok = True
    all_ok &= compare_luts("偏转相位增量", calc_steer, sv_steer, tol=0)
    all_ok &= compare_luts("聚焦相位", calc_focus, sv_focus, tol=1)
    for name in ("hann", "hamm", "blk"):
        all_ok &= compare_luts(f"窗口 ({name})", calc_windows[name], sv_windows[name], tol=0)

    if all_ok:
        print("\n[PASS] 所有 LUT 验证通过！")
    else:
        print("\n[FAIL] 存在 LUT 差异，请检查上述报告。")
        print("       提示：聚焦相位允许 ±1 的舍入误差；")
        print("       窗函数差异可能是由于 Verilog 使用了近似/自定义公式。")

    # 生成 SV 文件
    print("\n[4/4] 生成 generated_lut.sv ...")
    out_dir = Path(__file__).resolve().parent.parent / "analysis_outputs"
    out_dir.mkdir(exist_ok=True)
    generate_sv_lut(calc_steer, calc_focus, calc_windows, out_dir / "generated_lut.sv")

    print("\n" + "=" * 70)


if __name__ == "__main__":
    run()
