"""
N32 同心圆环阵列：±15° x轴偏转聚焦的参数公式推导

方法：
  1. 相位 = 理论聚焦相位 (复用 software_weighting_compare.phase_only_weights)
  2. 幅度 = RLS 约束优化 (复用 software_weighting_compare.rls_amplitude_weights)
  3. 对每个工况的幅度分布进行多项式拟合 A(ρ) = c0 + c1*ρ_n + c2*ρ_n^2 + c3*ρ_n^3
  4. 拟合系数 ci 随距离 r 和角度 θ 的变化公式

输出：
  - analysis_outputs/n32_focus_steer_formula_weights.csv
  - analysis_outputs/n32_focus_steer_formula_summary.txt
  - analysis_outputs/n32_focus_steer_formula.png
"""

import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import software_weighting_compare as sw
import design_n36_array as base

# ==================== 基础参数 ====================
C = 343.0
FREQ = 40000.0
LAMBDA = C / FREQ
K = 2.0 * np.pi / LAMBDA

# 扫描工况：距离 (mm) 和 偏转角度 (deg)
DISTANCES_MM = [100, 150, 200, 300, 400, 500, 700, 1000, 1500, 2000, 3000, 4000, 4999]
STEER_ANGLES_DEG = [-15.0, -7.5, 0.0, 7.5, 15.0]

ROOT = os.path.dirname(os.path.dirname(__file__))
OUT_DIR = os.path.join(ROOT, "analysis_outputs")
LAYOUT_CSV = os.path.join(os.path.dirname(__file__), "n32_array_coordinates.csv")


def ensure_out_dir():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)


def read_layout(path):
    points = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            points.append([float(row["X_mm"]) * 1e-3, float(row["Y_mm"]) * 1e-3])
    return np.asarray(points, dtype=float)


def make_case(distance_m, theta_x_deg):
    """构建 software_weighting_compare 所需的 case 字典"""
    sx = math.sin(math.radians(theta_x_deg))
    sz = math.cos(math.radians(theta_x_deg))
    x_m = distance_m * sx
    z_m = distance_m * sz
    return {
        "name": "focus_r%.0f_t%g" % (distance_m * 1000, theta_x_deg),
        "domain": "near",
        "x_m": x_m,
        "y_m": 0.0,
        "z_m": z_m,
    }


def fit_amplitude_per_case(points_mm, amp, rmax):
    """对每个工况的幅度分布进行多项式拟合，返回系数 [c0, c1, c2, c3]"""
    rho = np.sqrt(points_mm[:, 0] ** 2 + points_mm[:, 1] ** 2)
    rho_norm = rho / rmax
    X = np.column_stack([np.ones_like(rho_norm), rho_norm, rho_norm ** 2, rho_norm ** 3])
    c, _, _, _ = np.linalg.lstsq(X, amp, rcond=None)
    return c


def fit_coefficient_model(all_rows):
    """
    拟合多项式系数 c0, c1, c2, c3 随距离 r 和角度 θ 的变化。
    模型：ci = a0 + a1*r_m + a2*θ_n + a3*r_m^2 + a4*θ_n^2 + a5*r_m*θ_n
    """
    coeffs_list = np.array([row["poly_coeffs"] for row in all_rows])
    r_vals = np.array([row["distance_m"] * 1000.0 for row in all_rows])
    theta_vals = np.array([abs(row["theta_x_deg"]) for row in all_rows])

    r_m = r_vals / 1000.0
    th_n = theta_vals / 15.0

    X = np.column_stack([
        np.ones(len(r_m)),
        r_m,
        th_n,
        r_m ** 2,
        th_n ** 2,
        r_m * th_n,
    ])

    fitted_models = []
    for i in range(4):
        y = coeffs_list[:, i]
        a, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        fitted_models.append(a)

    def predict_coeffs(r_mm, theta_deg):
        r_m = r_mm / 1000.0
        th_n = abs(theta_deg) / 15.0
        x_vec = np.array([1.0, r_m, th_n, r_m ** 2, th_n ** 2, r_m * th_n])
        return np.array([np.dot(a, x_vec) for a in fitted_models])

    return fitted_models, predict_coeffs


def write_weights_csv(path, points_mm, rows):
    fields = [
        "case", "distance_mm", "theta_x_deg", "theta_y_deg",
        "element_id", "x_mm", "y_mm", "rho_mm",
        "amplitude", "phase_deg"
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            phase_w = row["weights"]
            for i, p in enumerate(points_mm):
                rho = math.hypot(p[0], p[1])
                writer.writerow({
                    "case": row["case"],
                    "distance_mm": row["distance_mm"],
                    "theta_x_deg": row["theta_x_deg"],
                    "theta_y_deg": row["theta_y_deg"],
                    "element_id": "E%02d" % i,
                    "x_mm": "%.3f" % (p[0] * 1000.0),
                    "y_mm": "%.3f" % (p[1] * 1000.0),
                    "rho_mm": "%.3f" % (rho * 1000.0),
                    "amplitude": "%.6f" % np.abs(phase_w[i]),
                    "phase_deg": "%.3f" % (np.rad2deg(np.angle(phase_w[i])) % 360.0),
                })


def write_summary(path, points_mm, rows, fitted_models, rmax, formula_metrics):
    rho = np.sqrt(points_mm[:, 0] ** 2 + points_mm[:, 1] ** 2)
    rho_norm = rho / rmax

    with open(path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("N32 超声阵列 ±15° 偏转聚焦参数公式推导报告\n")
        f.write("=" * 70 + "\n\n")
        f.write("阵列参数:\n")
        f.write("  频率: %.1f Hz\n" % FREQ)
        f.write("  声速: %.1f m/s\n" % C)
        f.write("  波长: %.3f mm\n" % (LAMBDA * 1000.0))
        f.write("  阵元数: %d\n" % len(points_mm))
        f.write("  阵元直径: %.1f mm\n" % (base.ELEMENT_DIAMETER_M * 1000.0))
        f.write("  孔径半径: %.1f mm\n" % (rmax * 1000.0))
        f.write("\n")

        f.write("扫描工况:\n")
        for row in rows:
            f.write("  %s: r=%4.0f mm, theta_x=%+5.1f°\n" % (
                row["case"], row["distance_mm"], row["theta_x_deg"]
            ))
        f.write("\n")

        f.write("优化指标 (RLS 幅度优化 + 理论聚焦相位):\n")
        f.write("%-12s %6s %6s %8s %8s %8s %8s\n" % (
            "case", "dist", "theta", "PSL_dB", "dBc", "gain_dB", "meanAmp"
        ))
        for row in rows:
            m = row["metrics"]
            f.write("%-12s %6.0f %6.1f %8.2f %8.2f %8.2f %8.3f\n" % (
                row["case"], row["distance_mm"], row["theta_x_deg"],
                m["psl_db"], m["dBc"],
                m["gain_db"], m["mean_amp"]
            ))
        f.write("\n")

        f.write("相位公式:\n")
        f.write("  φ_i(r, θ) = k * (R_i - r)   [rad]\n")
        f.write("  R_i = sqrt((r*sinθ - x_i)^2 + y_i^2 + (r*cosθ)^2)\n")
        f.write("  k = %.6f rad/m\n\n" % K)

        f.write("幅度参数公式:\n")
        f.write("  A_i(r, θ) = clip(c0 + c1*ρ_n + c2*ρ_n^2 + c3*ρ_n^3, 0, 1)\n")
        f.write("  ρ_n = ρ_i / %.3f  (归一化径向位置, ρ_i 单位 mm)\n\n" % (rmax * 1000.0))

        f.write("  系数 ci 随 r 和 θ 的变化:\n")
        f.write("    ci(r, θ) = a0 + a1*r_m + a2*θ_n + a3*r_m^2 + a4*θ_n^2 + a5*r_m*θ_n\n")
        f.write("    r_m = r / 1000   (距离，单位米)\n")
        f.write("    θ_n = |θ| / 15   (归一化偏转角，|θ| ≤ 15°)\n\n")

        labels = ["c0", "c1", "c2", "c3"]
        for i, lab in enumerate(labels):
            a = fitted_models[i]
            f.write("    %s: a0=% .6f, a1=% .6f, a2=% .6f, a3=% .6f, a4=% .6f, a5=% .6f\n" % (
                lab, a[0], a[1], a[2], a[3], a[4], a[5]
            ))
        f.write("\n")

        f.write("优化幅度按环统计 (均值 ± 标准差):\n")
        ring_names = ["Center", "Ring1", "Ring2", "Ring3"]
        ring_indices = [[0], list(range(1, 8)), list(range(8, 17)), list(range(17, 32))]
        for rname, idxs in zip(ring_names, ring_indices):
            amps = np.hstack([np.abs(row["weights"][idxs]) for row in rows])
            f.write("  %-8s: %.4f ± %.4f (range: %.4f ~ %.4f)\n" % (
                rname, np.mean(amps), np.std(amps), np.min(amps), np.max(amps)
            ))
        f.write("\n")

        f.write("公式验证 (拟合公式 vs RLS 优化):\n")
        f.write("%-12s %6s %6s %8s %8s %8s %8s\n" % (
            "case", "dist", "theta", "opt_PSL", "form_PSL", "ΔPSL", "Δgain"
        ))
        for row, fm in zip(rows, formula_metrics):
            f.write("%-12s %6.0f %6.1f %8.2f %8.2f %8.3f %8.3f\n" % (
                row["case"], row["distance_mm"], row["theta_x_deg"],
                row["metrics"]["psl_db"],
                fm["psl_db"],
                fm["psl_db"] - row["metrics"]["psl_db"],
                fm["gain_db"] - row["metrics"]["gain_db"],
            ))
        f.write("\n")

        f.write("=" * 70 + "\n")
        f.write("使用说明:\n")
        f.write("1. 对每个阵元 i (坐标 x_i, y_i [m], 径向位置 ρ_i [mm]):\n")
        f.write("   R_i = sqrt((r*sinθ - x_i)^2 + y_i^2 + (r*cosθ)^2)\n")
        f.write("   φ_i = k * (R_i - r)  [rad]\n")
        f.write("2. ρ_n = ρ_i / %.3f\n" % (rmax * 1000.0))
        f.write("   r_m = r / 1000,  θ_n = |θ| / 15\n")
        f.write("3. 计算系数:\n")
        for i, lab in enumerate(labels):
            a = fitted_models[i]
            f.write("     %s = % .4f + % .4f*r_m + % .4f*θ_n + % .4f*r_m^2 + % .4f*θ_n^2 + % .4f*r_m*θ_n\n" % (
                lab, a[0], a[1], a[2], a[3], a[4], a[5]
            ))
        f.write("4. A_i = max(0, min(1, c0 + c1*ρ_n + c2*ρ_n^2 + c3*ρ_n^3))\n")
        f.write("5. w_i = A_i * exp(-j * φ_i)\n")
        f.write("6. 归一化: w_i /= sum(|w_i|)\n")
        f.write("=" * 70 + "\n")


def plot_results(points_mm, rows, fitted_models, rmax, formula_metrics, path):
    rho = np.sqrt(points_mm[:, 0] ** 2 + points_mm[:, 1] ** 2)
    rho_norm = rho / rmax
    n_rows = len(rows)

    fig = plt.figure(figsize=(16, 10))
    gs = matplotlib.gridspec.GridSpec(3, 3, hspace=0.35, wspace=0.3)

    # 图1: 优化幅度 vs 径向位置
    ax1 = fig.add_subplot(gs[0, 0])
    for row in rows:
        color = "blue" if row["theta_x_deg"] > 0 else "red"
        ax1.scatter(rho * 1000.0, np.abs(row["weights"]), c=color, s=20, alpha=0.5)
    ax1.set_xlabel("ρ (mm)")
    ax1.set_ylabel("optimized amplitude")
    ax1.set_title("Amplitude vs radial position (RLS opt)")
    ax1.grid(True, alpha=0.25)

    # 图2: 幅度按环统计
    ax2 = fig.add_subplot(gs[0, 1])
    ring_names = ["Center", "Ring1", "Ring2", "Ring3"]
    ring_indices = [[0], list(range(1, 8)), list(range(8, 17)), list(range(17, 32))]
    ring_means = []
    ring_stds = []
    for idxs in ring_indices:
        amps = np.hstack([np.abs(row["weights"][idxs]) for row in rows])
        ring_means.append(np.mean(amps))
        ring_stds.append(np.std(amps))
    x_pos = np.arange(len(ring_names))
    ax2.bar(x_pos, ring_means, yerr=ring_stds, capsize=5, color="steelblue", alpha=0.8, edgecolor="black")
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(ring_names)
    ax2.set_ylabel("mean amplitude")
    ax2.set_title("Amplitude by ring")
    ax2.set_ylim(0, 1.05)
    ax2.grid(True, axis="y", alpha=0.25)

    # 图3: 系数随距离变化
    ax3 = fig.add_subplot(gs[0, 2])
    dists = sorted(set(row["distance_mm"] for row in rows))
    for angle in [-15.0, 15.0]:
        c0_vals = []
        c1_vals = []
        for d in dists:
            row = next((r for r in rows if r["distance_mm"] == d and abs(r["theta_x_deg"] - angle) < 0.1), None)
            if row:
                c0_vals.append(row["poly_coeffs"][0])
                c1_vals.append(row["poly_coeffs"][1])
        if c0_vals:
            ax3.plot(dists, c0_vals, marker="o", label="c0 θ=%+.0f°" % angle)
            ax3.plot(dists, c1_vals, marker="s", linestyle="--", label="c1 θ=%+.0f°" % angle)
    ax3.set_xlabel("distance (mm)")
    ax3.set_ylabel("coefficient value")
    ax3.set_title("Fitted coeffs vs distance")
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.25)

    # 图4-6: 典型工况的公式 vs 优化幅度对比
    example_cases = [0, 5, 9]
    for plot_idx, case_idx in enumerate(example_cases):
        if case_idx >= n_rows:
            continue
        row = rows[case_idx]
        ax = fig.add_subplot(gs[1, plot_idx])
        ax.scatter(rho * 1000.0, np.abs(row["weights"]), c="blue", s=40, alpha=0.7, label="RLS optimized")

        r_mm = row["distance_mm"]
        theta_deg = row["theta_x_deg"]
        c_pred = np.array([np.dot(a, np.array([1.0, r_mm/1000.0, abs(theta_deg)/15.0, (r_mm/1000.0)**2, (abs(theta_deg)/15.0)**2, (r_mm/1000.0)*(abs(theta_deg)/15.0)])) for a in fitted_models])
        pred = np.clip(c_pred[0] + c_pred[1]*rho_norm + c_pred[2]*rho_norm**2 + c_pred[3]*rho_norm**3, 0.0, 1.0)
        ax.scatter(rho * 1000.0, pred, c="red", s=40, alpha=0.7, marker="x", label="formula")
        ax.set_xlabel("ρ (mm)")
        ax.set_ylabel("amplitude")
        ax.set_title("r=%.0f mm, θ=%+.0f°" % (r_mm, theta_deg))
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.25)
        if plot_idx == 0:
            ax.legend()

    # 图7-9: 指标对比
    x_labels = ["%s\n%.0fmm" % (row["case"], row["distance_mm"]) for row in rows]
    x_pos = np.arange(len(rows))
    width = 0.35

    ax7 = fig.add_subplot(gs[2, 0])
    opt_psl = [row["metrics"]["psl_db"] for row in rows]
    form_psl = [fm["psl_db"] for fm in formula_metrics]
    ax7.bar(x_pos - width/2, opt_psl, width, label="RLS opt", color="steelblue")
    ax7.bar(x_pos + width/2, form_psl, width, label="formula", color="coral")
    ax7.set_xticks(x_pos)
    ax7.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=7)
    ax7.set_ylabel("PSL (dB)")
    ax7.set_title("Sidelobe comparison")
    ax7.legend()
    ax7.grid(True, axis="y", alpha=0.25)

    ax8 = fig.add_subplot(gs[2, 1])
    opt_gain = [row["metrics"]["gain_db"] for row in rows]
    form_gain = [fm["gain_db"] for fm in formula_metrics]
    ax8.bar(x_pos - width/2, opt_gain, width, label="RLS opt", color="steelblue")
    ax8.bar(x_pos + width/2, form_gain, width, label="formula", color="coral")
    ax8.set_xticks(x_pos)
    ax8.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=7)
    ax8.set_ylabel("gain loss (dB)")
    ax8.set_title("Gain loss comparison")
    ax8.legend()
    ax8.grid(True, axis="y", alpha=0.25)

    ax9 = fig.add_subplot(gs[2, 2])
    opt_dbc = [row["metrics"]["dBc"] for row in rows]
    form_dbc = [fm["dBc"] for fm in formula_metrics]
    ax9.bar(x_pos - width/2, opt_dbc, width, label="RLS opt", color="steelblue")
    ax9.bar(x_pos + width/2, form_dbc, width, label="formula", color="coral")
    ax9.set_xticks(x_pos)
    ax9.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=7)
    ax9.set_ylabel("dBc")
    ax9.set_title("Main-to-sidelobe ratio")
    ax9.legend()
    ax9.grid(True, axis="y", alpha=0.25)

    fig.suptitle("N32 Focus+Steer ±15°: RLS Optimization vs Fitted Formula", fontsize=14)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    ensure_out_dir()
    points_mm = read_layout(LAYOUT_CSV)
    n = len(points_mm)
    print("Loaded N32 layout: %d elements" % n)

    rho = np.sqrt(points_mm[:, 0] ** 2 + points_mm[:, 1] ** 2)
    rmax = np.max(rho)
    rho_norm = rho / rmax

    rows = []
    for dist_mm in DISTANCES_MM:
        for theta_x in STEER_ANGLES_DEG:
            case_name = "R%.0f_x%g" % (dist_mm, theta_x)
            print("Processing %s ..." % case_name)

            distance_m = dist_mm * 1e-3
            case = make_case(distance_m, theta_x)

            # Phase-only weights
            w_phase = sw.phase_only_weights(points_mm, case)
            # RLS amplitude optimization
            w_rls = sw.rls_amplitude_weights(points_mm, case)

            # Evaluate
            psl_phase, p99_phase, width_phase, target_amp_phase, err_phase, x_mm, y_mm, db_phase = sw.nearfield_eval(points_mm, case, w_phase)
            psl_rls, p99_rls, width_rls, target_amp_rls, err_rls, x_mm, y_mm, db_rls = sw.nearfield_eval(points_mm, case, w_rls)

            # Fit polynomial to amplitude distribution
            amp = np.abs(w_rls)
            poly_coeffs = fit_amplitude_per_case(points_mm, amp, rmax)

            rows.append({
                "case": case_name,
                "distance_mm": dist_mm,
                "distance_m": distance_m,
                "theta_x_deg": theta_x,
                "theta_y_deg": 0.0,
                "weights": w_rls,
                "poly_coeffs": poly_coeffs,
                "metrics": {
                    "psl_db": psl_rls,
                    "dBc": -psl_rls,
                    "gain_db": 20.0 * math.log10(target_amp_rls / max(target_amp_phase, 1e-18) + 1e-12),
                    "mean_amp": float(np.mean(amp)),
                },
            })
            print("  Phase-only PSL=%.2f dB | RLS PSL=%.2f dB, gain=%.2f dB, meanAmp=%.3f" % (
                psl_phase, psl_rls,
                20.0 * math.log10(target_amp_rls / max(target_amp_phase, 1e-18) + 1e-12),
                np.mean(amp)
            ))

    # 拟合系数随 r 和 θ 的变化模型
    print("\nFitting coefficient model ...")
    fitted_models, predict_coeffs = fit_coefficient_model(rows)
    print("  Models fitted for c0, c1, c2, c3")

    # 用公式验证每个工况
    print("\nValidating formula ...")
    formula_metrics = []
    for row in rows:
        dist_mm = row["distance_mm"]
        theta_x = row["theta_x_deg"]
        distance_m = dist_mm * 1e-3
        case = make_case(distance_m, theta_x)

        # 公式预测幅度
        c = predict_coeffs(dist_mm, theta_x)
        formula_amp = np.clip(c[0] + c[1]*rho_norm + c[2]*rho_norm**2 + c[3]*rho_norm**3, 0.0, 1.0)

        # 理论聚焦相位
        w_phase = sw.phase_only_weights(points_mm, case)
        formula_weights = formula_amp * (w_phase / np.abs(w_phase))

        psl_form, p99_form, width_form, target_amp_form, err_form, x_mm, y_mm, db_form = sw.nearfield_eval(points_mm, case, formula_weights)

        psl_phase = row["metrics"]["psl_db"]  # 实际上这里保存的是 RLS 的 PSL，需要重新计算 phase-only 的 gain
        # 重新计算 phase-only 目标增益用于对比
        w_phase_only = sw.phase_only_weights(points_mm, case)
        _, _, _, target_amp_phase, _, _, _, _ = sw.nearfield_eval(points_mm, case, w_phase_only)

        formula_metrics.append({
            "psl_db": psl_form,
            "dBc": -psl_form,
            "gain_db": 20.0 * math.log10(target_amp_form / max(target_amp_phase, 1e-18) + 1e-12),
        })
        print("  %s: formula PSL=%.2f dB (opt=%.2f), dBc=%.2f" % (
            row["case"], psl_form, row["metrics"]["psl_db"], -psl_form
        ))

    # 输出 CSV
    weights_csv = os.path.join(OUT_DIR, "n32_focus_steer_formula_weights.csv")
    write_weights_csv(weights_csv, points_mm, rows)
    print("\nWrote:", weights_csv)

    # 输出 Summary
    summary_txt = os.path.join(OUT_DIR, "n32_focus_steer_formula_summary.txt")
    write_summary(summary_txt, points_mm, rows, fitted_models, rmax, formula_metrics)
    print("Wrote:", summary_txt)

    # 输出 Plot
    plot_png = os.path.join(OUT_DIR, "n32_focus_steer_formula.png")
    plot_results(points_mm, rows, fitted_models, rmax, formula_metrics, plot_png)
    print("Wrote:", plot_png)


if __name__ == "__main__":
    main()
