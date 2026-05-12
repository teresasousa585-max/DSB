import csv
import math
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import j1


ROOT = os.path.dirname(os.path.dirname(__file__))
OUT_DIR = os.path.join(ROOT, "analysis_outputs")

C = 343.0
FREQ = 40000.0
LAMBDA = C / FREQ
K = 2.0 * np.pi / LAMBDA
ELEMENT_RADIUS = 8e-3
SEGMENT_WIDTH = 100e-3
RING_RADII = [18e-3, 38e-3]
RING_COUNTS = [4, 10]

TARGET_RANGE_M = 2.0
TARGET_ANGLES_DEG = [-60, -45, -30, -15, 0, 15, 30, 45, 60]
SEARCH_ANGLES_DEG = [0, 15, 30, 45, 60]
MAIN_EXCLUDE_MM = 120.0
DB_FLOOR = -60.0


def ensure_out_dir():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)


def piston_from_cos(cos_theta):
    cos_theta = np.asarray(cos_theta)
    front = cos_theta > 0.0
    sin_theta = np.sqrt(np.maximum(0.0, 1.0 - cos_theta * cos_theta))
    x = K * ELEMENT_RADIUS * sin_theta
    out = np.zeros_like(cos_theta, dtype=float)
    small = np.abs(x) < 1e-10
    out[front & small] = 1.0
    mask = front & (~small)
    out[mask] = np.abs(2.0 * j1(x[mask]) / x[mask])
    return out


def local_15_positions():
    pts = [[0.0, 0.0]]
    for r, n in zip(RING_RADII, RING_COUNTS):
        for i in range(n):
            a = 2.0 * np.pi * i / n
            pts.append([r * math.cos(a), r * math.sin(a)])
    return np.asarray(pts, dtype=float)


def connected_array(a1_deg, a2_deg):
    tilts = [float(a2_deg), float(a1_deg), 0.0, -float(a1_deg), -float(a2_deg)]
    local = local_15_positions()
    uxs = []
    normals_by_seg = []
    for tilt_deg in tilts:
        alpha = math.radians(tilt_deg)
        normals_by_seg.append(np.asarray([math.sin(alpha), 0.0, math.cos(alpha)]))
        uxs.append(np.asarray([math.cos(alpha), 0.0, -math.sin(alpha)]))

    centers = [None] * 5
    centers[2] = np.asarray([0.0, 0.0, 0.0])
    for seg in [1, 0]:
        hinge = centers[seg + 1] - 0.5 * SEGMENT_WIDTH * uxs[seg + 1]
        centers[seg] = hinge - 0.5 * SEGMENT_WIDTH * uxs[seg]
    for seg in [3, 4]:
        hinge = centers[seg - 1] + 0.5 * SEGMENT_WIDTH * uxs[seg - 1]
        centers[seg] = hinge + 0.5 * SEGMENT_WIDTH * uxs[seg]

    positions = []
    normals = []
    segment_ids = []
    for seg in range(5):
        uy = np.asarray([0.0, 1.0, 0.0])
        for xl, yl in local:
            positions.append(centers[seg] + xl * uxs[seg] + yl * uy)
            normals.append(normals_by_seg[seg])
            segment_ids.append(seg)
    return np.asarray(positions), np.asarray(normals), np.asarray(segment_ids), tilts, np.asarray(centers)


def target_for_angle(angle_deg, range_m=TARGET_RANGE_M):
    a = math.radians(angle_deg)
    return np.asarray([range_m * math.sin(a), 0.0, range_m * math.cos(a)])


def nearfield_matrix(positions, normals, x, y, z):
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    pts = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
    diff = pts[:, None, :] - positions[None, :, :]
    r = np.linalg.norm(diff, axis=2)
    dirs = diff / np.maximum(r[:, :, None], 1e-12)
    cosang = np.sum(dirs * normals[None, :, :], axis=2)
    elem = piston_from_cos(cosang)
    return elem * np.exp(1j * K * r) / np.maximum(r, 1e-12)


def phase_weights(h_target, amp=None):
    if amp is None:
        amp = np.ones_like(h_target, dtype=float)
    return amp * np.exp(-1j * np.angle(h_target))


def strategy_weights(strategy, positions, normals, segment_ids, target):
    ht = nearfield_matrix(positions, normals, np.asarray([target[0]]), np.asarray([target[1]]), np.asarray([target[2]]))[0]
    if strategy == "all":
        return phase_weights(ht)

    seg_score = []
    for seg in range(5):
        mask = segment_ids == seg
        seg_score.append(float(np.sum(np.abs(ht[mask]))))

    if strategy.startswith("best"):
        nseg = int(strategy[-1])
        keep = set(np.argsort(seg_score)[-nseg:])
        amp = np.asarray([1.0 if int(seg) in keep else 0.0 for seg in segment_ids])
        return phase_weights(ht, amp)

    if strategy == "soft_gate":
        amp_seg = np.asarray(seg_score) / max(max(seg_score), 1e-18)
        amp_seg = amp_seg ** 1.5
        amp_seg[amp_seg < 0.18] = 0.0
        amp = np.asarray([amp_seg[int(seg)] for seg in segment_ids])
        return phase_weights(ht, amp)

    raise ValueError(strategy)


def line_width(axis_mm, line_db, idx, threshold_db=-3.0):
    left = idx
    while left > 0 and line_db[left] >= threshold_db:
        left -= 1
    right = idx
    while right < len(line_db) - 1 and line_db[right] >= threshold_db:
        right += 1
    if left == 0 or right == len(line_db) - 1:
        return float("nan")

    def interp(i0, i1):
        if abs(line_db[i1] - line_db[i0]) < 1e-12:
            return axis_mm[i0]
        return axis_mm[i0] + (threshold_db - line_db[i0]) * (axis_mm[i1] - axis_mm[i0]) / (line_db[i1] - line_db[i0])

    return float(interp(right - 1, right) - interp(left, left + 1))


def eval_focus(positions, normals, segment_ids, target, strategy, span_mm=650.0, step_mm=20.0):
    w = strategy_weights(strategy, positions, normals, segment_ids, target)
    x0, y0, z0 = target
    x_mm = np.arange(x0 * 1000.0 - span_mm, x0 * 1000.0 + span_mm + 0.001, step_mm)
    y_mm = np.arange(y0 * 1000.0 - span_mm, y0 * 1000.0 + span_mm + 0.001, step_mm)
    xg, yg = np.meshgrid(x_mm * 1e-3, y_mm * 1e-3)
    zg = np.zeros_like(xg) + z0
    h = nearfield_matrix(positions, normals, xg, yg, zg)
    resp = np.dot(h, w).reshape(xg.shape)
    ht = nearfield_matrix(positions, normals, np.asarray([x0]), np.asarray([y0]), np.asarray([z0]))[0]
    target_amp = float(abs(np.dot(ht, w)))
    db = 20.0 * np.log10(np.abs(resp) / max(target_amp, 1e-18) + 1e-12)
    rr = np.sqrt((x_mm[None, :] - x0 * 1000.0) ** 2 + (y_mm[:, None] - y0 * 1000.0) ** 2)
    outside = rr > MAIN_EXCLUDE_MM
    idx = np.unravel_index(np.argmax(np.abs(resp)), resp.shape)
    ix = int(np.argmin(np.abs(x_mm - x0 * 1000.0)))
    iy = int(np.argmin(np.abs(y_mm - y0 * 1000.0)))
    focus_error = float(math.hypot(x_mm[idx[1]] - x0 * 1000.0, y_mm[idx[0]] - y0 * 1000.0))
    return {
        "pressure_psl_db": float(np.max(db[outside])),
        "pressure_p99_5_db": float(np.percentile(db[outside], 99.5)),
        "target_amp": target_amp,
        "focus_error_mm": focus_error,
        "fwhm_x_mm": line_width(x_mm, db[iy, :], ix),
        "fwhm_y_mm": line_width(y_mm, db[:, ix], iy),
    }, x_mm, y_mm, np.maximum(db, DB_FLOOR)


def score_row(row):
    psl_penalty = max(0.0, row["pressure_psl_db"] + 8.0) * 2.0
    gain_penalty = max(0.0, -10.0 - row["gain_vs_ref_db"])
    err_penalty = max(0.0, row["focus_error_mm"] - 40.0) * 0.05
    return psl_penalty + gain_penalty + err_penalty


def search_candidates():
    strategies = ["all", "soft_gate", "best1", "best2", "best3"]
    candidate_rows = []
    detail_rows = []
    ref_positions, ref_normals, ref_seg, _, _ = connected_array(30, 60)
    ref_target = target_for_angle(0)
    ref_metrics, _, _, _ = eval_focus(ref_positions, ref_normals, ref_seg, ref_target, "all", span_mm=650.0, step_mm=20.0)
    ref_amp = ref_metrics["target_amp"]

    for a1 in range(15, 46, 5):
        for a2 in range(max(40, a1 + 10), 76, 5):
            positions, normals, segment_ids, tilts, centers = connected_array(a1, a2)
            best_rows = []
            for angle in SEARCH_ANGLES_DEG:
                target = target_for_angle(angle)
                trial_rows = []
                for strategy in strategies:
                    metrics, _, _, _ = eval_focus(positions, normals, segment_ids, target, strategy, span_mm=650.0, step_mm=20.0)
                    row = {
                        "a1_deg": a1,
                        "a2_deg": a2,
                        "angle_deg": angle,
                        "strategy": strategy,
                        "gain_vs_ref_db": 20.0 * math.log10(metrics["target_amp"] / max(ref_amp, 1e-18) + 1e-12),
                    }
                    row.update(metrics)
                    row["score"] = score_row(row)
                    trial_rows.append(row)
                    detail_rows.append(row)
                best = min(trial_rows, key=lambda r: r["score"])
                best_rows.append(best)
            worst_psl = max(r["pressure_psl_db"] for r in best_rows)
            min_gain = min(r["gain_vs_ref_db"] for r in best_rows)
            max_err = max(r["focus_error_mm"] for r in best_rows)
            mean_score = float(np.mean([r["score"] for r in best_rows]))
            candidate_rows.append({
                "a1_deg": a1,
                "a2_deg": a2,
                "tilts_left_to_right": "+%g/+%g/0/-%g/-%g" % (a2, a1, a1, a2),
                "mean_score": mean_score,
                "worst_psl_db": worst_psl,
                "min_gain_vs_ref_db": min_gain,
                "max_focus_error_mm": max_err,
                "best_strategies_0_15_30_45_60": ",".join(r["strategy"] for r in best_rows),
            })
    candidate_rows.sort(key=lambda r: (r["mean_score"], r["worst_psl_db"], -r["min_gain_vs_ref_db"]))
    return candidate_rows, detail_rows, ref_amp


def refine_candidate(a1, a2, ref_amp):
    strategies = ["all", "soft_gate", "best1", "best2", "best3"]
    positions, normals, segment_ids, tilts, centers = connected_array(a1, a2)
    rows = []
    heatmaps = []
    for angle in TARGET_ANGLES_DEG:
        target = target_for_angle(angle)
        trial_rows = []
        for strategy in strategies:
            metrics, x_mm, y_mm, db = eval_focus(positions, normals, segment_ids, target, strategy, span_mm=700.0, step_mm=10.0)
            row = {
                "a1_deg": a1,
                "a2_deg": a2,
                "angle_deg": angle,
                "target_x_mm": target[0] * 1000.0,
                "target_z_mm": target[2] * 1000.0,
                "strategy": strategy,
                "gain_vs_ref_db": 20.0 * math.log10(metrics["target_amp"] / max(ref_amp, 1e-18) + 1e-12),
            }
            row.update(metrics)
            row["score"] = score_row(row)
            trial_rows.append((row, x_mm, y_mm, db))
        best = min(trial_rows, key=lambda item: item[0]["score"])
        rows.append(best[0])
        if angle in [-60, -30, 0, 30, 60]:
            heatmaps.append(best)
    return rows, heatmaps, positions, normals, segment_ids, tilts, centers


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_summary(rows, path):
    angles = [r["angle_deg"] for r in rows]
    psl = [r["pressure_psl_db"] for r in rows]
    gain = [r["gain_vs_ref_db"] for r in rows]
    err = [r["focus_error_mm"] for r in rows]
    fig, axes = plt.subplots(3, 1, figsize=(8.5, 9.0), sharex=True)
    axes[0].plot(angles, psl, "o-")
    axes[0].axhline(-8.0, color="k", ls=":", lw=1)
    axes[0].set_ylabel("PSL outside 120 mm (dB)")
    axes[0].grid(True, alpha=0.25)
    axes[1].plot(angles, gain, "o-", color="tab:orange")
    axes[1].axhline(-10.0, color="k", ls=":", lw=1)
    axes[1].set_ylabel("target gain vs ref (dB)")
    axes[1].grid(True, alpha=0.25)
    axes[2].plot(angles, err, "o-", color="tab:green")
    axes[2].axhline(40.0, color="k", ls=":", lw=1)
    axes[2].set_ylabel("focus error (mm)")
    axes[2].set_xlabel("target angle at R=2 m (deg)")
    axes[2].grid(True, alpha=0.25)
    fig.suptitle("Recommended connected 5-segment strategy at R=2 m")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_geometry(positions, normals, segment_ids, tilts, centers, path):
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    for seg in range(5):
        mask = segment_ids == seg
        ax.scatter(positions[mask, 0] * 1000.0, positions[mask, 2] * 1000.0, s=22, label="seg%d %.0f deg" % (seg, tilts[seg]))
        c = centers[seg]
        n = normals[np.where(mask)[0][0]]
        ax.arrow(c[0] * 1000.0, c[2] * 1000.0, n[0] * 55.0, n[2] * 55.0, head_width=8.0, length_includes_head=True)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title("Connected folded geometry")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_heatmaps(heatmaps, path):
    fig, axes = plt.subplots(1, len(heatmaps), figsize=(3.7 * len(heatmaps), 4.0), sharey=True)
    if len(heatmaps) == 1:
        axes = [axes]
    contour = None
    for ax, (row, x_mm, y_mm, db) in zip(axes, heatmaps):
        contour = ax.contourf(x_mm, y_mm, db, levels=np.linspace(DB_FLOOR, 0.0, 49), cmap="magma", extend="min")
        ax.plot([row["target_x_mm"]], [0.0], "cx", ms=8, mew=2)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title("%+g deg %s\nPSL %.1f dB" % (row["angle_deg"], row["strategy"], row["pressure_psl_db"]))
        ax.set_xlabel("x (mm)")
        ax.grid(True, color="white", alpha=0.12)
    axes[0].set_ylabel("y (mm)")
    cbar = fig.colorbar(contour, ax=axes, shrink=0.78)
    cbar.set_label("relative to target (dB)")
    fig.suptitle("Best strategy focus maps at R=2 m")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary(path, candidates, rows, tilts):
    with open(path, "w") as f:
        f.write("Connected segmented angle search for R=2 m, target angles +/-60 deg\n")
        f.write("Assumption: target range is radial distance from array origin, not fixed z=2 m.\n")
        f.write("Element model: 16 mm circular piston at 40 kHz. Metrics use pressure field in a constant-z target plane.\n")
        f.write("Main exclusion radius for PSL: %.0f mm.\n\n" % MAIN_EXCLUDE_MM)
        f.write("Top candidate angle sets:\n")
        for row in candidates[:10]:
            f.write("  %-18s mean_score=%5.2f worst_PSL=%6.2f dB min_gain=%6.2f dB strategies=%s\n" % (
                row["tilts_left_to_right"], row["mean_score"], row["worst_psl_db"], row["min_gain_vs_ref_db"], row["best_strategies_0_15_30_45_60"]
            ))
        f.write("\nRecommended tilts left-to-right: %s\n\n" % ("/".join("%+.1f" % t for t in tilts)))
        f.write("Recommended per-angle strategy:\n")
        for row in rows:
            f.write("  angle=%+5.1f deg x=%8.1f z=%8.1f strategy=%-9s PSL=%6.2f dB gain=%6.2f dB err=%5.1f mm FWHM=(%6.1f,%6.1f) mm\n" % (
                row["angle_deg"], row["target_x_mm"], row["target_z_mm"], row["strategy"], row["pressure_psl_db"],
                row["gain_vs_ref_db"], row["focus_error_mm"], row["fwhm_x_mm"], row["fwhm_y_mm"]
            ))


def main():
    ensure_out_dir()
    candidates, detail_rows, ref_amp = search_candidates()
    cand_csv = os.path.join(OUT_DIR, "seg75_2m_angle_candidates.csv")
    write_csv(cand_csv, candidates)
    detail_csv = os.path.join(OUT_DIR, "seg75_2m_angle_search_detail.csv")
    write_csv(detail_csv, detail_rows)
    best = candidates[0]
    rows, heatmaps, positions, normals, segment_ids, tilts, centers = refine_candidate(best["a1_deg"], best["a2_deg"], ref_amp)
    rec_csv = os.path.join(OUT_DIR, "seg75_2m_recommended_strategy.csv")
    write_csv(rec_csv, rows)
    summary_png = os.path.join(OUT_DIR, "seg75_2m_recommended_performance.png")
    plot_summary(rows, summary_png)
    geometry_png = os.path.join(OUT_DIR, "seg75_2m_recommended_geometry.png")
    plot_geometry(positions, normals, segment_ids, tilts, centers, geometry_png)
    heatmap_png = os.path.join(OUT_DIR, "seg75_2m_recommended_heatmaps.png")
    plot_heatmaps(heatmaps, heatmap_png)
    summary = os.path.join(OUT_DIR, "seg75_2m_angle_search_summary.txt")
    write_summary(summary, candidates, rows, tilts)
    for path in [summary, cand_csv, detail_csv, rec_csv, summary_png, geometry_png, heatmap_png]:
        print(path)


if __name__ == "__main__":
    main()
