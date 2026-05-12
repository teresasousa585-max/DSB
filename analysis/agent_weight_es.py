import csv
import math
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import differential_evolution

import design_n12_unit as base

ROOT = os.path.dirname(os.path.dirname(__file__))
OUT_DIR = os.path.join(ROOT, "analysis_outputs")
LAYOUT_CSV = os.path.join(OUT_DIR, "n12_unit_recommended_layout.csv")

CASES = [
    ("broadside", 0.0, 0.0),
    ("steer_x15", 15.0, 0.0),
    ("steer_x-15", -15.0, 0.0),
    ("steer_y15", 0.0, 15.0),
    ("steer_y-15", 0.0, -15.0),
]
GAIN_FLOOR = 0.8
MAIN_EXCLUDE_DEG = 11.0
SMOOTH_BETA = 25.0
REG_TO_UNIFORM = 1e-4
DB_FLOOR = -55.0


def ensure_out_dir():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)


def read_layout(path):
    points = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            points.append([float(row["X_mm"]), float(row["Y_mm"])])
    return base.centered(np.asarray(points, dtype=float))


def direction_from_xy(theta_x_deg, theta_y_deg):
    sx = math.sin(math.radians(theta_x_deg))
    sy = math.sin(math.radians(theta_y_deg))
    sz = math.sqrt(max(0.0, 1.0 - sx * sx - sy * sy))
    return np.asarray([sx, sy, sz], dtype=float)


def grid_response_matrix(points_mm, steer_x_deg, steer_y_deg):
    axis, tx, ty, theta_deg, dirs, elem, mask = base.GRID
    pts_m = np.column_stack([
        points_mm[:, 0] * 1e-3,
        points_mm[:, 1] * 1e-3,
        np.zeros(len(points_mm)),
    ])
    target_dir = direction_from_xy(steer_x_deg, steer_y_deg)
    phase_to_target = np.exp(-1j * base.K * np.dot(pts_m, target_dir))
    phase_grid = np.exp(1j * base.K * np.tensordot(dirs, pts_m.T, axes=([2], [0])))
    h = elem[..., None] * phase_grid * phase_to_target[None, None, :]
    target_elem = float(base.piston_from_cos(np.asarray([target_dir[2]]))[0])
    target_vec = np.ones(len(points_mm), dtype=complex) * target_elem
    sep = base.angular_separation_deg(tx, ty, steer_x_deg, steer_y_deg)
    outside = (sep > MAIN_EXCLUDE_DEG) & mask
    return h, target_vec, outside, mask, tx, ty


def smooth_max(values, beta=SMOOTH_BETA):
    m = float(np.max(values))
    return m + float(np.log(np.mean(np.exp(beta * (values - m)))) / beta)


def metrics_for_weights(points_mm, steer_x_deg, steer_y_deg, amp, phase_offset):
    h, target_vec, outside, mask, tx, ty = grid_response_matrix(points_mm, steer_x_deg, steer_y_deg)
    weights = amp * np.exp(1j * phase_offset)
    resp = np.dot(h, weights)
    target_amp = float(abs(np.dot(target_vec, weights)))
    full_target_amp = float(abs(np.dot(target_vec, np.ones(len(points_mm)))))
    db = 20.0 * np.log10(np.abs(resp) / max(target_amp, 1e-18) + 1e-12)
    db[~mask] = np.nan
    outside_db = db[outside]
    return {
        "max_sidelobe_db": float(np.nanmax(outside_db)),
        "p99_sidelobe_db": float(np.nanpercentile(outside_db, 99.0)),
        "p95_sidelobe_db": float(np.nanpercentile(outside_db, 95.0)),
        "target_gain_loss_db": 20.0 * math.log10(target_amp / max(full_target_amp, 1e-18) + 1e-12),
        "mean_amp": float(np.mean(amp)),
        "min_amp": float(np.min(amp)),
        "max_amp": float(np.max(amp)),
    }, db


def optimize_weights_de(points_mm, steer_x_deg, steer_y_deg, gain_floor=GAIN_FLOOR):
    h, target_vec, outside, mask, tx, ty = grid_response_matrix(points_mm, steer_x_deg, steer_y_deg)
    side = h[outside].reshape(-1, len(points_mm))
    target_real = np.real(target_vec)
    n = len(points_mm)

    uniform_target = abs(np.dot(target_vec, np.ones(n, dtype=complex)))

    def objective(params):
        amp = params[:n]
        phase_offset = params[n:]
        weights = amp * np.exp(1j * phase_offset)
        target = max(abs(np.dot(target_vec, weights)), 1e-18)
        rel = np.abs(side @ weights) / target
        sm = smooth_max(rel)
        reg = REG_TO_UNIFORM * float(np.mean((amp - 1.0) ** 2))
        # Penalty for gain below floor (relative to uniform excitation)
        gain_penalty = 0.0
        if target < gain_floor * uniform_target:
            gain_penalty = 10.0 * (gain_floor * uniform_target - target) / uniform_target
        return sm + reg + gain_penalty

    bounds = [(0.0, 1.0)] * n + [(-np.pi, np.pi)] * n
    result = differential_evolution(
        objective,
        bounds,
        maxiter=200,
        popsize=15,
        tol=1e-6,
        polish=True,
        seed=20260511 + int(round(steer_x_deg * 100 + steer_y_deg)),
    )
    amp = np.clip(result.x[:n], 0.0, 1.0)
    phase_offset = result.x[n:]
    return amp, phase_offset


def phase_degrees(points_mm, steer_x_deg, steer_y_deg):
    pts_m = np.column_stack([
        points_mm[:, 0] * 1e-3,
        points_mm[:, 1] * 1e-3,
        np.zeros(len(points_mm)),
    ])
    target_dir = direction_from_xy(steer_x_deg, steer_y_deg)
    phase = -base.K * np.dot(pts_m, target_dir)
    return (np.rad2deg(phase) + 180.0) % 360.0 - 180.0


def write_weights(path, points_mm, rows):
    with open(path, "w", newline="") as f:
        fields = [
            "case",
            "steer_x_deg",
            "steer_y_deg",
            "element_id",
            "x_mm",
            "y_mm",
            "amplitude",
            "phase_offset_deg",
            "total_phase_deg",
        ]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            steering_phase = phase_degrees(points_mm, row["steer_x_deg"], row["steer_y_deg"])
            for i, (p, amp, ph_off, ph_total) in enumerate(zip(points_mm, row["amp"], row["phase_offset"], steering_phase + np.rad2deg(row["phase_offset"]))):
                writer.writerow({
                    "case": row["case"],
                    "steer_x_deg": row["steer_x_deg"],
                    "steer_y_deg": row["steer_y_deg"],
                    "element_id": "U%02d" % i,
                    "x_mm": "%.3f" % p[0],
                    "y_mm": "%.3f" % p[1],
                    "amplitude": "%.6f" % amp,
                    "phase_offset_deg": "%.3f" % np.rad2deg(ph_off),
                    "total_phase_deg": "%.3f" % ((ph_total + 180.0) % 360.0 - 180.0),
                })


def write_metrics(path, rows):
    fields = [
        "case",
        "steer_x_deg",
        "steer_y_deg",
        "method",
        "max_sidelobe_db",
        "p99_sidelobe_db",
        "p95_sidelobe_db",
        "target_gain_loss_db",
        "mean_amp",
        "min_amp",
        "max_amp",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            uniform = row["uniform"]
            writer.writerow({
                "case": row["case"],
                "steer_x_deg": row["steer_x_deg"],
                "steer_y_deg": row["steer_y_deg"],
                "method": "uniform",
                **{k: "%.6f" % uniform[k] for k in fields if k in uniform},
            })
            opt = row["opt"]
            writer.writerow({
                "case": row["case"],
                "steer_x_deg": row["steer_x_deg"],
                "steer_y_deg": row["steer_y_deg"],
                "method": "phase_amp_de",
                **{k: "%.6f" % opt[k] for k in fields if k in opt},
            })


def plot_patterns(points_mm, rows, path):
    axis, tx, ty, theta_deg, dirs, elem, mask = base.GRID
    fig, axes = plt.subplots(len(rows), 2, figsize=(10.5, 3.5 * len(rows)))
    levels = np.linspace(DB_FLOOR, 0.0, 56)
    contour = None
    for r, row in enumerate(rows):
        for c, (title, db) in enumerate([
            ("uniform", row["uniform_db"]),
            ("phase+amp optimized", row["opt_db"]),
        ]):
            ax = axes[r, c] if len(rows) > 1 else axes[c]
            plot_db = np.maximum(db, DB_FLOOR)
            contour = ax.contourf(tx, ty, plot_db, levels=levels, cmap="viridis", extend="min")
            ax.plot([row["steer_x_deg"]], [row["steer_y_deg"]], "rx", ms=7, mew=1.6)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlim(-70, 70)
            ax.set_ylim(-70, 70)
            data = row["uniform"] if title == "uniform" else row["opt"]
            ax.set_title("%s %s\nSL %.1f dB, p99 %.1f dB, gain %.1f dB" % (
                row["case"], title,
                data["max_sidelobe_db"],
                data["p99_sidelobe_db"],
                data["target_gain_loss_db"],
            ))
            ax.set_xlabel("theta_x (deg)")
            ax.set_ylabel("theta_y (deg)")
            ax.grid(True, color="white", alpha=0.12)
    cbar = fig.colorbar(contour, ax=axes, shrink=0.75)
    cbar.set_label("relative pressure (dB)")
    fig.suptitle("N12 phase + amplitude optimization (DE)")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_weights(rows, path):
    fig, axes = plt.subplots(2, 1, figsize=(9.0, 8.0))
    x = np.arange(12)
    width = 0.14
    for idx, row in enumerate(rows):
        axes[0].bar(x + (idx - (len(rows) - 1) / 2.0) * width, row["amp"], width=width, label=row["case"])
        axes[1].bar(x + (idx - (len(rows) - 1) / 2.0) * width, np.rad2deg(row["phase_offset"]), width=width, label=row["case"])
    axes[0].set_ylim(0, 1.05)
    axes[0].set_ylabel("optimized amplitude")
    axes[1].set_ylabel("phase offset (deg)")
    for ax in axes:
        ax.set_xlabel("element id")
        ax.set_xticks(x)
        ax.set_xticklabels(["U%02d" % i for i in x])
        ax.grid(True, axis="y", alpha=0.25)
        ax.legend(ncol=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary(path, rows):
    with open(path, "w") as f:
        f.write("N12 phase + amplitude optimization report (DE)\n")
        f.write("Layout: analysis_outputs/n12_unit_recommended_layout.csv\n")
        f.write("Algorithm: steering phase + phase offset + amplitude, optimized via differential evolution.\n")
        f.write("Bounds: amp in [0,1], phase offset in [-pi, pi].\n\n")
        for row in rows:
            u = row["uniform"]
            o = row["opt"]
            f.write("  %-15s steer=(%5.1f,%5.1f): uniform SL=%6.2f p99=%6.2f; opt SL=%6.2f p99=%6.2f gain=%6.2f meanAmp=%4.2f\n" % (
                row["case"], row["steer_x_deg"], row["steer_y_deg"],
                u["max_sidelobe_db"], u["p99_sidelobe_db"],
                o["max_sidelobe_db"], o["p99_sidelobe_db"], o["target_gain_loss_db"], o["mean_amp"]
            ))
        f.write("\nOptimized amplitudes:\n")
        for row in rows:
            f.write("  %s: %s\n" % (row["case"], ", ".join("%.3f" % a for a in row["amp"])))
        f.write("\nOptimized phase offsets (deg):\n")
        for row in rows:
            f.write("  %s: %s\n" % (row["case"], ", ".join("%.2f" % np.rad2deg(p) for p in row["phase_offset"])))


def main():
    ensure_out_dir()
    points = read_layout(LAYOUT_CSV)
    rows = []
    for case, sx, sy in CASES:
        print("Optimizing %s (%.1f, %.1f) ..." % (case, sx, sy))
        uniform_amp = np.ones(len(points), dtype=float)
        uniform_phase_offset = np.zeros(len(points), dtype=float)
        uniform_metrics, uniform_db = metrics_for_weights(points, sx, sy, uniform_amp, uniform_phase_offset)
        amp, phase_offset = optimize_weights_de(points, sx, sy)
        opt_metrics, opt_db = metrics_for_weights(points, sx, sy, amp, phase_offset)
        rows.append({
            "case": case,
            "steer_x_deg": sx,
            "steer_y_deg": sy,
            "uniform": uniform_metrics,
            "uniform_db": uniform_db,
            "opt": opt_metrics,
            "opt_db": opt_db,
            "amp": amp,
            "phase_offset": phase_offset,
        })
        print("  uniform SL=%.2f p99=%.2f | opt SL=%.2f p99=%.2f gain=%.2f dB" % (
            uniform_metrics["max_sidelobe_db"], uniform_metrics["p99_sidelobe_db"],
            opt_metrics["max_sidelobe_db"], opt_metrics["p99_sidelobe_db"],
            opt_metrics["target_gain_loss_db"],
        ))

    weights_csv = os.path.join(OUT_DIR, "agent_weight_es_weights.csv")
    write_weights(weights_csv, points, rows)
    metrics_csv = os.path.join(OUT_DIR, "agent_weight_es_metrics.csv")
    write_metrics(metrics_csv, rows)
    patterns_png = os.path.join(OUT_DIR, "agent_weight_es_patterns.png")
    plot_patterns(points, rows, patterns_png)
    weights_png = os.path.join(OUT_DIR, "agent_weight_es_weights.png")
    plot_weights(rows, weights_png)
    summary = os.path.join(OUT_DIR, "agent_weight_es_summary.txt")
    write_summary(summary, rows)
    for path in [summary, weights_csv, metrics_csv, patterns_png, weights_png]:
        print(path)


if __name__ == "__main__":
    main()
