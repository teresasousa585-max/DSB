import csv
import math
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

import design_n12_unit as base


ROOT = os.path.dirname(os.path.dirname(__file__))
OUT_DIR = os.path.join(ROOT, "analysis_outputs")
LAYOUT_CSV = os.path.join(OUT_DIR, "n12_unit_recommended_layout.csv")

CASES = [
    ("broadside", 0.0, 0.0),
    ("steer_x15", 15.0, 0.0),
    ("steer_x30", 30.0, 0.0),
    ("steer_y30", 0.0, 30.0),
    ("steer_x30_y30", 30.0, 30.0),
]
GAIN_FLOORS = [0.9, 0.8, 0.7, 0.6]
MAIN_EXCLUDE_DEG = 11.0
IRLS_ITERS = 8
WEIGHT_POWER = 6.0
REG_TO_UNIFORM = 2e-4
SMOOTH_BETA = 25.0
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


def metrics_for_amplitudes(points_mm, steer_x_deg, steer_y_deg, amp):
    h, target_vec, outside, mask, tx, ty = grid_response_matrix(points_mm, steer_x_deg, steer_y_deg)
    resp = np.dot(h, amp)
    target_amp = float(abs(np.dot(target_vec, amp)))
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


def optimize_amplitudes(points_mm, steer_x_deg, steer_y_deg, gain_floor=0.8):
    h, target_vec, outside, mask, tx, ty = grid_response_matrix(points_mm, steer_x_deg, steer_y_deg)
    side = h[outside].reshape(-1, len(points_mm))
    target_real = np.real(target_vec)
    n = len(points_mm)
    amp = np.ones(n, dtype=float)
    history = []

    for iteration in range(IRLS_ITERS):
        resp = side @ amp
        rel = np.abs(resp) / max(float(np.dot(target_real, amp)), 1e-18)
        rel = rel / max(np.percentile(rel, 90.0), 1e-12)
        sample_weights = np.clip(rel, 0.2, 10.0) ** max(0.0, WEIGHT_POWER - 2.0)
        q = np.real(side.conj().T @ (side * sample_weights[:, None])) / max(len(side), 1)
        q += REG_TO_UNIFORM * np.eye(n)

        def obj(a):
            d = a - 1.0
            return float(a @ q @ a + REG_TO_UNIFORM * d @ d)

        def jac(a):
            return 2.0 * (q @ a + REG_TO_UNIFORM * (a - 1.0))

        cons = [{
            "type": "ineq",
            "fun": lambda a, tr=target_real: float(np.dot(tr, a) - gain_floor * np.dot(tr, np.ones_like(a))),
            "jac": lambda a, tr=target_real: tr,
        }]
        res = minimize(
            obj,
            amp,
            jac=jac,
            bounds=[(0.0, 1.0)] * n,
            constraints=cons,
            method="SLSQP",
            options={"maxiter": 250, "ftol": 1e-11, "disp": False},
        )
        if not res.success:
            amp = np.clip(res.x, 0.0, 1.0)
        else:
            amp = np.clip(res.x, 0.0, 1.0)
        m, _ = metrics_for_amplitudes(points_mm, steer_x_deg, steer_y_deg, amp)
        m["iteration"] = iteration + 1
        history.append(m)
    return amp, history


def optimize_minimax_amplitudes(points_mm, steer_x_deg, steer_y_deg, gain_floor=0.8):
    h, target_vec, outside, mask, tx, ty = grid_response_matrix(points_mm, steer_x_deg, steer_y_deg)
    side = h[outside].reshape(-1, len(points_mm))
    target_real = np.real(target_vec)
    n = len(points_mm)

    def obj(a):
        target = max(float(np.dot(target_real, a)), 1e-18)
        rel = np.abs(side @ a) / target
        m = float(np.max(rel))
        # Smooth max keeps SLSQP stable while still focusing on the hottest sidelobes.
        smooth = m + float(np.log(np.mean(np.exp(SMOOTH_BETA * (rel - m)))) / SMOOTH_BETA)
        return smooth + REG_TO_UNIFORM * float(np.mean((a - 1.0) ** 2))

    cons = [{
        "type": "ineq",
        "fun": lambda a, tr=target_real: float(np.dot(tr, a) - gain_floor * np.dot(tr, np.ones_like(a))),
        "jac": lambda a, tr=target_real: tr,
    }]

    rng = np.random.RandomState(20260511)
    starts = [
        np.ones(n, dtype=float),
        np.linspace(0.65, 1.0, n),
        rng.uniform(0.5, 1.0, n),
    ]
    best_amp = None
    best_metric = None
    for start in starts:
        res = minimize(
            obj,
            start,
            bounds=[(0.0, 1.0)] * n,
            constraints=cons,
            method="SLSQP",
            options={"maxiter": 320, "ftol": 1e-10, "disp": False},
        )
        amp = np.clip(res.x, 0.0, 1.0)
        metric, _ = metrics_for_amplitudes(points_mm, steer_x_deg, steer_y_deg, amp)
        if best_metric is None or metric["max_sidelobe_db"] < best_metric["max_sidelobe_db"]:
            best_metric = metric
            best_amp = amp
    return best_amp


def choose_best(points_mm, steer_x_deg, steer_y_deg):
    trials = []
    for floor in GAIN_FLOORS:
        amp = optimize_minimax_amplitudes(points_mm, steer_x_deg, steer_y_deg, floor)
        m, db = metrics_for_amplitudes(points_mm, steer_x_deg, steer_y_deg, amp)
        m["gain_floor"] = floor
        m["amp"] = amp
        m["db"] = db
        # Prefer real relative sidelobe improvement; do not buy tiny PSL gains with large target loss.
        m["selection_score"] = m["max_sidelobe_db"] + 0.25 * m["p99_sidelobe_db"] + max(0.0, -m["target_gain_loss_db"] - 2.2)
        trials.append(m)
    return min(trials, key=lambda row: row["selection_score"]), trials


def phase_degrees(points_mm, steer_x_deg, steer_y_deg):
    pts_m = np.column_stack([
        points_mm[:, 0] * 1e-3,
        points_mm[:, 1] * 1e-3,
        np.zeros(len(points_mm)),
    ])
    target_dir = direction_from_xy(steer_x_deg, steer_y_deg)
    phase = -base.K * np.dot(pts_m, target_dir)
    phase = (np.rad2deg(phase) + 180.0) % 360.0 - 180.0
    return phase


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
            "phase_deg",
        ]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            phase = phase_degrees(points_mm, row["steer_x_deg"], row["steer_y_deg"])
            for i, (p, amp, ph) in enumerate(zip(points_mm, row["best"]["amp"], phase)):
                writer.writerow({
                    "case": row["case"],
                    "steer_x_deg": row["steer_x_deg"],
                    "steer_y_deg": row["steer_y_deg"],
                    "element_id": "U%02d" % i,
                    "x_mm": "%.3f" % p[0],
                    "y_mm": "%.3f" % p[1],
                    "amplitude": "%.6f" % amp,
                    "phase_deg": "%.3f" % ph,
                })


def write_metrics(path, rows):
    fields = [
        "case",
        "steer_x_deg",
        "steer_y_deg",
        "method",
        "gain_floor",
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
                "gain_floor": "",
                **{k: "%.6f" % uniform[k] for k in fields if k in uniform},
            })
            for trial in row["trials"]:
                writer.writerow({
                    "case": row["case"],
                    "steer_x_deg": row["steer_x_deg"],
                    "steer_y_deg": row["steer_y_deg"],
                "method": "amp_minimax",
                    "gain_floor": "%.2f" % trial["gain_floor"],
                    **{k: "%.6f" % trial[k] for k in fields if k in trial},
                })


def plot_patterns(points_mm, rows, path):
    axis, tx, ty, theta_deg, dirs, elem, mask = base.GRID
    fig, axes = plt.subplots(len(rows), 2, figsize=(10.5, 3.5 * len(rows)))
    levels = np.linspace(DB_FLOOR, 0.0, 56)
    contour = None
    for r, row in enumerate(rows):
        for c, (title, db) in enumerate([
            ("uniform", row["uniform_db"]),
            ("amplitude optimized", row["best"]["db"]),
        ]):
            ax = axes[r, c] if len(rows) > 1 else axes[c]
            plot_db = np.maximum(db, DB_FLOOR)
            contour = ax.contourf(tx, ty, plot_db, levels=levels, cmap="viridis", extend="min")
            ax.plot([row["steer_x_deg"]], [row["steer_y_deg"]], "rx", ms=7, mew=1.6)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlim(-70, 70)
            ax.set_ylim(-70, 70)
            ax.set_title("%s %s\nSL %.1f dB, p99 %.1f dB, gain %.1f dB" % (
                row["case"], title, row[title.replace("amplitude optimized", "best").replace("uniform", "uniform")]["max_sidelobe_db"]
                if title == "uniform" else row["best"]["max_sidelobe_db"],
                row["uniform"]["p99_sidelobe_db"] if title == "uniform" else row["best"]["p99_sidelobe_db"],
                row["uniform"]["target_gain_loss_db"] if title == "uniform" else row["best"]["target_gain_loss_db"],
            ))
            ax.set_xlabel("theta_x (deg)")
            ax.set_ylabel("theta_y (deg)")
            ax.grid(True, color="white", alpha=0.12)
    cbar = fig.colorbar(contour, ax=axes, shrink=0.75)
    cbar.set_label("relative pressure (dB)")
    fig.suptitle("N12 amplitude optimization: uniform vs smooth-minimax")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_weights(rows, path):
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    x = np.arange(12)
    width = 0.14
    for idx, row in enumerate(rows):
        ax.bar(x + (idx - (len(rows) - 1) / 2.0) * width, row["best"]["amp"], width=width, label=row["case"])
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("element id")
    ax.set_ylabel("optimized amplitude")
    ax.set_xticks(x)
    ax.set_xticklabels(["U%02d" % i for i in x])
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(ncol=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary(path, rows):
    with open(path, "w") as f:
        f.write("N12 amplitude optimization report\n")
        f.write("Layout: analysis_outputs/n12_unit_recommended_layout.csv\n")
        f.write("Algorithm: fixed steering phase + non-negative amplitude smooth-minimax/SLSQP, bounds 0..1.\n")
        f.write("Sidelobe region excludes %.1f deg around target; gain floors tested: %s.\n\n" % (
            MAIN_EXCLUDE_DEG, ", ".join("%.2f" % g for g in GAIN_FLOORS)
        ))
        f.write("Per-case comparison:\n")
        for row in rows:
            u = row["uniform"]
            b = row["best"]
            f.write("  %-15s steer=(%5.1f,%5.1f): uniform SL=%6.2f p99=%6.2f; opt SL=%6.2f p99=%6.2f gain=%6.2f meanAmp=%4.2f floor=%.2f\n" % (
                row["case"], row["steer_x_deg"], row["steer_y_deg"],
                u["max_sidelobe_db"], u["p99_sidelobe_db"],
                b["max_sidelobe_db"], b["p99_sidelobe_db"], b["target_gain_loss_db"], b["mean_amp"], b["gain_floor"]
            ))
        f.write("\nOptimized amplitude rows:\n")
        for row in rows:
            f.write("  %s: %s\n" % (row["case"], ", ".join("%.3f" % a for a in row["best"]["amp"])))


def main():
    ensure_out_dir()
    points = read_layout(LAYOUT_CSV)
    rows = []
    for case, sx, sy in CASES:
        uniform_amp = np.ones(len(points), dtype=float)
        uniform_metrics, uniform_db = metrics_for_amplitudes(points, sx, sy, uniform_amp)
        best, trials = choose_best(points, sx, sy)
        rows.append({
            "case": case,
            "steer_x_deg": sx,
            "steer_y_deg": sy,
            "uniform": uniform_metrics,
            "uniform_db": uniform_db,
            "best": best,
            "trials": trials,
        })

    weights_csv = os.path.join(OUT_DIR, "n12_ampopt_weights.csv")
    write_weights(weights_csv, points, rows)
    metrics_csv = os.path.join(OUT_DIR, "n12_ampopt_metrics.csv")
    write_metrics(metrics_csv, rows)
    patterns_png = os.path.join(OUT_DIR, "n12_ampopt_patterns.png")
    plot_patterns(points, rows, patterns_png)
    weights_png = os.path.join(OUT_DIR, "n12_ampopt_weights.png")
    plot_weights(rows, weights_png)
    summary = os.path.join(OUT_DIR, "n12_ampopt_summary.txt")
    write_summary(summary, rows)
    for path in [summary, weights_csv, metrics_csv, patterns_png, weights_png]:
        print(path)


if __name__ == "__main__":
    main()
