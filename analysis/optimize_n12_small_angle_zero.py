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

CASES = [5.0, 10.0, 15.0]
GAIN_FLOORS = [0.9, 0.8, 0.7, 0.6]
TARGET_EXCLUDE_DEG = 6.0
ZERO_REGION_DEG = 2.5
ZERO_WEIGHT = 2.0
SIDE_WEIGHT = 0.55
SMOOTH_BETA = 30.0
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


def grid_response_matrix(points_mm, steer_x_deg, steer_y_deg=0.0):
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
    target_sep = base.angular_separation_deg(tx, ty, steer_x_deg, steer_y_deg)
    zero_sep = base.angular_separation_deg(tx, ty, 0.0, 0.0)
    outside_target = (target_sep > TARGET_EXCLUDE_DEG) & mask
    zero_region = (zero_sep <= ZERO_REGION_DEG) & mask
    return h, target_vec, outside_target, zero_region, mask, tx, ty


def response_at_angle(points_mm, steer_x_deg, amp, probe_x_deg, probe_y_deg=0.0):
    pts_m = np.column_stack([
        points_mm[:, 0] * 1e-3,
        points_mm[:, 1] * 1e-3,
        np.zeros(len(points_mm)),
    ])
    target_dir = direction_from_xy(steer_x_deg, 0.0)
    probe_dir = direction_from_xy(probe_x_deg, probe_y_deg)
    phase_to_target = np.exp(-1j * base.K * np.dot(pts_m, target_dir))
    phase_probe = np.exp(1j * base.K * np.dot(pts_m, probe_dir))
    elem = float(base.piston_from_cos(np.asarray([probe_dir[2]]))[0])
    target_elem = float(base.piston_from_cos(np.asarray([target_dir[2]]))[0])
    resp = elem * np.dot(phase_probe * phase_to_target, amp)
    target_resp = target_elem * np.dot(np.ones(len(points_mm)), amp)
    return 20.0 * math.log10(abs(resp) / max(abs(target_resp), 1e-18) + 1e-12)


def metrics(points_mm, steer_x_deg, amp):
    h, target_vec, outside_target, zero_region, mask, tx, ty = grid_response_matrix(points_mm, steer_x_deg)
    resp = np.dot(h, amp)
    target_amp = float(abs(np.dot(target_vec, amp)))
    full_target_amp = float(abs(np.dot(target_vec, np.ones(len(points_mm)))))
    db = 20.0 * np.log10(np.abs(resp) / max(target_amp, 1e-18) + 1e-12)
    db[~mask] = np.nan
    out_db = db[outside_target]
    zero_db = db[zero_region]
    return {
        "max_outside_target_db": float(np.nanmax(out_db)),
        "p99_outside_target_db": float(np.nanpercentile(out_db, 99.0)),
        "zero_point_db": response_at_angle(points_mm, steer_x_deg, amp, 0.0, 0.0),
        "zero_region_max_db": float(np.nanmax(zero_db)),
        "zero_region_p99_db": float(np.nanpercentile(zero_db, 99.0)),
        "target_gain_loss_db": 20.0 * math.log10(target_amp / max(full_target_amp, 1e-18) + 1e-12),
        "mean_amp": float(np.mean(amp)),
        "min_amp": float(np.min(amp)),
        "max_amp": float(np.max(amp)),
    }, db


def optimize_zero_suppression(points_mm, steer_x_deg, gain_floor):
    h, target_vec, outside_target, zero_region, mask, tx, ty = grid_response_matrix(points_mm, steer_x_deg)
    side = h[outside_target].reshape(-1, len(points_mm))
    zero = h[zero_region].reshape(-1, len(points_mm))
    target_real = np.real(target_vec)
    n = len(points_mm)

    def smooth_max(values):
        m = float(np.max(values))
        return m + float(np.log(np.mean(np.exp(SMOOTH_BETA * (values - m)))) / SMOOTH_BETA)

    def obj(a):
        target = max(float(np.dot(target_real, a)), 1e-18)
        side_rel = SIDE_WEIGHT * np.abs(side @ a) / target
        zero_rel = ZERO_WEIGHT * np.abs(zero @ a) / target
        vals = np.concatenate([side_rel, zero_rel])
        return smooth_max(vals) + REG_TO_UNIFORM * float(np.mean((a - 1.0) ** 2))

    cons = [{
        "type": "ineq",
        "fun": lambda a, tr=target_real: float(np.dot(tr, a) - gain_floor * np.dot(tr, np.ones_like(a))),
        "jac": lambda a, tr=target_real: tr,
    }]
    rng = np.random.RandomState(20260511 + int(round(steer_x_deg * 10)))
    starts = [
        np.ones(n, dtype=float),
        np.linspace(1.0, 0.6, n),
        np.linspace(0.6, 1.0, n),
        rng.uniform(0.45, 1.0, n),
    ]
    best_amp = None
    best_score = None
    for start in starts:
        res = minimize(
            obj,
            start,
            bounds=[(0.0, 1.0)] * n,
            constraints=cons,
            method="SLSQP",
            options={"maxiter": 360, "ftol": 1e-10, "disp": False},
        )
        amp = np.clip(res.x, 0.0, 1.0)
        m, _ = metrics(points_mm, steer_x_deg, amp)
        score = m["zero_region_max_db"] + 0.15 * m["max_outside_target_db"] + max(0.0, -m["target_gain_loss_db"] - 2.2)
        if best_score is None or score < best_score:
            best_score = score
            best_amp = amp
    return best_amp


def choose_best(points_mm, steer_x_deg):
    trials = []
    for floor in GAIN_FLOORS:
        amp = optimize_zero_suppression(points_mm, steer_x_deg, floor)
        m, db = metrics(points_mm, steer_x_deg, amp)
        m["gain_floor"] = floor
        m["amp"] = amp
        m["db"] = db
        # Keep a practical gain limit unless zero suppression becomes much better.
        m["selection_score"] = m["zero_region_max_db"] + 0.10 * m["max_outside_target_db"] + max(0.0, -m["target_gain_loss_db"] - 2.2)
        trials.append(m)
    return min(trials, key=lambda r: r["selection_score"]), trials


def phase_degrees(points_mm, steer_x_deg):
    pts_m = np.column_stack([
        points_mm[:, 0] * 1e-3,
        points_mm[:, 1] * 1e-3,
        np.zeros(len(points_mm)),
    ])
    target_dir = direction_from_xy(steer_x_deg, 0.0)
    phase = -base.K * np.dot(pts_m, target_dir)
    return (np.rad2deg(phase) + 180.0) % 360.0 - 180.0


def linecut(points_mm, steer_x_deg, amp):
    xs = np.linspace(-25.0, 30.0, 551)
    pts_m = np.column_stack([
        points_mm[:, 0] * 1e-3,
        points_mm[:, 1] * 1e-3,
        np.zeros(len(points_mm)),
    ])
    target_dir = direction_from_xy(steer_x_deg, 0.0)
    phase_to_target = np.exp(-1j * base.K * np.dot(pts_m, target_dir))
    target_elem = float(base.piston_from_cos(np.asarray([target_dir[2]]))[0])
    target_resp = target_elem * np.dot(np.ones(len(points_mm)), amp)
    vals = []
    for x in xs:
        probe = direction_from_xy(float(x), 0.0)
        elem = float(base.piston_from_cos(np.asarray([probe[2]]))[0])
        phase_probe = np.exp(1j * base.K * np.dot(pts_m, probe))
        vals.append(20.0 * math.log10(abs(elem * np.dot(phase_probe * phase_to_target, amp)) / max(abs(target_resp), 1e-18) + 1e-12))
    return xs, np.asarray(vals)


def write_metrics(path, rows):
    fields = [
        "case",
        "steer_x_deg",
        "method",
        "gain_floor",
        "zero_point_db",
        "zero_region_max_db",
        "zero_region_p99_db",
        "max_outside_target_db",
        "p99_outside_target_db",
        "target_gain_loss_db",
        "mean_amp",
        "min_amp",
        "max_amp",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            for method, data in [("uniform", row["uniform"])]:
                writer.writerow({
                    "case": row["case"],
                    "steer_x_deg": row["steer_x_deg"],
                    "method": method,
                    "gain_floor": "",
                    **{k: "%.6f" % data[k] for k in fields if k in data},
                })
            for trial in row["trials"]:
                writer.writerow({
                    "case": row["case"],
                    "steer_x_deg": row["steer_x_deg"],
                    "method": "zero_suppressed",
                    "gain_floor": "%.2f" % trial["gain_floor"],
                    **{k: "%.6f" % trial[k] for k in fields if k in trial},
                })


def write_weights(path, points_mm, rows):
    fields = ["case", "steer_x_deg", "element_id", "x_mm", "y_mm", "amplitude", "phase_deg"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            phase = phase_degrees(points_mm, row["steer_x_deg"])
            for i, (p, amp, ph) in enumerate(zip(points_mm, row["best"]["amp"], phase)):
                writer.writerow({
                    "case": row["case"],
                    "steer_x_deg": "%.1f" % row["steer_x_deg"],
                    "element_id": "U%02d" % i,
                    "x_mm": "%.3f" % p[0],
                    "y_mm": "%.3f" % p[1],
                    "amplitude": "%.6f" % amp,
                    "phase_deg": "%.3f" % ph,
                })


def plot_linecuts(points_mm, rows, path):
    fig, axes = plt.subplots(len(rows), 1, figsize=(8.0, 8.2), sharex=True)
    for ax, row in zip(axes, rows):
        xu, du = linecut(points_mm, row["steer_x_deg"], np.ones(len(points_mm)))
        xo, do = linecut(points_mm, row["steer_x_deg"], row["best"]["amp"])
        ax.plot(xu, du, label="uniform")
        ax.plot(xo, do, label="zero suppressed")
        ax.axvline(0.0, color="k", ls=":", lw=1.0)
        ax.axvline(row["steer_x_deg"], color="r", ls=":", lw=1.0)
        ax.set_ylim(-45, 5)
        ax.set_ylabel("rel. pressure (dB)")
        ax.set_title("%s: target %.1f deg, zero %.1f -> %.1f dB" % (
            row["case"],
            row["steer_x_deg"],
            row["uniform"]["zero_point_db"],
            row["best"]["zero_point_db"],
        ))
        ax.grid(True, alpha=0.25)
        ax.legend()
    axes[-1].set_xlabel("theta_x (deg), theta_y=0")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_patterns(rows, path):
    axis, tx, ty, theta_deg, dirs, elem, mask = base.GRID
    fig, axes = plt.subplots(len(rows), 2, figsize=(10.5, 3.8 * len(rows)))
    levels = np.linspace(DB_FLOOR, 0.0, 56)
    contour = None
    for r, row in enumerate(rows):
        for c, (label, db) in enumerate([("uniform", row["uniform_db"]), ("zero suppressed", row["best"]["db"])]):
            ax = axes[r, c]
            plot_db = np.maximum(db, DB_FLOOR)
            contour = ax.contourf(tx, ty, plot_db, levels=levels, cmap="viridis", extend="min")
            ax.plot([row["steer_x_deg"]], [0.0], "rx", ms=7, mew=1.6)
            ax.plot([0.0], [0.0], "wo", ms=4)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlim(-45, 45)
            ax.set_ylim(-45, 45)
            data = row["uniform"] if label == "uniform" else row["best"]
            ax.set_title("%s %s\nzero %.1f dB, max %.1f dB, gain %.1f dB" % (
                row["case"], label, data["zero_region_max_db"], data["max_outside_target_db"], data["target_gain_loss_db"]
            ))
            ax.set_xlabel("theta_x (deg)")
            ax.set_ylabel("theta_y (deg)")
            ax.grid(True, color="white", alpha=0.12)
    cbar = fig.colorbar(contour, ax=axes, shrink=0.75)
    cbar.set_label("relative pressure (dB)")
    fig.suptitle("N12 small-angle steering with 0 deg suppression")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary(path, rows):
    with open(path, "w") as f:
        f.write("N12 small-angle 0-degree suppression report\n")
        f.write("Cases: steering only to 5/10/15 deg in x direction.\n")
        f.write("Objective: fixed steering phase, non-negative amplitudes 0..1, suppress +/-%.1f deg around 0 deg.\n" % ZERO_REGION_DEG)
        f.write("Target main exclusion for general sidelobe metric: %.1f deg.\n\n" % TARGET_EXCLUDE_DEG)
        for row in rows:
            u = row["uniform"]
            b = row["best"]
            f.write("  %-10s: zero point %6.2f -> %6.2f dB, zero region %6.2f -> %6.2f dB, max outside target %6.2f -> %6.2f dB, gain=%6.2f dB, floor=%.2f\n" % (
                row["case"],
                u["zero_point_db"], b["zero_point_db"],
                u["zero_region_max_db"], b["zero_region_max_db"],
                u["max_outside_target_db"], b["max_outside_target_db"],
                b["target_gain_loss_db"], b["gain_floor"],
            ))
        f.write("\nSelected amplitudes:\n")
        for row in rows:
            f.write("  %s: %s\n" % (row["case"], ", ".join("%.3f" % a for a in row["best"]["amp"])))


def main():
    ensure_out_dir()
    points = read_layout(LAYOUT_CSV)
    rows = []
    for steer in CASES:
        uniform_amp = np.ones(len(points), dtype=float)
        uniform_metric, uniform_db = metrics(points, steer, uniform_amp)
        best, trials = choose_best(points, steer)
        rows.append({
            "case": "steer_x%.0f" % steer,
            "steer_x_deg": steer,
            "uniform": uniform_metric,
            "uniform_db": uniform_db,
            "best": best,
            "trials": trials,
        })
    metrics_csv = os.path.join(OUT_DIR, "n12_smallangle_zero_metrics.csv")
    write_metrics(metrics_csv, rows)
    weights_csv = os.path.join(OUT_DIR, "n12_smallangle_zero_weights.csv")
    write_weights(weights_csv, points, rows)
    line_png = os.path.join(OUT_DIR, "n12_smallangle_zero_linecuts.png")
    plot_linecuts(points, rows, line_png)
    patterns_png = os.path.join(OUT_DIR, "n12_smallangle_zero_patterns.png")
    plot_patterns(rows, patterns_png)
    summary = os.path.join(OUT_DIR, "n12_smallangle_zero_summary.txt")
    write_summary(summary, rows)
    for path in [summary, metrics_csv, weights_csv, line_png, patterns_png]:
        print(path)


if __name__ == "__main__":
    main()
