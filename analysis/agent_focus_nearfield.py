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

C = 343.0
FREQ = 40000.0
LAMBDA = C / FREQ
K = 2.0 * np.pi / LAMBDA
ELEMENT_DIAMETER_MM = 16.0
ELEMENT_RADIUS_M = ELEMENT_DIAMETER_MM * 0.5e-3
MIN_CENTER_SPACING_MM = 17.2
DB_FLOOR = -55.0

# Focus cases: (distance_m, theta_x_deg, theta_y_deg)
FOCUS_CASES = [
    (0.2, 0.0, 0.0),
    (0.2, 15.0, 0.0),
    (0.2, -15.0, 0.0),
    (0.2, 0.0, 15.0),
    (0.5, 0.0, 0.0),
    (0.5, 15.0, 0.0),
    (0.5, -15.0, 0.0),
    (0.5, 0.0, 15.0),
    (0.8, 0.0, 0.0),
    (0.8, 15.0, 0.0),
    (0.8, -15.0, 0.0),
    (0.8, 0.0, 15.0),
    (1.0, 0.0, 0.0),
    (1.0, 15.0, 0.0),
    (1.0, -15.0, 0.0),
    (1.0, 0.0, 15.0),
]

MAIN_EXCLUDE_DEG = 3.0
SMOOTH_BETA = 25.0
REG_TO_UNIFORM = 1e-4
GAIN_FLOOR = 0.7


def ensure_out_dir():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)


def read_layout(path):
    points = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            points.append([float(row["X_mm"]), float(row["Y_mm"])])
    return base.centered(np.asarray(points, dtype=float))


def focal_point_xyz(distance_m, theta_x_deg, theta_y_deg):
    sx = math.sin(math.radians(theta_x_deg))
    sy = math.sin(math.radians(theta_y_deg))
    sz2 = max(0.0, 1.0 - sx * sx - sy * sy)
    sz = math.sqrt(sz2)
    return np.asarray([distance_m * sx, distance_m * sy, distance_m * sz])


def probe_grid_sphere(distance_m, theta_x_deg, theta_y_deg, span_deg=20.0, n=41):
    """Generate probe points on a sphere at fixed distance around focal direction."""
    tx = np.linspace(-span_deg, span_deg, n)
    ty = np.linspace(-span_deg, span_deg, n)
    txs, tys = np.meshgrid(tx, ty)
    # Convert offsets to absolute direction angles (approximate for small offsets)
    abs_tx = theta_x_deg + txs
    abs_ty = theta_y_deg + tys
    sx = np.sin(np.deg2rad(abs_tx))
    sy = np.sin(np.deg2rad(abs_ty))
    sz = np.sqrt(np.maximum(0.0, 1.0 - sx * sx - sy * sy))
    x = distance_m * sx
    y = distance_m * sy
    z = distance_m * sz
    return x, y, z, txs, tys


def angular_separation_from_focal(txs_probe, tys_probe):
    """Approximate angular separation for small offsets."""
    return np.sqrt(txs_probe ** 2 + tys_probe ** 2)


def element_pattern_to_point(cos_theta):
    """Piston element pattern (far-field approx) vs cos(theta)."""
    return base.piston_from_cos(cos_theta)


def build_focus_matrix(points_mm, distance_m, theta_x_deg, theta_y_deg):
    """
    Build focused response matrix h_focused[probe, element].
    Pressure at probe p with weights w is: sum_n h_focused[p,n] * w[n]
    where w[n] = amp[n] * exp(1j * phase_offset[n]).
    The focusing phase (steering to focal point) is already folded in.
    """
    pts_m = np.column_stack([
        points_mm[:, 0] * 1e-3,
        points_mm[:, 1] * 1e-3,
        np.zeros(len(points_mm)),
    ])
    focal = focal_point_xyz(distance_m, theta_x_deg, theta_y_deg)

    # Probe points on sphere
    x, y, z, txs, tys = probe_grid_sphere(distance_m, theta_x_deg, theta_y_deg)
    probe_shape = x.shape
    n_probes = x.size
    x_flat = x.ravel()
    y_flat = y.ravel()
    z_flat = z.ravel()
    probes = np.column_stack([x_flat, y_flat, z_flat])

    # Distance from each element to each probe: [n_probes, n_elements]
    dx = probes[:, 0][:, None] - pts_m[:, 0][None, :]
    dy = probes[:, 1][:, None] - pts_m[:, 1][None, :]
    dz = probes[:, 2][:, None] - pts_m[:, 2][None, :]
    dist = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

    # Cosine of angle between element normal (+z) and direction to probe
    cos_theta = dz / dist
    elem = element_pattern_to_point(cos_theta)

    # Propagation phase and 1/r amplitude
    h = elem * np.exp(1j * K * dist) / dist

    # Focusing phase: cancel phase at focal point
    # dist_to_focal[n_elements]
    dx_f = focal[0] - pts_m[:, 0]
    dy_f = focal[1] - pts_m[:, 1]
    dz_f = focal[2] - pts_m[:, 2]
    dist_f = np.sqrt(dx_f ** 2 + dy_f ** 2 + dz_f ** 2)
    focus_phase = np.exp(-1j * K * dist_f)

    # Apply focusing phase to all probes
    h_focused = h * focus_phase[None, :]

    # Focal point response (for normalization)
    focal_cos = dz_f / dist_f
    focal_elem = element_pattern_to_point(np.asarray([focal_cos]))[0]
    focal_resp_per_elem = focal_elem * np.exp(1j * K * dist_f) / dist_f * focus_phase
    # This simplifies to: focal_elem / dist_f (real, positive)
    focal_resp = np.sum(focal_resp_per_elem)  # Should be real and positive

    # Angular separation for sidelobe mask
    sep = angular_separation_from_focal(txs, tys).ravel()
    outside = sep > MAIN_EXCLUDE_DEG

    return h_focused, focal_resp, outside, probe_shape, txs, tys


def smooth_max(values, beta=SMOOTH_BETA):
    m = float(np.max(values))
    return m + float(np.log(np.mean(np.exp(beta * (values - m)))) / beta)


def optimize_focus_weights(points_mm, distance_m, theta_x_deg, theta_y_deg):
    h_focused, focal_resp_uniform, outside, probe_shape, txs, tys = build_focus_matrix(
        points_mm, distance_m, theta_x_deg, theta_y_deg
    )
    n = len(points_mm)
    uniform_target = float(abs(focal_resp_uniform))

    def objective(params):
        amp = params[:n]
        phase_offset = params[n:]
        weights = amp * np.exp(1j * phase_offset)
        resp = h_focused @ weights
        center_idx = (probe_shape[0] // 2) * probe_shape[1] + probe_shape[1] // 2
        focal_resp = resp[center_idx]
        target = max(abs(focal_resp), 1e-18)

        side_resp = resp[outside]
        # Objective: maximize dBc = 20*log10(target / max(side))
        # DE minimizes, so minimize -dBc = 20*log10(max(side) / target)
        rel = np.abs(side_resp) / target
        sm = smooth_max(rel)
        reg = REG_TO_UNIFORM * float(np.mean((amp - 1.0) ** 2))
        return sm + reg

    bounds = [(0.0, 1.0)] * n + [(-np.pi, np.pi)] * n
    result = differential_evolution(
        objective,
        bounds,
        maxiter=120,
        popsize=12,
        tol=1e-6,
        polish=True,
        seed=20260511 + int(round(distance_m * 1000 + theta_x_deg * 10 + theta_y_deg)),
    )
    amp = np.clip(result.x[:n], 0.0, 1.0)
    phase_offset = result.x[n:]

    # Compute metrics
    weights = amp * np.exp(1j * phase_offset)
    resp = h_focused @ weights
    center_idx = (probe_shape[0] // 2) * probe_shape[1] + probe_shape[1] // 2
    focal_resp = resp[center_idx]
    target = max(abs(focal_resp), 1e-18)

    db = 20.0 * np.log10(np.abs(resp) / target + 1e-12)
    db = db.reshape(probe_shape)
    outside_db = db.ravel()[outside]

    # dBc = mainlobe_dB - max_sidelobe_dB (using current mainlobe as 0 dB reference)
    max_sl_db = float(np.nanmax(outside_db))
    dBc = -max_sl_db
    metrics = {
        "max_sidelobe_db": max_sl_db,
        "p99_sidelobe_db": float(np.nanpercentile(outside_db, 99.0)),
        "p95_sidelobe_db": float(np.nanpercentile(outside_db, 95.0)),
        "target_gain_loss_db": 20.0 * math.log10(target / max(uniform_target, 1e-18) + 1e-12),
        "dBc": dBc,
        "mean_amp": float(np.mean(amp)),
        "min_amp": float(np.min(amp)),
        "max_amp": float(np.max(amp)),
    }
    return amp, phase_offset, metrics, db, txs, tys


def write_focus_weights(path, points_mm, rows):
    with open(path, "w", newline="") as f:
        fields = [
            "case",
            "distance_m",
            "theta_x_deg",
            "theta_y_deg",
            "element_id",
            "x_mm",
            "y_mm",
            "amplitude",
            "phase_offset_deg",
        ]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            for i, (p, amp, ph_off) in enumerate(zip(points_mm, row["amp"], row["phase_offset"])):
                writer.writerow({
                    "case": row["case"],
                    "distance_m": row["distance_m"],
                    "theta_x_deg": row["theta_x_deg"],
                    "theta_y_deg": row["theta_y_deg"],
                    "element_id": "U%02d" % i,
                    "x_mm": "%.3f" % p[0],
                    "y_mm": "%.3f" % p[1],
                    "amplitude": "%.6f" % amp,
                    "phase_offset_deg": "%.3f" % np.rad2deg(ph_off),
                })


def write_focus_metrics(path, rows):
    fields = [
        "case",
        "distance_m",
        "theta_x_deg",
        "theta_y_deg",
        "method",
        "max_sidelobe_db",
        "p99_sidelobe_db",
        "p95_sidelobe_db",
        "target_gain_loss_db",
        "dBc",
        "mean_amp",
        "min_amp",
        "max_amp",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            for method, data in [("uniform", row["uniform"]), ("focus_opt", row["opt"])]:
                writer.writerow({
                    "case": row["case"],
                    "distance_m": row["distance_m"],
                    "theta_x_deg": row["theta_x_deg"],
                    "theta_y_deg": row["theta_y_deg"],
                    "method": method,
                    **{k: "%.6f" % data[k] for k in fields if k in data},
                })


def plot_focus_patterns(rows, path):
    fig, axes = plt.subplots(len(rows), 2, figsize=(11.0, 3.8 * len(rows)))
    levels = np.linspace(DB_FLOOR, 0.0, 56)
    contour = None
    for r, row in enumerate(rows):
        for c, (label, db) in enumerate([("uniform", row["uniform_db"]), ("focused", row["opt_db"])]):
            ax = axes[r, c] if len(rows) > 1 else axes[c]
            plot_db = np.maximum(db, DB_FLOOR)
            contour = ax.contourf(row["txs"], row["tys"], plot_db, levels=levels, cmap="viridis", extend="min")
            ax.plot([0.0], [0.0], "rx", ms=7, mew=1.6)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlim(-20, 20)
            ax.set_ylim(-20, 20)
            data = row["uniform"] if label == "uniform" else row["opt"]
            ax.set_title("%s %s @ %.1fm\nSL %.1f dB, p99 %.1f dB, gain %.1f dB" % (
                row["case"], label, row["distance_m"],
                data["max_sidelobe_db"], data["p99_sidelobe_db"], data["target_gain_loss_db"],
            ))
            ax.set_xlabel("offset theta_x (deg)")
            ax.set_ylabel("offset theta_y (deg)")
            ax.grid(True, color="white", alpha=0.12)
    cbar = fig.colorbar(contour, ax=axes, shrink=0.75)
    cbar.set_label("relative pressure (dB)")
    fig.suptitle("N12 near-field focus optimization (ring_1_4_7)")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_focus_summary(path, rows):
    with open(path, "w") as f:
        f.write("N12 near-field focus optimization report\n")
        f.write("Layout: ring_1_4_7 (agent_parametric_best.csv)\n")
        f.write("Algorithm: differential evolution on amplitude + phase offset.\n")
        f.write("Near-field model: exact spherical propagation + piston element pattern.\n")
        f.write("Focus region exclusion: %.1f deg around focal direction.\n\n" % MAIN_EXCLUDE_DEG)
        f.write("%-12s %6s %6s %6s %8s %8s %8s %8s %8s %8s\n" % (
            "case", "dist", "tx", "ty", "u_dBc", "u_gain", "o_dBc", "o_gain", "o_SL", "meanAmp"
        ))
        for row in rows:
            u = row["uniform"]
            o = row["opt"]
            f.write("%-12s %6.1f %6.1f %6.1f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n" % (
                row["case"], row["distance_m"], row["theta_x_deg"], row["theta_y_deg"],
                -u["max_sidelobe_db"], u["target_gain_loss_db"],
                o["dBc"], o["target_gain_loss_db"], o["max_sidelobe_db"], o["mean_amp"]
            ))


def main():
    ensure_out_dir()
    points = read_layout(os.path.join(OUT_DIR, "agent_parametric_best.csv"))
    rows = []

    for dist, tx, ty in FOCUS_CASES:
        case_name = "R%.1f_x%g_y%g" % (dist, tx, ty)
        print("Optimizing %s ..." % case_name)

        # Uniform weights baseline
        h_focused, focal_resp_uniform, outside, probe_shape, txs, tys = build_focus_matrix(
            points, dist, tx, ty
        )
        uniform_weights = np.ones(len(points), dtype=complex)
        resp_uniform = h_focused @ uniform_weights
        center_idx = (probe_shape[0] // 2) * probe_shape[1] + probe_shape[1] // 2
        target_uniform = max(abs(resp_uniform[center_idx]), 1e-18)
        db_uniform = 20.0 * np.log10(np.abs(resp_uniform) / target_uniform + 1e-12).reshape(probe_shape)
        outside_db_uniform = db_uniform.ravel()[outside]
        uniform_metrics = {
            "max_sidelobe_db": float(np.nanmax(outside_db_uniform)),
            "p99_sidelobe_db": float(np.nanpercentile(outside_db_uniform, 99.0)),
            "p95_sidelobe_db": float(np.nanpercentile(outside_db_uniform, 95.0)),
            "target_gain_loss_db": 0.0,
            "mean_amp": 1.0,
            "min_amp": 1.0,
            "max_amp": 1.0,
        }

        # Optimize
        amp, phase_offset, opt_metrics, opt_db, _, _ = optimize_focus_weights(points, dist, tx, ty)

        rows.append({
            "case": case_name,
            "distance_m": dist,
            "theta_x_deg": tx,
            "theta_y_deg": ty,
            "uniform": uniform_metrics,
            "uniform_db": db_uniform,
            "opt": opt_metrics,
            "opt_db": opt_db,
            "txs": txs,
            "tys": tys,
            "amp": amp,
            "phase_offset": phase_offset,
        })
        print("  uniform SL=%.2f p99=%.2f | opt SL=%.2f p99=%.2f gain=%.2f dB" % (
            uniform_metrics["max_sidelobe_db"], uniform_metrics["p99_sidelobe_db"],
            opt_metrics["max_sidelobe_db"], opt_metrics["p99_sidelobe_db"], opt_metrics["target_gain_loss_db"],
        ))

    weights_csv = os.path.join(OUT_DIR, "agent_focus_nearfield_weights.csv")
    write_focus_weights(weights_csv, points, rows)
    metrics_csv = os.path.join(OUT_DIR, "agent_focus_nearfield_metrics.csv")
    write_focus_metrics(metrics_csv, rows)
    patterns_png = os.path.join(OUT_DIR, "agent_focus_nearfield_patterns.png")
    plot_focus_patterns(rows, patterns_png)
    summary = os.path.join(OUT_DIR, "agent_focus_nearfield_summary.txt")
    write_focus_summary(summary, rows)
    for path in [summary, weights_csv, metrics_csv, patterns_png]:
        print(path)


if __name__ == "__main__":
    main()
