import csv
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(__file__))
import design_n36_array as base


ROOT = os.path.dirname(os.path.dirname(__file__))
OUT_DIR = os.path.join(ROOT, "analysis_outputs")
N38_CSV = os.path.join(OUT_DIR, "n38_refined_layout.csv")

THETA_MAX_DEG = 55.0
MAIN_EXCLUDE_DEG = 8.0
FOCUS_Z = 0.200
OFFSET_X = 0.040
DB_FLOOR = -60.0


def load_xy_csv(path):
    pts = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pts.append([float(row["X_mm"]) * 1e-3, float(row["Y_mm"]) * 1e-3])
    return np.asarray(pts, dtype=float)


def ensure_out_dir():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)


def farfield_vector(positions, theta, phi):
    theta = np.asarray(theta)
    phi = np.asarray(phi)
    kx = base.K * np.sin(theta) * np.cos(phi)
    ky = base.K * np.sin(theta) * np.sin(phi)
    phase = np.outer(kx.ravel(), positions[:, 0]) + np.outer(ky.ravel(), positions[:, 1])
    return base.element_pattern(theta).ravel()[:, None] * np.exp(-1j * phase)


def nearfield_vector(positions, x, y, z):
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()
    z = np.asarray(z).ravel()
    out = np.zeros((x.size, len(positions)), dtype=complex)
    for i, (x0, y0) in enumerate(positions):
        dx = x - x0
        dy = y - y0
        rr = np.sqrt(dx * dx + dy * dy + z * z)
        theta = np.arctan2(np.sqrt(dx * dx + dy * dy), z)
        out[:, i] = base.element_pattern(theta) * np.exp(1j * base.K * rr) / np.maximum(rr, 1e-9)
    return out


def target_vector(positions, case):
    if case["domain"] == "far":
        return farfield_vector(positions, np.deg2rad([case["theta_deg"]]), np.deg2rad([case["phi_deg"]]))[0]
    return nearfield_vector(
        positions,
        np.asarray([case["x_m"]]),
        np.asarray([case["y_m"]]),
        np.asarray([case["z_m"]]),
    )[0]


def phase_only_weights(positions, case):
    vt = target_vector(positions, case)
    return np.exp(-1j * np.angle(vt))


def radial_hamming_weights(positions, case):
    phase = phase_only_weights(positions, case)
    r = np.sqrt(np.sum(positions * positions, axis=1))
    rmax = np.max(r)
    amp = 0.54 + 0.46 * np.cos(np.pi * r / max(rmax, 1e-12))
    return amp * phase


def radial_tukey_weights(positions, case):
    phase = phase_only_weights(positions, case)
    r = np.sqrt(np.sum(positions * positions, axis=1))
    rmax = np.max(r)
    beta = 0.58
    edge = 0.22
    amp = np.ones_like(r)
    mask = r > beta * rmax
    t = (r[mask] - beta * rmax) / max((1.0 - beta) * rmax, 1e-12)
    amp[mask] = edge + (1.0 - edge) * 0.5 * (1.0 + np.cos(np.pi * t))
    return amp * phase


def farfield_sidelobe_matrix(positions, case, theta_step=1.0, phi_step=3.0):
    theta_deg = np.arange(0.0, THETA_MAX_DEG + 0.001, theta_step)
    phi_deg = np.arange(0.0, 360.0, phi_step)
    theta_grid_deg, phi_grid_deg = np.meshgrid(theta_deg, phi_deg)
    theta = np.deg2rad(theta_grid_deg)
    phi = np.deg2rad(phi_grid_deg)
    h = farfield_vector(positions, theta, phi)

    target = base.direction_vector(case["theta_deg"], case["phi_deg"])
    ux = np.sin(theta) * np.cos(phi)
    uy = np.sin(theta) * np.sin(phi)
    uz = np.cos(theta)
    sep = np.rad2deg(np.arccos(np.clip(ux * target[0] + uy * target[1] + uz * target[2], -1.0, 1.0)))
    mask = sep.ravel() > MAIN_EXCLUDE_DEG
    return h[mask]


def nearfield_sidelobe_matrix(positions, case, step_mm=3.0):
    x_mm = np.arange(-120.0, 120.0 + 0.001, step_mm)
    y_mm = np.arange(-120.0, 120.0 + 0.001, step_mm)
    xg, yg = np.meshgrid(x_mm * 1e-3, y_mm * 1e-3)
    zg = np.zeros_like(xg) + case["z_m"]
    h = nearfield_vector(positions, xg, yg, zg)
    rr_mm = np.sqrt((x_mm[None, :] - case["x_m"] * 1000.0) ** 2 + (y_mm[:, None] - case["y_m"] * 1000.0) ** 2)
    return h[(rr_mm.ravel() > 20.0)]


def sidelobe_matrix(positions, case):
    if case["domain"] == "far":
        return farfield_sidelobe_matrix(positions, case)
    return nearfield_sidelobe_matrix(positions, case)


def rls_amplitude_weights(positions, case, target_fraction=0.72, mu=1e-4):
    phase = phase_only_weights(positions, case)
    hs = sidelobe_matrix(positions, case)
    b = hs * phase[None, :]
    vt = target_vector(positions, case)
    g = np.real(vt * phase)
    g = np.maximum(g, 0.0)
    g0 = float(np.sum(g))
    c = np.real(np.dot(b.conj().T, b)) / max(len(b), 1)
    c += mu * np.eye(len(positions))

    def obj(a):
        return float(np.dot(a, np.dot(c, a)))

    def jac(a):
        return 2.0 * np.dot(c, a)

    cons = [{"type": "ineq", "fun": lambda a: float(np.dot(g, a) - target_fraction * g0),
             "jac": lambda a: g}]
    x0 = np.clip(radial_tukey_weights(positions, case) / phase, 0.0, 1.0).real
    res = minimize(obj, x0, jac=jac, bounds=[(0.0, 1.0)] * len(positions), constraints=cons, method="SLSQP", options={"maxiter": 160, "ftol": 1e-10, "disp": False})
    amp = np.clip(res.x if res.success else x0, 0.0, 1.0)
    return amp * phase


def mvdr_complex_weights(positions, case, diagonal_load=1e-3):
    hs = sidelobe_matrix(positions, case)
    vt = target_vector(positions, case)
    r = np.dot(hs.conj().T, hs) / max(len(hs), 1)
    tr = np.real(np.trace(r)) / len(positions)
    r += diagonal_load * max(tr, 1e-12) * np.eye(len(positions))
    rinv_v = np.linalg.solve(r, vt.conj())
    denom = np.dot(vt, rinv_v)
    w = rinv_v / denom
    m = np.max(np.abs(w))
    return w / max(m, 1e-12)


def make_weights(positions, case, method):
    if method == "phase_only":
        return phase_only_weights(positions, case)
    if method == "radial_tukey":
        return radial_tukey_weights(positions, case)
    if method == "radial_hamming":
        return radial_hamming_weights(positions, case)
    if method == "rls_amp":
        return rls_amplitude_weights(positions, case)
    if method == "mvdr_complex":
        return mvdr_complex_weights(positions, case)
    raise ValueError(method)


def response_db(response, target_amp):
    amp = np.abs(response)
    return 20.0 * np.log10(amp / max(target_amp, 1e-18) + 1e-12)


def farfield_eval(positions, case, weights):
    theta_deg = np.arange(0.0, THETA_MAX_DEG + 0.001, 0.5)
    phi_deg = np.arange(0.0, 360.0, 1.5)
    theta_grid_deg, phi_grid_deg = np.meshgrid(theta_deg, phi_deg)
    theta = np.deg2rad(theta_grid_deg)
    phi = np.deg2rad(phi_grid_deg)
    h = farfield_vector(positions, theta, phi)
    resp = np.dot(h, weights)
    vt = target_vector(positions, case)
    target_amp = float(abs(np.dot(vt, weights)))

    target = base.direction_vector(case["theta_deg"], case["phi_deg"])
    ux = np.sin(theta) * np.cos(phi)
    uy = np.sin(theta) * np.sin(phi)
    uz = np.cos(theta)
    sep = np.rad2deg(np.arccos(np.clip(ux * target[0] + uy * target[1] + uz * target[2], -1.0, 1.0)))
    outside = sep.ravel() > MAIN_EXCLUDE_DEG
    db = response_db(resp, target_amp)
    psl = float(np.max(db[outside]))
    p995 = float(np.percentile(db[outside], 99.5))

    alpha = np.linspace(-55.0, 55.0, 1101)
    hcut = farfield_vector(positions, np.abs(np.deg2rad(alpha)), np.where(alpha >= 0.0, 0.0, np.pi))
    cut_db = response_db(np.dot(hcut, weights), target_amp)
    peak_idx = int(np.argmin(np.abs(alpha - case["theta_deg"])))
    width = base.fwhm(alpha, np.maximum(cut_db, DB_FLOOR), peak_idx)
    return psl, p995, width, target_amp, alpha, np.maximum(cut_db, DB_FLOOR)


def nearfield_eval(positions, case, weights):
    x_mm = np.linspace(-120.0, 120.0, 241)
    y_mm = np.linspace(-120.0, 120.0, 241)
    xg, yg = np.meshgrid(x_mm * 1e-3, y_mm * 1e-3)
    zg = np.zeros_like(xg) + case["z_m"]
    h = nearfield_vector(positions, xg, yg, zg)
    resp = np.dot(h, weights)
    vt = target_vector(positions, case)
    target_amp = float(abs(np.dot(vt, weights)))
    db = response_db(resp, target_amp).reshape(xg.shape)
    rr = np.sqrt((x_mm[None, :] - case["x_m"] * 1000.0) ** 2 + (y_mm[:, None] - case["y_m"] * 1000.0) ** 2)
    outside = rr > 20.0
    psl = float(np.max(db[outside]))
    p995 = float(np.percentile(db[outside], 99.5))
    idx = np.unravel_index(np.argmax(np.abs(resp)), xg.shape)
    peak_x = float(x_mm[idx[1]])
    peak_y = float(y_mm[idx[0]])
    err = math.hypot(peak_x - case["x_m"] * 1000.0, peak_y - case["y_m"] * 1000.0)
    wx = base.fwhm(x_mm, np.maximum(db[idx[0], :], DB_FLOOR), idx[1])
    wy = base.fwhm(y_mm, np.maximum(db[:, idx[1]], DB_FLOOR), idx[0])
    return psl, p995, 0.5 * (wx + wy), target_amp, err, x_mm, y_mm, np.maximum(db, DB_FLOOR)


def summarize_weights(weights):
    amp = np.abs(weights)
    active_50 = int(np.sum(amp >= 0.5 * np.max(amp)))
    return float(np.min(amp)), float(np.mean(amp)), float(np.max(amp)), active_50


def run():
    ensure_out_dir()
    layouts = [
        ("N32_ring", base.load_current_layout()),
        ("N38_blue_noise", load_xy_csv(N38_CSV)),
    ]
    cases = [
        {"name": "broadside", "domain": "far", "theta_deg": 0.0, "phi_deg": 0.0},
        {"name": "steer20", "domain": "far", "theta_deg": 20.0, "phi_deg": 0.0},
        {"name": "focus200", "domain": "near", "x_m": 0.0, "y_m": 0.0, "z_m": FOCUS_Z},
        {"name": "offset40_200", "domain": "near", "x_m": OFFSET_X, "y_m": 0.0, "z_m": FOCUS_Z},
    ]
    methods = ["phase_only", "radial_tukey", "radial_hamming", "rls_amp", "mvdr_complex"]

    rows = []
    cut_data = {}
    heat_data = {}
    for layout_name, positions in layouts:
        for case in cases:
            phase_w = make_weights(positions, case, "phase_only")
            phase_target = abs(np.dot(target_vector(positions, case), phase_w))
            for method in methods:
                weights = make_weights(positions, case, method)
                amin, amean, amax, active_50 = summarize_weights(weights)
                if case["domain"] == "far":
                    psl, p995, width, target_amp, alpha, cut_db = farfield_eval(positions, case, weights)
                    err = 0.0
                    if case["name"] == "steer20":
                        cut_data[(layout_name, method)] = (alpha, cut_db)
                else:
                    psl, p995, width, target_amp, err, x_mm, y_mm, db = nearfield_eval(positions, case, weights)
                    if case["name"] == "offset40_200":
                        heat_data[(layout_name, method)] = (x_mm, y_mm, db)
                gain_drop = 20.0 * math.log10(target_amp / max(phase_target, 1e-18) + 1e-12)
                rows.append({
                    "layout": layout_name,
                    "case": case["name"],
                    "method": method,
                    "psl_db": psl,
                    "p99_5_db": p995,
                    "width": width,
                    "target_gain_drop_db": gain_drop,
                    "focus_error_mm": err,
                    "amp_min": amin,
                    "amp_mean": amean,
                    "amp_max": amax,
                    "active_ge_50pct": active_50,
                })

    csv_path = os.path.join(OUT_DIR, "software_weighting_summary.csv")
    with open(csv_path, "w", newline="") as f:
        fieldnames = [
            "layout", "case", "method", "psl_db", "p99_5_db", "width",
            "target_gain_drop_db", "focus_error_mm", "amp_min", "amp_mean",
            "amp_max", "active_ge_50pct",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    txt_path = os.path.join(OUT_DIR, "software_weighting_summary.txt")
    with open(txt_path, "w") as f:
        f.write("Software weighting comparison, 16 mm circular piston, 40 kHz\n")
        f.write("Methods: phase_only, radial_tukey, radial_hamming, rls_amp, mvdr_complex\n")
        f.write("PSL is relative to the target point/direction. Gain drop is relative to phase_only with max channel amplitude=1.\n\n")
        for layout_name, _ in layouts:
            f.write("[%s]\n" % layout_name)
            for case in cases:
                f.write("  %s\n" % case["name"])
                subset = [r for r in rows if r["layout"] == layout_name and r["case"] == case["name"]]
                for r in subset:
                    unit = "deg" if case["domain"] == "far" else "mm"
                    f.write(
                        "    %-14s PSL=%6.2f dB, p99.5=%6.2f dB, width=%6.2f %s, gain=%6.2f dB, err=%5.2f mm, amp_mean=%.2f, active=%d\n"
                        % (
                            r["method"], r["psl_db"], r["p99_5_db"], r["width"], unit,
                            r["target_gain_drop_db"], r["focus_error_mm"], r["amp_mean"], r["active_ge_50pct"],
                        )
                    )
            f.write("\n")

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), sharey=True)
    for ax, layout_name in zip(axes, ["N32_ring", "N38_blue_noise"]):
        for method in methods:
            alpha, cut_db = cut_data[(layout_name, method)]
            ax.plot(alpha, cut_db, lw=1.4, label=method)
        ax.axvline(20.0, color="k", lw=0.8, alpha=0.35)
        ax.set_title("%s, steer 20 deg" % layout_name)
        ax.set_xlabel("signed angle (deg)")
        ax.set_xlim(-55, 55)
        ax.set_ylim(-45, 4)
        ax.grid(True, alpha=0.25)
    axes[0].set_ylabel("level relative to target (dB)")
    axes[1].legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "software_weighting_steer20_cuts.png"), dpi=180)
    plt.close(fig)

    selected_methods = ["phase_only", "radial_tukey", "rls_amp", "mvdr_complex"]
    fig, axes = plt.subplots(2, 4, figsize=(15.0, 7.8))
    levels = np.linspace(-45.0, 3.0, 49)
    contour = None
    for row_idx, layout_name in enumerate(["N32_ring", "N38_blue_noise"]):
        for col_idx, method in enumerate(selected_methods):
            x_mm, y_mm, db = heat_data[(layout_name, method)]
            ax = axes[row_idx, col_idx]
            contour = ax.contourf(x_mm, y_mm, db, levels=levels, cmap="magma", extend="both")
            ax.plot([40.0], [0.0], "cx", ms=7, mew=1.5)
            ax.set_aspect("equal", adjustable="box")
            ax.set_title("%s\n%s" % (layout_name, method), fontsize=9)
            ax.set_xlabel("x (mm)")
            if col_idx == 0:
                ax.set_ylabel("y (mm)")
            ax.grid(True, color="white", alpha=0.10, lw=0.5)
    cbar = fig.colorbar(contour, ax=axes.ravel().tolist(), shrink=0.86, pad=0.02)
    cbar.set_label("level relative to target (dB)")
    fig.suptitle("Offset focus x=40 mm, z=200 mm: software weighting")
    fig.savefig(os.path.join(OUT_DIR, "software_weighting_offset_focus_heatmaps.png"), dpi=180)
    plt.close(fig)

    print(csv_path)
    print(txt_path)
    print(os.path.join(OUT_DIR, "software_weighting_steer20_cuts.png"))
    print(os.path.join(OUT_DIR, "software_weighting_offset_focus_heatmaps.png"))


if __name__ == "__main__":
    run()
