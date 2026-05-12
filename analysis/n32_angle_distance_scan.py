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
import software_weighting_compare as sw


ROOT = os.path.dirname(os.path.dirname(__file__))
OUT_DIR = os.path.join(ROOT, "analysis_outputs")


def ensure_out_dir():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)


def target_case(x_m, z_m):
    return {"name": "target", "domain": "near", "x_m": x_m, "y_m": 0.0, "z_m": z_m}


def local_nearfield_matrix(positions, target, span_mm, step_mm):
    x0 = target["x_m"] * 1000.0
    y0 = target["y_m"] * 1000.0
    x_mm = np.arange(x0 - span_mm, x0 + span_mm + 0.001, step_mm)
    y_mm = np.arange(y0 - span_mm, y0 + span_mm + 0.001, step_mm)
    xg, yg = np.meshgrid(x_mm * 1e-3, y_mm * 1e-3)
    zg = np.zeros_like(xg) + target["z_m"]
    h = sw.nearfield_vector(positions, xg, yg, zg)
    rr_mm = np.sqrt((x_mm[None, :] - x0) ** 2 + (y_mm[:, None] - y0) ** 2)
    return x_mm, y_mm, h, rr_mm


def local_rls_amp_weights(positions, target, span_mm=140.0, step_mm=4.0, target_fraction=0.72):
    phase = sw.phase_only_weights(positions, target)
    _, _, h, rr = local_nearfield_matrix(positions, target, span_mm, step_mm)
    hs = h[rr.ravel() > 20.0]
    b = hs * phase[None, :]
    vt = sw.target_vector(positions, target)
    g = np.real(vt * phase)
    g = np.maximum(g, 0.0)
    g0 = float(np.sum(g))
    c = np.real(np.dot(b.conj().T, b)) / max(len(b), 1)
    c += 1e-4 * np.eye(len(positions))

    def obj(a):
        return float(np.dot(a, np.dot(c, a)))

    def jac(a):
        return 2.0 * np.dot(c, a)

    cons = [{"type": "ineq", "fun": lambda a: float(np.dot(g, a) - target_fraction * g0), "jac": lambda a: g}]
    x0 = np.ones(len(positions)) * target_fraction
    res = minimize(obj, x0, jac=jac, bounds=[(0.0, 1.0)] * len(positions), constraints=cons, method="SLSQP", options={"maxiter": 120, "ftol": 1e-9})
    amp = np.clip(res.x if res.success else x0, 0.0, 1.0)
    return amp * phase


def line_width(axis_mm, line_db, idx, threshold_db):
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


def eval_target(positions, target, method, reference_amp, span_mm=160.0, step_mm=2.5):
    phase_w = sw.phase_only_weights(positions, target)
    if method == "phase_only":
        weights = phase_w
    elif method == "rls_amp":
        weights = local_rls_amp_weights(positions, target, span_mm=span_mm, step_mm=max(4.0, step_mm * 1.6))
    else:
        raise ValueError(method)

    x_mm, y_mm, h, rr = local_nearfield_matrix(positions, target, span_mm, step_mm)
    resp = np.dot(h, weights).reshape(rr.shape)
    vt = sw.target_vector(positions, target)
    target_amp = float(abs(np.dot(vt, weights)))
    db = 20.0 * np.log10(np.abs(resp) / max(target_amp, 1e-18) + 1e-12)
    outside = rr > 20.0
    psl = float(np.max(db[outside]))
    p995 = float(np.percentile(db[outside], 99.5))
    idx = np.unravel_index(np.argmax(np.abs(resp)), resp.shape)
    tx = target["x_m"] * 1000.0
    ty = target["y_m"] * 1000.0
    err = math.hypot(float(x_mm[idx[1]]) - tx, float(y_mm[idx[0]]) - ty)
    ix = int(np.argmin(np.abs(x_mm - tx)))
    iy = int(np.argmin(np.abs(y_mm - ty)))
    wx = line_width(x_mm, db[iy, :], ix, -3.0)
    wy = line_width(y_mm, db[:, ix], iy, -3.0)
    gain_ref_db = 20.0 * math.log10(target_amp / max(reference_amp, 1e-18) + 1e-12)
    return psl, p995, err, wx, wy, gain_ref_db, target_amp


def run_angle_scan(positions):
    ref_target = target_case(0.0, 0.200)
    ref_w = sw.phase_only_weights(positions, ref_target)
    ref_amp = float(abs(np.dot(sw.target_vector(positions, ref_target), ref_w)))
    rows = []
    for angle in range(-85, 86, 5):
        if abs(angle) >= 90:
            continue
        z_m = 0.200
        x_m = z_m * math.tan(math.radians(angle))
        target = target_case(x_m, z_m)
        elem_db = 20.0 * math.log10(float(base.element_pattern(np.asarray([abs(math.radians(angle))]))[0]) + 1e-12)
        for method in ["phase_only", "rls_amp"]:
            psl, p995, err, wx, wy, gain_ref_db, target_amp = eval_target(positions, target, method, ref_amp)
            rows.append({
                "angle_deg": angle,
                "x_mm": x_m * 1000.0,
                "z_mm": z_m * 1000.0,
                "method": method,
                "element_directivity_db": elem_db,
                "pressure_psl_db": psl,
                "audio_proxy_psl_db": 2.0 * psl,
                "pressure_p99_5_db": p995,
                "focus_error_mm": err,
                "fwhm_x_mm": wx,
                "fwhm_y_mm": wy,
                "target_gain_vs_onaxis200_db": gain_ref_db,
                "target_amp": target_amp,
            })
    return rows


def run_distance_scan(positions):
    ref_target = target_case(0.0, 0.200)
    ref_w = sw.phase_only_weights(positions, ref_target)
    ref_amp = float(abs(np.dot(sw.target_vector(positions, ref_target), ref_w)))
    distances = [60, 80, 100, 120, 150, 200, 250, 300, 350, 400, 500, 700, 1000, 1500, 2000]
    angles = [0, 10, 15]
    rows = []
    for z_mm in distances:
        span = max(140.0, 0.45 * z_mm)
        step = 3.0 if z_mm <= 500 else 6.0
        for angle in angles:
            z_m = z_mm * 1e-3
            x_m = z_m * math.tan(math.radians(angle))
            target = target_case(x_m, z_m)
            for method in ["phase_only", "rls_amp"]:
                psl, p995, err, wx, wy, gain_ref_db, target_amp = eval_target(positions, target, method, ref_amp, span_mm=span, step_mm=step)
                rows.append({
                    "z_mm": z_mm,
                    "angle_deg": angle,
                    "x_mm": x_m * 1000.0,
                    "method": method,
                    "pressure_psl_db": psl,
                    "audio_proxy_psl_db": 2.0 * psl,
                    "pressure_p99_5_db": p995,
                    "focus_error_mm": err,
                    "fwhm_x_mm": wx,
                    "fwhm_y_mm": wy,
                    "target_gain_vs_onaxis200_db": gain_ref_db,
                    "target_amp": target_amp,
                })
    return rows


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_angle(rows):
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 7.2), sharex=True)
    for method in ["phase_only", "rls_amp"]:
        subset = [r for r in rows if r["method"] == method]
        axes[0].plot([r["angle_deg"] for r in subset], [r["pressure_psl_db"] for r in subset], marker="o", lw=1.5, label=method)
        axes[1].plot([r["angle_deg"] for r in subset], [r["target_gain_vs_onaxis200_db"] for r in subset], marker="o", lw=1.5, label=method)
    elem = [r for r in rows if r["method"] == "phase_only"]
    axes[1].plot([r["angle_deg"] for r in elem], [r["element_directivity_db"] for r in elem], color="k", ls="--", lw=1.2, label="single element")
    axes[0].axhline(-10.0, color="k", ls=":", lw=1.0)
    axes[0].set_ylabel("pressure PSL (dB)")
    axes[0].set_title("N32 focus at z=200 mm vs steering angle")
    axes[0].grid(True, alpha=0.25)
    axes[1].set_ylabel("target gain vs on-axis 200 mm (dB)")
    axes[1].set_xlabel("steering/focus angle (deg)")
    axes[1].grid(True, alpha=0.25)
    axes[0].legend()
    axes[1].legend()
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "n32_angle_scan.png")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_distance(rows):
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 7.2), sharex=True)
    for angle in [0, 10, 15]:
        for method, ls in [("phase_only", "-"), ("rls_amp", "--")]:
            subset = [r for r in rows if r["angle_deg"] == angle and r["method"] == method]
            label = "%d deg %s" % (angle, method)
            axes[0].plot([r["z_mm"] for r in subset], [r["pressure_psl_db"] for r in subset], marker="o", lw=1.4, ls=ls, label=label)
            axes[1].plot([r["z_mm"] for r in subset], [0.5 * (r["fwhm_x_mm"] + r["fwhm_y_mm"]) for r in subset], marker="o", lw=1.4, ls=ls, label=label)
    axes[0].axhline(-10.0, color="k", ls=":", lw=1.0)
    axes[0].set_ylabel("pressure PSL (dB)")
    axes[0].set_title("N32 focus performance vs distance")
    axes[0].grid(True, alpha=0.25)
    axes[1].set_ylabel("mean pressure FWHM (mm)")
    axes[1].set_xlabel("focus z distance (mm)")
    axes[1].grid(True, alpha=0.25)
    axes[0].legend(ncol=2, fontsize=8)
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "n32_distance_scan.png")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def write_summary(angle_rows, distance_rows, path):
    with open(path, "w") as f:
        f.write("N32 angle and distance scan\n")
        f.write("Angle scan: targets are on z=200 mm plane, x=z*tan(angle).\n")
        f.write("Distance scan: targets use angles 0/10/15 deg.\n")
        f.write("PSL uses a local z-plane map and excludes a 20 mm radius around the target.\n\n")
        f.write("Angle scan representative phase-only rows:\n")
        for angle in [0, 10, 15, 20, 25, 30, 40, 50, 60]:
            row = next((r for r in angle_rows if r["angle_deg"] == angle and r["method"] == "phase_only"), None)
            if row:
                f.write(
                    "  angle=%3.0f deg x=%7.1f mm: PSL=%6.2f dB, gain=%7.2f dB, elem=%7.2f dB, FWHM=(%5.1f,%5.1f) mm\n"
                    % (row["angle_deg"], row["x_mm"], row["pressure_psl_db"], row["target_gain_vs_onaxis200_db"], row["element_directivity_db"], row["fwhm_x_mm"], row["fwhm_y_mm"])
                )
        f.write("\nAngle scan representative rls_amp rows:\n")
        for angle in [0, 10, 15, 20, 25, 30, 40, 50, 60]:
            row = next((r for r in angle_rows if r["angle_deg"] == angle and r["method"] == "rls_amp"), None)
            if row:
                f.write(
                    "  angle=%3.0f deg x=%7.1f mm: PSL=%6.2f dB, gain=%7.2f dB, FWHM=(%5.1f,%5.1f) mm\n"
                    % (row["angle_deg"], row["x_mm"], row["pressure_psl_db"], row["target_gain_vs_onaxis200_db"], row["fwhm_x_mm"], row["fwhm_y_mm"])
                )
        f.write("\nDistance scan, axial phase-only:\n")
        for z in [100, 150, 200, 300, 400, 500, 700, 1000, 1500, 2000]:
            row = next((r for r in distance_rows if r["z_mm"] == z and r["angle_deg"] == 0 and r["method"] == "phase_only"), None)
            if row:
                f.write(
                    "  z=%4.0f mm: PSL=%6.2f dB, gain=%7.2f dB, FWHM=(%6.1f,%6.1f) mm\n"
                    % (row["z_mm"], row["pressure_psl_db"], row["target_gain_vs_onaxis200_db"], row["fwhm_x_mm"], row["fwhm_y_mm"])
                )
        f.write("\nDistance scan, 10 deg rls_amp:\n")
        for z in [100, 150, 200, 300, 400, 500, 700, 1000]:
            row = next((r for r in distance_rows if r["z_mm"] == z and r["angle_deg"] == 10 and r["method"] == "rls_amp"), None)
            if row:
                f.write(
                    "  z=%4.0f mm x=%6.1f mm: PSL=%6.2f dB, gain=%7.2f dB, FWHM=(%6.1f,%6.1f) mm\n"
                    % (row["z_mm"], row["x_mm"], row["pressure_psl_db"], row["target_gain_vs_onaxis200_db"], row["fwhm_x_mm"], row["fwhm_y_mm"])
                )


def main():
    ensure_out_dir()
    positions = base.load_current_layout()
    angle_rows = run_angle_scan(positions)
    distance_rows = run_distance_scan(positions)
    angle_csv = os.path.join(OUT_DIR, "n32_angle_scan.csv")
    distance_csv = os.path.join(OUT_DIR, "n32_distance_scan.csv")
    write_csv(angle_csv, angle_rows)
    write_csv(distance_csv, distance_rows)
    angle_png = plot_angle(angle_rows)
    distance_png = plot_distance(distance_rows)
    summary = os.path.join(OUT_DIR, "n32_angle_distance_summary.txt")
    write_summary(angle_rows, distance_rows, summary)
    for path in [angle_csv, distance_csv, angle_png, distance_png, summary]:
        print(path)


if __name__ == "__main__":
    main()
