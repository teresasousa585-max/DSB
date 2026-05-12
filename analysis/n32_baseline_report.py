import csv
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import audio_focus_performance as audio
import design_n36_array as base
import software_weighting_compare as sw


ROOT = os.path.dirname(os.path.dirname(__file__))
OUT_DIR = os.path.join(ROOT, "analysis_outputs")

THETA_MAX_DEG = 55.0
MAIN_EXCLUDE_DEG = 8.0
FOCUS_Z_M = 0.200
OFFSET_X_M = 0.040
DB_FLOOR = -45.0


def ensure_out_dir():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)


def n32_positions():
    return base.load_current_layout()


def normalize_db_in_mask(field, mask):
    amp = np.abs(field)
    peak = np.nanmax(amp[mask])
    db = 20.0 * np.log10(amp / max(peak, 1e-18) + 1e-12)
    db = np.where(mask, np.maximum(db, DB_FLOOR), np.nan)
    return db


def circular_cases():
    return [
        {
            "name": "broadside",
            "title": "Broadside, phase-only",
            "domain": "far",
            "theta_deg": 0.0,
            "phi_deg": 0.0,
            "target_theta_deg": 0.0,
            "target_phi_deg": 0.0,
        },
        {
            "name": "steer20",
            "title": "Steer 20 deg, phase-only",
            "domain": "far",
            "theta_deg": 20.0,
            "phi_deg": 0.0,
            "target_theta_deg": 20.0,
            "target_phi_deg": 0.0,
        },
        {
            "name": "focus200",
            "title": "Focus z=200 mm",
            "domain": "near",
            "x_m": 0.0,
            "y_m": 0.0,
            "z_m": FOCUS_Z_M,
            "target_theta_deg": 0.0,
            "target_phi_deg": 0.0,
        },
        {
            "name": "offset40_200",
            "title": "Focus x=40, z=200 mm",
            "domain": "near",
            "x_m": OFFSET_X_M,
            "y_m": 0.0,
            "z_m": FOCUS_Z_M,
            "target_theta_deg": math.degrees(math.atan2(OFFSET_X_M, FOCUS_Z_M)),
            "target_phi_deg": 0.0,
        },
    ]


def circular_plot_and_metrics(positions):
    axis_deg = np.linspace(-THETA_MAX_DEG, THETA_MAX_DEG, 241)
    tx, ty = np.meshgrid(axis_deg, axis_deg)
    theta_deg = np.sqrt(tx * tx + ty * ty)
    phi_deg = (np.rad2deg(np.arctan2(ty, tx)) + 360.0) % 360.0
    mask = theta_deg <= THETA_MAX_DEG
    theta = np.deg2rad(np.minimum(theta_deg, THETA_MAX_DEG))
    phi = np.deg2rad(phi_deg)
    h = sw.farfield_vector(positions, theta, phi)

    ux = np.sin(theta) * np.cos(phi)
    uy = np.sin(theta) * np.sin(phi)
    uz = np.cos(theta)

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 10.2))
    levels = np.linspace(DB_FLOOR, 0.0, 46)
    rows = []
    contour = None
    for ax, case in zip(axes.ravel(), circular_cases()):
        weights = sw.phase_only_weights(positions, case)
        response = np.dot(h, weights).reshape(theta.shape)
        db = normalize_db_in_mask(response, mask)
        contour = ax.contourf(tx, ty, db, levels=levels, cmap="viridis", extend="min")
        target_theta = case["target_theta_deg"]
        target_phi = case["target_phi_deg"]
        mx = target_theta * math.cos(math.radians(target_phi))
        my = target_theta * math.sin(math.radians(target_phi))
        ax.plot([mx], [my], "rx", ms=8, mew=1.8)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(-THETA_MAX_DEG, THETA_MAX_DEG)
        ax.set_ylim(-THETA_MAX_DEG, THETA_MAX_DEG)
        ax.set_title(case["title"])
        ax.set_xlabel("theta_x (deg)")
        ax.set_ylabel("theta_y (deg)")
        ax.grid(True, color="white", alpha=0.12, lw=0.5)

        target_vec = base.direction_vector(target_theta, target_phi)
        sep = np.rad2deg(np.arccos(np.clip(ux * target_vec[0] + uy * target_vec[1] + uz * target_vec[2], -1.0, 1.0)))
        outside = (sep > MAIN_EXCLUDE_DEG) & mask
        peak_idx = np.nanargmax(db)
        peak_2d = np.unravel_index(peak_idx, db.shape)
        peak_theta = float(theta_deg[peak_2d])
        peak_phi = float(phi_deg[peak_2d])
        peak_vec = base.direction_vector(peak_theta, peak_phi)
        pointing_error = float(np.rad2deg(np.arccos(np.clip(np.dot(peak_vec, target_vec), -1.0, 1.0))))
        rows.append({
            "case": case["name"],
            "target_theta_deg": target_theta,
            "peak_theta_deg": peak_theta,
            "peak_phi_deg": peak_phi,
            "pointing_error_deg": pointing_error,
            "max_outside_8deg_db": float(np.nanmax(db[outside])),
        })
    cbar = fig.colorbar(contour, ax=axes.ravel().tolist(), shrink=0.86, pad=0.03)
    cbar.set_label("normalized level (dB)")
    fig.suptitle("N32 circular normalized 40 kHz carrier patterns")
    path = os.path.join(OUT_DIR, "n32_baseline_circular_patterns.png")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return rows, path


def focus_heatmaps(positions):
    x_mm = np.linspace(-120.0, 120.0, 241)
    y_mm = np.linspace(-120.0, 120.0, 241)
    xg, yg = np.meshgrid(x_mm * 1e-3, y_mm * 1e-3)
    zg = np.zeros_like(xg) + FOCUS_Z_M

    cases = [
        ("focus200", "Focus z=200 mm", {"name": "focus200", "domain": "near", "x_m": 0.0, "y_m": 0.0, "z_m": FOCUS_Z_M}, 0.0, 0.0),
        ("offset40_200", "Focus x=40, z=200 mm", {"name": "offset40_200", "domain": "near", "x_m": OFFSET_X_M, "y_m": 0.0, "z_m": FOCUS_Z_M}, 40.0, 0.0),
    ]
    methods = [("phase_only", sw.phase_only_weights), ("rls_amp", audio.rls_amp_weights)]
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 10.2))
    levels = np.linspace(DB_FLOOR, 3.0, 49)
    contour = None
    rows = []
    for row_idx, (case_name, case_title, case, marker_x, marker_y) in enumerate(cases):
        phase_w = sw.phase_only_weights(positions, case)
        phase_target_amp = abs(np.dot(sw.target_vector(positions, case), phase_w))
        for col_idx, (method_name, weight_fn) in enumerate(methods):
            weights = weight_fn(positions, case)
            h = sw.nearfield_vector(positions, xg, yg, zg)
            response = np.dot(h, weights).reshape(xg.shape)
            target_amp = abs(np.dot(sw.target_vector(positions, case), weights))
            db = 20.0 * np.log10(np.abs(response) / max(target_amp, 1e-18) + 1e-12)
            db = np.maximum(db, DB_FLOOR)
            ax = axes[row_idx, col_idx]
            contour = ax.contourf(x_mm, y_mm, db, levels=levels, cmap="magma", extend="both")
            ax.plot([marker_x], [marker_y], "cx", ms=8, mew=1.8)
            ax.set_aspect("equal", adjustable="box")
            ax.set_title("%s\n%s" % (case_title, method_name))
            ax.set_xlabel("x (mm)")
            ax.set_ylabel("y (mm)")
            ax.grid(True, color="white", alpha=0.12, lw=0.5)

            rr = np.sqrt((x_mm[None, :] - marker_x) ** 2 + (y_mm[:, None] - marker_y) ** 2)
            outside = rr > 20.0
            idx = np.unravel_index(np.argmax(np.abs(response)), response.shape)
            peak_x = float(x_mm[idx[1]])
            peak_y = float(y_mm[idx[0]])
            target_ix = int(np.argmin(np.abs(x_mm - marker_x)))
            target_iy = int(np.argmin(np.abs(y_mm - marker_y)))
            rows.append({
                "case": case_name,
                "method": method_name,
                "pressure_psl_db": float(np.max(db[outside])),
                "audio_source_psl_proxy_db": float(2.0 * np.max(db[outside])),
                "pressure_p99_5_db": float(np.percentile(db[outside], 99.5)),
                "focus_error_mm": float(math.hypot(peak_x - marker_x, peak_y - marker_y)),
                "pressure_fwhm_x_mm": audio.line_width(x_mm, db[target_iy, :], target_ix, -3.0),
                "pressure_fwhm_y_mm": audio.line_width(y_mm, db[:, target_ix], target_iy, -3.0),
                "audio_source_fwhm_x_mm": audio.line_width(x_mm, db[target_iy, :], target_ix, -1.5),
                "audio_source_fwhm_y_mm": audio.line_width(y_mm, db[:, target_ix], target_iy, -1.5),
                "target_gain_drop_db": float(20.0 * np.log10(target_amp / max(phase_target_amp, 1e-18) + 1e-12)),
            })
    fig.subplots_adjust(hspace=0.36, wspace=0.28, top=0.90, right=0.86)
    cbar_ax = fig.add_axes([0.89, 0.16, 0.025, 0.68])
    cbar = fig.colorbar(contour, cax=cbar_ax)
    cbar.set_label("level relative to target (dB)")
    fig.suptitle("N32 z=200 mm carrier focus maps", y=0.98)
    path = os.path.join(OUT_DIR, "n32_baseline_focus_heatmaps.png")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return rows, path


def scan_audio_focus(positions):
    rows = audio.scan_layout("N32_ring", positions)
    scan_path = os.path.join(OUT_DIR, "n32_audio_focus_scan.csv")
    audio.write_csv(scan_path, rows)
    usable = audio.usable_radius_by_z(rows)
    usable_path = os.path.join(OUT_DIR, "n32_audio_focus_usable_radius.csv")
    with open(usable_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["layout", "method", "z_mm", "usable_radius_mm", "usable_angle_deg"])
        writer.writeheader()
        for row in usable:
            writer.writerow(row)
    return rows, usable, scan_path, usable_path


def write_metrics_csv(path, circular_rows, focus_rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["section", "case", "method", "metric", "value"])
        for row in circular_rows:
            for key, value in row.items():
                if key != "case":
                    writer.writerow(["circular_farfield", row["case"], "phase_only", key, value])
        for row in focus_rows:
            for key, value in row.items():
                if key not in ["case", "method"]:
                    writer.writerow(["nearfield_z200", row["case"], row["method"], key, value])


def write_summary(path, positions, circular_rows, focus_rows, usable_rows):
    aperture_mm = 2.0 * base.aperture_radius(positions) * 1000.0
    min_dist_mm = base.min_distance(positions) * 1000.0
    with open(path, "w") as f:
        f.write("N32 baseline report, 40 kHz carrier, 16 mm circular piston elements\n")
        f.write("element_count: %d\n" % len(positions))
        f.write("wavelength_mm: %.3f\n" % (base.LAMBDA * 1000.0))
        f.write("element_diameter_mm: %.3f\n" % (base.ELEMENT_DIAMETER_M * 1000.0))
        f.write("aperture_diameter_mm: %.3f\n" % aperture_mm)
        f.write("aperture_in_lambda: %.3f\n" % (aperture_mm / (base.LAMBDA * 1000.0)))
        f.write("min_center_distance_mm: %.3f\n" % min_dist_mm)
        f.write("min_center_distance_in_lambda: %.3f\n" % (min_dist_mm / (base.LAMBDA * 1000.0)))
        f.write("\nCircular normalized carrier patterns:\n")
        for row in circular_rows:
            f.write(
                "  %-14s target=%5.2f deg, peak=%5.2f deg phi=%6.2f deg, err=%5.2f deg, outside8=%6.2f dB\n"
                % (
                    row["case"], row["target_theta_deg"], row["peak_theta_deg"],
                    row["peak_phi_deg"], row["pointing_error_deg"], row["max_outside_8deg_db"],
                )
            )
        f.write("\nZ=200 mm focus maps, level relative to target:\n")
        for row in focus_rows:
            f.write(
                "  %-14s %-10s pressure_PSL=%6.2f dB, audio_proxy=%6.2f dB, "
                "P_FWHM=(%5.1f,%5.1f)mm, audio_FWHM=(%5.1f,%5.1f)mm, gain=%6.2f dB, err=%4.1f mm\n"
                % (
                    row["case"], row["method"], row["pressure_psl_db"], row["audio_source_psl_proxy_db"],
                    row["pressure_fwhm_x_mm"], row["pressure_fwhm_y_mm"],
                    row["audio_source_fwhm_x_mm"], row["audio_source_fwhm_y_mm"],
                    row["target_gain_drop_db"], row["focus_error_mm"],
                )
            )
        f.write("\nUsable radius by z, criteria: pressure PSL <= -10 dB, gain drop >= -5 dB, focus error <= 5 mm:\n")
        for row in usable_rows:
            f.write(
                "  %-10s z=%3.0f mm: radius=%5s mm, angle=%5s deg\n"
                % (
                    row["method"], row["z_mm"],
                    "nan" if math.isnan(float(row["usable_radius_mm"])) else "%.1f" % row["usable_radius_mm"],
                    "nan" if math.isnan(float(row["usable_angle_deg"])) else "%.1f" % row["usable_angle_deg"],
                )
            )


def main():
    ensure_out_dir()
    positions = n32_positions()
    circular_rows, circular_path = circular_plot_and_metrics(positions)
    focus_rows, heatmap_path = focus_heatmaps(positions)
    scan_rows, usable_rows, scan_path, usable_path = scan_audio_focus(positions)
    metrics_path = os.path.join(OUT_DIR, "n32_baseline_metrics.csv")
    write_metrics_csv(metrics_path, circular_rows, focus_rows)
    summary_path = os.path.join(OUT_DIR, "n32_baseline_summary.txt")
    write_summary(summary_path, positions, circular_rows, focus_rows, usable_rows)
    for path in [summary_path, metrics_path, circular_path, heatmap_path, scan_path, usable_path]:
        print(path)


if __name__ == "__main__":
    main()
