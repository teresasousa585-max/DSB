import csv
import math
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from scipy.special import j1


C = 343.0
FREQ = 40000.0
LAMBDA = C / FREQ
K = 2.0 * np.pi / LAMBDA
ELEMENT_DIAMETER = 16e-3
ELEMENT_RADIUS = ELEMENT_DIAMETER / 2.0

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "analysis_outputs")
CSV_PATH = os.path.join(os.path.dirname(__file__), "n32_array_coordinates.csv")

STEER_DEG = 20.0
FOCUS_Z = 0.200
OFFSET_FOCUS_X = 0.040
PLOT_THETA_MAX_DEG = 60.0
MAIN_EXCLUDE_DEG = 8.0
DB_FLOOR = -40.0


def load_positions(csv_path):
    ids = []
    rings = []
    pos = []
    with open(csv_path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ids.append(row["Element_ID"])
            rings.append(row["Ring"])
            pos.append([float(row["X_mm"]) * 1e-3, float(row["Y_mm"]) * 1e-3])
    return ids, rings, np.asarray(pos, dtype=float)


def element_pattern(theta):
    x = K * ELEMENT_RADIUS * np.sin(theta)
    out = np.ones_like(x, dtype=float)
    mask = np.abs(x) > 1e-10
    out[mask] = np.abs(2.0 * j1(x[mask]) / x[mask])
    return out


def normalize_complex(weights):
    weights = np.asarray(weights, dtype=complex)
    return weights / np.sum(np.abs(weights))


def steer_weights(positions, theta_deg, phi_deg):
    theta = np.deg2rad(theta_deg)
    phi = np.deg2rad(phi_deg)
    phase = K * np.sin(theta) * (positions[:, 0] * np.cos(phi) + positions[:, 1] * np.sin(phi))
    return normalize_complex(np.exp(1j * phase))


def focus_weights(positions, x_focus, y_focus, z_focus):
    r = np.sqrt((x_focus - positions[:, 0]) ** 2 + (y_focus - positions[:, 1]) ** 2 + z_focus**2)
    r_ref = np.sqrt(x_focus**2 + y_focus**2 + z_focus**2)
    return normalize_complex(np.exp(-1j * K * (r - r_ref)))


def farfield_complex(positions, weights, theta, phi):
    kx = K * np.sin(theta) * np.cos(phi)
    ky = K * np.sin(theta) * np.sin(phi)
    phase = np.outer(kx.ravel(), positions[:, 0]) + np.outer(ky.ravel(), positions[:, 1])
    field = np.dot(np.exp(-1j * phase), weights).reshape(theta.shape)
    return field * element_pattern(theta)


def db_norm(field):
    amp = np.abs(field)
    peak = np.max(amp)
    if peak <= 0.0:
        return np.zeros_like(amp) + DB_FLOOR
    return np.maximum(20.0 * np.log10(amp / peak + 1e-12), DB_FLOOR)


def direction_vector(theta_deg, phi_deg):
    theta = np.deg2rad(theta_deg)
    phi = np.deg2rad(phi_deg)
    return np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])


def circular_pattern_metrics(db, theta_grid_deg, phi_grid_deg, target_theta_deg, target_phi_deg):
    target = direction_vector(target_theta_deg, target_phi_deg)
    u_x = np.sin(np.deg2rad(theta_grid_deg)) * np.cos(np.deg2rad(phi_grid_deg))
    u_y = np.sin(np.deg2rad(theta_grid_deg)) * np.sin(np.deg2rad(phi_grid_deg))
    u_z = np.cos(np.deg2rad(theta_grid_deg))
    dot = np.clip(u_x * target[0] + u_y * target[1] + u_z * target[2], -1.0, 1.0)
    sep_deg = np.rad2deg(np.arccos(dot))
    outside = sep_deg > MAIN_EXCLUDE_DEG
    max_outside = np.max(db[outside]) if np.any(outside) else np.nan

    max_idx = np.unravel_index(np.argmax(db), db.shape)
    peak_theta = theta_grid_deg[max_idx]
    peak_phi = phi_grid_deg[max_idx]
    peak_u = direction_vector(float(peak_theta), float(peak_phi))
    pointing_error = np.rad2deg(np.arccos(np.clip(np.dot(peak_u, target), -1.0, 1.0)))
    return float(max_outside), float(peak_theta), float(peak_phi % 360.0), float(pointing_error)


def hemispherical_directivity_db(positions, weights):
    theta = np.deg2rad(np.linspace(0.25, 89.75, 180))
    phi = np.deg2rad(np.linspace(0.0, 359.0, 360))
    theta_grid, phi_grid = np.meshgrid(theta, phi)
    field = farfield_complex(positions, weights, theta_grid, phi_grid)
    power = np.abs(field) ** 2
    sin_theta = np.sin(theta_grid)
    avg_power = np.sum(power * sin_theta) / np.sum(sin_theta)
    peak_power = np.max(power)
    return float(10.0 * np.log10(peak_power / avg_power + 1e-12))


def signed_cut(positions, weights, alpha_deg):
    alpha = np.deg2rad(alpha_deg)
    theta = np.abs(alpha)
    phase = K * np.outer(np.sin(alpha), positions[:, 0])
    field = np.dot(np.exp(-1j * phase), weights) * element_pattern(theta)
    return field


def interpolate_crossing(x0, y0, x1, y1, threshold):
    if abs(y1 - y0) < 1e-12:
        return x0
    return x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)


def beamwidth_3db(alpha_deg, db_cut, target_alpha_deg):
    window = np.abs(alpha_deg - target_alpha_deg) <= 6.0
    if not np.any(window):
        idx_peak = int(np.argmax(db_cut))
    else:
        idx_candidates = np.where(window)[0]
        idx_peak = int(idx_candidates[np.argmax(db_cut[idx_candidates])])
    threshold = db_cut[idx_peak] - 3.0

    left = idx_peak
    while left > 0 and db_cut[left] >= threshold:
        left -= 1
    right = idx_peak
    while right < len(db_cut) - 1 and db_cut[right] >= threshold:
        right += 1

    if left == 0 or right == len(db_cut) - 1:
        return np.nan

    x_left = interpolate_crossing(alpha_deg[left], db_cut[left], alpha_deg[left + 1], db_cut[left + 1], threshold)
    x_right = interpolate_crossing(alpha_deg[right - 1], db_cut[right - 1], alpha_deg[right], db_cut[right], threshold)
    return float(x_right - x_left)


def pressure_field(positions, weights, x_grid, y_grid, z_grid):
    field = np.zeros_like(x_grid, dtype=complex)
    for idx in range(positions.shape[0]):
        x0, y0 = positions[idx]
        dx = x_grid - x0
        dy = y_grid - y0
        r = np.sqrt(dx * dx + dy * dy + z_grid * z_grid)
        theta = np.arctan2(np.sqrt(dx * dx + dy * dy), z_grid)
        field += weights[idx] * element_pattern(theta) * np.exp(1j * K * r) / np.maximum(r, 1e-9)
    return field


def fwhm_mm(axis_mm, db_line, idx_peak):
    threshold = db_line[idx_peak] - 3.0
    left = idx_peak
    while left > 0 and db_line[left] >= threshold:
        left -= 1
    right = idx_peak
    while right < len(db_line) - 1 and db_line[right] >= threshold:
        right += 1
    if left == 0 or right == len(db_line) - 1:
        return np.nan
    x_left = interpolate_crossing(axis_mm[left], db_line[left], axis_mm[left + 1], db_line[left + 1], threshold)
    x_right = interpolate_crossing(axis_mm[right - 1], db_line[right - 1], axis_mm[right], db_line[right], threshold)
    return float(x_right - x_left)


def plot_layout(ids, rings, positions):
    fig, ax = plt.subplots(figsize=(7.0, 7.0))
    colors = {"Center": "#1f77b4", "Ring1": "#2ca02c", "Ring2": "#ff7f0e", "Ring3": "#d62728"}
    for ring in sorted(set(rings)):
        mask = np.asarray([r == ring for r in rings])
        ax.scatter(positions[mask, 0] * 1000.0, positions[mask, 1] * 1000.0, s=35, label=ring, color=colors.get(ring, None))
    for x, y in positions:
        ax.add_patch(Circle((x * 1000.0, y * 1000.0), ELEMENT_RADIUS * 1000.0, fill=False, lw=0.7, alpha=0.45))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("N32 array layout, element diameter 16 mm")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper right")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "n32_00_array_layout.png")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_circular_patterns(positions, cases):
    theta_deg = np.linspace(0.0, PLOT_THETA_MAX_DEG, 241)
    phi_deg = np.linspace(0.0, 360.0, 361)
    theta_grid_deg, phi_grid_deg = np.meshgrid(theta_deg, phi_deg)
    theta_grid = np.deg2rad(theta_grid_deg)
    phi_grid = np.deg2rad(phi_grid_deg)
    xx = theta_grid_deg * np.cos(phi_grid)
    yy = theta_grid_deg * np.sin(phi_grid)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 10.5))
    axes = axes.ravel()
    metrics = {}
    levels = np.linspace(DB_FLOOR, 0.0, 41)
    contour = None
    for ax, case in zip(axes, cases):
        field = farfield_complex(positions, case["weights"], theta_grid, phi_grid)
        db = db_norm(field)
        contour = ax.contourf(xx, yy, db, levels=levels, cmap="viridis", extend="min")
        t = case["target_theta_deg"]
        p = np.deg2rad(case["target_phi_deg"])
        ax.plot([t * np.cos(p)], [t * np.sin(p)], "rx", ms=8, mew=1.8)
        ax.add_patch(Circle((0.0, 0.0), PLOT_THETA_MAX_DEG, fill=False, color="white", lw=0.8, alpha=0.7))
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(-PLOT_THETA_MAX_DEG, PLOT_THETA_MAX_DEG)
        ax.set_ylim(-PLOT_THETA_MAX_DEG, PLOT_THETA_MAX_DEG)
        ax.set_title(case["title"])
        ax.set_xlabel("theta_x (deg)")
        ax.set_ylabel("theta_y (deg)")
        ax.grid(True, color="white", alpha=0.15, lw=0.5)
        metrics[case["name"]] = circular_pattern_metrics(
            db, theta_grid_deg, phi_grid_deg, case["target_theta_deg"], case["target_phi_deg"]
        )
    if contour is not None:
        cbar = fig.colorbar(contour, ax=axes.tolist(), shrink=0.86, pad=0.03)
        cbar.set_label("normalized level (dB)")
    fig.suptitle("Circular normalized far-field patterns, 40 kHz")
    fig.savefig(os.path.join(OUT_DIR, "n32_10_circular_normalized_patterns.png"), dpi=180)
    plt.close(fig)
    return metrics


def plot_farfield_cuts(positions, cases):
    alpha_deg = np.linspace(-60.0, 60.0, 1201)
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    cut_metrics = {}
    for case in cases:
        field = signed_cut(positions, case["weights"], alpha_deg)
        db = db_norm(field)
        ax.plot(alpha_deg, db, lw=1.8, label=case["title"])
        target = case["target_signed_deg"]
        ax.axvline(target, color="k", lw=0.6, alpha=0.20)
        bw = beamwidth_3db(alpha_deg, db, target)
        outside = np.abs(alpha_deg - target) > MAIN_EXCLUDE_DEG
        side = float(np.max(db[outside])) if np.any(outside) else np.nan
        cut_metrics[case["name"]] = (bw, side)
    ax.set_xlim(-60, 60)
    ax.set_ylim(DB_FLOOR, 2)
    ax.set_xlabel("signed angle in x-z plane (deg)")
    ax.set_ylabel("normalized level (dB)")
    ax.set_title("Far-field cuts in x-z plane")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower right")
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "n32_11_farfield_cuts.png")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return cut_metrics


def plot_xy_heatmaps(positions, cases):
    x_mm = np.linspace(-120.0, 120.0, 241)
    y_mm = np.linspace(-120.0, 120.0, 241)
    x, y = np.meshgrid(x_mm * 1e-3, y_mm * 1e-3)
    z = np.zeros_like(x) + FOCUS_Z

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 10.2))
    axes = axes.ravel()
    levels = np.linspace(DB_FLOOR, 0.0, 41)
    contour = None
    metrics = {}
    for ax, case in zip(axes, cases):
        field = pressure_field(positions, case["weights"], x, y, z)
        db = db_norm(field)
        contour = ax.contourf(x_mm, y_mm, db, levels=levels, cmap="magma", extend="min")
        target_x_mm = case["target_x_m"] * 1000.0
        target_y_mm = case["target_y_m"] * 1000.0
        ax.plot([target_x_mm], [target_y_mm], "cx", ms=8, mew=1.8)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(case["title"])
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.grid(True, color="white", alpha=0.12, lw=0.5)

        peak_idx = np.unravel_index(np.argmax(db), db.shape)
        peak_x_mm = float(x_mm[peak_idx[1]])
        peak_y_mm = float(y_mm[peak_idx[0]])
        err_mm = math.hypot(peak_x_mm - target_x_mm, peak_y_mm - target_y_mm)
        width_x = fwhm_mm(x_mm, db[peak_idx[0], :], peak_idx[1])
        width_y = fwhm_mm(y_mm, db[:, peak_idx[1]], peak_idx[0])
        metrics[case["name"]] = (peak_x_mm, peak_y_mm, err_mm, width_x, width_y)
    if contour is not None:
        cbar = fig.colorbar(contour, ax=axes.tolist(), shrink=0.86, pad=0.03)
        cbar.set_label("normalized pressure (dB)")
    fig.suptitle("XY pressure heatmaps at z = 200 mm")
    fig.savefig(os.path.join(OUT_DIR, "n32_20_xy_pressure_heatmaps.png"), dpi=180)
    plt.close(fig)
    return metrics


def plot_xz_heatmaps(positions, cases):
    x_mm = np.linspace(-120.0, 120.0, 241)
    z_mm = np.linspace(60.0, 450.0, 261)
    x, z = np.meshgrid(x_mm * 1e-3, z_mm * 1e-3)
    y = np.zeros_like(x)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 10.2))
    axes = axes.ravel()
    levels = np.linspace(DB_FLOOR, 0.0, 41)
    contour = None
    metrics = {}
    for ax, case in zip(axes, cases):
        field = pressure_field(positions, case["weights"], x, y, z)
        db = db_norm(field)
        contour = ax.contourf(x_mm, z_mm, db, levels=levels, cmap="magma", extend="min")
        target_x_mm = case["target_x_m"] * 1000.0
        target_z_mm = case["target_z_m"] * 1000.0
        ax.plot([target_x_mm], [target_z_mm], "cx", ms=8, mew=1.8)
        ax.set_title(case["title"])
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("z (mm)")
        ax.grid(True, color="white", alpha=0.12, lw=0.5)

        peak_idx = np.unravel_index(np.argmax(db), db.shape)
        peak_x_mm = float(x_mm[peak_idx[1]])
        peak_z_mm = float(z_mm[peak_idx[0]])
        err_mm = math.hypot(peak_x_mm - target_x_mm, peak_z_mm - target_z_mm)
        metrics[case["name"]] = (peak_x_mm, peak_z_mm, err_mm)
    if contour is not None:
        cbar = fig.colorbar(contour, ax=axes.tolist(), shrink=0.86, pad=0.03)
        cbar.set_label("normalized pressure (dB)")
    fig.suptitle("XZ pressure heatmaps, y = 0")
    fig.savefig(os.path.join(OUT_DIR, "n32_21_xz_pressure_heatmaps.png"), dpi=180)
    plt.close(fig)
    return metrics


def write_metrics(positions, cases, circular_metrics, cut_metrics, xy_metrics, xz_metrics):
    pairwise = []
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            pairwise.append(np.linalg.norm(positions[i] - positions[j]))
    pairwise = np.asarray(pairwise)
    aperture_m = 2.0 * np.max(np.sqrt(np.sum(positions * positions, axis=1)))
    element_first_null_deg = math.degrees(math.asin(min(1.0, 3.831705970 / (K * ELEMENT_RADIUS))))

    metrics_path = os.path.join(OUT_DIR, "n32_metrics.csv")
    with open(metrics_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "case",
                "target_theta_deg",
                "peak_theta_deg",
                "peak_phi_deg",
                "pointing_error_deg",
                "max_outside_8deg_db",
                "cut_3db_width_deg",
                "cut_max_outside_8deg_db",
                "hemisphere_directivity_db",
                "xy_peak_x_mm",
                "xy_peak_y_mm",
                "xy_target_error_mm",
                "xy_fwhm_x_mm",
                "xy_fwhm_y_mm",
                "xz_peak_x_mm",
                "xz_peak_z_mm",
                "xz_target_error_mm",
            ]
        )
        for case in cases:
            name = case["name"]
            c_side, peak_theta, peak_phi, pointing = circular_metrics[name]
            bw, cut_side = cut_metrics[name]
            xy_peak_x, xy_peak_y, xy_err, xy_wx, xy_wy = xy_metrics[name]
            xz_peak_x, xz_peak_z, xz_err = xz_metrics[name]
            writer.writerow(
                [
                    name,
                    case["target_theta_deg"],
                    peak_theta,
                    peak_phi,
                    pointing,
                    c_side,
                    bw,
                    cut_side,
                    hemispherical_directivity_db(positions, case["weights"]),
                    xy_peak_x,
                    xy_peak_y,
                    xy_err,
                    xy_wx,
                    xy_wy,
                    xz_peak_x,
                    xz_peak_z,
                    xz_err,
                ]
            )

    summary_path = os.path.join(OUT_DIR, "n32_summary.txt")
    with open(summary_path, "w") as f:
        f.write("N32 sparse circular-ring array summary\n")
        f.write("frequency_hz: %.1f\n" % FREQ)
        f.write("wavelength_mm: %.3f\n" % (LAMBDA * 1000.0))
        f.write("element_diameter_mm: %.3f\n" % (ELEMENT_DIAMETER * 1000.0))
        f.write("aperture_diameter_mm: %.3f\n" % (aperture_m * 1000.0))
        f.write("aperture_in_wavelengths: %.3f\n" % (aperture_m / LAMBDA))
        f.write("min_center_distance_mm: %.3f\n" % (np.min(pairwise) * 1000.0))
        f.write("min_center_distance_in_wavelengths: %.3f\n" % (np.min(pairwise) / LAMBDA))
        f.write("element_first_null_deg: %.2f\n" % element_first_null_deg)
        f.write("main_exclusion_for_sidelobe_deg: %.2f\n" % MAIN_EXCLUDE_DEG)
        f.write("\nSee n32_metrics.csv for per-case values.\n")
    return metrics_path, summary_path


def main():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)
    ids, rings, positions = load_positions(CSV_PATH)
    broadside = normalize_complex(np.ones(positions.shape[0], dtype=complex))
    steer = steer_weights(positions, STEER_DEG, 0.0)
    focus = focus_weights(positions, 0.0, 0.0, FOCUS_Z)
    offset_focus = focus_weights(positions, OFFSET_FOCUS_X, 0.0, FOCUS_Z)
    offset_focus_theta = math.degrees(math.atan2(OFFSET_FOCUS_X, FOCUS_Z))
    steer_x = FOCUS_Z * math.tan(math.radians(STEER_DEG))

    cases = [
        {
            "name": "broadside",
            "title": "Broadside",
            "weights": broadside,
            "target_theta_deg": 0.0,
            "target_phi_deg": 0.0,
            "target_signed_deg": 0.0,
            "target_x_m": 0.0,
            "target_y_m": 0.0,
            "target_z_m": FOCUS_Z,
        },
        {
            "name": "steer_20deg",
            "title": "Steer 20 deg",
            "weights": steer,
            "target_theta_deg": STEER_DEG,
            "target_phi_deg": 0.0,
            "target_signed_deg": STEER_DEG,
            "target_x_m": steer_x,
            "target_y_m": 0.0,
            "target_z_m": FOCUS_Z,
        },
        {
            "name": "focus_200mm",
            "title": "Focus z=200 mm",
            "weights": focus,
            "target_theta_deg": 0.0,
            "target_phi_deg": 0.0,
            "target_signed_deg": 0.0,
            "target_x_m": 0.0,
            "target_y_m": 0.0,
            "target_z_m": FOCUS_Z,
        },
        {
            "name": "offset_focus_40mm_200mm",
            "title": "Focus x=40, z=200 mm",
            "weights": offset_focus,
            "target_theta_deg": offset_focus_theta,
            "target_phi_deg": 0.0,
            "target_signed_deg": offset_focus_theta,
            "target_x_m": OFFSET_FOCUS_X,
            "target_y_m": 0.0,
            "target_z_m": FOCUS_Z,
        },
    ]

    layout_path = plot_layout(ids, rings, positions)
    circular_metrics = plot_circular_patterns(positions, cases)
    cut_metrics = plot_farfield_cuts(positions, cases)
    xy_metrics = plot_xy_heatmaps(positions, cases)
    xz_metrics = plot_xz_heatmaps(positions, cases)
    metrics_path, summary_path = write_metrics(positions, cases, circular_metrics, cut_metrics, xy_metrics, xz_metrics)

    print("Generated:")
    for path in [
        layout_path,
        os.path.join(OUT_DIR, "n32_10_circular_normalized_patterns.png"),
        os.path.join(OUT_DIR, "n32_11_farfield_cuts.png"),
        os.path.join(OUT_DIR, "n32_20_xy_pressure_heatmaps.png"),
        os.path.join(OUT_DIR, "n32_21_xz_pressure_heatmaps.png"),
        metrics_path,
        summary_path,
    ]:
        print(path)


if __name__ == "__main__":
    main()
