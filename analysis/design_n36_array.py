import csv
import math
import os
import random

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
ELEMENT_DIAMETER_M = 16e-3
ELEMENT_RADIUS_M = ELEMENT_DIAMETER_M / 2.0
MIN_CENTER_M = 16.5e-3

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
OUT_DIR = os.path.join(ROOT, "analysis_outputs")
CURRENT_CSV = os.path.join(HERE, "n32_array_coordinates.csv")

DB_FLOOR = -45.0
THETA_MAX_DEG = 55.0
MAIN_EXCLUDE_DEG = 8.0
FOCUS_Z_M = 0.200
OFFSET_X_M = 0.040


def ensure_out_dir():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)


def load_current_layout():
    positions = []
    with open(CURRENT_CSV, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            positions.append([float(row["X_mm"]) * 1e-3, float(row["Y_mm"]) * 1e-3])
    return np.asarray(positions, dtype=float)


def min_distance(positions):
    best = 1e9
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            d = float(np.linalg.norm(positions[i] - positions[j]))
            if d < best:
                best = d
    return best


def aperture_radius(positions):
    return float(np.max(np.sqrt(np.sum(positions * positions, axis=1))))


def element_pattern(theta):
    x = K * ELEMENT_RADIUS_M * np.sin(theta)
    out = np.ones_like(theta, dtype=float)
    mask = np.abs(x) > 1e-10
    out[mask] = np.abs(2.0 * j1(x[mask]) / x[mask])
    return out


def normalize_weights(weights):
    weights = np.asarray(weights, dtype=complex)
    return weights / np.sum(np.abs(weights))


def steer_weights(positions, theta_deg, phi_deg):
    theta = np.deg2rad(theta_deg)
    phi = np.deg2rad(phi_deg)
    phase = K * np.sin(theta) * (positions[:, 0] * np.cos(phi) + positions[:, 1] * np.sin(phi))
    return normalize_weights(np.exp(1j * phase))


def focus_weights(positions, x_focus, y_focus, z_focus):
    r = np.sqrt((x_focus - positions[:, 0]) ** 2 + (y_focus - positions[:, 1]) ** 2 + z_focus**2)
    return normalize_weights(np.exp(-1j * K * r))


def field_on_grid(positions, weights, theta_grid, phi_grid):
    kx = K * np.sin(theta_grid) * np.cos(phi_grid)
    ky = K * np.sin(theta_grid) * np.sin(phi_grid)
    phase = np.outer(kx.ravel(), positions[:, 0]) + np.outer(ky.ravel(), positions[:, 1])
    field = np.dot(np.exp(-1j * phase), weights).reshape(theta_grid.shape)
    return field * element_pattern(theta_grid)


def db_norm(field):
    amp = np.abs(field)
    peak = float(np.max(amp))
    if peak <= 0.0:
        return np.zeros_like(amp) + DB_FLOOR
    return np.maximum(20.0 * np.log10(amp / peak + 1e-12), DB_FLOOR)


def direction_vector(theta_deg, phi_deg):
    theta = np.deg2rad(theta_deg)
    phi = np.deg2rad(phi_deg)
    return np.asarray([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])


def prepare_farfield_grid(theta_step=1.0, phi_step=2.5):
    theta_deg = np.arange(0.0, THETA_MAX_DEG + 0.001, theta_step)
    phi_deg = np.arange(0.0, 360.0, phi_step)
    theta_grid_deg, phi_grid_deg = np.meshgrid(theta_deg, phi_deg)
    theta_grid = np.deg2rad(theta_grid_deg)
    phi_grid = np.deg2rad(phi_grid_deg)
    ux = np.sin(theta_grid) * np.cos(phi_grid)
    uy = np.sin(theta_grid) * np.sin(phi_grid)
    uz = np.cos(theta_grid)
    return theta_grid_deg, phi_grid_deg, theta_grid, phi_grid, ux, uy, uz


def max_outside_main(db, ux, uy, uz, target_theta_deg, target_phi_deg, exclude_deg=MAIN_EXCLUDE_DEG):
    target = direction_vector(target_theta_deg, target_phi_deg)
    sep = np.rad2deg(np.arccos(np.clip(ux * target[0] + uy * target[1] + uz * target[2], -1.0, 1.0)))
    mask = sep > exclude_deg
    return float(np.max(db[mask]))


def peak_pointing_error_deg(db, theta_grid_deg, phi_grid_deg, target_theta_deg, target_phi_deg):
    idx = np.unravel_index(np.argmax(db), db.shape)
    peak = direction_vector(float(theta_grid_deg[idx]), float(phi_grid_deg[idx]))
    target = direction_vector(target_theta_deg, target_phi_deg)
    return float(np.rad2deg(np.arccos(np.clip(np.dot(peak, target), -1.0, 1.0))))


def evaluate_farfield(positions, grid):
    theta_grid_deg, phi_grid_deg, theta_grid, phi_grid, ux, uy, uz = grid
    cases = [
        ("broadside", normalize_weights(np.ones(len(positions))), 0.0, 0.0, 1.00),
        ("steer15_x", steer_weights(positions, 15.0, 0.0), 15.0, 0.0, 1.10),
        ("steer20_x", steer_weights(positions, 20.0, 0.0), 20.0, 0.0, 1.35),
        ("steer20_diag", steer_weights(positions, 20.0, 45.0), 20.0, 45.0, 1.35),
        ("steer25_x", steer_weights(positions, 25.0, 0.0), 25.0, 0.0, 1.55),
        ("focus200_ff", focus_weights(positions, 0.0, 0.0, FOCUS_Z_M), 0.0, 0.0, 0.85),
        (
            "offset_focus_ff",
            focus_weights(positions, OFFSET_X_M, 0.0, FOCUS_Z_M),
            math.degrees(math.atan2(OFFSET_X_M, FOCUS_Z_M)),
            0.0,
            1.10,
        ),
    ]

    weighted = []
    raw = []
    pointing = []
    for name, weights, theta_t, phi_t, weight in cases:
        db = db_norm(field_on_grid(positions, weights, theta_grid, phi_grid))
        sll = max_outside_main(db, ux, uy, uz, theta_t, phi_t)
        err = peak_pointing_error_deg(db, theta_grid_deg, phi_grid_deg, theta_t, phi_t)
        weighted.append(weight * max(0.0, sll + 18.0))
        weighted.append(2.0 * max(0.0, err - 1.5))
        raw.append((name, sll))
        pointing.append((name, err))
    objective = max(weighted) + 0.08 * sum(weighted)
    return objective, raw, pointing


def jittered_ring_layout(counts, radii_m, jitter_frac, rng):
    positions = [[0.0, 0.0]]
    for n, r in zip(counts, radii_m):
        sector = 2.0 * np.pi / n
        offsets = []
        for i in range(n):
            offsets.append((i + rng.uniform(-jitter_frac, jitter_frac)) * sector + rng.uniform(0.0, sector))
        offsets = sorted(offsets)
        for th in offsets:
            positions.append([r * math.cos(th), r * math.sin(th)])
    return np.asarray(positions, dtype=float)


def fermat_layout(n, radius_m, inner_m, rng):
    positions = [[0.0, 0.0]]
    golden = math.pi * (3.0 - math.sqrt(5.0))
    rot = rng.uniform(0.0, 2.0 * math.pi)
    radial_power = rng.uniform(0.48, 0.62)
    angle_jitter = rng.uniform(0.00, 0.16)
    for j in range(n - 1):
        frac = (j + 0.5) / (n - 1)
        r = inner_m + (radius_m - inner_m) * (frac ** radial_power)
        th = rot + j * golden + rng.uniform(-angle_jitter, angle_jitter)
        positions.append([r * math.cos(th), r * math.sin(th)])
    return np.asarray(positions, dtype=float)


def poisson_disk_layout(n, radius_m, rng, max_attempts=25000):
    positions = [[0.0, 0.0]]
    attempts = 0
    while len(positions) < n and attempts < max_attempts:
        attempts += 1
        rr = radius_m * math.sqrt(rng.random())
        th = rng.random() * 2.0 * math.pi
        candidate = np.asarray([rr * math.cos(th), rr * math.sin(th)])
        ok = True
        for p in positions:
            if np.linalg.norm(candidate - np.asarray(p)) < MIN_CENTER_M:
                ok = False
                break
        if ok:
            positions.append(candidate.tolist())
    if len(positions) != n:
        return None
    return np.asarray(positions, dtype=float)


def valid_layout(positions, radius_limit_m):
    if positions is None:
        return False
    if np.max(np.sqrt(np.sum(positions * positions, axis=1))) > radius_limit_m + 1e-9:
        return False
    return min_distance(positions) >= MIN_CENTER_M - 1e-9


def search_layouts(samples_per_family=180, seed=20260510):
    rng = random.Random(seed)
    grid = prepare_farfield_grid(theta_step=2.0, phi_step=6.0)
    records = []

    current = load_current_layout()
    obj, raw, pointing = evaluate_farfield(current, grid)
    records.append(("current_n32", obj, current, raw, pointing))

    ring_specs = [
        ("ring36_6_11_18", [6, 11, 18], [(18.2, 20.2), (37.0, 42.0), (60.0, 66.0)], 36),
        ("ring37_6_10_20", [6, 10, 20], [(18.2, 20.5), (36.0, 41.0), (62.0, 67.0)], 37),
        ("ring37_7_11_18", [7, 11, 18], [(19.5, 22.0), (39.0, 44.0), (61.0, 68.0)], 37),
        ("ring40_7_12_20", [7, 12, 20], [(19.5, 22.5), (40.0, 46.0), (64.0, 72.0)], 40),
    ]
    for name, counts, ranges_mm, n_total in ring_specs:
        for _ in range(samples_per_family):
            radii_m = [rng.uniform(lo, hi) * 1e-3 for lo, hi in ranges_mm]
            layout = jittered_ring_layout(counts, radii_m, rng.uniform(0.05, 0.34), rng)
            if not valid_layout(layout, max(radii_m)):
                continue
            obj, raw, pointing = evaluate_farfield(layout, grid)
            records.append((name, obj, layout, raw, pointing))

    for n in [34, 36, 38, 40]:
        for radius_mm in [58.0, 62.0, 66.0, 70.0]:
            for _ in range(samples_per_family // 3):
                layout = fermat_layout(n, radius_mm * 1e-3, rng.uniform(15.5, 21.5) * 1e-3, rng)
                if not valid_layout(layout, radius_mm * 1e-3):
                    continue
                obj, raw, pointing = evaluate_farfield(layout, grid)
                records.append(("fermat_n%d_r%.0f" % (n, radius_mm), obj, layout, raw, pointing))

    for n in [36, 38, 40]:
        for radius_mm in [62.0, 66.0, 70.0]:
            for _ in range(max(36, samples_per_family // 5)):
                layout = poisson_disk_layout(n, radius_mm * 1e-3, rng)
                if not valid_layout(layout, radius_mm * 1e-3):
                    continue
                obj, raw, pointing = evaluate_farfield(layout, grid)
                records.append(("poisson_n%d_r%.0f" % (n, radius_mm), obj, layout, raw, pointing))

    records.sort(key=lambda item: item[1])
    return records


def signed_cut_db(positions, weights, alpha_deg):
    alpha = np.deg2rad(alpha_deg)
    phase = K * np.outer(np.sin(alpha), positions[:, 0])
    field = np.dot(np.exp(-1j * phase), weights) * element_pattern(np.abs(alpha))
    return db_norm(field)


def fwhm(axis, db, idx_peak):
    threshold = db[idx_peak] - 3.0
    left = idx_peak
    while left > 0 and db[left] >= threshold:
        left -= 1
    right = idx_peak
    while right < len(db) - 1 and db[right] >= threshold:
        right += 1
    if left == 0 or right == len(db) - 1:
        return float("nan")

    def interp(i0, i1):
        if abs(db[i1] - db[i0]) < 1e-12:
            return axis[i0]
        return axis[i0] + (threshold - db[i0]) * (axis[i1] - axis[i0]) / (db[i1] - db[i0])

    return float(interp(left, left + 1) - interp(right, right - 1)) * -1.0


def pressure_field(positions, weights, x_grid, y_grid, z_grid):
    field = np.zeros_like(x_grid, dtype=complex)
    for idx in range(len(positions)):
        x0, y0 = positions[idx]
        dx = x_grid - x0
        dy = y_grid - y0
        r = np.sqrt(dx * dx + dy * dy + z_grid * z_grid)
        theta = np.arctan2(np.sqrt(dx * dx + dy * dy), z_grid)
        field += weights[idx] * element_pattern(theta) * np.exp(1j * K * r) / np.maximum(r, 1e-9)
    return field


def write_layout_csv(positions, path):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Element_ID", "r_mm", "theta_deg", "X_mm", "Y_mm"])
        for i, (x, y) in enumerate(positions):
            r = math.hypot(x, y)
            th = math.degrees(math.atan2(y, x)) % 360.0
            writer.writerow(["E%02d" % i, "%.3f" % (r * 1000.0), "%.3f" % th, "%.3f" % (x * 1000.0), "%.3f" % (y * 1000.0)])


def plot_layout(positions, path):
    fig, ax = plt.subplots(figsize=(7.0, 7.0))
    ax.scatter(positions[:, 0] * 1000.0, positions[:, 1] * 1000.0, s=32, color="#1261a0")
    for i, (x, y) in enumerate(positions):
        ax.add_patch(Circle((x * 1000.0, y * 1000.0), ELEMENT_RADIUS_M * 1000.0, fill=False, lw=0.7, alpha=0.45))
        ax.text(x * 1000.0, y * 1000.0, str(i), fontsize=6, ha="center", va="center")
    ax.set_aspect("equal", adjustable="box")
    r = aperture_radius(positions) * 1000.0
    ax.set_xlim(-r - 12.0, r + 12.0)
    ax.set_ylim(-r - 12.0, r + 12.0)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("Recommended aperiodic planar array")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_report(positions, current, path_prefix):
    grid = prepare_farfield_grid(theta_step=0.5, phi_step=1.5)
    theta_grid_deg, phi_grid_deg, theta_grid, phi_grid, ux, uy, uz = grid
    axis_deg = np.linspace(-THETA_MAX_DEG, THETA_MAX_DEG, 241)
    xang, yang = np.meshgrid(axis_deg, axis_deg)
    theta_plot_deg = np.sqrt(xang * xang + yang * yang)
    phi_plot_deg = (np.rad2deg(np.arctan2(yang, xang)) + 360.0) % 360.0
    plot_mask = theta_plot_deg <= THETA_MAX_DEG
    theta_plot = np.deg2rad(np.minimum(theta_plot_deg, THETA_MAX_DEG))
    phi_plot = np.deg2rad(phi_plot_deg)

    cases = [
        ("Broadside", normalize_weights(np.ones(len(positions))), 0.0, 0.0),
        ("Steer 20 deg", steer_weights(positions, 20.0, 0.0), 20.0, 0.0),
        ("Focus z=200 mm", focus_weights(positions, 0.0, 0.0, FOCUS_Z_M), 0.0, 0.0),
        (
            "Focus x=40, z=200 mm",
            focus_weights(positions, OFFSET_X_M, 0.0, FOCUS_Z_M),
            math.degrees(math.atan2(OFFSET_X_M, FOCUS_Z_M)),
            0.0,
        ),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 10.4))
    levels = np.linspace(DB_FLOOR, 0.0, 46)
    contour = None
    for ax, (title, weights, theta_t, phi_t) in zip(axes.ravel(), cases):
        field = field_on_grid(positions, weights, theta_plot, phi_plot)
        field = np.where(plot_mask, field, 0.0)
        db = db_norm(field)
        db[~plot_mask] = np.nan
        contour = ax.contourf(xang, yang, db, levels=levels, cmap="viridis", extend="min")
        marker_x = theta_t * math.cos(math.radians(phi_t))
        marker_y = theta_t * math.sin(math.radians(phi_t))
        ax.plot([marker_x], [marker_y], "rx", ms=8, mew=1.8)
        ax.add_patch(Circle((0.0, 0.0), THETA_MAX_DEG, fill=False, color="white", lw=0.8, alpha=0.75))
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(-THETA_MAX_DEG, THETA_MAX_DEG)
        ax.set_ylim(-THETA_MAX_DEG, THETA_MAX_DEG)
        ax.set_title(title)
        ax.set_xlabel("theta_x (deg)")
        ax.set_ylabel("theta_y (deg)")
        ax.grid(True, color="white", alpha=0.12, lw=0.5)
    cbar = fig.colorbar(contour, ax=axes.ravel().tolist(), shrink=0.86, pad=0.03)
    cbar.set_label("normalized level (dB)")
    fig.suptitle("Recommended layout: circular normalized far-field")
    fig.savefig(path_prefix + "_circular_patterns.png", dpi=180)
    plt.close(fig)

    alpha = np.linspace(-55.0, 55.0, 1101)
    fig, ax = plt.subplots(figsize=(10.5, 6.0))
    cut_rows = []
    for title, weights, theta_t, phi_t in cases:
        db = signed_cut_db(positions, weights, alpha)
        ax.plot(alpha, db, lw=1.7, label=title)
        outside = np.abs(alpha - theta_t) > MAIN_EXCLUDE_DEG
        sll = float(np.max(db[outside]))
        width = fwhm(alpha, db, int(np.argmax(db)))
        cut_rows.append((title, sll, width))
    ax.set_xlim(-55, 55)
    ax.set_ylim(DB_FLOOR, 2)
    ax.set_xlabel("signed angle in x-z plane (deg)")
    ax.set_ylabel("normalized level (dB)")
    ax.set_title("Recommended layout: x-z far-field cuts")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(path_prefix + "_cuts.png", dpi=180)
    plt.close(fig)

    x_mm = np.linspace(-120.0, 120.0, 241)
    y_mm = np.linspace(-120.0, 120.0, 241)
    x, y = np.meshgrid(x_mm * 1e-3, y_mm * 1e-3)
    z = np.zeros_like(x) + FOCUS_Z_M
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.3))
    heat_cases = [
        ("Focus z=200 mm", focus_weights(positions, 0.0, 0.0, FOCUS_Z_M), 0.0, 0.0),
        ("Focus x=40, z=200 mm", focus_weights(positions, OFFSET_X_M, 0.0, FOCUS_Z_M), 40.0, 0.0),
    ]
    heat_rows = []
    contour = None
    for ax, (title, weights, tx, ty) in zip(axes.ravel(), heat_cases):
        db = db_norm(pressure_field(positions, weights, x, y, z))
        contour = ax.contourf(x_mm, y_mm, db, levels=levels, cmap="magma", extend="min")
        idx = np.unravel_index(np.argmax(db), db.shape)
        px = float(x_mm[idx[1]])
        py = float(y_mm[idx[0]])
        wx = fwhm(x_mm, db[idx[0], :], idx[1])
        wy = fwhm(y_mm, db[:, idx[1]], idx[0])
        heat_rows.append((title, px, py, math.hypot(px - tx, py - ty), wx, wy))
        ax.plot([tx], [ty], "cx", ms=8, mew=1.8)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(title)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.grid(True, color="white", alpha=0.12, lw=0.5)
    cbar = fig.colorbar(contour, ax=axes.ravel().tolist(), shrink=0.86, pad=0.03)
    cbar.set_label("normalized pressure (dB)")
    fig.suptitle("Recommended layout: XY pressure at z=200 mm")
    fig.savefig(path_prefix + "_xy_focus_heatmaps.png", dpi=180)
    plt.close(fig)

    current_obj, current_raw, current_point = evaluate_farfield(current, grid)
    best_obj, best_raw, best_point = evaluate_farfield(positions, grid)
    return current_obj, current_raw, current_point, best_obj, best_raw, best_point, cut_rows, heat_rows


def local_refine_layout(positions, iterations=450, seed=20260511, radius_limit_m=70e-3):
    rng = random.Random(seed)
    grid = prepare_farfield_grid(theta_step=2.0, phi_step=6.0)
    best = positions.copy()
    best_obj, best_raw, best_point = evaluate_farfield(best, grid)
    for it in range(iterations):
        cand = best.copy()
        idx = rng.randrange(1, len(cand))
        step = 0.0045 * (1.0 - float(it) / iterations) + 0.0010
        cand[idx, 0] += rng.gauss(0.0, step)
        cand[idx, 1] += rng.gauss(0.0, step)
        r = math.hypot(cand[idx, 0], cand[idx, 1])
        if r > radius_limit_m:
            cand[idx] *= radius_limit_m / r
        if min_distance(cand) < MIN_CENTER_M:
            continue
        obj, raw, point = evaluate_farfield(cand, grid)
        if obj < best_obj:
            best = cand
            best_obj = obj
            best_raw = raw
            best_point = point
    return best, best_obj, best_raw, best_point


def xy_nearfield_sidelobe_metrics(positions, target):
    x_mm = np.linspace(-130.0, 130.0, 261)
    y_mm = np.linspace(-130.0, 130.0, 261)
    x, y = np.meshgrid(x_mm * 1e-3, y_mm * 1e-3)
    z = np.zeros_like(x) + target[2]
    weights = focus_weights(positions, target[0], target[1], target[2])
    db = db_norm(pressure_field(positions, weights, x, y, z))
    idx = np.unravel_index(np.argmax(db), db.shape)
    tx_mm = target[0] * 1000.0
    ty_mm = target[1] * 1000.0
    rr = np.sqrt((x_mm[None, :] - tx_mm) ** 2 + (y_mm[:, None] - ty_mm) ** 2)
    outside = rr > 20.0
    peak_x = float(x_mm[idx[1]])
    peak_y = float(y_mm[idx[0]])
    err = math.hypot(peak_x - tx_mm, peak_y - ty_mm)
    return peak_x, peak_y, err, float(np.max(db[outside])), float(np.percentile(db[outside], 99.5))


def write_final_summary(summary_path, positions, current, current_raw, final_raw, final_point, cut_rows, heat_rows):
    with open(summary_path, "w") as f:
        f.write("recommended_layout: %d-element refined blue-noise / Poisson-disk planar aperture\n" % len(positions))
        f.write("frequency_hz: %.1f\n" % FREQ)
        f.write("wavelength_mm: %.3f\n" % (LAMBDA * 1000.0))
        f.write("element_diameter_mm: %.3f\n" % (ELEMENT_DIAMETER_M * 1000.0))
        f.write("element_count: %d\n" % len(positions))
        f.write("aperture_diameter_mm: %.3f\n" % (2.0 * aperture_radius(positions) * 1000.0))
        f.write("aperture_in_lambda: %.3f\n" % (2.0 * aperture_radius(positions) / LAMBDA))
        f.write("min_center_distance_mm: %.3f\n" % (min_distance(positions) * 1000.0))
        f.write("baseline_current32_aperture_mm: %.3f\n" % (2.0 * aperture_radius(current) * 1000.0))
        f.write("baseline_current32_min_distance_mm: %.3f\n" % (min_distance(current) * 1000.0))
        f.write("\nfarfield_worst_sidelobe_db_current32:\n")
        for name, val in current_raw:
            f.write("  %s: %.2f\n" % (name, val))
        f.write("\nfarfield_worst_sidelobe_db_refined%d:\n" % len(positions))
        for name, val in final_raw:
            f.write("  %s: %.2f\n" % (name, val))
        f.write("\nfarfield_pointing_error_deg_refined%d:\n" % len(positions))
        for name, val in final_point:
            f.write("  %s: %.2f\n" % (name, val))
        f.write("\nxz_cut_rows_title_sll_db_fwhm_deg:\n")
        for row in cut_rows:
            f.write("  %s: %.2f, %.2f\n" % row)
        f.write("\nxy_focus_rows_title_peakx_peaky_error_fwhmx_fwhmy_mm:\n")
        for row in heat_rows:
            f.write("  %s: %.2f, %.2f, %.2f, %.2f, %.2f\n" % row)
        f.write("\nxy_nearfield_psl_outside_20mm_current32:\n")
        for label, target in [("focus200", (0.0, 0.0, FOCUS_Z_M)), ("offset40_200", (OFFSET_X_M, 0.0, FOCUS_Z_M))]:
            peak_x, peak_y, err, psl, p995 = xy_nearfield_sidelobe_metrics(current, target)
            f.write("  %s: peak=(%.1f, %.1f)mm err=%.1fmm psl=%.2fdB p99.5=%.2fdB\n" % (
                label, peak_x, peak_y, err, psl, p995
            ))
        f.write("\nxy_nearfield_psl_outside_20mm_refined%d:\n" % len(positions))
        for label, target in [("focus200", (0.0, 0.0, FOCUS_Z_M)), ("offset40_200", (OFFSET_X_M, 0.0, FOCUS_Z_M))]:
            peak_x, peak_y, err, psl, p995 = xy_nearfield_sidelobe_metrics(positions, target)
            f.write("  %s: peak=(%.1f, %.1f)mm err=%.1fmm psl=%.2fdB p99.5=%.2fdB\n" % (
                label, peak_x, peak_y, err, psl, p995
            ))


def main():
    ensure_out_dir()
    records = search_layouts(samples_per_family=180)
    best_name, best_obj, best_positions, best_raw, best_pointing = records[0]
    current = load_current_layout()

    csv_path = os.path.join(OUT_DIR, "n36_recommended_layout.csv")
    layout_png = os.path.join(OUT_DIR, "n36_recommended_layout.png")
    write_layout_csv(best_positions, csv_path)
    plot_layout(best_positions, layout_png)

    prefix = os.path.join(OUT_DIR, "n36_recommended")
    report = plot_report(best_positions, current, prefix)
    current_obj, current_raw, current_point, final_obj, final_raw, final_point, cut_rows, heat_rows = report

    summary_path = os.path.join(OUT_DIR, "n36_recommended_summary.txt")
    with open(summary_path, "w") as f:
        f.write("recommended_family: %s\n" % best_name)
        f.write("element_count: %d\n" % len(best_positions))
        f.write("frequency_hz: %.1f\n" % FREQ)
        f.write("wavelength_mm: %.3f\n" % (LAMBDA * 1000.0))
        f.write("element_diameter_mm: %.3f\n" % (ELEMENT_DIAMETER_M * 1000.0))
        f.write("aperture_diameter_mm: %.3f\n" % (2.0 * aperture_radius(best_positions) * 1000.0))
        f.write("aperture_in_lambda: %.3f\n" % (2.0 * aperture_radius(best_positions) / LAMBDA))
        f.write("min_center_distance_mm: %.3f\n" % (min_distance(best_positions) * 1000.0))
        f.write("search_objective_current_n32: %.4f\n" % current_obj)
        f.write("search_objective_recommended: %.4f\n" % final_obj)
        f.write("\nfarfield_worst_sidelobe_db_recommended:\n")
        for name, val in final_raw:
            f.write("  %s: %.2f\n" % (name, val))
        f.write("\nfarfield_pointing_error_deg_recommended:\n")
        for name, val in final_point:
            f.write("  %s: %.2f\n" % (name, val))
        f.write("\nxz_cut_rows_title_sll_db_fwhm_deg:\n")
        for row in cut_rows:
            f.write("  %s: %.2f, %.2f\n" % row)
        f.write("\nxy_focus_rows_title_peakx_peaky_error_fwhmx_fwhmy_mm:\n")
        for row in heat_rows:
            f.write("  %s: %.2f, %.2f, %.2f, %.2f, %.2f\n" % row)
        f.write("\ntop_candidates:\n")
        for item in records[:20]:
            f.write("  %s obj=%.4f n=%d aperture=%.1fmm min_dist=%.2fmm\n" % (
                item[0],
                item[1],
                len(item[2]),
                2.0 * aperture_radius(item[2]) * 1000.0,
                min_distance(item[2]) * 1000.0,
            ))

    refined_positions, refined_obj, refined_raw, refined_point = local_refine_layout(best_positions)
    refined_prefix = os.path.join(OUT_DIR, "n%d_refined" % len(refined_positions))
    refined_csv_path = refined_prefix + "_layout.csv"
    refined_layout_png = refined_prefix + "_layout.png"
    write_layout_csv(refined_positions, refined_csv_path)
    plot_layout(refined_positions, refined_layout_png)
    refined_report = plot_report(refined_positions, current, refined_prefix)
    _, refined_current_raw, _, _, refined_final_raw, refined_final_point, refined_cut_rows, refined_heat_rows = refined_report
    refined_summary_path = refined_prefix + "_summary.txt"
    write_final_summary(
        refined_summary_path,
        refined_positions,
        current,
        refined_current_raw,
        refined_final_raw,
        refined_final_point,
        refined_cut_rows,
        refined_heat_rows,
    )

    print("Generated:")
    for path in [
        csv_path,
        layout_png,
        prefix + "_circular_patterns.png",
        prefix + "_cuts.png",
        prefix + "_xy_focus_heatmaps.png",
        summary_path,
        refined_csv_path,
        refined_layout_png,
        refined_prefix + "_circular_patterns.png",
        refined_prefix + "_cuts.png",
        refined_prefix + "_xy_focus_heatmaps.png",
        refined_summary_path,
    ]:
        print(path)


if __name__ == "__main__":
    main()
