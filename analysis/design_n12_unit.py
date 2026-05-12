import csv
import math
import os
import random

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
ELEMENT_DIAMETER_MM = 16.0
ELEMENT_RADIUS_M = ELEMENT_DIAMETER_MM * 0.5e-3
MIN_CENTER_SPACING_MM = 17.2
DB_FLOOR = -55.0


def ensure_out_dir():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)


def piston_from_cos(cos_theta):
    cos_theta = np.asarray(cos_theta)
    front = cos_theta > 0.0
    sin_theta = np.sqrt(np.maximum(0.0, 1.0 - cos_theta * cos_theta))
    x = K * ELEMENT_RADIUS_M * sin_theta
    out = np.zeros_like(cos_theta, dtype=float)
    small = np.abs(x) < 1e-10
    out[front & small] = 1.0
    mask = front & (~small)
    out[mask] = np.abs(2.0 * j1(x[mask]) / x[mask])
    return out


def ring_points(radius_mm, count, rotation_deg=0.0):
    pts = []
    for i in range(count):
        a = math.radians(rotation_deg + 360.0 * i / count)
        pts.append([radius_mm * math.cos(a), radius_mm * math.sin(a)])
    return pts


def min_spacing_mm(points):
    pts = np.asarray(points, dtype=float)
    best = 1e9
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            best = min(best, float(np.linalg.norm(pts[i] - pts[j])))
    return best


def aperture_mm(points):
    pts = np.asarray(points, dtype=float)
    return {
        "center_radius_mm": float(np.max(np.linalg.norm(pts, axis=1))),
        "width_centers_mm": float(np.max(pts[:, 0]) - np.min(pts[:, 0])),
        "height_centers_mm": float(np.max(pts[:, 1]) - np.min(pts[:, 1])),
    }


def centered(points):
    pts = np.asarray(points, dtype=float)
    return pts - np.mean(pts, axis=0)


def layout_rect_3x4(px_mm, py_mm):
    xs = [-1.5 * px_mm, -0.5 * px_mm, 0.5 * px_mm, 1.5 * px_mm]
    ys = [-py_mm, 0.0, py_mm]
    return centered([[x, y] for y in ys for x in xs])


def layout_ring_1_4_7(r1_mm, r2_mm, rot4_deg, rot7_deg):
    return centered([[0.0, 0.0]] + ring_points(r1_mm, 4, rot4_deg) + ring_points(r2_mm, 7, rot7_deg))


def layout_ring_5_7(r1_mm, r2_mm, rot5_deg, rot7_deg):
    return centered(ring_points(r1_mm, 5, rot5_deg) + ring_points(r2_mm, 7, rot7_deg))


def layout_sunflower(radius_mm, rotation_deg=0.0):
    golden = math.radians(137.50776405)
    pts = []
    for i in range(12):
        r = radius_mm * math.sqrt((i + 0.5) / 12.0)
        a = rotation_deg + golden * i
        pts.append([r * math.cos(a), r * math.sin(a)])
    return centered(pts)


def random_poisson_disk(radius_mm, min_dist_mm, n=12, attempts=4000):
    pts = []
    for _ in range(attempts):
        r = radius_mm * math.sqrt(random.random())
        a = 2.0 * math.pi * random.random()
        p = np.asarray([r * math.cos(a), r * math.sin(a)])
        if all(np.linalg.norm(p - np.asarray(q)) >= min_dist_mm for q in pts):
            pts.append(p)
            if len(pts) == n:
                return centered(np.asarray(pts))
    return None


def farfield_grid(axis_deg=70.0, n=91):
    axis = np.linspace(-axis_deg, axis_deg, n)
    tx, ty = np.meshgrid(axis, axis)
    theta_deg = np.sqrt(tx * tx + ty * ty)
    phi = np.arctan2(ty, tx)
    mask = theta_deg <= axis_deg
    theta = np.deg2rad(np.minimum(theta_deg, axis_deg))
    dirs = np.stack([
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta),
    ], axis=-1)
    elem = piston_from_cos(dirs[..., 2])
    return axis, tx, ty, theta_deg, dirs, elem, mask


GRID = farfield_grid()


def steering_vector(points_mm, target_theta_x_deg=0.0, target_theta_y_deg=0.0):
    tx = math.radians(target_theta_x_deg)
    ty = math.radians(target_theta_y_deg)
    sx = math.sin(tx)
    sy = math.sin(ty)
    sz2 = max(0.0, 1.0 - sx * sx - sy * sy)
    s = np.asarray([sx, sy, math.sqrt(sz2)])
    pts_m = np.column_stack([np.asarray(points_mm)[:, 0] * 1e-3, np.asarray(points_mm)[:, 1] * 1e-3, np.zeros(len(points_mm))])
    return np.exp(-1j * K * np.dot(pts_m, s))


def response_db(points_mm, steer_x_deg=0.0, steer_y_deg=0.0):
    axis, tx, ty, theta_deg, dirs, elem, mask = GRID
    pts_m = np.column_stack([np.asarray(points_mm)[:, 0] * 1e-3, np.asarray(points_mm)[:, 1] * 1e-3, np.zeros(len(points_mm))])
    w = steering_vector(points_mm, steer_x_deg, steer_y_deg)
    phase = np.exp(1j * K * np.tensordot(dirs, pts_m.T, axes=([2], [0])))
    resp = elem * np.dot(phase, w)
    target_amp = abs(np.dot(np.exp(1j * K * np.dot(pts_m, np.asarray([
        math.sin(math.radians(steer_x_deg)),
        math.sin(math.radians(steer_y_deg)),
        math.sqrt(max(0.0, 1.0 - math.sin(math.radians(steer_x_deg)) ** 2 - math.sin(math.radians(steer_y_deg)) ** 2)),
    ]))), w))
    db = 20.0 * np.log10(np.abs(resp) / max(target_amp, 1e-18) + 1e-12)
    db[~mask] = np.nan
    return db


def angular_separation_deg(tx, ty, steer_x_deg, steer_y_deg):
    sx0 = math.sin(math.radians(steer_x_deg))
    sy0 = math.sin(math.radians(steer_y_deg))
    sz0 = math.sqrt(max(0.0, 1.0 - sx0 * sx0 - sy0 * sy0))
    sx = np.sin(np.deg2rad(tx))
    sy = np.sin(np.deg2rad(ty))
    sz = np.sqrt(np.maximum(0.0, 1.0 - sx * sx - sy * sy))
    return np.rad2deg(np.arccos(np.clip(sx * sx0 + sy * sy0 + sz * sz0, -1.0, 1.0)))


def evaluate_layout(name, points_mm):
    pts = centered(points_mm)
    if len(pts) != 12:
        raise ValueError("layout must have 12 points")
    spacing = min_spacing_mm(pts)
    ap = aperture_mm(pts)
    axis, tx, ty, theta_deg, dirs, elem, mask = GRID
    cases = [(0, 0), (15, 0), (30, 0), (0, 30), (30, 30)]
    max_sidelobe = -1e9
    p99_sidelobe = -1e9
    case_rows = []
    for sx, sy in cases:
        db = response_db(pts, sx, sy)
        sep = angular_separation_deg(tx, ty, sx, sy)
        outside = (sep > 11.0) & mask
        sl = float(np.nanmax(db[outside]))
        p99 = float(np.nanpercentile(db[outside], 99.0))
        max_sidelobe = max(max_sidelobe, sl)
        p99_sidelobe = max(p99_sidelobe, p99)
        case_rows.append((sx, sy, sl, p99))
    unit_width = max(ap["width_centers_mm"], ap["height_centers_mm"]) + ELEMENT_DIAMETER_MM + 4.0
    score = max_sidelobe + 0.35 * p99_sidelobe
    if spacing < MIN_CENTER_SPACING_MM:
        score += 20.0 * (MIN_CENTER_SPACING_MM - spacing)
    return {
        "name": name,
        "points": pts,
        "score": score,
        "min_spacing_mm": spacing,
        "center_radius_mm": ap["center_radius_mm"],
        "width_centers_mm": ap["width_centers_mm"],
        "height_centers_mm": ap["height_centers_mm"],
        "recommended_unit_width_mm": unit_width,
        "worst_max_sidelobe_db": max_sidelobe,
        "worst_p99_sidelobe_db": p99_sidelobe,
        "case_rows": case_rows,
    }


def build_candidates():
    candidates = []
    for px in [18, 20, 22, 24]:
        for py in [18, 20, 22, 24]:
            pts = layout_rect_3x4(px, py)
            candidates.append(evaluate_layout("rect_3x4_%g_%g" % (px, py), pts))

    for r1 in np.arange(18.0, 26.1, 1.0):
        for r2 in np.arange(36.0, 47.1, 1.0):
            for rot4 in np.arange(0.0, 45.0, 15.0):
                for rot7 in np.arange(0.0, 360.0 / 7.0, 12.0):
                    pts = layout_ring_1_4_7(r1, r2, rot4, rot7)
                    if min_spacing_mm(pts) >= MIN_CENTER_SPACING_MM:
                        candidates.append(evaluate_layout("ring_1_4_7_r%.0f_%.0f" % (r1, r2), pts))

    for r1 in np.arange(21.0, 30.1, 1.0):
        for r2 in np.arange(37.0, 49.1, 1.0):
            for rot5 in np.arange(0.0, 72.0, 12.0):
                for rot7 in np.arange(0.0, 360.0 / 7.0, 12.0):
                    pts = layout_ring_5_7(r1, r2, rot5, rot7)
                    if min_spacing_mm(pts) >= MIN_CENTER_SPACING_MM:
                        candidates.append(evaluate_layout("ring_5_7_r%.0f_%.0f" % (r1, r2), pts))

    for radius in np.arange(39.0, 50.1, 1.0):
        for rot in np.arange(0.0, 2.0 * math.pi, math.pi / 10.0):
            pts = layout_sunflower(radius, rot)
            if min_spacing_mm(pts) >= MIN_CENTER_SPACING_MM:
                candidates.append(evaluate_layout("sunflower_r%.0f" % radius, pts))

    random.seed(20260511)
    for radius in [42.0, 44.0, 46.0, 48.0, 50.0]:
        kept = 0
        tries = 0
        while kept < 120 and tries < 4000:
            tries += 1
            pts = random_poisson_disk(radius, MIN_CENTER_SPACING_MM, attempts=1600)
            if pts is None:
                continue
            candidates.append(evaluate_layout("blue_noise_r%.0f_%03d" % (radius, kept), pts))
            kept += 1
    candidates.sort(key=lambda row: row["score"])
    return candidates


def write_layout_csv(path, points_mm):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Element_ID", "X_mm", "Y_mm"])
        for i, (x, y) in enumerate(points_mm):
            writer.writerow(["U%02d" % i, "%.3f" % x, "%.3f" % y])


def write_candidate_csv(path, candidates):
    fields = [
        "rank",
        "name",
        "score",
        "min_spacing_mm",
        "center_radius_mm",
        "width_centers_mm",
        "height_centers_mm",
        "recommended_unit_width_mm",
        "worst_max_sidelobe_db",
        "worst_p99_sidelobe_db",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for rank, row in enumerate(candidates[:200], start=1):
            out = {"rank": rank}
            for k in fields:
                if k == "rank":
                    continue
                out[k] = "%.6f" % row[k] if isinstance(row.get(k), float) else row.get(k)
            writer.writerow(out)


def plot_layout(points_mm, path):
    pts = np.asarray(points_mm)
    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    for i, (x, y) in enumerate(pts):
        circ = plt.Circle((x, y), ELEMENT_DIAMETER_MM * 0.5, fill=False, lw=1.8)
        ax.add_patch(circ)
        ax.text(x, y, str(i), ha="center", va="center", fontsize=8)
    half = math.ceil((max(np.max(np.abs(pts[:, 0])), np.max(np.abs(pts[:, 1]))) + ELEMENT_DIAMETER_MM * 0.5 + 3.0) / 5.0) * 5.0
    ax.plot([-half, half, half, -half, -half], [-half, -half, half, half, -half], "k--", lw=1.0, alpha=0.45)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-half - 8, half + 8)
    ax.set_ylim(-half - 8, half + 8)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("Recommended 12-element unit layout")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_patterns(points_mm, path):
    axis, tx, ty, theta_deg, dirs, elem, mask = GRID
    cases = [(0, 0), (15, 0), (30, 0), (0, 30)]
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 9.0))
    contour = None
    for ax, (sx, sy) in zip(axes.ravel(), cases):
        db = response_db(points_mm, sx, sy)
        db = np.maximum(db, DB_FLOOR)
        contour = ax.contourf(tx, ty, db, levels=np.linspace(DB_FLOOR, 0.0, 56), cmap="viridis", extend="min")
        ax.plot([sx], [sy], "rx", ms=8, mew=1.8)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(-70, 70)
        ax.set_ylim(-70, 70)
        ax.set_title("steer x=%g deg, y=%g deg" % (sx, sy))
        ax.set_xlabel("theta_x (deg)")
        ax.set_ylabel("theta_y (deg)")
        ax.grid(True, color="white", alpha=0.12)
    cbar = fig.colorbar(contour, ax=axes, shrink=0.82)
    cbar.set_label("normalized pressure (dB)")
    fig.suptitle("12-element unit circular far-field patterns")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary(path, candidates, chosen):
    with open(path, "w") as f:
        f.write("12-element ultrasonic unit layout search\n")
        f.write("Transducer: 16 mm diameter, 40 kHz, wavelength %.3f mm.\n" % (LAMBDA * 1000.0))
        f.write("Minimum center spacing used in search: %.1f mm.\n" % MIN_CENTER_SPACING_MM)
        f.write("Far-field score checks steering at (0,0), (15,0), (30,0), (0,30), (30,30) deg; sidelobe excludes 11 deg around target.\n\n")
        f.write("Top candidates:\n")
        for rank, row in enumerate(candidates[:12], start=1):
            f.write("  %2d %-22s score=%7.2f min_d=%5.1f mm unit=%5.1f mm worstSL=%6.2f dB p99=%6.2f dB\n" % (
                rank,
                row["name"],
                row["score"],
                row["min_spacing_mm"],
                row["recommended_unit_width_mm"],
                row["worst_max_sidelobe_db"],
                row["worst_p99_sidelobe_db"],
            ))
        f.write("\nRecommended layout: %s\n" % chosen["name"])
        f.write("Recommended square unit width: %.1f mm minimum; use 95-100 mm if mechanical margin is needed.\n" % chosen["recommended_unit_width_mm"])
        f.write("Center radius: %.1f mm, center-span: %.1f x %.1f mm, min center distance: %.1f mm.\n" % (
            chosen["center_radius_mm"],
            chosen["width_centers_mm"],
            chosen["height_centers_mm"],
            chosen["min_spacing_mm"],
        ))
        f.write("Worst max sidelobe across steering cases: %.2f dB; worst p99 sidelobe: %.2f dB.\n\n" % (
            chosen["worst_max_sidelobe_db"],
            chosen["worst_p99_sidelobe_db"],
        ))
        f.write("Coordinates in mm, unit centered at (0,0):\n")
        for i, (x, y) in enumerate(chosen["points"]):
            f.write("  U%02d  x=%8.3f  y=%8.3f\n" % (i, x, y))


def main():
    ensure_out_dir()
    candidates = build_candidates()
    manufacturable = [row for row in candidates if row["name"].startswith("ring_5_7")]
    chosen = manufacturable[0]
    layout_csv = os.path.join(OUT_DIR, "n12_unit_recommended_layout.csv")
    write_layout_csv(layout_csv, chosen["points"])
    cand_csv = os.path.join(OUT_DIR, "n12_unit_candidate_metrics.csv")
    write_candidate_csv(cand_csv, candidates)
    layout_png = os.path.join(OUT_DIR, "n12_unit_recommended_layout.png")
    plot_layout(chosen["points"], layout_png)
    patterns_png = os.path.join(OUT_DIR, "n12_unit_recommended_patterns.png")
    plot_patterns(chosen["points"], patterns_png)
    summary = os.path.join(OUT_DIR, "n12_unit_summary.txt")
    write_summary(summary, candidates, chosen)
    for path in [summary, layout_csv, cand_csv, layout_png, patterns_png]:
        print(path)


if __name__ == "__main__":
    main()
