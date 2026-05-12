"""
Agent parametric search for 12-element ultrasonic transducer array.
Extends design_n12_unit.py with finer resolution, local Nelder-Mead optimization,
additional evaluation cases, and extra blue-noise candidates.
"""

import csv
import math
import os
import random
import time

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

from design_n12_unit import (
    C,
    FREQ,
    LAMBDA,
    K,
    ELEMENT_DIAMETER_MM,
    ELEMENT_RADIUS_M,
    MIN_CENTER_SPACING_MM,
    DB_FLOOR,
    ensure_out_dir,
    piston_from_cos,
    ring_points,
    min_spacing_mm,
    aperture_mm,
    centered,
    layout_rect_3x4,
    layout_ring_1_4_7,
    layout_ring_5_7,
    layout_sunflower,
    random_poisson_disk,
    farfield_grid,
    GRID,
    steering_vector,
    response_db,
    angular_separation_deg,
    evaluate_layout,
    write_layout_csv,
    write_candidate_csv,
    plot_layout,
    plot_patterns,
)

OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "analysis_outputs")

# Extended evaluation cases per task requirement
CASES = [(0, 0), (15, 0), (30, 0), (0, 30), (30, 30), (0, 15), (15, 15)]


def evaluate_layout_extended(name, points_mm):
    """Same as evaluate_layout but with extended CASES list."""
    pts = centered(points_mm)
    if len(pts) != 12:
        raise ValueError("layout must have 12 points")
    spacing = min_spacing_mm(pts)
    ap = aperture_mm(pts)
    axis, tx, ty, theta_deg, dirs, elem, mask = GRID
    max_sidelobe = -1e9
    p99_sidelobe = -1e9
    case_rows = []
    for sx, sy in CASES:
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


# ---------------------------------------------------------------------------
# 1. Grid search over parametric families
# ---------------------------------------------------------------------------

def grid_search_candidates():
    candidates = []
    t0 = time.time()

    # Rect 3x4
    rect_count = 0
    for px in [18, 19, 20, 21, 22, 23, 24]:
        for py in [18, 19, 20, 21, 22, 23, 24]:
            pts = layout_rect_3x4(px, py)
            candidates.append(evaluate_layout_extended("rect_3x4_%g_%g" % (px, py), pts))
            rect_count += 1
    print("[grid] rect 3x4: %d candidates in %.1f s" % (rect_count, time.time() - t0))

    # Ring 1+4+7
    t1 = time.time()
    ring147_count = 0
    for r1 in np.arange(18.0, 26.01, 1.0):
        for r2 in np.arange(36.0, 47.01, 1.0):
            for rot4 in [0, 7.5, 15, 22.5, 30, 37.5]:
                # finer rot7: step 360/7/4 ~ 12.857/4 ~ 3.214, use 360/28 = 12.857? Wait:
                # 360/7 ~ 51.428. Finer step: 51.428/4 = 12.857. Let's step by 360/28 = 12.857 deg.
                # Actually we want finer than 12 deg. Use step = 360/7/4 = 12.857 deg.
                for rot7 in np.arange(0.0, 360.0 / 7.0, 360.0 / 7.0 / 4.0):
                    pts = layout_ring_1_4_7(r1, r2, rot4, rot7)
                    if min_spacing_mm(pts) >= MIN_CENTER_SPACING_MM:
                        candidates.append(
                            evaluate_layout_extended(
                                "ring_1_4_7_r%.1f_%.1f_rot4_%.1f_rot7_%.1f" % (r1, r2, rot4, rot7), pts
                            )
                        )
                        ring147_count += 1
    print("[grid] ring 1+4+7: %d candidates in %.1f s" % (ring147_count, time.time() - t1))

    # Ring 5+7
    t2 = time.time()
    ring57_count = 0
    for r1 in np.arange(21.0, 30.01, 1.0):
        for r2 in np.arange(37.0, 49.01, 1.0):
            for rot5 in np.arange(0.0, 72.0, 12.0):
                for rot7 in np.arange(0.0, 360.0 / 7.0, 360.0 / 7.0 / 4.0):
                    pts = layout_ring_5_7(r1, r2, rot5, rot7)
                    if min_spacing_mm(pts) >= MIN_CENTER_SPACING_MM:
                        candidates.append(
                            evaluate_layout_extended(
                                "ring_5_7_r%.1f_%.1f_rot5_%.1f_rot7_%.1f" % (r1, r2, rot5, rot7), pts
                            )
                        )
                        ring57_count += 1
    print("[grid] ring 5+7: %d candidates in %.1f s" % (ring57_count, time.time() - t2))

    # Sunflower
    t3 = time.time()
    sunflower_count = 0
    for radius in np.arange(39.0, 50.01, 1.0):
        # finer rotation: step pi/20 ~ 9 deg instead of 18 deg
        for rot in np.arange(0.0, 2.0 * math.pi, math.pi / 20.0):
            pts = layout_sunflower(radius, math.degrees(rot))
            if min_spacing_mm(pts) >= MIN_CENTER_SPACING_MM:
                candidates.append(
                    evaluate_layout_extended("sunflower_r%.1f_rot%.1f" % (radius, math.degrees(rot)), pts)
                )
                sunflower_count += 1
    print("[grid] sunflower: %d candidates in %.1f s" % (sunflower_count, time.time() - t3))

    return candidates


# ---------------------------------------------------------------------------
# 2. Local Nelder-Mead optimization on top candidates per family
# ---------------------------------------------------------------------------

def optimize_top(candidates, family_prefix, top_n=5):
    family = [c for c in candidates if c["name"].startswith(family_prefix)]
    family.sort(key=lambda x: x["score"])
    top = family[:top_n]
    optimized = []
    for cand in top:
        pts0 = cand["points"]
        # Choose parameterization based on family
        if family_prefix == "rect_3x4":
            # initial guess from name or from pts0
            # parse if possible, else estimate
            name = cand["name"]
            try:
                parts = name.split("_")
                px0 = float(parts[2])
                py0 = float(parts[3])
            except Exception:
                # estimate from span
                xs = pts0[:, 0]
                ys = pts0[:, 1]
                px0 = (np.max(xs) - np.min(xs)) / 3.0
                py0 = (np.max(ys) - np.min(ys)) / 2.0

            def obj_rect(params):
                px, py = params
                if px < 15 or py < 15 or px > 30 or py > 30:
                    return 1e6
                pts = layout_rect_3x4(px, py)
                if min_spacing_mm(pts) < MIN_CENTER_SPACING_MM:
                    return 1e6
                return evaluate_layout_extended("tmp", pts)["score"]

            res = minimize(obj_rect, [px0, py0], method="Nelder-Mead", options={"maxiter": 200, "xatol": 0.01, "fatol": 0.01})
            if res.success and res.fun < cand["score"] - 1e-6:
                pts = layout_rect_3x4(res.x[0], res.x[1])
                opt = evaluate_layout_extended(
                    "rect_3x4_%.3f_%.3f" % (res.x[0], res.x[1]), pts
                )
                optimized.append(opt)
            else:
                optimized.append(cand)

        elif family_prefix == "ring_1_4_7":
            name = cand["name"]
            try:
                parts = name.split("_")
                # ring_1_4_7_r18.0_36.0_rot4_0.0_rot7_0.0
                r1_0 = float(parts[4])
                r2_0 = float(parts[5])
                rot4_0 = float(parts[7])
                rot7_0 = float(parts[9])
            except Exception:
                r1_0 = 22.0
                r2_0 = 41.0
                rot4_0 = 0.0
                rot7_0 = 0.0

            def obj_ring147(params):
                r1, r2, rot4, rot7 = params
                if r1 < 10 or r2 < 20 or r1 > 35 or r2 > 55 or abs(rot4) > 90 or abs(rot7) > 360:
                    return 1e6
                pts = layout_ring_1_4_7(r1, r2, rot4, rot7)
                if min_spacing_mm(pts) < MIN_CENTER_SPACING_MM:
                    return 1e6
                return evaluate_layout_extended("tmp", pts)["score"]

            res = minimize(
                obj_ring147,
                [r1_0, r2_0, rot4_0, rot7_0],
                method="Nelder-Mead",
                options={"maxiter": 200, "xatol": 0.01, "fatol": 0.01},
            )
            if res.success and res.fun < cand["score"] - 1e-6:
                pts = layout_ring_1_4_7(res.x[0], res.x[1], res.x[2], res.x[3])
                opt = evaluate_layout_extended(
                    "ring_1_4_7_r%.3f_%.3f_rot4_%.3f_rot7_%.3f" % tuple(res.x), pts
                )
                optimized.append(opt)
            else:
                optimized.append(cand)

        elif family_prefix == "ring_5_7":
            name = cand["name"]
            try:
                parts = name.split("_")
                r1_0 = float(parts[3])
                r2_0 = float(parts[4])
                rot5_0 = float(parts[6])
                rot7_0 = float(parts[8])
            except Exception:
                r1_0 = 25.0
                r2_0 = 43.0
                rot5_0 = 0.0
                rot7_0 = 0.0

            def obj_ring57(params):
                r1, r2, rot5, rot7 = params
                if r1 < 10 or r2 < 20 or r1 > 35 or r2 > 55 or abs(rot5) > 90 or abs(rot7) > 360:
                    return 1e6
                pts = layout_ring_5_7(r1, r2, rot5, rot7)
                if min_spacing_mm(pts) < MIN_CENTER_SPACING_MM:
                    return 1e6
                return evaluate_layout_extended("tmp", pts)["score"]

            res = minimize(
                obj_ring57,
                [r1_0, r2_0, rot5_0, rot7_0],
                method="Nelder-Mead",
                options={"maxiter": 200, "xatol": 0.01, "fatol": 0.01},
            )
            if res.success and res.fun < cand["score"] - 1e-6:
                pts = layout_ring_5_7(res.x[0], res.x[1], res.x[2], res.x[3])
                opt = evaluate_layout_extended(
                    "ring_5_7_r%.3f_%.3f_rot5_%.3f_rot7_%.3f" % tuple(res.x), pts
                )
                optimized.append(opt)
            else:
                optimized.append(cand)

        elif family_prefix == "sunflower":
            name = cand["name"]
            try:
                parts = name.split("_")
                rad0 = float(parts[1][1:])
                rot0 = float(parts[2][3:])
            except Exception:
                rad0 = 45.0
                rot0 = 0.0

            def obj_sun(params):
                rad, rot = params
                if rad < 20 or rad > 60 or abs(rot) > 360:
                    return 1e6
                pts = layout_sunflower(rad, rot)
                if min_spacing_mm(pts) < MIN_CENTER_SPACING_MM:
                    return 1e6
                return evaluate_layout_extended("tmp", pts)["score"]

            res = minimize(
                obj_sun,
                [rad0, rot0],
                method="Nelder-Mead",
                options={"maxiter": 200, "xatol": 0.01, "fatol": 0.01},
            )
            if res.success and res.fun < cand["score"] - 1e-6:
                pts = layout_sunflower(res.x[0], res.x[1])
                opt = evaluate_layout_extended(
                    "sunflower_r%.3f_rot%.3f" % (res.x[0], res.x[1]), pts
                )
                optimized.append(opt)
            else:
                optimized.append(cand)
        else:
            optimized.append(cand)
    return optimized


# ---------------------------------------------------------------------------
# 3. Blue-noise / Poisson disk layouts
# ---------------------------------------------------------------------------

def generate_blue_noise(n=200, seed=20260512):
    random.seed(seed)
    candidates = []
    t0 = time.time()
    kept = 0
    tries = 0
    while kept < n and tries < 10000:
        tries += 1
        radius = 50.0
        pts = random_poisson_disk(radius, MIN_CENTER_SPACING_MM, attempts=2000)
        if pts is None:
            continue
        candidates.append(evaluate_layout_extended("blue_noise_%03d" % kept, pts))
        kept += 1
    print("[blue-noise] generated %d candidates in %.1f s (tries=%d)" % (kept, time.time() - t0, tries))
    return candidates


# ---------------------------------------------------------------------------
# 4. Plotting helpers for extended cases
# ---------------------------------------------------------------------------

def plot_patterns_extended(points_mm, path):
    axis, tx, ty, theta_deg, dirs, elem, mask = GRID
    cases = CASES  # 7 cases
    # Use 3 rows x 3 cols, leaving one empty
    fig, axes = plt.subplots(3, 3, figsize=(14.0, 13.0))
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
    # hide unused subplots
    for ax in axes.ravel()[len(cases):]:
        ax.axis("off")
    cbar = fig.colorbar(contour, ax=axes, shrink=0.82)
    cbar.set_label("normalized pressure (dB)")
    fig.suptitle("12-element unit circular far-field patterns (extended cases)")
    fig.savefig(path, dpi=180)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 5. Main
# ---------------------------------------------------------------------------

def main():
    ensure_out_dir()
    t_start = time.time()

    # Grid search
    candidates = grid_search_candidates()

    # Local optimization per family (top 5)
    families = ["rect_3x4", "ring_1_4_7", "ring_5_7", "sunflower"]
    for fam in families:
        t0 = time.time()
        opt = optimize_top(candidates, fam, top_n=5)
        # Replace old family entries with optimized ones
        candidates = [c for c in candidates if not c["name"].startswith(fam)]
        candidates.extend(opt)
        print("[optimize] %s: %.1f s" % (fam, time.time() - t0))

    # Blue noise
    blue = generate_blue_noise(n=200)
    candidates.extend(blue)

    # Sort globally
    candidates.sort(key=lambda row: row["score"])

    best = candidates[0]

    # Write outputs
    cand_csv = os.path.join(OUT_DIR, "agent_parametric_candidates.csv")
    write_candidate_csv(cand_csv, candidates)

    best_csv = os.path.join(OUT_DIR, "agent_parametric_best.csv")
    write_layout_csv(best_csv, best["points"])

    layout_png = os.path.join(OUT_DIR, "agent_parametric_layout.png")
    plot_layout(best["points"], layout_png)

    patterns_png = os.path.join(OUT_DIR, "agent_parametric_patterns.png")
    plot_patterns_extended(best["points"], patterns_png)

    print("\n===== SUMMARY =====")
    print("Best score: %.4f" % best["score"])
    print("Best family/name: %s" % best["name"])
    print("Min spacing: %.3f mm" % best["min_spacing_mm"])
    print("Recommended unit width: %.2f mm" % best["recommended_unit_width_mm"])
    print("Worst max sidelobe: %.2f dB" % best["worst_max_sidelobe_db"])
    print("Worst p99 sidelobe: %.2f dB" % best["worst_p99_sidelobe_db"])
    print("Total runtime: %.1f s" % (time.time() - t_start))
    print("===================\n")

    for path in [cand_csv, best_csv, layout_png, patterns_png]:
        print(path)


if __name__ == "__main__":
    main()
