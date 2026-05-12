import csv
import math
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

import design_n12_unit as base


ROOT = os.path.dirname(os.path.dirname(__file__))
OUT_DIR = os.path.join(ROOT, "analysis_outputs")


def ensure_out_dir():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)


def layout_1_5_7(r1_mm, r2_mm, rot5_deg, rot7_deg):
    return base.centered([[0.0, 0.0]] + base.ring_points(r1_mm, 5, rot5_deg) + base.ring_points(r2_mm, 7, rot7_deg))


def evaluate_layout(name, points_mm):
    pts = base.centered(points_mm)
    spacing = base.min_spacing_mm(pts)
    ap = base.aperture_mm(pts)
    axis, tx, ty, theta_deg, dirs, elem, mask = base.GRID
    cases = [(0, 0), (15, 0), (30, 0), (0, 30), (30, 30)]
    max_sidelobe = -1e9
    p99_sidelobe = -1e9
    case_rows = []
    for sx, sy in cases:
        db = base.response_db(pts, sx, sy)
        sep = base.angular_separation_deg(tx, ty, sx, sy)
        outside = (sep > 11.0) & mask
        sl = float(np.nanmax(db[outside]))
        p99 = float(np.nanpercentile(db[outside], 99.0))
        max_sidelobe = max(max_sidelobe, sl)
        p99_sidelobe = max(p99_sidelobe, p99)
        case_rows.append({"steer_x_deg": sx, "steer_y_deg": sy, "max_sidelobe_db": sl, "p99_sidelobe_db": p99})
    unit_width = max(ap["width_centers_mm"], ap["height_centers_mm"]) + base.ELEMENT_DIAMETER_MM + 4.0
    score = max_sidelobe + 0.35 * p99_sidelobe
    if spacing < base.MIN_CENTER_SPACING_MM:
        score += 20.0 * (base.MIN_CENTER_SPACING_MM - spacing)
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


def search_1_5_7():
    candidates = []
    for r1 in np.arange(21.0, 29.1, 1.0):
        for r2 in np.arange(38.0, 48.1, 1.0):
            for rot5 in np.arange(0.0, 72.0, 12.0):
                for rot7 in np.arange(0.0, 360.0 / 7.0, 12.0):
                    pts = layout_1_5_7(r1, r2, rot5, rot7)
                    if base.min_spacing_mm(pts) >= base.MIN_CENTER_SPACING_MM:
                        candidates.append(evaluate_layout("ring_1_5_7_r%.0f_%.0f" % (r1, r2), pts))
    candidates.sort(key=lambda row: row["score"])
    return candidates


def write_layout_csv(path, points_mm):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Element_ID", "X_mm", "Y_mm"])
        for i, (x, y) in enumerate(points_mm):
            writer.writerow(["U%02d" % i, "%.3f" % x, "%.3f" % y])


def write_case_csv(path, rows):
    fields = ["layout", "steer_x_deg", "steer_y_deg", "max_sidelobe_db", "p99_sidelobe_db"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_layout(points_mm, path):
    pts = np.asarray(points_mm)
    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    for i, (x, y) in enumerate(pts):
        circ = plt.Circle((x, y), base.ELEMENT_DIAMETER_MM * 0.5, fill=False, lw=1.8)
        ax.add_patch(circ)
        ax.text(x, y, str(i), ha="center", va="center", fontsize=8)
    half = math.ceil((max(np.max(np.abs(pts[:, 0])), np.max(np.abs(pts[:, 1]))) + base.ELEMENT_DIAMETER_MM * 0.5 + 3.0) / 5.0) * 5.0
    ax.plot([-half, half, half, -half, -half], [-half, -half, half, half, -half], "k--", lw=1.0, alpha=0.45)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-half - 8, half + 8)
    ax.set_ylim(-half - 8, half + 8)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("Recommended 13-element unit layout")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_patterns(points_mm, path):
    axis, tx, ty, theta_deg, dirs, elem, mask = base.GRID
    cases = [(0, 0), (15, 0), (30, 0), (0, 30)]
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 9.0))
    contour = None
    for ax, (sx, sy) in zip(axes.ravel(), cases):
        db = base.response_db(points_mm, sx, sy)
        db = np.maximum(db, base.DB_FLOOR)
        contour = ax.contourf(tx, ty, db, levels=np.linspace(base.DB_FLOOR, 0.0, 56), cmap="viridis", extend="min")
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
    fig.suptitle("13-element unit circular far-field patterns")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary(path, direct, optimized, top):
    gain_db = 20.0 * math.log10(13.0 / 12.0)
    with open(path, "w") as f:
        f.write("13-element ultrasonic unit report\n")
        f.write("Direct version: add one center transducer to the previous 5+7 ring unit.\n")
        f.write("Theoretical coherent on-target pressure gain vs 12 elements: %.2f dB.\n\n" % gain_db)
        f.write("12-element reference from previous report: worstSL=-6.16 dB, p99=-9.28 dB, unit ~=100 mm.\n\n")
        for label, row in [("direct_center_added", direct), ("optimized_1_5_7", optimized)]:
            f.write("%s:\n" % label)
            f.write("  layout=%s\n" % row["name"])
            f.write("  unit_width=%.1f mm, center_radius=%.1f mm, min_spacing=%.1f mm\n" % (
                row["recommended_unit_width_mm"], row["center_radius_mm"], row["min_spacing_mm"]
            ))
            f.write("  worst max sidelobe=%.2f dB, worst p99 sidelobe=%.2f dB, score=%.2f\n" % (
                row["worst_max_sidelobe_db"], row["worst_p99_sidelobe_db"], row["score"]
            ))
            f.write("  coordinates:\n")
            for i, (x, y) in enumerate(row["points"]):
                f.write("    U%02d x=%8.3f y=%8.3f\n" % (i, x, y))
            f.write("\n")
        f.write("Top optimized 1+5+7 candidates:\n")
        for rank, row in enumerate(top[:10], start=1):
            f.write("  %2d %-20s score=%7.2f min_d=%5.1f unit=%5.1f worstSL=%6.2f p99=%6.2f\n" % (
                rank, row["name"], row["score"], row["min_spacing_mm"], row["recommended_unit_width_mm"],
                row["worst_max_sidelobe_db"], row["worst_p99_sidelobe_db"]
            ))


def main():
    ensure_out_dir()
    direct_points = layout_1_5_7(22.0, 41.0, 60.0, 0.0)
    direct = evaluate_layout("direct_1_5_7_r22_41", direct_points)
    candidates = search_1_5_7()
    optimized = candidates[0]
    chosen = direct

    layout_csv = os.path.join(OUT_DIR, "n13_unit_center_added_layout.csv")
    write_layout_csv(layout_csv, chosen["points"])
    case_rows = []
    for label, row in [("direct_center_added", direct), ("optimized_1_5_7", optimized)]:
        for case in row["case_rows"]:
            case_rows.append({"layout": label, **case})
    case_csv = os.path.join(OUT_DIR, "n13_unit_case_metrics.csv")
    write_case_csv(case_csv, case_rows)
    layout_png = os.path.join(OUT_DIR, "n13_unit_center_added_layout.png")
    plot_layout(chosen["points"], layout_png)
    patterns_png = os.path.join(OUT_DIR, "n13_unit_center_added_patterns.png")
    plot_patterns(chosen["points"], patterns_png)
    summary = os.path.join(OUT_DIR, "n13_unit_summary.txt")
    write_summary(summary, direct, optimized, candidates)
    for path in [summary, layout_csv, case_csv, layout_png, patterns_png]:
        print(path)


if __name__ == "__main__":
    main()
