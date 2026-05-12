import csv
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import design_n36_array as base
import software_weighting_compare as sw


ROOT = os.path.dirname(os.path.dirname(__file__))
OUT_DIR = os.path.join(ROOT, "analysis_outputs")
N38_CSV = os.path.join(OUT_DIR, "n38_refined_layout.csv")

FOCUS_Z_LIST_MM = [100, 150, 200, 300, 400, 500]
OFFSET_LIST_MM = [0, 10, 20, 30, 40, 50, 60, 70, 80]
GRID_RANGE_MM = 140.0
GRID_STEP_MM = 2.0
MAIN_EXCLUDE_MM = 20.0


def ensure_out_dir():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)


def load_xy_csv(path):
    pts = []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pts.append([float(row["X_mm"]) * 1e-3, float(row["Y_mm"]) * 1e-3])
    return np.asarray(pts, dtype=float)


def line_width(axis_mm, line_db, center_idx, threshold_db):
    left = center_idx
    while left > 0 and line_db[left] >= threshold_db:
        left -= 1
    right = center_idx
    while right < len(line_db) - 1 and line_db[right] >= threshold_db:
        right += 1
    if left == 0 or right == len(line_db) - 1:
        return float("nan")

    def interp(i0, i1):
        if abs(line_db[i1] - line_db[i0]) < 1e-12:
            return axis_mm[i0]
        return axis_mm[i0] + (threshold_db - line_db[i0]) * (axis_mm[i1] - axis_mm[i0]) / (line_db[i1] - line_db[i0])

    return float(interp(right - 1, right) - interp(left, left + 1))


def phase_only_weights(positions, target):
    return sw.phase_only_weights(positions, target)


def rls_amp_weights(positions, target):
    return sw.rls_amplitude_weights(positions, target, target_fraction=0.72, mu=1e-4)


def focus_metrics(positions, target, weights, phase_target_amp):
    x_mm = np.arange(-GRID_RANGE_MM, GRID_RANGE_MM + 0.001, GRID_STEP_MM)
    y_mm = np.arange(-GRID_RANGE_MM, GRID_RANGE_MM + 0.001, GRID_STEP_MM)
    xg, yg = np.meshgrid(x_mm * 1e-3, y_mm * 1e-3)
    zg = np.zeros_like(xg) + target["z_m"]
    h = sw.nearfield_vector(positions, xg, yg, zg)
    resp = np.dot(h, weights).reshape(xg.shape)
    vt = sw.target_vector(positions, target)
    target_amp = float(abs(np.dot(vt, weights)))
    pressure_db = 20.0 * np.log10(np.abs(resp) / max(target_amp, 1e-18) + 1e-12)

    tx_mm = target["x_m"] * 1000.0
    ty_mm = target["y_m"] * 1000.0
    rr = np.sqrt((x_mm[None, :] - tx_mm) ** 2 + (y_mm[:, None] - ty_mm) ** 2)
    outside = rr > MAIN_EXCLUDE_MM
    psl_pressure_db = float(np.max(pressure_db[outside]))
    p995_pressure_db = float(np.percentile(pressure_db[outside], 99.5))

    peak_idx = np.unravel_index(np.argmax(np.abs(resp)), resp.shape)
    peak_x_mm = float(x_mm[peak_idx[1]])
    peak_y_mm = float(y_mm[peak_idx[0]])
    focus_error_mm = math.hypot(peak_x_mm - tx_mm, peak_y_mm - ty_mm)

    target_ix = int(np.argmin(np.abs(x_mm - tx_mm)))
    target_iy = int(np.argmin(np.abs(y_mm - ty_mm)))
    width_pressure_x_mm = line_width(x_mm, pressure_db[target_iy, :], target_ix, -3.0)
    width_pressure_y_mm = line_width(y_mm, pressure_db[:, target_ix], target_iy, -3.0)
    width_audio_x_mm = line_width(x_mm, pressure_db[target_iy, :], target_ix, -1.5)
    width_audio_y_mm = line_width(y_mm, pressure_db[:, target_ix], target_iy, -1.5)

    gain_drop_db = 20.0 * math.log10(target_amp / max(phase_target_amp, 1e-18) + 1e-12)
    angle_deg = math.degrees(math.atan2(target["x_m"], target["z_m"]))
    return {
        "z_mm": target["z_m"] * 1000.0,
        "x_mm": tx_mm,
        "angle_deg": angle_deg,
        "pressure_psl_db": psl_pressure_db,
        "pressure_p99_5_db": p995_pressure_db,
        "audio_source_psl_proxy_db": 2.0 * psl_pressure_db,
        "audio_source_p99_5_proxy_db": 2.0 * p995_pressure_db,
        "focus_error_mm": focus_error_mm,
        "pressure_fwhm_x_mm": width_pressure_x_mm,
        "pressure_fwhm_y_mm": width_pressure_y_mm,
        "audio_source_fwhm_x_mm": width_audio_x_mm,
        "audio_source_fwhm_y_mm": width_audio_y_mm,
        "target_gain_drop_db": gain_drop_db,
        "target_amp": target_amp,
        "peak_x_mm": peak_x_mm,
        "peak_y_mm": peak_y_mm,
    }


def target_case(x_mm, z_mm):
    return {
        "name": "x%d_z%d" % (x_mm, z_mm),
        "domain": "near",
        "x_m": x_mm * 1e-3,
        "y_m": 0.0,
        "z_m": z_mm * 1e-3,
    }


def scan_layout(layout_name, positions):
    rows = []
    for z_mm in FOCUS_Z_LIST_MM:
        for x_mm in OFFSET_LIST_MM:
            target = target_case(x_mm, z_mm)
            phase_w = phase_only_weights(positions, target)
            phase_target_amp = float(abs(np.dot(sw.target_vector(positions, target), phase_w)))
            for method, fn in [("phase_only", phase_only_weights), ("rls_amp", rls_amp_weights)]:
                weights = fn(positions, target)
                metrics = focus_metrics(positions, target, weights, phase_target_amp)
                metrics["layout"] = layout_name
                metrics["method"] = method
                metrics["element_count"] = len(positions)
                metrics["aperture_mm"] = 2.0 * base.aperture_radius(positions) * 1000.0
                rows.append(metrics)
    return rows


def usable_radius_by_z(rows, pressure_psl_limit=-10.0, gain_limit=-5.0, err_limit=5.0):
    out = []
    for layout in sorted(set(r["layout"] for r in rows)):
        for method in sorted(set(r["method"] for r in rows)):
            for z_mm in FOCUS_Z_LIST_MM:
                subset = [
                    r for r in rows
                    if r["layout"] == layout and r["method"] == method and abs(r["z_mm"] - z_mm) < 1e-6
                ]
                usable = [
                    r for r in subset
                    if r["pressure_psl_db"] <= pressure_psl_limit
                    and r["target_gain_drop_db"] >= gain_limit
                    and r["focus_error_mm"] <= err_limit
                ]
                radius = max([r["x_mm"] for r in usable], default=float("nan"))
                angle = math.degrees(math.atan2(radius, z_mm)) if not math.isnan(radius) else float("nan")
                out.append({
                    "layout": layout,
                    "method": method,
                    "z_mm": z_mm,
                    "usable_radius_mm": radius,
                    "usable_angle_deg": angle,
                })
    return out


def write_csv(path, rows):
    fields = [
        "layout", "method", "element_count", "aperture_mm", "z_mm", "x_mm", "angle_deg",
        "pressure_psl_db", "pressure_p99_5_db", "audio_source_psl_proxy_db", "audio_source_p99_5_proxy_db",
        "focus_error_mm", "pressure_fwhm_x_mm", "pressure_fwhm_y_mm",
        "audio_source_fwhm_x_mm", "audio_source_fwhm_y_mm",
        "target_gain_drop_db", "target_amp", "peak_x_mm", "peak_y_mm",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_usable_radius(usable_rows):
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    for layout in sorted(set(r["layout"] for r in usable_rows)):
        for method in ["phase_only", "rls_amp"]:
            subset = [r for r in usable_rows if r["layout"] == layout and r["method"] == method]
            ax.plot(
                [r["z_mm"] for r in subset],
                [r["usable_radius_mm"] for r in subset],
                marker="o",
                lw=1.8,
                label="%s %s" % (layout, method),
            )
    ax.set_xlabel("focus distance z (mm)")
    ax.set_ylabel("max usable lateral offset (mm)")
    ax.set_title("Usable focus range, pressure PSL <= -10 dB, gain drop >= -5 dB")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    path = os.path.join(OUT_DIR, "audio_focus_usable_radius.png")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def write_summary(path, rows, usable_rows):
    def find(layout, method, z, x):
        for r in rows:
            if r["layout"] == layout and r["method"] == method and abs(r["z_mm"] - z) < 1e-6 and abs(r["x_mm"] - x) < 1e-6:
                return r
        return None

    with open(path, "w") as f:
        f.write("Audio transmission focus performance, carrier model: 40 kHz 16 mm circular piston\n")
        f.write("Audio source proxy assumes local demodulated source strength proportional to |P_ultra|^2.\n")
        f.write("Therefore audio_source_psl_proxy_db = 2 * pressure_psl_db; audible propagation is not included.\n\n")
        for layout in ["N32_ring", "N38_blue_noise"]:
            f.write("[%s]\n" % layout)
            for method in ["phase_only", "rls_amp"]:
                f.write("  %s representative points\n" % method)
                for z, x in [(150, 0), (150, 30), (200, 0), (200, 40), (300, 0), (300, 60), (500, 0), (500, 80)]:
                    r = find(layout, method, z, x)
                    if r is None:
                        continue
                    f.write(
                        "    z=%3.0f x=%2.0f angle=%5.1f deg: pressure_PSL=%6.2f dB, audio_proxy_PSL=%6.2f dB, "
                        "P_FWHM=(%5.1f,%5.1f)mm, audio_proxy_FWHM=(%5.1f,%5.1f)mm, gain=%6.2f dB, err=%4.1fmm\n"
                        % (
                            r["z_mm"], r["x_mm"], r["angle_deg"], r["pressure_psl_db"],
                            r["audio_source_psl_proxy_db"], r["pressure_fwhm_x_mm"], r["pressure_fwhm_y_mm"],
                            r["audio_source_fwhm_x_mm"], r["audio_source_fwhm_y_mm"],
                            r["target_gain_drop_db"], r["focus_error_mm"],
                        )
                    )
            f.write("\n")
        f.write("Usable radius by z with pressure PSL <= -10 dB, gain drop >= -5 dB, focus error <= 5 mm:\n")
        for r in usable_rows:
            f.write(
                "  %-15s %-10s z=%3.0f mm: radius=%5.1f mm, angle=%5.1f deg\n"
                % (r["layout"], r["method"], r["z_mm"], r["usable_radius_mm"], r["usable_angle_deg"])
            )


def main():
    ensure_out_dir()
    layouts = [
        ("N32_ring", base.load_current_layout()),
        ("N38_blue_noise", load_xy_csv(N38_CSV)),
    ]
    all_rows = []
    for layout_name, positions in layouts:
        all_rows.extend(scan_layout(layout_name, positions))

    csv_path = os.path.join(OUT_DIR, "audio_focus_scan.csv")
    write_csv(csv_path, all_rows)
    usable = usable_radius_by_z(all_rows)
    usable_csv = os.path.join(OUT_DIR, "audio_focus_usable_radius.csv")
    with open(usable_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["layout", "method", "z_mm", "usable_radius_mm", "usable_angle_deg"])
        writer.writeheader()
        for row in usable:
            writer.writerow(row)
    plot_path = plot_usable_radius(usable)
    summary_path = os.path.join(OUT_DIR, "audio_focus_summary.txt")
    write_summary(summary_path, all_rows, usable)
    print(csv_path)
    print(usable_csv)
    print(summary_path)
    print(plot_path)


if __name__ == "__main__":
    main()
