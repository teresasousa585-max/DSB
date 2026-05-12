import csv
import math
import os

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
ELEMENT_RADIUS = 8e-3
SEGMENT_WIDTH = 100e-3
RING_RADII = [18e-3, 38e-3]
RING_COUNTS = [4, 10]

INNER_TILT_DEG = 45.0
OUTER_TILT_DEG = 65.0
EAR_SPACING_M = 0.18
LISTENER_Z_M = [0.5, 0.8, 1.0]
REGION_HALF_MM = 35.0
REGION_STEP_MM = 10.0


def ensure_out_dir():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)


def piston_from_cos(cos_theta):
    cos_theta = np.asarray(cos_theta)
    front = cos_theta > 0.0
    sin_theta = np.sqrt(np.maximum(0.0, 1.0 - cos_theta * cos_theta))
    x = K * ELEMENT_RADIUS * sin_theta
    out = np.zeros_like(cos_theta, dtype=float)
    small = np.abs(x) < 1e-10
    out[front & small] = 1.0
    mask = front & (~small)
    out[mask] = np.abs(2.0 * j1(x[mask]) / x[mask])
    return out


def local_15_positions():
    pts = [[0.0, 0.0]]
    for r, n in zip(RING_RADII, RING_COUNTS):
        for i in range(n):
            a = 2.0 * np.pi * i / n
            pts.append([r * math.cos(a), r * math.sin(a)])
    return np.asarray(pts, dtype=float)


def connected_array(a1_deg=INNER_TILT_DEG, a2_deg=OUTER_TILT_DEG):
    tilts = [float(a2_deg), float(a1_deg), 0.0, -float(a1_deg), -float(a2_deg)]
    local = local_15_positions()
    uxs = []
    normals_by_seg = []
    for tilt_deg in tilts:
        alpha = math.radians(tilt_deg)
        normals_by_seg.append(np.asarray([math.sin(alpha), 0.0, math.cos(alpha)]))
        uxs.append(np.asarray([math.cos(alpha), 0.0, -math.sin(alpha)]))

    centers = [None] * 5
    centers[2] = np.asarray([0.0, 0.0, 0.0])
    for seg in [1, 0]:
        hinge = centers[seg + 1] - 0.5 * SEGMENT_WIDTH * uxs[seg + 1]
        centers[seg] = hinge - 0.5 * SEGMENT_WIDTH * uxs[seg]
    for seg in [3, 4]:
        hinge = centers[seg - 1] + 0.5 * SEGMENT_WIDTH * uxs[seg - 1]
        centers[seg] = hinge + 0.5 * SEGMENT_WIDTH * uxs[seg]

    positions = []
    normals = []
    segment_ids = []
    for seg in range(5):
        uy = np.asarray([0.0, 1.0, 0.0])
        for xl, yl in local:
            positions.append(centers[seg] + xl * uxs[seg] + yl * uy)
            normals.append(normals_by_seg[seg])
            segment_ids.append(seg)
    return np.asarray(positions), np.asarray(normals), np.asarray(segment_ids), tilts, np.asarray(centers)


def nearfield_matrix(positions, normals, x, y, z):
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    pts = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
    diff = pts[:, None, :] - positions[None, :, :]
    r = np.linalg.norm(diff, axis=2)
    dirs = diff / np.maximum(r[:, :, None], 1e-12)
    cosang = np.sum(dirs * normals[None, :, :], axis=2)
    elem = piston_from_cos(cosang)
    return elem * np.exp(1j * K * r) / np.maximum(r, 1e-12)


def h_at(positions, normals, point):
    return nearfield_matrix(
        positions,
        normals,
        np.asarray([point[0]]),
        np.asarray([point[1]]),
        np.asarray([point[2]]),
    )[0]


def phase_weights_for_point(positions, normals, point):
    h = h_at(positions, normals, point)
    return np.exp(-1j * np.angle(h))


def wrong_region_points(wrong_point):
    xs = np.arange(wrong_point[0] * 1000.0 - REGION_HALF_MM, wrong_point[0] * 1000.0 + REGION_HALF_MM + 0.001, REGION_STEP_MM)
    ys = np.arange(wrong_point[1] * 1000.0 - REGION_HALF_MM, wrong_point[1] * 1000.0 + REGION_HALF_MM + 0.001, REGION_STEP_MM)
    xg, yg = np.meshgrid(xs * 1e-3, ys * 1e-3)
    zg = np.zeros_like(xg) + wrong_point[2]
    return xg.ravel(), yg.ravel(), zg.ravel()


def suppress_wrong_region_weights(positions, normals, target_point, wrong_point, reg=5e-4):
    ht = h_at(positions, normals, target_point)
    x, y, z = wrong_region_points(wrong_point)
    hs = nearfield_matrix(positions, normals, x, y, z)
    rmat = np.real(hs.conj().T @ hs) / len(hs)
    rmat = rmat + reg * np.eye(len(positions))
    v = np.linalg.solve(rmat, ht.conj())
    denom = ht @ v
    if abs(denom) < 1e-18:
        w = phase_weights_for_point(positions, normals, target_point)
    else:
        w = v / denom
    max_amp = np.max(np.abs(w))
    if max_amp > 1e-18:
        w = w / max_amp
    return w


def line_width(axis_mm, line_db, idx, threshold_db=-3.0):
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


def focus_width_x(positions, normals, target_point, weights, span_mm=240.0, step_mm=2.0):
    x0, y0, z0 = target_point
    xs = np.arange(x0 * 1000.0 - span_mm, x0 * 1000.0 + span_mm + 0.001, step_mm)
    ys = np.zeros_like(xs) + y0
    zs = np.zeros_like(xs) + z0
    h = nearfield_matrix(positions, normals, xs * 1e-3, ys, zs)
    resp = h @ weights
    target_amp = abs(h_at(positions, normals, target_point) @ weights)
    db = 20.0 * np.log10(np.abs(resp) / max(target_amp, 1e-18) + 1e-12)
    ix = int(np.argmin(np.abs(xs - x0 * 1000.0)))
    return line_width(xs, db, ix), xs, np.maximum(db, -60.0)


def evaluate_channel(positions, normals, target_point, wrong_point, weights):
    target_resp = abs(h_at(positions, normals, target_point) @ weights)
    wrong_center_resp = abs(h_at(positions, normals, wrong_point) @ weights)
    x, y, z = wrong_region_points(wrong_point)
    wrong_region_resp = np.abs(nearfield_matrix(positions, normals, x, y, z) @ weights)
    fwhm_x, _, _ = focus_width_x(positions, normals, target_point, weights)
    center_xtalk_db = 20.0 * math.log10(wrong_center_resp / max(target_resp, 1e-18) + 1e-12)
    region_xtalk_db = 20.0 * math.log10(float(np.max(wrong_region_resp)) / max(target_resp, 1e-18) + 1e-12)
    return {
        "target_pressure_amp": float(target_resp),
        "point_crosstalk_pressure_db": center_xtalk_db,
        "region_crosstalk_pressure_db": region_xtalk_db,
        "point_crosstalk_audio_proxy_db": 2.0 * center_xtalk_db,
        "region_crosstalk_audio_proxy_db": 2.0 * region_xtalk_db,
        "fwhm_x_mm": fwhm_x,
        "max_weight_amp": float(np.max(np.abs(weights))),
        "mean_weight_amp": float(np.mean(np.abs(weights))),
    }


def run_scan():
    positions, normals, segment_ids, tilts, centers = connected_array()
    rows = []
    line_plots = []
    for z_m in LISTENER_Z_M:
        left = np.asarray([-EAR_SPACING_M / 2.0, 0.0, z_m])
        right = np.asarray([EAR_SPACING_M / 2.0, 0.0, z_m])
        for side, target, wrong in [("L", left, right), ("R", right, left)]:
            for method in ["phase", "wrong_region_null"]:
                if method == "phase":
                    weights = phase_weights_for_point(positions, normals, target)
                else:
                    weights = suppress_wrong_region_weights(positions, normals, target, wrong)
                metrics = evaluate_channel(positions, normals, target, wrong, weights)
                metrics.update({
                    "listener_z_m": z_m,
                    "ear_spacing_mm": EAR_SPACING_M * 1000.0,
                    "side": side,
                    "method": method,
                    "target_x_mm": target[0] * 1000.0,
                    "target_z_mm": target[2] * 1000.0,
                    "target_angle_deg": math.degrees(math.atan2(target[0], target[2])),
                })
                rows.append(metrics)
                if side == "R":
                    fwhm_x, xs, db = focus_width_x(positions, normals, target, weights)
                    line_plots.append((z_m, method, xs, db, fwhm_x))
    return rows, line_plots, positions, normals, segment_ids, tilts, centers


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_summary(rows, path):
    fig, axes = plt.subplots(3, 1, figsize=(8.0, 8.5), sharex=True)
    for method, style in [("phase", "o-"), ("wrong_region_null", "s--")]:
        subset = [r for r in rows if r["side"] == "R" and r["method"] == method]
        zs = [r["listener_z_m"] for r in subset]
        xtalk = [r["region_crosstalk_pressure_db"] for r in subset]
        audio = [r["region_crosstalk_audio_proxy_db"] for r in subset]
        fwhm = [r["fwhm_x_mm"] for r in subset]
        axes[0].plot(zs, xtalk, style, label=method)
        axes[1].plot(zs, audio, style, label=method)
        axes[2].plot(zs, fwhm, style, label=method)
    axes[0].set_ylabel("region pressure crosstalk (dB)")
    axes[0].axhline(-10.0, color="k", ls=":", lw=1)
    axes[1].set_ylabel("audio proxy crosstalk (dB)")
    axes[1].axhline(-20.0, color="k", ls=":", lw=1)
    axes[2].set_ylabel("horizontal FWHM (mm)")
    axes[2].set_xlabel("listener z distance (m)")
    for ax in axes:
        ax.grid(True, alpha=0.25)
        ax.legend()
    fig.suptitle("Stereo feasibility, ear spacing 180 mm, +65/+45/0/-45/-65")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_linecuts(line_plots, path):
    fig, axes = plt.subplots(len(LISTENER_Z_M), 1, figsize=(8.5, 8.5), sharex=False)
    for ax, z_m in zip(axes, LISTENER_Z_M):
        for zz, method, xs, db, fwhm in line_plots:
            if abs(zz - z_m) < 1e-9:
                ax.plot(xs, db, label="%s FWHM %.1f mm" % (method, fwhm))
        ax.axvline(EAR_SPACING_M * 500.0, color="k", lw=0.8, alpha=0.4)
        ax.axvline(-EAR_SPACING_M * 500.0, color="k", lw=0.8, alpha=0.4)
        ax.set_ylim(-45, 3)
        ax.set_title("right-ear beam line cut at z=%.1f m" % z_m)
        ax.set_ylabel("relative pressure (dB)")
        ax.grid(True, alpha=0.25)
        ax.legend()
    axes[-1].set_xlabel("x (mm)")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary(path, rows, tilts):
    with open(path, "w") as f:
        f.write("Stereo feasibility below 1 m\n")
        f.write("Array: connected folded 75-element geometry, tilts left-to-right %s deg.\n" % ("/".join("%+.1f" % t for t in tilts)))
        f.write("Listener model: ears at x=+/-%.0f mm, y=0, z in %.1f..%.1f m.\n" % (EAR_SPACING_M * 500.0, min(LISTENER_Z_M), max(LISTENER_Z_M)))
        f.write("Wrong-ear region: +/-%.0f mm square around the opposite ear.\n\n" % REGION_HALF_MM)
        f.write("Per-side rows are symmetric; R-side values shown below:\n")
        for row in rows:
            if row["side"] != "R":
                continue
            f.write("  z=%.1f m method=%-17s target_angle=%5.2f deg FWHMx=%6.1f mm pointXT=%7.2f dB regionXT=%7.2f dB audioRegionXT=%7.2f dB meanAmp=%5.2f\n" % (
                row["listener_z_m"], row["method"], row["target_angle_deg"], row["fwhm_x_mm"],
                row["point_crosstalk_pressure_db"], row["region_crosstalk_pressure_db"],
                row["region_crosstalk_audio_proxy_db"], row["mean_weight_amp"],
            ))


def main():
    ensure_out_dir()
    rows, line_plots, positions, normals, segment_ids, tilts, centers = run_scan()
    csv_path = os.path.join(OUT_DIR, "stereo_under1m_metrics.csv")
    write_csv(csv_path, rows)
    summary_png = os.path.join(OUT_DIR, "stereo_under1m_summary.png")
    plot_summary(rows, summary_png)
    line_png = os.path.join(OUT_DIR, "stereo_under1m_linecuts.png")
    plot_linecuts(line_plots, line_png)
    summary = os.path.join(OUT_DIR, "stereo_under1m_summary.txt")
    write_summary(summary, rows, tilts)
    for path in [summary, csv_path, summary_png, line_png]:
        print(path)


if __name__ == "__main__":
    main()
