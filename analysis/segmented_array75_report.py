import csv
import math
import os

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.special import j1


ROOT = os.path.dirname(os.path.dirname(__file__))
OUT_DIR = os.path.join(ROOT, "analysis_outputs")

C = 343.0
FREQ = 40000.0
LAMBDA = C / FREQ
K = 2.0 * np.pi / LAMBDA
ELEMENT_DIAMETER = 16e-3
ELEMENT_RADIUS = ELEMENT_DIAMETER / 2.0

SEGMENT_SPACING = 100e-3
SEGMENT_WIDTH = 100e-3
RING_RADII = [18e-3, 38e-3]
RING_COUNTS = [4, 10]
TILT_LIST_DEG = [-45.0, -22.5, 0.0, 22.5, 45.0]
THETA_MAX_DEG = 75.0
DB_FLOOR = -50.0
MAIN_EXCLUDE_MM = 20.0


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


def generate_segmented_array(mode):
    if mode in ("toe_in", "connected_toe_in"):
        tilts = [45.0, 22.5, 0.0, -22.5, -45.0]
    elif mode in ("listed_fan", "connected_fan"):
        tilts = TILT_LIST_DEG[:]
    else:
        raise ValueError(mode)

    local = local_15_positions()
    uxs = []
    normals_by_seg = []
    for tilt_deg in tilts:
        alpha = math.radians(tilt_deg)
        normals_by_seg.append(np.asarray([math.sin(alpha), 0.0, math.cos(alpha)]))
        uxs.append(np.asarray([math.cos(alpha), 0.0, -math.sin(alpha)]))

    centers = [None] * len(tilts)
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
    local_ids = []
    for seg, tilt_deg in enumerate(tilts):
        center = centers[seg]
        normal = normals_by_seg[seg]
        ux = uxs[seg]
        uy = np.asarray([0.0, 1.0, 0.0])
        for j, (xl, yl) in enumerate(local):
            positions.append(center + ux * xl + uy * yl)
            normals.append(normal)
            segment_ids.append(seg)
            local_ids.append(j)
    return np.asarray(positions), np.asarray(normals), np.asarray(segment_ids), np.asarray(local_ids), tilts


def farfield_matrix(positions, normals, theta, phi):
    theta = np.asarray(theta)
    phi = np.asarray(phi)
    sx = np.sin(theta) * np.cos(phi)
    sy = np.sin(theta) * np.sin(phi)
    sz = np.cos(theta)
    dirs = np.column_stack([sx.ravel(), sy.ravel(), sz.ravel()])
    cosang = np.dot(dirs, normals.T)
    elem = piston_from_cos(cosang)
    phase = np.exp(-1j * K * np.dot(dirs, positions.T))
    return elem * phase


def nearfield_matrix(positions, normals, x, y, z):
    points = np.column_stack([np.asarray(x).ravel(), np.asarray(y).ravel(), np.asarray(z).ravel()])
    out = np.zeros((len(points), len(positions)), dtype=complex)
    for i, (p, n) in enumerate(zip(positions, normals)):
        diff = points - p
        rr = np.linalg.norm(diff, axis=1)
        dirs = diff / np.maximum(rr[:, None], 1e-12)
        elem = piston_from_cos(np.dot(dirs, n))
        out[:, i] = elem * np.exp(1j * K * rr) / np.maximum(rr, 1e-9)
    return out


def phase_weights_for_farfield(positions, normals, theta_deg, phi_deg):
    h = farfield_matrix(positions, normals, np.deg2rad([theta_deg]), np.deg2rad([phi_deg]))[0]
    return np.exp(-1j * np.angle(h))


def phase_weights_for_target(positions, normals, target):
    h = nearfield_matrix(
        positions,
        normals,
        np.asarray([target[0]]),
        np.asarray([target[1]]),
        np.asarray([target[2]]),
    )[0]
    return np.exp(-1j * np.angle(h))


def response_db(resp, target_amp=None):
    amp = np.abs(resp)
    if target_amp is None:
        target_amp = float(np.max(amp))
    return np.maximum(20.0 * np.log10(amp / max(target_amp, 1e-18) + 1e-12), DB_FLOOR)


def write_layout_csv(path, positions, normals, seg_ids, local_ids, tilts):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Element_ID", "Segment", "Local_ID", "Segment_Tilt_deg", "X_mm", "Y_mm", "Z_mm", "Nx", "Ny", "Nz"])
        for i, (p, n, seg, local_id) in enumerate(zip(positions, normals, seg_ids, local_ids)):
            writer.writerow([
                "E%02d" % i,
                int(seg),
                int(local_id),
                "%.3f" % tilts[int(seg)],
                "%.3f" % (p[0] * 1000.0),
                "%.3f" % (p[1] * 1000.0),
                "%.3f" % (p[2] * 1000.0),
                "%.6f" % n[0],
                "%.6f" % n[1],
                "%.6f" % n[2],
            ])


def plot_geometry(positions, normals, seg_ids, path):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))
    ax = axes[0]
    for seg in sorted(set(seg_ids)):
        mask = seg_ids == seg
        ax.scatter(positions[mask, 0] * 1000.0, positions[mask, 1] * 1000.0, s=24, label="seg %d" % seg)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("Front view")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=3, fontsize=8)

    ax = axes[1]
    for seg in sorted(set(seg_ids)):
        mask = seg_ids == seg
        ax.scatter(positions[mask, 0] * 1000.0, positions[mask, 2] * 1000.0, s=24)
        c = np.mean(positions[mask], axis=0)
        n = normals[np.where(mask)[0][0]]
        ax.arrow(c[0] * 1000.0, c[2] * 1000.0, n[0] * 40.0, n[2] * 40.0, head_width=8.0, length_includes_head=True)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title("Side view and segment normals")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_orientation_compare():
    theta_deg = np.linspace(-90.0, 90.0, 1441)
    theta_abs = np.abs(np.deg2rad(theta_deg))
    phi = np.where(theta_deg >= 0.0, 0.0, np.pi)
    fig, ax = plt.subplots(figsize=(10.5, 5.4))
    rows = []
    for mode, label in [("listed_fan", "connected outward fan"), ("toe_in", "connected toe-in")]:
        pos, norm, seg, local, tilts = generate_segmented_array(mode)
        h = farfield_matrix(pos, norm, theta_abs, phi)
        w = np.ones(len(pos), dtype=complex)
        db = response_db(np.dot(h, w))
        ax.plot(theta_deg, db, lw=1.8, label=label)
        rows.append({
            "mode": mode,
            "level_0deg_db": float(db[np.argmin(np.abs(theta_deg - 0.0))]),
            "level_20deg_db": float(db[np.argmin(np.abs(theta_deg - 20.0))]),
            "level_45deg_db": float(db[np.argmin(np.abs(theta_deg - 45.0))]),
            "peak_angle_deg": float(theta_deg[int(np.argmax(db))]),
        })
    ax.set_xlim(-90, 90)
    ax.set_ylim(-50, 3)
    ax.set_xlabel("angle in x-z plane (deg)")
    ax.set_ylabel("normalized level (dB)")
    ax.set_title("75-element connected segmented array, all channels in phase")
    ax.grid(True, alpha=0.25)
    ax.legend()
    path = os.path.join(OUT_DIR, "seg75_orientation_compare.png")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return rows, path


def plot_circular_patterns(positions, normals):
    axis = np.linspace(-THETA_MAX_DEG, THETA_MAX_DEG, 241)
    tx, ty = np.meshgrid(axis, axis)
    theta_deg = np.sqrt(tx * tx + ty * ty)
    phi_deg = (np.rad2deg(np.arctan2(ty, tx)) + 360.0) % 360.0
    mask = theta_deg <= THETA_MAX_DEG
    theta = np.deg2rad(np.minimum(theta_deg, THETA_MAX_DEG))
    phi = np.deg2rad(phi_deg)
    h = farfield_matrix(positions, normals, theta, phi)
    cases = [
        ("All in phase", np.ones(len(positions), dtype=complex), 0.0, 0.0),
        ("Steer 0 deg", phase_weights_for_farfield(positions, normals, 0.0, 0.0), 0.0, 0.0),
        ("Steer 20 deg", phase_weights_for_farfield(positions, normals, 20.0, 0.0), 20.0, 0.0),
        ("Steer 45 deg", phase_weights_for_farfield(positions, normals, 45.0, 0.0), 45.0, 0.0),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 10.0))
    levels = np.linspace(DB_FLOOR, 0.0, 51)
    rows = []
    contour = None
    for ax, (title, w, target_theta, target_phi) in zip(axes.ravel(), cases):
        resp = np.dot(h, w).reshape(theta.shape)
        db = response_db(resp)
        db[~mask] = np.nan
        contour = ax.contourf(tx, ty, db, levels=levels, cmap="viridis", extend="min")
        ax.plot([target_theta * math.cos(math.radians(target_phi))], [target_theta * math.sin(math.radians(target_phi))], "rx", ms=8, mew=1.8)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(-THETA_MAX_DEG, THETA_MAX_DEG)
        ax.set_ylim(-THETA_MAX_DEG, THETA_MAX_DEG)
        ax.set_title(title)
        ax.set_xlabel("theta_x (deg)")
        ax.set_ylabel("theta_y (deg)")
        ax.grid(True, color="white", alpha=0.12)
        target_vec = np.asarray([
            math.sin(math.radians(target_theta)) * math.cos(math.radians(target_phi)),
            math.sin(math.radians(target_theta)) * math.sin(math.radians(target_phi)),
            math.cos(math.radians(target_theta)),
        ])
        ux = np.sin(theta) * np.cos(phi)
        uy = np.sin(theta) * np.sin(phi)
        uz = np.cos(theta)
        sep = np.rad2deg(np.arccos(np.clip(ux * target_vec[0] + uy * target_vec[1] + uz * target_vec[2], -1.0, 1.0)))
        outside = (sep > 8.0) & mask
        rows.append({"case": title, "max_outside_8deg_db": float(np.nanmax(db[outside]))})
    cbar = fig.colorbar(contour, ax=axes.ravel().tolist(), shrink=0.86, pad=0.03)
    cbar.set_label("normalized level (dB)")
    fig.suptitle("75-element connected toe-in segmented array, circular far-field patterns")
    path = os.path.join(OUT_DIR, "seg75_circular_patterns.png")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return rows, path


def local_focus_eval(positions, normals, target, span_mm=180.0, step_mm=3.0, method="phase"):
    x0, y0, z0 = target
    if method == "phase":
        w = phase_weights_for_target(positions, normals, target)
    elif method == "rls_amp":
        phase = phase_weights_for_target(positions, normals, target)
        x_mm = np.arange(x0 * 1000.0 - span_mm, x0 * 1000.0 + span_mm + 0.001, max(4.0, step_mm * 1.5))
        y_mm = np.arange(y0 * 1000.0 - span_mm, y0 * 1000.0 + span_mm + 0.001, max(4.0, step_mm * 1.5))
        xg, yg = np.meshgrid(x_mm * 1e-3, y_mm * 1e-3)
        zg = np.zeros_like(xg) + z0
        h = nearfield_matrix(positions, normals, xg, yg, zg)
        rr = np.sqrt((x_mm[None, :] - x0 * 1000.0) ** 2 + (y_mm[:, None] - y0 * 1000.0) ** 2)
        hs = h[rr.ravel() > MAIN_EXCLUDE_MM]
        b = hs * phase[None, :]
        ht = nearfield_matrix(positions, normals, np.asarray([x0]), np.asarray([y0]), np.asarray([z0]))[0]
        g = np.maximum(np.real(ht * phase), 0.0)
        g0 = float(np.sum(g))
        c = np.real(np.dot(b.conj().T, b)) / max(len(b), 1)
        c += 1e-4 * np.eye(len(positions))

        def obj(a):
            return float(np.dot(a, np.dot(c, a)))

        def jac(a):
            return 2.0 * np.dot(c, a)

        cons = [{"type": "ineq", "fun": lambda a: float(np.dot(g, a) - 0.72 * g0), "jac": lambda a: g}]
        res = minimize(obj, np.ones(len(positions)) * 0.72, jac=jac, bounds=[(0.0, 1.0)] * len(positions), constraints=cons, method="SLSQP", options={"maxiter": 120, "ftol": 1e-9})
        amp = np.clip(res.x if res.success else np.ones(len(positions)) * 0.72, 0.0, 1.0)
        w = amp * phase
    else:
        raise ValueError(method)

    x_mm = np.arange(x0 * 1000.0 - span_mm, x0 * 1000.0 + span_mm + 0.001, step_mm)
    y_mm = np.arange(y0 * 1000.0 - span_mm, y0 * 1000.0 + span_mm + 0.001, step_mm)
    xg, yg = np.meshgrid(x_mm * 1e-3, y_mm * 1e-3)
    zg = np.zeros_like(xg) + z0
    h = nearfield_matrix(positions, normals, xg, yg, zg)
    resp = np.dot(h, w).reshape(xg.shape)
    ht = nearfield_matrix(positions, normals, np.asarray([x0]), np.asarray([y0]), np.asarray([z0]))[0]
    target_amp = float(abs(np.dot(ht, w)))
    db = 20.0 * np.log10(np.abs(resp) / max(target_amp, 1e-18) + 1e-12)
    rr = np.sqrt((x_mm[None, :] - x0 * 1000.0) ** 2 + (y_mm[:, None] - y0 * 1000.0) ** 2)
    outside = rr > MAIN_EXCLUDE_MM
    idx = np.unravel_index(np.argmax(np.abs(resp)), resp.shape)
    ix = int(np.argmin(np.abs(x_mm - x0 * 1000.0)))
    iy = int(np.argmin(np.abs(y_mm - y0 * 1000.0)))
    return {
        "pressure_psl_db": float(np.max(db[outside])),
        "pressure_p99_5_db": float(np.percentile(db[outside], 99.5)),
        "focus_error_mm": float(math.hypot(x_mm[idx[1]] - x0 * 1000.0, y_mm[idx[0]] - y0 * 1000.0)),
        "fwhm_x_mm": line_width(x_mm, db[iy, :], ix, -3.0),
        "fwhm_y_mm": line_width(y_mm, db[:, ix], iy, -3.0),
        "audio_proxy_psl_db": float(2.0 * np.max(db[outside])),
        "target_amp": target_amp,
    }, x_mm, y_mm, np.maximum(db, DB_FLOOR)


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


def focus_scan(positions, normals):
    ref, _, _, _ = local_focus_eval(positions, normals, (0.0, 0.0, 0.220), span_mm=160.0, step_mm=4.0, method="phase")
    ref_amp = ref["target_amp"]
    z_values = [100, 150, 200, 220, 250, 300, 400, 500, 700, 1000]
    angles = [-45, -30, -22.5, -15, 0, 15, 22.5, 30, 45]
    rows = []
    for z_mm in z_values:
        for angle in angles:
            z = z_mm * 1e-3
            x = z * math.tan(math.radians(angle))
            span = max(160.0, min(320.0, 0.5 * z_mm))
            for method in ["phase", "rls_amp"]:
                metrics, _, _, _ = local_focus_eval(positions, normals, (x, 0.0, z), span_mm=span, step_mm=4.0 if z_mm <= 500 else 8.0, method=method)
                metrics.update({
                    "method": method,
                    "z_mm": z_mm,
                    "angle_deg": angle,
                    "x_mm": x * 1000.0,
                    "target_gain_vs_z220_onaxis_db": 20.0 * math.log10(metrics["target_amp"] / max(ref_amp, 1e-18) + 1e-12),
                })
                rows.append(metrics)
    return rows


def orientation_focus_compare():
    ref_positions, ref_normals, *_ = generate_segmented_array("toe_in")
    ref_metrics, _, _, _ = local_focus_eval(ref_positions, ref_normals, (0.0, 0.0, 0.220), span_mm=160.0, step_mm=4.0, method="phase")
    ref_amp = ref_metrics["target_amp"]
    targets = [(150, 0.0), (220, 0.0), (220, 22.5), (250, 30.0), (300, 0.0)]
    rows = []
    for mode in ["toe_in", "listed_fan"]:
        positions, normals, *_ = generate_segmented_array(mode)
        for z_mm, angle in targets:
            z = z_mm * 1e-3
            x = z * math.tan(math.radians(angle))
            span = max(160.0, min(320.0, 0.5 * z_mm))
            metrics, _, _, _ = local_focus_eval(positions, normals, (x, 0.0, z), span_mm=span, step_mm=4.0, method="phase")
            metrics.update({
                "mode": mode,
                "z_mm": z_mm,
                "angle_deg": angle,
                "x_mm": x * 1000.0,
                "target_gain_vs_toe_in_z220_onaxis_db": 20.0 * math.log10(metrics["target_amp"] / max(ref_amp, 1e-18) + 1e-12),
            })
            rows.append(metrics)
    return rows


def plot_focus_heatmaps(positions, normals):
    cases = [
        ("z=220 x=0 phase", (0.0, 0.0, 0.220), "phase"),
        ("z=220 x=0 rls", (0.0, 0.0, 0.220), "rls_amp"),
        ("z=250 angle=30 phase", (0.250 * math.tan(math.radians(30.0)), 0.0, 0.250), "phase"),
        ("z=250 angle=30 rls", (0.250 * math.tan(math.radians(30.0)), 0.0, 0.250), "rls_amp"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.8))
    levels = np.linspace(DB_FLOOR, 3.0, 54)
    contour = None
    for ax, (title, target, method) in zip(axes.ravel(), cases):
        _, x_mm, y_mm, db = local_focus_eval(positions, normals, target, span_mm=200.0, step_mm=3.0, method=method)
        contour = ax.contourf(x_mm, y_mm, db, levels=levels, cmap="magma", extend="both")
        ax.plot([target[0] * 1000.0], [target[1] * 1000.0], "cx", ms=8, mew=1.8)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(title)
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.grid(True, color="white", alpha=0.12)
    fig.subplots_adjust(right=0.86, hspace=0.32, wspace=0.22, top=0.90)
    cax = fig.add_axes([0.89, 0.16, 0.025, 0.68])
    fig.colorbar(contour, cax=cax).set_label("level relative to target (dB)")
    fig.suptitle("75-element connected toe-in segmented array, focus maps")
    path = os.path.join(OUT_DIR, "seg75_focus_heatmaps.png")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def usable_summary(rows, psl_limit=-10.0, gain_limit=-8.0, err_limit=8.0):
    out = []
    for method in ["phase", "rls_amp"]:
        for z_mm in sorted(set(r["z_mm"] for r in rows)):
            subset = [r for r in rows if r["method"] == method and r["z_mm"] == z_mm]
            usable = [
                r for r in subset
                if r["pressure_psl_db"] <= psl_limit
                and r["target_gain_vs_z220_onaxis_db"] >= gain_limit
                and r["focus_error_mm"] <= err_limit
            ]
            if usable:
                max_angle = max(abs(r["angle_deg"]) for r in usable)
                max_abs_x = max(abs(r["x_mm"]) for r in usable)
            else:
                max_angle = float("nan")
                max_abs_x = float("nan")
            out.append({"method": method, "z_mm": z_mm, "max_usable_angle_deg": max_angle, "max_usable_x_mm": max_abs_x})
    return out


def plot_focus_scan(rows, usable):
    fig, axes = plt.subplots(2, 1, figsize=(10.0, 7.6), sharex=True)
    for method, ls in [("phase", "-"), ("rls_amp", "--")]:
        for z in [150, 220, 300, 500]:
            subset = [r for r in rows if r["method"] == method and r["z_mm"] == z]
            axes[0].plot([r["angle_deg"] for r in subset], [r["pressure_psl_db"] for r in subset], marker="o", ls=ls, lw=1.4, label="%s z=%d" % (method, z))
            axes[1].plot([r["angle_deg"] for r in subset], [r["target_gain_vs_z220_onaxis_db"] for r in subset], marker="o", ls=ls, lw=1.4)
    axes[0].axhline(-10.0, color="k", ls=":", lw=1.0)
    axes[0].set_ylabel("pressure PSL (dB)")
    axes[0].set_title("Focus quality vs angle")
    axes[0].grid(True, alpha=0.25)
    axes[1].set_ylabel("target gain vs z=220 mm on-axis (dB)")
    axes[1].set_xlabel("focus angle (deg)")
    axes[1].grid(True, alpha=0.25)
    axes[0].legend(ncol=2, fontsize=8)
    fig.tight_layout()
    path1 = os.path.join(OUT_DIR, "seg75_focus_angle_scan.png")
    fig.savefig(path1, dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    for method in ["phase", "rls_amp"]:
        subset = [r for r in usable if r["method"] == method]
        ax.plot([r["z_mm"] for r in subset], [r["max_usable_angle_deg"] for r in subset], marker="o", lw=1.8, label=method)
    ax.set_xlabel("focus distance z (mm)")
    ax.set_ylabel("max usable absolute angle (deg)")
    ax.set_title("Usable focus angle, PSL <= -10 dB and gain >= -8 dB")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    path2 = os.path.join(OUT_DIR, "seg75_usable_angle_vs_distance.png")
    fig.savefig(path2, dpi=180)
    plt.close(fig)
    return path1, path2


def write_dict_csv(path, rows):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_summary(path, orientation_rows, circular_rows, focus_rows, usable, orientation_focus_rows=None):
    def find(z, angle, method):
        return next((r for r in focus_rows if r["z_mm"] == z and abs(r["angle_deg"] - angle) < 1e-6 and r["method"] == method), None)

    with open(path, "w") as f:
        f.write("75-element segmented array report\n")
        f.write("Assumption: 5 adjacent 100 mm segments connected edge-to-edge as a folded surface; each segment is center + 4@18 mm + 10@38 mm.\n")
        f.write("The simulated orientation is connected toe-in: left-to-right tilts +45/+22.5/0/-22.5/-45 deg, with the center segment at z=0.\n")
        f.write("Element model: 16 mm circular piston at 40 kHz, front hemisphere only.\n\n")
        f.write("All-in-phase orientation comparison:\n")
        for row in orientation_rows:
            f.write("  %s: peak=%5.1f deg, L0=%6.2f dB, L20=%6.2f dB, L45=%6.2f dB\n" % (
                row["mode"], row["peak_angle_deg"], row["level_0deg_db"], row["level_20deg_db"], row["level_45deg_db"]
            ))
        f.write("\nCircular far-field max outside 8 deg:\n")
        for row in circular_rows:
            f.write("  %-16s %7.2f dB\n" % (row["case"], row["max_outside_8deg_db"]))
        if orientation_focus_rows:
            f.write("\nNear-field orientation comparison, phase-only:\n")
            for row in orientation_focus_rows:
                f.write("  %-10s z=%3d angle=%5.1f x=%7.1f: PSL=%6.2f dB, gain=%6.2f dB, FWHM=(%5.1f,%5.1f)mm\n" % (
                    row["mode"], row["z_mm"], row["angle_deg"], row["x_mm"], row["pressure_psl_db"],
                    row["target_gain_vs_toe_in_z220_onaxis_db"], row["fwhm_x_mm"], row["fwhm_y_mm"]
                ))
        f.write("\nRepresentative focus points, phase-only:\n")
        for z, angle in [(150, 0), (220, 0), (220, 22.5), (250, 30), (300, 30), (500, 0), (500, 22.5)]:
            row = find(z, angle, "phase")
            if row:
                f.write("  z=%3d angle=%5.1f x=%7.1f: PSL=%6.2f dB, audio_proxy=%6.2f dB, gain=%6.2f dB, FWHM=(%5.1f,%5.1f)mm, err=%4.1fmm\n" % (
                    z, angle, row["x_mm"], row["pressure_psl_db"], row["audio_proxy_psl_db"], row["target_gain_vs_z220_onaxis_db"], row["fwhm_x_mm"], row["fwhm_y_mm"], row["focus_error_mm"]
                ))
        f.write("\nRepresentative focus points, rls_amp:\n")
        for z, angle in [(150, 0), (220, 0), (220, 22.5), (250, 30), (300, 30), (500, 0), (500, 22.5)]:
            row = find(z, angle, "rls_amp")
            if row:
                f.write("  z=%3d angle=%5.1f x=%7.1f: PSL=%6.2f dB, audio_proxy=%6.2f dB, gain=%6.2f dB, FWHM=(%5.1f,%5.1f)mm, err=%4.1fmm\n" % (
                    z, angle, row["x_mm"], row["pressure_psl_db"], row["audio_proxy_psl_db"], row["target_gain_vs_z220_onaxis_db"], row["fwhm_x_mm"], row["fwhm_y_mm"], row["focus_error_mm"]
                ))
        f.write("\nUsable angle by distance, criteria PSL <= -10 dB, gain >= -8 dB, err <= 8 mm:\n")
        for row in usable:
            f.write("  %-7s z=%4d mm: angle=%5s deg, x=%7s mm\n" % (
                row["method"], row["z_mm"],
                "nan" if math.isnan(row["max_usable_angle_deg"]) else "%.1f" % row["max_usable_angle_deg"],
                "nan" if math.isnan(row["max_usable_x_mm"]) else "%.1f" % row["max_usable_x_mm"],
            ))


def main():
    ensure_out_dir()
    orientation_rows, orientation_png = plot_orientation_compare()
    positions, normals, seg_ids, local_ids, tilts = generate_segmented_array("toe_in")
    layout_csv = os.path.join(OUT_DIR, "seg75_toe_in_layout.csv")
    write_layout_csv(layout_csv, positions, normals, seg_ids, local_ids, tilts)
    geometry_png = os.path.join(OUT_DIR, "seg75_geometry.png")
    plot_geometry(positions, normals, seg_ids, geometry_png)
    circular_rows, circular_png = plot_circular_patterns(positions, normals)
    orientation_focus_rows = orientation_focus_compare()
    orientation_focus_csv = os.path.join(OUT_DIR, "seg75_orientation_focus_compare.csv")
    write_dict_csv(orientation_focus_csv, orientation_focus_rows)
    focus_rows = focus_scan(positions, normals)
    focus_csv = os.path.join(OUT_DIR, "seg75_focus_scan.csv")
    write_dict_csv(focus_csv, focus_rows)
    usable = usable_summary(focus_rows)
    usable_csv = os.path.join(OUT_DIR, "seg75_usable_focus.csv")
    write_dict_csv(usable_csv, usable)
    heatmap_png = plot_focus_heatmaps(positions, normals)
    focus_angle_png, usable_png = plot_focus_scan(focus_rows, usable)
    summary = os.path.join(OUT_DIR, "seg75_summary.txt")
    write_summary(summary, orientation_rows, circular_rows, focus_rows, usable, orientation_focus_rows)
    for path in [summary, layout_csv, orientation_png, geometry_png, circular_png, orientation_focus_csv, focus_csv, usable_csv, heatmap_png, focus_angle_png, usable_png]:
        print(path)


if __name__ == "__main__":
    main()
