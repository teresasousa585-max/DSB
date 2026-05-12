import csv
import math
import os
import random
import time

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

import design_n12_unit as base

ROOT = os.path.dirname(os.path.dirname(__file__))
OUT_DIR = os.path.join(ROOT, "analysis_outputs")
LAYOUT_CSV = os.path.join(OUT_DIR, "agent_parametric_best.csv")

C = 343.0
FREQ = 40000.0
LAMBDA = C / FREQ
K = 2.0 * np.pi / LAMBDA
ELEMENT_DIAMETER_MM = 16.0
ELEMENT_RADIUS_M = ELEMENT_DIAMETER_MM * 0.5e-3
DB_FLOOR = -55.0

UNIT_ANGLES_DEG = [-45.0, -22.5, 0.0, 22.5, 45.0]

# Probe grid
PROBE_SPAN_DEG = 18.0
PROBE_N = 25
MAIN_EXCLUDE_DEG = 8.0

# GA parameters (amplitude only)
AMP_POP_SIZE = 50
AMP_N_GEN = 80
AMP_STAGNATION = 25
AMP_CROSSOVER_PROB = 0.9
AMP_MUTATION_PROB = 0.15
AMP_MUTATION_SIGMA = 0.08
AMP_ELITE = 2
AMP_REG_UNIFORM = 1e-3

# GA parameters (amplitude + phase, single case)
AP_POP_SIZE = 40
AP_N_GEN = 60
AP_STAGNATION = 20
AP_CROSSOVER_PROB = 0.9
AP_MUTATION_PROB = 0.15
AP_MUTATION_SIGMA_AMP = 0.08
AP_MUTATION_SIGMA_PHASE = 0.15
AP_ELITE = 2
AP_REG_UNIFORM = 1e-3

# Focal cases for shared-amplitude optimization
FOCUS_CASES = [
    (0.5, 0.0, 0.0),
    (0.5, 15.0, 0.0),
    (0.5, 30.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 15.0, 0.0),
    (1.0, 30.0, 0.0),
    (2.0, 0.0, 0.0),
    (2.0, 15.0, 0.0),
    (2.0, 30.0, 0.0),
    (3.0, 0.0, 0.0),
    (3.0, 15.0, 0.0),
    (3.0, 30.0, 0.0),
    (5.0, 0.0, 0.0),
    (5.0, 15.0, 0.0),
    (5.0, 30.0, 0.0),
]


def ensure_out_dir():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)


def read_layout(path):
    points = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            points.append([float(row["X_mm"]), float(row["Y_mm"])])
    return np.asarray(points, dtype=float)


def build_array_geometry(layout_mm):
    """
    Bowl-shaped (concave) array: units placed as tangent planes on a circular arc.
    Adjacent units meet at their edges, forming a continuous bowl.
    Arc center is at origin (zero point); normals point outward (away from origin)
    so the reverse extensions intersect at the origin.
    """
    unit_local = np.asarray(layout_mm, dtype=float)
    angles_rad = np.deg2rad(UNIT_ANGLES_DEG)
    n_units = len(UNIT_ANGLES_DEG)

    # Half-width of unit in x-direction (m)
    w_m = float(np.max(unit_local[:, 0]) - np.min(unit_local[:, 0])) * 0.5e-3

    # Angular spacing between adjacent units
    delta_rad = abs(angles_rad[1] - angles_rad[0])

    # Radius for edge-touching tangent planes: R = w / tan(delta/2)
    R_m = w_m / math.tan(delta_rad * 0.5)

    all_positions = []
    all_normals = []
    all_unit_ids = []
    all_element_ids = []

    for u, alpha in enumerate(angles_rad):
        # Center of unit lies on the lower circular arc in the x-z plane
        # Bowl opening faces +z; normals point toward +z (inward to the bowl center at origin)
        Cx = R_m * math.sin(alpha)
        Cy = 0.0
        Cz = -R_m * math.cos(alpha)
        center = np.array([Cx, Cy, Cz])

        # Local x-axis: tangent to arc, pointing in +alpha direction
        x_local = np.array([math.cos(alpha), 0.0, math.sin(alpha)])
        y_local = np.array([0.0, 1.0, 0.0])

        # Normal points toward origin (upward/outward from bowl surface toward +z)
        normal = np.array([-math.sin(alpha), 0.0, math.cos(alpha)])

        for e, (x_loc, y_loc) in enumerate(unit_local):
            x_loc_m = x_loc * 1e-3
            y_loc_m = y_loc * 1e-3
            P = center + x_loc_m * x_local + y_loc_m * y_local
            all_positions.append(P)
            all_normals.append(normal)
            all_unit_ids.append(u)
            all_element_ids.append(e)

    return {
        "positions": np.asarray(all_positions, dtype=float),
        "normals": np.asarray(all_normals, dtype=float),
        "unit_ids": np.asarray(all_unit_ids, dtype=int),
        "element_ids": np.asarray(all_element_ids, dtype=int),
        "n_units": n_units,
        "n_elements_per_unit": len(unit_local),
        "R_m": R_m,
    }


def probe_grid(distance_m, theta_x_deg, theta_y_deg, span_deg=PROBE_SPAN_DEG, n=PROBE_N):
    tx = np.linspace(-span_deg, span_deg, n)
    ty = np.linspace(-span_deg, span_deg, n)
    txs, tys = np.meshgrid(tx, ty)
    abs_tx = theta_x_deg + txs
    abs_ty = theta_y_deg + tys
    sx = np.sin(np.deg2rad(abs_tx))
    sy = np.sin(np.deg2rad(abs_ty))
    sz = np.sqrt(np.maximum(0.0, 1.0 - sx * sx - sy * sy))
    x = distance_m * sx
    y = distance_m * sy
    z = distance_m * sz
    return x, y, z, txs, tys


def build_case_data(geometry, distance_m, theta_x_deg, theta_y_deg):
    positions = geometry["positions"]
    normals = geometry["normals"]
    n_elements = len(positions)

    x, y, z, txs, tys = probe_grid(distance_m, theta_x_deg, theta_y_deg)
    probe_shape = x.shape
    n_probes = x.size
    probes = np.column_stack([x.ravel(), y.ravel(), z.ravel()])

    sx = math.sin(math.radians(theta_x_deg))
    sy = math.sin(math.radians(theta_y_deg))
    sz = math.sqrt(max(0.0, 1.0 - sx * sx - sy * sy))
    focal = distance_m * np.array([sx, sy, sz])

    dx = probes[:, 0][:, None] - positions[:, 0][None, :]
    dy = probes[:, 1][:, None] - positions[:, 1][None, :]
    dz = probes[:, 2][:, None] - positions[:, 2][None, :]
    dist = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    dist = np.maximum(dist, 1e-6)

    dir_x = dx / dist
    dir_y = dy / dist
    dir_z = dz / dist
    cos_theta = (dir_x * normals[None, :, 0] + dir_y * normals[None, :, 1] + dir_z * normals[None, :, 2])
    elem = base.piston_from_cos(cos_theta)

    dx_f = focal[0] - positions[:, 0]
    dy_f = focal[1] - positions[:, 1]
    dz_f = focal[2] - positions[:, 2]
    dist_f = np.sqrt(dx_f ** 2 + dy_f ** 2 + dz_f ** 2)
    dist_f = np.maximum(dist_f, 1e-6)

    h = elem * np.exp(1j * K * dist) / dist * np.exp(-1j * K * dist_f[None, :])

    # Focal point element pattern (for normalization)
    dir_fx = dx_f / dist_f
    dir_fy = dy_f / dist_f
    dir_fz = dz_f / dist_f
    cos_theta_f = dir_fx * normals[:, 0] + dir_fy * normals[:, 1] + dir_fz * normals[:, 2]
    elem_f = base.piston_from_cos(np.asarray([cos_theta_f]))[0]
    focal_resp_per_elem = elem_f / dist_f

    # Angular separation from focal direction
    target_dir = np.array([sx, sy, sz])
    probe_dirs = probes / distance_m
    cos_sep = np.dot(probe_dirs, target_dir)
    cos_sep = np.clip(cos_sep, -1.0, 1.0)
    sep_deg = np.rad2deg(np.arccos(cos_sep))
    outside = sep_deg > MAIN_EXCLUDE_DEG

    return {
        "h": h,
        "focal_resp_per_elem": focal_resp_per_elem,
        "outside": outside,
        "probe_shape": probe_shape,
        "txs": txs,
        "tys": tys,
        "distance_m": distance_m,
        "theta_x_deg": theta_x_deg,
        "theta_y_deg": theta_y_deg,
    }


def evaluate_amplitudes(amp, cases_data):
    amp = np.asarray(amp, dtype=float)
    max_sl = -1e9
    for case in cases_data:
        resp = np.dot(case["h"], amp)
        focal_resp = np.dot(case["focal_resp_per_elem"], amp)
        db = 20.0 * np.log10(np.abs(resp) / max(abs(focal_resp), 1e-18) + 1e-12)
        sl = float(np.max(db[case["outside"]]))
        max_sl = max(max_sl, sl)
    reg = AMP_REG_UNIFORM * float(np.mean((amp - 1.0) ** 2))
    return max_sl + reg


def evaluate_amp_phase(params, case_data):
    n = len(case_data["focal_resp_per_elem"])
    amp = np.clip(np.asarray(params[:n]), 0.0, 1.0)
    phase = np.asarray(params[n:])
    weights = amp * np.exp(1j * phase)
    resp = np.dot(case_data["h"], weights)
    center_idx = (case_data["probe_shape"][0] // 2) * case_data["probe_shape"][1] + case_data["probe_shape"][1] // 2
    focal_resp = resp[center_idx]
    target = max(abs(focal_resp), 1e-18)
    db = 20.0 * np.log10(np.abs(resp) / target + 1e-12)
    sl = float(np.max(db[case_data["outside"]]))
    reg = AP_REG_UNIFORM * float(np.mean((amp - 1.0) ** 2))
    return sl + reg


def tournament_select(pop, fitnesses, k=2):
    best = random.randrange(len(pop))
    for _ in range(k - 1):
        contender = random.randrange(len(pop))
        if fitnesses[contender] < fitnesses[best]:
            best = contender
    return pop[best].copy()


def blend_crossover(p1, p2, alpha=0.5):
    c1, c2 = p1.copy(), p2.copy()
    for i in range(len(p1)):
        if random.random() < 0.5:
            d = abs(p1[i] - p2[i])
            lo = min(p1[i], p2[i]) - alpha * d
            hi = max(p1[i], p2[i]) + alpha * d
            c1[i] = random.uniform(lo, hi)
            c2[i] = random.uniform(lo, hi)
    return c1, c2


def gaussian_mutate(individual, prob, sigma, low, high):
    mutant = individual.copy()
    for i in range(len(individual)):
        if random.random() < prob:
            mutant[i] += random.gauss(0.0, sigma)
    return np.clip(mutant, low, high)


def run_ga_amp_only(cases_data, n_elements, seed=20260511):
    random.seed(seed)
    np.random.seed(seed)

    pop = [np.random.uniform(0.0, 1.0, n_elements).tolist() for _ in range(AMP_POP_SIZE)]
    fitnesses = [None] * AMP_POP_SIZE
    evaluated = [False] * AMP_POP_SIZE

    best_hist = []
    mean_hist = []
    best_ever = None
    best_fitness_ever = float("inf")
    stagnation = 0

    start = time.time()
    for gen in range(AMP_N_GEN):
        for i in range(AMP_POP_SIZE):
            if not evaluated[i]:
                fitnesses[i] = evaluate_amplitudes(pop[i], cases_data)
                evaluated[i] = True

        sorted_idx = np.argsort(fitnesses)
        new_pop = [pop[idx].copy() for idx in sorted_idx[:AMP_ELITE]]
        new_fitness = [fitnesses[idx] for idx in sorted_idx[:AMP_ELITE]]
        new_evaluated = [True] * AMP_ELITE

        if fitnesses[sorted_idx[0]] < best_fitness_ever:
            best_fitness_ever = fitnesses[sorted_idx[0]]
            best_ever = pop[sorted_idx[0]].copy()
            stagnation = 0
        else:
            stagnation += 1

        best_hist.append(best_fitness_ever)
        mean_hist.append(float(np.mean(fitnesses)))

        print("AMP Gen %3d | best=%.4f | mean=%.4f | stag=%d" % (gen, best_fitness_ever, mean_hist[-1], stagnation))

        if stagnation >= AMP_STAGNATION:
            print("AMP GA stagnated at gen %d" % gen)
            break

        while len(new_pop) < AMP_POP_SIZE:
            p1 = tournament_select(pop, fitnesses)
            p2 = tournament_select(pop, fitnesses)
            if random.random() < AMP_CROSSOVER_PROB:
                c1, c2 = blend_crossover(p1, p2)
            else:
                c1, c2 = p1.copy(), p2.copy()
            c1 = gaussian_mutate(c1, AMP_MUTATION_PROB, AMP_MUTATION_SIGMA, 0.0, 1.0)
            c2 = gaussian_mutate(c2, AMP_MUTATION_PROB, AMP_MUTATION_SIGMA, 0.0, 1.0)
            new_pop.append(c1)
            new_fitness.append(None)
            new_evaluated.append(False)
            if len(new_pop) < AMP_POP_SIZE:
                new_pop.append(c2)
                new_fitness.append(None)
                new_evaluated.append(False)

        pop = new_pop
        fitnesses = new_fitness
        evaluated = new_evaluated

    elapsed = time.time() - start
    print("AMP GA finished in %.1f s" % elapsed)
    return np.clip(best_ever, 0.0, 1.0), best_hist, mean_hist


def run_ga_amp_phase(case_data, n_elements, seed=20260511):
    random.seed(seed)
    np.random.seed(seed)

    n_params = 2 * n_elements
    pop = []
    for _ in range(AP_POP_SIZE):
        ind = np.random.uniform(0.0, 1.0, n_elements).tolist()
        ind += np.random.uniform(-np.pi, np.pi, n_elements).tolist()
        pop.append(ind)

    fitnesses = [None] * AP_POP_SIZE
    evaluated = [False] * AP_POP_SIZE

    best_hist = []
    mean_hist = []
    best_ever = None
    best_fitness_ever = float("inf")
    stagnation = 0

    start = time.time()
    for gen in range(AP_N_GEN):
        for i in range(AP_POP_SIZE):
            if not evaluated[i]:
                fitnesses[i] = evaluate_amp_phase(pop[i], case_data)
                evaluated[i] = True

        sorted_idx = np.argsort(fitnesses)
        new_pop = [pop[idx].copy() for idx in sorted_idx[:AP_ELITE]]
        new_fitness = [fitnesses[idx] for idx in sorted_idx[:AP_ELITE]]
        new_evaluated = [True] * AP_ELITE

        if fitnesses[sorted_idx[0]] < best_fitness_ever:
            best_fitness_ever = fitnesses[sorted_idx[0]]
            best_ever = pop[sorted_idx[0]].copy()
            stagnation = 0
        else:
            stagnation += 1

        best_hist.append(best_fitness_ever)
        mean_hist.append(float(np.mean(fitnesses)))

        print("AP Gen %3d | best=%.4f | mean=%.4f | stag=%d" % (gen, best_fitness_ever, mean_hist[-1], stagnation))

        if stagnation >= AP_STAGNATION:
            print("AP GA stagnated at gen %d" % gen)
            break

        while len(new_pop) < AP_POP_SIZE:
            p1 = tournament_select(pop, fitnesses)
            p2 = tournament_select(pop, fitnesses)
            if random.random() < AP_CROSSOVER_PROB:
                c1, c2 = blend_crossover(p1, p2)
            else:
                c1, c2 = p1.copy(), p2.copy()
            c1 = gaussian_mutate(c1[:n_elements], AP_MUTATION_PROB, AP_MUTATION_SIGMA_AMP, 0.0, 1.0).tolist()
            c1 += gaussian_mutate(np.asarray(best_ever[n_elements:]) if best_ever is not None else np.zeros(n_elements), AP_MUTATION_PROB, AP_MUTATION_SIGMA_PHASE, -np.pi, np.pi).tolist()
            # Actually the mutate should apply to the full individual
            c1_full = np.array(c1)
            c1_full[:n_elements] = gaussian_mutate(c1_full[:n_elements], AP_MUTATION_PROB, AP_MUTATION_SIGMA_AMP, 0.0, 1.0)
            c1_full[n_elements:] = gaussian_mutate(c1_full[n_elements:], AP_MUTATION_PROB, AP_MUTATION_SIGMA_PHASE, -np.pi, np.pi)
            c2_full = np.array(c2)
            c2_full[:n_elements] = gaussian_mutate(c2_full[:n_elements], AP_MUTATION_PROB, AP_MUTATION_SIGMA_AMP, 0.0, 1.0)
            c2_full[n_elements:] = gaussian_mutate(c2_full[n_elements:], AP_MUTATION_PROB, AP_MUTATION_SIGMA_PHASE, -np.pi, np.pi)
            new_pop.append(c1_full.tolist())
            new_fitness.append(None)
            new_evaluated.append(False)
            if len(new_pop) < AP_POP_SIZE:
                new_pop.append(c2_full.tolist())
                new_fitness.append(None)
                new_evaluated.append(False)

        pop = new_pop
        fitnesses = new_fitness
        evaluated = new_evaluated

    elapsed = time.time() - start
    print("AP GA finished in %.1f s" % elapsed)
    best = np.asarray(best_ever)
    amp = np.clip(best[:n_elements], 0.0, 1.0)
    phase = best[n_elements:]
    return amp, phase, best_hist, mean_hist


def metrics_for_case(amp, case_data):
    amp = np.asarray(amp, dtype=float)
    resp = np.dot(case_data["h"], amp)
    focal_resp = np.dot(case_data["focal_resp_per_elem"], amp)
    target = max(abs(focal_resp), 1e-18)
    db = 20.0 * np.log10(np.abs(resp) / target + 1e-12)
    db = db.reshape(case_data["probe_shape"])
    outside_db = db.ravel()[case_data["outside"]]
    return {
        "max_sidelobe_db": float(np.max(outside_db)),
        "p99_sidelobe_db": float(np.percentile(outside_db, 99.0)),
        "p95_sidelobe_db": float(np.percentile(outside_db, 95.0)),
        "mean_amp": float(np.mean(amp)),
        "min_amp": float(np.min(amp)),
        "max_amp": float(np.max(amp)),
        "db": db,
    }


def metrics_for_case_amp_phase(amp, phase, case_data):
    weights = amp * np.exp(1j * phase)
    resp = np.dot(case_data["h"], weights)
    center_idx = (case_data["probe_shape"][0] // 2) * case_data["probe_shape"][1] + case_data["probe_shape"][1] // 2
    focal_resp = resp[center_idx]
    target = max(abs(focal_resp), 1e-18)
    db = 20.0 * np.log10(np.abs(resp) / target + 1e-12)
    db = db.reshape(case_data["probe_shape"])
    outside_db = db.ravel()[case_data["outside"]]
    return {
        "max_sidelobe_db": float(np.max(outside_db)),
        "p99_sidelobe_db": float(np.percentile(outside_db, 99.0)),
        "p95_sidelobe_db": float(np.percentile(outside_db, 95.0)),
        "mean_amp": float(np.mean(amp)),
        "min_amp": float(np.min(amp)),
        "max_amp": float(np.max(amp)),
        "db": db,
    }


def plot_geometry(geometry, path):
    positions = geometry["positions"]
    unit_ids = geometry["unit_ids"]
    n_units = geometry["n_units"]

    fig, ax = plt.subplots(figsize=(10.0, 7.0))
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
    for u in range(n_units):
        mask = unit_ids == u
        pts = positions[mask]
        ax.scatter(pts[:, 0] * 1000.0, pts[:, 2] * 1000.0, s=40, color=colors[u], label="Unit %d (%.1f deg)" % (u, UNIT_ANGLES_DEG[u]), zorder=3)
        for p in pts:
            circ = plt.Circle((p[0] * 1000.0, p[2] * 1000.0), ELEMENT_DIAMETER_MM * 0.5, fill=False, lw=0.8, color=colors[u], alpha=0.5)
            ax.add_patch(circ)
        # Draw normal arrow from unit center
        center = np.mean(pts, axis=0)
        normal = geometry["normals"][mask][0]
        ax.arrow(center[0] * 1000.0, center[2] * 1000.0,
                 normal[0] * 120.0, normal[2] * 120.0,
                 head_width=8.0, head_length=10.0, fc=colors[u], ec=colors[u], alpha=0.7, zorder=2)

    ax.plot([0], [0], "k+", ms=12, mew=2.0, label="Zero point (normal intersection)", zorder=4)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title("5-unit bowl-shaped array geometry (x-z view), R=%.3f m" % geometry["R_m"])
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_patterns(rows, path, title):
    n = len(rows)
    fig, axes = plt.subplots(n, 2, figsize=(11.0, 3.8 * n))
    if n == 1:
        axes = axes.reshape(1, 2)
    fig.subplots_adjust(hspace=0.45, wspace=0.35)
    levels = np.linspace(DB_FLOOR, 0.0, 56)
    contour = None
    for r, row in enumerate(rows):
        for c, (label, db, data) in enumerate([
            ("uniform", row["uniform_db"], row["uniform"]),
            ("optimized", row["opt_db"], row["opt"]),
        ]):
            ax = axes[r, c]
            plot_db = np.maximum(db, DB_FLOOR)
            contour = ax.contourf(row["txs"], row["tys"], plot_db, levels=levels, cmap="viridis", extend="min")
            ax.plot([0], [0], "rx", ms=7, mew=1.6)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlim(-PROBE_SPAN_DEG, PROBE_SPAN_DEG)
            ax.set_ylim(-PROBE_SPAN_DEG, PROBE_SPAN_DEG)
            ax.set_title("%s %s @ %.1fm, %.0f deg\nSL %.1f dB, p99 %.1f dB" % (
                row["case"], label, row["distance_m"], row["theta_x_deg"],
                data["max_sidelobe_db"], data["p99_sidelobe_db"],
            ))
            ax.set_xlabel("offset theta_x (deg)")
            ax.set_ylabel("offset theta_y (deg)")
            ax.grid(True, color="white", alpha=0.12)
    cbar = fig.colorbar(contour, ax=axes, shrink=0.75)
    cbar.set_label("relative pressure (dB)")
    fig.suptitle(title)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_weights(amp, geometry, path, title):
    n_units = geometry["n_units"]
    n_per = geometry["n_elements_per_unit"]
    fig, ax = plt.subplots(figsize=(10.0, 5.5))
    x = np.arange(n_per)
    width = 0.14
    for u in range(n_units):
        vals = amp[u * n_per:(u + 1) * n_per]
        ax.bar(x + (u - (n_units - 1) / 2.0) * width, vals, width=width, label="Unit %d (%.1f deg)" % (u, UNIT_ANGLES_DEG[u]))
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("element index within unit")
    ax.set_ylabel("optimized amplitude")
    ax.set_xticks(x)
    ax.set_xticklabels(["E%02d" % i for i in x])
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(ncol=3, fontsize=8)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_history(best_hist, mean_hist, path, title):
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    ax.plot(best_hist, lw=1.8, label="best fitness")
    ax.plot(mean_hist, lw=1.5, alpha=0.7, label="mean fitness")
    ax.set_xlabel("generation")
    ax.set_ylabel("fitness (dB)")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_amp_weights(path, geometry, amp):
    with open(path, "w", newline="") as f:
        fields = ["global_id", "unit_id", "unit_angle_deg", "element_id", "x_m", "y_m", "z_m", "amplitude"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for i in range(len(amp)):
            p = geometry["positions"][i]
            writer.writerow({
                "global_id": i,
                "unit_id": geometry["unit_ids"][i],
                "unit_angle_deg": UNIT_ANGLES_DEG[geometry["unit_ids"][i]],
                "element_id": geometry["element_ids"][i],
                "x_m": "%.6f" % p[0],
                "y_m": "%.6f" % p[1],
                "z_m": "%.6f" % p[2],
                "amplitude": "%.6f" % amp[i],
            })


def write_amp_phase_weights(path, geometry, amp, phase):
    with open(path, "w", newline="") as f:
        fields = ["global_id", "unit_id", "unit_angle_deg", "element_id", "x_m", "y_m", "z_m", "amplitude", "phase_rad", "phase_deg"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for i in range(len(amp)):
            p = geometry["positions"][i]
            writer.writerow({
                "global_id": i,
                "unit_id": geometry["unit_ids"][i],
                "unit_angle_deg": UNIT_ANGLES_DEG[geometry["unit_ids"][i]],
                "element_id": geometry["element_ids"][i],
                "x_m": "%.6f" % p[0],
                "y_m": "%.6f" % p[1],
                "z_m": "%.6f" % p[2],
                "amplitude": "%.6f" % amp[i],
                "phase_rad": "%.6f" % phase[i],
                "phase_deg": "%.3f" % np.rad2deg(phase[i]),
            })


def write_metrics(path, rows):
    fields = [
        "case", "distance_m", "theta_x_deg", "theta_y_deg", "method",
        "max_sidelobe_db", "p99_sidelobe_db", "p95_sidelobe_db",
        "mean_amp", "min_amp", "max_amp",
    ]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            for method, data in [("uniform", row["uniform"]), ("optimized", row["opt"])]:
                out = {
                    "case": row["case"],
                    "distance_m": row["distance_m"],
                    "theta_x_deg": row["theta_x_deg"],
                    "theta_y_deg": row["theta_y_deg"],
                    "method": method,
                }
                for k in fields:
                    if k in data:
                        out[k] = "%.6f" % data[k]
                writer.writerow(out)


def main():
    ensure_out_dir()
    t0 = time.time()

    print("Loading layout from %s ..." % LAYOUT_CSV)
    layout = read_layout(LAYOUT_CSV)
    geometry = build_array_geometry(layout)
    n_elements = len(geometry["positions"])
    print("Array: %d units x %d elements = %d total elements" % (
        geometry["n_units"], geometry["n_elements_per_unit"], n_elements))

    print("Building case data for %d focal cases ..." % len(FOCUS_CASES))
    cases_data = []
    for dist, tx, ty in FOCUS_CASES:
        cases_data.append(build_case_data(geometry, dist, tx, ty))

    # --- GA 1: shared amplitude optimization across all focal cases ---
    print("\n=== GA: shared amplitude optimization (%d cases) ===" % len(FOCUS_CASES))
    best_amp, best_hist, mean_hist = run_ga_amp_only(cases_data, n_elements)

    # Evaluate shared amplitude on all cases
    shared_rows = []
    for case_data in cases_data:
        uniform_amp = np.ones(n_elements, dtype=float)
        uniform_metrics = metrics_for_case(uniform_amp, case_data)
        opt_metrics = metrics_for_case(best_amp, case_data)
        shared_rows.append({
            "case": "R%.1f_x%.0f_y%.0f" % (case_data["distance_m"], case_data["theta_x_deg"], case_data["theta_y_deg"]),
            "distance_m": case_data["distance_m"],
            "theta_x_deg": case_data["theta_x_deg"],
            "theta_y_deg": case_data["theta_y_deg"],
            "uniform": uniform_metrics,
            "uniform_db": uniform_metrics.pop("db"),
            "opt": opt_metrics,
            "opt_db": opt_metrics.pop("db"),
            "txs": case_data["txs"],
            "tys": case_data["tys"],
        })

    # --- GA 2: amplitude + phase for center case at 2m ---
    center_case = None
    for cd in cases_data:
        if cd["distance_m"] == 2.0 and cd["theta_x_deg"] == 0.0 and cd["theta_y_deg"] == 0.0:
            center_case = cd
            break

    if center_case is not None:
        print("\n=== GA: amplitude + phase optimization (center 2m case) ===")
        best_amp_ap, best_phase_ap, best_hist_ap, mean_hist_ap = run_ga_amp_phase(center_case, n_elements)
        ap_uniform = metrics_for_case(np.ones(n_elements), center_case)
        ap_opt = metrics_for_case_amp_phase(best_amp_ap, best_phase_ap, center_case)
        ap_rows = [{
            "case": "center_2m_amp_phase",
            "distance_m": 2.0,
            "theta_x_deg": 0.0,
            "theta_y_deg": 0.0,
            "uniform": ap_uniform,
            "uniform_db": ap_uniform.pop("db"),
            "opt": ap_opt,
            "opt_db": ap_opt.pop("db"),
            "txs": center_case["txs"],
            "tys": center_case["tys"],
        }]
    else:
        ap_rows = []
        best_amp_ap = None
        best_phase_ap = None
        best_hist_ap = []
        mean_hist_ap = []

    # --- Write outputs ---
    print("\nWriting outputs ...")

    geo_png = os.path.join(OUT_DIR, "multi_unit_geometry.png")
    plot_geometry(geometry, geo_png)

    shared_patterns_png = os.path.join(OUT_DIR, "multi_unit_shared_amp_patterns.png")
    plot_patterns(shared_rows[:6], shared_patterns_png,
                  "Shared amplitude optimization: uniform vs optimized (first 6 cases)")

    shared_weights_png = os.path.join(OUT_DIR, "multi_unit_shared_amp_weights.png")
    plot_weights(best_amp, geometry, shared_weights_png,
                 "Shared optimized amplitudes across %d focal cases" % len(FOCUS_CASES))

    shared_history_png = os.path.join(OUT_DIR, "multi_unit_shared_amp_history.png")
    plot_history(best_hist, mean_hist, shared_history_png,
                 "Shared amplitude GA convergence (%d cases)" % len(FOCUS_CASES))

    shared_weights_csv = os.path.join(OUT_DIR, "multi_unit_shared_amp_weights.csv")
    write_amp_weights(shared_weights_csv, geometry, best_amp)

    shared_metrics_csv = os.path.join(OUT_DIR, "multi_unit_shared_amp_metrics.csv")
    write_metrics(shared_metrics_csv, shared_rows)

    if ap_rows:
        ap_patterns_png = os.path.join(OUT_DIR, "multi_unit_amp_phase_patterns.png")
        plot_patterns(ap_rows, ap_patterns_png,
                      "Amplitude + phase optimization: center 2m case")

        ap_weights_png = os.path.join(OUT_DIR, "multi_unit_amp_phase_weights.png")
        plot_weights(best_amp_ap, geometry, ap_weights_png,
                     "Amplitude + phase optimized weights (center 2m)")

        ap_history_png = os.path.join(OUT_DIR, "multi_unit_amp_phase_history.png")
        plot_history(best_hist_ap, mean_hist_ap, ap_history_png,
                     "Amplitude + phase GA convergence (center 2m)")

        ap_weights_csv = os.path.join(OUT_DIR, "multi_unit_amp_phase_weights.csv")
        write_amp_phase_weights(ap_weights_csv, geometry, best_amp_ap, best_phase_ap)

    summary_path = os.path.join(OUT_DIR, "multi_unit_focus_summary.txt")
    with open(summary_path, "w") as f:
        f.write("Multi-unit bowl-shaped focused array optimization report\n")
        f.write("=========================================================\n")
        f.write("Units: %d x %d = %d elements\n" % (geometry["n_units"], geometry["n_elements_per_unit"], n_elements))
        f.write("Unit angles (deg): %s\n" % ", ".join("%.1f" % a for a in UNIT_ANGLES_DEG))
        f.write("Bowl radius R (edge-touching tangent planes): %.3f m\n" % geometry["R_m"])
        f.write("Focus cases: %d cases (distances %.1f-%.1f m, angles 0-30 deg)\n" % (
            len(FOCUS_CASES), min(d for d, _, _ in FOCUS_CASES), max(d for d, _, _ in FOCUS_CASES)))
        f.write("Probe grid: %d x %d over +/-%.1f deg, main exclusion %.1f deg\n" % (
            PROBE_N, PROBE_N, PROBE_SPAN_DEG, MAIN_EXCLUDE_DEG))
        f.write("GA shared amp: pop=%d, gen=%d (with stagnation limit %d)\n" % (AMP_POP_SIZE, AMP_N_GEN, AMP_STAGNATION))
        f.write("GA amp+phase: pop=%d, gen=%d (with stagnation limit %d)\n\n" % (AP_POP_SIZE, AP_N_GEN, AP_STAGNATION))

        f.write("Shared amplitude optimization results:\n")
        f.write("%-20s %6s %6s %6s %10s %10s %10s %10s\n" % (
            "case", "dist", "tx", "ty", "u_SL_dB", "u_p99_dB", "o_SL_dB", "o_p99_dB"))
        for row in shared_rows:
            u = row["uniform"]
            o = row["opt"]
            f.write("%-20s %6.1f %6.1f %6.1f %10.2f %10.2f %10.2f %10.2f\n" % (
                row["case"], row["distance_m"], row["theta_x_deg"], row["theta_y_deg"],
                u["max_sidelobe_db"], u["p99_sidelobe_db"],
                o["max_sidelobe_db"], o["p99_sidelobe_db"],
            ))

        if ap_rows:
            f.write("\nAmplitude + phase optimization (center 2m case):\n")
            u = ap_rows[0]["uniform"]
            o = ap_rows[0]["opt"]
            f.write("  uniform SL=%.2f p99=%.2f | opt SL=%.2f p99=%.2f\n" % (
                u["max_sidelobe_db"], u["p99_sidelobe_db"],
                o["max_sidelobe_db"], o["p99_sidelobe_db"],
            ))

        f.write("\nOptimized shared amplitudes per element:\n")
        for u in range(geometry["n_units"]):
            start = u * geometry["n_elements_per_unit"]
            end = start + geometry["n_elements_per_unit"]
            vals = best_amp[start:end]
            f.write("  Unit %d (%.1f deg): %s\n" % (u, UNIT_ANGLES_DEG[u], ", ".join("%.3f" % v for v in vals)))

    elapsed = time.time() - t0
    print("\nTotal runtime: %.1f s" % elapsed)
    print("Outputs:")
    for p in [geo_png, shared_patterns_png, shared_weights_png, shared_history_png,
              shared_weights_csv, shared_metrics_csv, summary_path]:
        print("  " + p)
    if ap_rows:
        for p in [ap_patterns_png, ap_weights_png, ap_history_png, ap_weights_csv]:
            print("  " + p)


if __name__ == "__main__":
    main()
