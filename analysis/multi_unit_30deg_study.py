import csv
import math
import os
import random

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
DB_FLOOR = -55.0

UNIT_ANGLES_DEG = [-45.0, -22.5, 0.0, 22.5, 45.0]

PROBE_SPAN_DEG = 18.0
PROBE_N = 25
MAIN_EXCLUDE_DEG = 8.0

AP_POP_SIZE = 50
AP_N_GEN = 80
AP_STAGNATION = 25
AP_CROSSOVER_PROB = 0.9
AP_MUTATION_PROB = 0.15
AP_MUTATION_SIGMA_AMP = 0.08
AP_MUTATION_SIGMA_PHASE = 0.12
AP_ELITE = 2
AP_REG_UNIFORM = 1e-3

THETA_X = 30.0
THETA_Y = 0.0
DISTANCES = [0.5, 1.0, 2.0, 3.0, 5.0]


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
    unit_local = np.asarray(layout_mm, dtype=float)
    angles_rad = np.deg2rad(UNIT_ANGLES_DEG)
    n_units = len(UNIT_ANGLES_DEG)
    w_m = float(np.max(unit_local[:, 0]) - np.min(unit_local[:, 0])) * 0.5e-3
    delta_rad = abs(angles_rad[1] - angles_rad[0])
    R_m = w_m / math.tan(delta_rad * 0.5)

    all_positions = []
    all_normals = []
    all_unit_ids = []
    all_element_ids = []

    for u, alpha in enumerate(angles_rad):
        Cx = R_m * math.sin(alpha)
        Cy = 0.0
        Cz = -R_m * math.cos(alpha)
        center = np.array([Cx, Cy, Cz])
        x_local = np.array([math.cos(alpha), 0.0, math.sin(alpha)])
        y_local = np.array([0.0, 1.0, 0.0])
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

    dir_fx = dx_f / dist_f
    dir_fy = dy_f / dist_f
    dir_fz = dz_f / dist_f
    cos_theta_f = dir_fx * normals[:, 0] + dir_fy * normals[:, 1] + dir_fz * normals[:, 2]
    elem_f = base.piston_from_cos(np.asarray([cos_theta_f]))[0]
    focal_resp_per_elem = elem_f / dist_f

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

    best = np.asarray(best_ever)
    amp = np.clip(best[:n_elements], 0.0, 1.0)
    phase = best[n_elements:]
    return amp, phase, best_hist, mean_hist


def metrics_for_case(amp, phase, case_data):
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


def compute_cos_theta_to_focus(geometry, distance_m, theta_x_deg, theta_y_deg):
    """Compute cos(theta) between each element normal and focal direction."""
    sx = math.sin(math.radians(theta_x_deg))
    sy = math.sin(math.radians(theta_y_deg))
    sz = math.sqrt(max(0.0, 1.0 - sx * sx - sy * sy))
    focal_dir = np.array([sx, sy, sz])
    normals = geometry["normals"]
    cos_vals = np.dot(normals, focal_dir)
    return cos_vals


def heuristic_amp_for_30deg(geometry, power=1.0):
    """Directional heuristic: weight elements by cos(theta)^power to 30 deg focus."""
    # Use a representative distance (focus direction independent of distance)
    cos_vals = compute_cos_theta_to_focus(geometry, 2.0, THETA_X, THETA_Y)
    amp = np.clip(cos_vals ** power, 0.0, 1.0)
    return amp


def plot_comparison(rows, path):
    n = len(rows)
    fig, axes = plt.subplots(n, 3, figsize=(13.5, 3.8 * n))
    if n == 1:
        axes = axes.reshape(1, 3)
    fig.subplots_adjust(hspace=0.45, wspace=0.35)
    levels = np.linspace(DB_FLOOR, 0.0, 56)
    contour = None
    for r, row in enumerate(rows):
        for c, (label, db, data) in enumerate([
            ("uniform", row["uniform_db"], row["uniform"]),
            ("heuristic", row["heuristic_db"], row["heuristic"]),
            ("GA opt", row["opt_db"], row["opt"]),
        ]):
            ax = axes[r, c]
            plot_db = np.maximum(db, DB_FLOOR)
            contour = ax.contourf(row["txs"], row["tys"], plot_db, levels=levels, cmap="viridis", extend="min")
            ax.plot([0], [0], "rx", ms=7, mew=1.6)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlim(-PROBE_SPAN_DEG, PROBE_SPAN_DEG)
            ax.set_ylim(-PROBE_SPAN_DEG, PROBE_SPAN_DEG)
            ax.set_title("%s %s\nSL %.1f dB, p99 %.1f dB" % (
                row["case"], label, data["max_sidelobe_db"], data["p99_sidelobe_db"],
            ))
            ax.set_xlabel("offset theta_x (deg)")
            ax.set_ylabel("offset theta_y (deg)")
            ax.grid(True, color="white", alpha=0.12)
    cbar = fig.colorbar(contour, ax=axes, shrink=0.75)
    cbar.set_label("relative pressure (dB)")
    fig.suptitle("30-degree steering: uniform vs heuristic vs GA optimized")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_weights_comparison(rows, geometry, path):
    n_units = geometry["n_units"]
    n_per = geometry["n_elements_per_unit"]
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.5))
    methods = [("uniform", np.ones(n_per * n_units)),
               ("heuristic", rows[0]["heuristic_amp"]),
               ("GA opt", rows[0]["opt_amp"])]
    x = np.arange(n_per)
    width = 0.14
    for ax, (label, amp) in zip(axes, methods):
        for u in range(n_units):
            vals = amp[u * n_per:(u + 1) * n_per]
            ax.bar(x + (u - (n_units - 1) / 2.0) * width, vals, width=width,
                   label="Unit %d (%.1f deg)" % (u, UNIT_ANGLES_DEG[u]))
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("element index")
        ax.set_ylabel("amplitude")
        ax.set_xticks(x)
        ax.set_xticklabels(["E%02d" % i for i in x])
        ax.grid(True, axis="y", alpha=0.25)
        ax.set_title(label)
        if ax == axes[0]:
            ax.legend(ncol=3, fontsize=7)
    fig.suptitle("Amplitude weights for 30-degree steering (%s)" % rows[0]["case"])
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_summary(path, rows, geometry):
    with open(path, "w") as f:
        f.write("30-degree steering optimization study\n")
        f.write("=====================================\n")
        f.write("Bowl radius R: %.3f m\n" % geometry["R_m"])
        f.write("Focus direction: theta_x = %.1f deg, theta_y = %.1f deg\n\n" % (THETA_X, THETA_Y))
        f.write("%-12s %6s %10s %10s %10s %10s %10s %10s\n" % (
            "case", "dist", "u_SL", "u_p99", "h_SL", "h_p99", "o_SL", "o_p99"
        ))
        for row in rows:
            f.write("%-12s %6.1f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n" % (
                row["case"], row["distance_m"],
                row["uniform"]["max_sidelobe_db"], row["uniform"]["p99_sidelobe_db"],
                row["heuristic"]["max_sidelobe_db"], row["heuristic"]["p99_sidelobe_db"],
                row["opt"]["max_sidelobe_db"], row["opt"]["p99_sidelobe_db"],
            ))
        f.write("\n")
        f.write("Unit cos(theta) to 30-deg direction:\n")
        cos_vals = compute_cos_theta_to_focus(geometry, 2.0, THETA_X, THETA_Y)
        for u in range(geometry["n_units"]):
            start = u * geometry["n_elements_per_unit"]
            c = cos_vals[start]
            f.write("  Unit %d (%.1f deg): cos_theta = %.3f\n" % (u, UNIT_ANGLES_DEG[u], c))


def main():
    ensure_out_dir()
    layout = read_layout(LAYOUT_CSV)
    geometry = build_array_geometry(layout)
    n_elements = len(geometry["positions"])
    print("Array: %d units x %d elements = %d total" % (
        geometry["n_units"], geometry["n_elements_per_unit"], n_elements))

    # Precompute heuristic amplitude (directional weighting)
    heuristic_amp = heuristic_amp_for_30deg(geometry, power=1.0)

    rows = []
    for dist in DISTANCES:
        case_name = "R%.1f_x%.0f_y%.0f" % (dist, THETA_X, THETA_Y)
        print("\n=== %s ===" % case_name)
        case_data = build_case_data(geometry, dist, THETA_X, THETA_Y)

        # Uniform baseline
        uniform_metrics = metrics_for_case(np.ones(n_elements), np.zeros(n_elements), case_data)

        # Heuristic baseline
        heuristic_metrics = metrics_for_case(heuristic_amp, np.zeros(n_elements), case_data)

        # GA optimization
        print("Running GA amp+phase ...")
        best_amp, best_phase, best_hist, mean_hist = run_ga_amp_phase(
            case_data, n_elements,
            seed=20260511 + int(round(dist * 1000))
        )
        opt_metrics = metrics_for_case(best_amp, best_phase, case_data)

        print("  uniform  SL=%.2f p99=%.2f" % (uniform_metrics["max_sidelobe_db"], uniform_metrics["p99_sidelobe_db"]))
        print("  heuristic SL=%.2f p99=%.2f" % (heuristic_metrics["max_sidelobe_db"], heuristic_metrics["p99_sidelobe_db"]))
        print("  GA opt   SL=%.2f p99=%.2f" % (opt_metrics["max_sidelobe_db"], opt_metrics["p99_sidelobe_db"]))

        rows.append({
            "case": case_name,
            "distance_m": dist,
            "uniform": uniform_metrics,
            "uniform_db": uniform_metrics.pop("db"),
            "heuristic": heuristic_metrics,
            "heuristic_db": heuristic_metrics.pop("db"),
            "heuristic_amp": heuristic_amp,
            "opt": opt_metrics,
            "opt_db": opt_metrics.pop("db"),
            "opt_amp": best_amp,
            "opt_phase": best_phase,
            "txs": case_data["txs"],
            "tys": case_data["tys"],
        })

    # Plot all cases
    patterns_png = os.path.join(OUT_DIR, "multi_unit_30deg_patterns.png")
    plot_comparison(rows, patterns_png)

    # Plot weights for first case
    weights_png = os.path.join(OUT_DIR, "multi_unit_30deg_weights.png")
    plot_weights_comparison(rows, geometry, weights_png)

    # Write summary
    summary_path = os.path.join(OUT_DIR, "multi_unit_30deg_summary.txt")
    write_summary(summary_path, rows, geometry)

    print("\nOutputs:")
    for p in [patterns_png, weights_png, summary_path]:
        print("  " + p)


if __name__ == "__main__":
    main()
