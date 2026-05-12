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

C = 343.0
FREQ = 40000.0
LAMBDA = C / FREQ
K = 2.0 * np.pi / LAMBDA
ELEMENT_DIAMETER_MM = 16.0
ELEMENT_RADIUS_M = ELEMENT_DIAMETER_MM * 0.5e-3
MIN_CENTER_SPACING_MM = 17.2
DB_FLOOR = -55.0

# GA parameters
POP_SIZE = 40
N_GEN = 80
STAGNATION_LIMIT = 25
SBX_ETA = 15.0
PM_ETA = 20.0
SBX_PROB = 0.9
PM_PROB = 0.15
ELITE_COUNT = 2
PENALTY_PER_VIOLATION_DB = 100.0

# Steering cases for fitness evaluation
STEERING_CASES = [
    (0, 0),
    (5, 0),
    (-5, 0),
    (10, 0),
    (-10, 0),
    (15, 0),
    (-15, 0),
    (0, 15),
    (0, -15),
]


def ensure_out_dir():
    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)


def min_spacing_mm(points):
    pts = np.asarray(points, dtype=float)
    best = 1e9
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            best = min(best, float(np.linalg.norm(pts[i] - pts[j])))
    return best


def count_violations(points, min_dist=MIN_CENTER_SPACING_MM):
    pts = np.asarray(points, dtype=float)
    count = 0
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if float(np.linalg.norm(pts[i] - pts[j])) < min_dist:
                count += 1
    return count


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


def response_db(points_mm, steer_x_deg=0.0, steer_y_deg=0.0):
    axis, tx, ty, theta_deg, dirs, elem, mask = base.GRID
    pts_m = np.column_stack([
        np.asarray(points_mm)[:, 0] * 1e-3,
        np.asarray(points_mm)[:, 1] * 1e-3,
        np.zeros(len(points_mm)),
    ])
    w = base.steering_vector(points_mm, steer_x_deg, steer_y_deg)
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


def evaluate_candidate(points_mm):
    pts = centered(points_mm)
    if len(pts) != 12:
        raise ValueError("layout must have 12 points")
    spacing = min_spacing_mm(pts)
    violations = count_violations(pts)
    axis, tx, ty, theta_deg, dirs, elem, mask = base.GRID
    max_sidelobe = -1e9
    p99_sidelobe = -1e9
    for sx, sy in STEERING_CASES:
        db = response_db(pts, sx, sy)
        sep = angular_separation_deg(tx, ty, sx, sy)
        outside = (sep > 11.0) & mask
        sl = float(np.nanmax(db[outside]))
        p99 = float(np.nanpercentile(db[outside], 99.0))
        max_sidelobe = max(max_sidelobe, sl)
        p99_sidelobe = max(p99_sidelobe, p99)
    fitness = max_sidelobe + 0.35 * p99_sidelobe
    if violations > 0:
        fitness += PENALTY_PER_VIOLATION_DB * violations
    return {
        "points": pts,
        "fitness": fitness,
        "max_sidelobe_db": max_sidelobe,
        "p99_sidelobe_db": p99_sidelobe,
        "min_spacing_mm": spacing,
        "violations": violations,
    }


def decode(individual):
    xs = individual[:12]
    ys = individual[12:]
    return np.column_stack([xs, ys])


def encode(points_mm):
    pts = np.asarray(points_mm, dtype=float)
    return np.concatenate([pts[:, 0], pts[:, 1]]).tolist()


def create_individual():
    return [random.uniform(-50.0, 50.0) for _ in range(24)]


def sbx_crossover(p1, p2, eta=SBX_ETA):
    c1 = p1.copy()
    c2 = p2.copy()
    for i in range(24):
        if random.random() <= 0.5:
            if abs(p1[i] - p2[i]) > 1e-14:
                if p1[i] < p2[i]:
                    y1, y2 = p1[i], p2[i]
                else:
                    y1, y2 = p2[i], p1[i]
                beta = 1.0 + (2.0 * (y1 - (-50.0)) / (y2 - y1))
                alpha = 2.0 - beta ** (-(eta + 1.0))
                rand = random.random()
                if rand <= 1.0 / alpha:
                    beta_q = (rand * alpha) ** (1.0 / (eta + 1.0))
                else:
                    beta_q = (1.0 / (2.0 - rand * alpha)) ** (1.0 / (eta + 1.0))
                c1[i] = 0.5 * ((y1 + y2) - beta_q * (y2 - y1))
                beta = 1.0 + (2.0 * (50.0 - y2) / (y2 - y1))
                alpha = 2.0 - beta ** (-(eta + 1.0))
                rand = random.random()
                if rand <= 1.0 / alpha:
                    beta_q = (rand * alpha) ** (1.0 / (eta + 1.0))
                else:
                    beta_q = (1.0 / (2.0 - rand * alpha)) ** (1.0 / (eta + 1.0))
                c2[i] = 0.5 * ((y1 + y2) + beta_q * (y2 - y1))
                c1[i] = np.clip(c1[i], -50.0, 50.0)
                c2[i] = np.clip(c2[i], -50.0, 50.0)
    return c1, c2


def polynomial_mutation(individual, eta=PM_ETA):
    mutant = individual.copy()
    for i in range(24):
        if random.random() < PM_PROB:
            x = individual[i]
            delta1 = (x - (-50.0)) / (50.0 - (-50.0))
            delta2 = (50.0 - x) / (50.0 - (-50.0))
            rand = random.random()
            mut_pow = 1.0 / (eta + 1.0)
            if rand <= 0.5:
                xy = 1.0 - delta1
                val = 2.0 * rand + (1.0 - 2.0 * rand) * (xy ** (eta + 1.0))
                delta_q = val ** mut_pow - 1.0
            else:
                xy = 1.0 - delta2
                val = 2.0 * (1.0 - rand) + 2.0 * (rand - 0.5) * (xy ** (eta + 1.0))
                delta_q = 1.0 - val ** mut_pow
            mutant[i] = x + delta_q * (50.0 - (-50.0))
            mutant[i] = np.clip(mutant[i], -50.0, 50.0)
    return mutant


def tournament_select(pop, fitnesses, k=2):
    best = random.randrange(len(pop))
    for _ in range(k - 1):
        contender = random.randrange(len(pop))
        if fitnesses[contender] < fitnesses[best]:
            best = contender
    return pop[best]


def run_ga(seed=20260511):
    random.seed(seed)
    np.random.seed(seed)
    population = [create_individual() for _ in range(POP_SIZE)]
    fitnesses = [None] * POP_SIZE
    evaluated = [False] * POP_SIZE

    best_fitness_history = []
    mean_fitness_history = []

    best_ever = None
    best_fitness_ever = float('inf')
    stagnation_count = 0

    start_time = time.time()

    for gen in range(N_GEN):
        # Evaluate unevaluated individuals
        for i in range(POP_SIZE):
            if not evaluated[i]:
                pts = decode(population[i])
                res = evaluate_candidate(pts)
                fitnesses[i] = res["fitness"]
                evaluated[i] = True

        # Elitism: sort by fitness
        sorted_indices = np.argsort(fitnesses)
        new_population = [population[idx].copy() for idx in sorted_indices[:ELITE_COUNT]]
        new_fitnesses = [fitnesses[idx] for idx in sorted_indices[:ELITE_COUNT]]
        new_evaluated = [True] * ELITE_COUNT

        # Update best ever
        if fitnesses[sorted_indices[0]] < best_fitness_ever:
            best_fitness_ever = fitnesses[sorted_indices[0]]
            best_ever = decode(population[sorted_indices[0]]).copy()
            stagnation_count = 0
        else:
            stagnation_count += 1

        best_fitness_history.append(best_fitness_ever)
        mean_fitness_history.append(np.mean(fitnesses))

        print("Gen %3d | best=%.4f | mean=%.4f | stag=%d" % (
            gen, best_fitness_ever, np.mean(fitnesses), stagnation_count
        ))

        if stagnation_count >= STAGNATION_LIMIT:
            print("Stagnation reached after %d generations." % gen)
            break

        # Generate offspring
        while len(new_population) < POP_SIZE:
            p1 = tournament_select(population, fitnesses)
            p2 = tournament_select(population, fitnesses)
            if random.random() < SBX_PROB:
                c1, c2 = sbx_crossover(p1, p2)
            else:
                c1, c2 = p1.copy(), p2.copy()
            c1 = polynomial_mutation(c1)
            c2 = polynomial_mutation(c2)
            new_population.append(c1)
            new_fitnesses.append(None)
            new_evaluated.append(False)
            if len(new_population) < POP_SIZE:
                new_population.append(c2)
                new_fitnesses.append(None)
                new_evaluated.append(False)

        population = new_population
        fitnesses = new_fitnesses
        evaluated = new_evaluated

    elapsed = time.time() - start_time
    print("GA finished in %.1f seconds." % elapsed)

    # Final evaluation of best
    best_res = evaluate_candidate(best_ever)
    return best_res, best_fitness_history, mean_fitness_history


def write_layout_csv(path, points_mm):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Element_ID", "X_mm", "Y_mm"])
        for i, (x, y) in enumerate(points_mm):
            writer.writerow(["U%02d" % i, "%.6f" % x, "%.6f" % y])


def write_history_csv(path, best_hist, mean_hist):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["generation", "best_fitness", "mean_fitness"])
        for i, (b, m) in enumerate(zip(best_hist, mean_hist)):
            writer.writerow([i, "%.6f" % b, "%.6f" % m])


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
    ax.set_title("GA-optimized 12-element layout")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    ensure_out_dir()
    best_res, best_hist, mean_hist = run_ga()
    pts = best_res["points"]
    ap = aperture_mm(pts)
    unit_width = max(ap["width_centers_mm"], ap["height_centers_mm"]) + ELEMENT_DIAMETER_MM + 4.0

    layout_csv = os.path.join(OUT_DIR, "agent_layout_ga_best.csv")
    write_layout_csv(layout_csv, pts)
    history_csv = os.path.join(OUT_DIR, "agent_layout_ga_history.csv")
    write_history_csv(history_csv, best_hist, mean_hist)
    layout_png = os.path.join(OUT_DIR, "agent_layout_ga_layout.png")
    plot_layout(pts, layout_png)

    print("\n=== GA Best Layout Summary ===")
    print("Best fitness: %.4f" % best_res["fitness"])
    print("Worst max sidelobe: %.2f dB" % best_res["max_sidelobe_db"])
    print("Worst p99 sidelobe: %.2f dB" % best_res["p99_sidelobe_db"])
    print("Min center spacing: %.3f mm" % best_res["min_spacing_mm"])
    print("Approx unit width: %.1f mm" % unit_width)
    print("Violations: %d" % best_res["violations"])
    print("Saved layout to: %s" % layout_csv)
    print("Saved history to: %s" % history_csv)
    print("Saved plot to: %s" % layout_png)


if __name__ == "__main__":
    main()
