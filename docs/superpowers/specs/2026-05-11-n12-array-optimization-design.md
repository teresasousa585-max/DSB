# N12 Ultrasonic Array Multi-Agent Optimization Design

## Objective
1. Find better 12-transducer layouts within 10x10 cm with >=17.2 mm center spacing.
2. Jointly optimize phase and amplitude for best directivity within +/-15 deg.

## Architecture
- Orchestrator launches 3 agents in parallel via multiprocessing.
- Each agent writes checkpoint.json periodically.
- Orchestrator monitors progress and terminates stuck agents.

## Agent A: Free-Form Layout GA (`agent_layout_ga.py`)
- Encoding: 24 real variables (x, y for 12 elements), bounded to +/-50 mm.
- Constraint: collision penalty for center spacing < 17.2 mm.
- Algorithm: real-coded GA with elitism, SBX crossover, polynomial mutation.
- Fitness: worst max-sidelobe + 0.35*p99-sidelobe over 9 steering cases:
  (0,0), (+/-5,0), (+/-10,0), (+/-15,0), (0,+/-15).
- Fast mode: use array factor only; verify top-10 with full piston pattern.
- Termination: 300 generations or 50 generations without improvement.

## Agent B: Complex Weight ES (`agent_weight_es.py`)
- Encoding: 24 real variables (12 amplitudes 0..1, 12 phase offsets -pi..pi).
- Actual phase = steering_phase + offset.
- Algorithm: CMA-ES (prefer) or Differential Evolution.
- Fitness: smooth minimax sidelobe over the same 9 steering cases.
- Constraint: mainlobe gain >= floor (penalty if violated).
- Runs on current best layout; switches if Agent A finds a better one.
- Termination: 50k function evaluations or CMA-ES convergence.

## Agent C: Parametric Refinement (`agent_parametric_search.py`)
- Fine-grained search over ring-5-7, sunflower, rect-3-4 families.
- Uses denser grid + random perturbation + local Nelder-Mead polish.
- Baseline comparison for Agent A.

## Synchronization
- Shared directory: `analysis_outputs/agent_checkpoints/`
- Each agent writes `{agent_name}_checkpoint.json` every 30 seconds.
- Orchestrator reads checkpoints and maintains `global_best.json`.
- Stuck threshold: 5 minutes without checkpoint update -> terminate.

## Deliverables
- `n12_layout_ga_best.csv` + plots
- `n12_weight_es_best.csv` + plots
- `n12_parametric_best.csv` + plots
- `n12_multiagent_report.txt` comparing all agents
