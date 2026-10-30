# STL operator parameter study

Measures how the SubtreeLeap (STL) operator's mixing efficiency depends on its
acceptance-probability target and its weight, on two datasets (small-2 and tall-1).

## Study design

Three arms, each with the single-STL default as a shared control:

| Arm | What varies | Fixed |
|-----|-------------|-------|
| `target_sweep` | STL `targetAcceptance` ∈ {0.05, 0.10, 0.15, **0.234**, 0.35, 0.50} | default weight (nTaxa) |
| `weight_sweep` | STL weight factor ∈ {0.1, 0.25, 0.5, **1.0**, 2.0, 4.0} × nTaxa | target 0.234 |
| `two_stl` | 2nd STL target ∈ {0.03, 0.05, 0.10, 0.15, 0.20} × factor ∈ {0.1, 0.3} × nTaxa | 1st STL at default |

5 seeds for the default control, 3 seeds for each non-default condition.
Datasets: **small-2** (165 taxa, 40 M states) and **tall-1** (1441 taxa, 60 M states),
both BEAGLE and MAT likelihoods.

Estimated scope: **260 distinct runs**, **575 core-hours**, ~48 hours wall time at 12 workers.

## Reusing stl2_study runs

The `two_stl` arm's run IDs use the same format as `stl2_study` for the overlapping conditions
(`stl2t{target}w{factor}`).  Before starting, hard-link completed stl2_study runs to save time:

```sh
cd stl_study
mkdir -p runs
for src in ../stl2_study/runs/*/; do
    id=$(basename "$src")
    dst="runs/$id"
    [ -d "$dst" ] && continue
    cp -rl "$src" "$dst"   # hard-link the whole directory tree
done
```

## Running

```sh
cd stl_study

# Generate conditions (all three arms at once)
python3 design.py all

# Optional: hard-link stl2_study runs (see above) before starting

# Run (add --arm target_sweep to start with just one arm)
python3 run_grid.py

# Analyse when enough runs are done
python3 analyse.py all
```

To run just one arm: `python3 run_grid.py --arm target_sweep`.

## Outputs

- `analysis/summary.csv` — per cell: ESS/h ratio vs control (bootstrap CI), stuck-run count
- `analysis/target_sweep_*.png` — line plots: STL acceptance rate vs efficiency
- `analysis/weight_sweep_*.png` — line plots: STL weight factor vs efficiency
- `analysis/two_stl_*.png` — heatmaps: (stl2 target × stl2 weight factor) vs efficiency
