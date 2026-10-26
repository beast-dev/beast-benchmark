# Second SubtreeLeap study

Tests whether adding a second SubtreeLeap, with its own target acceptance, improves mixing. For each dataset it looks for the best combination of two things:
- the second SubtreeLeap's target acceptance;
- its weight, relative to the existing 0.234 SubtreeLeap and FHSPR.

Every run keeps the benchmark's operators as they are:
- SubtreeLeap: weight = number of taxa, adapting to 0.234;
- FHSPR: weight = taxa / 10;
- the parameter operators (plus history Gibbs, under MAT).

A grid run adds a second SubtreeLeap (XML id `stl2`). It copies the first, including the MAT resampler hook, and has:
- its own `targetAcceptance`;
- weight = factor × taxa.

The control is the unchanged default.

A low target makes the second move adapt to long leaps, for topology. A high target makes it adapt to short steps, for local mixing.

**Requires a BEAST build that includes operator ids in names** (beast-mcmc `operator_study`, commit `67b414ab1` or later). With it, the two moves report as `subtreeLeap:subtreeLeap(treeModel)` and `stl2:subtreeLeap(treeModel)`. `run_grid.py --check` fails if the build is older.

## Design

| | Levels |
|---|---|
| stl2 target acceptance | 0.02, 0.05, 0.10, 0.50 |
| stl2 weight (× number of taxa) | 0.03, 0.1, 0.3, 1.0 |
| Seeds | 3 per grid cell, 5 for the control |

| Dataset | Shape | Likelihoods | Chain length |
|---|---|---|---|
| small-2 | 165 taxa, 1698 sites | BEAGLE, MAT | 40M |
| tall-1 | 1441 taxa, 987 sites | BEAGLE (2×2 grid: targets 0.05/0.1, weights 0.1/0.3; 3 control seeds), MAT | 60M |
| wide-1 | 64 taxa, 185,578 sites; clock rate fixed | BEAGLE, MAT | set from the pilot |
| big-1 | 976 taxa, 29,409 sites | MAT only | set from the pilot |

That's 333 runs, about 355 core-hours with the placeholder chain lengths for wide-1 and big-1. Levels can be overridden per dataset and likelihood in `study.json` (`levels`, `reps`, `reps_default`); tall-1 BEAGLE uses this.

**Reference posteriors.** Each dataset is scored against one reference posterior. small-2 and tall-1 reuse the weights study's (`../operator_study/reference`, committed on this branch). wide-1 needs a new BEAGLE reference. big-1 needs a MAT one (`reference_likelihood: "mat"`); MAT has been validated against BEAGLE on small-2. `design.py reference` makes runs only for datasets that have no reference yet.

## Workflow

```sh
cd stl2_study
# build BEAST (operator_study, 67b414ab1 or later): `ant dist`; the first run freezes build/dist/beast.jar

# 1. pilot wide-1 and big-1 (their "states" are null): default runs, 20M states, 3 seeds
python3 design.py pilot && python3 run_grid.py --check && caffeinate -i python3 run_grid.py
python3 analyse.py runs
#    set datasets.wide-1/big-1 "states" (min ESS ~300 for the default) and "minutes_per_mstate"
#    from runs/*/result.json; then regenerate (no --append):

# 2. reference (wide-1 BEAGLE, big-1 MAT) and the grid, in ONE run_grid.py
python3 design.py reference grid
python3 run_grid.py --check
caffeinate -i python3 run_grid.py
python3 analyse.py reference                 # check max pairwise ASDSF < 0.01
python3 analyse.py all                       # analysis/summary.csv, grid_DATASET-LIKELIHOOD.png
```

To split the work across two machines, give each its own datasets, for example:
- `--dataset tall-1` on one machine (about 170 core-hours);
- `--dataset small-2 --dataset wide-1 --dataset big-1` on the other.

Each machine then holds the reference for its own datasets. Run `analyse.py` where the runs are, or bring the `runs/` together first.

Stopping, restarting, locking, and the run order (longest first) work as in `operator_study/`. `run_grid.py --check` runs a short BEAST chain of each kind of condition, including one with the second SubtreeLeap under each likelihood.

## Output

- **`analysis/summary.csv`:** one row per (dataset, likelihood, cell), with `stl2_target`, `stl2_factor` and:
  - ESS per hour (`eff_*`) and per million states (`mst_*`) relative to the control, with 95% bootstrap intervals, for the lowest continuous ESS (`min_cont`), the Robinson–Foulds distance (`rf`) and the worst uncertain clade (`split_min`);
  - both SubtreeLeaps' acceptance and step size (`sl_*`, `sl_stl2_*`), and each operator's share of run time.
- **`grid_DATASET-LIKELIHOOD.png`:** heatmaps of ESS per hour relative to the control, one panel per metric. Blue is better than without the second move, red is worse.

ESS per hour is the main measure here. The second move changes the cost per state, especially under MAT, where tree moves are cheap. That makes it important to run on an otherwise idle machine.
