# Operator study

Tests two BEAST defaults using ESS per hour:

- **A: the 0.234 target acceptance**, especially for SubtreeLeap.
- **B: the relative operator weights.** The default gives SubtreeLeap weight = nTaxa, FHSPR = SubtreeLeap/10, and each scalar parameter 3.

**Study A now runs separately**, on the `subtreeleap_acceptance` branch in `acceptance_study/`, so that it can use a second machine. This machine runs B, the reference posteriors (which both studies use) and the pilot. `design.py` can still generate the A arms.

Datasets: small-1, small-2 and tall-1. BEAGLE (`-min`) is the primary arm. The MAT variants of small-2 and tall-1 form a secondary arm.

Requires the BEAST `operator_study` branch, which adds:
- per-operator timing in the `.ops` table;
- a working `-adaptation_target`;
- `-operator_weights KEY=FACTOR,...`.

## Arms

| Arm | What varies | How |
|---|---|---|
| reference | nothing (default settings, 10× chain length) | gives the reference clade frequencies and the focal tree |
| pilot | nothing (default) | sets chain lengths and the A1 size grid |
| A1 | SubtreeLeap `size`, with adaptation off | XML |
| A2 | SubtreeLeap `targetAcceptance` | XML |
| A3 | target for every adaptive operator | `-adaptation_target` |
| B1 | SubtreeLeap and FHSPR weight × m | `-operator_weights` |
| B2 | FHSPR weight as a fraction of SubtreeLeap's | `-operator_weights` |
| confirm | default against the best settings in `study.json` `"confirm"` | either |
| mat | A2 and B1 on the MAT XMLs | as above |

Run ids are built from the effective settings. A cell that is just the default (A2 at 0.234, B1 at ×1, B2 at 0.1, the pilot) therefore reuses the same default runs for each seed.

## Workflow

The pilot is done: its results are in `pilot/`, and the chain lengths and A1 sizes in `study.json` come from it. This machine runs the reference posteriors and the B arms in **one** `run_grid.py`, with 12 workers on the M2 Ultra's 16 performance cores. Don't start a second grid alongside it: that is what overloaded the first attempt.

```sh
cd operator_study
# build BEAST first (BEAST branch operator_study, `ant dist`); the first run freezes build/dist/beast.jar
python3 design.py reference B1 B2 mat        # -> conditions.csv (159 runs to run, below)
python3 run_grid.py --check --arm reference --arm B1 --arm B2 --arm mat-B1    # pre-flight, about 15 s
python3 run_grid.py --arm reference --arm B1 --arm B2 --arm mat-B1            # prints core-hours and ETA
python3 analyse.py reference                 # when the reference is done; check max pairwise ASDSF < 0.01
git add -f reference/*.json && git commit -m "Reference posteriors" && git push   # the acceptance study reads these
python3 analyse.py all                       # analysis/runs.csv, summary.csv, *.png
#    fill study.json "confirm" (e.g. {"sl_target": 0.1, "operator_weights": "..."})
python3 design.py confirm --append && python3 run_grid.py --arm confirm && python3 analyse.py summary
```

Runs start longest first, so the long reference runs go first. Ctrl-C, or `kill` of the `run_grid.py` process, stops every BEAST run it started. Interrupted runs are redone on restart, and finished runs are kept. Each run locks its directory while it runs, so a second grid started by mistake skips runs that are already in progress.

To start completely afresh, stop the grid, then delete `runs/`, `conditions.csv`, `beast-frozen.jar`, `reference/` and `analysis/`.

Run ids include the chain length (e.g. `small-1-min_30M_default_s1`). Changing `states` therefore creates new runs; it never silently reuses runs of the wrong length.

For a quick smoke test, point `OPERATOR_STUDY_CONFIG` at a copy of `study.json` with a different `work_dir` and short chains.

Per-dataset settings in `study.json`: `minutes_per_mstate` (pilot speeds, used for ordering and the ETA), `arms` (which arms include the dataset), `reps` ({arm: n}) and `reference_multiplier`. tall-1 uses them because it mixes slowly (pilot: about 36 ESS/h). It runs at 60M states with 3 replicates of A2, B1 and B2 only, and its reference is 5× long.

## Metrics (analyse.py)

- **ESS** uses the same estimator as Tracer/`loganalyser` (`TraceCorrelation`). It is computed after 10% burn-in, and reported per hour of BEAST run time (`eff_*`) and per million states (`mst_*`). The per-state figure doesn't depend on how loaded the machine was.
- **min_cont:** the minimum ESS over joint, prior, likelihood, rootHeight, treeLength, popSize, clock.rate and kappa.
- **rf:** the ESS of the Robinson–Foulds distance to the reference focal tree.
- **split_min:** the minimum ESS over the indicator traces of the 50 most uncertain reference clades (frequency between 0.05 and 0.95). A clade that is never or always present counts as an ESS of 1.
- **asdsf:** the mean absolute difference in clade frequency from the reference, over clades at ≥ 1%. `asdsf_x_sqrt_hours` is the matching accuracy-per-cost measure.
- **Operator diagnostics,** taken from the `.ops` table and the file log:
  - SubtreeLeap acceptance and adapted size;
  - μs per operation;
  - share of run time;
  - the rolling `topologyChanged` fraction and `mean_delta`.

`summary.csv` gives each cell's geometric-mean efficiency relative to the default runs on the same dataset and likelihood, with 95% bootstrap intervals.
