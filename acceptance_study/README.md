# SubtreeLeap acceptance study

Tests BEAST's 0.234 target acceptance, especially for SubtreeLeap's adaptive step size, using ESS per hour. It is split off from the operator weights study (`operator_study/` on the `operator_study` branch) so it can run on a second machine.

The two studies use the same harness. This directory has its own copy (`design.py`, `make_variants.py`, `run_grid.py`, `analyse.py`, `common.py`), its own `study.json`, and its own `runs/`. Requires the BEAST `operator_study` branch, built with `ant dist`. That branch provides per-operator timing in `.ops` and a working `-adaptation_target`.

## Arms

| Arm | What varies | Levels | Seeds (small / tall-1) |
|---|---|---|---|
| A1 | SubtreeLeap `size` with adaptation off (pilot's adapted size × m) | m = 1/16 … 8 | 3 / 2 |
| A2 | SubtreeLeap `targetAcceptance` | 0.05, 0.10, 0.15, 0.234, 0.35, 0.50 | 5 / 3 |
| A3 | `-adaptation_target` for every adaptive operator | 0.10, 0.234, 0.44 | 5 / 3 |
| A4 | SubtreeLeap with a Cauchy distance kernel | 0.10, 0.234, 0.35 | 3 / – |
| mat | A2 on the MAT XMLs (small-2, tall-1) | as A2 | 3 |
| confirm | default against the best settings in `study.json` `"confirm"` | | 10 / 3 |

The chain lengths (small-1 30M, small-2 40M, tall-1 60M) and the A1 sizes (44.6, 0.416 and 0.36) come from the pilot in the weights study. They do not need re-piloting unless this machine is a different model. The design is 222 runs, about 315 core-hours, which is roughly 26 h on 12 workers. The MAT timings are guesses.

## Workflow

```sh
cd acceptance_study
python3 design.py A1 A2 A3 A4 mat            # -> conditions.csv (222 runs)
python3 run_grid.py                          # resumable, shuffled, 12 workers; freezes ../../BEAST/build/dist/beast.jar
```

The split-frequency metrics (rf, split_min, asdsf) need the reference posteriors. The weights study computes them, and `analyse.py` reads them from `../operator_study/reference/` (`"reference_dir"` in `study.json`). Once the weights machine has run `analyse.py reference` and committed `operator_study/reference/*.json` (with `git add -f`; the directory is ignored), fetch them here:

```sh
git fetch origin && git checkout origin/operator_study -- operator_study/reference
python3 analyse.py all                       # analysis/runs.csv, summary.csv, *.png
#    fill study.json "confirm" (e.g. {"best_target": {"sl_target": 0.1}})
python3 design.py confirm --append && python3 run_grid.py --arm confirm && python3 analyse.py summary
```

`analyse.py runs` works without the reference and gives the continuous-parameter ESS.

## What to look at

- **`A1_acceptance.png`:** efficiency against the realised SubtreeLeap acceptance, with adaptation off. This is the actual efficiency curve. If its peak is far from 0.234, the default target is wrong for this operator.
- **`A2.png`:** whether adapting to each target reaches the A1 optimum. Robbins–Monro adapts on min(1, r) rather than on the accept/reject outcome, so the two may differ.
- **`A3.png`:** whether the scalar operators also prefer a different target. One-dimensional theory gives about 0.44.
- **`sl_topologyChanged` in `summary.csv`:** the fraction of SubtreeLeap proposals that change topology. At the default target the pilot gave about 0.5 (small-1) to 0.8 (tall-1).
