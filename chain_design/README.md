# Chain-design study: how best to use 16 cores on tall-1 (BEAGLE)

This study compares ways of spending the 16 performance cores of Machine 1 (Mac Studio M2
Ultra, 16 performance + 8 efficiency cores) on a single analysis:
- **one chain** with several threads;
- **several independent chains** sharing the threads;
- **parallel-tempering sets.**

Each configuration is scored by how soon, in wall-clock time, it reaches a good effective
sample size (ESS) for continuous parameters and for the tree topology, with **burn-in
counted**. Burn-in is also measured in its own right: its time is a cost that every extra
chain pays again.

The dataset is tall-1 (1441 taxa, 987 sites, `benchmark-tall-1-min.xml`). It is the
hardest benchmark to mix, and single chains sometimes get stuck on a wrong clade. The
likelihood is BEAGLE only.

**Requires** a BEAST build of the `multichain_runner` branch (`-chains`, `-mc3_*`,
`-thread_scheme`), checked out next to this repository as `../BEAST`. Build it with
`ant dist`; the first run freezes `build/dist/beast.jar` into `beast-frozen.jar`.

## Configurations

All configurations use `-thread_scheme serial -threads 0` (one partition, so all of a chain's
threads go to BEAGLE) and an explicit `-beagle_threads`. None uses more than 16 threads.

| Id | Configuration | Threads per chain |
|---|---|---|
| I1 | 1 chain | t* (the best single-chain thread count, from pilot P1) |
| I2 | 2 chains | 8 |
| I4 | 4 chains | 4 |
| I8 | 8 chains | 2 |
| I16 | 16 chains | 1 |
| T1 | 1 tempering set of 1 cold + 3 heated chains | 4 |
| T2 | 2 tempering sets | 2 |
| T4 | 4 tempering sets | 1 |

**Tempering:**
- non-reversible swap scheme (`-mc3_scheme deterministic_multiple`);
- a geometric ladder βᵢ = rⁱ, with the ratio r chosen in pilot P2;
- swaps attempted so that each block between swaps takes about 0.5 s, so that waiting at the
  swap barrier is a small cost.

Only the cold chain of each set is logged and scored.

## Design

1. **Pilots** (about 5 h; see `study.json` → `pilot`):
   - **P1, throughput:** states per second per chain for chains × threads (K ∈ {1,2,4,8,16},
     t ∈ {1,2,4,8,16}, K·t ≤ 16), plus T1, T2 and T4. Five minutes each, two repeats. This
     gives t* (the fewest threads within 3% of the fastest single chain) and the speeds that
     set every main run's chain length and log interval.
   - **P2, tempering:** T1 with ladder ratios 0.99, 0.98 and 0.96 (the joint's posterior SD is
     about 59, so neighbouring rungs need Δβ around 0.017). The hottest ladder whose
     neighbouring chains still swap at least 20% of the time is used. The pilot's time is a
     cost of tempering.
   - **P3, JVM:** 16 chains in one JVM (ParallelGC and G1) against 16 separate JVMs. This
     checks that running chains in one JVM costs nothing and that GC pauses are small.
2. **Main runs:**
   - 8 configurations × 4 seeds × 4 h = 128 h;
   - one run at a time with the whole machine;
   - configurations in random order within each seed block, so drift in the machine's state
     is spread over them;
   - each run's chain length is set from P1 so that it ends at about 4 h. A run that hasn't
     finished by 5 h is stopped. BEAST's logs are flushed line by line, so a stopped run's
     logs are complete.
3. **Logging:** every file log gets a `<timeLogger/>` column, `secondsElapsed`, so ESS can be
   measured against wall-clock time. `run_grid.py` also records when each chain's first
   sample appears, which is the start-up cost: parsing K models and creating BEAGLE.

## Chain length needed (Phase 0, from earlier runs)

13 default single-threaded tall-1 BEAGLE runs of 60M states from the second-SubtreeLeap study
(`stl2_study/runs/tall-1-min_60M_default_s*`), at about 4.7 minutes per million states:

| Quantity | Median (range) |
|---|---|
| Burn-in: joint reaches the reference's 1% quantile | 1.1M states (0.8–1.6M), 3–8 minutes |
| min_cont ESS per million states | 1.55 (0.9–2.3) |
| Worst-clade ESS per million states (50 most uncertain clades) | 0.19 (0.02–0.33) |
| Stuck (worst-clade ESS below 2) | 2 of 13 (seeds 3 and 8) |

**Total post-burn-in states needed:**
- continuous ESS 200: about 130M;
- clade ESS 100: about 500M, about 39 hours of one single-threaded chain.

So 16 single-threaded chains need about 31M states each, roughly 2.5 h plus burn-in, while
one chain needs about 10 h even if threading makes it about 4× faster. A 4 h run separates
the configurations. The few-chain configurations are expected to run out of time on the clade
target; for those, the time to target is extrapolated from the fitted ESS rate
(`predicted_time_to_target`).

## Metrics (`analyse.py`)

The reference posterior is `operator_study/reference/tall-1-min.json.gz` (four 300M-state
runs, pairwise ASDSF 0.006). Its 1% joint quantile, means and SDs are in
`study.json` → `reference_stats`.

**Per chain** (`analysis/main/chains.csv`; for tempered runs, each set's cold chain):
- **Burn-in:** the first sample whose joint density reaches the reference's 1% quantile, in
  states and in minutes from process start (start-up included). `configs.csv` also gives the
  burn-in in core-minutes (minutes × chains × threads), which is the cost of a configuration's
  burn-ins.
- **Speed:** states per second after the first minute.
- **Post-burn-in ESS per million states and per hour:**
  - min_cont: the minimum over joint, prior, likelihood, rootHeight, treeLength, popSize and
    clock.rate;
  - clade: the minimum over the indicators of the reference's 50 most uncertain clades;
  - RF: the RF distance to the reference tree.
- **Stuck:** a clade more than 0.3 away from the reference frequency in the second half of the
  chain. **Escaped:** stuck in the first half but not the second.

**Per run, as curves over wall-clock time** (`analysis/main/curves.csv`, every 5 minutes; chains
truncated at t; each chain's burn-in removed):
- **Continuous ESS:** rank-normalised multi-chain bulk ESS (Vehtari et al. 2021), minimum over
  the min_cont parameters. Chains that disagree lower it, so a stuck chain is not rewarded.
- **Clade ESS:** multi-chain ESS of the 50 clade indicators (worst clade).
- **Summed per-chain ESS** (Tracer's estimator): raw throughput, for comparison.
- **Accuracy:** every uncertain clade within 0.05 of the reference frequency, and every
  min_cont mean within 0.1 reference SD.
- **Success:** accurate, continuous ESS ≥ 200 and clade ESS ≥ 100.

**Per run** (`analysis/main/runs.csv`):
- time to target, or "censored" if it was not reached;
- ESS rate and effective burn-in, from a straight-line fit ESS(t) = r·(t − b);
- the predicted time to target;
- stuck and escaped chains;
- the practical burn-in: the split-R-hat rule on a 2% grid, applied to a single chain's halves
  or across chains, for comparison with the reference-based burn-in;
- GC pause fraction and machine load.

**Per configuration** (`analysis/main/configs.csv`):
- Kaplan–Meier median time to target (runs that did not reach it are censored);
- median predicted time to target;
- median ESS per hour and its coefficient of variation across seeds;
- burn-in minutes, and burn-in core-minutes;
- runs with a stuck chain, and runs accurate at the end;
- **diagnostic power:** P(K chains include both a stuck and a good chain), using the stuck
  rate over all chains.

Plots: `analysis/main/ess_curves.png` and `analysis/main/burnin.png`.

## Workflow (on Machine 1)

```sh
cd chain_design
# 0. build BEAST: in ../BEAST, `git checkout multichain_runner && ant dist`

# 1. pilots P1 (throughput) and P3 (JVM), about 3 h
python3 design.py throughput
python3 run_grid.py --check                     # load must be low; nothing else running
python3 run_grid.py --phase throughput --phase jvm
python3 analyse.py pilot                        # -> pilot/summary.json: speeds, t*, JVM/GC

# 2. pilot P2 (tempering ladder), about 2 h
python3 design.py tempering
python3 run_grid.py --phase tempering
python3 analyse.py pilot                        # -> chosen ladder ratio
#    if the JVM pilot shows a real cost or GC pauses > 2%, change java_args in study.json

# 3. main runs, about 128 h
python3 design.py main
python3 run_grid.py --check
python3 run_grid.py --phase main
python3 analyse.py main                         # -> analysis/main/; can be run at any time on the finished runs
```

Ctrl-C (or `kill` of `run_grid.py`) stops the current run, which is redone next time; finished
runs are kept. Each run locks its directory, so starting a second grid by mistake skips runs
that are in progress.

**Stall watchdog.** If a run's chain logs stop growing for `stall_minutes` (10), `run_grid.py`
saves a Java thread dump (`stall_*.txt`), stops the run, and moves its directory to
`runs/<run_id>.stalled_<time>` so that the run is redone. During testing, one tempered smoke run
(1×(1+3), 4 BEAGLE threads per chain) stopped logging after 2.5 minutes. Six identical runs did
not stall, so the cause is unknown; a later smoke stall turned out to be the laptop going to
sleep, and that one may have been too. If it happens on Machine 1, the thread dump will show
where the chains are blocked.

**Sleep.** On macOS `run_grid.py` runs `caffeinate -i` for as long as it runs, so the Mac does
not idle-sleep. A sleep it cannot prevent (a closed lid, a thermal emergency) is detected from
the gap between the wall clock and the monotonic clock: a run during which the machine slept for
more than a minute is moved to `runs/<run_id>.slept_<time>` and redone, because its timings are
meaningless.

After each run is analysed, its trees are condensed into per-sample traces
(`*.condensed.json.gz`: the joint and min_cont parameters, the RF distance to the reference
tree, and the 50 clade indicators) and the tree file is gzipped.

For a quick end-to-end test, point `CHAIN_DESIGN_CONFIG` at a copy of `study.json` with a
`work_dir`, short times and a few configurations.

---

# The big-2 study on darwin (CPU cores and GPUs)

The same question for a different machine, with `study_darwin.json` (work directory
`darwin/`). Every command below needs `CHAIN_DESIGN_CONFIG=study_darwin.json`.

**darwin:**
- 2 × Xeon E5-2687W v3: 20 physical cores (40 hyper-threads) on 2 NUMA nodes;
- 2 × Tesla P100 (16 GB);
- Linux, with CUDA BEAGLE.

**Dataset:** big-2 (1610 taxa, 7927 patterns, `benchmark-big-2-min.xml`). It costs about 4× tall-1 per state; its mixing has not been measured, and it has no reference posterior.

## Configurations

| Id | Device | Configuration |
|---|---|---|
| C1 | CPU | 1 chain × t* threads (t* from the pilot, up to 40) |
| C2 | CPU | 2 × 10 threads (one JVM per socket if the pilot shows binding helps by at least 10%) |
| C4, C10, C20 | CPU | 4 × 5, 10 × 2, 20 × 1 |
| C40 | CPU | 40 × 1 on the hyper-threads, **only if** the pilot shows 40 × 1 gives at least 10% more than 20 × 1 |
| CT | CPU | 1 tempering set of 1 + 3 chains, 5 threads each |
| G1 | GPU | 1 chain on one GPU |
| G2 | GPU | 2 chains, one per GPU |
| G4 | GPU | 4 chains, two per GPU, **only if** the pilot shows that gives at least 10% more than G2 |
| GT | GPU | 1 tempering set of 1 + 3 chains, two per GPU |
| M | GPU + CPU | 2 GPU chains plus 4 CPU chains × 4 threads |

**How chains are placed:**
- GPU chains use `-beagle_GPU -beagle_order` with resources 1 and 2 alternating (CPU chains in M get resource 0). Check the resource numbers with `beast -beagle_info` and set `machine.gpu_resources` if they differ.
- BEAGLE instances are created in chain order, but only chain 1's BEAGLE messages are printed. So the placement is confirmed from `nvidia-smi`, which `run_grid.py` samples into each GPU run's `gpu.csv`.

## Reference

There is no big-2 reference, so the study uses its own pooled posterior (`"reference": "pooled"`):
- **Pool:** the second halves of all chains of the main and mixing runs.
- **Scoring:** each run is scored against the pool with its own samples removed.
- **What the pool gives:** the 50 most uncertain clades, a majority-rule focal tree, the joint's 1% quantile (for burn-in), and parameter means and SDs.
- **Convergence check:** `analyse.py main` prints the ASDSF between the pools of seeds {1, 2} and {3}, and warns if it is 0.01 or more. Then the pool is not converged, and accuracy cannot be judged.

## Workflow (on darwin)

```sh
cd chain_design
export CHAIN_DESIGN_CONFIG=study_darwin.json
# 0. build BEAST (multichain_runner) in ../BEAST; check the GPUs and the BEAGLE path
java -Djava.library.path=/usr/local/lib -cp ../../BEAST/build/dist/beast.jar dr.app.beast.BeastMain -beagle_info
nvidia-smi -L
#    adjust java_args (-Djava.library.path) and machine.gpu_resources in study_darwin.json if needed
#    use tmux or nohup for the long phases (there is no caffeinate on Linux)

# 1. throughput pilot (about 5 h): CPU chains x threads, hyper-threads, NUMA binding, GPU sharing, M, CT, GT
python3 design.py throughput && python3 run_grid.py --check
python3 run_grid.py --phase throughput
python3 analyse.py pilot               # t*, totals, hyper-threading/GPU-sharing/NUMA gains, GPU %

# 2. mixing pilot (24 h on G2): burn-in, ESS per million states, the joint's SD
python3 design.py mixing && python3 run_grid.py --phase mixing
python3 analyse.py pilot               # joint SD -> the tempering ladders
python3 analyse.py mixing              # mixing of the pilot against its own pool -> choose main.hours

# 3. tempering ladders (about 4 h): GT and CT at ratios 1 - f / SD, f in {0.5, 1, 2}
python3 design.py tempering && python3 run_grid.py --phase tempering
python3 analyse.py pilot               # the chosen ladder ratio

# 4. main grid (10-12 configurations x 3 seeds x main.hours; about 250 h at 8 h)
python3 design.py main && python3 run_grid.py --check
python3 run_grid.py --phase main
python3 analyse.py main                # -> darwin/analysis/main/
```

After the mixing pilot, set `main.hours` so that the fastest configurations (on the pilot's rates) reach the targets with room to spare, keeping the main grid within about 250 h.
