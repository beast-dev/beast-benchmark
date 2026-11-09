# BEAST Benchmarking

This folder contains the benchmarking infrastructure for measuring and validating BEAST performance.

## Folder layout

```
beast-benchmark/
├── benchmark.py          # The benchmark runner script
├── benchmarks.py         # Finds the XMLs; shared by the scripts below
├── normalise.py          # Turns a BEAUti XML into a benchmark -min.xml
├── make_variants.py      # Derives the mat, gamma and skygrid variants
├── calibrate.py          # Measures speed and proposes chain lengths
├── apply_lengths.py      # Writes those chain lengths into the XMLs
├── apply_densities.py    # Writes each run's final joint density into its XML
├── run_all.sh            # The `min` set only
├── run_all_models.sh     # Every XML, every model variant
├── run_all_gpu.sh        # Every BEAGLE XML on the GPU (no mat)
├── README.md             # This file
├── Benchmarks/           # Benchmark XML files, organised by type
│   ├── Benchmark-big/
│   │   └── benchmark-big-1-min.xml
│   ├── Benchmark-small/
│   │   └── benchmark-small-1-min.xml
│   ├── Benchmark-tall/
│   ├── Benchmark-wide/
│   └── Benchmark-medium/
└── Results/              # Created automatically; one sub-folder per run
    └── {machine}_{type}_{number}_{model}_{timestamp}/
        ├── {prefix}benchmark-{type}-{number}-{model}.stdout
        ├── benchmark-{type}-{number}-{model}.log
        ├── benchmark-{type}-{number}-{model}.trees
        ├── {run_name}_report.md
        └── {run_name}_results.json
```

## Benchmark XML format

Each XML file begins with a structured header comment that records the dataset properties, expected MCMC settings, and the expected log-probability values for a canonical run at a fixed random seed.  The script parses this header to validate results.

```xml
<!-->
Benchmark: benchmark-small-1
Source data: YFV.fasta
Taxa: 71
Sites: 654
Site patterns: 280
Model: JC69, strict clock, coalescent constant population tree prior
Priors: rate ~ ApproximateReferencePrior, popSize ~ Gamma(0.001, 1000)
Operators:
  - SubtreeLeap on tree with weight 71.0
  ...
MCMC length: 10000000
Logging every: 10000
Seed: 1234
Initial Joint density: -29529.1683
Final Joint density: -6981.7346
<-->
```

The **Initial Joint density** is the log-probability of the initial state drawn at the given seed. The **Final Joint density** is the log-probability at the end of the chain. For a given seed, both values are fully deterministic and must match exactly across platforms if the implementation is correct. The benchmark runner reports `PASS` or `FAIL` for each.

### XML file naming convention

```
benchmark-{type}-{number}-{model}.xml
```

| Field    | Description |
| --- | --- |
| `type`   | Dataset category: `big`, `small`, `tall`, `wide`, `medium` |
| `number` | Dataset index within a type (1, 2, …) |
| `model`  | Model label — e.g. `min` for the minimum/baseline model configuration |

Model labels in use:

| `model` | Likelihood | Model |
| --- | --- | --- |
| `min` / `mat` | BEAGLE / MAT | the base model: strict clock, constant-size coalescent |
| `gamma` / `gamma-mat` | BEAGLE / MAT | the base model plus a 4-category discrete gamma site rate model |
| `skygrid` / `skygrid-mat` | BEAGLE / MAT | the base model with a Bayesian skygrid in place of the constant-size coalescent |

### Model variants

The `gamma` and `skygrid` XMLs are generated from each dataset's `min` and `mat` XMLs by `make_variants.py`:

```bash
python3 make_variants.py --beast "java -Djava.library.path=/opt/homebrew/lib -cp beast.jar dr.app.beast.BeastMain"
```

- **gamma:** the shape `alpha` has an Exponential(mean 0.5) prior and a scale operator (weight 3). The MAT
  version also samples each variable site's rate category with a `siteRateCategoryGibbsOperator` at 1% of
  the SubtreeLeap weight.
- **skygrid:** about one grid interval per 20 coalescent events (taxa − 1), at least 10 and at most 100.
  The cut-off is the 97.5% quantile of the root height from earlier BEAGLE and MAT runs, rounded up to 2
  significant figures (see `CUT_OFFS` in the script). The constant-size model is kept only to simulate the
  starting tree.

With `--beast`, the script runs each XML to state 0 to fill in `Initial Joint density`. It leaves out
`Final Joint density`, which needs a full run (step 4 of [Adding a new benchmark](#adding-a-new-benchmark)).

## Running a benchmark

```bash
python3 benchmark.py <type> <number> <model> [options]
```

**Example — run the small-1 minimum-model benchmark:**

```bash
python3 benchmark.py small 1 min
```

**Example — tag the run with a machine name:**

```bash
python3 benchmark.py small 1 min --machine bivalve
```

**Full usage:**

```
usage: benchmark.py [-h] [--seed SEED] [--beast BEAST] [--states STATES]
                    [--prefix PREFIX] [--benchmarks DIR] [--results DIR]
                    [--machine NAME] [--no-overwrite]
                    type number model

positional arguments:
  type              Benchmark type: big, small, tall, wide, medium
  number            Benchmark number (e.g. 1)
  model             Model label used in XML filename (e.g. min)

optional arguments:
  --seed SEED       Random seed (default: value from XML header, or 1234)
  --beast BEAST     Path to beast executable (default: 'beast' on PATH)
  --states STATES   Chain length to scale the progress bar against; does
                    NOT shorten the run (the length comes from the XML)
  --prefix PREFIX   Override the run directory name
                    (default: {machine_}{type}_{number}_{model}_{timestamp})
  --benchmarks DIR  Path to Benchmarks folder
                    (default: Benchmarking/Benchmarks)
  --results DIR     Path to results parent folder
                    (default: Benchmarking/Results)
  --machine NAME    Machine identifier — prefixed to the run folder and
                    output file names
  --no-overwrite    Do not pass -overwrite to BEAST
```

## Setting the chain lengths

The benchmarks measure wall-clock time and states per second, so they are most useful when
they all take about the same time — a run that takes 40 minutes tells you no more than one
that takes 10, and makes a full sweep cost a day. `calibrate.py` works out the chain length
that hits a target time on the machine you are on.

It does not run a chain to the end. It starts BEAST, waits for the screen log, times a
stretch of states, and kills it:

```bash
python3 calibrate.py --target 600 --beast /path/to/beast --machine machine1
```

That writes `calibration.csv` with, per benchmark, the startup time (JVM, XML parsing,
BEAGLE setup, first likelihood), the measured minutes per million states, and the chain
length those imply:

```
states = (target - startup) / seconds_per_state
```

rounded to two significant figures and a whole multiple of 1000. Then:

```bash
python3 apply_lengths.py --dry-run     # show what would change
python3 apply_lengths.py               # change it
```

### Correcting them from a completed sweep

`calibrate.py` only times the opening stretch of a chain, and that stretch is not
representative: on the M2 Ultra the BEAGLE benchmarks all ran faster over a full chain than
over their first minute, by 20% on the small alignments and by eight-fold on `wide-1`. The
MAT benchmarks, which do not use BEAGLE, came out within 20% throughout.

So treat calibration as the first guess and the sweep as the measurement. After a sweep,
rescale each length by how far its run landed from the target:

```bash
python3 apply_lengths.py --from-results --dry-run
python3 apply_lengths.py --from-results
```

It reads `results.csv`, takes the newest run of each benchmark, holds back the startup time
that no chain length can shorten, and leaves alone anything already within `--tolerance`
(20% by default) so a length only moves when it has to. One pass gets most benchmarks
inside the tolerance; repeat it after the next sweep for any that are still out.

## Recording the final joint density

Each XML header carries the joint density its chain ends on, and `benchmark.py` compares
every run against it: a change in BEAST that alters the chain shows up as a deviation. The
value belongs to one exact chain, so it is only written when the run's chain length and
seed still match the XML:

```bash
python3 apply_densities.py --dry-run   # show what would be written
python3 apply_densities.py             # write it
```

A benchmark whose length has moved since its last run is reported as `STALE` and left
without a density until the next sweep. CPU runs win over GPU ones, which end on a
different value because the arithmetic differs.

`--beast` takes the BEAST launcher script, or a `.jar`, or a bare name found on `PATH`; a
relative path is fine, as both runners make it absolute before running BEAST in the run
directory. A `.jar` is run with `java -jar`, which passes no JVM options, so BEAGLE has to be
somewhere Java already looks — otherwise BEAST exits with "No acceptable BEAGLE library
plugins found" and only the `mat` benchmarks run. Point `JDK_JAVA_OPTIONS` at it if needed:

```bash
JDK_JAVA_OPTIONS=-Djava.library.path=/usr/local/lib python3 calibrate.py --beast ./beast.jar
```

`calibrate.py` says `NO BEAGLE` against any non-`mat` benchmark that ran without it, and
lists them again at the end: such a run measures the Java likelihood, so its rate — and the
chain length taken from it — would not be comparable with the others.

`apply_lengths.py` sets all five places that have to agree — `chainLength` on `<mcmc>`,
`logEvery` on the screen, file and tree logs, and the `MCMC length` and `Logging every`
header lines — keeping 1000 samples per run. It also drops the `Final Joint density` line,
because that value belongs to the old chain length; run the sweep again to replace it.

**Chain lengths are machine-specific.** Calibrate on the machine the benchmarks are for.
The GPU runs are not separately calibrated: one XML holds one chain length, and it is set
for the CPU sweep, so a GPU run of the same XML will not take the target time.

## Output files

Each run creates a sub-folder inside `Results/`.  The folder name is:

```
{machine_}{type}_{number}_{model}_{timestamp}
```

where `{machine_}` is omitted if `--machine` is not supplied.

| File | Description |
| --- | --- |
| `{prefix}benchmark-{type}-{number}-{model}.stdout` | Full captured stdout from BEAST |
| `benchmark-{type}-{number}-{model}.log` | BEAST MCMC log (written by BEAST) |
| `benchmark-{type}-{number}-{model}.trees` | Tree log (written by BEAST) |
| `{run_name}_report.md` | Human-readable markdown report (see below) |
| `{run_name}_results.json` | Machine-readable JSON summary (see below) |

### Markdown report sections

- **Header table** — date, run directory, full command line
- **Dataset** — benchmark name, type, number, model, taxa, sites, site patterns
- **MCMC Settings** — chain length, logging frequency, seed
- **System** — macOS version, platform string, CPU, RAM, core counts, Python version
- **Software** — BEAST version, BEAGLE version, BEAGLE resource, BEAGLE thread count
- **Results / Joint Density** — table with observed vs expected initial and final joint log-density, deviation, and `✓ PASS` / `✗ FAIL` status
- **Results / Performance** — BEAST-reported runtime, wall-clock time, sample rate (min / million states)
- **Operator Analysis** — verbatim operator acceptance rate block from BEAST stdout

### JSON results file

```json
{
  "benchmark":        "benchmark-small-1-min",
  "type":             "small",
  "number":           1,
  "model":            "min",
  "taxa":             71,
  "sites":            654,
  "patterns":         280,
  "chain_length":     10000000,
  "sample_frequency": 10000,
  "seed":             1234,
  "timestamp":        "2026-05-14 12:00:00",
  "platform":         "macOS-26.3.1-arm64-arm-64bit",
  "cpu":              "Apple M3 Max",
  "ram_gb":           36,
  "cores_physical":   14,
  "cores_logical":    14,
  "beast_version":    "v10.5.0",
  "beagle_version":   "v4.0.1 (PRE-RELEASE)",
  "beagle_resource":  "CPU (arm64)",
  "beagle_threads":   14,
  "initial_joint":    -29529.1683,
  "final_joint":      -6981.7346,
  "runtime":          3.55,
  "wall_clock_time":  213.0,
  "sample_rate":      0.355
}
```

`runtime` is the BEAST-reported execution time in minutes. `wall_clock_time` is the wall-clock elapsed time in seconds. `sample_rate` is minutes per million MCMC states.

## Validation logic

Both the initial and final joint log-densities are **fully deterministic** for a given XML and seed. The runner checks both against the expected values recorded in the XML header:

- Deviation < 0.01 → **PASS**
- Deviation ≥ 0.01 → **FAIL** (exit code 2)

A failure indicates that either the BEAST implementation has changed in a way that alters the likelihood or MCMC trajectory, or that the benchmark XML header values need updating.

## Adding a new benchmark

The runners find the XMLs through `benchmarks.py`, so a new file is picked up by
`run_all_models.sh` and `run_all_gpu.sh` without editing them; a removed one disappears from
both. The name has to follow the convention above for that to work.

1. Create a folder `Benchmarks/Benchmark-{type}/`.
2. Place the XML file as `benchmark-{type}-{number}-{model}.xml`.
3. Add the structured header comment (see format above).
4. Set the chain length so the run takes about as long as the others:
   ```bash
   python3 calibrate.py {type} {number} {model} --target 600
   python3 apply_lengths.py
   ```
5. Run once with the canonical seed to obtain the initial and final joint densities:
   ```bash
   python3 benchmark.py {type} {number} {model} --seed 1234
   ```
6. Copy the reported joint densities into the XML header.
7. Re-run to confirm both values show `✓ PASS`.

## Adding a benchmark

A new benchmark arrives as BEAUti output and becomes six XMLs: the BEAGLE `min` one, the
MAT one, and the gamma and skygrid variants of each. Put the BEAUti file at
`Benchmarks/Benchmark-{type}/benchmark-{type}-{number}-min.xml`, then:

```bash
BEAST="java -Djava.library.path=/usr/local/lib -jar /path/to/beast.jar"

python3 normalise.py big 3 --source SARS-CoV-2.fasta --beast "$BEAST"
python3 make_variants.py --datasets big-3 --variants mat --beast "$BEAST"
python3 make_variants.py --datasets big-3 --variants gamma,skygrid --beast "$BEAST"
```

`normalise.py` rewrites the header, the `beast` version, the coalescent, the optimization
schedule and the output names, and drops the checkpoint log. `make_variants.py` needs the
dataset's skygrid cut-off in its `CUT_OFFS` first: the depth the grid reaches back to, which
should cover the root. Estimate it to start with and replace it with the 97.5% quantile of
the sampled root height once the benchmark has been run.

The chain length that `normalise.py` writes is a guess. Run the sweep, then
`apply_lengths.py --from-results` to set it from the measured wall clock, and
`apply_densities.py` to record the density each chain ends on.

Everything else picks the new XMLs up on its own: `benchmarks.py` finds them, so the
`run_all*.sh` scripts, `calibrate.py`, `apply_lengths.py` and `apply_densities.py` all
include them without being edited.
