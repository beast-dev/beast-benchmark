"""
Shared helpers for the chain-design study: configuration, run conditions, XML variants
and BEAST command lines. The log and tree parsers and the Tracer-equivalent ESS come from
stl2_study/common.py.
"""

import csv
import gzip
import json
import os
import re
from pathlib import Path

import importlib.util

STUDY_DIR = Path(__file__).resolve().parent
REPO_DIR = STUDY_DIR.parent

# the parsers and Tracer-equivalent ESS of the second SubtreeLeap study
_spec = importlib.util.spec_from_file_location("stl2_common", REPO_DIR / "stl2_study" / "common.py")
stl2_common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(stl2_common)
ess = stl2_common.ess
read_log = stl2_common.read_log
read_trees = stl2_common.read_trees
clades_from_json = stl2_common.clades_from_json

CONFIG_PATH = Path(os.environ.get("CHAIN_DESIGN_CONFIG", STUDY_DIR / "study.json"))

RUNS_DIR = STUDY_DIR / "runs"
PILOT_DIR = STUDY_DIR / "pilot"
ANALYSIS_DIR = STUDY_DIR / "analysis"
CONDITIONS_PATH = STUDY_DIR / "conditions.csv"
PILOT_SUMMARY_PATH = PILOT_DIR / "summary.json"

CONTINUOUS_PARAMETERS = ["joint", "prior", "likelihood", "rootHeight", "treeLength",
                         "constant.popSize", "clock.rate"]

CONDITION_FIELDS = ["run_id", "phase", "config", "chains", "tempered", "threads", "seed",
                    "minutes", "states", "log_every", "ladder_ratio", "swap_every",
                    "processes", "gc", "device", "gpus", "gpu_chains", "numa"]

CPU_DEVICE_ARGS = ["-beagle_CPU", "-beagle_SSE"]

# chain lengths for runs that are stopped at their time limit rather than finishing
UNBOUNDED_STATES = 2_000_000_000


def load_config() -> dict:
    with open(CONFIG_PATH) as f:
        config = json.load(f)
    if "work_dir" in config:
        global RUNS_DIR, PILOT_DIR, ANALYSIS_DIR, CONDITIONS_PATH, PILOT_SUMMARY_PATH
        base = (STUDY_DIR / config["work_dir"]).resolve()
        RUNS_DIR = base / "runs"
        PILOT_DIR = base / "pilot"
        ANALYSIS_DIR = base / "analysis"
        CONDITIONS_PATH = base / "conditions.csv"
        PILOT_SUMMARY_PATH = PILOT_DIR / "summary.json"
    return config


def paths():
    return dict(runs=RUNS_DIR, pilot=PILOT_DIR, analysis=ANALYSIS_DIR,
                conditions=CONDITIONS_PATH, pilot_summary=PILOT_SUMMARY_PATH)


def read_conditions():
    with open(CONDITIONS_PATH) as f:
        return list(csv.DictReader(f))


def write_conditions(rows, append=False):
    CONDITIONS_PATH.parent.mkdir(parents=True, exist_ok=True)
    existing = read_conditions() if append and CONDITIONS_PATH.exists() else []
    seen = {r["run_id"] for r in existing}
    rows = existing + [r for r in rows if r["run_id"] not in seen]
    with open(CONDITIONS_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CONDITION_FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in CONDITION_FIELDS})
    return rows


def load_pilot_summary():
    if not PILOT_SUMMARY_PATH.exists():
        return None
    with open(PILOT_SUMMARY_PATH) as f:
        return json.load(f)


def load_reference(config) -> dict:
    with gzip.open(REPO_DIR / config["reference"], "rt") as f:
        ref = json.load(f)
    freqs = clades_from_json(ref["clade_frequencies"])
    focal = {int(c, 16) for c in ref["focal_clades"]}
    uncertain = sorted((c for c, f in freqs.items() if 0.05 <= f <= 0.95),
                       key=lambda c: abs(freqs[c] - 0.5))[:config["focal_clades"]]
    return dict(frequencies=freqs, focal=focal, uncertain=uncertain,
                uncertain_freq=[freqs[c] for c in uncertain])


def ladder(ratio: float, chains: int):
    """Inverse temperatures of the heated chains of a set: ratio, ratio^2, ..."""
    return [round(ratio ** i, 6) for i in range(1, chains)]


def make_xml(config, cond, run_dir: Path) -> Path:
    """Writes the run's XML: per-sample wall-clock time in the file log, the run's log
    interval, its chain length and output file names."""
    text = (REPO_DIR / config["xml"]).read_text()
    run_id = cond["run_id"]
    log_every = int(cond["log_every"])
    states = int(cond["states"]) if cond.get("states") else UNBOUNDED_STATES

    text = re.sub(r'(<mcmc id="mcmc" chainLength=")\d+"', rf'\g<1>{states}"', text)
    text = re.sub(r'\s+operatorAnalysis="[^"]*"', "", text)
    screen_every = min(log_every * 20, 2**31 - 1)  # BEAST reads logEvery as an int
    text = re.sub(r'(<log id="screenLog" logEvery=")\d+"', rf'\g<1>{screen_every}"', text)
    text = re.sub(r'(<log id="fileLog" logEvery=")\d+"', rf'\g<1>{log_every}"', text)
    text = re.sub(r'(<logTree id="treeFileLog" logEvery=")\d+"', rf'\g<1>{log_every}"', text)
    text = re.sub(r'fileName="[^"]*\.log"', f'fileName="{run_id}.log"', text)
    text = re.sub(r'fileName="[^"]*\.trees"', f'fileName="{run_id}.trees"', text)
    # seconds since the chain's first log line, for ESS against wall-clock time
    text, n = re.subn(r'(<log id="fileLog"[^>]*>\n)', r'\1\t\t\t<timeLogger/>\n', text)
    if n != 1:
        raise ValueError("could not add timeLogger to the file log")

    path = run_dir / f"{run_id}.xml"
    path.write_text(text)
    return path


def beagle_order(config, cond):
    """BEAGLE resource for each instance, in creation order (chain 1's cold chain, its heated
    chains, then chain 2, ...): GPU chains alternate over the GPUs used, CPU chains (in a
    mixed configuration, after the GPU chains) get resource 0."""
    device = cond.get("device") or "cpu"
    if device == "cpu":
        return None
    gpus = config["machine"]["gpu_resources"][:int(cond.get("gpus") or 1)]
    chains = int(cond["chains"])
    tempered = int(cond.get("tempered") or 1)
    gpu_chains = chains if device == "gpu" else int(cond["gpu_chains"])
    order = []
    for k in range(chains):
        for _ in range(tempered):
            if k < gpu_chains:
                order.append(gpus[len(order) % len(gpus)])
            else:
                order.append(0)
    return order


def beast_command(config, jar: Path, cond, xml_name: str, seed=None, process_index=0) -> list:
    java = list(config["java_args"])
    if int(cond.get("processes") or 1) > 1:
        # one small JVM per single-threaded chain
        java = [a for a in java if not a.startswith(("-Xms", "-Xmx"))] + ["-Xms1g", "-Xmx1g"]
    if cond.get("gc"):
        java = [a for a in java if not a.startswith("-XX:+Use") or "GC" not in a] + [cond["gc"]]
    beast = list(config["beast_args"])
    order = beagle_order(config, cond)
    if order is None:
        beast += [a for a in CPU_DEVICE_ARGS if a not in beast]
    else:
        beast = [a for a in beast if a not in CPU_DEVICE_ARGS]
        beast += ["-beagle_GPU", "-beagle_order", ",".join(str(r) for r in order)]
    cmd = ["java"] + java + ["-cp", str(jar), "dr.app.beast.BeastMain"] + beast
    cmd += ["-seed", str(seed if seed is not None else cond["seed"]),
            "-beagle_threads", str(cond["threads"])]

    chains = int(cond["chains"])
    tempered = int(cond.get("tempered") or 1)
    if chains > 1 or tempered > 1:
        cmd += ["-chains", str(chains), "-chain_threads", str(chains * tempered),
                "-convergence_every", "0"]
    if tempered > 1:
        temps = ladder(float(cond["ladder_ratio"]), tempered)
        cmd += ["-mc3_chains", str(tempered),
                "-mc3_temperatures", ",".join(str(t) for t in temps),
                "-mc3_swap", str(cond["swap_every"]),
                "-mc3_scheme", config["tempering"]["scheme"]]
    cmd += [xml_name]
    if cond.get("numa") and int(cond["numa"]):
        # one JVM per NUMA node, bound to its cores and memory
        node = process_index % config["machine"]["numa_nodes"]
        cmd = ["numactl", f"--cpunodebind={node}", f"--membind={node}"] + cmd
    return cmd


def chain_files(run_dir: Path, run_id: str, cond) -> list:
    """The (log, trees) files of each chain whose samples count: one per independent chain,
    or one per tempering set (its cold chain)."""
    chains = int(cond["chains"])
    processes = int(cond.get("processes") or 1)
    dirs = [run_dir / f"p{i + 1:02d}" for i in range(processes)] if processes > 1 else [run_dir]
    files = []
    for d in dirs:
        if chains == 1:
            files.append((d / f"{run_id}.log", d / f"{run_id}.trees"))
        else:
            files += [(d / f"{run_id}.chain{k}.log", d / f"{run_id}.chain{k}.trees")
                      for k in range(1, chains + 1)]
    return files
