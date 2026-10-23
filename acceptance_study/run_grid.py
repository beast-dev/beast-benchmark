#!/usr/bin/env python3
"""
Runs the distinct runs in a conditions table through a pool of BEAST workers.

Each run gets runs/RUN_ID/ holding the variant XML, BEAST's stdout, the .log,
.trees and .ops files, and result.json. Runs with a result.json whose exit code
is 0 are skipped, so the grid can be stopped and restarted.

Runs start longest first (by study.json minutes_per_mstate x states), so the
long reference and tall-1 runs do not finish last; within a dataset the order
is shuffled so machine load and thermal drift are spread across conditions.

Ctrl-C, or `kill PID` of this script, stops every BEAST process it started and
launches no more; the interrupted runs are redone on restart. Each run holds a lock on its directory
while BEAST runs, so a second run_grid.py on the same runs skips those in progress.

The BEAST jar is copied to beast-frozen.jar on first use so a rebuild in the
BEAST tree cannot change the study part way through; delete it to pick up a
new build deliberately.

Usage:
  run_grid.py --check [--arm ...]     pre-flight: short BEAST runs of each kind of condition
  run_grid.py [--arm A2 --dataset small-1 --workers 12 --dry-run] [conditions.csv]
"""

import argparse
import csv
import fcntl
import hashlib
import json
import random
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path

from common import REPO_DIR, load_config, read_log, read_ops, read_trees, work_dir
from make_variants import write_variant

sys.path.insert(0, str(REPO_DIR))
from benchmark import get_system_info, parse_stdout  # noqa: E402

# set from the config in main()
RUNS_DIR = FROZEN_JAR = None

print_lock = threading.Lock()
stopping = threading.Event()
running = set()  # BEAST processes in flight, so Ctrl-C can stop them
running_lock = threading.Lock()


def log(message):
    with print_lock:
        print(f"[{datetime.now():%H:%M:%S}] {message}", flush=True)


def frozen_jar(config) -> tuple:
    if not FROZEN_JAR.exists():
        source = (REPO_DIR / config["beast_jar"]).resolve()
        shutil.copy2(source, FROZEN_JAR)
        log(f"Froze {source} -> {FROZEN_JAR.name}")
    digest = hashlib.sha256(FROZEN_JAR.read_bytes()).hexdigest()[:16]
    return FROZEN_JAR, digest


def is_complete(run_dir) -> bool:
    # a missing or unreadable result (e.g. cut off by Ctrl-C) means the run is redone
    try:
        return json.loads((run_dir / "result.json").read_text()).get("exit_code") == 0
    except (OSError, ValueError):
        return False


def expected_minutes(config, cond) -> float:
    rate = config["datasets"][cond["dataset"]].get("minutes_per_mstate", {}).get(cond["likelihood"], 1.0)
    return rate * int(cond["states"]) / 1e6


def beast_command(config, jar, cond, xml_name=None) -> list:
    cmd = ["java", *config["java_args"], "-cp", str(jar), "dr.app.beast.BeastMain",
           *config["beast_args"], "-seed", str(cond["seed"])]
    if cond["adaptation_target"]:
        cmd += ["-adaptation_target", str(cond["adaptation_target"])]
    if cond["operator_weights"]:
        cmd += ["-operator_weights", cond["operator_weights"]]
    cmd.append(xml_name or f"{cond['run_id']}.xml")
    return cmd


def run_one(config, jar, jar_hash, cond, system):
    """Returns (run_id, exit code or None if skipped or interrupted, wall seconds)."""
    if stopping.is_set():
        return cond["run_id"], None, 0.0
    run_dir = RUNS_DIR / cond["run_id"]
    run_dir.mkdir(parents=True, exist_ok=True)

    # the lock is released by the OS if this process dies, so it can never go stale
    lock = open(run_dir / "run.lock", "w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        log(f"{cond['run_id']}: already running in another run_grid.py, skipped")
        return cond["run_id"], None, 0.0
    try:
        if is_complete(run_dir):
            return cond["run_id"], None, 0.0
        return run_locked(config, jar, jar_hash, cond, system, run_dir)
    finally:
        lock.close()


def run_beast(cmd, cwd, stdout_path) -> int:
    with open(stdout_path, "w") as out:
        proc = subprocess.Popen(cmd, cwd=cwd, stdout=out, stderr=subprocess.STDOUT)
        with running_lock:
            running.add(proc)
        try:
            return proc.wait()
        finally:
            with running_lock:
                running.discard(proc)


def run_locked(config, jar, jar_hash, cond, system, run_dir):
    write_variant(config, cond, run_dir / f"{cond['run_id']}.xml")
    cmd = beast_command(config, jar, cond)

    start = time.time()
    exit_code = run_beast(cmd, run_dir, run_dir / "beast.stdout")
    wall = time.time() - start
    if stopping.is_set():
        return cond["run_id"], None, wall  # interrupted: no result, so it is redone on restart

    parsed = parse_stdout((run_dir / "beast.stdout").read_text())
    settings = {k: v for k, v in cond.items() if k not in ("arm", "level")}
    result = {
        **settings,
        "exit_code": exit_code,
        "wall_seconds": wall,
        "beast_minutes": parsed["beast_minutes"],
        "final_rate_min_per_mstate": parsed["final_rate"],
        "command": cmd,
        "jar_sha256": jar_hash,
        "finished": datetime.now().isoformat(timespec="seconds"),
        "system": system,
    }
    # write then rename, so an interrupted write never leaves a partial result.json
    tmp = run_dir / "result.json.tmp"
    tmp.write_text(json.dumps(result, indent=2))
    tmp.replace(run_dir / "result.json")
    return cond["run_id"], exit_code, wall


def schedule(config, conds, seed):
    """Longest expected runs first; shuffled within each (dataset, likelihood, length)."""
    rng = random.Random(seed)
    conds = list(conds)
    rng.shuffle(conds)
    return sorted(conds, key=lambda c: -expected_minutes(config, c))


# ---------------------------------------------------------------------------
# pre-flight check
# ---------------------------------------------------------------------------

def performance_cores():
    try:
        return int(subprocess.check_output(["sysctl", "-n", "hw.perflevel0.physicalcpu"], text=True))
    except (OSError, subprocess.CalledProcessError, ValueError):
        return None


def check(config, conds, workers) -> bool:
    ok = True

    def report(passed, message):
        nonlocal ok
        ok &= passed
        log(("ok    " if passed else "FAIL  ") + message)

    others = subprocess.run(["pgrep", "-f", "BeastMain"], capture_output=True, text=True).stdout.split()
    report(not others, f"no other BEAST processes running ({len(others)} found)")
    cores = performance_cores()
    if cores:
        report(workers <= cores - 2, f"{workers} workers for {cores} performance cores (keep at least 2 free)")

    source = (REPO_DIR / config["beast_jar"]).resolve()
    jar = FROZEN_JAR if FROZEN_JAR.exists() else source
    report(jar.exists(), f"BEAST jar {jar}" + (" (frozen)" if jar == FROZEN_JAR else " (not yet frozen)"))
    if not jar.exists():
        return False

    # one short run per kind of condition: each dataset/likelihood and each kind of setting
    kinds = {}
    for c in conds:
        for key in [("data", c["dataset"], c["likelihood"]), ("weights", bool(c["operator_weights"])),
                    ("target", bool(c["adaptation_target"])), ("sl_target", bool(c["sl_target"])),
                    ("fixed_size", c["sl_adapt"] == "false"), ("kernel", c["sl_kernel"])]:
            if key not in kinds and (key[0] == "data" or key[1]):
                kinds[key] = c
    samples_config = dict(config, samples=20)
    for cond in {id(c): c for c in kinds.values()}.values():
        short = dict(cond, states=str(20000))
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            name = f"{cond['run_id']}.xml"
            write_variant(samples_config, short, tmp / name)
            start = time.time()
            exit_code = run_beast(beast_command(config, jar, short, name), tmp, tmp / "beast.stdout")
            label = f"{cond['run_id']} (20k states, {time.time() - start:.0f} s)"
            if exit_code != 0:
                tail = (tmp / "beast.stdout").read_text().strip().splitlines()[-5:]
                report(False, f"{label}: exit {exit_code}\n        " + "\n        ".join(tail))
                continue
            try:
                columns, data = read_log(tmp / f"{cond['run_id']}.log")
                _, trees = read_trees(tmp / f"{cond['run_id']}.trees")
                ops = read_ops(tmp / f"{cond['run_id']}.ops")
                sl = next((o for n, o in ops.items() if n.startswith("subtreeLeap")), None)
                problems = []
                if "topologyChanged" not in columns:
                    problems.append("no SubtreeLeap columns in the log")
                if not trees:
                    problems.append("no trees")
                if not sl or sl["time_per_op_us"] <= 0:
                    problems.append("no SubtreeLeap timing in the .ops (old BEAST build?)")
                if cond["operator_weights"] and "Operator weight override" not in (tmp / "beast.stdout").read_text():
                    problems.append("-operator_weights not applied")
                report(not problems, label + (": " + "; ".join(problems) if problems else
                                              f": {len(data)} log rows, {len(trees)} trees, "
                                              f"SubtreeLeap accept {sl['accept']:.2f}"))
            except Exception as e:  # report any parse failure rather than stopping the check
                report(False, f"{label}: could not read the output ({e})")
    return ok


def interrupt(signum, frame):
    raise KeyboardInterrupt


def main():
    # handle these explicitly: a shell starts background jobs with SIGINT ignored
    signal.signal(signal.SIGINT, interrupt)
    signal.signal(signal.SIGTERM, interrupt)

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("conditions", nargs="?", help="default: WORK_DIR/conditions.csv")
    parser.add_argument("--arm", action="append", help="only runs belonging to this arm (repeatable)")
    parser.add_argument("--dataset", action="append", help="only runs on this dataset (repeatable)")
    parser.add_argument("--workers", type=int)
    parser.add_argument("--shuffle-seed", type=int, default=1)
    parser.add_argument("--check", action="store_true", help="pre-flight check, then exit")
    parser.add_argument("--dry-run", action="store_true", help="list what would run, and write the XMLs")
    args = parser.parse_args()

    global RUNS_DIR, FROZEN_JAR
    config = load_config()
    RUNS_DIR = work_dir(config) / "runs"
    FROZEN_JAR = work_dir(config) / "beast-frozen.jar"
    with open(args.conditions or work_dir(config) / "conditions.csv") as f:
        rows = list(csv.DictReader(f))
    if args.arm:
        rows = [r for r in rows if r["arm"] in args.arm]
    if args.dataset:
        rows = [r for r in rows if r["dataset"] in args.dataset]

    runs = {}
    for r in rows:
        runs.setdefault(r["run_id"], r)
    pending = schedule(config, [c for rid, c in runs.items() if not is_complete(RUNS_DIR / rid)],
                       args.shuffle_seed)
    workers = args.workers or config["workers"]

    if args.check:
        sys.exit(0 if check(config, list(runs.values()), workers) else 1)

    minutes = sum(expected_minutes(config, c) for c in pending)
    longest = max((expected_minutes(config, c) for c in pending), default=0)
    estimate = max(minutes / workers, longest)
    log(f"{len(runs)} distinct runs selected, {len(runs) - len(pending)} complete, {len(pending)} to run "
        f"on {workers} workers: about {minutes / 60:.0f} core-hours, finishing around "
        f"{datetime.now() + timedelta(minutes=estimate):%a %H:%M}")

    if args.dry_run:
        for cond in pending:
            run_dir = RUNS_DIR / cond["run_id"]
            run_dir.mkdir(parents=True, exist_ok=True)
            write_variant(config, cond, run_dir / f"{cond['run_id']}.xml")
            print(" ".join(beast_command(config, "beast.jar", cond)))
        return

    jar, jar_hash = frozen_jar(config)
    system = get_system_info()
    failures = 0
    done = 0
    pool = ThreadPoolExecutor(max_workers=workers)
    try:
        futures = [pool.submit(run_one, config, jar, jar_hash, cond, system) for cond in pending]
        for future in as_completed(futures):
            rid, exit_code, wall = future.result()
            if exit_code is None:
                continue
            done += 1
            failures += exit_code != 0
            log(f"{done}/{len(pending)} {rid}: exit {exit_code}, {wall / 60:.1f} min")
    except KeyboardInterrupt:
        stopping.set()
        with running_lock:
            procs = list(running)
        log(f"Interrupted: stopping {len(procs)} BEAST runs; they will be redone on restart")
        for proc in procs:
            proc.terminate()
        pool.shutdown(wait=True, cancel_futures=True)
        sys.exit(130)
    pool.shutdown()
    log(f"Done, {failures} failed")


if __name__ == "__main__":
    main()
