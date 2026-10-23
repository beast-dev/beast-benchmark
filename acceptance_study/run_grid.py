#!/usr/bin/env python3
"""
Runs the distinct runs in a conditions table through a pool of BEAST workers.

Each run gets runs/RUN_ID/ holding the variant XML, BEAST's stdout, the .log,
.trees and .ops files, and result.json. Runs with a result.json whose exit code
is 0 are skipped, so the grid can be stopped and restarted. Order is shuffled so
that machine load and thermal drift are spread across conditions.

The BEAST jar is copied to beast-frozen.jar on first use so a rebuild in the
BEAST tree cannot change the study part way through; delete it to pick up a
new build deliberately.

Usage:
  run_grid.py [conditions.csv] [--arm A2 --dataset small-1 --workers 12 --dry-run]
"""

import argparse
import csv
import hashlib
import json
import random
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime

from common import REPO_DIR, load_config, work_dir
from make_variants import write_variant

sys.path.insert(0, str(REPO_DIR))
from benchmark import get_system_info, parse_stdout  # noqa: E402

# set from the config in main()
RUNS_DIR = FROZEN_JAR = None

print_lock = threading.Lock()


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


def beast_command(config, jar, cond) -> list:
    cmd = ["java", *config["java_args"], "-cp", str(jar), "dr.app.beast.BeastMain",
           *config["beast_args"], "-seed", str(cond["seed"])]
    if cond["adaptation_target"]:
        cmd += ["-adaptation_target", str(cond["adaptation_target"])]
    if cond["operator_weights"]:
        cmd += ["-operator_weights", cond["operator_weights"]]
    cmd.append(f"{cond['run_id']}.xml")
    return cmd


def run_one(config, jar, jar_hash, cond, system):
    run_dir = RUNS_DIR / cond["run_id"]
    run_dir.mkdir(parents=True, exist_ok=True)
    write_variant(config, cond, run_dir / f"{cond['run_id']}.xml")
    cmd = beast_command(config, jar, cond)

    start = time.time()
    with open(run_dir / "beast.stdout", "w") as out:
        exit_code = subprocess.call(cmd, cwd=run_dir, stdout=out, stderr=subprocess.STDOUT)
    wall = time.time() - start

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


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("conditions", nargs="?", help="default: WORK_DIR/conditions.csv")
    parser.add_argument("--arm", action="append", help="only runs belonging to this arm (repeatable)")
    parser.add_argument("--dataset", action="append", help="only runs on this dataset (repeatable)")
    parser.add_argument("--workers", type=int)
    parser.add_argument("--shuffle-seed", type=int, default=1)
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
    pending = [c for rid, c in runs.items() if not is_complete(RUNS_DIR / rid)]
    random.Random(args.shuffle_seed).shuffle(pending)

    workers = args.workers or config["workers"]
    log(f"{len(runs)} distinct runs selected, {len(runs) - len(pending)} complete, {len(pending)} to run "
        f"on {workers} workers")

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
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(run_one, config, jar, jar_hash, cond, system) for cond in pending]
        for i, future in enumerate(as_completed(futures), 1):
            rid, exit_code, wall = future.result()
            failures += exit_code != 0
            log(f"{i}/{len(pending)} {rid}: exit {exit_code}, {wall / 60:.1f} min")
    log(f"Done, {failures} failed")


if __name__ == "__main__":
    main()
