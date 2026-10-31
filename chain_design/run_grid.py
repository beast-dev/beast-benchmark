#!/usr/bin/env python3
"""
Runs the conditions in conditions.csv one at a time, each with the whole machine.

  python3 run_grid.py [--phase PHASE] [--check] [--dry-run]

Each run is stopped at its time limit (minutes in conditions.csv) if it has not finished.
BEAST's file logs are flushed line by line, so a stopped run's logs are complete up to the
last sample. Finished runs are skipped on restart; an interrupted run is redone. Ctrl-C
stops the current run and exits.
"""

import argparse
import fcntl
import hashlib
import json
import os
import shutil
import signal
import socket
import subprocess
import sys
import threading
import time
from pathlib import Path

import cd_common as cc

stopping = threading.Event()


def log(message):
    print(time.strftime("%Y-%m-%d %H:%M:%S"), message, flush=True)


def frozen_jar(config) -> tuple:
    """The first run copies the BEAST jar into the study directory; every run uses that copy."""
    frozen = cc.STUDY_DIR / "beast-frozen.jar"
    if not frozen.exists():
        source = (cc.REPO_DIR / config["beast_jar"]).resolve()
        if not source.exists():
            sys.exit(f"BEAST jar not found: {source} (build the multichain_runner branch with `ant dist`)")
        shutil.copy2(source, frozen)
        log(f"froze {source} -> {frozen}")
    digest = hashlib.sha256(frozen.read_bytes()).hexdigest()[:12]
    return frozen, digest


def load_average():
    return list(os.getloadavg())


def first_log_times(files, start, deadline, times):
    """Records, for each chain's log, the seconds from process start to its first sample."""
    pending = {str(f) for f in files}
    while pending and time.time() < deadline and not stopping.is_set():
        for name in list(pending):
            p = Path(name)
            if p.exists():
                with open(p) as f:
                    data = [line for line in f if line[:1].isdigit()]
                if data:
                    times[name] = time.time() - start
                    pending.discard(name)
        time.sleep(2)


def log_size(paths):
    return sum(p.stat().st_size for p in map(Path, paths) if p.exists())


def run_processes(cmds, cwds, stdout_paths, limit_seconds, chain_logs, stall_seconds, dump_path):
    """Starts the processes and stops them at the time limit, or if the chain logs stop growing
    for stall_seconds (a thread dump is saved first). Returns (exit codes, stopped, stalled,
    start, end, first log times)."""
    start = time.time()
    procs = []
    for cmd, cwd, out in zip(cmds, cwds, stdout_paths):
        f = open(out, "w")
        procs.append((subprocess.Popen(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT,
                                       start_new_session=True), f))
    times = {}
    watcher = threading.Thread(target=first_log_times,
                               args=(chain_logs, start, start + limit_seconds, times), daemon=True)
    watcher.start()

    stopped = False
    stalled = False
    deadline = start + limit_seconds
    last_size, last_growth = -1, time.time()
    while any(p.poll() is None for p, _ in procs):
        size = log_size(chain_logs)
        if size != last_size:
            last_size, last_growth = size, time.time()
        elif size > 0 and time.time() - last_growth > stall_seconds:
            stalled = True
            with open(dump_path, "w") as dump:
                for p, _ in procs:
                    if p.poll() is None:
                        subprocess.run(["jcmd", str(p.pid), "Thread.print"], stdout=dump,
                                       stderr=subprocess.STDOUT, timeout=60)
            log(f"STALLED: no log output for {stall_seconds / 60:.0f} min; thread dump in {dump_path}")
        if stopping.is_set() or stalled or time.time() >= deadline:
            stopped = not stopping.is_set() and not stalled
            for p, _ in procs:
                if p.poll() is None:
                    os.killpg(p.pid, signal.SIGTERM)
            for p, _ in procs:
                try:
                    p.wait(timeout=60)
                except subprocess.TimeoutExpired:
                    os.killpg(p.pid, signal.SIGKILL)
            break
        time.sleep(2)
    end = time.time()
    for p, f in procs:
        p.wait()
        f.close()
    watcher.join(timeout=5)
    return [p.returncode for p, _ in procs], stopped, stalled, start, end, times


def run_one(config, jar, digest, cond):
    run_id = cond["run_id"]
    run_dir = cc.paths()["runs"] / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    lock = open(run_dir / "run.lock", "w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        log(f"{run_id}: running in another run_grid.py, skipped")
        return
    try:
        if (run_dir / "result.json").exists():
            return
        processes = int(cond.get("processes") or 1)
        cmds, cwds, outs = [], [], []
        if processes > 1:
            # separate JVMs, one single-threaded chain each (pilot P3)
            for i in range(processes):
                d = run_dir / f"p{i + 1:02d}"
                d.mkdir(exist_ok=True)
                xml = cc.make_xml(config, cond, d)
                cmds.append(cc.beast_command(config, jar, cond, xml.name, seed=int(cond["seed"]) + i))
                cwds.append(d)
                outs.append(d / "beast.out")
        else:
            xml = cc.make_xml(config, cond, run_dir)
            cmd = cc.beast_command(config, jar, cond, xml.name)
            cmd.insert(1, f"-Xlog:gc:file={run_dir / 'gc.log'}")  # GC pause fraction
            cmds, cwds, outs = [cmd], [run_dir], [run_dir / "beast.out"]

        chain_logs = [lg for lg, _ in cc.chain_files(run_dir, run_id, cond)]
        limit = float(cond["minutes"]) * 60
        # runs with a chain length finish by themselves; the limit is then only a backstop
        if cond.get("states"):
            limit *= 1.25
        load_before = load_average()
        log(f"{run_id}: start ({cond['config']}, {processes} process(es), limit {limit / 60:.0f} min)")
        codes, stopped, stalled, start, end, first_logs = run_processes(
            cmds, cwds, outs, limit, chain_logs, config.get("stall_minutes", 10) * 60,
            run_dir / f"stall_{time.strftime('%Y%m%d_%H%M%S')}.txt")
        if stopping.is_set():
            log(f"{run_id}: interrupted, will be redone")
            return
        if stalled:
            # keep the evidence, then redo the run from scratch next time
            keep = run_dir.parent / f"{run_id}.stalled_{time.strftime('%Y%m%d_%H%M%S')}"
            run_dir.rename(keep)
            log(f"{run_id}: stalled; moved to {keep.name}, will be redone")
            return
        result = dict(cond)
        result.update(dict(
            exit_codes=codes, stopped_at_limit=stopped, start=start, end=end,
            wall_seconds=end - start, first_log_seconds=first_logs,
            load_before=load_before, load_after=load_average(),
            host=socket.gethostname(), jar=str(jar), jar_sha256=digest, command=cmds[0]))
        ok = all(c == 0 for c in codes) or stopped
        if not ok:
            log(f"{run_id}: FAILED (exit codes {codes}); see {outs[0]}")
            result["failed"] = True
        with open(run_dir / "result.json", "w") as f:
            json.dump(result, f, indent=1)
        log(f"{run_id}: done in {(end - start) / 60:.1f} min{' (stopped at limit)' if stopped else ''}")
    finally:
        lock.close()


def check(config, conds):
    jar = (cc.REPO_DIR / config["beast_jar"]).resolve()
    problems = []
    if not jar.exists() and not (cc.STUDY_DIR / "beast-frozen.jar").exists():
        problems.append(f"no BEAST jar at {jar}")
    if not (cc.REPO_DIR / config["xml"]).exists():
        problems.append(f"no XML {config['xml']}")
    try:
        cc.load_reference(config)
    except Exception as e:  # noqa: BLE001
        problems.append(f"reference: {e}")
    tmp = cc.STUDY_DIR / "check_tmp"
    tmp.mkdir(exist_ok=True)
    try:
        cc.make_xml(config, conds[0], tmp)
    except Exception as e:  # noqa: BLE001
        problems.append(f"XML variant: {e}")
    shutil.rmtree(tmp, ignore_errors=True)
    load = load_average()[0]
    if load > 2:
        problems.append(f"load average is {load:.1f}: something else is running")
    for p in problems:
        log("PROBLEM: " + p)
    if not problems:
        log(f"check OK: {len(conds)} runs, load {load:.1f}")
    return not problems


def interrupt(signum, frame):
    log("stopping (the current run will be redone next time)")
    stopping.set()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", action="append", help="only these phases (default: all)")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="print the commands")
    args = ap.parse_args()

    config = cc.load_config()
    conds = cc.read_conditions()
    if args.phase:
        conds = [c for c in conds if c["phase"] in args.phase]
    if not conds:
        sys.exit("no conditions (run design.py)")
    if args.check:
        sys.exit(0 if check(config, conds) else 1)

    if args.dry_run:
        for c in conds:
            print(c["run_id"], " ".join(cc.beast_command(config, Path("beast.jar"), c, c["run_id"] + ".xml")))
        return

    signal.signal(signal.SIGINT, interrupt)
    signal.signal(signal.SIGTERM, interrupt)
    jar, digest = frozen_jar(config)
    todo = [c for c in conds if not (cc.paths()["runs"] / c["run_id"] / "result.json").exists()]
    hours = sum(float(c["minutes"]) for c in todo) / 60
    log(f"{len(todo)} of {len(conds)} runs to do, about {hours:.1f} hours")
    for cond in todo:
        if stopping.is_set():
            break
        run_one(config, jar, digest, cond)


if __name__ == "__main__":
    main()
