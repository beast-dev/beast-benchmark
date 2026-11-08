#!/usr/bin/env python3
"""
Measure how fast each benchmark XML runs on this machine, and propose the chain
length that would make it take a target wall-clock time.

The benchmarks exist to compare wall-clock time and states per second, so they are
most useful when they all take about the same time. This script does not run a
chain to the end: it starts BEAST, waits for the screen log, times a stretch of
states, and kills it.

    python3 calibrate.py                       # all 42, 10 minute target
    python3 calibrate.py --target 300          # 5 minute target
    python3 calibrate.py small 1 min           # just one
    python3 calibrate.py --beast /path/to/beast --machine machine1

Unrecognised options are passed through to BEAST, so --beagle_gpu works here as it
does for benchmark.py.

Writes calibration.csv and prints a table. It changes no XML: feed the proposed
lengths to apply_lengths.py.
"""

import argparse
import csv
import re
import subprocess
import sys
import tempfile
import threading
import time
from datetime import datetime
from pathlib import Path

from benchmarks import discover, xml_path

STATE_HEADER_RE = re.compile(r'^state\t')
STATE_ROW_RE = re.compile(r'^(\d+)\t')


def round_states(states: float) -> int:
    """Two significant figures, and a whole multiple of 1000 so logEvery stays integral."""
    if states < 1000:
        return 1000
    magnitude = 10 ** (len(str(int(states))) - 2)
    rounded = round(states / magnitude) * magnitude
    return max(1000, int(round(rounded / 1000)) * 1000)


def time_one(xml: Path, beast: str, seed: int, measure: float, max_wait: float,
             extra_args: list) -> dict:
    """Run BEAST briefly and return startup time and seconds per state."""
    cmd = [beast, "-seed", str(seed), "-citations_off", "-tests", "0", "-overwrite"]
    cmd += extra_args
    cmd.append(str(xml))

    samples = []          # (state, elapsed seconds since launch)
    header_found = False
    started = time.monotonic()
    startup = None
    error = None

    with tempfile.TemporaryDirectory(prefix="calibrate_") as workdir:
        try:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    text=True, bufsize=1, cwd=workdir)
        except FileNotFoundError:
            return {"error": f"beast not found: {beast}"}

        # If BEAST ever stops producing output the read below would block for good,
        # so give it a hard stop as well as the per-line check.
        watchdog = threading.Timer(max_wait, proc.kill)
        watchdog.start()
        try:
            for line in proc.stdout:
                elapsed = time.monotonic() - started
                stripped = line.rstrip("\n")

                if elapsed > max_wait:
                    break

                if STATE_HEADER_RE.match(stripped):
                    header_found = True
                    continue

                if header_found and STATE_ROW_RE.match(stripped):
                    state = int(STATE_ROW_RE.match(stripped).group(1))
                    samples.append((state, elapsed))
                    if startup is None:
                        startup = elapsed
                    # enough rows, and enough chain timed since the first one
                    if len(samples) >= 4 and elapsed - samples[0][1] >= measure:
                        break
        finally:
            watchdog.cancel()
            proc.kill()
            proc.wait()

    if not samples:
        error = f"no state rows within {max_wait:.0f}s"

    if error:
        return {"error": error}
    if len(samples) < 2:
        return {"error": f"only {len(samples)} state row(s) seen"}

    # Drop the state 0 row where we can: the first MCMC interval carries warm-up.
    first = samples[1] if len(samples) >= 3 else samples[0]
    last = samples[-1]
    d_states = last[0] - first[0]
    d_time = last[1] - first[1]
    if d_states <= 0 or d_time <= 0:
        return {"error": "no measurable progress"}

    return {
        "startup_s": round(startup, 2),
        "sec_per_state": d_time / d_states,
        "rows_used": len(samples),
        "states_timed": d_states,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Propose chain lengths that hit a target wall-clock time.")
    parser.add_argument("type", nargs="?", help="Benchmark type (default: all)")
    parser.add_argument("number", nargs="?", type=int, help="Benchmark number")
    parser.add_argument("model", nargs="?", help="Model label")
    parser.add_argument("--target", type=float, default=600.0,
                        help="Target wall-clock seconds per run (default: 600)")
    parser.add_argument("--measure", type=float, default=60.0,
                        help="Seconds of chain to time, after the first row (default: 60)")
    parser.add_argument("--max-wait", type=float, default=300.0,
                        help="Give up on a run after this many seconds (default: 300)")
    parser.add_argument("--seed", type=int, default=1234, help="Random seed (default: 1234)")
    parser.add_argument("--beast", default="beast", help="Path to beast (default: on PATH)")
    parser.add_argument("--benchmarks", default=None, metavar="DIR")
    parser.add_argument("--out", default="calibration.csv", metavar="FILE")
    parser.add_argument("--machine", default=None, metavar="NAME",
                        help="Machine identifier, recorded in the csv")
    args, extra_beast_args = parser.parse_known_args()

    here = Path(__file__).resolve().parent
    benchmarks = Path(args.benchmarks) if args.benchmarks else here / "Benchmarks"
    if not benchmarks.is_dir():
        sys.exit(f"No Benchmarks folder at {benchmarks}")

    if args.type:
        if args.number is None or args.model is None:
            sys.exit("Give all three of type, number and model, or none of them.")
        targets = [(args.type.lower(), args.number, args.model)]
    else:
        targets = discover(benchmarks)

    print(f"Calibrating {len(targets)} benchmark(s) for a {args.target:.0f}s target")
    print(f"  beast: {args.beast}" + (f"  extra: {' '.join(extra_beast_args)}" if extra_beast_args else ""))
    print()

    rows = []
    for btype, number, model in targets:
        xml = xml_path(benchmarks, btype, number, model)
        name = xml.stem
        if not xml.exists():
            print(f"  {name:36s} MISSING")
            continue

        print(f"  {name:36s} ", end="", flush=True)
        result = time_one(xml, args.beast, args.seed, args.measure, args.max_wait,
                          extra_beast_args)

        if "error" in result:
            print(f"FAILED: {result['error']}")
            rows.append({"benchmark": name, "type": btype, "number": number, "model": model,
                         "error": result["error"]})
            continue

        usable = max(args.target - result["startup_s"], 1.0)
        proposed = round_states(usable / result["sec_per_state"])
        min_per_mstate = result["sec_per_state"] * 1e6 / 60.0

        print(f"{min_per_mstate:8.3f} min/Mstate  startup {result['startup_s']:6.1f}s"
              f"  ->{proposed:>12,} states")

        rows.append({
            "benchmark": name, "type": btype, "number": number, "model": model,
            "machine": args.machine or "",
            "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "startup_s": result["startup_s"],
            "min_per_mstate": round(min_per_mstate, 4),
            "states_timed": result["states_timed"],
            "rows_used": result["rows_used"],
            "target_s": args.target,
            "proposed_states": proposed,
            "proposed_log_every": proposed // 1000,
            "error": "",
        })

    out = Path(args.out)
    fields = ["benchmark", "type", "number", "model", "machine", "timestamp", "startup_s",
              "min_per_mstate", "states_timed", "rows_used", "target_s", "proposed_states",
              "proposed_log_every", "error"]
    with open(out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    ok = [r for r in rows if not r.get("error")]
    print(f"\nWrote {out.name}: {len(ok)} measured, {len(rows) - len(ok)} failed")
    if ok:
        print("Apply with:  python3 apply_lengths.py")


if __name__ == "__main__":
    main()
