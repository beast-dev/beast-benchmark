#!/usr/bin/env python3
"""
Write the Final Joint density of a completed run into the benchmark XML header.

benchmark.py compares the density a run ends on against the one recorded in the
header, which is how a change in BEAST that alters the chain shows up. The value
is therefore only meaningful for the exact chain it came from: a different chain
length, or a different seed, ends somewhere else entirely.

So a row is only used when its chain length and seed still match the XML. A
benchmark whose length has moved since its last run is reported and left alone,
to be filled in after the next sweep.

    python3 apply_densities.py --dry-run     # show what would be written
    python3 apply_densities.py               # write it
"""

import argparse
import csv
import re
import sys
from pathlib import Path

from benchmarks import xml_path
from apply_lengths import HEADER_BLOCK_RE, CHAIN_LENGTH_RE, FINAL_JOINT_LINE_RE

HEADER_SEED_RE = re.compile(r'^Seed:\s*(\d+)', re.MULTILINE)
INITIAL_LINE_RE = re.compile(r'^Initial Joint density:[^\n]*', re.MULTILINE)


def newest_rows(path: Path) -> dict:
    """The run to believe for each benchmark: the most recent, preferring CPU over GPU.

    A GPU run reaches a different final density -- the arithmetic differs -- so it
    must not overwrite the CPU value that the header is there to protect.
    """
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    best = {}
    for row in rows:
        name = row["benchmark"]
        on_gpu = "GPU" in (row.get("beagle_resource") or "")
        key = (not on_gpu, row["timestamp"])
        if name not in best or key > best[name][0]:
            best[name] = (key, row)
    return {name: row for name, (_, row) in best.items()}


def set_density(text: str, density: str):
    """Put the Final Joint density line in the header, replacing any line already there."""
    block_match = HEADER_BLOCK_RE.search(text)
    if not block_match:
        raise ValueError("no header comment block")

    block = block_match.group(2)
    line = f"Final Joint density: {density}"
    if FINAL_JOINT_LINE_RE.search(block):
        block = FINAL_JOINT_LINE_RE.sub("\n" + line, block)
    elif INITIAL_LINE_RE.search(block):
        block = INITIAL_LINE_RE.sub(lambda m: m.group(0) + "\n" + line, block)
    else:
        raise ValueError("header has no Initial Joint density line to follow")

    return text[:block_match.start(2)] + block + text[block_match.end(2):]


def main():
    parser = argparse.ArgumentParser(description="Record each run's final joint density in its XML.")
    parser.add_argument("--results", default="results.csv", metavar="FILE")
    parser.add_argument("--benchmarks", default=None, metavar="DIR")
    parser.add_argument("--dry-run", action="store_true", help="Show the changes, write nothing")
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    benchmarks = Path(args.benchmarks) if args.benchmarks else here / "Benchmarks"

    results = Path(args.results)
    if not results.exists():
        sys.exit(f"No {results} -- run the sweep first")

    written = stale = failed = 0
    for name, row in sorted(newest_rows(results).items()):
        xml = xml_path(benchmarks, row["type"], row["number"], row["model"])
        if not xml.exists():
            print(f"  {name:32s} MISSING {xml}")
            failed += 1
            continue

        text = xml.read_text(encoding="utf-8")
        raw = row.get("final_joint") or ""
        if not raw:
            print(f"  {name:32s} NO DENSITY in the run")
            failed += 1
            continue
        # BEAST prints four decimal places; the CSV drops trailing zeros
        density = f"{float(raw):.4f}"

        xml_states = int(CHAIN_LENGTH_RE.search(text).group(2))
        run_states = int(row["chain_length"])
        if xml_states != run_states:
            print(f"  {name:32s} STALE: ran {run_states:,}, XML is now {xml_states:,}")
            stale += 1
            continue

        seed_match = HEADER_SEED_RE.search(text)
        if seed_match and seed_match.group(1) != row["seed"]:
            print(f"  {name:32s} STALE: ran seed {row['seed']}, header says {seed_match.group(1)}")
            stale += 1
            continue

        try:
            new_text = set_density(text, density)
        except ValueError as exc:
            print(f"  {name:32s} FAILED: {exc}")
            failed += 1
            continue

        arrow = "would be" if args.dry_run else "->"
        print(f"  {name:32s} {arrow} {density}")
        if not args.dry_run:
            xml.write_text(new_text, encoding="utf-8")
        written += 1

    verb = "would write" if args.dry_run else "wrote"
    print(f"\n{verb} {written}, stale {stale}, failed {failed}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
