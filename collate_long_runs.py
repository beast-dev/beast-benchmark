#!/usr/bin/env python3
"""
Collates the longest default-settings runs of each benchmark XML into archive/.

Searches the study run directories (*_study/runs in this checkout and in sibling
beast-benchmark* worktrees, or the directories given) for completed runs with
default settings: the reference runs (`..._ref_sN`) and the default runs
(`..._default_sN`). For each dataset and likelihood it keeps every seed at the
longest chain length found and hard-links (or copies, across file systems) the
run's XML, .log, .trees, .ops and result.json to

  archive/{dataset}-{likelihood}/{run_id}/

and writes archive/MANIFEST.tsv. Hard links take no extra disk space and keep
the files if the study directories are later cleared. archive/ is not committed.
Re-running adds new runs (e.g. references still in progress) and skips existing ones.

Usage:
  collate_long_runs.py [RUNS_DIR ...] [--archive DIR] [--dry-run]
"""

import argparse
import json
import os
import re
import shutil
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DEFAULT_RUN = re.compile(r"_(ref|default)_s\d+$")
SUFFIXES = (".xml", ".log", ".trees", ".ops")


def default_runs_dirs():
    dirs = []
    for repo in [ROOT] + sorted(p for p in ROOT.parent.glob("beast-benchmark*") if p.resolve() != ROOT):
        dirs.extend(sorted(repo.glob("*_study/runs")))
    return dirs


def completed_runs(runs_dir):
    for run_dir in sorted(runs_dir.iterdir()):
        if not DEFAULT_RUN.search(run_dir.name):
            continue
        try:
            result = json.loads((run_dir / "result.json").read_text())
        except (OSError, ValueError):
            continue
        if result.get("exit_code") != 0 or not (run_dir / f"{run_dir.name}.trees").exists():
            continue
        yield run_dir, result


def link_or_copy(src, dst):
    if dst.exists():
        return
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs_dirs", nargs="*", type=Path)
    parser.add_argument("--archive", type=Path, default=ROOT / "archive")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    runs_dirs = args.runs_dirs or default_runs_dirs()
    for d in runs_dirs:
        print(f"searching {d}")

    # (dataset, likelihood) -> states -> {run_id: (run_dir, result)}; a run id seen in two
    # checkouts (e.g. a reference copied between studies) is kept once
    found = defaultdict(lambda: defaultdict(dict))
    for runs_dir in runs_dirs:
        for run_dir, result in completed_runs(runs_dir):
            key = (result["dataset"], result["likelihood"])
            found[key][int(result["states"])].setdefault(run_dir.name, (run_dir, result))

    rows = []
    for (dataset, likelihood), by_states in sorted(found.items()):
        states = max(by_states)
        runs = by_states[states]
        print(f"{dataset}-{likelihood}: {len(runs)} run(s) at {states // 1_000_000}M states")
        for run_id, (run_dir, result) in sorted(runs.items()):
            dest = args.archive / f"{dataset}-{likelihood}" / run_id
            files = [run_dir / f"{run_id}{s}" for s in SUFFIXES] + [run_dir / "result.json"]
            files = [f for f in files if f.exists()]
            size = sum(f.stat().st_size for f in files)
            print(f"  {run_id}  {size / 1e6:,.0f} MB  <- {run_dir}")
            if not args.dry_run:
                dest.mkdir(parents=True, exist_ok=True)
                for f in files:
                    link_or_copy(f, dest / f.name)
            rows.append([dataset, likelihood, run_id, result["seed"], str(states),
                         f"{result.get('beast_minutes', 0) / 60:.1f}", result.get("finished", ""),
                         result.get("jar_sha256", ""), str(run_dir)])

    if not args.dry_run and rows:
        header = ["dataset", "likelihood", "run_id", "seed", "states", "hours", "finished", "jar_sha256", "source"]
        args.archive.mkdir(parents=True, exist_ok=True)
        (args.archive / "MANIFEST.tsv").write_text("\n".join("\t".join(r) for r in [header] + rows) + "\n")
        print(f"{len(rows)} runs -> {args.archive}")


if __name__ == "__main__":
    main()
