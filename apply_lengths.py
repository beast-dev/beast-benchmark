#!/usr/bin/env python3
"""
Write the chain lengths proposed by calibrate.py into the benchmark XMLs.

Each XML carries the length in five places that have to agree: chainLength on
<mcmc>, logEvery on the screen log, the file log and the tree log, and the
"MCMC length" and "Logging every" lines of the header comment. This sets all of
them, keeping 1000 samples per run.

It also drops the "Final Joint density" line, because that value belongs to the
old chain length. Run the sweep again to replace it.

    python3 apply_lengths.py --dry-run     # show what would change
    python3 apply_lengths.py               # change it

calibrate.py only times the opening stretch of a chain, which turns out to
understate the speed of the BEAGLE benchmarks badly on the larger alignments.
--from-results replaces that estimate with the wall clock of a completed run,
rescaling each length by how far it landed from the target:

    python3 apply_lengths.py --from-results results.csv --dry-run
"""

import argparse
import csv
import re
import sys
from pathlib import Path

from benchmarks import xml_path
from calibrate import round_states

HEADER_BLOCK_RE = re.compile(r'(<!-->\s*)(.*?)(\s*<-->)', re.DOTALL)
CHAIN_LENGTH_RE = re.compile(r'(chainLength=")(\d+)(")')
LOG_EVERY_RE = re.compile(r'(logEvery=")(\d+)(")')
HEADER_LENGTH_RE = re.compile(r'^(MCMC length:\s*)([\d,]+)', re.MULTILINE)
HEADER_EVERY_RE = re.compile(r'^(Logging every:\s*)([\d,]+)', re.MULTILINE)
# Eat the newline before the line too, so removing the last line of the header
# does not leave a blank one behind.
FINAL_JOINT_LINE_RE = re.compile(r'\n?^Final Joint density:[^\n]*', re.MULTILINE)


def rewrite(text: str, states: int, log_every: int, drop_final: bool):
    """Return (new text, old chain length) or raise ValueError if the XML is not as expected."""
    if len(CHAIN_LENGTH_RE.findall(text)) != 1:
        raise ValueError("expected exactly one chainLength")
    if len(LOG_EVERY_RE.findall(text)) != 3:
        raise ValueError(f"expected three logEvery, found {len(LOG_EVERY_RE.findall(text))}")

    old_states = int(CHAIN_LENGTH_RE.search(text).group(2))

    text = CHAIN_LENGTH_RE.sub(lambda m: f"{m.group(1)}{states}{m.group(3)}", text)
    text = LOG_EVERY_RE.sub(lambda m: f"{m.group(1)}{log_every}{m.group(3)}", text)

    block_match = HEADER_BLOCK_RE.search(text)
    if not block_match:
        raise ValueError("no header comment block")

    block = block_match.group(2)
    if not HEADER_LENGTH_RE.search(block) or not HEADER_EVERY_RE.search(block):
        raise ValueError("header has no MCMC length / Logging every line")
    block = HEADER_LENGTH_RE.sub(lambda m: f"{m.group(1)}{states}", block)
    block = HEADER_EVERY_RE.sub(lambda m: f"{m.group(1)}{log_every}", block)
    if drop_final:
        block = FINAL_JOINT_LINE_RE.sub("", block)

    text = text[:block_match.start(2)] + block + text[block_match.end(2):]
    return text, old_states


def verify(text: str, states: int, log_every: int):
    """Every one of the five places agrees, or a list of complaints."""
    problems = []
    chain = [int(m[1]) for m in CHAIN_LENGTH_RE.findall(text)]
    every = [int(m[1]) for m in LOG_EVERY_RE.findall(text)]
    if chain != [states]:
        problems.append(f"chainLength {chain} != [{states}]")
    if every != [log_every] * 3:
        problems.append(f"logEvery {every} != {[log_every] * 3}")

    block_match = HEADER_BLOCK_RE.search(text)
    if not block_match:
        # An XML straight out of BEAUti, not yet given a benchmark header
        problems.append("no benchmark header block")
        return problems

    block = block_match.group(2)
    hl = HEADER_LENGTH_RE.search(block)
    he = HEADER_EVERY_RE.search(block)
    if not hl or int(hl.group(2).replace(",", "")) != states:
        problems.append("header MCMC length disagrees")
    if not he or int(he.group(2).replace(",", "")) != log_every:
        problems.append("header Logging every disagrees")
    return problems


def proposals_from_results(path: Path, target: float, tolerance: float) -> list:
    """Rescale each chain length by how far its last completed run was from the target.

    Returns rows shaped like calibrate.py's, so the caller does not care which
    source it got them from. The newest run of each benchmark wins; a run already
    within the tolerance is left alone, so a length only moves when it has to.
    """
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    newest = {}
    for row in rows:
        name = row["benchmark"]
        if name not in newest or row["timestamp"] > newest[name]["timestamp"]:
            newest[name] = row

    proposals = []
    for name, row in sorted(newest.items()):
        proposal = {"benchmark": name, "type": row["type"],
                    "number": row["number"], "model": row["model"], "error": ""}
        try:
            wall = float(row["wall_clock_time"])
            states = int(row["chain_length"])
            # BEAST's own runtime covers the chain alone, so the difference is the
            # startup that no chain length can shorten: hold it back from both sides.
            startup = max(0.0, wall - float(row["runtime"]) * 60.0)
        except (ValueError, KeyError):
            proposal["error"] = "no usable timing in the row"
            proposals.append(proposal)
            continue

        chain_time = wall - startup
        if chain_time <= 0:
            proposal["error"] = f"ran for {wall:.0f}s, all of it startup"
            proposals.append(proposal)
            continue

        proposal["ratio"] = wall / target
        if abs(wall - target) <= tolerance * target:
            proposal["error"] = f"within tolerance at {wall:.0f}s"
            proposals.append(proposal)
            continue

        scaled = round_states(states * (target - startup) / chain_time)
        proposal["proposed_states"] = scaled
        proposal["proposed_log_every"] = scaled // 1000
        proposals.append(proposal)

    return proposals


def check_all(benchmarks: Path) -> int:
    """Confirm every XML's five length fields agree and give 1000 samples."""
    xmls = sorted(benchmarks.glob("Benchmark-*/benchmark-*.xml"))
    if not xmls:
        print(f"No XMLs under {benchmarks}")
        return 1

    bad = 0
    for xml in xmls:
        text = xml.read_text(encoding="utf-8")
        chain = [int(m[1]) for m in CHAIN_LENGTH_RE.findall(text)]
        every = [int(m[1]) for m in LOG_EVERY_RE.findall(text)]
        problems = []
        if len(chain) != 1:
            problems.append(f"{len(chain)} chainLength")
        if len(every) != 3 or len(set(every)) != 1:
            problems.append(f"logEvery {every}")

        if not problems:
            states, log_every = chain[0], every[0]
            problems = verify(text, states, log_every)
            samples = states / log_every if log_every else 0
            if samples != 1000:
                problems.append(f"{samples:g} samples, not 1000")
            print(f"  {xml.stem:36s} {states:>12,}  every {log_every:>8,}"
                  f"  {'OK' if not problems else '; '.join(problems)}")
        else:
            print(f"  {xml.stem:36s} {'; '.join(problems)}")
        bad += bool(problems)

    print(f"\nChecked {len(xmls)}, {bad} with problems")
    return 1 if bad else 0


def main():
    parser = argparse.ArgumentParser(description="Apply calibrated chain lengths to the XMLs.")
    parser.add_argument("--check", action="store_true",
                        help="Only check that every XML's length fields agree; change nothing")
    parser.add_argument("--calibration", default="calibration.csv", metavar="FILE")
    parser.add_argument("--from-results", nargs="?", const="results.csv", metavar="FILE",
                        help="Rescale from the wall clock of completed runs in results.csv, "
                             "rather than from a calibrate.py estimate")
    parser.add_argument("--target", type=float, default=600.0, metavar="SECONDS",
                        help="Wall clock each benchmark should take (default 600)")
    parser.add_argument("--tolerance", type=float, default=0.2, metavar="FRACTION",
                        help="Leave a length alone if its run was this close to the target "
                             "(default 0.2)")
    parser.add_argument("--benchmarks", default=None, metavar="DIR")
    parser.add_argument("--dry-run", action="store_true", help="Show the changes, write nothing")
    parser.add_argument("--keep-final", action="store_true",
                        help="Leave the Final Joint density line alone (it will be stale)")
    args = parser.parse_args()

    here = Path(__file__).resolve().parent
    benchmarks = Path(args.benchmarks) if args.benchmarks else here / "Benchmarks"
    if args.check:
        sys.exit(check_all(benchmarks))

    if args.from_results:
        source = Path(args.from_results)
        if not source.exists():
            sys.exit(f"No {source} -- run the sweep first")
        rows = proposals_from_results(source, args.target, args.tolerance)
    else:
        source = Path(args.calibration)
        if not source.exists():
            sys.exit(f"No {source} -- run calibrate.py first")
        with open(source, encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))

    changed = skipped = failed = 0
    for row in rows:
        name = row["benchmark"]
        if row.get("error"):
            print(f"  {name:36s} SKIPPED ({row['error']})")
            skipped += 1
            continue

        states = int(row["proposed_states"])
        log_every = int(row["proposed_log_every"])
        xml = xml_path(benchmarks, row["type"], row["number"], row["model"])
        if not xml.exists():
            print(f"  {name:36s} MISSING {xml}")
            failed += 1
            continue

        text = xml.read_text(encoding="utf-8")
        try:
            new_text, old_states = rewrite(text, states, log_every, not args.keep_final)
        except ValueError as exc:
            print(f"  {name:36s} FAILED: {exc}")
            failed += 1
            continue

        problems = verify(new_text, states, log_every)
        if problems:
            print(f"  {name:36s} FAILED: {'; '.join(problems)}")
            failed += 1
            continue

        arrow = "would be" if args.dry_run else "->"
        print(f"  {name:36s} {old_states:>12,} {arrow} {states:>12,}  (logEvery {log_every:,})")
        if not args.dry_run:
            xml.write_text(new_text, encoding="utf-8")
        changed += 1

    verb = "would change" if args.dry_run else "changed"
    print(f"\n{verb} {changed}, skipped {skipped}, failed {failed}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
