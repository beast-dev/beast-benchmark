#!/usr/bin/env python3
"""
Where the benchmark XMLs are and what they are called.

Shared by calibrate.py and apply_lengths.py, and run as a script by the run_all
shell scripts so that adding or removing a benchmark XML needs no edit to them:

    python3 benchmarks.py          # one "type number model" per line
"""

import re
import sys
from pathlib import Path

# Dataset categories, and model labels in the order the runners use them.
TYPES = ["small", "big", "tall", "wide", "medium"]
MODELS = ["min", "mat", "gamma", "gamma-mat", "skygrid", "skygrid-mat"]

BENCHMARKS = Path(__file__).resolve().parent / "Benchmarks"


def xml_path(benchmarks: Path, btype: str, number, model: str) -> Path:
    return benchmarks / f"Benchmark-{btype}" / f"benchmark-{btype}-{number}-{model}.xml"


def discover(benchmarks: Path = None):
    """Every (type, number, model) that has an XML, grouped by dataset."""
    benchmarks = benchmarks or BENCHMARKS
    found = []
    for btype in TYPES:
        numbers = set()
        for path in (benchmarks / f"Benchmark-{btype}").glob("benchmark-*.xml"):
            match = re.match(rf'benchmark-{btype}-(\d+)-(.+)\.xml$', path.name)
            if match:
                numbers.add(int(match.group(1)))
        for number in sorted(numbers):
            for model in MODELS:
                if xml_path(benchmarks, btype, number, model).exists():
                    found.append((btype, number, model))
    return found


if __name__ == "__main__":
    benchmarks = Path(sys.argv[1]) if len(sys.argv) > 1 else BENCHMARKS
    for btype, number, model in discover(benchmarks):
        print(btype, number, model)
