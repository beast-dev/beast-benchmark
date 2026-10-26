#!/usr/bin/env python3
"""
make_report.py — build a one-page HTML summary of results.csv

Usage:
    python3 make_report.py [results.csv] [-o report.html]

The CSV rows are embedded in the page as JSON; all derived quantities
(throughput, speed-ups, groupings) are computed in the browser so the page
is self-contained apart from Plotly, which is loaded from jsDelivr.
"""

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
TEMPLATE    = _SCRIPT_DIR / "report_template.html"


def main() -> None:
    parser = argparse.ArgumentParser(description="Build an HTML report from results.csv")
    parser.add_argument("csv", nargs="?", default=str(_SCRIPT_DIR / "results.csv"),
                        help="Results CSV (default: results.csv beside this script)")
    parser.add_argument("-o", "--output", default=str(_SCRIPT_DIR / "report.html"),
                        help="Output HTML file (default: report.html beside this script)")
    args = parser.parse_args()

    with open(args.csv, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))

    # Keep '</' out of the inline script
    data = json.dumps(rows).replace("</", "<\\/")
    html = (TEMPLATE.read_text(encoding="utf-8")
            .replace("/*__DATA__*/[]", data)
            .replace("__GENERATED__", datetime.now().strftime("%Y-%m-%d %H:%M"))
            .replace("__SOURCE__", Path(args.csv).name))
    Path(args.output).write_text(html, encoding="utf-8")
    print(f"Wrote {args.output} ({len(rows)} runs)")


if __name__ == "__main__":
    main()
