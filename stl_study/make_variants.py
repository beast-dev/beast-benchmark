#!/usr/bin/env python3
"""
Writes a study variant of a benchmark XML for one condition.

The edits are textual so the benchmark header and formatting survive:
  - chainLength, and logEvery so each run logs `samples` states
  - log/trees/ops file names set to the run id
  - subtreeLeap: id, size, weight, targetAcceptance, autoOptimize, distanceKernel
  - subtreeLeap added to the file log (topologyChanged, mean_delta columns)
  - optionally a second subtreeLeap (id stl2): a copy of the first, children (the MAT
    resampler hook) included, with its own weight and target acceptance

Usage (to inspect a variant):
  make_variants.py conditions.csv RUN_ID [-o out.xml]
"""

import argparse
import csv
import re
import sys

from common import base_xml_path, load_config

SL_ID = "subtreeLeap"
SL2_ID = "stl2"


def _sub_once(pattern, repl, text, what):
    new_text, n = re.subn(pattern, repl, text, count=1)
    if n != 1:
        raise ValueError(f"Could not find {what} in the base XML")
    return new_text


def _set_attribute(tag: str, name: str, value) -> str:
    """Sets (or adds) an attribute on an opening tag string."""
    attr = f'{name}="{value}"'
    if re.search(rf'\s{name}="[^"]*"', tag):
        return re.sub(rf'\s{name}="[^"]*"', " " + attr, tag)
    return re.sub(r"\s*(/?>)$", f" {attr}\\1", tag)


def make_variant(text: str, cond: dict, run_id: str, samples: int) -> str:
    states = int(cond["states"])
    log_every = max(1, states // samples)
    screen_every = max(1, states // 100)

    text = _sub_once(r'(<mcmc id="mcmc"[^>]*?chainLength=")\d+(")', rf"\g<1>{states}\2", text, "mcmc chainLength")
    text = _sub_once(r'(<mcmc id="mcmc"[^>]*?operatorAnalysis=")[^"]+(")', rf"\g<1>{run_id}.ops\2", text,
                     "mcmc operatorAnalysis")
    text = _sub_once(r'(<log id="screenLog" logEvery=")\d+(")', rf"\g<1>{screen_every}\2", text, "screenLog")
    text = _sub_once(r'(<log id="fileLog" logEvery=")\d+(" fileName=")[^"]+(")',
                     rf"\g<1>{log_every}\g<2>{run_id}.log\3", text, "fileLog")
    text = _sub_once(r'(<logTree id="treeFileLog" logEvery=")\d+("[^>]*?fileName=")[^"]+(")',
                     rf"\g<1>{log_every}\g<2>{run_id}.trees\3", text, "treeFileLog")

    m = re.search(r"<subtreeLeap\b[^>]*>", text)
    if not m:
        raise ValueError("Could not find subtreeLeap in the base XML")
    tag = _set_attribute(m.group(0), "id", SL_ID)
    if cond.get("sl_size"):
        tag = _set_attribute(tag, "size", cond["sl_size"])
    if cond.get("sl_weight"):
        tag = _set_attribute(tag, "weight", cond["sl_weight"])
    if cond.get("sl_target"):
        tag = _set_attribute(tag, "targetAcceptance", cond["sl_target"])
    if cond.get("sl_adapt") == "false":
        tag = _set_attribute(tag, "autoOptimize", "false")
    if cond.get("sl_kernel"):
        tag = _set_attribute(tag, "distanceKernel", cond["sl_kernel"])
    text = text[:m.start()] + tag + text[m.end():]

    if cond.get("stl2_target"):
        # a second SubtreeLeap: the same element, children included, under its own id, weight and
        # target; it is not added to the file log, whose topologyChanged columns would then clash
        start = m.start()
        end = text.index("</subtreeLeap>", start) + len("</subtreeLeap>")
        block = text[start:end]
        tag2 = re.match(r"<subtreeLeap\b[^>]*>", block).group(0)
        tag2 = _set_attribute(tag2, "id", SL2_ID)
        tag2 = _set_attribute(tag2, "weight", cond["stl2_weight"])
        tag2 = _set_attribute(tag2, "targetAcceptance", cond["stl2_target"])
        indent = text[text.rindex("\n", 0, start) + 1:start]
        text = text[:end] + "\n" + indent + tag2 + block[len(tag):] + text[end:]

    # log the operator's rolling topology-change fraction and mean leap distance
    start = text.index('<log id="fileLog"')
    end = text.index("</log>", start)
    text = text[:end] + f'\t<subtreeLeap idref="{SL_ID}"/>\n\t\t' + text[end:]
    return text


def write_variant(config: dict, cond: dict, out_path) -> None:
    base = base_xml_path(config, cond["dataset"], cond["likelihood"]).read_text()
    out_path.write_text(make_variant(base, cond, cond["run_id"], config["samples"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("conditions")
    parser.add_argument("run_id")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()

    config = load_config()
    with open(args.conditions) as f:
        cond = next((r for r in csv.DictReader(f) if r["run_id"] == args.run_id), None)
    if cond is None:
        sys.exit(f"No condition with run_id {args.run_id}")
    base = base_xml_path(config, cond["dataset"], cond["likelihood"]).read_text()
    text = make_variant(base, cond, cond["run_id"], config["samples"])
    if args.output:
        with open(args.output, "w") as f:
            f.write(text)
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
