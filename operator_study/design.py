#!/usr/bin/env python3
"""
Generates the conditions table for one or more study arms.

Each row is one (arm, run) membership. The run id is derived from the effective
settings, so the default configuration shared by several arms (A2 at 0.234,
A3 at 0.234, B1 at x1, B2 at 0.1, the pilot, confirmation) is run once per seed.

Arms:
  reference  long default runs giving the reference split frequencies
  pilot      default runs used to set chain lengths and the A1 size grid
  A1         SubtreeLeap fixed size (adaptation off), around the pilot's adapted size
  A2         SubtreeLeap targetAcceptance
  A3         -adaptation_target for every adaptive operator
  B1         -operator_weights multiplier on SubtreeLeap and FHSPR
  B2         FHSPR weight as a fraction of SubtreeLeap's
  confirm    default against the configurations in study.json "confirm"
  mat        A2 and B1 on the MAT variants

Usage:
  design.py ARM [ARM ...] [-o conditions.csv]
"""

import argparse
import csv
import sys

from common import format_level, load_config, work_dir

FIELDS = ["arm", "level", "run_id", "dataset", "likelihood", "rep", "seed", "states",
          "sl_size", "sl_target", "sl_adapt", "sl_kernel", "adaptation_target", "operator_weights"]

DEFAULT_TARGET = 0.234
FHSPR_DEFAULT_FRACTION = 0.1  # the benchmark XMLs give FHSPR a tenth of SubtreeLeap's weight
TREE_OPERATORS = ("subtreeLeap", "fixedHeightSubtreePruneRegraft")


def weights_spec(factors: dict) -> str:
    return ",".join(f"{k}={v:g}" for k, v in factors.items() if v != 1)


def run_id(row: dict) -> str:
    parts = []
    if row["states_tag"]:
        parts.append(row["states_tag"])
    if row["sl_adapt"] == "false":
        parts.append(f"slsize{format_level(row['sl_size'])}")
    if row["sl_target"]:
        parts.append(f"sltarget{format_level(row['sl_target'])}")
    if row["sl_kernel"]:
        parts.append(f"slkernel{row['sl_kernel']}")
    if row["adaptation_target"]:
        parts.append(f"target{format_level(row['adaptation_target'])}")
    if row["operator_weights"]:
        parts.append(f"w{format_level(row['operator_weights'])}")
    sig = "_".join(parts) or "default"
    states = int(row["states"])
    length = f"{states // 1000000}M" if states % 1000000 == 0 else f"{states // 1000}k"
    return f"{row['dataset']}-{row['likelihood']}_{length}_{sig}_s{row['seed']}"


def make_row(config, arm, level, dataset, likelihood, rep, states_tag="", states=None, **settings):
    row = {
        "arm": arm, "level": level, "dataset": dataset, "likelihood": likelihood,
        "rep": rep, "seed": rep, "states": states or config["datasets"][dataset]["states"],
        "states_tag": states_tag,
        "sl_size": "", "sl_target": "", "sl_adapt": "", "sl_kernel": "",
        "adaptation_target": "", "operator_weights": "",
    }
    for key, value in settings.items():
        row[key] = "" if value is None else value
    # normalise settings that equal the defaults so identical runs share an id
    if row["sl_target"] != "" and float(row["sl_target"]) == DEFAULT_TARGET:
        row["sl_target"] = ""
    if row["adaptation_target"] != "" and float(row["adaptation_target"]) == DEFAULT_TARGET:
        row["adaptation_target"] = ""
    row["run_id"] = run_id(row)
    del row["states_tag"]
    return row


def arm_rows(config, arm):
    levels = config["levels"]
    reps = config["reps"]
    datasets = config["study_datasets"]
    rows = []

    if arm == "reference":
        for ds in datasets:
            states = config["datasets"][ds]["states"] * config["reference_multiplier"]
            for rep in range(1, config["reference_runs"] + 1):
                rows.append(make_row(config, arm, "default", ds, "min", 1000 + rep, states_tag="ref", states=states))

    elif arm == "pilot":
        for ds in datasets:
            for rep in range(1, config["pilot_runs"] + 1):
                rows.append(make_row(config, arm, "default", ds, "min", rep))

    elif arm == "A1":
        for ds in datasets:
            size = config["datasets"][ds]["sl_size"]
            if not size:
                sys.exit(f"A1 needs datasets.{ds}.sl_size in study.json (the pilot's adapted SubtreeLeap size)")
            for mult in levels["A1_size_multipliers"]:
                for rep in range(1, reps["A1"] + 1):
                    rows.append(make_row(config, arm, mult, ds, "min", rep,
                                         sl_size=f"{size * mult:.6g}", sl_adapt="false"))

    elif arm == "A2":
        for ds in datasets:
            for target in levels["A2_targets"]:
                for rep in range(1, reps["A2"] + 1):
                    rows.append(make_row(config, arm, target, ds, "min", rep, sl_target=target))

    elif arm == "A3":
        for ds in datasets:
            for target in levels["A3_targets"]:
                for rep in range(1, reps["A3"] + 1):
                    rows.append(make_row(config, arm, target, ds, "min", rep, adaptation_target=target))

    elif arm in ("B1", "mat"):
        likelihood = "mat" if arm == "mat" else "min"
        for ds in (config["mat_datasets"] if arm == "mat" else datasets):
            for mult in levels["B1_tree_multipliers"]:
                for rep in range(1, reps[arm] + 1):
                    rows.append(make_row(config, "B1" if arm == "B1" else "mat-B1", mult, ds, likelihood, rep,
                                         operator_weights=weights_spec({op: mult for op in TREE_OPERATORS})))
            if arm == "mat":
                for target in levels["A2_targets"]:
                    for rep in range(1, reps["mat"] + 1):
                        rows.append(make_row(config, "mat-A2", target, ds, likelihood, rep, sl_target=target))

    elif arm == "B2":
        for ds in datasets:
            for fraction in levels["B2_fhspr_fractions"]:
                factor = fraction / FHSPR_DEFAULT_FRACTION
                for rep in range(1, reps["B2"] + 1):
                    rows.append(make_row(config, arm, fraction, ds, "min", rep,
                                         operator_weights=weights_spec({TREE_OPERATORS[1]: factor})))

    elif arm == "confirm":
        configs = {"default": {}}
        configs.update({k: v for k, v in config["confirm"].items() if v})
        if len(configs) == 1:
            sys.exit('Fill in study.json "confirm" with the best settings from A and B first')
        for ds in datasets:
            for name, settings in configs.items():
                for rep in range(1, reps["confirm"] + 1):
                    rows.append(make_row(config, arm, name, ds, "min", rep, **settings))
    else:
        sys.exit(f"Unknown arm {arm}")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("arms", nargs="+")
    parser.add_argument("-o", "--output", help="default: WORK_DIR/conditions.csv")
    parser.add_argument("--append", action="store_true", help="add to an existing conditions file")
    args = parser.parse_args()

    config = load_config()
    if args.output is None:
        work_dir(config).mkdir(parents=True, exist_ok=True)
        args.output = str(work_dir(config) / "conditions.csv")
    rows = [r for arm in args.arms for r in arm_rows(config, arm)]

    existing = []
    if args.append:
        try:
            with open(args.output) as f:
                existing = list(csv.DictReader(f))
        except FileNotFoundError:
            pass
    seen = {(r["arm"], r["run_id"]) for r in existing}
    rows = existing + [r for r in rows if (r["arm"], r["run_id"]) not in seen]

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    unique = len({r["run_id"] for r in rows})
    print(f"{len(rows)} arm memberships, {unique} distinct runs -> {args.output}")


if __name__ == "__main__":
    main()
