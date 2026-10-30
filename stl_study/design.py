#!/usr/bin/env python3
"""
Generates the conditions table for the STL operator parameter study.

Three arms study the SubtreeLeap (STL) operator's effect on mixing efficiency:

  target_sweep  single STL, vary targetAcceptance; weight fixed at default (nTaxa)
  weight_sweep  single STL, vary weight relative to default (factor × nTaxa); target fixed at 0.234
  two_stl       two STLs: first at default, second with varied target and weight factor

The default run (single STL, target 0.234, default weight) is the shared control.
Reference runs are long default runs for comparison.

Arms:
  pilot          default runs for datasets with no chain length yet
  reference      long default runs giving the reference posterior
  target_sweep   the control plus each non-default target acceptance rate
  weight_sweep   the control plus each non-default weight factor
  two_stl        the control plus every (stl2_target, stl2_factor) pair
  all            target_sweep + weight_sweep + two_stl

Usage:
  design.py ARM [ARM ...] [-o conditions.csv] [--append]
"""

import argparse
import csv
import sys

from common import find_reference, format_level, load_config, work_dir

DEFAULT_TARGET = 0.234
DEFAULT_WEIGHT_FACTOR = 1.0

FIELDS = ["arm", "level", "run_id", "dataset", "likelihood", "rep", "seed", "states",
          "sl_target", "sl_weight_factor", "sl_weight",
          "stl2_target", "stl2_factor", "stl2_weight",
          # consumed by make_variants.py / run_grid.py but not set here
          "sl_size", "sl_adapt", "sl_kernel"]


def length_tag(states: int) -> str:
    return f"{states // 1000000}M" if states % 1000000 == 0 else f"{states // 1000}k"


def make_row(config, arm, level, dataset, likelihood, seed, states,
             sl_target=None, sl_weight_factor=None,
             stl2_target=None, stl2_factor=None):
    taxa = config["datasets"][dataset]["taxa"]
    row = {k: "" for k in FIELDS}
    row.update(arm=arm, level=level, dataset=dataset, likelihood=likelihood, rep=seed, seed=seed, states=states)

    parts = []
    if sl_target is not None and sl_target != DEFAULT_TARGET:
        row["sl_target"] = sl_target
        parts.append(f"t{format_level(sl_target)}")
    if sl_weight_factor is not None and sl_weight_factor != DEFAULT_WEIGHT_FACTOR:
        row["sl_weight_factor"] = sl_weight_factor
        row["sl_weight"] = f"{sl_weight_factor * taxa:.6g}"
        parts.append(f"w{format_level(sl_weight_factor)}")
    if stl2_target is not None:
        row.update(stl2_target=stl2_target, stl2_factor=stl2_factor,
                   stl2_weight=f"{stl2_factor * taxa:.6g}")
        parts.append(f"stl2t{format_level(stl2_target)}w{format_level(stl2_factor)}")

    sig = "_".join(parts) if parts else "default"
    row["run_id"] = f"{dataset}-{likelihood}_{length_tag(int(states))}_{sig}_s{seed}"
    return row


def study_states(config, dataset):
    states = config["datasets"][dataset].get("states")
    if not states:
        sys.exit(f"datasets.{dataset}.states is not set: run the pilot first")
    return states


def arm_rows(config, arm):
    datasets = config["study_datasets"]
    reps = config["reps"]
    rows = []

    if arm == "pilot":
        for ds in datasets:
            d = config["datasets"][ds]
            if d.get("states"):
                continue
            for lik in d["likelihoods"]:
                for seed in range(1, config["pilot_runs"] + 1):
                    rows.append(make_row(config, arm, "default", ds, lik, seed, config["pilot_states"]))

    elif arm == "reference":
        for ds in datasets:
            d = config["datasets"][ds]
            lik = d.get("reference_likelihood", "min")
            if find_reference(config, f"{ds}-{lik}"):
                print(f"{ds}: reference already found, skipping reference runs")
                continue
            states = study_states(config, ds) * d.get("reference_multiplier", config["reference_multiplier"])
            for seed in range(1001, 1001 + config["reference_runs"]):
                rows.append(make_row(config, arm, "default", ds, lik, seed, states))

    elif arm == "target_sweep":
        for ds in datasets:
            d = config["datasets"][ds]
            states = study_states(config, ds)
            targets = config["levels"]["stl_targets"]
            for lik in d["likelihoods"]:
                for seed in range(1, reps["default"] + 1):
                    rows.append(make_row(config, arm, "default", ds, lik, seed, states))
                for target in targets:
                    if target == DEFAULT_TARGET:
                        continue
                    for seed in range(1, reps["grid"] + 1):
                        rows.append(make_row(config, arm, f"t{format_level(target)}", ds, lik, seed, states,
                                             sl_target=target))

    elif arm == "weight_sweep":
        for ds in datasets:
            d = config["datasets"][ds]
            states = study_states(config, ds)
            factors = config["levels"]["stl_weight_factors"]
            for lik in d["likelihoods"]:
                for seed in range(1, reps["default"] + 1):
                    rows.append(make_row(config, arm, "default", ds, lik, seed, states))
                for factor in factors:
                    if factor == DEFAULT_WEIGHT_FACTOR:
                        continue
                    for seed in range(1, reps["grid"] + 1):
                        rows.append(make_row(config, arm, f"w{format_level(factor)}", ds, lik, seed, states,
                                             sl_weight_factor=factor))

    elif arm == "two_stl":
        for ds in datasets:
            d = config["datasets"][ds]
            states = study_states(config, ds)
            stl2_targets = config["levels"]["stl2_targets"]
            stl2_factors = config["levels"]["stl2_factors"]
            for lik in d["likelihoods"]:
                for seed in range(1, reps["default"] + 1):
                    rows.append(make_row(config, arm, "default", ds, lik, seed, states))
                for target in stl2_targets:
                    for factor in stl2_factors:
                        for seed in range(1, reps["grid"] + 1):
                            rows.append(make_row(config, arm,
                                                 f"stl2t{format_level(target)}w{format_level(factor)}",
                                                 ds, lik, seed, states,
                                                 stl2_target=target, stl2_factor=factor))

    elif arm == "all":
        for a in ["target_sweep", "weight_sweep", "two_stl"]:
            rows.extend(arm_rows(config, a))

    else:
        sys.exit(f"Unknown arm {arm!r}. Valid: pilot reference target_sweep weight_sweep two_stl all")

    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("arms", nargs="+")
    parser.add_argument("-o", "--output", help="default: WORK_DIR/conditions.csv")
    parser.add_argument("--append", action="store_true")
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
