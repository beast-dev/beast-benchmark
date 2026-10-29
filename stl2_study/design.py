#!/usr/bin/env python3
"""
Generates the conditions table for the second-SubtreeLeap study.

Every run keeps the benchmark's operators as they are: SubtreeLeap (weight = number
of taxa, adapting to 0.234) and FHSPR (weight = taxa / 10). A grid run adds a second
SubtreeLeap (id stl2) with its own target acceptance and a weight of factor x taxa.
The control is the unchanged default.

Each row is one (arm, run) membership. The run id is derived from the effective
settings, so the default runs carry the same ids as in the other studies and can be
copied in from them when they were made with the same jar.

Arms:
  pilot      default runs for datasets with no chain length yet (study.json "states": null)
  reference  long default runs giving the reference posterior, for datasets that have none
  grid       the control plus every (stl2 target, stl2 factor) pair, per dataset and likelihood
  confirm    the control against the settings in study.json "confirm", on the datasets and
             likelihoods in "confirm_datasets" (default: every study dataset and likelihood)

Usage:
  design.py ARM [ARM ...] [-o conditions.csv] [--append]
"""

import argparse
import csv
import sys

from common import find_reference, format_level, load_config, work_dir

FIELDS = ["arm", "level", "run_id", "dataset", "likelihood", "rep", "seed", "states",
          "stl2_target", "stl2_factor", "stl2_weight",
          # read by make_variants.py / run_grid.py; unused here
          "sl_size", "sl_target", "sl_adapt", "sl_kernel", "adaptation_target", "operator_weights"]


def length_tag(states: int) -> str:
    return f"{states // 1000000}M" if states % 1000000 == 0 else f"{states // 1000}k"


def make_row(config, arm, level, dataset, likelihood, seed, states, tag="", stl2_target=None, stl2_factor=None):
    row = {k: "" for k in FIELDS}
    row.update(arm=arm, level=level, dataset=dataset, likelihood=likelihood, rep=seed, seed=seed, states=states)
    if stl2_target is not None:
        taxa = config["datasets"][dataset]["taxa"]
        row.update(stl2_target=stl2_target, stl2_factor=stl2_factor, stl2_weight=f"{stl2_factor * taxa:.6g}")
        sig = f"stl2t{format_level(stl2_target)}w{format_level(stl2_factor)}"
    else:
        sig = tag or "default"
    row["run_id"] = f"{dataset}-{likelihood}_{length_tag(int(states))}_{sig}_s{seed}"
    return row


def study_states(config, dataset):
    states = config["datasets"][dataset].get("states")
    if not states:
        sys.exit(f"datasets.{dataset}.states is not set: run the pilot and set it in study.json first")
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
                print(f"{ds}: reference {ds}-{lik} already exists, no reference runs needed")
                continue
            states = study_states(config, ds) * d.get("reference_multiplier", config["reference_multiplier"])
            for seed in range(1001, 1001 + config["reference_runs"]):
                rows.append(make_row(config, arm, "default", ds, lik, seed, states, tag="ref"))

    elif arm == "grid":
        for ds in datasets:
            d = config["datasets"][ds]
            states = study_states(config, ds)
            for lik in d["likelihoods"]:
                levels = dict(config["levels"], **d.get("levels", {}).get(lik, {}))
                n = d.get("reps", {}).get(lik, reps["grid"])
                for seed in range(1, d.get("reps_default", {}).get(lik, reps["default"]) + 1):
                    rows.append(make_row(config, arm, "default", ds, lik, seed, states))
                for target in levels["stl2_targets"]:
                    for factor in levels["stl2_factors"]:
                        for seed in range(1, n + 1):
                            rows.append(make_row(config, arm, f"t{target}_w{factor}", ds, lik, seed, states,
                                                 stl2_target=target, stl2_factor=factor))

    elif arm == "confirm":
        if not config.get("confirm"):
            sys.exit('Fill in study.json "confirm" ({"name": {"stl2_target": t, "stl2_factor": f}}) first')
        only = config.get("confirm_datasets")
        for ds in datasets:
            if only is not None and ds not in only:
                continue
            d = config["datasets"][ds]
            states = study_states(config, ds)
            for lik in (only[ds] if only is not None else d["likelihoods"]):
                for seed in range(1, reps["confirm"] + 1):
                    rows.append(make_row(config, arm, "default", ds, lik, seed, states))
                    for name, s in config["confirm"].items():
                        rows.append(make_row(config, arm, name, ds, lik, seed, states,
                                             stl2_target=s["stl2_target"], stl2_factor=s["stl2_factor"]))
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
