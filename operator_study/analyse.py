#!/usr/bin/env python3
"""
Analyses completed study runs.

  analyse.py reference   pool the reference runs per dataset into reference/DATASET-LIK.json:
                         clade frequencies, a focal tree (the sampled tree with the highest
                         product of clade frequencies) and between-run ASDSF
  analyse.py runs        per-run metrics -> runs/RUN_ID/metrics.json and analysis/runs.csv
  analyse.py summary     per (arm, dataset, level) efficiency relative to the default
                         configuration -> analysis/summary.csv and analysis/*.png
  analyse.py all         all three

Efficiency is ESS per hour of BEAST run time (eff_*), and ESS per million states
(mst_*), which does not depend on how loaded the machine was, for:
  min_cont   the smallest ESS over the headline continuous parameters
  rf         Robinson-Foulds distance to the reference focal tree
  split_min  the smallest ESS over indicator traces of the uncertain clades
             (reference frequency in [0.05, 0.95], up to the 50 nearest 0.5)
and accuracy per cost as ASDSF against the reference x sqrt(hours).
"""

import argparse
import csv
import gzip
import json
import math
import random
import sys

import numpy as np

from common import (NON_PARAMETER_COLUMNS, clades_from_json, clades_to_json, ess,
                    geometric_mean, load_config, read_log, read_ops, read_trees, split_deviation,
                    split_frequencies, work_dir)

# set from the config in main()
RUNS_DIR = REFERENCE_DIR = ANALYSIS_DIR = None

HEADLINE_PARAMETERS = ["joint", "prior", "likelihood", "rootHeight", "treeLength",
                       "constant.popSize", "clock.rate", "kappa"]
MAX_SPLIT_INDICATORS = 50
# clades rarer than this are left out of the stored reference: the metrics only use clades at
# >= 1% (ASDSF) or 5-95% (split ESS), and keeping every clade of a diffuse tree posterior
# (tall-1) made the file over 500 MB
REFERENCE_MIN_FREQ = 0.001
METRICS = ["min_cont", "rf", "split_min"]


def load_conditions(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def load_result(run_id):
    try:
        result = json.loads((RUNS_DIR / run_id / "result.json").read_text())
    except (OSError, ValueError):
        return None
    return result if result.get("exit_code") == 0 else None


def post_burnin(seq, burnin):
    return seq[int(len(seq) * burnin):]


def reference_key(dataset, likelihood):
    # the MAT arm is scored against the BEAGLE reference: both target the same posterior
    return f"{dataset}-min"


# ---------------------------------------------------------------------------
# reference
# ---------------------------------------------------------------------------

def build_reference(config, conditions):
    REFERENCE_DIR.mkdir(exist_ok=True)
    by_dataset = {}
    for row in conditions:
        if row["arm"] == "reference":
            by_dataset.setdefault(row["dataset"], set()).add(row["run_id"])

    for dataset, run_ids in sorted(by_dataset.items()):
        done = [r for r in sorted(run_ids) if load_result(r)]
        if len(done) < 2:
            print(f"{dataset}: {len(done)} reference runs complete, need at least 2")
            continue
        clade_sets = []
        all_trees = []
        params = {}
        for rid in done:
            _, trees = read_trees(RUNS_DIR / rid / f"{rid}.trees")
            trees = post_burnin([c for _, c in trees], config["burnin"])
            clade_sets.append(trees)
            all_trees.extend(trees)
            columns, data = read_log(RUNS_DIR / rid / f"{rid}.log")
            data = post_burnin(data, config["burnin"])
            for i, col in enumerate(columns):
                if col not in NON_PARAMETER_COLUMNS:
                    params.setdefault(col, []).append(float(data[:, i].mean()))

        pooled = split_frequencies(clade_sets)
        per_run = [split_frequencies([t]) for t in clade_sets]
        pairwise = [split_deviation(per_run[i], per_run[j])[0]
                    for i in range(len(per_run)) for j in range(i + 1, len(per_run))]

        def credibility(clades):
            return sum(math.log(pooled[c]) for c in clades)
        focal = max(all_trees, key=credibility)

        out = {
            "dataset": dataset,
            "runs": done,
            "trees": len(all_trees),
            "max_pairwise_asdsf": max(pairwise),
            "parameter_means": {k: float(np.mean(v)) for k, v in params.items()},
            "parameter_run_spread": {k: float(np.std(v)) for k, v in params.items()},
            "focal_clades": [format(c, "x") for c in focal],
            "min_stored_frequency": REFERENCE_MIN_FREQ,
            "clades_seen": len(pooled),
            "clade_frequencies": clades_to_json({c: f for c, f in pooled.items() if f >= REFERENCE_MIN_FREQ}),
        }
        path = REFERENCE_DIR / f"{dataset}-min.json.gz"
        with gzip.open(path, "wt") as f:
            json.dump(out, f)
        (REFERENCE_DIR / f"{dataset}-min.json").unlink(missing_ok=True)  # the old uncompressed form
        flag = "" if out["max_pairwise_asdsf"] < 0.01 else "  <-- reference runs disagree, run longer"
        print(f"{dataset}: {len(done)} runs, {len(all_trees)} trees, "
              f"{len(out['clade_frequencies'])} of {len(pooled)} clades stored, "
              f"{path.stat().st_size / 1e6:.1f} MB, max pairwise ASDSF {out['max_pairwise_asdsf']:.4f}{flag}")


def load_reference(key):
    path = REFERENCE_DIR / f"{key}.json.gz"
    if not path.exists():
        return None
    with gzip.open(path, "rt") as f:
        data = json.load(f)
    data["clade_frequencies"] = clades_from_json(data["clade_frequencies"])
    data["focal_clades"] = {int(c, 16) for c in data["focal_clades"]}
    return data


# ---------------------------------------------------------------------------
# per-run metrics
# ---------------------------------------------------------------------------

def run_metrics(config, run_id, result, reference):
    run_dir = RUNS_DIR / run_id
    burnin = config["burnin"]
    hours = (result["beast_minutes"] or result["wall_seconds"] / 60) / 60
    m = {"run_id": run_id, "hours": hours}

    columns, data = read_log(run_dir / f"{run_id}.log")
    data = post_burnin(data, burnin)
    step = int(data[1, 0] - data[0, 0]) if len(data) > 1 else 1
    m["samples"] = len(data)
    for i, col in enumerate(columns):
        if col == "state":
            continue
        if col in ("topologyChanged", "mean_delta"):
            m[f"sl_{col}"] = float(np.nanmean(data[:, i]))
        else:
            m[f"ess_{col}"] = ess(data[:, i], step)
    headline = [m[f"ess_{p}"] for p in HEADLINE_PARAMETERS
                if f"ess_{p}" in m and not math.isnan(m[f"ess_{p}"])]
    m["ess_min_cont"] = min(headline) if headline else float("nan")
    m["ess_min_cont_parameter"] = min((p for p in HEADLINE_PARAMETERS if f"ess_{p}" in m
                                       and not math.isnan(m[f"ess_{p}"])),
                                      key=lambda p: m[f"ess_{p}"], default="")

    if reference:
        _, trees = read_trees(run_dir / f"{run_id}.trees")
        clade_sets = post_burnin([c for _, c in trees], burnin)
        focal = reference["focal_clades"]
        rf = [len(c ^ focal) for c in clade_sets]
        m["ess_rf"] = ess(rf, step)
        m["rf_mean"] = float(np.mean(rf))

        ref_freqs = reference["clade_frequencies"]
        uncertain = sorted((c for c, f in ref_freqs.items() if 0.05 <= f <= 0.95),
                           key=lambda c: abs(ref_freqs[c] - 0.5))[:MAX_SPLIT_INDICATORS]
        split_ess = []
        for clade in uncertain:
            indicator = [1.0 if clade in c else 0.0 for c in clade_sets]
            value = ess(indicator, step)
            # a clade never (or always) visited has not mixed: count it as a single sample
            split_ess.append(1.0 if math.isnan(value) else value)
        m["split_indicators"] = len(split_ess)
        m["ess_split_min"] = min(split_ess) if split_ess else float("nan")
        m["ess_split_median"] = float(np.median(split_ess)) if split_ess else float("nan")

        asdsf, max_dev = split_deviation(split_frequencies([clade_sets]), ref_freqs)
        m["asdsf"] = asdsf
        m["max_split_deviation"] = max_dev
        m["asdsf_x_sqrt_hours"] = asdsf * math.sqrt(hours)

    mstates = int(result["states"]) / 1e6
    for metric in METRICS:
        value = m.get(f"ess_{metric}")
        valid = value is not None and not math.isnan(value)
        m[f"eff_{metric}"] = value / hours if valid else float("nan")
        m[f"mst_{metric}"] = value / mstates if valid else float("nan")

    ops = read_ops(run_dir / f"{run_id}.ops")
    total_time = sum(o["time_s"] for o in ops.values())
    for name, o in ops.items():
        short = name.split("(")[0]
        if short in ("subtreeLeap", "fixedHeightSubtreePruneRegraft"):
            key = "sl" if short == "subtreeLeap" else "fhspr"
            m[f"{key}_accept"] = o["accept"]
            m[f"{key}_tuning"] = o["tuning"]
            m[f"{key}_us_per_op"] = o["time_per_op_us"]
            m[f"{key}_time_fraction"] = o["time_s"] / total_time if total_time else float("nan")
    m["operators"] = ops
    return m


def analyse_runs(config, conditions, force=False):
    rows = []
    references = {}
    run_info = {}
    for row in conditions:
        run_info.setdefault(row["run_id"], row)
    for run_id, cond in sorted(run_info.items()):
        result = load_result(run_id)
        if result is None or cond["arm"] == "reference":
            continue
        cache = RUNS_DIR / run_id / "metrics.json"
        key = reference_key(cond["dataset"], cond["likelihood"])
        if key not in references:
            references[key] = load_reference(key)
        reference = references[key]
        if cache.exists() and not force:
            m = json.loads(cache.read_text())
            if (reference and "ess_rf" not in m) or "mst_min_cont" not in m:
                m = run_metrics(config, run_id, result, reference)
        else:
            m = run_metrics(config, run_id, result, reference)
        cache.write_text(json.dumps(m, indent=1))
        rows.append(m)
        print(f"{run_id}: min ESS/h {m['eff_min_cont']:.0f} ({m['ess_min_cont_parameter']})"
              + (f", split-min ESS/h {m['eff_split_min']:.0f}" if "eff_split_min" in m else ""))

    ANALYSIS_DIR.mkdir(exist_ok=True)
    if rows:
        fields = sorted({k for r in rows for k in r if k != "operators"},
                        key=lambda k: (k != "run_id", k))
        with open(ANALYSIS_DIR / "runs.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for r in rows:
                writer.writerow({k: v for k, v in r.items() if k != "operators"})
    return {r["run_id"]: r for r in rows}


# ---------------------------------------------------------------------------
# summary
# ---------------------------------------------------------------------------

def bootstrap_ratio(values, baseline, n=2000, rng=random.Random(1)):
    """Geometric-mean ratio and 95% bootstrap interval, resampling runs within each cell."""
    values = [v for v in values if v and v > 0 and math.isfinite(v)]
    baseline = [v for v in baseline if v and v > 0 and math.isfinite(v)]
    if not values or not baseline:
        return float("nan"), float("nan"), float("nan")
    ratio = geometric_mean(values) / geometric_mean(baseline)
    draws = sorted(geometric_mean(rng.choices(values, k=len(values))) /
                   geometric_mean(rng.choices(baseline, k=len(baseline))) for _ in range(n))
    return ratio, draws[int(0.025 * n)], draws[int(0.975 * n) - 1]


def is_default(run_id):
    return "_default_s" in run_id


def run_prefix(run_id):
    return "_".join(run_id.split("_")[:2])


def summarise(conditions, metrics):
    cells = {}
    for row in conditions:
        if row["arm"] in ("reference", "pilot") or row["run_id"] not in metrics:
            continue
        cells.setdefault((row["arm"], row["dataset"], row["likelihood"], row["level"]), []).append(
            metrics[row["run_id"]])
    # baseline: default runs on the same dataset, likelihood and chain length ("small-1-min_20M")
    defaults = {}
    for run_id, m in metrics.items():
        if is_default(run_id):
            defaults.setdefault(run_prefix(run_id), []).append(m)

    rows = []
    for (arm, dataset, likelihood, level), runs in sorted(cells.items(), key=lambda kv: (
            kv[0][0], kv[0][1], kv[0][2], _numeric(kv[0][3]))):
        base = defaults.get(run_prefix(runs[0]["run_id"]), [])
        row = {"arm": arm, "dataset": dataset, "likelihood": likelihood, "level": level, "n": len(runs),
               "n_default": len(base)}
        for metric in METRICS:
            for key in (f"eff_{metric}", f"mst_{metric}"):
                row[f"{key}_gmean"] = geometric_mean([r.get(key) for r in runs])
                ratio, lo, hi = bootstrap_ratio([r.get(key) for r in runs], [r.get(key) for r in base])
                row[f"{key}_ratio"], row[f"{key}_lo"], row[f"{key}_hi"] = ratio, lo, hi
        for key in ("asdsf", "asdsf_x_sqrt_hours", "sl_accept", "sl_tuning", "sl_topologyChanged",
                    "sl_us_per_op", "fhspr_accept", "hours"):
            vals = [r[key] for r in runs if r.get(key) is not None and math.isfinite(r[key])]
            row[key] = float(np.mean(vals)) if vals else float("nan")
        rows.append(row)

    ANALYSIS_DIR.mkdir(exist_ok=True)
    if rows:
        with open(ANALYSIS_DIR / "summary.csv", "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"{len(rows)} cells -> {ANALYSIS_DIR / 'summary.csv'}")
        plot(rows, conditions, metrics)
    return rows


def _numeric(level):
    try:
        return (0, float(level), "")
    except ValueError:
        return (1, 0.0, level)


def plot(rows, conditions, metrics):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping plots")
        return

    for arm in sorted({r["arm"] for r in rows}):
        arm_rows = [r for r in rows if r["arm"] == arm]
        fig, axes = plt.subplots(1, len(METRICS), figsize=(5 * len(METRICS), 4), squeeze=False)
        for ax, metric in zip(axes[0], METRICS):
            for (dataset, likelihood) in sorted({(r["dataset"], r["likelihood"]) for r in arm_rows}):
                cell = [r for r in arm_rows if r["dataset"] == dataset and r["likelihood"] == likelihood
                        and _numeric(r["level"])[0] == 0]
                if not cell:
                    continue
                x = [float(r["level"]) for r in cell]
                y = [r[f"eff_{metric}_ratio"] for r in cell]
                err = [[max(0, r[f"eff_{metric}_ratio"] - r[f"eff_{metric}_lo"]) for r in cell],
                       [max(0, r[f"eff_{metric}_hi"] - r[f"eff_{metric}_ratio"]) for r in cell]]
                ax.errorbar(x, y, yerr=err, marker="o", capsize=3, label=f"{dataset}-{likelihood}")
            ax.axhline(1.0, color="grey", lw=0.8, ls="--")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel(f"{arm} level")
            ax.set_ylabel(f"ESS/h relative to default ({metric})")
            ax.set_title(metric)
        axes[0][0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(ANALYSIS_DIR / f"{arm}.png", dpi=120)
        plt.close(fig)

    # A1: efficiency against the realised SubtreeLeap acceptance, run by run
    a1 = [(row, metrics[row["run_id"]]) for row in conditions if row["arm"] == "A1" and row["run_id"] in metrics]
    if a1:
        fig, axes = plt.subplots(1, len(METRICS), figsize=(5 * len(METRICS), 4), squeeze=False)
        for ax, metric in zip(axes[0], METRICS):
            for dataset in sorted({row["dataset"] for row, _ in a1}):
                pts = [(m.get("sl_accept"), m.get(f"eff_{metric}")) for row, m in a1 if row["dataset"] == dataset]
                pts = [(a, e) for a, e in pts if a and e and math.isfinite(e)]
                if pts:
                    ax.scatter(*zip(*pts), label=dataset, s=14)
            ax.axvline(0.234, color="grey", lw=0.8, ls="--")
            ax.set_yscale("log")
            ax.set_xlabel("SubtreeLeap acceptance")
            ax.set_ylabel(f"ESS/h ({metric})")
        axes[0][0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(ANALYSIS_DIR / "A1_acceptance.png", dpi=120)
        plt.close(fig)
    print(f"plots -> {ANALYSIS_DIR}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["reference", "runs", "summary", "all"])
    parser.add_argument("conditions", nargs="?", help="default: WORK_DIR/conditions.csv")
    parser.add_argument("--force", action="store_true", help="recompute cached per-run metrics")
    args = parser.parse_args()

    global RUNS_DIR, REFERENCE_DIR, ANALYSIS_DIR
    config = load_config()
    RUNS_DIR = work_dir(config) / "runs"
    # the reference posteriors can be shared with another study ("reference_dir" in study.json)
    REFERENCE_DIR = work_dir(config) / config.get("reference_dir", "reference")
    ANALYSIS_DIR = work_dir(config) / "analysis"
    conditions = load_conditions(args.conditions or work_dir(config) / "conditions.csv")
    if args.command in ("reference", "all"):
        build_reference(config, conditions)
    if args.command in ("runs", "summary", "all"):
        metrics = analyse_runs(config, conditions, force=args.force)
        if args.command != "runs":
            summarise(conditions, metrics)


if __name__ == "__main__":
    sys.exit(main())
