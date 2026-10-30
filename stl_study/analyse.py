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
A run is stuck when some uncertain clade's indicator never changes (split ESS < 2); each cell
counts its stuck runs (n_stuck, n_stuck_default) with a one-sided Fisher exact p (fewer_stuck_p).
"""

import argparse
import csv
import gzip
import json
import math
import random
import sys

import numpy as np

from common import (NON_PARAMETER_COLUMNS, clades_from_json, clades_to_json, ess, find_reference,
                    geometric_mean, load_config, operator_kind, read_log, read_ops, read_trees, split_deviation,
                    split_frequencies, split_operator_name, work_dir)
from make_variants import SL_ID  # the id make_variants gives the main SubtreeLeap

# set from the config in main()
RUNS_DIR = REFERENCE_DIR = ANALYSIS_DIR = None

HEADLINE_PARAMETERS = ["joint", "prior", "likelihood", "rootHeight", "treeLength",
                       "constant.popSize", "clock.rate", "kappa"]
MAX_SPLIT_INDICATORS = 50
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


def reference_key(config, dataset):
    return f"{dataset}-{config['datasets'][dataset].get('reference_likelihood', 'min')}"


# ---------------------------------------------------------------------------
# reference
# ---------------------------------------------------------------------------

def build_reference(config, conditions):
    REFERENCE_DIR.mkdir(exist_ok=True)
    by_dataset = {}
    for row in conditions:
        if row["arm"] == "reference":
            by_dataset.setdefault((row["dataset"], row["likelihood"]), set()).add(row["run_id"])

    for (dataset, likelihood), run_ids in sorted(by_dataset.items()):
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
        path = REFERENCE_DIR / f"{dataset}-{likelihood}.json.gz"
        with gzip.open(path, "wt") as f:
            json.dump(out, f)
        flag = "" if out["max_pairwise_asdsf"] < 0.01 else "  <-- reference runs disagree, run longer"
        print(f"{dataset}: {len(done)} runs, {len(all_trees)} trees, "
              f"{len(out['clade_frequencies'])} of {len(pooled)} clades stored, "
              f"{path.stat().st_size / 1e6:.1f} MB, max pairwise ASDSF {out['max_pairwise_asdsf']:.4f}{flag}")


def load_reference(config, key):
    path = find_reference(config, key)
    if path is None:
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
        op_id, _ = split_operator_name(name)
        short = operator_kind(name)
        if short in ("subtreeLeap", "fixedHeightSubtreePruneRegraft"):
            key = "sl" if short == "subtreeLeap" else "fhspr"
            if op_id and op_id != SL_ID:
                key += f"_{op_id}"  # a second SubtreeLeap (or FHSPR) gets its own columns
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
        key = reference_key(config, cond["dataset"])
        if key not in references:
            references[key] = load_reference(config, key)
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
    values = [v for v in values if v and v > 0 and math.isfinite(v)]
    baseline = [v for v in baseline if v and v > 0 and math.isfinite(v)]
    if not values or not baseline:
        return float("nan"), float("nan"), float("nan")
    ratio = geometric_mean(values) / geometric_mean(baseline)
    draws = sorted(geometric_mean(rng.choices(values, k=len(values))) /
                   geometric_mean(rng.choices(baseline, k=len(baseline))) for _ in range(n))
    return ratio, draws[int(0.025 * n)], draws[int(0.975 * n) - 1]


STUCK_SPLIT_ESS = 2.0


def is_stuck(m):
    v = m.get("ess_split_min")
    return isinstance(v, (int, float)) and math.isfinite(v) and v < STUCK_SPLIT_ESS


def fewer_stuck_p(stuck, n, stuck_base, n_base):
    total, k = n + n_base, stuck + stuck_base
    if k == 0:
        return float("nan")
    return sum(math.comb(k, i) * math.comb(total - k, n - i) for i in range(0, stuck + 1)) / math.comb(total, n)


def is_default(run_id):
    return "_default_s" in run_id


def run_prefix(run_id):
    return "_".join(run_id.split("_")[:2])


def summarise(conditions, metrics):
    cells = {}
    settings = {}
    for row in conditions:
        if row["arm"] in ("reference", "pilot") or row["run_id"] not in metrics:
            continue
        cells.setdefault((row["arm"], row["dataset"], row["likelihood"], row["level"]), []).append(
            metrics[row["run_id"]])
        settings[(row["arm"], row["dataset"], row["likelihood"], row["level"])] = row
    defaults = {}
    for run_id, m in metrics.items():
        if is_default(run_id):
            defaults.setdefault(run_prefix(run_id), []).append(m)

    rows = []
    for (arm, dataset, likelihood, level), runs in sorted(cells.items(), key=lambda kv: (
            kv[0][0], kv[0][1], kv[0][2], _numeric(kv[0][3]))):
        base = defaults.get(run_prefix(runs[0]["run_id"]), [])
        cond = settings[(arm, dataset, likelihood, level)]
        row = {"arm": arm, "dataset": dataset, "likelihood": likelihood, "level": level,
               "sl_target": cond.get("sl_target", ""),
               "sl_weight_factor": cond.get("sl_weight_factor", ""),
               "stl2_target": cond.get("stl2_target", ""),
               "stl2_factor": cond.get("stl2_factor", ""),
               "n": len(runs), "n_default": len(base)}
        row["n_stuck"] = sum(map(is_stuck, runs))
        row["n_stuck_default"] = sum(map(is_stuck, base))
        row["fewer_stuck_p"] = fewer_stuck_p(row["n_stuck"], len(runs), row["n_stuck_default"], len(base))
        for metric in METRICS:
            for key in (f"eff_{metric}", f"mst_{metric}"):
                row[f"{key}_gmean"] = geometric_mean([r.get(key) for r in runs])
                ratio, lo, hi = bootstrap_ratio([r.get(key) for r in runs], [r.get(key) for r in base])
                row[f"{key}_ratio"], row[f"{key}_lo"], row[f"{key}_hi"] = ratio, lo, hi
        for key in ("asdsf", "asdsf_x_sqrt_hours", "sl_accept", "sl_tuning", "sl_topologyChanged",
                    "sl_us_per_op", "sl_time_fraction", "sl_stl2_accept", "sl_stl2_tuning",
                    "sl_stl2_us_per_op", "sl_stl2_time_fraction", "fhspr_accept", "fhspr_time_fraction", "hours"):
            vals = [r[key] for r in runs if isinstance(r.get(key), (int, float)) and math.isfinite(r[key])]
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

    def line_plot(arm_rows, x_field, x_label, filename_prefix):
        for dataset, likelihood in sorted({(r["dataset"], r["likelihood"]) for r in arm_rows}):
            cell = [r for r in arm_rows if r["dataset"] == dataset and r["likelihood"] == likelihood]
            xs_raw = sorted({r[x_field] for r in cell if r[x_field] != ""})
            try:
                xs = sorted(float(x) for x in xs_raw)
            except (ValueError, TypeError):
                xs = xs_raw
            fig, axes = plt.subplots(1, len(METRICS), figsize=(4.2 * len(METRICS), 3.5), squeeze=False)
            for ax, metric in zip(axes[0], METRICS):
                by_x = {}
                for r in cell:
                    try:
                        x = float(r[x_field]) if r[x_field] != "" else None
                    except (ValueError, TypeError):
                        x = None
                    if x is not None:
                        by_x.setdefault(x, []).append(r[f"eff_{metric}_ratio"])
                ys = [geometric_mean(by_x.get(x, [])) for x in xs]
                ax.axhline(1.0, color="gray", linewidth=0.8, linestyle="--")
                ax.plot(xs, ys, marker="o")
                ax.set_xlabel(x_label)
                ax.set_ylabel("ESS/h ratio vs control")
                ax.set_title(f"{metric}")
                if x_field == "sl_target":
                    ax.axvline(0.234, color="gray", linewidth=0.8, linestyle=":")
            fig.suptitle(f"{dataset} ({likelihood}): {filename_prefix.replace('_', ' ')}", y=1.02)
            fig.tight_layout()
            fig.savefig(ANALYSIS_DIR / f"{filename_prefix}_{dataset}-{likelihood}.png", dpi=120, bbox_inches="tight")
            plt.close(fig)

    target_rows = [r for r in rows if r["arm"] == "target_sweep" and r["level"] != "default"]
    if target_rows:
        line_plot(target_rows, "sl_target", "STL target acceptance", "target_sweep")

    weight_rows = [r for r in rows if r["arm"] == "weight_sweep" and r["level"] != "default"]
    if weight_rows:
        line_plot(weight_rows, "sl_weight_factor", "STL weight factor (× nTaxa)", "weight_sweep")

    two_stl_rows = [r for r in rows if r["arm"] == "two_stl" and r["stl2_target"] != ""]
    for dataset, likelihood in sorted({(r["dataset"], r["likelihood"]) for r in two_stl_rows}):
        cell = [r for r in two_stl_rows if r["dataset"] == dataset and r["likelihood"] == likelihood]
        targets = sorted({float(r["stl2_target"]) for r in cell})
        factors = sorted({float(r["stl2_factor"]) for r in cell})
        if not targets or not factors:
            continue
        fig, axes = plt.subplots(1, len(METRICS), figsize=(4.2 * len(METRICS), 0.6 * len(targets) + 2),
                                 squeeze=False)
        for k, (ax, metric) in enumerate(zip(axes[0], METRICS)):
            z = np.full((len(targets), len(factors)), np.nan)
            for r in cell:
                try:
                    ti = targets.index(float(r["stl2_target"]))
                    fi = factors.index(float(r["stl2_factor"]))
                    z[ti, fi] = r[f"eff_{metric}_ratio"]
                except (ValueError, IndexError):
                    pass
            image = ax.imshow(np.log2(z), cmap="RdBu", vmin=-1.5, vmax=1.5, aspect="auto")
            for i in range(len(targets)):
                for j in range(len(factors)):
                    if math.isfinite(z[i, j]):
                        ax.text(j, i, f"{z[i, j]:.2f}", ha="center", va="center", fontsize=9,
                                color="white" if abs(math.log2(z[i, j])) > 0.9 else "black")
            ax.set_xticks(range(len(factors)), [f"{f:g}" for f in factors])
            ax.set_yticks(range(len(targets)), [f"{t:g}" for t in targets] if k == 0 else [])
            ax.set_xlabel("stl2 weight factor (× nTaxa)")
            if k == 0:
                ax.set_ylabel("stl2 target acceptance")
            ax.set_title(f"{metric}: ESS/h vs control")
        fig.colorbar(image, ax=axes[0].tolist(), label="log2 ratio to control", shrink=0.8)
        fig.suptitle(f"{dataset} ({likelihood}): second SubtreeLeap", y=1.06)
        fig.savefig(ANALYSIS_DIR / f"two_stl_{dataset}-{likelihood}.png", dpi=120, bbox_inches="tight")
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
    REFERENCE_DIR = work_dir(config) / "reference"
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
