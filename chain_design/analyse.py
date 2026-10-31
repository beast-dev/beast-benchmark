#!/usr/bin/env python3
"""
Analysis for the chain-design study.

  python3 analyse.py pilot     # pilot/summary.json: speeds, best threads, gains, ladder, JVM/GC, joint SD
  python3 analyse.py mixing    # the mixing pilot scored against its own pooled posterior
  python3 analyse.py main      # analysis/chains.csv, runs.csv, curves.csv, configs.csv, plots

Every chain whose samples count (an independent chain, or a tempering set's cold chain) is
condensed in two stages, both cached next to the run:
  A. its log traces, and the clade counts over its second half (for a pooled reference);
  B. its trees as per-sample traces against the reference: the RF distance to the reference's
     focal tree and indicators of its most uncertain clades.
The reference is either fixed (study.json "reference": a reference file, as for tall-1) or
"pooled": built from the second halves of all the chains of the main and mixing runs, with
each run scored against the pool without its own samples. Tree files are gzipped once read.
"""

import csv
import gzip
import hashlib
import json
import math
import multiprocessing
import os
import re
import sys
from pathlib import Path

import numpy as np

import cd_common as cc

WARMUP_SECONDS = 60  # JIT warm-up excluded from speed measurements (still charged to time)


# ---------------------------------------------------------------------------
# Multi-chain ESS and R-hat (Vehtari et al. 2021, as in ArviZ/Stan)
# ---------------------------------------------------------------------------

def _norm_ppf(p):
    """Inverse standard normal CDF (Acklam's rational approximation), vectorised."""
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00]
    p = np.asarray(p, dtype=float)
    out = np.empty_like(p)
    lo, hi = p < 0.02425, p > 1 - 0.02425
    mid = ~(lo | hi)
    q = np.sqrt(-2 * np.log(p[lo]))
    out[lo] = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
              ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    q = np.sqrt(-2 * np.log(1 - p[hi]))
    out[hi] = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    q = p[mid] - 0.5
    r = q * q
    out[mid] = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
               (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)
    return out


def _rank_normalise(x):
    flat = x.ravel()
    order = flat.argsort(kind="mergesort")
    ranks = np.empty(len(flat))
    ranks[order] = np.arange(1, len(flat) + 1)
    # average the ranks of ties
    sorted_vals = flat[order]
    start = 0
    for i in range(1, len(flat) + 1):
        if i == len(flat) or sorted_vals[i] != sorted_vals[start]:
            if i - start > 1:
                ranks[order[start:i]] = (start + 1 + i) / 2.0
            start = i
    return _norm_ppf((ranks - 0.375) / (len(flat) + 0.25)).reshape(x.shape)


def _split(x):
    n = x.shape[1] // 2
    return np.concatenate([x[:, :n], x[:, x.shape[1] - n:]], axis=0)


def multichain_ess(chains, rank=True):
    """Bulk ESS of an (m chains x n samples) array, split chains, Geyer's initial monotone
    sequence. Disagreement between chains lowers it."""
    x = np.asarray(chains, dtype=float)
    if x.ndim != 2 or x.shape[1] < 8:
        return float("nan")
    if np.all(x == x.flat[0]):
        return float("nan")
    x = _split(x)
    if rank:
        x = _rank_normalise(x)
    m, n = x.shape
    centred = x - x.mean(axis=1, keepdims=True)
    n_fft = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(centred, n_fft, axis=1)
    acov = np.fft.irfft(f * np.conjugate(f), n_fft, axis=1)[:, :n] / n
    mean_var = acov[:, 0].mean() * n / (n - 1)
    var_plus = mean_var * (n - 1) / n
    if m > 1:
        var_plus += x.mean(axis=1).var(ddof=1)
    if var_plus <= 0:
        return float("nan")
    rho = np.zeros(n)
    rho_even = 1.0
    rho[0] = 1.0
    rho_odd = 1.0 - (mean_var - acov[:, 1].mean()) / var_plus
    rho[1] = rho_odd
    t = 1
    while t < n - 3 and rho_even + rho_odd > 0:
        rho_even = 1.0 - (mean_var - acov[:, t + 1].mean()) / var_plus
        rho_odd = 1.0 - (mean_var - acov[:, t + 2].mean()) / var_plus
        if rho_even + rho_odd >= 0:
            rho[t + 1] = rho_even
            rho[t + 2] = rho_odd
        t += 2
    max_t = t - 2
    if rho_even > 0 and max_t + 1 < n:
        rho[max_t + 1] = rho_even
    t = 1
    while t <= max_t - 2:
        if rho[t + 1] + rho[t + 2] > rho[t - 1] + rho[t]:
            rho[t + 1] = (rho[t - 1] + rho[t]) / 2.0
            rho[t + 2] = rho[t + 1]
        t += 2
    total = m * n
    tau = -1.0 + 2.0 * rho[:max_t + 1].sum() + rho[max_t + 1:max_t + 2].sum()
    tau = max(tau, 1.0 / math.log10(total))
    return total / tau


def split_rhat(chains):
    """Split R-hat (not rank-normalised) of an (m x n) array; with m = 1 it compares the
    chain's two halves."""
    x = np.asarray(chains, dtype=float)
    if x.shape[1] < 8:
        return float("nan")
    x = _split(x)
    n = x.shape[1]
    w = x.var(axis=1, ddof=1).mean()
    b = n * x.mean(axis=1).var(ddof=1)
    if w <= 0:
        return float("nan")
    return math.sqrt(((n - 1) / n * w + b / n) / w)


# ---------------------------------------------------------------------------
# Condensing a chain
# ---------------------------------------------------------------------------

def _trees(trees_path: Path):
    """(state, clade set) of each complete tree, from the tree file or its gzipped copy."""
    plain = trees_path
    zipped = Path(str(trees_path) + ".gz")
    if plain.exists():
        _, trees = cc.read_trees(plain)
        return trees
    if zipped.exists():
        tmp = trees_path.with_suffix(f".trees.tmp{os.getpid()}")
        with gzip.open(zipped, "rt") as fin, open(tmp, "w") as fout:
            fout.write(fin.read())
        try:
            _, trees = cc.read_trees(tmp)
        finally:
            tmp.unlink()
        return trees
    return []


def _gzip_trees(trees_path: Path):
    if trees_path.exists():
        with open(trees_path, "rb") as fin, gzip.open(str(trees_path) + ".gz", "wb") as fout:
            fout.write(fin.read())
        trees_path.unlink()


def stage_a(log_path: Path, trees_path: Path, count_clades: bool):
    """Log traces of a chain (states that have a complete tree) and, if asked, the counts of
    each clade over the chain's second half."""
    out = log_path.with_suffix(".stageA.json.gz")
    if out.exists():
        with gzip.open(out, "rt") as f:
            chain = json.load(f)
        if chain.get("counted") or not count_clades:
            return chain
    columns, data = cc.read_log(log_path)
    if data.size == 0:
        return None
    col = {c: i for i, c in enumerate(columns)}
    trees = _trees(trees_path)
    by_state = {s: clades for s, clades in trees}
    states = data[:, col["state"]].astype(np.int64)
    keep = [i for i, s in enumerate(states) if int(s) in by_state]  # a stopped run's last tree may be cut short
    data = data[keep]
    states = states[keep]
    chain = dict(
        state=[int(s) for s in states],
        seconds=[float(v) for v in data[:, col["secondsElapsed"]]],
        params={p: [float(v) for v in data[:, col[p]]] for p in cc.CONTINUOUS_PARAMETERS if p in col},
        counted=count_clades,
    )
    if count_clades:
        half = len(states) // 2
        counts = {}
        for s in states[half:]:
            for c in by_state[int(s)]:
                counts[c] = counts.get(c, 0) + 1
        chain["second_half"] = len(states) - half
        chain["clade_counts"] = {format(c, "x"): n for c, n in counts.items()}
    with gzip.open(out, "wt") as f:
        json.dump(chain, f)
    return chain


def stage_b(log_path: Path, trees_path: Path, reference, ref_key: str):
    """RF distances to the reference's focal tree and indicators of its uncertain clades, for
    each sample kept in stage A."""
    out = log_path.with_suffix(".stageB.json.gz")
    if out.exists():
        with gzip.open(out, "rt") as f:
            traces = json.load(f)
        if traces.get("ref_key") == ref_key:
            return traces
    a = stage_a(log_path, trees_path, False)
    by_state = {s: clades for s, clades in _trees(trees_path)}
    focal = reference["focal"]
    uncertain = reference["uncertain"]
    sets = [by_state[s] for s in a["state"]]
    traces = dict(ref_key=ref_key,
                  rf=[len(c ^ focal) for c in sets],
                  clades="".join("".join("1" if u in c else "0" for u in uncertain) for c in sets))
    with gzip.open(out, "wt") as f:
        json.dump(traces, f)
    return traces


def _stage_a_job(args):
    log_path, trees_path, count = args
    return str(log_path), stage_a(Path(log_path), Path(trees_path), count)


def _stage_b_job(args):
    log_path, trees_path, reference, ref_key = args
    traces = stage_b(Path(log_path), Path(trees_path), reference, ref_key)
    _gzip_trees(Path(trees_path))
    return str(log_path), traces


def parallel_map(job, items):
    items = list(items)
    if not items:
        return []
    workers = max(1, min(len(items), (os.cpu_count() or 2) // 2))
    if workers == 1:
        return [job(i) for i in items]
    with multiprocessing.get_context("spawn").Pool(workers) as pool:
        return pool.map(job, items)


def clade_matrix(traces, n_clades):
    a = np.frombuffer(traces["clades"].encode(), dtype=np.uint8) - ord("0")
    return a.reshape(-1, n_clades).astype(float)


# ---------------------------------------------------------------------------
# References
# ---------------------------------------------------------------------------

def reference_key(reference):
    h = hashlib.sha256()
    h.update(",".join(format(c, "x") for c in reference["uncertain"]).encode())
    h.update(",".join(sorted(format(c, "x") for c in reference["focal"])).encode())
    return h.hexdigest()[:16]


def fixed_reference(config):
    ref = cc.load_reference(config)
    stats = config["reference_stats"]
    ref.update(joint_q01=stats["joint_q01"], mean=stats["mean"], sd=stats["sd"])
    return ref


def pooled_reference(config, runs_chains):
    """runs_chains: {run_id: [stage A chain, ...]}. Returns the reference built from every run
    plus, per run, what to subtract to leave the run out."""
    total = {}
    n_total = 0
    per_run = {}
    param_sums, param_sq, param_n = {}, {}, 0
    joints = []
    for run_id, chains in runs_chains.items():
        counts = {}
        n = 0
        sums = {}
        m = 0
        for ch in chains:
            for k, v in ch["clade_counts"].items():
                c = int(k, 16)
                counts[c] = counts.get(c, 0) + v
            n += ch["second_half"]
            half = len(ch["state"]) // 2
            for p, values in ch["params"].items():
                tail = values[half:]
                sums[p] = sums.get(p, 0.0) + sum(tail)
                param_sq[p] = param_sq.get(p, 0.0) + sum(v * v for v in tail)
                if p == "joint":
                    joints += tail
            m += len(ch["state"]) - half
        per_run[run_id] = dict(counts=counts, n=n, sums=sums, m=m)
        for c, v in counts.items():
            total[c] = total.get(c, 0) + v
        n_total += n
        for p, v in sums.items():
            param_sums[p] = param_sums.get(p, 0.0) + v
        param_n += m
    if n_total == 0:
        sys.exit("no samples for the pooled reference")
    freqs = {c: v / n_total for c, v in total.items()}
    uncertain = sorted((c for c, f in freqs.items() if 0.05 <= f <= 0.95),
                       key=lambda c: abs(freqs[c] - 0.5))[:config["focal_clades"]]
    focal = {c for c, f in freqs.items() if f > 0.5}  # majority-rule tree
    mean = {p: v / param_n for p, v in param_sums.items()}
    sd = {p: math.sqrt(max(0.0, param_sq[p] / param_n - mean[p] ** 2)) for p in mean}
    joints.sort()
    ref = dict(frequencies=freqs, uncertain=uncertain, uncertain_freq=[freqs[c] for c in uncertain],
               focal=focal, joint_q01=joints[int(0.01 * (len(joints) - 1))], mean=mean, sd=sd,
               pool=dict(total=total, n_total=n_total, per_run=per_run, param_sums=param_sums,
                         param_n=param_n))
    return ref


def run_reference(reference, run_id):
    """The reference for scoring a run: a pooled reference without the run's own samples."""
    pool = reference.get("pool")
    if not pool or run_id not in pool["per_run"]:
        return reference
    own = pool["per_run"][run_id]
    n = pool["n_total"] - own["n"]
    m = pool["param_n"] - own["m"]
    if n <= 0 or m <= 0:
        return reference
    out = dict(reference)
    out["uncertain_freq"] = [(pool["total"].get(c, 0) - own["counts"].get(c, 0)) / n
                             for c in reference["uncertain"]]
    out["mean"] = {p: (pool["param_sums"][p] - own["sums"].get(p, 0.0)) / m for p in pool["param_sums"]}
    return out


def pool_check(reference, conds, groups):
    """ASDSF between pools built from two groups of seeds (clades of at least 1% in either)."""
    per_run = reference["pool"]["per_run"]
    seed_of = {c["run_id"]: int(c["seed"]) for c in conds}
    freqs = []
    for group in groups:
        counts, n = {}, 0
        for run_id, own in per_run.items():
            if seed_of.get(run_id) in group:
                for c, v in own["counts"].items():
                    counts[c] = counts.get(c, 0) + v
                n += own["n"]
        if n == 0:
            return float("nan"), float("nan")
        freqs.append({c: v / n for c, v in counts.items()})
    keys = {c for f in freqs for c, v in f.items() if v >= 0.01}
    devs = [abs(freqs[0].get(c, 0.0) - freqs[1].get(c, 0.0)) for c in keys]
    return float(np.mean(devs)) if devs else 0.0, float(np.max(devs)) if devs else 0.0


# ---------------------------------------------------------------------------
# Pilot
# ---------------------------------------------------------------------------

def chain_speed(chain):
    sec = np.array(chain["seconds"])
    st = np.array(chain["state"], dtype=float)
    after = sec >= WARMUP_SECONDS
    if after.sum() < 2:
        after = np.ones_like(sec, dtype=bool)
    s, t = st[after], sec[after]
    return float((s[-1] - s[0]) / (t[-1] - t[0])) if t[-1] > t[0] else float("nan")


def log_speed(log_path):
    columns, data = cc.read_log(log_path)
    if data.size == 0:
        return float("nan")
    col = {c: i for i, c in enumerate(columns)}
    return chain_speed(dict(seconds=list(data[:, col["secondsElapsed"]]), state=list(data[:, col["state"]])))


def gc_pause_fraction(run_dir: Path, wall):
    path = run_dir / "gc.log"
    if not path.exists() or not wall:
        return float("nan")
    total_ms = 0.0
    for line in path.read_text().splitlines():
        m = re.search(r"Pause.*?([\d.]+)ms\s*$", line)
        if m:
            total_ms += float(m.group(1))
    return total_ms / 1000.0 / wall


def gpu_utilisation(run_dir: Path):
    """Mean utilisation (%) of each GPU sampled during the run."""
    path = run_dir / "gpu.csv"
    if not path.exists():
        return {}
    util = {}
    for line in path.read_text().splitlines():
        parts = [p.strip() for p in line.split(",")]
        if len(parts) >= 3:
            try:
                util.setdefault(parts[1], []).append(float(parts[2]))
            except ValueError:
                continue
    return {k: float(np.mean(v)) for k, v in util.items()}


SWAP_ROW = re.compile(r"^\s*(\d+)\s*<->\s*(\d+)\s+\S+\s*<->\s*\S+\s+(\d+)\s+(\d+)\s+([\d.]+)")


def swap_table(out_path: Path):
    rows = []
    if out_path.exists():
        for line in out_path.read_text().splitlines():
            m = SWAP_ROW.match(line)
            if m:
                i, j, n, acc, rate = m.groups()
                rows.append(dict(i=int(i), j=int(j), attempts=int(n), accepted=int(acc), rate=float(rate)))
    return rows


def analyse_pilot(config):
    conds = [c for c in cc.read_conditions() if c["phase"] in ("throughput", "tempering", "jvm", "mixing")]
    speeds, totals = {}, {}
    jvm, tempering, gpu = {}, {}, {}
    mixing_joint = []
    for cond in conds:
        run_dir = cc.paths()["runs"] / cond["run_id"]
        result_path = run_dir / "result.json"
        if not result_path.exists():
            continue
        result = json.load(open(result_path))
        logs = [lg for lg, _ in cc.chain_files(run_dir, cond["run_id"], cond) if lg.exists()]
        if cond["phase"] == "mixing":
            for lg in logs:
                columns, data = cc.read_log(lg)
                if data.size:
                    joint = data[:, columns.index("joint")]
                    mixing_joint += list(joint[len(joint) // 2:])
            continue
        chain_speeds = [s for s in (log_speed(lg) for lg in logs) if not math.isnan(s)]
        if not chain_speeds:
            continue
        # a mixed configuration's chain length is set by its fastest (GPU) chains
        per_chain = float(max(chain_speeds) if cond.get("device") == "mixed" else np.mean(chain_speeds))
        total = float(np.sum(chain_speeds))
        if cond["phase"] == "throughput":
            speeds.setdefault(cond["config"], []).append(per_chain)
            totals.setdefault(cond["config"], []).append(total)
            if (cond.get("device") or "cpu") != "cpu":
                gpu.setdefault(cond["config"], []).append(gpu_utilisation(run_dir))
        elif cond["phase"] == "jvm":
            jvm[cond["run_id"]] = dict(per_chain=float(np.mean(chain_speeds)), total=total,
                                       gc_pause_fraction=gc_pause_fraction(run_dir, result["wall_seconds"]))
        else:
            table = swap_table(run_dir / "beast.out")
            neighbours = [r["rate"] for r in table if abs(r["i"] - r["j"]) == 1]
            tempering.setdefault(cond["ladder_ratio"], []).append(dict(
                config=cond["config"], per_chain=per_chain, swaps=table,
                min_neighbour_rate=min(neighbours) if neighbours else None))

    summary = dict(speeds={}, speed_spread={}, totals={}, jvm=jvm, tempering=tempering, gpu_utilisation=gpu)
    for key, values in speeds.items():
        summary["speeds"][key] = float(np.mean(values))
        summary["speed_spread"][key] = float((max(values) - min(values)) / np.mean(values)) if len(values) > 1 else None
        summary["totals"][key] = float(np.mean(totals[key]))
    singles = {}
    for key, v in summary["speeds"].items():
        m = re.fullmatch(r"K1xt(\d+)", key)
        if m:
            singles[int(m.group(1))] = v
    if singles:
        best = max(singles.values())
        # the fewest threads within 3% of the best single-chain speed
        summary["t_star"] = min(t for t, v in singles.items() if v >= 0.97 * best)
    summary["include"] = {}
    for name, (a, b, factor) in config["pilot"].get("gains", {}).items():
        if a in summary["totals"] and b in summary["totals"]:
            ratio = summary["totals"][a] / summary["totals"][b]
            summary["include"][name] = bool(ratio >= factor)
            summary.setdefault("gain_ratios", {})[name] = ratio
    if mixing_joint:
        summary["joint_sd"] = float(np.std(mixing_joint))
    # a ladder must work for every tempered configuration tried with it
    rated = {r: min(v["min_neighbour_rate"] for v in vs) for r, vs in tempering.items()
             if all(v["min_neighbour_rate"] is not None for v in vs)}
    usable = [r for r, rate in rated.items() if rate >= 0.2]
    if usable:
        # the hottest ladder whose neighbouring chains still swap at least 20% of the time
        summary["ladder_ratio"] = float(min(usable, key=float))
    elif rated:
        best = max(rated, key=rated.get)
        summary["ladder_ratio"] = float(best)
        print(f"WARNING: no ladder reached a 20% swap rate between every pair of neighbours; "
              f"using {best}, the best ({rated[best]:.3f}). Consider a gentler ratio.")
    cc.paths()["pilot"].mkdir(parents=True, exist_ok=True)
    with open(cc.paths()["pilot_summary"], "w") as f:
        json.dump(summary, f, indent=1)

    print("\nStates per second (per chain; total over chains; spread = (max - min) / mean over repeats)")
    for key in sorted(summary["speeds"]):
        spread = summary["speed_spread"][key]
        util = gpu.get(key)
        util_text = f"  GPU % {util[-1]}" if util and util[-1] else ""
        print(f"  {key:14} {summary['speeds'][key]:9.0f} {summary['totals'][key]:10.0f}  "
              f"spread {spread if spread is None else round(spread, 3)}{util_text}")
    if "t_star" in summary:
        print(f"best threads for one chain (t*): {summary['t_star']}")
    for name, ratio in summary.get("gain_ratios", {}).items():
        print(f"{name}: x{ratio:.2f} -> {'include' if summary['include'][name] else 'leave out'}")
    if "joint_sd" in summary:
        print(f"joint SD (mixing pilot, second halves): {summary['joint_sd']:.1f}")
    for r, vs in sorted(tempering.items()):
        print(f"ladder ratio {r}: " + ", ".join(f"{v['config']} min neighbour swap rate {v['min_neighbour_rate']}" for v in vs))
    if "ladder_ratio" in summary:
        print(f"chosen ladder ratio: {summary['ladder_ratio']}")
    for k, v in jvm.items():
        print(f"{k}: total {v['total']:.0f} states/s, GC pauses {v['gc_pause_fraction']:.2%}")


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def chain_metrics(chain, traces, startup, config, ref):
    joint = np.array(chain["params"]["joint"])
    reached = np.nonzero(joint >= ref["joint_q01"])[0]
    burnin = int(reached[0]) if len(reached) else len(joint)
    sec = np.array(chain["seconds"])
    n = len(joint)
    clades = clade_matrix(traces, len(ref["uncertain"]))
    ref_f = np.array(ref["uncertain_freq"])

    def max_dev(rows):
        return float(np.max(np.abs(clades[rows].mean(axis=0) - ref_f))) if len(rows) else float("nan")

    post = np.arange(burnin, n)
    halves = (post[:len(post) // 2], post[len(post) // 2:])
    m = dict(samples=n, burnin_index=burnin,
             burnin_states=chain["state"][burnin] if burnin < n else float("nan"),
             burnin_minutes=(startup + sec[burnin]) / 60 if burnin < n else float("nan"),
             speed=chain_speed(chain),
             max_dev_first_half=max_dev(halves[0]), max_dev_second_half=max_dev(halves[1]))
    dev = config["stuck"]["clade_deviation"]
    m["stuck_at_end"] = bool(m["max_dev_second_half"] > dev)
    m["escaped"] = bool(m["max_dev_first_half"] > dev and not m["stuck_at_end"])
    if len(post) >= 8:
        step = chain["state"][1] - chain["state"][0]
        hours = (sec[-1] - sec[burnin]) / 3600
        mstates = (chain["state"][-1] - chain["state"][burnin]) / 1e6
        cont = min(cc.ess(np.array(chain["params"][p])[post], step) for p in chain["params"])
        clade = min((v if not math.isnan(v) else 1.0) for v in
                    (cc.ess(clades[post, j], step) for j in range(clades.shape[1])))
        rf = cc.ess(np.array(traces["rf"])[post], step)
        for name, v in (("cont", cont), ("clade", clade), ("rf", rf)):
            m[f"ess_{name}"] = v
            m[f"ess_{name}_per_mstate"] = v / mstates if mstates else float("nan")
            m[f"ess_{name}_per_hour"] = v / hours if hours else float("nan")
    return m


def curve_point(chains, traces, burnins, startups, t, config, ref):
    """Metrics of a run using the samples available t seconds after the process started."""
    targets = config["targets"]
    n_clades = len(ref["uncertain"])
    usable = []
    for chain, tr, b, s in zip(chains, traces, burnins, startups):
        sec = np.array(chain["seconds"]) + s
        idx = np.nonzero(sec <= t)[0]
        idx = idx[idx >= b]
        if len(idx) >= 8:
            usable.append((chain, tr, idx))
    point = dict(t_minutes=t / 60, chains_used=len(usable), cont_ess=0.0, clade_ess=0.0,
                 summed_cont_ess=0.0, summed_clade_ess=0.0, accurate=False, success=False)
    if not usable:
        return point
    length = min(len(idx) for _, _, idx in usable)
    params = list(usable[0][0]["params"])
    cont, summed_cont = [], []
    for p in params:
        arr = np.array([np.array(c["params"][p])[idx[-length:]] for c, _, idx in usable])
        cont.append(multichain_ess(arr, rank=True))
        summed_cont.append(sum(cc.ess(np.array(c["params"][p])[idx], c["state"][1] - c["state"][0])
                               for c, _, idx in usable))
    mats = [(clade_matrix(tr, n_clades), idx) for _, tr, idx in usable]
    clade, summed_clade = [], []
    for j in range(n_clades):
        arr = np.array([mat[idx[-length:], j] for mat, idx in mats])
        v = multichain_ess(arr, rank=False)
        clade.append(1.0 if math.isnan(v) else v)
        s = 0.0
        for (c, _, idx), (mat, _) in zip(usable, mats):
            e = cc.ess(mat[idx, j], c["state"][1] - c["state"][0])
            s += 1.0 if math.isnan(e) else e
        summed_clade.append(s)
    pooled = np.concatenate([mat[idx] for mat, idx in mats])
    clade_ok = np.max(np.abs(pooled.mean(axis=0) - np.array(ref["uncertain_freq"]))) <= targets["clade_tolerance"]
    means_ok = all(abs(np.mean(np.concatenate([np.array(c["params"][p])[idx] for c, _, idx in usable]))
                       - ref["mean"][p]) <= targets["mean_tolerance_sd"] * ref["sd"][p]
                   for p in params if p in ref["mean"] and ref["sd"].get(p))
    point.update(cont_ess=float(np.nanmin(cont)), clade_ess=float(np.min(clade)),
                 summed_cont_ess=float(np.nanmin(summed_cont)), summed_clade_ess=float(np.min(summed_clade)),
                 accurate=bool(clade_ok and means_ok))
    point["success"] = bool(point["accurate"] and point["cont_ess"] >= targets["continuous_ess"]
                            and point["clade_ess"] >= targets["clade_ess"])
    return point


def practical_burnin(chains, burnin_max=0.5, step=0.02, threshold=1.01):
    """Smallest common burn-in fraction after which the split R-hat of the joint density stays
    below the threshold for every larger fraction tried (2% grid)."""
    joints = [np.array(c["params"]["joint"]) for c in chains]
    n = min(len(j) for j in joints)
    fractions = np.arange(0, burnin_max + 1e-9, step)
    ok = []
    for f in fractions:
        start = int(math.ceil(n * f))
        arr = np.array([j[start:n] for j in joints])
        r = split_rhat(arr)
        ok.append(not math.isnan(r) and r < threshold)
    for i in range(len(fractions)):
        if all(ok[i:]):
            return float(fractions[i])
    return float("nan")


def fit_rate(curve, key):
    """Least-squares fit of ESS(t) = r (t - b) over the points after ESS first becomes positive."""
    pts = [(p["t_minutes"] / 60, p[key]) for p in curve if p[key] > 0]
    if len(pts) < 3:
        return float("nan"), float("nan")
    t, e = np.array(pts).T
    r, c = np.polyfit(t, e, 1)
    if r <= 0:
        return float("nan"), float("nan")  # no growth: the rate is not estimable
    return float(r), float(-c / r)


def predict(b_cont, r_cont, b_clade, r_clade, targets, limit_minutes):
    """Minutes to reach both ESS targets from the fitted rates; inf if not estimable or more than
    20 times the run length (an extrapolation that far means only "far beyond the run")."""
    if not (r_cont > 0 and r_clade > 0):
        return float("inf")
    minutes = max(b_cont + targets["continuous_ess"] / r_cont, b_clade + targets["clade_ess"] / r_clade) * 60
    return minutes if minutes <= 20 * limit_minutes else float("inf")


def km_median(times, censored):
    """Kaplan-Meier median of times with right-censoring (nan if not reached)."""
    order = sorted(zip(times, censored))
    at_risk = len(order)
    surv = 1.0
    for t, cen in order:
        if not cen:
            surv *= 1 - 1 / at_risk
            if surv <= 0.5:
                return t
        at_risk -= 1
    return float("nan")


def run_chain_files(cond):
    run_dir = cc.paths()["runs"] / cond["run_id"]
    return [(lg, tr) for lg, tr in cc.chain_files(run_dir, cond["run_id"], cond) if lg.exists()]


def finished(cond):
    path = cc.paths()["runs"] / cond["run_id"] / "result.json"
    return path.exists() and not json.load(open(path)).get("failed")


def analyse_main(config, score_phases):
    pooled = config["reference"] == "pooled"
    all_conds = cc.read_conditions()
    pool_conds = [c for c in all_conds if c["phase"] in (("main", "mixing") if "main" in score_phases else score_phases)
                  and finished(c)]
    conds = [c for c in all_conds if c["phase"] in score_phases and finished(c)]
    targets = config["targets"]
    out_dir = cc.paths()["analysis"] / ("_".join(score_phases))
    out_dir.mkdir(parents=True, exist_ok=True)

    # stage A for every chain that counts (clade counts only for a pooled reference)
    jobs = [(str(lg), str(tr), pooled) for c in (pool_conds if pooled else conds) for lg, tr in run_chain_files(c)]
    stage = dict(parallel_map(_stage_a_job, jobs))
    if pooled:
        runs_chains = {c["run_id"]: [stage[str(lg)] for lg, _ in run_chain_files(c) if stage.get(str(lg))]
                       for c in pool_conds}
        reference = pooled_reference(config, runs_chains)
        groups = config.get("reference_check_seeds")
        check = pool_check(reference, pool_conds, groups) if groups else (float("nan"), float("nan"))
        print(f"pooled reference: {len(pool_conds)} runs, {reference['pool']['n_total']} trees; "
              f"seeds {groups} ASDSF {check[0]:.4f} (max {check[1]:.3f})"
              + ("" if check[0] < 0.01 else "  WARNING: the pooled reference is not converged"))
        with open(out_dir / "reference.json", "w") as f:
            json.dump(dict(runs=len(pool_conds), trees=reference["pool"]["n_total"], asdsf_check=check[0],
                           max_dev_check=check[1], joint_q01=reference["joint_q01"], mean=reference["mean"],
                           sd=reference["sd"], uncertain_freq=reference["uncertain_freq"]), f, indent=1)
    else:
        reference = fixed_reference(config)
    ref_key = reference_key(reference)
    light = {k: v for k, v in reference.items() if k not in ("frequencies", "pool")}
    jobs = [(str(lg), str(tr), light, ref_key) for c in conds for lg, tr in run_chain_files(c)]
    traces_of = dict(parallel_map(_stage_b_job, jobs))

    chain_rows, run_rows, curve_rows = [], [], []
    step = config["curve_step_minutes"] * 60
    for cond in conds:
        run_id = cond["run_id"]
        run_dir = cc.paths()["runs"] / run_id
        result = json.load(open(run_dir / "result.json"))
        ref = run_reference(reference, run_id)
        chains, traces, startups, metrics = [], [], [], []
        for k, (lg, tr) in enumerate(run_chain_files(cond), start=1):
            chain = stage.get(str(lg)) or stage_a(lg, tr, False)
            if chain is None or str(lg) not in traces_of:
                continue
            startup = result["first_log_seconds"].get(str(lg), 0.0)
            m = chain_metrics(chain, traces_of[str(lg)], startup, config, ref)
            m.update(run_id=run_id, config=cond["config"], seed=cond["seed"], chain=k,
                     threads=cond["threads"], startup_minutes=startup / 60)
            chain_rows.append(m)
            chains.append(chain)
            traces.append(traces_of[str(lg)])
            startups.append(startup)
            metrics.append(m)
        if not chains:
            continue
        burnins = [m["burnin_index"] for m in metrics]
        limit = float(cond["minutes"]) * 60
        curve = []
        t = step
        while t <= limit + 1e-6:
            p = curve_point(chains, traces, burnins, startups, t, config, ref)
            p.update(run_id=run_id, config=cond["config"], seed=cond["seed"])
            curve.append(p)
            t += step
        curve_rows += curve
        reached = [p["t_minutes"] for p in curve if p["success"]]
        r_cont, b_cont = fit_rate(curve, "cont_ess")
        r_clade, b_clade = fit_rate(curve, "clade_ess")
        final = curve[-1]
        util = gpu_utilisation(run_dir)
        run = dict(run_id=run_id, config=cond["config"], seed=cond["seed"], chains=len(chains),
                   threads=cond["threads"], device=cond.get("device") or "cpu",
                   wall_minutes=result["wall_seconds"] / 60,
                   time_to_target=reached[0] if reached else float("nan"), censored=not reached,
                   final_cont_ess=final["cont_ess"], final_clade_ess=final["clade_ess"],
                   final_summed_cont_ess=final["summed_cont_ess"], final_summed_clade_ess=final["summed_clade_ess"],
                   final_accurate=final["accurate"],
                   cont_ess_per_hour=r_cont, cont_burnin_hours=b_cont,
                   clade_ess_per_hour=r_clade, clade_burnin_hours=b_clade,
                   predicted_time_to_target=predict(b_cont, r_cont, b_clade, r_clade, targets,
                                                    float(cond["minutes"])),
                   stuck_chains=sum(m["stuck_at_end"] for m in metrics),
                   escaped_chains=sum(m["escaped"] for m in metrics),
                   median_burnin_minutes=float(np.nanmedian([m["burnin_minutes"] for m in metrics])),
                   max_burnin_minutes=float(np.nanmax([m["burnin_minutes"] for m in metrics])),
                   practical_burnin_fraction=practical_burnin(chains),
                   total_states_per_second=float(np.nansum([m["speed"] for m in metrics])),
                   mean_gpu_utilisation=float(np.mean(list(util.values()))) if util else float("nan"),
                   gc_pause_fraction=gc_pause_fraction(run_dir, result["wall_seconds"]),
                   load_before=result["load_before"][0], load_after=result["load_after"][0])
        run_rows.append(run)
        print(f"{run_id}: {len(chains)} chains, target {'reached at ' + str(run['time_to_target']) + ' min' if reached else 'not reached'}, "
              f"final cont ESS {final['cont_ess']:.0f}, clade ESS {final['clade_ess']:.0f}, stuck {run['stuck_chains']}")

    write_csv(out_dir / "chains.csv", chain_rows)
    write_csv(out_dir / "runs.csv", run_rows)
    write_csv(out_dir / "curves.csv", curve_rows)

    p_stuck = np.mean([m["stuck_at_end"] for m in chain_rows]) if chain_rows else float("nan")
    config_rows = []
    config_ids = list(dict.fromkeys(c["config"] for c in conds))
    for cid in config_ids:
        runs = [r for r in run_rows if r["config"] == cid]
        if not runs:
            continue
        chains_here = [m for m in chain_rows if m["config"] == cid]
        k = runs[0]["chains"]
        limit = float(next(c for c in conds if c["config"] == cid)["minutes"])
        spec = config["main"]["configs"].get(cid, {})
        gpu_chains = k if spec.get("device") == "gpu" else int(spec.get("gpu_chains") or 0)
        cpu_threads = (k - gpu_chains) * int(runs[0]["threads"]) * (spec.get("tempered") or 1)

        def cv(values):
            v = [x for x in values if not math.isnan(x)]
            return float(np.std(v) / np.mean(v)) if len(v) > 1 and np.mean(v) else float("nan")

        median_burnin = float(np.nanmedian([m["burnin_minutes"] for m in chains_here])) \
            if any(not math.isnan(m["burnin_minutes"]) for m in chains_here) else float("nan")
        config_rows.append(dict(
            config=cid, chains=k, threads=runs[0]["threads"], device=runs[0]["device"], runs=len(runs),
            reached=sum(not r["censored"] for r in runs),
            km_median_time_to_target=km_median([r["time_to_target"] if not r["censored"] else limit for r in runs],
                                               [r["censored"] for r in runs]),
            median_predicted_time_to_target=float(np.nanmedian([r["predicted_time_to_target"] for r in runs])),
            median_cont_ess_per_hour=float(np.nanmedian([r["cont_ess_per_hour"] for r in runs])),
            cv_cont_ess_per_hour=cv([r["cont_ess_per_hour"] for r in runs]),
            median_clade_ess_per_hour=float(np.nanmedian([r["clade_ess_per_hour"] for r in runs])),
            cv_clade_ess_per_hour=cv([r["clade_ess_per_hour"] for r in runs]),
            median_final_cont_ess=float(np.median([r["final_cont_ess"] for r in runs])),
            median_final_clade_ess=float(np.median([r["final_clade_ess"] for r in runs])),
            median_burnin_minutes=median_burnin,
            cv_burnin_minutes=cv([m["burnin_minutes"] for m in chains_here]),
            burnin_cpu_core_minutes=median_burnin * cpu_threads,
            burnin_gpu_chain_minutes=median_burnin * gpu_chains * (spec.get("tempered") or 1),
            median_speed_per_chain=float(np.nanmedian([m["speed"] for m in chains_here])),
            median_total_states_per_second=float(np.nanmedian([r["total_states_per_second"] for r in runs])),
            mean_gpu_utilisation=float(np.nanmean([r["mean_gpu_utilisation"] for r in runs]))
            if any(not math.isnan(r["mean_gpu_utilisation"]) for r in runs) else float("nan"),
            runs_with_stuck_chain=sum(r["stuck_chains"] > 0 for r in runs),
            runs_accurate_at_end=sum(r["final_accurate"] for r in runs),
            diagnostic_power=1 - (1 - p_stuck) ** k - p_stuck ** k if k > 1 else 0.0,
            median_practical_burnin_fraction=float(np.nanmedian([r["practical_burnin_fraction"] for r in runs])),
            max_gc_pause_fraction=float(np.nanmax([r["gc_pause_fraction"] for r in runs])),
        ))
    write_csv(out_dir / "configs.csv", config_rows)
    print(f"\nchain stuck rate over all chains: {p_stuck:.3f}")
    for r in config_rows:
        print(f"{r['config']:4} reached {r['reached']}/{r['runs']}  KM median {r['km_median_time_to_target']:.0f} min  "
              f"predicted {r['median_predicted_time_to_target']:.0f} min  cont ESS/h {r['median_cont_ess_per_hour']:.0f}  "
              f"clade ESS/h {r['median_clade_ess_per_hour']:.0f}  burn-in {r['median_burnin_minutes']:.1f} min  "
              f"{r['median_total_states_per_second']:.0f} states/s")
    plot(curve_rows, config_rows, chain_rows, out_dir)


def write_csv(path, rows):
    if not rows:
        return
    keys = list(rows[0])
    for r in rows[1:]:
        keys += [k for k in r if k not in keys]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def plot(curves, configs, chains, out_dir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available: no plots")
        return
    ids = [c["config"] for c in configs]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, key, target in ((axes[0], "cont_ess", 200), (axes[1], "clade_ess", 100)):
        for i, cid in enumerate(ids):
            colour = plt.cm.tab20(i % 20)
            runs = sorted({p["run_id"] for p in curves if p["config"] == cid})
            for j, run in enumerate(runs):
                pts = [p for p in curves if p["run_id"] == run]
                ax.plot([p["t_minutes"] / 60 for p in pts], [p[key] for p in pts], color=colour, alpha=0.7,
                        label=cid if j == 0 else None)
        ax.axhline(target, color="grey", ls="--", lw=1)
        ax.set_xlabel("wall-clock hours")
        ax.set_ylabel("multi-chain ESS (" + ("continuous, min" if key == "cont_ess" else "worst of 50 clades") + ")")
        ax.set_yscale("symlog")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "ess_curves.png", dpi=130)

    fig, ax = plt.subplots(figsize=(8, 4))
    data = [[m["burnin_minutes"] for m in chains if m["config"] == cid and not math.isnan(m["burnin_minutes"])]
            for cid in ids]
    try:
        ax.boxplot(data, tick_labels=ids)
    except TypeError:  # matplotlib < 3.9
        ax.boxplot(data, labels=ids)
    ax.set_ylabel("burn-in per chain (minutes from start)")
    fig.tight_layout()
    fig.savefig(out_dir / "burnin.png", dpi=130)


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ("pilot", "mixing", "main"):
        sys.exit(__doc__)
    config = cc.load_config()
    if sys.argv[1] == "pilot":
        analyse_pilot(config)
    else:
        analyse_main(config, [sys.argv[1]])


if __name__ == "__main__":
    main()
