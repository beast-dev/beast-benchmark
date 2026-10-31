#!/usr/bin/env python3
"""
Writes conditions.csv for a phase of the chain-design study.

  python3 design.py throughput    # pilot P1: states/s for chains x threads (and P3, JVM/GC)
  python3 design.py tempering     # pilot P2: ladders (needs the P1 speeds)
  python3 design.py main          # the main grid (needs pilot/summary.json)

Phases are appended to conditions.csv; existing rows are kept.
"""

import random
import sys

import cd_common as cc


def row(**kw):
    r = {k: "" for k in cc.CONDITION_FIELDS}
    r.update({k: str(v) for k, v in kw.items()})
    return r


def throughput_rows(config):
    p = config["pilot"]
    rows = []
    for rep in range(1, p["repeats"] + 1):
        for chains, threads in p["throughput_grid"]:
            rows.append(row(run_id=f"p1_K{chains}_t{threads}_r{rep}", phase="throughput",
                            config=f"K{chains}xt{threads}", chains=chains, tempered=1, threads=threads,
                            seed=100 + rep, minutes=p["minutes"], log_every=p["log_every"]))
        # tempered sets at the main design's thread counts (speed of a set, barrier included)
        for cid, c in config["main"]["configs"].items():
            if c.get("tempered"):
                rows.append(row(run_id=f"p1_{cid}_r{rep}", phase="throughput", config=cid,
                                chains=c["chains"], tempered=c["tempered"], threads=c["threads"],
                                seed=100 + rep, minutes=p["minutes"], log_every=p["log_every"],
                                ladder_ratio=config["tempering"]["pilot_ladder_ratios"][1],
                                swap_every=1000))
    # P3: 16 chains in one JVM (two collectors) against 16 separate JVMs
    jvm = p["jvm_minutes"]
    rows.append(row(run_id="p3_onejvm_parallelgc", phase="jvm", config="K16xt1", chains=16, tempered=1,
                    threads=1, seed=300, minutes=jvm, log_every=p["log_every"], gc="-XX:+UseParallelGC"))
    rows.append(row(run_id="p3_onejvm_g1", phase="jvm", config="K16xt1", chains=16, tempered=1,
                    threads=1, seed=300, minutes=jvm, log_every=p["log_every"], gc="-XX:+UseG1GC"))
    rows.append(row(run_id="p3_16jvms", phase="jvm", config="K1xt1", chains=1, tempered=1,
                    threads=1, seed=300, minutes=jvm, log_every=p["log_every"], processes=16))
    return rows


def speed(summary, key):
    s = summary["speeds"].get(key)
    if s is None:
        sys.exit(f"pilot/summary.json has no speed for {key}: run the throughput pilot and analyse.py pilot")
    return s


def tempering_rows(config, summary):
    t = config["tempering"]
    minutes = config["pilot"]["tempering_minutes"]
    rows = []
    # a single tempering set at 4 threads per chain (T1): ladders compared at the same cost
    c = config["main"]["configs"]["T1"]
    per_chain = speed(summary, "T1")
    swap = max(100, int(round(per_chain * t["swap_block_seconds"], -2)))
    for ratio in t["pilot_ladder_ratios"]:
        rows.append(row(run_id=f"p2_T1_r{ratio}", phase="tempering", config="T1", chains=c["chains"],
                        tempered=c["tempered"], threads=c["threads"], seed=200, minutes=minutes,
                        states=int(per_chain * minutes * 60), log_every=config["pilot"]["log_every"],
                        ladder_ratio=ratio, swap_every=swap))
    return rows


def main_rows(config, summary):
    m = config["main"]
    ratio = config["tempering"]["ladder_ratio"] or summary.get("ladder_ratio")
    if not ratio:
        sys.exit('set "tempering": {"ladder_ratio": ...} in study.json (from the tempering pilot)')
    t_star = summary["t_star"]
    seconds = m["hours"] * 3600
    blocks = []
    for seed in m["seeds"]:
        block = []
        for cid, c in m["configs"].items():
            threads = t_star if c["threads"] == "t_star" else c["threads"]
            key = cid if c.get("tempered") else f"K{c['chains']}xt{threads}"
            per_chain = speed(summary, key)
            states = int(per_chain * seconds * 1.02)  # finishes at about the time limit
            log_every = max(1000, int(round(states / m["samples_per_chain"], -3)))
            r = row(run_id=f"{cid}_s{seed}", phase="main", config=cid, chains=c["chains"],
                    tempered=c.get("tempered", 1), threads=threads, seed=seed,
                    minutes=m["hours"] * 60, states=states, log_every=log_every)
            if c.get("tempered"):
                r["ladder_ratio"] = str(ratio)
                r["swap_every"] = str(max(100, int(round(per_chain * config["tempering"]["swap_block_seconds"], -2))))
            block.append(r)
        # configurations in random order within each seed block (drift in machine state
        # is spread over the configurations rather than confounded with them)
        random.Random(seed).shuffle(block)
        blocks += block
    return blocks


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ("throughput", "tempering", "main"):
        sys.exit(__doc__)
    config = cc.load_config()
    phase = sys.argv[1]
    if phase == "throughput":
        rows = throughput_rows(config)
    else:
        summary = cc.load_pilot_summary()
        if summary is None:
            sys.exit("no pilot/summary.json: run the throughput pilot, then python3 analyse.py pilot")
        rows = tempering_rows(config, summary) if phase == "tempering" else main_rows(config, summary)
    all_rows = cc.write_conditions(rows, append=True)
    added = [r for r in all_rows if r["phase"] == phase]
    hours = sum(float(r["minutes"]) for r in added) / 60  # each run has the machine to itself
    print(f"{phase}: {len(added)} runs, about {hours:.1f} hours of machine time -> {cc.paths()['conditions']}")


if __name__ == "__main__":
    main()
