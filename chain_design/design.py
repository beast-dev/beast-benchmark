#!/usr/bin/env python3
"""
Writes conditions.csv for a phase of the chain-design study.

  python3 design.py throughput    # pilot: states/s for each configuration (and JVM/GC on tall-1)
  python3 design.py mixing        # pilot (darwin): one long run to measure mixing and the joint's SD
  python3 design.py tempering     # pilot: ladders (needs the throughput speeds; on darwin the joint SD)
  python3 design.py main          # the main grid (needs pilot/summary.json)

Phases are appended to conditions.csv; existing rows are kept.
"""

import random
import sys

import cd_common as cc

DEVICE_FIELDS = ("device", "gpus", "gpu_chains", "processes", "numa")


def row(**kw):
    r = {k: "" for k in cc.CONDITION_FIELDS}
    r.update({k: str(v) for k, v in kw.items() if v is not None})
    return r


def throughput_entries(config):
    """The pilot's throughput configurations as dicts with an id (tall-1 lists [chains, threads]
    pairs; darwin names its entries)."""
    p = config["pilot"]
    entries = []
    for e in p.get("throughput", p.get("throughput_grid", [])):
        if isinstance(e, list):
            chains, threads = e
            entries.append(dict(id=f"K{chains}xt{threads}", chains=chains, threads=threads))
        else:
            entries.append(dict(e))
    if "throughput_grid" in p:
        # tall-1: the tempered sets at the main design's thread counts
        for cid, c in config["main"]["configs"].items():
            if c.get("tempered"):
                entries.append(dict(id=cid, chains=c["chains"], tempered=c["tempered"], threads=c["threads"]))
    return entries


def throughput_rows(config):
    p = config["pilot"]
    rows = []
    for rep in range(1, p["repeats"] + 1):
        for e in throughput_entries(config):
            r = row(run_id=f"p1_{e['id']}_r{rep}", phase="throughput", config=e["id"],
                    chains=e["chains"], tempered=e.get("tempered", 1), threads=e["threads"],
                    seed=100 + rep, minutes=p["minutes"], log_every=p["log_every"],
                    **{k: e.get(k) for k in DEVICE_FIELDS})
            if e.get("tempered"):
                # a provisional ladder; only the speed matters here
                r["ladder_ratio"] = "0.99"
                r["swap_every"] = "1000"
            rows.append(r)
    if "jvm_minutes" in p:
        # tall-1 P3: 16 chains in one JVM (two collectors) against 16 separate JVMs
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


def config_fields(c):
    return {k: c.get(k) for k in DEVICE_FIELDS}


def speed_key(cid, c, summary):
    if c.get("speed_key"):
        return c["speed_key"]
    threads = summary["t_star"] if c["threads"] == "t_star" else c["threads"]
    return cid if c.get("tempered") else f"K{c['chains']}xt{threads}"


def mixing_rows(config, summary):
    m = config["pilot"]["mixing"]
    cid = m["config"]
    c = config["main"]["configs"][cid]
    per_chain = speed(summary, speed_key(cid, c, summary))
    states = int(per_chain * m["hours"] * 3600 * 1.02)
    return [row(run_id=f"mixing_{cid}", phase="mixing", config=cid, chains=c["chains"],
                tempered=c.get("tempered", 1), threads=c["threads"], seed=m["seed"],
                minutes=m["hours"] * 60, states=states,
                log_every=max(1000, int(round(states / config["main"]["samples_per_chain"], -3))),
                **config_fields(c))]


def ladder_ratios(config, summary):
    t = config["tempering"]
    if "ladder_sd_factors" in t:
        sd = summary.get("joint_sd")
        if not sd:
            sys.exit("no joint SD in pilot/summary.json: run the mixing pilot and analyse.py pilot first")
        # neighbouring rungs swap well when the change in beta times the joint's SD is about 1
        return [round(1 - f / sd, 5) for f in t["ladder_sd_factors"]]
    return t["pilot_ladder_ratios"]


def tempering_rows(config, summary):
    t = config["tempering"]
    p = config["pilot"]
    minutes = p["tempering_minutes"]
    rows = []
    for cid in p.get("tempering_configs", ["T1"]):
        c = config["main"]["configs"][cid]
        per_chain = speed(summary, speed_key(cid, c, summary))
        swap = max(100, int(round(per_chain * t["swap_block_seconds"], -2)))
        for ratio in ladder_ratios(config, summary):
            rows.append(row(run_id=f"p2_{cid}_r{ratio}", phase="tempering", config=cid, chains=c["chains"],
                            tempered=c["tempered"], threads=c["threads"], seed=200, minutes=minutes,
                            states=int(per_chain * minutes * 60), log_every=p["log_every"],
                            ladder_ratio=ratio, swap_every=swap, **config_fields(c)))
    return rows


def main_rows(config, summary):
    m = config["main"]
    ratio = config["tempering"]["ladder_ratio"] or summary.get("ladder_ratio")
    if any(c.get("tempered") for c in m["configs"].values()) and not ratio:
        sys.exit('set "tempering": {"ladder_ratio": ...} in study.json (from the tempering pilot)')
    seconds = m["hours"] * 3600
    included = summary.get("include", {})
    blocks = []
    for seed in m["seeds"]:
        block = []
        for cid, c in m["configs"].items():
            if c.get("include_if") and not included.get(c["include_if"], False):
                continue  # the pilot did not justify this configuration
            threads = summary["t_star"] if c["threads"] == "t_star" else c["threads"]
            per_chain = speed(summary, speed_key(cid, c, summary))
            states = int(per_chain * seconds * 1.02)  # finishes at about the time limit
            log_every = max(1000, int(round(states / m["samples_per_chain"], -3)))
            fields = config_fields(c)
            chains = c["chains"]
            if c.get("numa_split_if") and included.get(c["numa_split_if"], False):
                # one JVM per NUMA node, each with its share of the chains
                nodes = config["machine"]["numa_nodes"]
                fields.update(processes=nodes, numa=1)
                chains = chains // nodes
            r = row(run_id=f"{cid}_s{seed}", phase="main", config=cid, chains=chains,
                    tempered=c.get("tempered", 1), threads=threads, seed=seed,
                    minutes=m["hours"] * 60, states=states, log_every=log_every, **fields)
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
    phases = ("throughput", "mixing", "tempering", "main")
    if len(sys.argv) < 2 or sys.argv[1] not in phases:
        sys.exit(__doc__)
    config = cc.load_config()
    phase = sys.argv[1]
    if phase == "throughput":
        rows = throughput_rows(config)
    else:
        summary = cc.load_pilot_summary()
        if summary is None:
            sys.exit("no pilot/summary.json: run the throughput pilot, then python3 analyse.py pilot")
        rows = {"mixing": mixing_rows, "tempering": tempering_rows, "main": main_rows}[phase](config, summary)
    all_rows = cc.write_conditions(rows, append=True)
    added = [r for r in all_rows if r["phase"] == phase]
    hours = sum(float(r["minutes"]) for r in added) / 60  # each run has the machine to itself
    print(f"{phase}: {len(added)} runs, about {hours:.1f} hours of machine time -> {cc.paths()['conditions']}")


if __name__ == "__main__":
    main()
