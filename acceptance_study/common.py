"""
Shared helpers for the operator study: configuration, run naming, and the
parsers and statistics used by run_grid.py and analyse.py.

ESS follows dr.inference.trace.TraceCorrelation exactly (Tracer's estimator),
so numbers here are comparable with what users see in Tracer/loganalyser.
"""

import json
import math
import os
import re
from pathlib import Path

import numpy as np

STUDY_DIR = Path(__file__).resolve().parent
REPO_DIR = STUDY_DIR.parent
CONFIG_PATH = STUDY_DIR / "study.json"

# Log columns that are not model parameters (operator diagnostics, state)
NON_PARAMETER_COLUMNS = {"state", "topologyChanged", "mean_delta"}


def load_config() -> dict:
    """study.json, or the file named by $OPERATOR_STUDY_CONFIG (e.g. for a smoke test)."""
    path = Path(os.environ.get("OPERATOR_STUDY_CONFIG", CONFIG_PATH))
    with open(path) as f:
        return json.load(f)


def work_dir(config: dict) -> Path:
    """Where runs/, reference/, analysis/ and conditions.csv live (study.json "work_dir")."""
    return STUDY_DIR / config.get("work_dir", ".")


def base_xml_path(config: dict, dataset: str, likelihood: str) -> Path:
    return REPO_DIR / config["datasets"][dataset]["xml"].format(likelihood=likelihood)


def format_level(value) -> str:
    """A filesystem-safe rendering of a factor level for run ids."""
    if value is None or value == "":
        return "default"
    text = str(value)
    return re.sub(r"[^A-Za-z0-9.]+", "-", text).strip("-")


# ---------------------------------------------------------------------------
# ESS (TraceCorrelation.analyseCorrelationNumeric)
# ---------------------------------------------------------------------------

MAX_LAG = 2000


def ess(values, step_size: int = 1) -> float:
    x = np.asarray(values, dtype=float)
    samples = len(x)
    if samples < 2 or np.all(x == x[0]):
        return float("nan")
    x = x - x.mean()
    max_lag = min(samples - 1, MAX_LAG)

    # autocovariance at every lag via FFT, normalised by (samples - lag) as Tracer does
    n_fft = 1 << (2 * samples - 1).bit_length()
    f = np.fft.rfft(x, n_fft)
    acov = np.fft.irfft(f * np.conjugate(f), n_fft)[:max_lag]
    gamma = acov / (samples - np.arange(max_lag))

    var_stat = gamma[0]
    for lag in range(2, max_lag, 2):
        pair = gamma[lag - 1] + gamma[lag]
        if pair > 0:
            var_stat += 2.0 * pair
        else:
            break

    if gamma[0] == 0:
        return 1.0
    act = step_size * var_stat / gamma[0]
    return (step_size * samples) / act if act != 0 else 1.0


# ---------------------------------------------------------------------------
# BEAST .log files
# ---------------------------------------------------------------------------

def read_log(path: Path):
    """Returns (columns, 2-D numpy array) from a tab-separated BEAST log."""
    columns = None
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if columns is None:
                columns = [p.strip() for p in parts]
                continue
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    return columns, np.array(rows)


# ---------------------------------------------------------------------------
# Operator analysis (.ops) files
# ---------------------------------------------------------------------------

def _is_number(token: str) -> bool:
    try:
        float(token)
        return True
    except ValueError:
        return False


def read_ops(path: Path) -> dict:
    """Parses the plain-text operator analysis table, keyed by operator name.

    Rows end in Count, Time(s), Time/Op(us), Pr(accept); before those come the
    weight and, for adaptive operators, the tuning value.
    """
    ops = {}
    if not path.exists():
        return ops
    for line in path.read_text().splitlines():
        tokens = line.split()
        if len(tokens) < 6 or tokens[0] == "Operator" or not all(_is_number(t) for t in tokens[-4:]):
            continue
        count, time_s, per_op, accept = tokens[-4:]
        head = tokens[:-4]
        numbers = []
        while head and _is_number(head[-1]) and len(numbers) < 2:
            numbers.insert(0, float(head.pop()))
        if not head or not numbers:
            continue
        ops[" ".join(head)] = {
            "weight": numbers[0],
            "tuning": numbers[1] if len(numbers) > 1 else None,
            "count": int(count),
            "time_s": float(time_s),
            "time_per_op_us": float(per_op),
            "accept": float(accept),
        }
    return ops


# ---------------------------------------------------------------------------
# NEXUS tree files -> clade sets
# ---------------------------------------------------------------------------

TRANSLATE_RE = re.compile(r"^\s*(\d+)\s+([^,;]+?)\s*[,;]?\s*$")
TREE_RE = re.compile(r"^\s*tree\s+STATE_(\d+)\s*(?:\[[^\]]*\]\s*)?=\s*(.*;)\s*$", re.IGNORECASE)
COMMENT_RE = re.compile(r"\[[^\]]*\]")
TOKEN_RE = re.compile(r"\(|\)|,|;|:[^(),;]*|[^(),;:]+")


def newick_clades(newick: str, taxon_index: dict) -> set:
    """Non-trivial clades (as integer bitmasks over taxon indices) of a rooted tree."""
    clades = set()
    stack = []
    current = 0
    for token in TOKEN_RE.findall(COMMENT_RE.sub("", newick)):
        if token == "(":
            stack.append(current)
            current = 0
        elif token == ")":
            child = current
            if bin(child).count("1") > 1:
                clades.add(child)
            current = stack.pop() | child
        elif token in (",", ";") or token.startswith(":") or not token.strip():
            continue
        else:
            current |= 1 << taxon_index[token.strip()]
    return clades


def read_trees(path: Path):
    """Yields (state, clade set) for each tree, plus taxa names via the first return value.

    Returns (taxa, [(state, clades)]) where taxa is ordered by translate number.
    """
    translate = {}
    taxa = []
    trees = []
    in_translate = False
    taxon_index = None
    with open(path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.lower() == "translate":
                in_translate = True
                continue
            if in_translate:
                m = TRANSLATE_RE.match(line)
                if m:
                    translate[m.group(1)] = m.group(2).strip("'\"")
                if stripped.endswith(";"):
                    in_translate = False
                    keys = sorted(translate, key=int)
                    taxa = [translate[k] for k in keys]
                    taxon_index = {k: i for i, k in enumerate(keys)}
                continue
            m = TREE_RE.match(line)
            if m:
                if taxon_index is None:
                    raise ValueError(f"{path}: trees without a translate block are not supported")
                trees.append((int(m.group(1)), newick_clades(m.group(2), taxon_index)))
    return taxa, trees


def split_frequencies(tree_sets) -> dict:
    """Pooled clade frequencies, each list of clade sets weighted equally."""
    freqs = {}
    for clade_sets in tree_sets:
        w = 1.0 / (len(clade_sets) * len(tree_sets))
        for clades in clade_sets:
            for c in clades:
                freqs[c] = freqs.get(c, 0.0) + w
    return freqs


def split_deviation(freqs: dict, reference: dict, min_freq: float = 0.01):
    """(ASDSF, max deviation) between two clade frequency tables."""
    keys = {c for c, f in freqs.items() if f >= min_freq} | {c for c, f in reference.items() if f >= min_freq}
    if not keys:
        return 0.0, 0.0
    devs = [abs(freqs.get(c, 0.0) - reference.get(c, 0.0)) for c in keys]
    return float(np.mean(devs)), float(np.max(devs))


def clades_to_json(freqs: dict) -> dict:
    return {format(c, "x"): f for c, f in freqs.items()}


def clades_from_json(data: dict) -> dict:
    return {int(k, 16): f for k, f in data.items()}


def split_operator_name(name: str):
    """'stl1:subtreeLeap(treeModel)' -> ('stl1', 'subtreeLeap(treeModel)'); no id -> (None, name).

    BEAST prefixes an operator's name with its XML id when it has one."""
    head = name.split("(")[0]
    if ":" in head:
        op_id, base = name.split(":", 1)
        return op_id, base
    return None, name


def operator_kind(name: str) -> str:
    """'stl1:subtreeLeap(treeModel)' -> 'subtreeLeap'."""
    return split_operator_name(name)[1].split("(")[0].strip()


def geometric_mean(values):
    v = [x for x in values if x is not None and x > 0 and math.isfinite(x)]
    return float(np.exp(np.mean(np.log(v)))) if v else float("nan")
