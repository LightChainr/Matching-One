#!/usr/bin/env python3
"""Canonical order crossings of twelve already-solved fixed physical clocks.

Exact Bernstein sign/root certificates consume saved safe-set coefficients.
No network reconstruction, reliability solve, Monte Carlo, or new selection.
"""
from fractions import Fraction as F
from functools import cmp_to_key
from hashlib import sha256
from itertools import combinations
import json
from math import comb, gcd, lcm
from pathlib import Path
import subprocess

import numpy as np
from scipy.optimize import brentq
from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "bd95f2a048d5780568b689bd42e0a684daf74315"
SOURCE_PATH = "results/p334-twelve-prefix-clocks/full_clocks.json"
OUT = ROOT/"results/p334-twelve-canonical-crossings"
N, K0, D = 425, 252, 173
P_REF = .59274605079
OFFSET = 43042500000


def signs(values):
    return [1 if x > 0 else -1 for x in values if x]


def variations(values):
    s = signs(values)
    return sum(a != b for a, b in zip(s, s[1:]))


def sign_runs(values):
    runs = []
    for i, x in enumerate(values):
        s = 1 if x > 0 else -1 if x < 0 else 0
        if runs and runs[-1]["sign"] == s:
            runs[-1]["last"] = i
        else:
            runs.append({"first": i, "last": i, "sign": s})
    return runs


def primitive_integer(values):
    factor = 0
    for x in values:
        factor = gcd(factor, x)
    return [x//factor for x in values] if factor else values


def split_half(values):
    row, degree = list(values), len(values)-1
    left, right = [], []
    for k in range(degree+1):
        left.append(row[0] << (degree-k))
        right.append(row[-1] << (degree-k))
        row = [a+b for a, b in zip(row, row[1:])]
    return primitive_integer(left), primitive_integer(right[::-1])


def root_certificate(coefficients):
    """0-variation leaves exclude roots; 1-variation leaves have exactly one.

    Exact double/dyadic-zero exceptional cases are not silently classified.
    None occurs in these frozen 66 comparisons.
    """
    pending = [(coefficients, F(0), F(1), 0)]
    leaves, intervals = [], []
    while pending:
        c, a, b, depth = pending.pop()
        v = variations(c)
        if v <= 1:
            leaves.append({"interval": [str(a), str(b)], "variations": v,
                           "first_nonzero_sign": signs(c)[0] if signs(c) else 0,
                           "last_nonzero_sign": signs(c)[-1] if signs(c) else 0})
            if v == 1:
                intervals.append((a, b))
            continue
        if depth >= 40:
            raise ValueError("Unresolved exact root multiplicity; do not claim a crossing count")
        left, right = split_half(c)
        if left[-1] == 0:
            raise ValueError("Exact dyadic zero requires a separate multiplicity certificate")
        mid = (a+b)/2
        pending.extend([(right, mid, b, depth+1), (left, a, mid, depth+1)])
    return sorted(intervals), sorted(leaves, key=lambda row: F(row["interval"][0]))


def thermal_survival_coefficients(row):
    survival = [F(x, comb(D, k)) for k, x in enumerate(row["true_safe_counts"])]
    return np.r_[np.ones(K0), [float(x) for x in survival]]


def evaluate(c, p):
    return float(c@binom.pmf(np.arange(len(c)), len(c)-1, p))


def derivative(c, p):
    return evaluate((len(c)-1)*np.diff(c), p)


def reduced_difference(a, b):
    # Delta=f2_a-f2_b has coefficients S_b-S_a.
    difference = [y-x for x, y in zip(a["true_safe_counts"], b["true_safe_counts"])]
    nonzero = [k for k, x in enumerate(difference) if x]
    if not nonzero:
        raise ValueError("Identical clocks need an equality label, not root isolation")
    low, high = nonzero[0], nonzero[-1]
    degree = high-low
    reduced = [F(difference[k]*comb(N, K0+k), comb(D, k)*comb(degree, k-low))
               for k in range(low, high+1)]
    denominator = lcm(*(x.denominator for x in reduced))
    integer = primitive_integer([int(x*denominator) for x in reduced])
    full = np.r_[np.zeros(K0), [float(F(x, comb(D, k))) for k, x in enumerate(difference)]]
    return difference, low, high, integer, full


def pair_score(a, b, curve_a, curve_b, survival_a, survival_b):
    difference, low, high, integer, full = reduced_difference(a, b)
    intervals, leaves = root_certificate(integer)
    magnitude = max(map(abs, integer))
    normalized = np.array([float(F(x, magnitude)) for x in integer])
    roots = []
    leading = signs(integer)[0]
    for number, (lo, hi) in enumerate(intervals):
        p = float(brentq(lambda x: evaluate(normalized, x), float(lo), float(hi), xtol=1e-14))
        roots.append({"p": p, "isolation_interval": [str(lo), str(hi)],
                      "delta_derivative": derivative(full, p),
                      "first_F2": evaluate(curve_a, p), "second_F2": evaluate(curve_b, p),
                      "first_one_minus_F2": evaluate(survival_a, p),
                      "second_one_minus_F2": evaluate(survival_b, p),
                      "delta_sign_below": leading*(-1)**number,
                      "delta_sign_above": -leading*(-1)**number})
    area_delta = -(F(a["mean_true_birth_step"]["exact"])-F(b["mean_true_birth_step"]["exact"]))/426
    return {"first_counter": a["counter"], "second_counter": b["counter"],
            "delta_definition": "F2_first-F2_second", "survival_first_minus_second_sign_runs": sign_runs([-x for x in difference]),
            "rank_sign_variations": variations(difference), "thermal_root_count": len(roots), "roots": roots,
            "reduction": {"positive_endpoint_factor": f"p^{K0+low}*(1-p)^{N-K0-high}",
                          "remaining_Bernstein_degree": high-low,
                          "primitive_integer_coefficient_sha256": sha256(json.dumps(integer, separators=(",", ":")).encode()).hexdigest(),
                          "initial_sign_variations": variations(integer), "complete_dyadic_partition": leaves},
            "p_ref": {"p": P_REF, "delta": evaluate(full, P_REF), "delta_derivative": derivative(full, P_REF),
                      "first_F2": evaluate(curve_a, P_REF), "second_F2": evaluate(curve_b, P_REF),
                      "first_derivative": derivative(curve_a, P_REF), "second_derivative": derivative(curve_b, P_REF)},
            "integrated_delta": {"exact": str(area_delta), "float": float(area_delta)},
            "initial_delta_sign": leading}


def main():
    raw = subprocess.check_output(["git", "show", SOURCE_COMMIT+":"+SOURCE_PATH], cwd=ROOT)
    records = sorted(json.loads(raw)["records"], key=lambda row: row["counter"])
    assert len(records) == 12 and all(row["status"] == "solved_full_physical" for row in records)
    assert all(int(row["original_row"]["k0"]) == K0 for row in records)
    survivals = {row["counter"]: thermal_survival_coefficients(row) for row in records}
    curves = {counter: 1-survival for counter, survival in survivals.items()}
    pairs = [pair_score(a, b, curves[a["counter"]], curves[b["counter"]],
                        survivals[a["counter"]], survivals[b["counter"]])
             for a, b in combinations(records, 2)]
    pair_map = {(row["first_counter"], row["second_counter"]): row for row in pairs}
    representative = pair_map[(OFFSET+83, OFFSET+1006)]
    crossing = [row for row in pairs if row["rank_sign_variations"]]
    events = sorted([{"p": root["p"], "pair": [row["first_counter"], row["second_counter"]]}
                     for row in crossing for root in row["roots"]], key=lambda row: row["p"])
    def compare(a, b):
        if a == b:
            return 0
        if a < b:
            return -pair_map[(a, b)]["initial_delta_sign"]
        return pair_map[(b, a)]["initial_delta_sign"]
    order = sorted(curves, key=cmp_to_key(compare))
    cells, lower = [], 0.
    for event in events:
        cells.append({"p_open_interval": [lower, event["p"]], "F2_descending_counters": list(order)})
        a, b = map(order.index, event["pair"])
        if abs(a-b) != 1:
            raise ValueError("A nonadjacent order swap would require unresolved simultaneous crossings")
        order[a], order[b] = order[b], order[a]
        lower = event["p"]
    cells.append({"p_open_interval": [lower, 1.], "F2_descending_counters": order})
    result = {"schema": "matching-one/p334-twelve-canonical-clock-crossings/v1",
              "source_commit": SOURCE_COMMIT, "source_path": SOURCE_PATH, "source_sha256": sha256(raw).hexdigest(),
              "new_prefixes": 0, "new_MC": 0, "network_or_reliability_reruns": 0,
              "fixed_parameters": {"N": N, "k0": K0, "vacant_sites": D, "p_ref": P_REF},
              "selected_count": 12, "fixed_pair_count": len(pairs),
              "summary": {"rank_crossing_pairs": len(crossing),
                          "thermal_crossing_pairs": sum(row["thermal_root_count"] > 0 for row in pairs),
                          "disappeared_rank_crossing_pairs": sum(row["thermal_root_count"] == 0 for row in crossing),
                          "total_thermal_crossings": len(events), "order_cells": len(cells)},
              "representative_83_minus_1006": representative,
              "pairs": pairs, "canonical_order_cells": cells,
              "boundary": "Twelve fixed old-source prefixes, not a population crossing rate. Exact Bernstein certificates determine root counts; root positions and slopes are numerical evaluations. Very-late crossings can have negligible absolute signal. All p values reuse the same exact conditional laws; no new independent evidence or per-prefix covariance rescore is claimed."}
    OUT.mkdir(exist_ok=True)
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# The twelve physical clocks retain their thermal order reversals", "",
             f"83/1006 has exactly one canonical crossing at p={representative['roots'][0]['p']:.13g}.",
             f"For Delta=F2_83-F2_1006, Delta(p_ref)={representative['p_ref']['delta']:.10g}, Delta'(p_ref)={representative['p_ref']['delta_derivative']:.10g}.", "",
             "| Prefix pair (counter suffixes) | Rank sign changes | Canonical roots |", "|---|---:|---|"]
    for row in crossing:
        lines.append(f"| {row['first_counter']-OFFSET}/{row['second_counter']-OFFSET} | {row['rank_sign_variations']} | "+", ".join(f"{r['p']:.10g}" for r in row['roots'])+" |")
    lines += ["", f"{result['summary']}", "", result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
