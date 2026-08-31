#!/usr/bin/env python3
"""Sharp fixed-bipartite-capacity overlap envelopes on saved trigger graphs."""
from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import json
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/local-20260831/P334-cooperative-closure/trigger_graph_structure_bounded.json"
OUTPUT = ROOT / "results/p334-sharp-overlap-envelope"


def frac(value):
    value = Fraction(value)
    return {"exact": str(value), "value": float(value)}


def balanced_choose_sum(vertices, edges):
    quotient, remainder = divmod(edges, vertices)
    return vertices * comb(quotient, 2) + remainder * quotient


@lru_cache(None)
def sharp_wedges(L, R, m):
    """All simple bipartite graphs on fixed L/R slots, with exactly m edges.

    Isolates/disconnection are allowed. Row compression proves a maximum is
    Ferrers; the recurrence searches only its integer partition, not graphs.
    """
    if not 0 <= m <= L * R or min(L, R) < 1:
        raise ValueError("infeasible capacity")
    if L > R:
        result = sharp_wedges(R, L, m)
        return {**result, "L": L, "R": R, "maximum_partition_side": "right"}
    minimum = balanced_choose_sum(L, m) + balanced_choose_sum(R, m)

    @lru_cache(None)
    def optimum(row, cap, remaining):
        if row == L:
            return (0, ()) if remaining == 0 else (-10**30, ())
        if remaining < 0 or remaining > (L - row) * cap:
            return -10**30, ()
        best = (-10**30, ())
        lower = (remaining + L - row - 1) // (L - row)
        for degree in range(lower, min(cap, remaining) + 1):
            score, suffix = optimum(row + 1, degree, remaining - degree)
            candidate = (comb(degree, 2) + row * degree + score, (degree,) + suffix)
            if candidate > best:
                best = candidate
        return best

    maximum, partition = optimum(0, R, m)
    return {"L": L, "R": R, "m": m, "minimum": minimum, "maximum": maximum,
            "maximum_ferrers_partition": list(partition), "maximum_partition_side": "left"}


def coop_numerator(a, m, W):
    return Fraction(2 * m + 2 * W) - Fraction(4 * m * m, a)


def score_record(record):
    structure, source = record["structure"], record["source_row"]
    a, m, observed = (structure[key] for key in ("safe_vertices", "trigger_edges", "trigger_wedges"))
    d = int(source["n"]) - int(source["k0"])
    denominator = d * (d - 1) ** 2
    components = []
    balanced_real = Fraction(0)
    within_side = Fraction(0)
    for component in structure["components_with_edges"]:
        L, R = map(len, component["sides"])
        edges = component["edges"]
        envelope = sharp_wedges(L, R, edges)
        W = sum(comb(degree, 2) for side in component["side_degree_sequences"] for degree in side)
        if not envelope["minimum"] <= W <= envelope["maximum"]:
            raise AssertionError("observed graph outside sharp envelope")
        balanced_real += Fraction(edges * edges, L) + Fraction(edges * edges, R)
        within = sum(Fraction((size * degree - edges) ** 2, size * size)
                     for size, side in zip((L, R), component["side_degree_sequences"])
                     for degree in side)
        within_side += within
        components.append({**envelope, "observed": W, "within_side_degree_variation": frac(within)})
    lower = sum(c["minimum"] for c in components)
    upper = sum(c["maximum"] for c in components)
    actual_num = coop_numerator(a, m, observed)
    structural_num = balanced_real - Fraction(4 * m * m, a)
    if actual_num != structural_num + within_side:
        raise AssertionError("capacity/within-side decomposition mismatch")
    safe_triples_pair_only = comb(a, 3) - m * (a - 2) + observed
    return {"selection": record["selection"], "source_graph": record["graph_artifact"],
            "N": int(source["n"]), "orientation": source["orientation"],
            "counter": int(source["replica"]), "d": d, "a_safe": a, "m_trigger": m,
            "components": components,
            "W2": {"minimum": lower, "observed": observed, "maximum": upper,
                   "fraction_through_sharp_range": frac(Fraction(observed - lower, upper - lower)) if upper > lower else None},
            "Delta_coop": {"minimum": frac(coop_numerator(a, m, lower) / denominator),
                           "observed": frac(actual_num / denominator),
                           "maximum": frac(coop_numerator(a, m, upper) / denominator),
                           "sharp_minimum_fraction_of_observed": frac(coop_numerator(a, m, lower) / actual_num) if actual_num else None},
            "degree_variance_decomposition": {
                "capacity_floor_real": frac(structural_num), "within_side": frac(within_side),
                "total": frac(actual_num), "capacity_fraction": frac(structural_num / actual_num) if actual_num else None},
            "pair_only_safe_triples": {"minimum": comb(a, 3) - m * (a - 2) + lower,
                                      "observed": safe_triples_pair_only,
                                      "maximum": comb(a, 3) - m * (a - 2) + upper},
            "pair_only_s3_upper": frac(Fraction(safe_triples_pair_only, comb(d, 3)))}


def build():
    source = json.loads(SOURCE.read_text())
    rows = [score_record(r) for r in source["records"]]
    A, B = rows[:2]
    a, d, m = A["a_safe"], A["d"], A["m_trigger"]
    assert (a, d, m) == (B["a_safe"], B["d"], B["m_trigger"])
    triple_path = ROOT / "results/local-20260831/P334-cooperative-closure/safe_triple_census.json"
    triple = {r["replica_counter"]: r for r in json.loads(triple_path.read_text())["checkpoints"]}
    cA, cB = (triple[r["counter"]]["minimal_nonfaces_size3"] for r in (A, B))
    safe_pairs = comb(a, 2) - m
    for row in (A, B):
        count = triple[row["counter"]]["actual_safe_triples"]
        row["actual_s3"] = frac(Fraction(count, comb(d, 3)))
        row["actual_h3_conditional_on_two_safe"] = frac(1 - Fraction(3 * count, (d - 2) * safe_pairs))
    possible_difference = {
        "minimum_W2_B_minus_A": B["W2"]["minimum"] - A["W2"]["maximum"],
        "observed_W2_B_minus_A": B["W2"]["observed"] - A["W2"]["observed"],
        "maximum_W2_B_minus_A": B["W2"]["maximum"] - A["W2"]["minimum"],
    }
    minimum_gap = possible_difference["minimum_W2_B_minus_A"]
    possible_difference.update({
        "branch_survival_gap_minimum": frac(Fraction(2 * minimum_gap, d * (d - 1)**2)),
        "actual_s3_gap_minimum_if_c3_held_fixed": frac(Fraction(minimum_gap + cA - cB, comb(d, 3))),
        "h3_gap_B_minus_A": frac(-Fraction(3 * (possible_difference["observed_W2_B_minus_A"] + cA - cB), (d - 2) * safe_pairs)),
        "h3_gap_upper_if_c3_held_fixed": frac(-Fraction(3 * (minimum_gap + cA - cB), (d - 2) * safe_pairs)),
        "observed_s3_gap": frac(Fraction(B["W2"]["observed"] - A["W2"]["observed"] + cA - cB, comb(d, 3))),
        "c3_A": cA, "c3_B": cB,
        "fraction_of_observed_W2_gap_forced_for_every_rewiring": frac(Fraction(minimum_gap, possible_difference["observed_W2_B_minus_A"]))})
    floor_fractions = [Fraction(r["Delta_coop"]["sharp_minimum_fraction_of_observed"]["exact"]) for r in rows]
    return {"parent_commit": "119cb5fc87a80285413e9bf019b2a6d5aa257c38", "new_samples": 0,
            "selection": source["selection_rule"], "records": rows, "saved_N425_pair": possible_difference,
            "selected_summary": {"graphs": len(rows), "zero_width_envelopes": sum(r["W2"]["minimum"] == r["W2"]["maximum"] for r in rows),
                                 "minimum_floor_fraction": frac(min(floor_fractions)), "maximum_floor_fraction": frac(max(floor_fractions))},
            "envelope_scope": "Sharp over simple bipartite graphs with fixed observed component vertex slots, L/R capacities and edge counts; inter-component edges are forbidden, but isolates/disconnection within each fixed block are allowed. Extremizers are not asserted to be realizable square-lattice checkpoints.",
            "population_boundary": "Only the existing two witnesses plus first five eligible counters per environment; 22 selected configurations, not all checkpoints or independent evidence."}


if __name__ == "__main__":
    result = build()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "sharp_envelopes.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    for row in result["records"]:
        print(row["selection"], row["W2"], row["Delta_coop"]["sharp_minimum_fraction_of_observed"])
    print("N425", json.dumps(result["saved_N425_pair"], sort_keys=True))
