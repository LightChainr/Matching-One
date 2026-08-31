#!/usr/bin/env python3
"""Exact pair+triple truncated clocks for the two archived N425 checkpoints.

No new random samples. False-twin site groups are weighted by (1+z)^s-1;
weighted deletion/contraction counts every allowed subset exactly once.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
import json
from math import comb
from pathlib import Path
import time

from p334_pair_only_survival import contiguous, frac, multiply

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/p334-pair-triple-clock"
OLD = ROOT / "results/local-20260831/P334-cooperative-closure"


def add(p, q):
    out = [0] * max(len(p), len(q))
    for i, value in enumerate(p):
        out[i] += value
    for i, value in enumerate(q):
        out[i] += value
    return tuple(out)


def binomial(n, nonempty=False):
    return tuple(0 if nonempty and k == 0 else comb(n, k) for k in range(n + 1))


def strip_one_plus_z(coefficients):
    """Exact repeated division; returns a short dominance certificate."""
    p = list(coefficients)
    while p and p[-1] == 0:
        p.pop()
    power = 0
    while len(p) > 1 and sum(v * (-1)**i for i, v in enumerate(p)) == 0:
        q = [0] * (len(p) - 1)
        q[-1] = p[-1]
        for i in range(len(q) - 1, 0, -1):
            q[i - 1] = p[i] - q[i]
        p, power = q, power + 1
    return power, p


def polynomial(graph, triples):
    edges = [frozenset(e) for e in graph["minimal_trigger_pairs"] + triples]
    sites = set(graph["safe_sites"])
    active = set().union(*edges)
    # Equal links imply no edge contains two group members. Any nonempty
    # subset of such a group imposes the same constraints on other groups.
    links = defaultdict(list)
    for v in sorted(active):
        link = tuple(sorted(tuple(sorted(e - {v})) for e in edges if v in e))
        links[link].append(v)
    groups = list(links.values())
    index = {v: j for j, group in enumerate(groups) for v in group}
    compressed = tuple(sorted({sum(1 << index[v] for v in e) for e in edges}))
    weights = [binomial(len(group), nonempty=True) for group in groups]
    started = time.monotonic()
    calls = 0

    @lru_cache(None)
    def solve(vertices, constraints):
        nonlocal calls
        calls += 1
        if calls > 200_000 or time.monotonic() - started > 90:
            raise RuntimeError("bounded exact-state budget exceeded")
        if 0 in constraints:
            return (0,)
        forced_absent = 0
        for edge in constraints:
            if edge.bit_count() == 1:
                forced_absent |= edge
        if forced_absent:
            return solve(vertices & ~forced_absent,
                         tuple(e for e in constraints if not e & forced_absent))
        # A pair forbids its triple supersets too; remove those redundancies.
        pairs = {e for e in constraints if e.bit_count() == 2}
        reduced = []
        for edge in constraints:
            if edge.bit_count() == 3:
                bits, remaining = [], edge
                while remaining:
                    bit = remaining & -remaining
                    bits.append(bit)
                    remaining ^= bit
                if any(edge ^ bit in pairs for bit in bits):
                    continue
            reduced.append(edge)
        reduced = tuple(reduced)
        if reduced != constraints:
            return solve(vertices, reduced)
        touched = 0
        for edge in constraints:
            touched |= edge
        isolated = vertices & ~touched
        if isolated:
            degree = sum(len(groups[j]) for j in range(len(groups)) if isolated & (1 << j))
            return tuple(multiply(binomial(degree), solve(touched, constraints)))
        if not constraints:
            return (1,)
        # Independent hypergraph components multiply exactly.
        pending = vertices
        components = []
        while pending:
            component = pending & -pending
            while True:
                grown = component
                for edge in constraints:
                    if component & edge:
                        grown |= edge
                if grown == component:
                    break
                component = grown
            components.append(component)
            pending &= ~component
        if len(components) > 1:
            result = (1,)
            for component in components:
                result = tuple(multiply(result, solve(component, tuple(e for e in constraints if e & component))))
            return result
        # High incidence reduces the two bounded real hypergraphs rapidly.
        chosen = max((j for j in range(len(groups)) if vertices & (1 << j)),
                     key=lambda j: sum(4 - e.bit_count() for e in constraints if e & (1 << j)))
        bit = 1 << chosen
        absent = solve(vertices ^ bit, tuple(e for e in constraints if not e & bit))
        present = solve(vertices ^ bit, tuple(sorted({e & ~bit for e in constraints})))
        return add(absent, multiply(weights[chosen], present))

    core = solve((1 << len(groups)) - 1, compressed)
    isolates = len(sites - active)
    counts = multiply(core, binomial(isolates))
    counts += [0] * (len(sites) + 1 - len(counts))
    return counts, {"active_sites": len(active), "isolated_sites": isolates,
                    "false_twin_groups": groups, "compressed_hyperedges": len(compressed),
                    "memoized_states": calls, "elapsed_seconds": time.monotonic() - started,
                    "core_independent_counts": list(core)}


def main():
    source = json.loads((OUTPUT / "full_triples.json").read_text())
    pair_only = {r["counter"]: r for r in json.loads((ROOT / "results/p334-pair-only-clock/pair_only_survival.json").read_text())["records"]}
    records = []
    for triple_row in source["checkpoints"]:
        counter = triple_row["replica_counter"]
        artifact = OLD / "trigger_graph_raw" / f"N425-second-{counter}.json"
        graph = json.loads(artifact.read_text())
        counts, structure = polynomial(graph, triple_row["all_minimal_nonfaces"])
        d = len(graph["safe_sites"])
        if counts[:4] != [1, d, triple_row["b2_safe_pairs"], triple_row["actual_safe_triples"]]:
            raise ValueError("new polynomial does not start with the archived physical counts")
        survival = [Fraction(counts[k], comb(d, k)) for k in range(d + 1)]
        hazard = [None] + [1 - survival[k] / survival[k - 1] if survival[k - 1] else None for k in range(1, d + 1)]
        pair = pair_only[counter]
        mean = sum(survival[:-1])
        row = {"counter": counter, "source_graph": str(artifact.relative_to(ROOT)),
               "seed": graph["seed"], "N": graph["N"], "k0": graph["k0"], "candidate_sites": d,
               "genuine_triples": len(triple_row["all_minimal_nonfaces"]), **structure,
               "independent_counts": counts, "survival": [frac(s) for s in survival],
               "hazard": [frac(h) if h is not None else None for h in hazard],
               "maximum_pair_triple_safe_k": max(k for k, count in enumerate(counts) if count),
               "mean_first_pair_triple_trigger_step": frac(mean),
               "mean_clock_shortening_from_pairs": frac(Fraction(pair["mean_first_pair_trigger_step"]["exact"]) - mean),
               "first_trigger_quantiles": {str(q): next(k for k in range(d + 1) if survival[k] <= 1 - q)
                                           for q in (Fraction(1, 10), Fraction(1, 2), Fraction(9, 10))},
               "removed_pair_safe_counts": [pair["independent_counts"][k] - counts[k] for k in range(d + 1)]}
        records.append(row)
        print(counter, "states", row["memoized_states"], "groups", len(row["false_twin_groups"]),
              "seconds", row["elapsed_seconds"], "mean", row["mean_first_pair_triple_trigger_step"],
              "maximum", row["maximum_pair_triple_safe_k"], flush=True)
    A, B = records
    difference = [b - a for a, b in zip(A["independent_counts"], B["independent_counts"])]
    crosses = {k: A["independent_counts"][k] * B["independent_counts"][k - 1]
                 - B["independent_counts"][k] * A["independent_counts"][k - 1]
               for k in range(1, d + 1) if A["independent_counts"][k - 1] and B["independent_counts"][k - 1]}
    factors = [strip_one_plus_z(row["independent_counts"]) for row in records]
    common = min(power for power, _ in factors)
    residuals = [multiply(p, binomial(power - common)) for power, p in factors]
    width = max(map(len, residuals))
    for p in residuals:
        p += [0] * (width - len(p))
    residual_difference = [b - a for a, b in zip(*residuals)]
    z_power = next(k for k, value in enumerate(residual_difference) if value)
    loss_difference = [b - a for a, b in zip(A["removed_pair_safe_counts"], B["removed_pair_safe_counts"])]
    new_gap = (Fraction(B["mean_first_pair_triple_trigger_step"]["exact"])
               - Fraction(A["mean_first_pair_triple_trigger_step"]["exact"]))
    old_gap = (Fraction(pair_only[B["counter"]]["mean_first_pair_trigger_step"]["exact"])
               - Fraction(pair_only[A["counter"]]["mean_first_pair_trigger_step"]["exact"]))
    comparison = {"independent_count_B_minus_A": difference,
                  "survival_B_above": contiguous([k for k, v in enumerate(difference) if v > 0]),
                  "survival_B_below": contiguous([k for k, v in enumerate(difference) if v < 0]),
                  "survival_equal": contiguous([k for k, v in enumerate(difference) if v == 0]),
                  "hazard_crossproducts_B_minus_A": crosses,
                  "hazard_B_above": contiguous([k for k, v in crosses.items() if v > 0]),
                  "hazard_B_below": contiguous([k for k, v in crosses.items() if v < 0]),
                  "hazard_equal": contiguous([k for k, v in crosses.items() if v == 0]),
                  "mean_first_trigger_B_minus_A": frac(new_gap),
                  "fraction_of_pair_only_mean_gap_retained": frac(new_gap / old_gap),
                  "triples_shortening_B_minus_A": frac(old_gap - new_gap),
                  "triple_loss_count_B_minus_A": loss_difference,
                  "triple_survival_loss_B_above": contiguous([k for k, v in enumerate(loss_difference) if v > 0]),
                  "triple_survival_loss_B_below": contiguous([k for k, v in enumerate(loss_difference) if v < 0]),
                  "dominance_factor": {"common_1_plus_z_power": common, "z_power": z_power,
                                       "remaining_coefficients": residual_difference[z_power:],
                                       "all_remaining_coefficients_strictly_positive": all(v > 0 for v in residual_difference[z_power:])}}
    result = {"parent_commit": "ad6c595a70c66ea4421c816b4c65b1cfe3d9c803", "new_samples": 0,
              "records": records, "comparison": comparison,
              "claim_boundary": "Exact truncated pair+genuine-triple trigger clocks of two fixed real checkpoints. All k<=3 match the full physical birth clock. At k>=4 omitted larger minimal triggers can shorten further. Neither hazard bounds nor full-birth A/B ordering follow from survival upper bounds."}
    (OUTPUT / "pair_triple_survival.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print({k: v for k, v in comparison.items() if k not in (
        "independent_count_B_minus_A", "hazard_crossproducts_B_minus_A",
        "triple_loss_count_B_minus_A", "dominance_factor")})


if __name__ == "__main__":
    main()
