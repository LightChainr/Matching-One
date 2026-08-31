#!/usr/bin/env python3
"""Exact independent-set clocks for the two saved N425 pair-trigger graphs."""
from __future__ import annotations

from collections import Counter
from fractions import Fraction
import json
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/local-20260831/P334-cooperative-closure"
OUTPUT = ROOT / "results/p334-pair-only-clock"


def frac(value):
    value = Fraction(value)
    return {"exact": str(value), "value": float(value)}


def multiply(p, q):
    result = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            result[i + j] += a * b
    return result


def independent_polynomial(graph, sides):
    left, right = sorted(sides, key=len)
    index = {v: i for i, v in enumerate(right)}
    adjacency = {v: set() for v in graph["safe_sites"]}
    for u, v in graph["minimal_trigger_pairs"]:
        adjacency[u].add(v)
        adjacency[v].add(u)
    masks = [sum(1 << index[v] for v in adjacency[u]) for u in left]
    coverage = [0] * (1 << len(left))
    histogram = Counter()
    core = [0] * (len(left) + len(right) + 1)
    for mask in range(1 << len(left)):
        if mask:
            bit = mask & -mask
            coverage[mask] = coverage[mask ^ bit] | masks[bit.bit_length() - 1]
        selected, free = mask.bit_count(), len(right) - coverage[mask].bit_count()
        histogram[selected, free] += 1
        for k in range(free + 1):
            core[selected + k] += comb(free, k)
    isolated = len(graph["safe_sites"]) - len(left) - len(right)
    polynomial = multiply(core, [comb(isolated, k) for k in range(isolated + 1)])
    universal = [v for v in adjacency if len(adjacency[v]) and
                 any(v in side and len(adjacency[v]) == len(other)
                     for side, other in (sides, list(reversed(sides))))]
    # Extract the real graph's optional universal-vertex/disjoint-Klr form.
    factorizations = []
    for u in universal:
        remaining = {v for v in adjacency if adjacency[v]} - {u}
        components = []
        while remaining:
            pending, component = [min(remaining)], set()
            while pending:
                v = pending.pop()
                if v in component:
                    continue
                component.add(v)
                pending.extend((adjacency[v] - {u}) & remaining)
            remaining -= component
            L, R = ([v for v in side if v in component] for side in sides)
            if all(adjacency[v] - {u} == set(R) for v in L) and all(adjacency[v] - {u} == set(L) for v in R):
                components.append({"left": L, "right": R, "shape": [len(L), len(R)]})
            else:
                break
        else:
            factorizations.append({"universal_site": u, "components_after_deletion": components,
                                   "same_side_other_vertices": next(len(side) - 1 for side in sides if u in side)})
    return polynomial, {"enumerated_side_subsets": 1 << len(left), "isolated_vertices": isolated,
                        "coverage_histogram": [{"selected_left": s, "free_right": f, "multiplicity": count}
                                               for (s, f), count in sorted(histogram.items())],
                        "core_independent_counts": core,
                        "universal_vertex_factorizations": factorizations}


def contiguous(values):
    output = []
    for k in values:
        if output and output[-1][-1] + 1 == k:
            output[-1].append(k)
        else:
            output.append([k])
    return [[row[0], row[-1]] for row in output]


def build():
    selected = json.loads((SOURCE / "trigger_graph_structure_bounded.json").read_text())["records"][:2]
    triple = {r["replica_counter"]: r for r in json.loads((SOURCE / "safe_triple_census.json").read_text())["checkpoints"]}
    records = []
    for source in selected:
        graph = json.loads((ROOT / source["graph_artifact"]).read_text())
        sides = source["structure"]["components_with_edges"][0]["sides"]
        if len(source["structure"]["components_with_edges"]) != 1:
            raise ValueError("this bounded scorer expects the saved one-component witnesses")
        counts, structure = independent_polynomial(graph, sides)
        a, d = len(graph["safe_sites"]), graph["N"] - graph["k0"]
        assert a == d == 173
        counts += [0] * (d + 1 - len(counts))
        survival = [Fraction(counts[k], comb(d, k)) for k in range(d + 1)]
        hazard = [None] + [(1 - survival[k] / survival[k - 1]) if survival[k - 1] else None
                          for k in range(1, d + 1)]
        # First coefficients reuse the existing exact counts as bookkeeping;
        # no homology, trigger, or triple census is rerun.
        exact3 = triple[graph["replica_counter"]]
        assert counts[:4] == [1, d, exact3["b2_safe_pairs"], exact3["safe_graph_triangles"]]
        record = {"selection": source["selection"], "source_graph": source["graph_artifact"],
                  "counter": graph["replica_counter"], "seed": graph["seed"],
                  "N": graph["N"], "h12": graph["h12"], "k0": graph["k0"], "candidate_sites": d,
                  "independent_counts": counts, **structure,
                  "survival": [frac(v) for v in survival],
                  "hazard": [frac(v) if v is not None else None for v in hazard],
                  "maximum_pair_safe_k": max(k for k, value in enumerate(counts) if value),
                  "mean_first_pair_trigger_step": frac(sum(survival[:-1])),
                  "first_trigger_quantiles": {str(q): next(k for k in range(d + 1) if survival[k] <= 1 - q)
                                              for q in (Fraction(1, 10), Fraction(1, 2), Fraction(9, 10))},
                  "true_birth_k3": {"safe_triples": exact3["actual_safe_triples"],
                                    "minimal_triple_nonfaces": exact3["minimal_nonfaces_size3"],
                                    "survival": frac(Fraction(exact3["actual_safe_triples"], comb(d, 3)))}}
        records.append(record)
    A, B = records
    # Both full polynomials contain (1+z)^139.  Remove that common isolate
    # factor; this yields a short positive-coefficient dominance certificate.
    common_isolates = min(A["isolated_vertices"], B["isolated_vertices"])
    residuals = []
    for row in (A, B):
        extra = row["isolated_vertices"] - common_isolates
        residuals.append(multiply(row["core_independent_counts"], [comb(extra, k) for k in range(extra + 1)]))
    width = max(map(len, residuals))
    for polynomial in residuals:
        polynomial += [0] * (width - len(polynomial))
    reduced_difference = [b - a for a, b in zip(*residuals)]
    while reduced_difference and reduced_difference[-1] == 0:
        reduced_difference.pop()
    first_nonzero = next(i for i, value in enumerate(reduced_difference) if value)
    positive_factor = reduced_difference[first_nonzero:]
    survival_difference = [B["independent_counts"][k] - A["independent_counts"][k] for k in range(d + 1)]
    hazard_integer_cross = {}
    for k in range(1, d + 1):
        if A["independent_counts"][k - 1] and B["independent_counts"][k - 1]:
            # Sign is the sign of h_B(k)-h_A(k); denominators are positive.
            hazard_integer_cross[k] = (A["independent_counts"][k] * B["independent_counts"][k - 1]
                                      - B["independent_counts"][k] * A["independent_counts"][k - 1])
    comparison = {
        "independent_count_difference_B_minus_A": survival_difference,
        "survival_B_above": contiguous([k for k, value in enumerate(survival_difference) if value > 0]),
        "survival_B_below": contiguous([k for k, value in enumerate(survival_difference) if value < 0]),
        "survival_equal": contiguous([k for k, value in enumerate(survival_difference) if value == 0]),
        "hazard_sign_integer_crossproducts": hazard_integer_cross,
        "hazard_B_above": contiguous([k for k, value in hazard_integer_cross.items() if value > 0]),
        "hazard_B_below": contiguous([k for k, value in hazard_integer_cross.items() if value < 0]),
        "hazard_equal": contiguous([k for k, value in hazard_integer_cross.items() if value == 0]),
        "dominance_factor": {"common_1_plus_z_power": common_isolates, "z_power": first_nonzero,
                             "remaining_coefficients": positive_factor,
                             "all_remaining_coefficients_strictly_positive": all(v > 0 for v in positive_factor)},
        "mean_first_trigger_difference_B_minus_A": frac(Fraction(B["mean_first_pair_trigger_step"]["exact"])
                                                        - Fraction(A["mean_first_pair_trigger_step"]["exact"]))}
    return {"parent_commit": "c1fbcc43beefacd889f215f6d05bc3c608e65652", "new_samples": 0,
            "records": records, "comparison": comparison,
            "claim_boundary": "Pair-only clock on two saved real trigger graphs, not the complete rank-two birth clock. Genuine triples and larger minimal triggers can absorb earlier; only true k<=3 is available here. No population or continuum inference."}


if __name__ == "__main__":
    result = build()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "pair_only_survival.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    for row in result["records"]:
        print(row["selection"], "max_k", row["maximum_pair_safe_k"], "mean", row["mean_first_pair_trigger_step"]["value"],
              "quantiles", row["first_trigger_quantiles"], "factors", row["universal_vertex_factorizations"])
    print({key: value for key, value in result["comparison"].items()
           if key not in ("independent_count_difference_B_minus_A", "hazard_sign_integer_crossproducts")})
