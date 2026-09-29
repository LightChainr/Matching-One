#!/usr/bin/env python3
"""Exact L4 successor geometry, not a new stochastic production block.

Reuse the existing lifted topology routine. Enumerate 2^16 subsets per
lattice, count prefix entry histories by integer DP (not 16! permutations),
and resolve the two-insertion law at matched direction/completion count.
"""
from collections import Counter, defaultdict
from fractions import Fraction
import importlib.util
from itertools import combinations
import json
import math
from pathlib import Path
import platform
import time


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "birth-completion-geometry-20260929/analyze.py"
SPEC = importlib.util.spec_from_file_location("existing_topology", SOURCE)
TOPO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(TOPO)
# Only topology() is called; its original L3 acquisition/DP is not run.
L, N = 4, 16
TOPO.L, TOPO.N = L, N
FACT = [math.factorial(k) for k in range(N+1)]


def occupied(mask):
    return [v for v in range(N) if mask >> v & 1]


def solve(lattice):
    started = time.perf_counter()
    steps = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    if lattice == "triangular":
        steps += [(1, 1), (-1, -1)]
    neighbors = [[(((v % L+dx) % L)+L*((v//L+dy) % L), dx, dy)
                  for dx, dy in steps] for v in range(N)]
    topo = [TOPO.topology(mask, neighbors) for mask in range(1 << N)]
    rank_counts = Counter((mask.bit_count(), r) for mask, (r, _) in enumerate(topo))
    # This is new L4 geometry, not a rerun of the old L3 complete table.
    assert topo[0] == (0, None) and topo[-1] == (2, None)
    assert topo[(1 << L)-1] == (1, (1, 0))  # one complete horizontal row
    completion = {}
    for mask, (rank, direction) in enumerate(topo):
        if rank == 1:
            completion[mask] = tuple(v for v in range(N) if not mask >> v & 1
                                     and topo[mask | (1 << v)][0] == 2)
    ways = {}
    cells = {}
    geometry = {}
    checked = 0
    for mask, complete_sites in completion.items():
        k, m = mask.bit_count(), N-mask.bit_count()
        direction = topo[mask][1]
        birth = Counter()
        for v in occupied(mask):
            old = mask ^ (1 << v)
            if topo[old][0] == 0:
                birth[k] += FACT[k-1]
            else:
                assert topo[old] == (1, direction)
                birth.update(ways[old])
        assert sum(birth.values()) == FACT[k]
        ways[mask] = birth
        safe = [v for v in range(N) if not mask >> v & 1 and v not in complete_sites]
        edges = [(v, w) for v, w in combinations(safe, 2)
                 if topo[mask | (1 << v) | (1 << w)][0] == 2]
        degrees = Counter(v for edge in edges for v in edge)
        c, e = len(complete_sites), len(edges)
        next_counts = Counter()
        for v in safe:
            new = mask | (1 << v)
            assert topo[new] == (1, direction)
            assert len(completion[new]) == c+degrees[v]
            next_counts[len(completion[new])] += 1
        # Direct unordered-pair outcome versus graph formula, one arithmetic
        # check inside the new enumeration, not a second topology census.
        vacancies = [v for v in range(N) if not mask >> v & 1]
        surviving_ordered_pairs = 2*sum(topo[mask | (1 << v) | (1 << w)][0] == 1
                                      for v, w in combinations(vacancies, 2))
        assert surviving_ordered_pairs == (m-c)*(m-c-1)-2*e
        checked += 1
        signature = tuple(sorted(next_counts.items()))
        key = k, direction, c
        if key not in cells:
            cells[key] = {"masks": 0, "edge_counts": Counter(), "signatures": {},
                          "birth": defaultdict(lambda: {"prefix_mass": 0, "edges_sum": 0,
                                                        "survival2_numerator_sum": 0,
                                                        "successor_nu_counts": Counter()})}
        cell = cells[key]
        cell["masks"] += 1
        cell["edge_counts"][e] += 1
        cell["signatures"].setdefault(signature, mask)
        for j, count in birth.items():
            entry = cell["birth"][j]
            entry["prefix_mass"] += count
            entry["edges_sum"] += count*e
            entry["survival2_numerator_sum"] += count*surviving_ordered_pairs
            for nu, number in next_counts.items():
                entry["successor_nu_counts"][nu] += count*number
        geometry[mask] = {"mask": mask, "occupied": occupied(mask), "k": k,
                          "direction": list(direction), "nu2": c,
                          "completion_sites": list(complete_sites), "safe_sites": safe,
                          "synergy_edges": edges, "edge_count": e,
                          "successor_nu_counts": dict(sorted(next_counts.items())),
                          "birth_prefix_counts": dict(sorted(birth.items())),
                          "survival2": str(Fraction(surviving_ordered_pairs, m*(m-1))) if m >= 2 else None}
    strong_witness = None
    weak_witness = None
    summaries = []
    variable_edge_cells = variable_transition_cells = 0
    for (k, direction, c), cell in sorted(cells.items()):
        m = N-k
        edge_variation = len(cell["edge_counts"]) > 1
        transition_variation = len(cell["signatures"]) > 1
        variable_edge_cells += edge_variation
        variable_transition_cells += transition_variation
        if edge_variation and strong_witness is None:
            candidate_masks = [mask for mask in cell["signatures"].values()]
            first = min(candidate_masks)
            second = min(mask for mask in candidate_masks
                         if geometry[mask]["edge_count"] != geometry[first]["edge_count"])
            strong_witness = {"selection": "first lexicographic cell with distinct edge counts",
                              "preparations": [geometry[first], geometry[second]]}
        birth_rows = []
        for j, stats in sorted(cell["birth"].items()):
            birth_rows.append({"j1": j, **stats,
                               "successor_nu_counts": dict(sorted(stats["successor_nu_counts"].items()))})
        summaries.append({"k": k, "direction": list(direction), "nu2": c,
                          "configurations": cell["masks"], "edge_count_histogram": dict(sorted(cell["edge_counts"].items())),
                          "distinct_successor_signatures": len(cell["signatures"]),
                          "birth_sufficient_statistics": birth_rows})
        if weak_witness is not None or m < 2:
            continue
        for a in range(1, k):
            early = [r for r in birth_rows if r["j1"] <= a]
            late = [r for r in birth_rows if r["j1"] > a]
            if not early or not late:
                continue
            def group(rows):
                mass = sum(r["prefix_mass"] for r in rows)
                edges_sum = sum(r["edges_sum"] for r in rows)
                num = sum(r["survival2_numerator_sum"] for r in rows)
                next_counts = Counter()
                for r in rows:
                    next_counts.update(r["successor_nu_counts"])
                return {"prefix_mass": mass, "edges_sum": edges_sum,
                        "mean_edges": str(Fraction(edges_sum, mass)),
                        "survival2": str(Fraction(num, mass*m*(m-1))),
                        "survivor_successor_law": {nu: str(Fraction(count, mass*(m-c)))
                                                   for nu, count in sorted(next_counts.items())} if m > c else {},
                        "mean_nu_next_given_survival": str(Fraction(sum(nu*count for nu, count in next_counts.items()), mass*(m-c))) if m > c else None}
            eg, lg = group(early), group(late)
            delta_e = Fraction(eg["mean_edges"])-Fraction(lg["mean_edges"])
            if not delta_e:
                continue
            delta_z = Fraction(eg["survival2"])-Fraction(lg["survival2"])
            assert delta_z == -2*delta_e/(m*(m-1))
            weak_witness = {"selection": "first lexicographic (k,direction,nu2,a) with nonzero exact mean-edge difference; structural existence, no p-value",
                            "k": k, "a": a, "direction": list(direction), "nu2": c,
                            "early": eg, "late": lg,
                            "mean_edge_difference": str(delta_e), "survival2_difference": str(delta_z),
                            "one_step_survival_both": str(Fraction(m-c, m))}
            break
    return {"lattice": lattice, "L": L, "N": N, "occupied_subsets": 1 << N,
            "rank_counts": [[k, rank, count] for (k, rank), count in sorted(rank_counts.items())],
            "rank_one_configurations": len(completion), "coarse_cells": len(cells),
            "cells_with_edge_variation": variable_edge_cells,
            "cells_with_successor_variation": variable_transition_cells,
            "pair_formula_rows": checked,
            "strong_preparation_witness": strong_witness,
            "uniform_start_birth_history_witness": weak_witness,
            "exact_cell_sufficient_statistics": summaries,
            "wall_seconds": time.perf_counter()-started}


def main():
    result = {"schema": "matching-one.completion-pair-synergy.v1", "date": "2026-09-29",
              "source": "uniform site permutation; square NN or triangular +diagonal occupied graph on L4 torus",
              "clock": "insertion count, not iid-label time", "python": platform.python_version(),
              "topology_source": "analysis/birth-completion-geometry-20260929/analyze.py:topology",
              "topology_source_commit": "8b37c20e5f2b1946cefbed98b09280a0cad25ed2",
              "scope": "Exact finite L4 computation, no Monte Carlo, 16! enumeration, cloud or continuum inference",
              "models": [solve(lattice) for lattice in ("square", "triangular")]}
    (HERE / "result.json").write_text(json.dumps(result, indent=2)+"\n")
    lines = ["# Completion-pair synergy: exact L4 successor calculation", "",
             "Uniform-permutation count clock. Exhaustive occupied subsets and integer prefix DP, not 16! paths.", ""]
    for model in result["models"]:
        lines += [f"## {model['lattice']}", "",
                  f"{model['rank_one_configurations']} rank-one configurations; {model['coarse_cells']} (k,D,nu2) cells.",
                  f"{model['cells_with_edge_variation']} cells have different pair counts; {model['cells_with_successor_variation']} different successor laws.", ""]
        witness = model["uniform_start_birth_history_witness"]
        if witness:
            lines += [f"First weak-history witness: k={witness['k']}, a={witness['a']}, D={witness['direction']}, nu={witness['nu2']}.",
                      f"Identical one-step survival: {witness['one_step_survival_both']}.",
                      f"Two-step early-minus-late survival: {witness['survival2_difference']}.",
                      f"Mean synergy-edge early-minus-late: {witness['mean_edge_difference']}.", ""]
        else:
            lines += ["No birth-cutoff mean-pair witness found; not a full weak-Markov proof.", ""]
        lines += [f"Wall time: {model['wall_seconds']:.3f}s.", ""]
    lines += ["Complete per-cell integer birth masses, edge sums and successor counts are in result.json.",
              "Existence witnesses are selected lexicographically from an exact finite census, not significance-ranked samples.", ""]
    (HERE / "RESULT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
