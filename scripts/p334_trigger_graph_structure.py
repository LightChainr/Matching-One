#!/usr/bin/env python3
"""Only bipartite and componentwise chain/Ferrers checks on archived graphs."""
from __future__ import annotations

import argparse
from collections import deque
import csv
from itertools import combinations
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT/"results/local-20260831/P334-cooperative-closure"


def canonical_cycle(cycle):
    return min(tuple(seq[i:]+seq[:i]) for seq in (cycle, list(reversed(cycle)))
               for i in range(len(seq)))


def shortest_odd_cycle(adj):
    best = None
    for root in adj:
        parent, distance, queue = {root: None}, {root: 0}, deque([root])
        while queue:
            v = queue.popleft()
            for w in sorted(adj[v]):
                if w not in distance:
                    distance[w], parent[w] = distance[v]+1, v
                    queue.append(w)
                elif distance[w] % 2 == distance[v] % 2:
                    left, x = [], v
                    while x is not None:
                        left.append(x); x = parent[x]
                    right, x = [], w
                    while x not in left:
                        right.append(x); x = parent[x]
                    cycle = canonical_cycle(left[:left.index(x)+1]+list(reversed(right)))
                    if best is None or (len(cycle), cycle) < (len(best), best):
                        best = cycle
    return list(best) if best is not None else None


def shortest_path(adj, start, target):
    queue, seen = deque([[start]]), {start}
    while queue:
        path = queue.popleft()
        if path[-1] == target: return path
        for v in sorted(adj[path[-1]]):
            if v not in seen:
                seen.add(v); queue.append(path+[v])
    raise ValueError("claimed component is disconnected")


def analyze(graph):
    adj = {v: set() for v in graph["safe_sites"]}
    for v, w in graph["minimal_trigger_pairs"]:
        adj[v].add(w); adj[w].add(v)
    components, seen = [], set()
    first_non_nested = None
    bipartite = True
    for root in sorted(adj):
        if root in seen: continue
        color, queue, valid = {root: 0}, deque([root]), True
        seen.add(root)
        while queue:
            v = queue.popleft()
            for w in sorted(adj[v]):
                if w not in color:
                    color[w] = 1-color[v]; seen.add(w); queue.append(w)
                elif color[w] == color[v]: valid = False
        bipartite &= valid
        sides = [sorted(v for v in color if color[v] == c) for c in (0, 1)]
        witness = None
        if valid:
            for side in sides:
                for a, b in combinations(side, 2):
                    left, right = adj[a]-adj[b], adj[b]-adj[a]
                    if left and right:
                        candidate = (a, b, min(left), min(right))
                        if witness is None or candidate < witness: witness = candidate
        if witness is not None:
            a, b, x, y = witness
            item = {"same_side_sites": [a, b], "exclusive_neighbors": [x, y],
                    "trigger_pairs": [[a, x], [b, y]],
                    "safe_cross_pairs": [[a, y], [b, x]],
                    "neighborhoods": {str(a): sorted(adj[a]), str(b): sorted(adj[b])},
                    "same_component_shortest_path": shortest_path(adj, a, b),
                    "component_vertices": sorted(color),
                    "minimal_witness_size": 4}
            if first_non_nested is None or witness < tuple(first_non_nested["same_side_sites"]+first_non_nested["exclusive_neighbors"]):
                first_non_nested = item
        components.append({"vertices": sorted(color),
                           "edges": sum(len(adj[v]) for v in color)//2,
                           "bipartite": valid, "sides": sides,
                           "side_degree_sequences": [sorted((len(adj[v]) for v in s), reverse=True) for s in sides],
                           "chain_ferrers": valid and witness is None})
    degrees = [len(adj[v]) for v in adj]
    nontrivial = [c for c in components if c["edges"]]
    return {
        "safe_vertices": len(adj), "trigger_edges": len(graph["minimal_trigger_pairs"]),
        "trigger_wedges": sum(comb(d, 2) for d in degrees),
        "trigger_triangles": sum(len(adj[v]&adj[w]) for v, w in graph["minimal_trigger_pairs"])//3,
        "isolated_vertices": sum(not adj[v] for v in adj),
        "nonisolated_components": len(nontrivial), "components_with_edges": nontrivial,
        "bipartite": bipartite,
        "every_connected_component_chain_ferrers": all(c["chain_ferrers"] for c in components),
        "shortest_odd_cycle": None if bipartite else shortest_odd_cycle(adj),
        "non_nested_neighborhood_witness": first_non_nested,
    }


def selected_rows(scope):
    saved = json.loads((RESULT/"scalar_state_collisions.json").read_text())
    witness = saved["environments"]["N425_second"]["witness"]
    rows = [("existing_N425_"+name, witness[name]["original_checkpoint_row"])
            for name in ("A", "B")]
    if scope == "bounded":
        for n in (325, 425):
            with (RESULT/f"raw/N{n}.geometry_pilot.csv").open(newline="") as stream:
                raw = list(csv.DictReader(stream))
            for o in ("first", "second"):
                eligible = [r for r in raw if r["orientation"] == o and
                            comb(int(r["checkpoint_b1_safe_count"]), 2) > int(r["checkpoint_b2_safe_pairs"])]
                for r in sorted(eligible, key=lambda r: int(r["replica"]))[:5]:
                    rows.append((f"first5_N{n}_{o}_{r['replica']}", r))
    return rows


def fraction_record(value):
    return {"numerator": value.numerator, "denominator": value.denominator,
            "text": str(value), "decimal": float(value)}


def saved_bipartite_baseline(records):
    output = {}
    expectations = []
    for record in records[:2]:
        structure = record["structure"]
        expected = Fraction(0)
        for component in structure["components_with_edges"]:
            L, R = map(len, component["sides"])
            m, capacity = component["edges"], L*R
            if m >= 2:
                expected += Fraction((L*comb(R, 2)+R*comb(L, 2))*m*(m-1),
                                     capacity*(capacity-1))
        observed = structure["trigger_wedges"]
        output[record["selection"]] = {
            "observed_wedges": observed,
            "componentwise_uniform_bipartite_expected": fraction_record(expected),
            "excess": fraction_record(Fraction(observed)-expected)}
        expectations.append(expected)
    difference = records[1]["structure"]["trigger_wedges"]-records[0]["structure"]["trigger_wedges"]
    expected_difference = expectations[1]-expectations[0]
    return {"rows": output, "observed_difference_B_minus_A": difference,
            "expected_difference_B_minus_A": fraction_record(expected_difference),
            "excess_difference_B_minus_A": fraction_record(Fraction(difference)-expected_difference),
            "expected_fraction_of_observed_difference": fraction_record(expected_difference/difference),
            "interpretation": "Posthoc structural arithmetic on the two selected graphs. Fix the observed component vertex sets, L/R capacities and m, but the uniform bipartite benchmark may be disconnected or contain isolates: it is not conditioned on preserving connectedness or the complete component decomposition. Capacity includes both component size and side imbalance. Not a new independent test, population percentage or causal attribution."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scope", choices=("existing", "bounded"), default="bounded")
    args = parser.parse_args()
    raw_dir = RESULT/"trigger_graph_raw"
    raw_dir.mkdir(exist_ok=True)
    records = []
    for label, row in selected_rows(args.scope):
        n, o, counter = int(row["n"]), row["orientation"], int(row["replica"])
        meta = json.loads((RESULT/f"raw/N{n}.metadata.json").read_text())
        h12 = meta["designs"][0][f"{o}_period_matrix"][0][1]
        path = raw_dir/f"N{n}-{o}-{counter}.json"
        if not path.exists():
            subprocess.run([str(args.binary), str(n), str(h12), str(row["k0"]),
                            str(meta["seed"]), str(counter), str(path)], check=True)
        graph = json.loads(path.read_text())
        structure = analyze(graph)
        a, m = structure["safe_vertices"], structure["trigger_edges"]
        trigger_degree_square_sum = 2*structure["trigger_wedges"]+2*m
        safe_degree_square_sum = a*(a-1)**2-4*(a-1)*m+trigger_degree_square_sum
        if (a != int(row["checkpoint_b1_safe_count"]) or
                comb(a, 2)-m != int(row["checkpoint_b2_safe_pairs"]) or
                safe_degree_square_sum != int(row["checkpoint_sum_child_b1_sq"])):
            raise ValueError("replayed trigger graph disagrees with archived counts")
        records.append({"selection": label, "source_row": row,
                        "graph_artifact": str(path.relative_to(ROOT)), "structure": structure})
        print(label, "bipartite", structure["bipartite"], "component_chain",
              structure["every_connected_component_chain_ferrers"],
              "edges/wedges/triangles", m, structure["trigger_wedges"], structure["trigger_triangles"],
              "odd_cycle", structure["shortest_odd_cycle"],
              "nonnested", structure["non_nested_neighborhood_witness"])
    payload = {
        "schema": "matching-one/p334-trigger-graph-two-class-check/v1", "new_samples": 0,
        "scope": args.scope,
        "selection_rule": "Two saved N425 witnesses, then (if bounded) the five lowest eligible counters in each of N325/N425 first/second. Eligible means an archived rank-one checkpoint with at least one minimal trigger pair; no graph-property screening.",
        "tested_classes_only": ["bipartite", "each connected component is chain/Ferrers"],
        "chain_definition": "Within each connected bipartite component, neighborhoods of every same-side vertex pair must be nested; disjoint 2K2 across different components is allowed.",
        "records": records,
        "summary": {"graphs": len(records),
                    "bipartite_pass": sum(r["structure"]["bipartite"] for r in records),
                    "component_chain_pass": sum(r["structure"]["every_connected_component_chain_ferrers"] for r in records)},
        "posthoc_two_saved_graph_bipartite_baseline": saved_bipartite_baseline(records),
        "claim_boundary": "A fixed small real-checkpoint census, not an all-configuration theorem, universality law, classifier search, or population prevalence estimate. All graphs reuse existing production checkpoints."
    }
    (RESULT/f"trigger_graph_structure_{args.scope}.json").write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")


if __name__ == "__main__":
    main()
