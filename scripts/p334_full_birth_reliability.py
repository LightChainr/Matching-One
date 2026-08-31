#!/usr/bin/env python3
"""Exact full physical birth clocks on the two saved two-port networks.

Tree-decomposition connectivity states keep both terminals in every bag.
Selected-site weights are counted once, when a site is forgotten. Integer
Kronecker encoding multiplies the exact generating polynomials without FFT.
"""
from collections import defaultdict
from fractions import Fraction
import itertools
import json
from math import comb
from pathlib import Path
import time

import networkx as nx

from p334_pair_only_survival import contiguous, frac
from p334_pair_triple_survival import strip_one_plus_z, binomial, multiply

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/p334-contracted-full-clock"


def canonical(labels):
    mapping = {}
    return tuple(-1 if a < 0 else mapping.setdefault(a, len(mapping)) for a in labels)


def safety_polynomial(network, all_sites):
    graph = nx.Graph()
    graph.add_nodes_from(network["vertices"])
    graph.add_edges_from(network["edges"])
    terminals = tuple(network["terminals"])
    random = set(network["vacant_sites"])
    fixed = set(network["fixed_components"]) | set(terminals)
    width, tree = nx.approximation.treewidth_min_degree(graph)
    root = min(tree.nodes, key=lambda bag: (-len(bag), tuple(sorted(bag))))
    bits, base = len(all_sites) + 2, 1 << (len(all_sites) + 2)
    started = time.monotonic()
    stats = {"treewidth_upper_bound": width, "bags": len(tree), "maximum_states": 0,
             "join_pairs": 0, "bag_state_counts": []}

    def guard(table):
        stats["maximum_states"] = max(stats["maximum_states"], len(table))
        if len(table) > 200_000 or time.monotonic() - started > 240:
            raise RuntimeError("bounded two-network reliability budget exceeded")

    def initial(bag):
        positions = [i for i, v in enumerate(bag) if v in random]
        result = {}
        for flags in itertools.product((False, True), repeat=len(positions)):
            selected = {bag[i] for i, yes in zip(positions, flags) if yes} | fixed
            state = canonical([i if v in selected else -1 for i, v in enumerate(bag)])
            result[state] = 1
        return result

    def forget(table, bag, vertex):
        index = bag.index(vertex)
        result = defaultdict(int)
        remaining = tuple(v for v in bag if v != vertex)
        neighbor_indices = [i for i, v in enumerate(bag) if graph.has_edge(vertex, v)]
        for state, weight in table.items():
            labels = list(state)
            if labels[index] >= 0:
                merge = {labels[i] for i in neighbor_indices if labels[i] >= 0} | {labels[index]}
                target = min(merge)
                labels = [target if label in merge else label for label in labels]
                if labels[bag.index(terminals[0])] == labels[bag.index(terminals[1])]:
                    continue
                if vertex in random:
                    weight <<= bits
            del labels[index]
            result[canonical(labels)] += weight
        return dict(result), remaining

    def project(table, bag, target):
        for vertex in sorted(set(bag) - set(target)):
            table, bag = forget(table, bag, vertex)
        for vertex in sorted(set(target) - set(bag)):
            expanded = tuple(sorted((*bag, vertex)))
            position = expanded.index(vertex)
            result = {}
            for state, weight in table.items():
                if vertex in random:
                    labels = list(state)
                    labels.insert(position, -1)
                    result[canonical(labels)] = weight
                labels = list(state)
                labels.insert(position, max(state, default=-1) + 1)
                result[canonical(labels)] = weight
            table, bag = result, expanded
        guard(table)
        return table

    def occupancy(state):
        return tuple(i for i, value in enumerate(state) if value >= 0)

    def join(left, right):
        groups = defaultdict(list)
        for state, weight in right.items():
            groups[occupancy(state)].append((state, weight))
        result = defaultdict(int)
        for state, weight in left.items():
            active = occupancy(state)
            for other, other_weight in groups[active]:
                stats["join_pairs"] += 1
                labels = list(state)
                # Add the second child's connectivity equivalences.
                for value in set(other) - {-1}:
                    merge = {labels[i] for i in active if other[i] == value}
                    target = min(merge)
                    labels = [target if label in merge else label for label in labels]
                # Terminals sort first in every bag and must remain separate.
                if labels[0] == labels[1]:
                    continue
                result[canonical(labels)] += weight * other_weight
        guard(result)
        return dict(result)

    def visit(node, parent=None):
        bag = tuple(sorted(set(node) | set(terminals)))
        table = None
        for child in sorted(tree.neighbors(node), key=lambda b: tuple(sorted(b))):
            if child == parent:
                continue
            child_table, child_bag = visit(child, node)
            message = project(child_table, child_bag, bag)
            table = message if table is None else join(table, message)
        if table is None:
            table = initial(bag)
        stats["bag_state_counts"].append({"bag": list(bag), "states": len(table)})
        guard(table)
        return table, bag

    table, bag = visit(root)
    result = project(table, bag, tuple(sorted(terminals)))
    packed = sum(result.values())
    irrelevant = len(set(all_sites) - random)
    packed *= pow(base + 1, irrelevant)
    mask = base - 1
    coefficients = [(packed >> (bits * k)) & mask for k in range(len(all_sites) + 1)]
    stats.update({"elapsed_seconds": time.monotonic() - started,
                  "random_network_vertices": len(random), "irrelevant_random_sites": irrelevant,
                  "coefficient_bit_stride": bits})
    return coefficients, stats


def main():
    source = json.loads((OUTPUT / "whole_event_networks.json").read_text())["records"]
    quartic = {r["counter"]: r for r in json.loads((ROOT / "results/p334-quartic-clock/quartic_survival.json").read_text())["records"]}
    records = []
    for mapped in source:
        networks = [c["two_terminal_network"] for c in mapped["port_components"] if "two_terminal_network" in c]
        if len(networks) != 1 or any(len(c["addresses"]) > 2 for c in mapped["port_components"]):
            raise ValueError("bounded scorer requires the observed unique two-port component")
        counts, stats = safety_polynomial(networks[0], mapped["vacant_sites"])
        old = quartic[mapped["counter"]]
        if counts[:5] != old["independent_counts"][:5]:
            raise ValueError("full physical polynomial disagrees with exact known k<=4")
        d = len(mapped["vacant_sites"])
        survival = [Fraction(counts[k], comb(d, k)) for k in range(d + 1)]
        hazard = [None] + [1 - survival[k] / survival[k-1] if survival[k-1] else None for k in range(1, d+1)]
        mean = sum(survival[:-1])
        row = {"counter": mapped["counter"], "seed": mapped["seed"], "N": mapped["N"], "k0": mapped["k0"],
               "true_safe_counts": counts, "true_survival": [frac(s) for s in survival],
               "true_hazard": [frac(h) if h is not None else None for h in hazard],
               "mean_true_birth_step": frac(mean), "maximum_true_safe_k": max(k for k, value in enumerate(counts) if value),
               "remaining_mean_shortening_after_quartics": frac(Fraction(old["mean_first_trigger_step"]["exact"]) - mean),
               "minimal_quintic_count_without_subset_enumeration": old["independent_counts"][5] - counts[5],
               "true_first_trigger_quantiles": {str(q): next(k for k in range(d+1) if survival[k] <= 1-q)
                                                for q in (Fraction(1,10), Fraction(1,2), Fraction(9,10))},
               **stats}
        records.append(row)
        print(mapped["counter"], "mean", row["mean_true_birth_step"], "states", stats["maximum_states"],
              "joins", stats["join_pairs"], "seconds", stats["elapsed_seconds"], flush=True)
    A, B = records
    difference = [b-a for a,b in zip(A["true_safe_counts"], B["true_safe_counts"])]
    crosses = {k: A["true_safe_counts"][k]*B["true_safe_counts"][k-1] - B["true_safe_counts"][k]*A["true_safe_counts"][k-1]
               for k in range(1,d+1) if A["true_safe_counts"][k-1] and B["true_safe_counts"][k-1]}
    factors = [strip_one_plus_z(row["true_safe_counts"]) for row in records]
    common = min(power for power, _ in factors)
    residuals = [multiply(p, binomial(power - common)) for power, p in factors]
    width = max(map(len, residuals))
    for p in residuals:
        p += [0] * (width - len(p))
    residual_difference = [b - a for a,b in zip(*residuals)]
    z_power = next(k for k, value in enumerate(residual_difference) if value)
    comparison = {"true_safe_count_B_minus_A": difference,
                  "survival_B_above": contiguous([k for k,v in enumerate(difference) if v>0]),
                  "survival_B_below": contiguous([k for k,v in enumerate(difference) if v<0]),
                  "survival_equal": contiguous([k for k,v in enumerate(difference) if v==0]),
                  "hazard_crossproducts_B_minus_A": crosses,
                  "hazard_B_above": contiguous([k for k,v in crosses.items() if v>0]),
                  "hazard_B_below": contiguous([k for k,v in crosses.items() if v<0]),
                  "hazard_equal": contiguous([k for k,v in crosses.items() if v==0]),
                  "dominance_factor": {"common_1_plus_z_power": common, "z_power": z_power,
                                       "remaining_coefficients": residual_difference[z_power:],
                                       "all_remaining_coefficients_strictly_positive": all(v > 0 for v in residual_difference[z_power:])},
                  "mean_true_birth_B_minus_A": frac(Fraction(B["mean_true_birth_step"]["exact"]) - Fraction(A["mean_true_birth_step"]["exact"]))}
    result = {"parent_commit": "1614a17e10997656fdf2d5520846fff2a228a5cd", "new_samples": 0,
              "records": records, "comparison": comparison,
              "scope": "Complete physical rank-two first-birth distribution for uniform fresh insertion order conditional on each of two fixed real N425 checkpoints. No truncation, no population law or continuum claim."}
    (OUTPUT / "full_physical_birth_clock.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print({k:v for k,v in comparison.items() if k not in ("true_safe_count_B_minus_A", "hazard_crossproducts_B_minus_A", "dominance_factor")})


if __name__ == "__main__":
    main()
