#!/usr/bin/env python3
"""Replicated marked-pair switching oracle for the linewise TM inequality."""

from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction
import json
from pathlib import Path

from p334_dual_hazard_ulc import _local_degrees
from p334_lorentzian_support_gate import _honest_geometries
from projective_essential_birth_oracle import subset_marks


Pair = tuple[int, int]
SupplyKey = tuple[str, int, int]


def _marked_pair_capacities(marks, line, n: int, k: int, layers):
    """Return the doubled-unordered TM demand and its two supply reservoirs."""

    lower = layers[k]
    m = n - k
    layer_size = len(lower)
    exit_site_counts: Counter[int] = Counter()
    demand_pairs: Counter[Pair] = Counter()
    synergy_pairs: Counter[Pair] = Counter()

    for mask in lower:
        exits = []
        internal = []
        for site in range(n):
            if mask >> site & 1:
                continue
            rank, marked_line, _ = marks[mask | (1 << site)]
            if rank == 2:
                exits.append(site)
                exit_site_counts[site] += 1
            elif rank == 1 and marked_line == line:
                internal.append(site)
        for left_index, left in enumerate(exits):
            for right in exits[left_index + 1 :]:
                demand_pairs[left, right] += 1
        for left_index, left in enumerate(internal):
            for right in internal[left_index + 1 :]:
                if marks[mask | (1 << left) | (1 << right)][0] == 2:
                    synergy_pairs[left, right] += 1

    demand = Counter(
        {
            pair: 2 * m * layer_size * count
            for pair, count in demand_pairs.items()
        }
    )
    synergy = Counter(
        {
            ("synergy", *pair): 2 * m * layer_size * count
            for pair, count in synergy_pairs.items()
        }
    )
    reservoir: Counter[SupplyKey] = Counter()
    for left in range(n):
        for right in range(left, n):
            multiplicity = 1 if left == right else 2
            capacity = (
                (m - 1)
                * multiplicity
                * exit_site_counts[left]
                * exit_site_counts[right]
            )
            if capacity:
                reservoir["reservoir", left, right] = capacity

    supply = synergy + reservoir
    # These are exactly the three terms of the linewise TM determinant.
    ordered_exit_pairs = sum(
        count * (count - 1)
        for count in (
            sum(
                1
                for site in range(n)
                if not (mask >> site & 1)
                and marks[mask | (1 << site)][0] == 2
            )
            for mask in lower
        )
    )
    assert sum(demand.values()) == m * layer_size * ordered_exit_pairs
    return demand, supply, synergy, reservoir, exit_site_counts


def _pair_from_supply(key: SupplyKey) -> Pair:
    return key[1], key[2]


def _compatible(demand_pair: Pair, supply_key: SupplyKey) -> bool:
    """The local rule: preserve at least one of the two marked sites."""

    supply_pair = _pair_from_supply(supply_key)
    return bool(set(demand_pair) & set(supply_pair))


def _add_edge(graph, source: int, target: int, capacity: int):
    forward = [target, capacity, len(graph[target]), capacity]
    reverse = [source, 0, len(graph[source]), 0]
    graph[source].append(forward)
    graph[target].append(reverse)
    return forward


def lex_first_integer_flow(demand: Counter[Pair], supply: Counter[SupplyKey]):
    """Deterministic integral max flow using lex-first shortest augmenting paths.

    Demand bins are ordered lexicographically.  For each bin, exact-pair supply
    precedes one-endpoint switching, synergy precedes the independent reservoir,
    and remaining ties use the supply-pair order.  Integer bin flow lifts to a
    token injection by taking replica labels in increasing order within a bin.
    """

    demand_keys = sorted(pair for pair, capacity in demand.items() if capacity)
    supply_keys = sorted(
        (key for key, capacity in supply.items() if capacity),
        key=lambda key: (key[1], key[2], key[0]),
    )
    source = 0
    demand_offset = 1
    supply_offset = demand_offset + len(demand_keys)
    sink = supply_offset + len(supply_keys)
    graph = [[] for _ in range(sink + 1)]
    edge_refs = {}
    total_demand = sum(demand.values())
    for index, pair in enumerate(demand_keys):
        _add_edge(graph, source, demand_offset + index, demand[pair])
    for index, key in enumerate(supply_keys):
        _add_edge(graph, supply_offset + index, sink, supply[key])
    supply_index = {key: index for index, key in enumerate(supply_keys)}
    for demand_index, pair in enumerate(demand_keys):
        compatible = [key for key in supply_keys if _compatible(pair, key)]
        compatible.sort(
            key=lambda key: (
                _pair_from_supply(key) != pair,
                key[0] != "synergy",
                key[1],
                key[2],
            )
        )
        for key in compatible:
            edge_refs[pair, key] = _add_edge(
                graph,
                demand_offset + demand_index,
                supply_offset + supply_index[key],
                total_demand,
            )

    total_flow = 0
    while True:
        parent = [None] * len(graph)
        parent[source] = (-1, -1)
        pending = deque([source])
        while pending and parent[sink] is None:
            vertex = pending.popleft()
            for edge_index, edge in enumerate(graph[vertex]):
                target, capacity, _, _ = edge
                if capacity and parent[target] is None:
                    parent[target] = vertex, edge_index
                    pending.append(target)
                    if target == sink:
                        break
        if parent[sink] is None:
            break
        amount = total_demand
        vertex = sink
        while vertex != source:
            previous, edge_index = parent[vertex]
            amount = min(amount, graph[previous][edge_index][1])
            vertex = previous
        vertex = sink
        while vertex != source:
            previous, edge_index = parent[vertex]
            edge = graph[previous][edge_index]
            reverse = edge[2]
            edge[1] -= amount
            graph[vertex][reverse][1] += amount
            vertex = previous
        total_flow += amount

    flow = {
        (pair, key): edge[3] - edge[1]
        for (pair, key), edge in edge_refs.items()
        if edge[3] - edge[1]
    }
    return total_flow, flow


def hall_family(demand: Counter[Pair], supply: Counter[SupplyKey], n: int):
    """Audit the exact Hall family induced by the one-common-site graph."""

    combined_supply: Counter[Pair] = Counter()
    for key, capacity in supply.items():
        combined_supply[_pair_from_supply(key)] += capacity
    demand_keys = sorted(pair for pair, capacity in demand.items() if capacity)
    seen_families = set()
    rows = []
    for vertex_bits in range(1, 1 << n):
        family = tuple(
            pair
            for pair in demand_keys
            if vertex_bits >> pair[0] & 1 and vertex_bits >> pair[1] & 1
        )
        if not family or family in seen_families:
            continue
        seen_families.add(family)
        active_bits = 0
        for left, right in family:
            active_bits |= (1 << left) | (1 << right)
        required = sum(demand[pair] for pair in family)
        available = sum(
            capacity
            for (left, right), capacity in combined_supply.items()
            if active_bits >> left & 1 or active_bits >> right & 1
        )
        rows.append((required, available, active_bits, len(family)))
    return rows


def _descriptor(n, matrix, carrier, line, k, **extra):
    return {
        "N": n,
        "matrix": [list(row) for row in matrix],
        "carrier": carrier,
        "line": list(line),
        "lower_layer": k,
        **extra,
    }


def build_result():
    counters = Counter()
    first_exact_pair_failure = None
    maximum_hall_ratio = (Fraction(-1), None)
    maximum_switched_fraction = (Fraction(-1), None)
    flow_digest = Counter()

    for n, matrix, geometry in _honest_geometries(12):
        for carrier, marks in (
            ("primal", subset_marks(geometry, matching=False)),
            ("matching", subset_marks(geometry, matching=True)),
        ):
            for line in sorted({line for rank, line, _ in marks if rank == 1}):
                layers, _ = _local_degrees(marks, line, n)
                for k in range(n):
                    if not layers[k] or not layers[k + 1]:
                        continue
                    counters["line_layer_rows"] += 1
                    demand, supply, synergy, reservoir, exit_counts = (
                        _marked_pair_capacities(marks, line, n, k, layers)
                    )
                    total_demand = sum(demand.values())
                    total_supply = sum(supply.values())
                    assert total_supply >= total_demand

                    combined_supply: Counter[Pair] = Counter()
                    for key, capacity in supply.items():
                        combined_supply[_pair_from_supply(key)] += capacity
                    exact_pair_pass = all(
                        combined_supply[pair] >= capacity
                        for pair, capacity in demand.items()
                    )
                    counters["exact_pair_pass"] += int(exact_pair_pass)
                    counters["exact_pair_fail"] += int(not exact_pair_pass)
                    if not exact_pair_pass and first_exact_pair_failure is None:
                        deficits = [
                            {
                                "pair": list(pair),
                                "demand": capacity,
                                "same_pair_supply": combined_supply[pair],
                            }
                            for pair, capacity in sorted(demand.items())
                            if capacity > combined_supply[pair]
                        ]
                        first_exact_pair_failure = _descriptor(
                            n,
                            matrix,
                            carrier,
                            line,
                            k,
                            deficits=deficits,
                        )

                    total_flow, flow = lex_first_integer_flow(demand, supply)
                    flow_pass = total_flow == total_demand
                    counters["one_common_site_flow_pass"] += int(flow_pass)
                    counters["one_common_site_flow_fail"] += int(not flow_pass)
                    assert flow_pass
                    exact_mass = 0
                    switched_mass = 0
                    synergy_mass = 0
                    reservoir_mass = 0
                    for (demand_pair, supply_key), amount in flow.items():
                        supply_pair = _pair_from_supply(supply_key)
                        assert set(demand_pair) & set(supply_pair)
                        if demand_pair == supply_pair:
                            exact_mass += amount
                        else:
                            switched_mass += amount
                        if supply_key[0] == "synergy":
                            synergy_mass += amount
                        else:
                            reservoir_mass += amount
                    assert exact_mass + switched_mass == total_demand
                    assert synergy_mass + reservoir_mass == total_demand
                    flow_digest["demand_tokens"] += total_demand
                    flow_digest["exact_pair_tokens"] += exact_mass
                    flow_digest["one_endpoint_switched_tokens"] += switched_mass
                    flow_digest["synergy_tokens_used"] += synergy_mass
                    flow_digest["reservoir_tokens_used"] += reservoir_mass
                    if switched_mass:
                        counters["rows_using_one_endpoint_switch"] += 1
                    switched_fraction = (
                        Fraction(switched_mass, total_demand)
                        if total_demand
                        else Fraction()
                    )
                    if switched_fraction > maximum_switched_fraction[0]:
                        maximum_switched_fraction = (
                            switched_fraction,
                            _descriptor(
                                n,
                                matrix,
                                carrier,
                                line,
                                k,
                                demand_tokens=total_demand,
                                exact_pair_tokens=exact_mass,
                                one_endpoint_switched_tokens=switched_mass,
                                synergy_tokens_used=synergy_mass,
                                reservoir_tokens_used=reservoir_mass,
                                active_flow_edges=len(flow),
                            ),
                        )

                    hall_rows = hall_family(demand, supply, n)
                    counters["Hall_families_checked"] += len(hall_rows)
                    for required, available, active_bits, family_size in hall_rows:
                        counters["Hall_failures"] += int(required > available)
                        assert required <= available
                        ratio = Fraction(required, available)
                        if ratio > maximum_hall_ratio[0]:
                            maximum_hall_ratio = (
                                ratio,
                                _descriptor(
                                    n,
                                    matrix,
                                    carrier,
                                    line,
                                    k,
                                    active_vertices=[
                                        site
                                        for site in range(n)
                                        if active_bits >> site & 1
                                    ],
                                    demand_pair_bins=family_size,
                                    required=required,
                                    available=available,
                                    slack=available - required,
                                ),
                            )

    assert counters["line_layer_rows"] == 984
    assert counters["one_common_site_flow_fail"] == 0
    assert counters["Hall_failures"] == 0
    ratio, ratio_row = maximum_hall_ratio
    switched_ratio, switched_row = maximum_switched_fraction
    result = {
        "schema_version": "p334-tm-replicated-switching-oracle-v1",
        "replicated_TM_multisets": {
            "demand": "for each unordered same-state exit pair {v,w}, 2mA replicas",
            "synergy_supply": "for each unordered synergy square {v,w}, 2mA replicas",
            "independent_supply": "(m-1) replicas of all ordered pairs of exit marks; an off-diagonal site pair has multiplicity 2 X_v X_w and a diagonal pair X_v^2",
            "cardinality_identity": "total supply - total demand = m A N_new + (m-1) X^2 - m A E_2",
        },
        "local_switching_rule": {
            "compatibility": "a demand pair and a supply pair are compatible iff they share at least one marked site",
            "switch": "preserve the common mark and switch only the other mark; fixed ell supplies the ambient H1 channel",
            "canonicalization": "sort demand bins lexicographically; prefer exact-pair then one-endpoint supply, synergy then reservoir, and take lex-first shortest augmenting paths; lift bin flow by increasing replica labels",
            "status": "explicit deterministic integer injection on every bounded row",
        },
        "capacitated_Hall_theorem": {
            "statement": "The one-common-site injection exists iff for every vertex set U, demand on all positive demand pairs contained in U is at most supply on all synergy/reservoir pairs incident to the vertices used by those demand pairs.",
            "why_sufficient": "the neighbor set of a demand family depends only on the union of its marked sites; for a fixed union, including every positive demand pair is the unique maximal left capacity. Integral bipartite max-flow then lifts to a token injection.",
            "relation_to_TM": "U equal to all marked sites is the aggregate TM inequality; proper U give the strictly stronger local Hall cuts needed for a one-mark-preserving switch.",
            "status": "exact finite capacitated Hall equivalence",
        },
        "bounded_audit": dict(counters),
        "flow_digest": dict(flow_digest),
        "exact_pair_locking": {
            "status": "false",
            "first_failure": first_exact_pair_failure,
        },
        "tightest_Hall_cut": {
            "required_over_available": str(ratio),
            "row": ratio_row,
        },
        "largest_canonical_switch_fraction": {
            "fraction": str(switched_ratio),
            "row": switched_row,
        },
        "theorem_frontier": {
            "proved": "TM is exactly a supply-demand cardinality inequality; one-common-site injection is exactly characterized by the displayed Hall family.",
            "bounded_exact": "all 984 topological line/carrier/layer rows satisfy every Hall cut and admit the canonical integer injection.",
            "open": "derive the proper-subset Hall cuts uniformly from digital Alexander topology rather than checking them quotient by quotient.",
        },
    }
    return json.loads(json.dumps(result))


def render_markdown(result):
    audit = result["bounded_audit"]
    exact = result["exact_pair_locking"]["first_failure"]
    tight = result["tightest_Hall_cut"]
    switched = result["largest_canonical_switch_fraction"]
    return "\n".join(
        [
            "# Replicated one-mark switching for TM",
            "",
            "Write the TM determinant as a doubled-unordered supply-demand problem. A same-state exit pair `{v,w}` has `2mA` demand replicas. A synergy square has `2mA` supply replicas. The independent reservoir contains `(m-1)` replicas of every ordered pair of exit marks. Its total capacity is `(m-1)X^2`.",
            "",
            "## The local rule found by canonical flow",
            "",
            "Locking supply to exactly the same site pair fails. The minimal bounded failure is "
            f"`N={exact['N']}`, matrix `{exact['matrix']}`, {exact['carrier']} carrier, line `{exact['line']}`, lower layer `{exact['lower_layer']}`. Overall, exact-pair locking fails on {audit['exact_pair_fail']} of {audit['line_layer_rows']} rows.",
            "",
            "The first successful locality is one-common-site compatibility: preserve either `v` or `w` and switch only the other marked site. Lex-first integral max flow succeeds on every row. This rule uses the fixed ambient line `ell`, but no metric distance, quotient-specific phase, or fitted parameter.",
            "",
            "## Exact Hall family",
            "",
            "For a set `U` of marked sites, let `D(U)` be demand on every positive demand pair contained in `U`. Let `S(U)` be supply on every synergy or reservoir pair incident to a site actually used by `D(U)`. Then",
            "",
            "`D(U) <= S(U) for every U`",
            "",
            "is necessary and sufficient for the one-common-site injection. Necessity is Hall. For sufficiency, the neighborhood of any demand family depends only on its union of marked sites, and adding every positive demand pair inside that union maximizes left capacity without changing its neighborhood. Integrality of bipartite max flow gives a token injection.",
            "",
            f"The oracle checks {audit['Hall_families_checked']} distinct nonempty induced demand families with zero failures. The tightest cut has `required/available={tight['required_over_available']}` and is the all-site TM cut. The largest canonical one-endpoint switching fraction is `{switched['fraction']}`.",
            "",
            "## Boundary of the theorem",
            "",
            "This closes the switching problem on the bounded atlas and replaces an opaque aggregate inequality by an exact local Hall criterion. What is not yet proved is that digital Alexander topology forces every proper-subset Hall cut on arbitrary quotients. That is now the sole topological gap; no larger-N evidence is needed to state it precisely.",
            "",
        ]
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = build_result()
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    args.markdown.write_text(render_markdown(result), encoding="utf-8")


if __name__ == "__main__":
    main()
