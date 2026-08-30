#!/usr/bin/env python3
"""Classify the matching defect against named surface polynomials.

The exact positive result is a site-state rank-image quotient.  The exact
negative result is local: one ordinary ribbon edge, or one transition
polynomial smoothing, cannot realize the four-terminal junction required by
one occupied site of the typed incidence spine.  Derivatives can reweight
states but cannot repair a missing local partition state.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from itertools import combinations
from math import comb
from pathlib import Path
from typing import Iterable

from integer_period_torus import classify_configuration, integer_torus_geometry
from p144_typed_incidence_spine import build_oracle as build_spine_oracle


Matrix = tuple[tuple[int, int], tuple[int, int]]
Partition = tuple[tuple[int, ...], ...]


def canonical_partition(blocks: Iterable[Iterable[int]]) -> Partition:
    normalized = []
    for block in blocks:
        value = tuple(sorted(block))
        if value:
            normalized.append(value)
    return tuple(sorted(normalized))


def endpoint_partition(edge: tuple[int, int] | None) -> Partition:
    if edge is None:
        return canonical_partition(((0,), (1,), (2,), (3,)))
    first, second = edge
    return canonical_partition(
        ((first, second),) + tuple((port,) for port in range(4) if port not in edge)
    )


def connected_on_four_ports(edges: Iterable[tuple[int, int]]) -> bool:
    adjacency = {port: set() for port in range(4)}
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    seen = {0}
    frontier = [0]
    while frontier:
        current = frontier.pop()
        for neighbor in adjacency[current] - seen:
            seen.add(neighbor)
            frontier.append(neighbor)
    return len(seen) == 4


def local_partition_nogo() -> dict[str, object]:
    ports = tuple(range(4))
    target_j4 = canonical_partition((ports,))
    all_edges = tuple(combinations(ports, 2))
    single_edge_states = [endpoint_partition(None)] + [
        endpoint_partition(edge) for edge in all_edges
    ]
    transition_pairings = [
        canonical_partition(pairing)
        for pairing in (
            ((0, 1), (2, 3)),
            ((0, 3), (1, 2)),
            ((0, 2), (1, 3)),
        )
    ]
    connected_edge_sets = [
        edges
        for size in range(len(all_edges) + 1)
        for edges in combinations(all_edges, size)
        if connected_on_four_ports(edges)
    ]
    minimum_edges = min(map(len, connected_edge_sets))
    minimum_witnesses = [
        [list(edge) for edge in edges]
        for edges in connected_edge_sets if len(edges) == minimum_edges
    ]
    return {
        "terminal_ports": ["north", "east", "south", "west"],
        "required_partition": [list(block) for block in target_j4],
        "ordinary_single_ribbon_edge_partitions": [
            [list(block) for block in partition] for partition in single_edge_states
        ],
        "ordinary_single_edge_contains_J4": target_j4 in single_edge_states,
        "topological_transition_pairings": [
            [list(block) for block in partition] for partition in transition_pairings
        ],
        "transition_pairings_contain_J4": target_j4 in transition_pairings,
        "minimum_independent_edges_needed_for_J4": minimum_edges,
        "minimum_J4_edge_gadgets": minimum_witnesses,
        "independent_subsets_of_minimum_gadget": 2**minimum_edges,
        "unwanted_mixed_subsets_if_only_all_off_and_all_on_are_kept": 2**minimum_edges - 2,
        "minimal_extension": (
            "promote one site to a typed four-terminal partition element; its two "
            "states are J_edge on edge-midpoint ports and J_face on face-center ports"
        ),
        "derivative_boundary": (
            "a specialization or derivative changes weights of existing states and "
            "cannot create the absent J4 terminal partition"
        ),
    }


def genus_flags(rank_black: int, rank_white: int) -> tuple[int, int]:
    """Rank-determined torus-side genus flags.

    On an honest quotient these are the actual genera of the complementary
    carriers.  On a degenerate quotient they are the finite-cover/rank-image
    extension and do not assert that a chosen quotient CW presentation is a
    Krushkal spanning ribbon subgraph.
    """
    if (rank_black, rank_white) not in {(0, 2), (1, 1), (2, 0)}:
        raise AssertionError("rank pair escaped unrestricted digital Alexander duality")
    return int(rank_black == 2), int(rank_white == 2)


def site_rank_image_quotient() -> dict[str, object]:
    exact = build_spine_oracle()["honest_torus_exact_oracle"]
    N = int(exact["N"])
    rows: Counter[tuple[int, int, int, int, int, int]] = Counter()
    from_rank = [0] * (N + 1)
    from_relative = [0] * (N + 1)
    from_genus = [0] * (N + 1)
    unrefined = [0] * (N + 1)
    rank_pairs = set()
    genus_pairs = set()

    for source in exact["coefficient_rows"]:
        occupied = int(source["occupied"])
        rank_black = int(source["rank_black"])
        rank_white = int(source["rank_white"])
        count = int(source["count"])
        if rank_black + rank_white != 2:
            raise AssertionError("unrestricted rank sum failed")
        q = rank_black - 1
        if q != (rank_black - rank_white) // 2:
            raise AssertionError("one-sided and relative charges disagree")
        genus_black, genus_white = genus_flags(rank_black, rank_white)
        if q != genus_black - genus_white:
            raise AssertionError("Krushkal genus-side charge disagrees")
        rows[(occupied, rank_black, rank_white, genus_black, genus_white, q)] += count
        from_rank[occupied] += (rank_black - 1) * count
        from_relative[occupied] += q * count
        from_genus[occupied] += (genus_black - genus_white) * count
        unrefined[occupied] += count
        rank_pairs.add((rank_black, rank_white))
        genus_pairs.add((genus_black, genus_white))

    direct = list(exact["matching_Bernstein_coefficients_direct"])
    expected_unrefined = [comb(N, occupied) for occupied in range(N + 1)]
    if not from_rank == from_relative == from_genus == direct:
        raise AssertionError("rank/genus quotient did not reproduce matching coefficients")
    if unrefined != expected_unrefined:
        raise AssertionError("occupation-only quotient is not binomial")

    return {
        "geometry": exact["geometry"],
        "N": N,
        "configurations": int(exact["configurations"]),
        "coefficient_rows": [
            {
                "occupied": key[0],
                "rank_black": key[1],
                "rank_white": key[2],
                "genus_black": key[3],
                "genus_white": key[4],
                "q": key[5],
                "count": count,
            }
            for key, count in sorted(rows.items())
        ],
        "rank_pairs": [list(pair) for pair in sorted(rank_pairs)],
        "carrier_genus_pairs": [list(pair) for pair in sorted(genus_pairs)],
        "matching_Bernstein_coefficients": direct,
        "matching_from_black_rank_derivative": from_rank,
        "matching_from_relative_rank_derivative": from_relative,
        "matching_from_Krushkal_genus_derivative": from_genus,
        "occupation_only_Bernstein_counts": unrefined,
        "exact_quotients": {
            "black_rank": (
                "R_site(a,b;z)=sum_B a^|B| b^(N-|B|) z^r_B; "
                "M=(z*d_z-1)R_site|z=1"
            ),
            "relative_rank": (
                "Z_rel(a,b;Q)=sum_B a^|B| b^(N-|B|) Q^q; "
                "M=Q*d_Q Z_rel|Q=1"
            ),
            "site_Krushkal": (
                "K_site(a,b;X,Y)=sum_B a^|B| b^(N-|B|) "
                "X^g_B Y^g_W; M=(X*d_X-Y*d_Y)K_site|X=Y=1"
            ),
        },
        "independent_output_sources": {
            "homogeneous_form": ["a_for_B", "b_for_W", "Q_relative_topology"],
            "fixed_N_reduced_form": ["t=a/b", "Q_relative_topology"],
            "extra_topology_sources_beyond_occupation": 1,
            "equivalent_topology_source": "z^r_B with the affine derivative z*d_z-1",
        },
        "all_derivative_routes_agree": True,
    }


def hnf_matrices_of_order(order: int) -> list[Matrix]:
    rows = []
    for first in range(1, order + 1):
        if order % first:
            continue
        second = order // first
        for shear in range(first):
            rows.append(((first, shear), (0, second)))
    return rows


def rank_defect_bernstein(matrix: Matrix) -> list[int]:
    geometry = integer_torus_geometry(matrix, name="p144-minimal-source-witness")
    coefficients = [0] * (geometry.n + 1)
    for mask in range(1 << geometry.n):
        active = tuple(bool(mask & (1 << vertex)) for vertex in range(geometry.n))
        rank_black = classify_configuration(geometry, active)[0].max_rank
        rank_white = classify_configuration(
            geometry, tuple(not value for value in active), matching=True
        )[0].max_rank
        if rank_black + rank_white != 2:
            raise AssertionError("unrestricted rank sum failed in minimal witness")
        coefficients[sum(active)] += rank_black - 1
    return coefficients


def minimal_topology_source_witness(maximum_order: int = 4) -> dict[str, object]:
    first = None
    audits = []
    for order in range(1, maximum_order + 1):
        classes: dict[tuple[int, ...], list[Matrix]] = {}
        for matrix in hnf_matrices_of_order(order):
            coefficients = tuple(rank_defect_bernstein(matrix))
            classes.setdefault(coefficients, []).append(matrix)
        row = {
            "order": order,
            "HNF_representatives": sum(len(group) for group in classes.values()),
            "distinct_matching_polynomials": len(classes),
            "matching_classes": [
                {
                    "bernstein": list(coefficients),
                    "matrices": [[list(part) for part in matrix] for matrix in matrices],
                }
                for coefficients, matrices in sorted(classes.items())
            ],
        }
        audits.append(row)
        if first is None and len(classes) > 1:
            groups = list(classes.items())
            first = {
                "order": order,
                "matrix_A": [list(part) for part in groups[0][1][0]],
                "matching_A": list(groups[0][0]),
                "matrix_B": [list(part) for part in groups[1][1][0]],
                "matching_B": list(groups[1][0]),
                "shared_occupation_only_counts": [comb(order, k) for k in range(order + 1)],
            }
    if first is None:
        raise AssertionError("no topology-source witness found")
    return {
        "search_orders": [1, maximum_order],
        "first_distinguishing_order": first["order"],
        "witness": first,
        "audits": audits,
        "conclusion": (
            "the occupation-only site state sum is identical for every geometry of "
            "fixed order, but the matching defect already differs at order two; at "
            "least one topology source is necessary"
        ),
    }


def named_family_classification() -> list[dict[str, object]]:
    common = (
        "fails the one-ground-element-per-site local-state contract: an edge has "
        "two ends, whereas a selected site requires a typed four-terminal J4"
    )
    return [
        {
            "family": "multivariate Tutte / random-cluster",
            "state_index": "subsets of edges of a fixed graph",
            "topological_readout": "component/matroid rank only",
            "direct_site_specialization": False,
            "reason": common,
        },
        {
            "family": "Bollobas-Riordan",
            "state_index": "spanning ribbon subgraphs indexed by edge subsets",
            "topological_readout": "rank, nullity and ribbon boundary/genus data",
            "direct_site_specialization": False,
            "reason": common,
        },
        {
            "family": "Las Vergnas",
            "state_index": "edge subsets in a cycle-to-bond matroid perspective",
            "topological_readout": "matroid-perspective ranks",
            "direct_site_specialization": False,
            "reason": common,
        },
        {
            "family": "Krushkal",
            "state_index": "spanning embedded subgraphs indexed by edge subsets",
            "topological_readout": (
                "subgraph and complement genus sources are sufficient for q on a torus"
            ),
            "direct_site_specialization": False,
            "reason": (
                "closest named readout, but " + common + "; K_site is therefore a "
                "site-state Krushkal analogue/rank-image quotient, not K_G itself"
            ),
        },
        {
            "family": "topological transition / 2025 embedded-graph vertex polynomial",
            "state_index": "ribbon-edge smoothings, partial dual/twist edge states",
            "topological_readout": "boundary-component or twisted-dual vertex counts",
            "direct_site_specialization": False,
            "reason": (
                "its local four-valent states are pairwise smoothings; none is the J4 "
                "partition, and a boundary-count derivative cannot create that state"
            ),
        },
    ]


def build_artifact() -> dict[str, object]:
    local = local_partition_nogo()
    quotient = site_rank_image_quotient()
    source_witness = minimal_topology_source_witness()
    classifications = named_family_classification()
    all_pass = (
        not local["ordinary_single_edge_contains_J4"]
        and not local["transition_pairings_contain_J4"]
        and local["minimum_independent_edges_needed_for_J4"] == 3
        and quotient["all_derivative_routes_agree"]
        and quotient["independent_output_sources"]["extra_topology_sources_beyond_occupation"] == 1
        and source_witness["first_distinguishing_order"] == 2
        and all(not row["direct_site_specialization"] for row in classifications)
    )
    return {
        "schema": "matching-one.p144-surface-polynomial-classification.v1",
        "issues": [144, 269],
        "pull_requests": [229, 267],
        "status": "named_edge_subset_nogo_and_exact_site_rank_image_quotient",
        "decision": {
            "named_polynomial": (
                "no direct specialization or derivative under the declared natural "
                "one-ground-element-per-site, product-local state-sum contract"
            ),
            "closest_named_readout": (
                "Krushkal subgraph/complement genus variables; the readout is exact "
                "after replacing edge subsets by typed site junction states"
            ),
            "exact_positive_object": (
                "R_site or equivalently Z_rel/K_site, the rank-image quotient of the "
                "typed incidence-spine partition state sum"
            ),
            "minimal_extension": (
                "one typed J4 site element locally and one relative topology source Q "
                "at terminal output, in addition to occupation"
            ),
        },
        "local_partition_no_go": local,
        "named_family_classification": classifications,
        "exact_site_rank_image_quotient": quotient,
        "minimal_topology_source_witness": source_witness,
        "claim_contract": {
            "excluded": (
                "arbitrary graph gadgets with coefficient extraction, cancellations or "
                "a graph reverse-engineered separately from each final scalar polynomial"
            ),
            "allowed": (
                "a fixed embedded object, one independent ground element per site, "
                "product local weights, specialization and source derivatives"
            ),
            "gadget_escape": (
                "three or more ordinary edges can realize J4 only if their intermediate "
                "subsets are removed by a block projector; that projector is precisely "
                "the additional typed site-state structure"
            ),
        },
        "dependencies": {
            "matching_defect_PR229_head": "83a14f3",
            "typed_incidence_spine": "e6953ba",
            "relative_source_collapse": "011d332",
            "unrestricted_degenerate_quotient_theorem": "73d4960",
        },
        "all_machine_gates_pass": all_pass,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--json",
        type=Path,
        default=Path("results/p144-surface-polynomial-classification/latest.json"),
    )
    args = parser.parse_args()
    artifact = build_artifact()
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(artifact["decision"], indent=2, sort_keys=True))
    if not artifact["all_machine_gates_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
