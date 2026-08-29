#!/usr/bin/env python3
"""Freeze and exactly normalize the N650 C2 x C5 mixed-join score."""

from __future__ import annotations

import argparse
from collections import Counter
from decimal import Decimal, getcontext
from fractions import Fraction
import json
from pathlib import Path
from typing import Iterable, Sequence


P_REF = Fraction(592746050790, 10**12)
FIBER_SIZE = 10
SOURCE_FIBERS = 65


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _decimal_text(value: Fraction, digits: int = 40) -> str:
    getcontext().prec = digits
    return str(Decimal(value.numerator) / Decimal(value.denominator))


def partition_rank(number_of_vertices: int, edges: Iterable[tuple[int, int]]) -> int:
    parent = list(range(number_of_vertices))

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for first, second in edges:
        root_first = find(first)
        root_second = find(second)
        if root_first != root_second:
            parent[max(root_first, root_second)] = min(root_first, root_second)
    return number_of_vertices - len({find(vertex) for vertex in range(number_of_vertices)})


def mixed_join_redundancy(
    number_of_vertices: int,
    base_edges: Sequence[tuple[int, int]],
    factor2_edges: Sequence[tuple[int, int]],
    factor5_edges: Sequence[tuple[int, int]],
) -> int:
    """The nonnegative partition-rank modular defect J_25."""

    base = list(base_edges)
    rank_base = partition_rank(number_of_vertices, base)
    rank_2 = partition_rank(number_of_vertices, base + list(factor2_edges))
    rank_5 = partition_rank(number_of_vertices, base + list(factor5_edges))
    rank_25 = partition_rank(
        number_of_vertices, base + list(factor2_edges) + list(factor5_edges)
    )
    value = rank_2 + rank_5 - rank_25 - rank_base
    if value < 0:
        raise AssertionError("partition-rank submodularity was violated")
    return value


def _selected_edges(selected: Sequence[bool]) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    if len(selected) != FIBER_SIZE:
        raise ValueError("a CRT fiber has ten cells")
    factor2 = []
    factor5 = []
    # Index a*5+b is the cell (a mod 2,b mod 5).
    for b in range(5):
        if selected[b] and selected[5 + b]:
            factor2.append((b, 5 + b))
    for a in range(2):
        row = [a * 5 + b for b in range(5) if selected[a * 5 + b]]
        if row:
            factor5.extend((row[0], vertex) for vertex in row[1:])
    return factor2, factor5


def isolated_fiber_join(selected: Sequence[bool]) -> int:
    """J_25 on one fiber with no pre-existing connectivity."""

    occupied = [index for index, keep in enumerate(selected) if keep]
    if not occupied:
        return 0
    relabel = {old: new for new, old in enumerate(occupied)}
    factor2, factor5 = _selected_edges(selected)
    restricted2 = [(relabel[a], relabel[b]) for a, b in factor2]
    restricted5 = [(relabel[a], relabel[b]) for a, b in factor5]
    return mixed_join_redundancy(len(occupied), [], restricted2, restricted5)


def incidence_cycle_rank(selected: Sequence[bool]) -> int:
    """First Betti number of the C2-row/C5-column incidence graph."""

    edges = [(a, b) for a in range(2) for b in range(5) if selected[a * 5 + b]]
    if not edges:
        return 0
    rows = {a for a, _ in edges}
    columns = {b for _, b in edges}
    adjacency = {("r", a): [] for a in rows} | {("c", b): [] for b in columns}
    for a, b in edges:
        adjacency[("r", a)].append(("c", b))
        adjacency[("c", b)].append(("r", a))
    components = 0
    seen: set[tuple[str, int]] = set()
    for start in adjacency:
        if start in seen:
            continue
        components += 1
        stack = [start]
        seen.add(start)
        while stack:
            vertex = stack.pop()
            for neighbor in adjacency[vertex]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
    return len(edges) - len(adjacency) + components


def local_color_record(bits: int) -> tuple[int, int, int, int, int]:
    black = [bool((bits >> index) & 1) for index in range(FIBER_SIZE)]
    white = [not value for value in black]
    join_black = isolated_fiber_join(black)
    join_white = isolated_fiber_join(white)
    if join_black != incidence_cycle_rank(black):
        raise AssertionError("black join/incidence identity failed")
    if join_white != incidence_cycle_rank(white):
        raise AssertionError("white join/incidence identity failed")
    common_black = sum(black[b] and black[5 + b] for b in range(5))
    common_white = sum(white[b] and white[5 + b] for b in range(5))
    if join_black != max(common_black - 1, 0):
        raise AssertionError("two-row black closed form failed")
    if join_white != max(common_white - 1, 0):
        raise AssertionError("two-row white closed form failed")
    return join_black, join_white, join_black - join_white, join_black + join_white, bits.bit_count()


def weighted_moments(p: Fraction, include_rationals: bool = False) -> dict:
    records = [local_color_record(bits) for bits in range(1 << FIBER_SIZE)]
    means = [Fraction() for _ in range(4)]
    seconds = [[Fraction() for _ in range(4)] for _ in range(4)]
    for record in records:
        values = record[:4]
        black_count = record[4]
        weight = p**black_count * (1 - p) ** (FIBER_SIZE - black_count)
        for row in range(4):
            means[row] += weight * values[row]
            for column in range(4):
                seconds[row][column] += weight * values[row] * values[column]
    covariance = [
        [seconds[row][column] - means[row] * means[column] for column in range(4)]
        for row in range(4)
    ]
    output = {
        "order": ["J_black", "J_white", "C_odd", "C_even"],
        "p": _decimal_text(p),
        "p_exact": _fraction_text(p),
        "mean_decimal": [_decimal_text(value) for value in means],
        "covariance_decimal": [
            [_decimal_text(value) for value in row] for row in covariance
        ],
    }
    if include_rationals:
        output["mean_exact"] = [_fraction_text(value) for value in means]
        output["covariance_exact"] = [
            [_fraction_text(value) for value in row] for row in covariance
        ]
    return output


def connected_two_fiber_witness() -> dict:
    """A four-vertex example where base connectivity creates one mixed cycle."""

    # Occupied global labels: fiber 0 cells (0,0),(0,1),(1,0), and fiber 1
    # cell (0,0).  Relabel them in this listed order.
    occupied = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0)]
    factor2 = [(0, 2)]
    factor5 = [(0, 1)]
    # Two ordinary connectivity edges from the two arms to the next fiber.
    base = [(1, 3), (2, 3)]
    full = mixed_join_redundancy(4, base, factor2, factor5)
    isolated = mixed_join_redundancy(4, [], factor2, factor5)
    if (isolated, full, full - isolated) != (0, 1, 1):
        raise AssertionError("connected mixed-interaction witness failed")
    return {
        "occupied_vertices_as_(fiber,a,b)": [list(value) for value in occupied],
        "base_connectivity_edges_relabelled": [list(value) for value in base],
        "factor2_edges_relabelled": [list(value) for value in factor2],
        "factor5_edges_relabelled": [list(value) for value in factor5],
        "J_isolated": isolated,
        "J_full": full,
        "R_connected": full - isolated,
    }


def render() -> dict:
    records = [local_color_record(bits) for bits in range(1 << FIBER_SIZE)]
    half = weighted_moments(Fraction(1, 2), include_rationals=True)
    p_ref = weighted_moments(P_REF)
    odd_distribution = Counter(record[2] for record in records)
    half_mean = half["mean_exact"]
    half_covariance = half["covariance_exact"]
    if half_mean[:2] != ["499/1024", "499/1024"]:
        raise AssertionError("unexpected p=1/2 join mean")
    if half_mean[2] != "0" or half_covariance[2][2] != "681/512":
        raise AssertionError("unexpected p=1/2 odd normalization")

    # Algebraic null and alternative in a four-corner response table.
    h0, u2, u5 = 7, -3, 11
    additive = {"00": h0, "10": h0 + u2, "01": h0 + u5, "11": h0 + u2 + u5}
    additive_defect = additive["11"] - additive["10"] - additive["01"] + additive["00"]
    interaction_strength = 4
    interacting = dict(additive)
    interacting["11"] += interaction_strength
    interaction_defect = (
        interacting["11"] - interacting["10"] - interacting["01"] + interacting["00"]
    )
    if additive_defect != 0 or interaction_defect != interaction_strength:
        raise AssertionError("four-corner mechanism gate failed")

    return {
        "schema": "matching-one.p200-n650-mixed-join-phaseB.v1",
        "issue": 200,
        "status": "frozen_phaseB_observable_before_N650_acquisition",
        "factors": {"K2": "C2", "K5": "C5", "fiber": "C2 x C5"},
        "primary_observable": {
            "per_color_full": "J_c=r(Pi_c join R2)+r(Pi_c join R5)-r(Pi_c join R2 join R5)-r(Pi_c)",
            "partition_rank": "r(Pi)=number_of_colored_vertices-number_of_blocks",
            "per_color_local": "J_c_local=sum_over_65_fibers b1(C2-C5 occupied-cell incidence graph)",
            "two_row_closed_form": "b1=max(k_c-1,0), k_c=#C5 columns with both C2 lifts color c",
            "connected_residual": "R_c=J_c_full-J_c_local",
            "color_channels": {"even": "R_black+R_white", "odd": "R_black-R_white"},
            "orientation_channels": {
                "ES": "(even_(23,11)+even_(17,19))/2",
                "ED": "(even_(23,11)-even_(17,19))/2",
                "OS": "(odd_(23,11)+odd_(17,19))/2",
                "OD": "(odd_(23,11)-odd_(17,19))/2",
            },
            "primary_state_order": ["ES", "ED", "OS", "OD"],
        },
        "null_and_alternative": {
            "exact_algebra": {
                "factor_additive_table": additive,
                "mixed_defect": additive_defect,
                "statement": "any h_ab=h_00+a*u2+b*u5 has zero four-corner mixed defect",
            },
            "synthetic_interaction": {
                "table": interacting,
                "inserted_lambda": interaction_strength,
                "recovered_mixed_defect": interaction_defect,
            },
            "radial_clock_boundary": "q2 and Jordan alone constrain only unmarked endpoint states; zero is the prediction of their minimal factor-additive topology extension, not a theorem of either scalar clock",
            "primary_null": "mean vector (ES,ED,OS,OD)=0 after the configurationwise isolated-fiber subtraction",
            "nonzero_mean_identifies": "a nonlocal mixed C2-by-C5 connectivity interaction",
            "nonzero_does_not_alone_identify": "continuum RG memory, Jordan structure, or path noncommutativity",
        },
        "acquisition": {
            "probability": _decimal_text(P_REF),
            "endpoint_orientations": [[23, 11], [17, 19]],
            "source_fibers_per_orientation": SOURCE_FIBERS,
            "canonical_sampling": "for each counter permutation draw K~Binomial(650,p_ref), use the prefix configuration, and use the same K/permutation labels for both orientations",
            "same_configuration_rule": "evaluate Pi, Pi join R2, Pi join R5, and Pi join R2 join R5 on the same configuration; never treat four corners or two paths as replicas",
            "full_lift_rule": "retain NN-black and NN+NNN-white displacement-potential connectivity; primary partition rank ignores potentials, while ambient-H1 ranks for all four corners are archived as a typed secondary vector",
            "batches": "at least 100 synchronized counter batches",
        },
        "score": {
            "primary": "joint GLS chi-square of the zero 4-vector (ES,ED,OS,OD)",
            "covariance": "full 4x4 delete-one covariance from synchronized same-flow batches",
            "degrees_of_freedom": "numerical rank of the frozen covariance pseudoinverse",
            "report_also": [
                "four marginal signed z values",
                "raw J_black/J_white and exact local baselines",
                "secondary four-corner ambient-H1 mixed defects without folding them into the primary score",
            ],
            "decision": {
                "pass": "no detected connected mixed-factor interaction at frozen precision",
                "fail": "factor-additive topology bridge is false; classify the nonzero ES/ED/OS/OD channel before proposing an RG-memory mechanism",
            },
        },
        "toy_exact_normalization": {
            "all_binary_fiber_configurations": 1 << FIBER_SIZE,
            "isolated_residual_R_c_zero_configurationwise": True,
            "p_half": half,
            "p_ref": p_ref,
            "p_half_C_odd_distribution_counts": {
                str(key): odd_distribution[key] for key in sorted(odd_distribution)
            },
            "recommended_dimensionless_reference": "divide an orientation's odd connected residual by sqrt(65*Var_p_ref[J_black-J_white]); empirical delete-one covariance remains primary",
        },
        "connected_toy_witness": connected_two_fiber_witness(),
        "claim_layers": {
            "exact": [
                "partition-rank J is nonnegative",
                "isolated J equals the C2-C5 incidence first Betti number",
                "the closed form max(k-1,0) and all toy moments",
                "factor-additive four-corner tables have zero mixed defect",
            ],
            "mechanism_null": "the minimal factor-additive topology extension of either radial clock has zero connected residual mean",
            "exploratory": "a persistent OD component would be a chiral mixed-factor interaction worth comparing with the P57 odd direction",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = render()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
