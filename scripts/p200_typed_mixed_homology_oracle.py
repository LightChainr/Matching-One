#!/usr/bin/env python3
"""Tiny exact typed-join oracle for Issue #200 Phase C."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gaussian_crt_commutator import ALPHA, BETA, multiplication_matrix  # noqa: E402
from integer_period_torus import (  # noqa: E402
    IntegerHomologyUnionFind,
    IntegerPeriods,
    _rank,
    integer_torus_geometry,
)
from p200_n650_hnf_maps import column_hnf  # noqa: E402


Matrix = tuple[tuple[int, int], tuple[int, int]]
PRODUCT = (3, 1)
FINAL_MATRIX = multiplication_matrix(PRODUCT)
GEOMETRY = integer_torus_geometry(FINAL_MATRIX, name="p200-N1-to-N10")
PERIOD_2 = IntegerPeriods(multiplication_matrix(ALPHA))
PERIOD_5 = IntegerPeriods(multiplication_matrix(BETA))
RELATIONS = {
    "0": (),
    "2": (PERIOD_5,),  # degree-two kernel: final -> norm-five quotient
    "5": (PERIOD_2,),  # degree-five kernel: final -> norm-two quotient
    "25": (PERIOD_5, PERIOD_2),
}


def fraction_record(value: Fraction) -> dict[str, str | float]:
    return {"exact": str(value), "decimal": float(value)}


def mask_from_integer(mask: int) -> tuple[bool, ...]:
    return tuple(bool(mask & (1 << index)) for index in range(10))


def connectivity_stats(active: tuple[bool, ...], edges, relations) -> tuple[int, int]:
    union_find = IntegerHomologyUnionFind(10, GEOMETRY.periods)
    for edge in edges:
        if active[edge.i] and active[edge.j]:
            union_find.add_edge(edge.i, edge.j, edge.dx, edge.dy)
    for quotient in relations:
        fibers: dict[tuple[int, int], list[int]] = {}
        for vertex, (enabled, coordinate) in enumerate(zip(active, GEOMETRY.coordinates)):
            if enabled:
                fibers.setdefault(quotient.quotient_key(coordinate), []).append(vertex)
        for fiber in fibers.values():
            anchor = fiber[0]
            anchor_x, anchor_y = GEOMETRY.coordinates[anchor]
            for vertex in fiber[1:]:
                x, y = GEOMETRY.coordinates[vertex]
                union_find.add_edge(anchor, vertex, x - anchor_x, y - anchor_y)
    roots = {union_find.find(i)[0] for i, enabled in enumerate(active) if enabled}
    generators = [vector for root in roots for vector in union_find.generators[root]]
    return sum(active) - len(roots), _rank(generators)


def local_incidence_cycle_rank(active: tuple[bool, ...]) -> int:
    """b1 of the occupied-edge incidence graph C2--C5 in one source fiber."""

    parent: dict[tuple[str, tuple[int, int]], tuple[str, tuple[int, int]]] = {}

    def find(node):
        parent.setdefault(node, node)
        if parent[node] != node:
            parent[node] = find(parent[node])
        return parent[node]

    def union(left, right):
        left, right = find(left), find(right)
        if left != right:
            parent[max(left, right)] = min(left, right)

    edges = 0
    for vertex, enabled in enumerate(active):
        if not enabled:
            continue
        coordinate = GEOMETRY.coordinates[vertex]
        node_2 = ("N2", PERIOD_2.quotient_key(coordinate))
        node_5 = ("N5", PERIOD_5.quotient_key(coordinate))
        union(node_2, node_5)
        edges += 1
    components = len({find(node) for node in parent})
    return edges - len(parent) + components


def color_rows(active: tuple[bool, ...], edges) -> dict:
    partition_ranks, ambient_ranks = {}, {}
    order_reverse_equal = True
    for name, relations in RELATIONS.items():
        partition_ranks[name], ambient_ranks[name] = connectivity_stats(active, edges, relations)
        if name == "25":
            reverse = connectivity_stats(active, edges, tuple(reversed(relations)))
            order_reverse_equal &= reverse == (partition_ranks[name], ambient_ranks[name])
    join_full = partition_ranks["2"] + partition_ranks["5"] - partition_ranks["25"] - partition_ranks["0"]
    join_local = local_incidence_cycle_rank(active)
    ambient_delta = ambient_ranks["25"] - ambient_ranks["2"] - ambient_ranks["5"] + ambient_ranks["0"]
    return {
        "partition_ranks": partition_ranks,
        "ambient_H1_ranks": ambient_ranks,
        "J_full": join_full,
        "J_local": join_local,
        "R_nonlocal": join_full - join_local,
        "Delta25_ambient_rank": ambient_delta,
        "join_order_equal": order_reverse_equal,
    }


def configuration_record(mask: int) -> dict:
    black = mask_from_integer(mask)
    white = tuple(not value for value in black)
    black_rows = color_rows(black, GEOMETRY.primal_edges)
    white_rows = color_rows(white, GEOMETRY.matching_edges)
    delta_black = black_rows["Delta25_ambient_rank"]
    delta_white = white_rows["Delta25_ambient_rank"]
    return {
        "mask": mask,
        "occupied_sites": [index for index, value in enumerate(black) if value],
        "black_NN": black_rows,
        "white_matching": white_rows,
        "ambient_color_even": Fraction(delta_black + delta_white, 2),
        "ambient_color_odd": Fraction(delta_black - delta_white, 2),
        "residual_color_even": Fraction(black_rows["R_nonlocal"] + white_rows["R_nonlocal"], 2),
        "residual_color_odd": Fraction(black_rows["R_nonlocal"] - white_rows["R_nonlocal"], 2),
    }


def pair_histogram(records: list[dict], field: str) -> dict[str, int]:
    counts = Counter(
        (record["black_NN"][field], record["white_matching"][field])
        for record in records
    )
    return {f"{left},{right}": count for (left, right), count in sorted(counts.items())}


def render() -> dict:
    records = [configuration_record(mask) for mask in range(1 << 10)]
    if not all(
        record["black_NN"]["join_order_equal"] and record["white_matching"]["join_order_equal"]
        for record in records
    ):
        raise AssertionError("typed joins failed path symmetry")
    # The matching involution swaps the two typed layers, rather than merely
    # complementing a mask while leaving the graph types fixed.
    if not all(
        Fraction(
            record["white_matching"]["Delta25_ambient_rank"]
            - record["black_NN"]["Delta25_ambient_rank"], 2
        ) == -record["ambient_color_odd"]
        and Fraction(
            record["white_matching"]["R_nonlocal"]
            - record["black_NN"]["R_nonlocal"], 2
        ) == -record["residual_color_odd"]
        for record in records
    ):
        raise AssertionError("typed matching involution did not negate the odd row")
    ambient_witness = next(
        record for record in records
        if len(record["occupied_sites"]) == 5 and record["ambient_color_odd"] != 0
    )
    residual_witness = next(
        record for record in records
        if len(record["occupied_sites"]) == 5 and record["residual_color_odd"] != 0
    )

    local_black = [record["black_NN"]["J_local"] for record in records]
    local_white = [record["white_matching"]["J_local"] for record in records]
    local_mean = sum(map(Fraction, local_black)) / 1024
    local_differences = [left - right for left, right in zip(local_black, local_white)]
    local_difference_mean = sum(map(Fraction, local_differences)) / 1024
    local_difference_variance = sum(
        (Fraction(value) - local_difference_mean) ** 2 for value in local_differences
    ) / 1024
    if local_mean != Fraction(499, 1024) or local_difference_variance != Fraction(681, 512):
        raise AssertionError("local incidence moments changed")

    final_hnf, _ = column_hnf(FINAL_MATRIX)
    hnf_2, _ = column_hnf(PERIOD_2.matrix)
    hnf_5, _ = column_hnf(PERIOD_5.matrix)
    return {
        "schema": "matching-one.p200-typed-mixed-homology-oracle.v1",
        "issue": 200,
        "status": "exact_tiny_configuration_oracle",
        "geometry": {
            "lineage": "N1 --(1+i)/(2-i)--> N2/N5 --other factor--> N10",
            "final_gaussian": [3, 1],
            "final_period_matrix": [list(row) for row in FINAL_MATRIX],
            "final_column_HNF": [list(row) for row in final_hnf],
            "N2_column_HNF": [list(row) for row in hnf_2],
            "N5_column_HNF": [list(row) for row in hnf_5],
            "coordinates_in_reference_engine_order": [list(value) for value in GEOMETRY.coordinates],
        },
        "definitions": {
            "typed_layers": "black uses NN connectivity; white uses NN+NNN matching connectivity on the same final lift",
            "relations": "R2 identifies same-colour sites in degree-two kernel fibers; R5 analogously for degree five; edges carry exact coordinate displacement",
            "ambient_rank": "rational rank of the span of all component winding generators in final-period coordinates",
            "ambient_Delta25_per_colour": "r25-r2-r5+r0",
            "ambient_colour_rows": {"even": "(Delta_B+Delta_W)/2", "odd": "(Delta_B-Delta_W)/2"},
            "partition_J_full": "rank2+rank5-rank25-rank0",
            "partition_J_local": "b1 of occupied C2--C5 incidence graph within the source fiber",
            "partition_R_nonlocal": "J_full-J_local",
        },
        "exhaustive_checks": {
            "configurations": 1024,
            "join_path_symmetry_every_configuration_and_colour": True,
            "typed_matching_layer_swap_negates_odd_rows": True,
            "ordinary_mask_complement_warning": "do not leave NN/matching graph types fixed; the exact involution swaps the typed layers",
            "ambient_delta_pair_histogram_B_W": pair_histogram(records, "Delta25_ambient_rank"),
            "nonlocal_residual_pair_histogram_B_W": pair_histogram(records, "R_nonlocal"),
        },
        "exact_local_normalization": {
            "E_J_local_black_at_p_half": fraction_record(local_mean),
            "E_J_local_white_at_p_half": fraction_record(sum(map(Fraction, local_white)) / 1024),
            "E_JB_minus_JW": fraction_record(local_difference_mean),
            "Var_JB_minus_JW": fraction_record(local_difference_variance),
        },
        "witnesses": {
            "nonzero_ambient_colour_odd": {
                key: (str(value) if isinstance(value, Fraction) else value)
                for key, value in ambient_witness.items()
            },
            "nonzero_partition_residual_colour_odd": {
                key: (str(value) if isinstance(value, Fraction) else value)
                for key, value in residual_witness.items()
            },
        },
        "production_interface": {
            "per_orientation_primary_rows": ["R_even=(R_black+R_white)/2", "R_odd=(R_black-R_white)/2"],
            "two_orientation_primary_state_order": ["R_S_even", "R_D_even", "R_S_odd", "R_D_odd"],
            "primary_candidate": "locally-subtracted nonlocal mixed incidence response R",
            "primary_null": "zero four-vector; q2/Jordan imply it only with an explicit factor-additive topology bridge",
            "secondary_mechanism_rows": ["ambient_even", "ambient_odd"],
            "two_orientation_transform": "form S/D across the two N650 orientations inside every same-stream batch",
            "covariance": "joint delete-one 4x4 covariance of the primary state; ambient rows remain correlated secondary diagnostics",
            "next_action": "bounded C++ fixture only; production remains stopped",
        },
        "evidence_boundary": {
            "exact": "path symmetry, typed involution oddness, nonzero tiny witnesses, and local moments",
            "mechanism_inference": "R_nonlocal removes isolated-fiber incidence cycles and targets connectivity-mediated C2xC5 interaction",
            "not_claimed": "a continuum memory field or q2/Jordan violation",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = render()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
