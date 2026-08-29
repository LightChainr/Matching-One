#!/usr/bin/env python3
"""Exact doubled-space intertwiner identity for the Q->1 matching tangent."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from q1_matching_derivative_defect import component_count, induced_component_count


Edge = tuple[int, int]


def text(value: Fraction | int) -> str:
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def complement(mask: int, vertex_count: int) -> int:
    return ((1 << vertex_count) - 1) ^ mask


def cluster_vector(vertex_count: int, edges: list[Edge]) -> list[int]:
    return [
        induced_component_count(vertex_count, edges, mask)
        for mask in range(1 << vertex_count)
    ]


def doubled_tangent(left: list[int], right: list[int]) -> list[Fraction]:
    return [Fraction(value) for value in left + right]


def exchange_map(index: int, vertex_count: int) -> int:
    size = 1 << vertex_count
    if index < size:
        return size + complement(index, vertex_count)
    return complement(index - size, vertex_count)


def conjugate_diagonal(vector: list[Fraction], vertex_count: int) -> list[Fraction]:
    return [vector[exchange_map(i, vertex_count)] for i in range(len(vector))]


def even_odd_tangents(vector: list[Fraction], vertex_count: int) -> tuple[list[Fraction], list[Fraction]]:
    swapped = conjugate_diagonal(vector, vertex_count)
    even = [(a + b) / 2 for a, b in zip(vector, swapped)]
    odd = [(a - b) / 2 for a, b in zip(vector, swapped)]
    return even, odd


def commutator_coefficients(vector: list[Fraction], vertex_count: int) -> list[Fraction]:
    """Coefficients of [K,J]e_i along e_{J(i)} for diagonal K and exchange J."""
    return [vector[exchange_map(i, vertex_count)] - vector[i] for i in range(len(vector))]


def occupancy_buckets(left: list[int], right: list[int], vertex_count: int) -> dict[int, list[int]]:
    buckets: dict[int, set[int]] = defaultdict(set)
    for mask, value in enumerate(left):
        buckets[mask.bit_count()].add(value - right[complement(mask, vertex_count)])
    return {n: sorted(values) for n, values in sorted(buckets.items())}


def subset(mask: int, vertex_count: int) -> list[int]:
    return [i for i in range(vertex_count) if (mask >> i) & 1]


def site_c4_k4_oracle() -> dict:
    vertex_count = 4
    cycle4 = [(0, 1), (1, 2), (2, 3), (3, 0)]
    complete4 = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    left = cluster_vector(vertex_count, cycle4)
    right = cluster_vector(vertex_count, complete4)
    tangent = doubled_tangent(left, right)
    even, odd = even_odd_tangents(tangent, vertex_count)
    commutator = commutator_coefficients(tangent, vertex_count)
    twice_odd_J = [2 * odd[exchange_map(i, vertex_count)] for i in range(len(tangent))]
    buckets = occupancy_buckets(left, right, vertex_count)

    adjacent, opposite = 0b0011, 0b0101
    score_at_half = lambda mask: 4 * mask.bit_count() - 2 * vertex_count

    return {
        "graphs": {"left": "C4 occupied", "right": "K4 vacant matching"},
        "configuration_space_dimension": len(tangent),
        "bare_exchange": {
            "definition": "J|G,A>=|Ghat,A^c>; J|Ghat,B>=|G,B^c>",
            "J_squared_is_identity": all(
                exchange_map(exchange_map(i, vertex_count), vertex_count) == i
                for i in range(len(tangent))
            ),
        },
        "operator_identity": {
            "K": "diag(k_C4(A), k_K4(B))",
            "S": "(K+J K J)/2",
            "D": "(K-J K J)/2",
            "Ad_J_S_equals_S": conjugate_diagonal(even, vertex_count) == even,
            "Ad_J_D_equals_minus_D": conjugate_diagonal(odd, vertex_count) == [-x for x in odd],
            "d_logQ_pullthrough_residual_at_1": "[K,J]=2 D J",
            "identity_verified": commutator == twice_odd_J,
            "nonzero_residual_rank": sum(value != 0 for value in commutator),
            "max_abs_residual": text(max(abs(value) for value in commutator)),
            "all_order_identity": "exp(tK)J-Jexp(tK)=2 exp(tS)sinh(tD)J, t=log Q",
        },
        "occupancy_counterterm_no_go": {
            "candidate_class": "any scalar ell(|A|) added to the right tangent",
            "defect_values_by_occupied_count": {str(n): values for n, values in buckets.items()},
            "exists": all(len(values) == 1 for values in buckets.values()),
            "witness_count": 2,
            "witnesses": [
                {
                    "shape": "adjacent pair",
                    "occupied": subset(adjacent, vertex_count),
                    "defect": left[adjacent] - right[complement(adjacent, vertex_count)],
                    "Bernoulli_measure_score_at_p_half": score_at_half(adjacent),
                },
                {
                    "shape": "opposite pair",
                    "occupied": subset(opposite, vertex_count),
                    "defect": left[opposite] - right[complement(opposite, vertex_count)],
                    "Bernoulli_measure_score_at_p_half": score_at_half(opposite),
                },
            ],
            "consequence": "the Q matching tangent is not an occupancy-only measure derivative",
        },
        "trace_selection": {
            "uniform_doubled_trace_D": text(sum(odd)),
            "uniform_doubled_trace_S_times_D": text(sum(a * b for a, b in zip(even, odd))),
            "meaning": "J-even and J-odd diagonal operators are exactly orthogonal in the J-invariant doubled trace",
        },
    }


def edge_dual_control() -> dict:
    triangle = [(0, 1), (1, 2), (2, 0)]
    triple_edge = [(0, 1), (0, 1), (0, 1)]
    differences = []
    adjusted = []
    for mask in range(8):
        occupied = [edge for i, edge in enumerate(triangle) if (mask >> i) & 1]
        vacant_dual = [edge for i, edge in enumerate(triple_edge) if not ((mask >> i) & 1)]
        value = component_count(3, occupied) - component_count(2, vacant_dual)
        differences.append(value)
        adjusted.append(value - (2 - mask.bit_count()))
    return {
        "pair": "C3 and its three-parallel-edge planar dual",
        "raw_defect_values_by_occupied_count": {
            str(n): sorted({differences[m] for m in range(8) if m.bit_count() == n})
            for n in range(4)
        },
        "local_Euler_counterterm": "ell(n)=2-n",
        "adjusted_pullthrough_tangent_zero_configurationwise": all(value == 0 for value in adjusted),
    }


def build_oracle() -> dict:
    return {
        "schema": "matching-one.q1-matching-intertwiner.v1",
        "issue": 233,
        "exact_result": {
            "statement": (
                "S and D are exact +/- eigenoperators of the exchange superoperator Ad_J on the doubled "
                "configuration space. D is also exactly half the first log-Q pull-through obstruction: "
                "[K,J]=2DJ. It is not thereby a parity field of either single theory."
            ),
            "interface_status": (
                "J is an involutive configuration bijection and an intertwiner at Q=1 because all Q weights are one; "
                "its generic-Q intertwining fails at first order precisely when D is nonzero."
            ),
        },
        "edge_FK_positive_control": edge_dual_control(),
        "site_matching_square_face": site_c4_k4_oracle(),
        "connection_to_existing_derivative_ledger": {
            "issue_258_measure_derivative": (
                "A Bernoulli measure score depends only on |A|. The equal-size adjacent/opposite witness has equal "
                "measure score but different matching tangent, so the latter requires an explicit topological/operator insertion."
            ),
            "issue_257_colour_selection": (
                "Ad_J graph-complement parity and S_Q colour charge are distinct. The matching D operator is colour "
                "singlet, so the [2] null for V_(2,+/-2) does not remove this D channel; thermal Q4 remains eligible."
            ),
            "three_term_Q_derivative": [
                "critical-manifold/measure score",
                "projector or defect derivative",
                "explicit bare-observable derivative",
            ],
        },
        "claim_boundary": {
            "proved": [
                "the doubled-space involution and Ad_J grading",
                "the exact first-order pull-through identity",
                "the C3 planar-dual Euler repair",
                "the C4/K4 no-go for every occupancy-only counterterm",
                "the separation from a Bernoulli measure score",
            ],
            "not_proved": [
                "a bounded-local transfer-matrix seam for infinite square-site percolation",
                "identification of empirical P4[S']/P4[D'] with these exact finite-space operators",
                "a continuum Jordan or OPE interpretation",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    payload = json.dumps(build_oracle(), indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
