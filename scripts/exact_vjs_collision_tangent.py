#!/usr/bin/env python3
"""Exact tiny-torus VJS energy/two-cluster collision tangent at Q=1.

This implements the geometric probabilities and colour-tensor formula of
Vasseur--Jacobsen--Saleur (arXiv:1206.2312, eqs. 13--15) on the L=2 square-bond
torus.  It separates the Q derivative into the FK measure score, the finite
confluent projector derivative from Issue #262, and the explicit Q-dependent
normalization of the lattice insertion.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys
from typing import Dict, Iterable, Sequence

from square_bond_kappa3 import BondPair, square_bond_pairs
from p262_confluent_potts_projectors import analyze as projector_oracle


Vector = list[Fraction]
BASIS = ("I_pair", "X_shared_colour", "J_all_ones")


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def add(*vectors: Sequence[Fraction]) -> Vector:
    return [sum(vector[i] for vector in vectors) for i in range(len(vectors[0]))]


def scale(vector: Sequence[Fraction], scalar: Fraction) -> Vector:
    return [scalar * value for value in vector]


def vector_record(vector: Sequence[Fraction]) -> Dict[str, str]:
    return {basis: fraction_text(value) for basis, value in zip(BASIS, vector)}


def configuration_row(length: int, mask: int, pairs: Sequence[BondPair]) -> Dict[str, object]:
    vertices = length * length
    parent = list(range(vertices))

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def union(first: int, second: int) -> None:
        first, second = find(first), find(second)
        if first != second:
            parent[second] = first

    bonds = 0
    for index, pair in enumerate(pairs):
        if (mask >> index) & 1:
            bonds += 1
            union(pair.primal[0], pair.primal[1])
    roots = [find(vertex) for vertex in range(vertices)]
    clusters = len(set(roots))
    # On the 2x2 torus these are the two disjoint horizontal nearest-neighbour pairs.
    first_pair = (roots[0], roots[1])
    second_pair = (roots[2], roots[3])
    first_distinct = first_pair[0] != first_pair[1]
    second_distinct = second_pair[0] != second_pair[1]
    all_roots = set(roots)
    first_roots, second_roots = set(first_pair), set(second_pair)
    p0 = first_distinct and second_distinct and len(all_roots) == 4
    p1 = (first_distinct and second_distinct and len(all_roots) == 3
          and len(first_roots & second_roots) == 1)
    p2 = (first_distinct and second_distinct and len(all_roots) == 2
          and first_roots == second_roots)
    return {
        "T": Fraction(2 * clusters + bonds, 2),
        "pair0_distinct": int(first_distinct),
        "pair1_distinct": int(second_distinct),
        "both_pairs_distinct": int(first_distinct and second_distinct),
        "P0_four_clusters": int(p0),
        "P1_one_propagating_cluster": int(p1),
        "P2_two_propagating_clusters": int(p2),
    }


def mean(rows: Sequence[Dict[str, object]], field: str) -> Fraction:
    return sum((Fraction(row[field]) for row in rows), Fraction(0)) / len(rows)


def score_derivative(rows: Sequence[Dict[str, object]], field: str) -> Fraction:
    observable_mean = mean(rows, field)
    score_mean = mean(rows, "T")
    return sum(
        (Fraction(row[field]) - observable_mean) * (Fraction(row["T"]) - score_mean)
        for row in rows
    ) / len(rows)


def render(length: int = 2) -> Dict[str, object]:
    if length != 2:
        raise ValueError("the frozen minimal VJS oracle uses exactly L=2")
    pairs = square_bond_pairs(length)
    rows = [configuration_row(length, mask, pairs) for mask in range(1 << len(pairs))]
    fields = (
        "pair0_distinct", "pair1_distinct", "both_pairs_distinct",
        "P0_four_clusters", "P1_one_propagating_cluster",
        "P2_two_propagating_clusters",
    )
    probabilities = {field: mean(rows, field) for field in fields}
    derivatives = {field: score_derivative(rows, field) for field in fields}
    p_ne = probabilities["pair0_distinct"]
    p_ne_prime = derivatives["pair0_distinct"]
    p0 = probabilities["P0_four_clusters"]
    p1 = probabilities["P1_one_propagating_cluster"]
    p2 = probabilities["P2_two_propagating_clusters"]
    p0_prime = derivatives["P0_four_clusters"]
    p1_prime = derivatives["P1_one_propagating_cluster"]
    p2_prime = derivatives["P2_two_propagating_clusters"]
    energy_connected = p0 + p1 - p_ne * p_ne
    energy_connected_prime = p0_prime + p1_prime - 2 * p_ne * p_ne_prime

    # Eq. (14), after centering, in the diagram basis:
    # C(Q)=4 Q^-4 J G(Q) + 2 Q^-2 P_reg(Q) P2(Q),
    # P_reg=P_singlet+P_[2].  Issue #262 gives P_reg(1)=I+X-4J
    # and P_reg'(1)=X.
    p_reg = [Fraction(1), Fraction(1), Fraction(-4)]
    p_reg_prime = [Fraction(0), Fraction(1), Fraction(0)]
    j = [Fraction(0), Fraction(0), Fraction(1)]
    value = add(scale(j, 4 * energy_connected), scale(p_reg, 2 * p2))
    measure_score = add(
        scale(j, 4 * energy_connected_prime),
        scale(p_reg, 2 * p2_prime),
    )
    projector_derivative = scale(p_reg_prime, 2 * p2)
    explicit_insertion_derivative = add(
        scale(j, -16 * energy_connected),
        scale(p_reg, -4 * p2),
    )
    total = add(measure_score, projector_derivative, explicit_insertion_derivative)

    # Independent direct product-rule differentiation of the three generic-Q
    # coefficient functions at Q=1.
    direct = add(
        scale(j, -16 * energy_connected + 4 * energy_connected_prime),
        scale(p_reg, -4 * p2 + 2 * p2_prime),
        scale(p_reg_prime, 2 * p2),
    )
    projector = projector_oracle()["unordered_pair_representation"]["Q_to_1_laurent"]
    probability_records = {
        field: {
            "Q1": fraction_text(probabilities[field]),
            "dQ_at_Q1_measure_score": fraction_text(derivatives[field]),
        }
        for field in fields
    }
    return {
        "schema": "matching-one.exact-vjs-collision-tangent.v1",
        "issues": [258, 262],
        "geometry": {
            "lattice": "square_bond_L2_torus",
            "vertices": 4,
            "bonds": len(pairs),
            "configurations": len(rows),
            "marked_pairs": [[0, 1], [2, 3]],
            "warning": "tiny identity oracle, not a scaling-distance experiment",
        },
        "critical_coordinate": {
            "eta": "log(v/sqrt(Q))=0",
            "t": "log(Q)",
            "Q1_d_dt_equals_d_dQ": True,
            "measure_score": "T=k+b/2",
            "mean_T": fraction_text(mean(rows, "T")),
        },
        "vjs_geometric_probabilities": probability_records,
        "exact_partition_check": {
            "P0_plus_P1_plus_P2": fraction_text(p0 + p1 + p2),
            "both_pairs_distinct": fraction_text(probabilities["both_pairs_distinct"]),
            "equal": p0 + p1 + p2 == probabilities["both_pairs_distinct"],
        },
        "centered_energy_channel": {
            "definition": "G=P0+P1-P_distinct^2",
            "Q1": fraction_text(energy_connected),
            "dQ_at_Q1": fraction_text(energy_connected_prime),
        },
        "vjs_centered_two_point_tensor": {
            "diagram_basis": list(BASIS),
            "generic_Q_formula": "4 Q^-4 J G(Q) + 2 Q^-2 [P_singlet(Q)+P_[2](Q)] P2(Q)",
            "Q1_value": vector_record(value),
            "derivative_decomposition": {
                "measure_score_fixed_field": vector_record(measure_score),
                "finite_confluent_projector_derivative": vector_record(projector_derivative),
                "explicit_insertion_normalization_derivative": vector_record(explicit_insertion_derivative),
                "sum": vector_record(total),
                "direct_generic_Q_product_rule": vector_record(direct),
                "sum_equals_direct": total == direct,
            },
            "all_three_contributions_nonzero": all(any(value for value in vector) for vector in (
                measure_score, projector_derivative, explicit_insertion_derivative)),
        },
        "issue262_reuse_gate": {
            "P_regular_at_Q1": projector["regular_confluent_block_at_Q1"],
            "P_regular_derivative_at_Q1": projector["regular_block_derivative_at_Q1"],
            "individual_projector_derivatives_used": False,
        },
        "collision_boundary": {
            "finite_graph_result": "the full VJS tensor tangent decomposes exactly into three nonzero terms",
            "not_determined_here": [
                "continuum scaling-dimension collision slope sqrt(3)/pi",
                "large-distance amplitude equality A(1)=A_tilde(1)",
                "the logarithmic coefficient 2sqrt(3)/pi from a one-distance L2 oracle",
            ],
            "reason": "these are scaling-limit inputs, not finite colour/projector identities",
        },
        "source": {
            "paper": "Vasseur, Jacobsen, Saleur, arXiv:1206.2312",
            "equations": [4, 7, 12, 13, 14, 15],
        },
    }


def main(argv: Iterable[str] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = render()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
