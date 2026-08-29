#!/usr/bin/env python3
"""Exact CRT/commutator oracle for the N65 -> N650 Gaussian cover square.

The point is structural: linear fiber averages and character projectors for
the coprime norm-two and norm-five factors commute before any probability or
scaling limit is introduced.  A small partition-lattice calculation records
the first nonlinear mixed object that can remain nonzero.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import reduce
import json
from math import gcd
from pathlib import Path
from typing import Callable, Iterable


Gaussian = tuple[int, int]
Matrix = tuple[tuple[int, int], tuple[int, int]]
DeckPoint = tuple[int, int]
Partition = tuple[tuple[DeckPoint, ...], ...]


ALPHA: Gaussian = (1, 1)   # 1+i, norm 2
BETA: Gaussian = (2, -1)   # 2-i, norm 5
SOURCES: tuple[Gaussian, ...] = ((8, 1), (7, 4))


def gaussian_multiply(left: Gaussian, right: Gaussian) -> Gaussian:
    a, b = left
    c, d = right
    return a * c - b * d, a * d + b * c


def gaussian_norm(value: Gaussian) -> int:
    return value[0] * value[0] + value[1] * value[1]


def multiplication_matrix(value: Gaussian) -> Matrix:
    a, b = value
    return ((a, -b), (b, a))


def determinant(matrix: Matrix) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def smith_invariants(matrix: Matrix) -> tuple[int, int]:
    d1 = reduce(gcd, (abs(item) for row in matrix for item in row))
    d2 = abs(determinant(matrix)) // d1
    if d1 == 0 or d2 % d1:
        raise ValueError("matrix does not define a full-rank Smith quotient")
    return d1, d2


def gaussian_text(value: Gaussian) -> str:
    a, b = value
    if b == 0:
        return str(a)
    sign = "+" if b > 0 else "-"
    magnitude = "" if abs(b) == 1 else str(abs(b))
    return f"{a}{sign}{magnitude}i"


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def deck_points() -> tuple[DeckPoint, ...]:
    """CRT coordinates for R/(alpha beta) = Z/2 x Z/5."""

    return tuple((a, b) for a in range(2) for b in range(5))


def average_over_factor(
    values: dict[DeckPoint, Fraction], factor: int
) -> dict[DeckPoint, Fraction]:
    """Average over the order-2 or order-5 translation subgroup."""

    if factor == 2:
        return {
            point: sum(
                (values[(a, point[1])] for a in range(2)), Fraction()
            )
            / 2
            for point in deck_points()
        }
    if factor == 5:
        return {
            point: sum(
                (values[(point[0], b)] for b in range(5)), Fraction()
            )
            / 5
            for point in deck_points()
        }
    raise ValueError("factor must be 2 or 5")


def canonical_partition(blocks: Iterable[Iterable[DeckPoint]]) -> Partition:
    normalized = [tuple(sorted(block)) for block in blocks]
    return tuple(sorted(normalized, key=lambda block: block[0]))


def discrete_partition() -> Partition:
    return canonical_partition(((point,) for point in deck_points()))


def factor_relation(factor: int) -> Partition:
    if factor == 2:
        return canonical_partition(
            (((a, b) for a in range(2)) for b in range(5))
        )
    if factor == 5:
        return canonical_partition(
            (((a, b) for b in range(5)) for a in range(2))
        )
    raise ValueError("factor must be 2 or 5")


def join_partitions(left: Partition, right: Partition) -> Partition:
    parent = {point: point for point in deck_points()}

    def find(point: DeckPoint) -> DeckPoint:
        while parent[point] != point:
            parent[point] = parent[parent[point]]
            point = parent[point]
        return point

    def union(first: DeckPoint, second: DeckPoint) -> None:
        root_first = find(first)
        root_second = find(second)
        if root_first != root_second:
            parent[max(root_first, root_second)] = min(root_first, root_second)

    for partition in (left, right):
        for block in partition:
            for point in block[1:]:
                union(block[0], point)

    output: dict[DeckPoint, list[DeckPoint]] = {}
    for point in deck_points():
        output.setdefault(find(point), []).append(point)
    return canonical_partition(output.values())


def partition_rank(partition: Partition) -> int:
    """Graphic/equivalence rank |D|-number_of_blocks."""

    return len(deck_points()) - len(partition)


def _lattice_record(value: Gaussian) -> dict:
    matrix = multiplication_matrix(value)
    return {
        "gaussian_period": gaussian_text(value),
        "coordinates": list(value),
        "index": gaussian_norm(value),
        "basis_matrix": [list(row) for row in matrix],
        "smith_invariants": list(smith_invariants(matrix)),
    }


def render() -> dict:
    gamma = gaussian_multiply(ALPHA, BETA)
    if gamma != gaussian_multiply(BETA, ALPHA) or gamma != (3, 1):
        raise AssertionError("Gaussian factor product is not 3+i")

    # Bezout certificate: 1 = i beta - alpha^2.
    i_beta = gaussian_multiply((0, 1), BETA)
    alpha_squared = gaussian_multiply(ALPHA, ALPHA)
    bezout_value = (i_beta[0] - alpha_squared[0], i_beta[1] - alpha_squared[1])
    if bezout_value != (1, 0):
        raise AssertionError("coprimality certificate failed")

    paths = []
    for source in SOURCES:
        after_alpha = gaussian_multiply(source, ALPHA)
        after_beta = gaussian_multiply(source, BETA)
        final = gaussian_multiply(source, gamma)
        paths.append(
            {
                "source": _lattice_record(source),
                "path_B_norm2_then_norm5": {
                    "intermediate": _lattice_record(after_alpha),
                    "final": _lattice_record(final),
                },
                "path_A_norm5_then_norm2": {
                    "intermediate": _lattice_record(after_beta),
                    "final": _lattice_record(final),
                },
            }
        )

    # An exact arbitrary fiber field; equality for this field is an oracle,
    # while the theorem in the note proves equality for every field.
    values = {
        (a, b): Fraction(17 * a + 3 * b * b - 5 * a * b - 2 * b, 7)
        for a, b in deck_points()
    }
    average_2_then_5 = average_over_factor(average_over_factor(values, 2), 5)
    average_5_then_2 = average_over_factor(average_over_factor(values, 5), 2)
    global_mean = sum(values.values(), Fraction()) / len(values)
    if average_2_then_5 != average_5_then_2:
        raise AssertionError("fiber Fubini identity failed")
    if any(value != global_mean for value in average_2_then_5.values()):
        raise AssertionError("two-stage average is not the full fiber average")

    # Fourier projectors are exact diagonal masks on character labels (r,s).
    # The two masks commute for all 2*5 joint characters.
    projector_checks = []
    for target_2 in range(2):
        for target_5 in range(5):
            for mode_2 in range(2):
                for mode_5 in range(5):
                    p2 = int(mode_2 == target_2)
                    p5 = int(mode_5 == target_5)
                    if p2 * p5 != p5 * p2:
                        raise AssertionError("character masks did not commute")
            projector_checks.append(
                {
                    "target_character": [target_2, target_5],
                    "commutator_zero_exact": True,
                }
            )

    base = discrete_partition()
    relation_2 = factor_relation(2)
    relation_5 = factor_relation(5)
    join_2_then_5 = join_partitions(join_partitions(base, relation_2), relation_5)
    join_5_then_2 = join_partitions(join_partitions(base, relation_5), relation_2)
    if join_2_then_5 != join_5_then_2:
        raise AssertionError("partition joins must commute")
    mixed_rank_defect = (
        partition_rank(join_2_then_5)
        - partition_rank(join_partitions(base, relation_2))
        - partition_rank(join_partitions(base, relation_5))
        + partition_rank(base)
    )

    return {
        "schema": "matching-one.gaussian-crt-commutator.v1",
        "issue": 200,
        "status": "exact_finite_group_result",
        "factors": {
            "alpha_norm2": _lattice_record(ALPHA),
            "beta_norm5": _lattice_record(BETA),
            "gamma_product": _lattice_record(gamma),
            "bezout_certificate": {
                "identity": "1 = i*(2-i) - (1+i)^2",
                "value": list(bezout_value),
                "crt_idempotent_mod_alpha": "i*(2-i)",
                "crt_idempotent_mod_beta": "-(1+i)^2",
            },
        },
        "n650_lattice_paths": paths,
        "quotient_square": {
            "rings": "R=Z[i]",
            "inclusions": "z*alpha*beta*R subset z*alpha*R,z*beta*R subset z*R",
            "meet": "z*alpha*R intersect z*beta*R = z*alpha*beta*R",
            "join": "z*alpha*R + z*beta*R = z*R",
            "relative_deck_group": "zR/(z*alpha*beta*R) ~= R/(alpha*beta) ~= C2 x C5 ~= C10",
            "kernel_to_N325": "z*beta*R/(z*alpha*beta*R) ~= R/(alpha) ~= C2",
            "kernel_to_N130": "z*alpha*R/(z*alpha*beta*R) ~= R/(beta) ~= C5",
            "cartesian_and_cocartesian": True,
        },
        "linear_fiber_oracle": {
            "test_field_values": {
                f"{a},{b}": _fraction_text(value)
                for (a, b), value in values.items()
            },
            "global_mean": _fraction_text(global_mean),
            "E5_E2_equals_E2_E5_equals_full_average": True,
            "ordered_linear_mark": {
                "definition": "C_mark=(E5 E2-E2 E5)f",
                "prediction": "identically_zero_for_every_f",
                "evidence_level": "exact_Fubini_CRT_identity",
            },
        },
        "character_projector_oracle": {
            "basis": "characters indexed by (r mod 2,s mod 5)",
            "factor_projectors_are_diagonal_masks": True,
            "all_joint_targets": projector_checks,
            "commutator_zero_exact": True,
        },
        "join_oracle": {
            "full_partition_join_commutator_zero": True,
            "rank_values_for_discrete_input": {
                "h_base": partition_rank(base),
                "h_join_K2": partition_rank(join_partitions(base, relation_2)),
                "h_join_K5": partition_rank(join_partitions(base, relation_5)),
                "h_join_K2_K5": partition_rank(join_2_then_5),
            },
            "mixed_rank_defect": mixed_rank_defect,
            "mixed_defect_formula": "h(Pi join R2 join R5)-h(Pi join R2)-h(Pi join R5)+h(Pi)",
            "interpretation": "minimal symmetric nonlinear interaction; not an order commutator",
        },
        "scientific_boundary": {
            "proved": [
                "the coprime Gaussian quotient square is CRT",
                "linear fiber averages and character projectors commute pointwise",
                "complete equivalence-relation joins commute",
            ],
            "excluded_as_discriminator": "any antisymmetric mark made only from honest linear fiber projections",
            "allowed_high_risk_observable": "a symmetric mixed topology/join defect after retaining the full typed connectivity state",
            "extra_structure_needed_for_antisymmetry": "an explicitly nonfunctorial intermediate truncation, representative choice, pruning rule, or orientation-sensitive local coarse graining",
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
