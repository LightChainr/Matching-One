#!/usr/bin/env python3
"""Exact configuration-level oracle for the two N650 Gaussian paths.

The proposed paths factor the same Gaussian cover, ``3+i=(2-i)(1+i)``.
This script decides what can be attached to one final binary configuration
without adding a coarse-graining convention.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from gaussian_cover_character_modes import (  # noqa: E402
    coset_key,
    matrix_multiply,
    multiplication_matrix,
    quotient_data,
    smith_invariants,
)


Vector = tuple[int, int]
Matrix = tuple[tuple[int, int], tuple[int, int]]

FACTOR_5: Vector = (2, -1)
FACTOR_2: Vector = (1, 1)
PRODUCT: Vector = (3, 1)
PARENTS: tuple[Vector, ...] = ((8, 1), (7, 4))
EXPECTED_FINAL: tuple[Vector, ...] = ((23, 11), (17, 19))


def group_summary(group: dict) -> dict:
    return {
        "label": group["label"],
        "multiplier": group["multiplier"],
        "multiplication_matrix": group["multiplication_matrix"],
        "norm_degree": group["norm_degree"],
        "smith_invariants": group["smith_invariants"],
        "additive_group": group["additive_group"],
        "group_exponent": group["group_exponent"],
        "element_representatives": group["element_representatives"],
    }


def vector_from_matrix(matrix: Matrix) -> Vector:
    return matrix[0][0], matrix[1][0]


def matrix_record(matrix: Matrix) -> list[list[int]]:
    return [list(row) for row in matrix]


def factor_projection(full: dict, factor: dict) -> list[int]:
    index = {
        coset_key(factor["_matrix"], representative): position
        for position, representative in enumerate(factor["_elements"])
    }
    return [
        index[coset_key(factor["_matrix"], representative)]
        for representative in full["_elements"]
    ]


def phase_exponent_10(phase: Fraction) -> int:
    scaled = phase * 10
    if scaled.denominator != 1:
        raise ArithmeticError("factor character did not embed in tenth roots")
    return scaled.numerator % 10


def coefficient_table(
    active: tuple[bool, ...],
    factor5: dict,
    factor2: dict,
    projection5: list[int],
    projection2: list[int],
    order: str,
) -> list[list[list[int]]]:
    """Nested exact DFT, stored as coefficients of tenth roots of unity."""

    values = {
        (projection5[index], projection2[index]): int(state)
        for index, state in enumerate(active)
    }
    output: list[list[list[int]]] = []
    for character5 in range(5):
        row: list[list[int]] = []
        for character2 in range(2):
            counts: Counter[int] = Counter()
            if order == "5_then_2":
                for element2 in range(2):
                    inner: Counter[int] = Counter()
                    for element5 in range(5):
                        if values[(element5, element2)]:
                            phase = factor5["_phase_table"][character5][element5]
                            inner[phase_exponent_10(phase)] += 1
                    shift = phase_exponent_10(
                        factor2["_phase_table"][character2][element2]
                    )
                    for exponent, coefficient in inner.items():
                        counts[(exponent + shift) % 10] += coefficient
            elif order == "2_then_5":
                for element5 in range(5):
                    inner = Counter()
                    for element2 in range(2):
                        if values[(element5, element2)]:
                            phase = factor2["_phase_table"][character2][element2]
                            inner[phase_exponent_10(phase)] += 1
                    shift = phase_exponent_10(
                        factor5["_phase_table"][character5][element5]
                    )
                    for exponent, coefficient in inner.items():
                        counts[(exponent + shift) % 10] += coefficient
            else:
                raise ValueError(f"unknown transform order: {order}")
            row.append([counts[index] for index in range(10)])
        output.append(row)
    return output


def descends(active: tuple[bool, ...], projection: list[int]) -> bool:
    states: dict[int, bool] = {}
    for value, target in zip(active, projection):
        if target in states and states[target] != value:
            return False
        states[target] = value
    return True


def nested_boolean_reduce(
    active: tuple[bool, ...],
    projection5: list[int],
    projection2: list[int],
    operation: str,
    order: str,
) -> bool:
    values = {
        (projection5[index], projection2[index]): state
        for index, state in enumerate(active)
    }
    reducer = any if operation == "OR" else all
    if operation not in {"OR", "AND"}:
        raise ValueError(operation)
    if order == "5_then_2":
        stage = [reducer(values[(element5, element2)] for element5 in range(5)) for element2 in range(2)]
    elif order == "2_then_5":
        stage = [reducer(values[(element5, element2)] for element2 in range(2)) for element5 in range(5)]
    else:
        raise ValueError(order)
    return reducer(stage)


def join_partitions(*partitions: list[int]) -> list[int]:
    size = len(partitions[0])
    parent = list(range(size))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    for partition in partitions:
        first_by_block: dict[int, int] = {}
        for index, block in enumerate(partition):
            if block in first_by_block:
                union(index, first_by_block[block])
            else:
                first_by_block[block] = index
    block_index: dict[int, int] = {}
    result = []
    for index in range(size):
        root = find(index)
        block_index.setdefault(root, len(block_index))
        result.append(block_index[root])
    return result


def partition_rank(partition: list[int]) -> int:
    return len(partition) - len(set(partition))


def exhaustive_fiber_oracle(
    full: dict, factor5: dict, factor2: dict
) -> dict:
    projection5 = factor_projection(full, factor5)
    projection2 = factor_projection(full, factor2)
    crt_pairs = list(zip(projection5, projection2))
    if len(set(crt_pairs)) != 10:
        raise AssertionError("CRT map is not bijective")

    descend5 = 0
    descend2 = 0
    descend_both = 0
    nonzero_commutators = 0
    direct_image_or_nonzero_commutators = 0
    direct_image_and_nonzero_commutators = 0
    ambiguity_example = None
    for mask in range(1 << 10):
        active = tuple(bool(mask & (1 << index)) for index in range(10))
        is5 = descends(active, projection5)
        is2 = descends(active, projection2)
        descend5 += int(is5)
        descend2 += int(is2)
        descend_both += int(is5 and is2)
        forward = coefficient_table(
            active, factor5, factor2, projection5, projection2, "5_then_2"
        )
        reverse = coefficient_table(
            active, factor5, factor2, projection5, projection2, "2_then_5"
        )
        nonzero_commutators += int(forward != reverse)
        direct_image_or_nonzero_commutators += int(
            nested_boolean_reduce(active, projection5, projection2, "OR", "5_then_2")
            != nested_boolean_reduce(active, projection5, projection2, "OR", "2_then_5")
        )
        direct_image_and_nonzero_commutators += int(
            nested_boolean_reduce(active, projection5, projection2, "AND", "5_then_2")
            != nested_boolean_reduce(active, projection5, projection2, "AND", "2_then_5")
        )
        if ambiguity_example is None and sum(active) == 1:
            occupied = active.index(True)
            ambiguity_example = {
                "full_element_index": occupied,
                "projects_to_Z5_cell": projection5[occupied],
                "projects_to_Z2_cell": projection2[occupied],
                "fiber_OR_value": 1,
                "fiber_AND_value": 0,
                "meaning": "a binary pushdown requires an extra rule",
            }

    expected = {"to_Z5_intermediate": 32, "to_Z2_intermediate": 4, "to_both": 2}
    observed = {
        "to_Z5_intermediate": descend5,
        "to_Z2_intermediate": descend2,
        "to_both": descend_both,
    }
    if (
        observed != expected
        or nonzero_commutators
        or direct_image_or_nonzero_commutators
        or direct_image_and_nonzero_commutators
    ):
        raise AssertionError("fiber exhaustive gate failed")

    discrete = list(range(10))
    join_5_then_2 = join_partitions(join_partitions(discrete, projection5), projection2)
    join_2_then_5 = join_partitions(join_partitions(discrete, projection2), projection5)
    if join_5_then_2 != join_2_then_5:
        raise AssertionError("partition joins did not commute")
    ranks = {
        "h_Pi": partition_rank(discrete),
        "h_Pi_join_R2": partition_rank(join_partitions(discrete, projection5)),
        "h_Pi_join_R5": partition_rank(join_partitions(discrete, projection2)),
        "h_Pi_join_R2_join_R5": partition_rank(join_5_then_2),
    }
    mixed_defect = (
        ranks["h_Pi_join_R2_join_R5"]
        - ranks["h_Pi_join_R2"]
        - ranks["h_Pi_join_R5"]
        + ranks["h_Pi"]
    )
    if ranks != {
        "h_Pi": 0,
        "h_Pi_join_R2": 5,
        "h_Pi_join_R5": 8,
        "h_Pi_join_R2_join_R5": 9,
    } or mixed_defect != -4:
        raise AssertionError("mixed partition-rank witness changed")

    return {
        "configuration_count": 1 << 10,
        "full_to_factor_indices": {
            "Z5": projection5,
            "Z2": projection2,
        },
        "crt_pairs_in_full_representative_order": [list(pair) for pair in crt_pairs],
        "crt_is_bijection": True,
        "descent_counts": observed,
        "expected_descent_counts": expected,
        "sequential_character_transform_nonzero_commutators": nonzero_commutators,
        "all_linear_character_projectors_commute_per_configuration": True,
        "functorial_boolean_direct_image": {
            "OR_nonzero_order_commutators": direct_image_or_nonzero_commutators,
            "AND_nonzero_order_commutators": direct_image_and_nonzero_commutators,
            "interpretation": "set image (OR), universal image (AND), and connectivity-equivalence joins compose in either order",
        },
        "minimal_binary_pushdown_ambiguity": ambiguity_example,
        "symmetric_mixed_partition_witness": {
            "partition": "discrete partition on C2 x C5",
            "rank_definition": "h(Pi)=10-number_of_blocks(Pi)",
            "join_orders_equal": True,
            "ranks": ranks,
            "Delta25_h": mixed_defect,
            "meaning": "the antisymmetric commutator is zero while a symmetric mixed interaction need not vanish",
        },
    }


def lineage_record(parent: Vector, expected: Vector) -> dict:
    parent_matrix = multiplication_matrix(parent)
    factor5_matrix = multiplication_matrix(FACTOR_5)
    factor2_matrix = multiplication_matrix(FACTOR_2)
    path_a_intermediate = matrix_multiply(parent_matrix, factor5_matrix)
    path_b_intermediate = matrix_multiply(parent_matrix, factor2_matrix)
    path_a_final = matrix_multiply(path_a_intermediate, factor2_matrix)
    path_b_final = matrix_multiply(path_b_intermediate, factor5_matrix)
    if path_a_final != path_b_final:
        raise AssertionError("Gaussian paths did not commute")
    if vector_from_matrix(path_a_final) != expected:
        raise AssertionError("unexpected N650 endpoint")
    if smith_invariants(path_a_final) != (1, 650):
        raise AssertionError("N650 endpoint should be cyclic")
    return {
        "parent_gaussian": list(parent),
        "parent_period_matrix": matrix_record(parent_matrix),
        "path_A": {
            "multipliers": ["2-i", "1+i"],
            "intermediate_gaussian": list(vector_from_matrix(path_a_intermediate)),
            "intermediate_period_matrix": matrix_record(path_a_intermediate),
        },
        "path_B": {
            "multipliers": ["1+i", "2-i"],
            "intermediate_gaussian": list(vector_from_matrix(path_b_intermediate)),
            "intermediate_period_matrix": matrix_record(path_b_intermediate),
        },
        "common_final_gaussian": list(expected),
        "common_final_period_matrix": matrix_record(path_a_final),
        "final_smith_invariants": [1, 650],
        "same_integer_period_graph_not_merely_isomorphic": True,
    }


def render() -> dict:
    full = quotient_data("3+i", PRODUCT, True)
    factor5 = quotient_data("2-i", FACTOR_5, True)
    factor2 = quotient_data("1+i", FACTOR_2, True)
    fiber = exhaustive_fiber_oracle(full, factor5, factor2)
    lineages = [
        lineage_record(parent, expected)
        for parent, expected in zip(PARENTS, EXPECTED_FINAL)
    ]
    return {
        "schema": "matching-one.p200-n650-path-operationalization.v1",
        "issue": 200,
        "status": "exact_obstruction_and_revised_acquisition_semantics",
        "gaussian_identity": "(2-i)(1+i)=(1+i)(2-i)=3+i",
        "deck_groups": {
            "final_over_parent": group_summary(full),
            "path_A_first_stage": group_summary(factor5),
            "path_B_first_stage": group_summary(factor2),
            "crt": "Z/10 is canonically represented here by the declared Gaussian basis as Z/5 x Z/2",
        },
        "n650_lineages": lineages,
        "single_parent_fiber_exhaustive_oracle": fiber,
        "exact_conclusions": {
            "endpoint_histogram_count": 1,
            "linear_deck_path_commutator": "identically zero on every configuration",
            "generic_binary_configuration_descends_to_intermediate": False,
            "intermediate_homology_flag_without_pushdown_rule": "undefined as a Bernoulli intermediate mask; functorial direct-image connectivity is defined and commutes",
            "old_C_mark_semantics": "exact zero for character projectors or honest functorial direct images; otherwise an extra nonfunctorial rule is missing",
        },
        "revised_acquisition": {
            "do_not_duplicate_endpoint_stream_by_path": True,
            "store_full_deck_label": "Z/10 label, equivalently the exact (Z/5,Z/2) CRT pair",
            "exact_zero_control": "Pi_Z5 Pi_Z2 - Pi_Z2 Pi_Z5 = 0 per configuration",
            "allowed_nonzero_linear_channel": "declared marked/covariant O_chibar times S_chi with full covariance",
            "allowed_invariant_quadratic_channel": "Bernoulli Hessian in chi tensor chibar, including diagonal correction",
            "minimal_nonlinear_mixed_candidate": "Delta25 h = h(Pi join R2 join R5)-h(Pi join R2)-h(Pi join R5)+h(Pi); symmetric, not an order commutator",
            "noncanonical_alternative": "freeze OR, AND, majority, or another pushdown as a new observable; never call its difference factorization memory",
            "production_decision": "do not extend the N650 C++ runner for the old path flag",
        },
        "evidence_boundary": {
            "exact": [
                "period matrices coincide pathwise",
                "Z10-to-Z5xZ2 CRT is bijective",
                "all 1024 fiber configurations have zero sequential-character commutator",
                "only 32/1024 and 4/1024 configurations descend to the two intermediates",
            ],
            "inference": "a nonzero morphism-memory experiment needs an explicitly new nonlinear or marked observable",
            "not_claimed": "absence of all possible nonlinear coarse-graining memory",
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
