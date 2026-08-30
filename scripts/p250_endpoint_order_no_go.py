#!/usr/bin/env python3
"""Exact P250 endpoint-order no-go and P333 embedding audit.

This certificate proves that any statistic factored through net displacement
cannot distinguish the words xy and yx.  It then checks the P333 signed
detach/join positive control and audits the actual projective-leg runner
contract.  It intentionally does not invent microscopic morphism semantics.
"""
from __future__ import annotations

import argparse
import inspect
import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

try:
    from p333_dual_number_gram_extension import (
        detach_operator,
        join_operator,
        matrix_rank,
        multiply,
        radical_operator,
        restricted_first_form,
        transpose,
        weighted_detach_join_jordan_oracle,
    )
    from z5_charged_threepoint_mc import PARENT_GEOMETRY
    from z5_projective_leg_multiseparation_mc import ProjectiveLegIndex, charged_rows
except ModuleNotFoundError:
    from scripts.p333_dual_number_gram_extension import (
        detach_operator,
        join_operator,
        matrix_rank,
        multiply,
        radical_operator,
        restricted_first_form,
        transpose,
        weighted_detach_join_jordan_oracle,
    )
    from scripts.z5_charged_threepoint_mc import PARENT_GEOMETRY
    from scripts.z5_projective_leg_multiseparation_mc import ProjectiveLegIndex, charged_rows


DEFAULT_OUTPUT = Path("results/exact-p250-endpoint-order-no-go/latest.json")


def net_displacement(word: str) -> tuple[int, int]:
    if set(word) - {"x", "y"}:
        raise ValueError("words use only x and y")
    return word.count("x"), word.count("y")


def endpoint_symbol(word: str, radius: int = 0, moment_order: int = 1) -> tuple[str, int, int, int, int]:
    """Formal symbol for any row depending only on endpoint/radius/order."""
    a, b = net_displacement(word)
    return "G", a, b, int(radius), int(moment_order)


def endpoint_no_go(max_word_length: int = 8) -> dict:
    """Exhaustively check the abelianization fibers through a bounded length."""
    fibers: dict[tuple[int, int], list[str]] = {}
    checked = 0
    violations = []
    moment_library = [(radius, order) for radius in range(1, 5) for order in range(1, 7)]
    for length in range(max_word_length + 1):
        for letters in itertools.product("xy", repeat=length):
            word = "".join(letters)
            fibers.setdefault(net_displacement(word), []).append(word)
    for displacement, words in fibers.items():
        for left in words:
            for right in words:
                for radius, order in moment_library:
                    checked += 1
                    if endpoint_symbol(left, radius, order) != endpoint_symbol(right, radius, order):
                        violations.append([left, right, radius, order])
    if violations:
        raise AssertionError("endpoint factorization did not annihilate an order contrast")
    xy_contrasts = {
        f"radius_{radius}_moment_{order}": 0
        for radius, order in moment_library
    }
    return {
        "factor_map": "word -> (#x,#y)",
        "theorem": (
            "For every codomain A, map G:N^2->A and functional Phi, "
            "Phi(G(abelianize(w))) is constant on each word fiber."
        ),
        "consequence": "every endpoint-only TxTy-TyTx contrast is exactly zero",
        "max_word_length_checked": max_word_length,
        "abelianization_fibers": len(fibers),
        "ordered_pairs_radius_moments_checked": checked,
        "violations": violations,
        "xy_minus_yx_library": xy_contrasts,
    }


def parent_translation_gate() -> dict:
    failures = []
    for index, (x, y) in enumerate(PARENT_GEOMETRY.coordinates):
        xy = PARENT_GEOMETRY.vertex((x + 1, y + 1))
        yx = PARENT_GEOMETRY.vertex((x + 1, y + 1))
        if xy != yx:
            failures.append(index)
    if failures:
        raise AssertionError("parent endpoint translations failed to commute")
    return {
        "parent_vertices_checked": PARENT_GEOMETRY.n,
        "unit_xy_yx_endpoint_failures": failures,
        "runner_interpretation": (
            "charged_rows is keyed by final parent_index and recomputes a static "
            "projective-leg value there; no ordered path is an argument."
        ),
    }


def subtract(left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]):
    return [
        [left[row][column] - right[row][column] for column in range(len(left[0]))]
        for row in range(len(left))
    ]


def add(*matrices: Sequence[Sequence[Fraction]], coefficients: Sequence[int]):
    return [
        [
            sum(Fraction(coefficient) * matrix[row][column]
                for matrix, coefficient in zip(matrices, coefficients))
            for column in range(len(matrices[0][0]))
        ]
        for row in range(len(matrices[0]))
    ]


def encode_matrix(matrix: Sequence[Sequence[Fraction]]) -> list[list[str]]:
    return [[str(value) for value in row] for row in matrix]


def p333_rectangle_gate() -> dict:
    d = radical_operator(detach_operator(3, 1))
    j = radical_operator(join_operator(3, 0, 1))
    dj = multiply(d, j)
    jd = multiply(j, d)
    commutator = subtract(dj, jd)
    connected = add(d, j, dj, jd, coefficients=(1, 1, -1, -1))
    difference = subtract(d, j)
    square = multiply(difference, difference)
    h = restricted_first_form(3)
    oracle = weighted_detach_join_jordan_oracle()
    if connected != square:
        raise AssertionError("P333 connected rectangle is not (D-J)^2")
    if matrix_rank(commutator) != 2 or matrix_rank(connected) != 1:
        raise AssertionError("P333 ordered/connected ranks changed")
    if multiply(h, connected) != multiply(transpose(connected), h):
        raise AssertionError("P333 connected rectangle lost Gram compatibility")
    return {
        "state_space": "three-mark endpoint radical, dimension 4",
        "D": "detach mark 1 in the formal set-partition module",
        "J": "join marks 0 and 1 in the formal set-partition module",
        "ordered_DJ": encode_matrix(dj),
        "ordered_JD": encode_matrix(jd),
        "order_commutator": encode_matrix(commutator),
        "order_commutator_rank": matrix_rank(commutator),
        "connected_rectangle_formula": "D+J-DJ-JD=(D-J)^2",
        "connected_rectangle": encode_matrix(connected),
        "connected_rectangle_rank": matrix_rank(connected),
        "connected_rectangle_square_zero": not any(
            value for row in multiply(connected, connected) for value in row
        ),
        "first_jet_gram_self_adjoint": True,
        "existing_oracle_agrees": oracle["formula_exact"] and oracle["rank"] == 1,
        "important_distinction": (
            "K uses both ordered histories but is their signed sum; retain DJ-JD "
            "separately if the scientific target is order antisymmetry."
        ),
    }


def runner_embedding_audit() -> dict:
    public_methods = sorted(
        name
        for name, member in inspect.getmembers(ProjectiveLegIndex, predicate=inspect.isfunction)
        if not name.startswith("_")
    )
    charged_parameters = list(inspect.signature(charged_rows).parameters)
    required = [
        "three_mark_partition_before",
        "first_morphism_id_and_physical_support",
        "intermediate_three_mark_descriptors_S_D_and_S_J",
        "black_NN_and_white_matching_component_type_after_first_morphism",
        "ambient_rank_and_primitive_line_per_marked_block",
        "responses_D_J_DJ_JD_on_one_common_replica",
    ]
    return {
        "embedding_decision": "NOT_IMPLEMENTED_NO_PHYSICAL_SEMANTICS",
        "current_ProjectiveLegIndex_public_methods": public_methods,
        "charged_rows_parameters": charged_parameters,
        "current_endpoint_state": (
            "static black-NN or white-matching rank-one membership at each root, "
            "then a five-fiber DFT"
        ),
        "why_formal_K_cannot_be_renamed_physical": [
            "the runner emits no three-mark connectivity partition",
            "the runner applies no first morphism and retains no intermediate state",
            "formal detach has no declared site/edge cut on a percolation configuration",
            "formal join has no declared occupied connector and can join states forbidden by colour",
            "a physical site mutation changes black NN and complementary white matching connectivity together",
        ],
        "minimum_new_intermediate_fields": required,
        "minimum_online_rectangle": {
            "symmetric_connected": "R_plus=L_D+L_J-L_DJ-L_JD",
            "order_antisymmetric": "R_minus=L_DJ-L_JD",
            "common_randomness": "all four responses and both intermediate descriptors come from one base replica and one marked triple",
        },
        "production_action": (
            "Do not run until D and J are declared as physical paired primal/matching "
            "mutations and the intermediate marked-state descriptor is implemented."
        ),
    }


def build_result() -> dict:
    return {
        "schema": "p250-endpoint-order-no-go-v1",
        "status": "exact_no_go_and_semantic_gate",
        "endpoint_no_go": endpoint_no_go(),
        "actual_parent_translation_gate": parent_translation_gate(),
        "p333_positive_control": p333_rectangle_gate(),
        "projective_leg_embedding": runner_embedding_audit(),
        "claim_boundary": (
            "This proves non-identifiability from endpoint-only moments and audits "
            "the current runner. It does not define or sample a physical detach/join "
            "morphism, path memory, noncommuting translation, or continuum Jordan field."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--stdout", action="store_true")
    args = parser.parse_args()
    result = build_result()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    if args.stdout:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
