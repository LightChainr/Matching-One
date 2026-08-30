#!/usr/bin/env python3
"""Exact two-module and confluent-extension certificate for P321.

This continues ``p321_homology_trace_certificate`` by keeping crossed and
ordinary closures in separate module blocks.  It tests central grading,
blockwise pull-through, and whether the Q derivative forces a unique Jordan
off-diagonal at widths 2--4.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Sequence

from noncrossing_connectivity_codec import noncrossing_states
from p321_homology_trace_certificate import (
    Matrix,
    action_matrix,
    exact_rank,
    identity,
    join_adjacent,
    matrix_add,
    matrix_multiply,
    matrix_stats,
    rotate_state,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE = ROOT / "analysis" / "p321_graded_closure_extension.json"
SCHEMA = "matching-one/p321-graded-closure-extension/v1"


def zero(rows: int, columns: int | None = None) -> Matrix:
    columns = rows if columns is None else columns
    return tuple(tuple(0 for _ in range(columns)) for _ in range(rows))


def block_matrix(
    top_left: Matrix,
    top_right: Matrix,
    bottom_left: Matrix,
    bottom_right: Matrix,
) -> Matrix:
    top = tuple(left + right for left, right in zip(top_left, top_right))
    bottom = tuple(left + right for left, right in zip(bottom_left, bottom_right))
    return top + bottom


def block_diagonal(top_left: Matrix, bottom_right: Matrix) -> Matrix:
    return block_matrix(
        top_left,
        zero(len(top_left), len(bottom_right)),
        zero(len(bottom_right), len(top_left)),
        bottom_right,
    )


def scalar_multiply(matrix: Matrix, scalar: int) -> Matrix:
    return tuple(tuple(scalar * value for value in row) for row in matrix)


def intertwiner_constraints(left_actions: Sequence[Matrix], right_actions: Sequence[Matrix]) -> Matrix:
    """Equations for X A_k = B_k X, with row-major X variables."""

    if len(left_actions) != len(right_actions):
        raise ValueError("left and right action lists must have equal length")
    size = len(left_actions[0])
    equations: list[tuple[int, ...]] = []
    for left, right in zip(left_actions, right_actions):
        for row in range(size):
            for column in range(size):
                equation = [0] * (size * size)
                for pivot in range(size):
                    equation[row * size + pivot] += left[pivot][column]
                    equation[pivot * size + column] -= right[row][pivot]
                equations.append(tuple(equation))
    return tuple(equations)


def exact_nullity(equations: Matrix) -> int:
    variables = len(equations[0]) if equations else 0
    return variables - exact_rank(equations)


def _fraction_rank(matrix: Sequence[Sequence[Fraction]]) -> int:
    work = [list(row) for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, rows) if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [
                value - scale * pivot_value
                for value, pivot_value in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
    return pivot_row


def verify_singular_confluence(size: int, off_diagonal: Matrix) -> dict[str, Any]:
    """Check one exact rational point of the universal singular-basis identity."""

    delta = Fraction(2, 3)
    unit = tuple(tuple(Fraction(int(row == column)) for column in range(size)) for row in range(size))
    null = tuple(tuple(Fraction(0) for _ in range(size)) for _ in range(size))
    x = tuple(tuple(Fraction(value) for value in row) for row in off_diagonal)

    def blocks(a: Sequence[Sequence[Fraction]], b: Sequence[Sequence[Fraction]],
               c: Sequence[Sequence[Fraction]], d: Sequence[Sequence[Fraction]]) -> tuple[tuple[Fraction, ...], ...]:
        return tuple(tuple(left) + tuple(right) for left, right in zip(a, b)) + tuple(
            tuple(left) + tuple(right) for left, right in zip(c, d)
        )

    def multiply(left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]) -> tuple[tuple[Fraction, ...], ...]:
        columns = tuple(zip(*right))
        return tuple(
            tuple(sum(a * b for a, b in zip(row, column)) for column in columns)
            for row in left
        )

    s_delta = blocks(unit, x, null, tuple(tuple(delta * value for value in row) for row in unit))
    s_inverse = blocks(
        unit,
        tuple(tuple(-value / delta for value in row) for row in x),
        null,
        tuple(tuple(value / delta for value in row) for row in unit),
    )
    weight = blocks(
        unit,
        null,
        null,
        tuple(tuple((1 + delta) * value for value in row) for row in unit),
    )
    expected = blocks(
        unit,
        x,
        null,
        tuple(tuple((1 + delta) * value for value in row) for row in unit),
    )
    actual = multiply(multiply(s_delta, weight), s_inverse)
    residual = [[actual[row][column] - expected[row][column] for column in range(2 * size)] for row in range(2 * size)]
    return {
        "identity": "S_X(delta) diag(I,(1+delta)I) S_X(delta)^-1=[[I,X],[0,(1+delta)I]]",
        "rational_check_delta": "2/3",
        "residual_rank": _fraction_rank(residual),
        "limit": "I+N_X with N_X=[[0,X],[0,0]] and N_X^2=0",
    }


def width_certificate(width: int) -> dict[str, Any]:
    states = noncrossing_states(width)
    size = len(states)
    unit = identity(size)
    null = zero(size)
    translation = action_matrix(width, lambda state: rotate_state(state, 1))
    joins = tuple(
        action_matrix(width, lambda state, site=site: join_adjacent(state, site))
        for site in range(width)
    )

    doubled_translation = block_diagonal(translation, translation)
    closure_seam = block_diagonal(translation, unit)
    # J_Q=J_0+Q J_1=diag(I,-Q I).  W_Q=Gamma J_Q=diag(I,QI).
    j_constant = block_diagonal(unit, null)
    j_q_coefficient = block_diagonal(null, scalar_multiply(unit, -1))
    w_derivative = block_diagonal(null, unit)
    j_derivative_squared = matrix_multiply(j_q_coefficient, j_q_coefficient)
    w_derivative_idempotency = matrix_add(matrix_multiply(w_derivative, w_derivative), w_derivative, -1)

    residual_rows = []
    for site, generator in enumerate(joins):
        shifted = joins[(site + 1) % width]
        doubled_generator = block_diagonal(generator, generator)
        target_generator = block_diagonal(shifted, generator)
        d_constant = matrix_multiply(j_constant, closure_seam)
        d_q_coefficient = matrix_multiply(j_q_coefficient, closure_seam)
        residual_rows.append(
            {
                "J_constant_commutator": matrix_stats(
                    matrix_add(
                        matrix_multiply(j_constant, doubled_generator),
                        matrix_multiply(doubled_generator, j_constant),
                        -1,
                    )
                ),
                "J_Q_coefficient_commutator": matrix_stats(
                    matrix_add(
                        matrix_multiply(j_q_coefficient, doubled_generator),
                        matrix_multiply(doubled_generator, j_q_coefficient),
                        -1,
                    )
                ),
                "closure_seam_pull_through": matrix_stats(
                    matrix_add(
                        matrix_multiply(closure_seam, doubled_generator),
                        matrix_multiply(target_generator, closure_seam),
                        -1,
                    )
                ),
                "graded_defect_constant_pull_through": matrix_stats(
                    matrix_add(
                        matrix_multiply(d_constant, doubled_generator),
                        matrix_multiply(target_generator, d_constant),
                        -1,
                    )
                ),
                "graded_defect_Q_coefficient_pull_through": matrix_stats(
                    matrix_add(
                        matrix_multiply(d_q_coefficient, doubled_generator),
                        matrix_multiply(target_generator, d_q_coefficient),
                        -1,
                    )
                ),
            }
        )
    if not all(row == residual_rows[0] for row in residual_rows):
        raise AssertionError("translation-related block residuals must agree")

    canonical_nilpotent = block_matrix(null, translation, null, null)
    canonical_nilpotent_pull_rows = []
    for site, generator in enumerate(joins):
        shifted = joins[(site + 1) % width]
        doubled_generator = block_diagonal(generator, generator)
        target_generator = block_diagonal(shifted, generator)
        canonical_nilpotent_pull_rows.append(
            matrix_stats(
                matrix_add(
                    matrix_multiply(canonical_nilpotent, doubled_generator),
                    matrix_multiply(target_generator, canonical_nilpotent),
                    -1,
                )
            )
        )
    if not all(row == canonical_nilpotent_pull_rows[0] for row in canonical_nilpotent_pull_rows):
        raise AssertionError("translation-related nilpotent residuals must agree")

    sigma_constraints = intertwiner_constraints(joins, joins[1:] + joins[:1])
    affine_sigma_constraints = sigma_constraints + intertwiner_constraints(
        (translation,), (translation,)
    )
    commutant_constraints = intertwiner_constraints(joins, joins)
    affine_commutant_constraints = commutant_constraints + intertwiner_constraints(
        (translation,), (translation,)
    )

    return {
        "width": width,
        "module_dimension": size,
        "doubled_dimension": 2 * size,
        "block_pull_through": {
            "sites_checked": width,
            "all_sites_equal": True,
            **residual_rows[0],
        },
        "canonical_off_diagonal": {
            "X": "T",
            "rank_N_X": exact_rank(canonical_nilpotent),
            "N_X_squared": matrix_stats(matrix_multiply(canonical_nilpotent, canonical_nilpotent)),
            "pull_through": canonical_nilpotent_pull_rows[0],
            "translation_commutator": matrix_stats(
                matrix_add(
                    matrix_multiply(canonical_nilpotent, doubled_translation),
                    matrix_multiply(doubled_translation, canonical_nilpotent),
                    -1,
                )
            ),
        },
        "intertwiner_dimensions_over_Q": {
            "join_sigma_intertwiners_Xe_i=e_(i+1)X": exact_nullity(sigma_constraints),
            "affine_sigma_intertwiners_also_XT=TX": exact_nullity(affine_sigma_constraints),
            "join_commutant": exact_nullity(commutant_constraints),
            "affine_commutant": exact_nullity(affine_commutant_constraints),
        },
        "regular_Q_derivative": {
            "J_prime": "diag(0,-I)",
            "rank_J_prime": exact_rank(j_q_coefficient),
            "rank_J_prime_squared": exact_rank(j_derivative_squared),
            "W_prime": "diag(0,I)",
            "W_prime_idempotency_residual": matrix_stats(w_derivative_idempotency),
            "decision": "diagonal/idempotent, not nilpotent",
        },
        "singular_confluence": verify_singular_confluence(size, translation),
    }


def build_certificate(widths: Sequence[int] = (2, 3, 4)) -> dict[str, Any]:
    rows = [width_certificate(width) for width in widths]
    return {
        "schema": SCHEMA,
        "data_class": "exact integer/rational block matrices; symbolic Q split into coefficients",
        "construction": {
            "module": "M_cross direct_sum M_trivial",
            "grading": "Gamma=diag(I,-I)",
            "weight": "W_Q=diag(I,QI)",
            "supertrace_label": "J_Q=Gamma W_Q=diag(I,-Q I)",
            "closure_seam": "C=diag(T,I)",
            "graded_closure_defect": "D_Q=J_Q C=diag(T,-Q I)",
            "evaluation": "Tr(D_Q diag(A,A))=Tr(TA)-QTr(A)",
        },
        "widths": rows,
        "decision": {
            "block_pull_through": (
                "yes: separating the closures converts the failed single-space sum into an exact "
                "blockwise intertwiner from diag(e_i,e_i) to diag(e_(i+1),e_i)"
            ),
            "regular_Q_derivative": (
                "does not generate a nilpotent; it is the diagonal sector projector with an explicit sign"
            ),
            "Jordan_uniqueness": (
                "no: a Jordan limit requires singular confluent identification, and the admissible "
                "off-diagonal intertwiner space has dimension greater than one at every checked width"
            ),
            "F_t_classification": (
                "the exact homology trace first requires independent scalar amplitudes in the crossed and "
                "trivial modules; a Jordan matrix element is additional extension data, not forced by Q differentiation"
            ),
        },
        "q_lift_crosswalk": (
            "The L_hom to L_CP tangent shift -pi_0D is the diagonal derivative of the sector weight. "
            "It is not evidence for an off-diagonal logarithmic partner."
        ),
        "claim_boundary": (
            "The intertwiner multiplicities are for the checked periodic connectivity-join quotient. "
            "The block pull-through identity is representation-independent for any affine-TL action with "
            "T e_i T^-1=e_(i+1), but a fuller detach/loop-weight algebra may reduce the Hom spaces."
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--print", action="store_true", dest="print_certificate")
    args = parser.parse_args(argv)
    certificate = build_certificate()
    expected = json.loads(args.check.read_text(encoding="utf-8"))
    if expected != certificate:
        raise SystemExit(f"certificate drifted from {args.check}")
    if args.print_certificate:
        print(json.dumps(certificate, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
