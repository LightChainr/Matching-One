#!/usr/bin/env python3
"""Exact finite-width certificate for the P321 homology trace candidate.

The calculation uses the periodic FK connectivity representation on circular
noncrossing partitions.  ``e_i`` is the idempotent which joins sites ``i`` and
``i+1`` and ``T`` translates every site by one step.  This is deliberately a
minimal representation-theoretic test, not a tube-algebra implementation.

Two exact closure functionals are available on a transfer word ``A``::

    Z_trivial(A) = Tr(A)
    Z_cross(A)   = Tr(T A)

Their formal homology contrast is represented inside this matrix module by
``D_hom = T-Q I``.  The certificate tests whether this representative is an
ordinary central element or one transparent crossed defect.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

from noncrossing_connectivity_codec import canonical_rgs, noncrossing_states


Matrix = tuple[tuple[int, ...], ...]
State = tuple[int, ...]

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE = ROOT / "analysis" / "p321_homology_trace_certificate.json"
SCHEMA = "matching-one/p321-homology-trace-functional/v1"


def identity(size: int) -> Matrix:
    return tuple(tuple(int(row == column) for column in range(size)) for row in range(size))


def matrix_add(left: Matrix, right: Matrix, right_scale: int = 1) -> Matrix:
    return tuple(
        tuple(a + right_scale * b for a, b in zip(left_row, right_row))
        for left_row, right_row in zip(left, right)
    )


def matrix_multiply(left: Matrix, right: Matrix) -> Matrix:
    columns = tuple(zip(*right))
    return tuple(
        tuple(sum(a * b for a, b in zip(row, column)) for column in columns)
        for row in left
    )


def matrix_power(matrix: Matrix, exponent: int) -> Matrix:
    if exponent < 0:
        raise ValueError("matrix exponent must be nonnegative")
    result = identity(len(matrix))
    for _ in range(exponent):
        result = matrix_multiply(result, matrix)
    return result


def matrix_trace(matrix: Matrix) -> int:
    return sum(matrix[index][index] for index in range(len(matrix)))


def exact_rank(matrix: Matrix) -> int:
    work = [[Fraction(value) for value in row] for row in matrix]
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


def matrix_stats(matrix: Matrix) -> dict[str, int]:
    values = [value for row in matrix for value in row]
    return {
        "nonzero_entries": sum(value != 0 for value in values),
        "rank_over_Q": exact_rank(matrix),
        "frobenius_norm_squared": sum(value * value for value in values),
    }


def rotate_state(state: State, steps: int = 1) -> State:
    width = len(state)
    rotated = [0] * width
    for site, label in enumerate(state):
        rotated[(site + steps) % width] = label
    return canonical_rgs(rotated)


def join_adjacent(state: State, site: int) -> State:
    width = len(state)
    other = (site + 1) % width
    left_label = state[site]
    right_label = state[other]
    if left_label == right_label:
        return state
    return canonical_rgs(
        left_label if label == right_label else label for label in state
    )


def action_matrix(width: int, action: Callable[[State], State]) -> Matrix:
    states = noncrossing_states(width)
    index = {state: rank for rank, state in enumerate(states)}
    output = [[0] * len(states) for _ in states]
    for column, state in enumerate(states):
        output[index[action(state)]][column] = 1
    return tuple(tuple(row) for row in output)


def _trace_difference_coefficients(translation: Matrix, matrix: Matrix) -> tuple[int, int]:
    """Return constant and Q coefficients of Tr((T-QI) matrix)."""

    return matrix_trace(matrix_multiply(translation, matrix)), -matrix_trace(matrix)


def _word(generators: dict[str, Matrix], names: Iterable[str]) -> Matrix:
    output = identity(len(next(iter(generators.values()))))
    for name in names:
        output = matrix_multiply(output, generators[name])
    return output


def width_certificate(width: int) -> dict[str, Any]:
    if width < 2:
        raise ValueError("the periodic certificate starts at width two")
    states = noncrossing_states(width)
    dimension = len(states)
    unit = identity(dimension)
    translation = action_matrix(width, lambda state: rotate_state(state, 1))
    inverse_translation = action_matrix(width, lambda state: rotate_state(state, -1))
    joins = tuple(
        action_matrix(width, lambda state, site=site: join_adjacent(state, site))
        for site in range(width)
    )

    local_rows: list[dict[str, Any]] = []
    for site, generator in enumerate(joins):
        translated_site = (site + 1) % width
        translated_generator = joins[translated_site]
        seam_covariance = matrix_add(
            matrix_multiply(matrix_multiply(translation, generator), inverse_translation),
            translated_generator,
            -1,
        )
        seam_pull_through = matrix_add(
            matrix_multiply(translation, generator),
            matrix_multiply(translated_generator, translation),
            -1,
        )
        d_commutator = matrix_add(
            matrix_multiply(translation, generator),
            matrix_multiply(generator, translation),
            -1,
        )
        # (T-QI)e_i-e_(i+1)(T-QI) = Q(e_(i+1)-e_i).
        d_pull_through_q_coefficient = matrix_add(translated_generator, generator, -1)
        d_translation_commutator = matrix_add(
            matrix_multiply(translation, translation),
            matrix_multiply(translation, translation),
            -1,
        )
        local_rows.append(
            {
                "seam_covariance_residual": matrix_stats(seam_covariance),
                "seam_pull_through_residual": matrix_stats(seam_pull_through),
                "D_hom_commutator": {
                    "polynomial_factor": "1",
                    **matrix_stats(d_commutator),
                },
                "D_hom_translation_commutator": matrix_stats(d_translation_commutator),
                "D_hom_pull_through_residual": {
                    "polynomial_factor": "Q",
                    **matrix_stats(d_pull_through_q_coefficient),
                },
            }
        )

    if not all(row == local_rows[0] for row in local_rows):
        raise AssertionError("translation-related local residuals must agree exactly")

    generators = {"T": translation, **{f"e{site}": matrix for site, matrix in enumerate(joins)}}
    crossed_twisted_residuals = []
    for left in generators.values():
        for right in generators.values():
            sigma_inverse_right = matrix_multiply(
                matrix_multiply(inverse_translation, right), translation
            )
            ordinary_product = matrix_multiply(left, right)
            twisted_product = matrix_multiply(sigma_inverse_right, left)
            residual = matrix_trace(matrix_multiply(translation, ordinary_product)) - matrix_trace(
                matrix_multiply(translation, twisted_product)
            )
            crossed_twisted_residuals.append(residual)

    # A minimal ordinary-cyclicity witness for Tr((T-QI) .).  Width two is
    # translation-degenerate; widths three and four already expose the seam.
    if width == 2:
        cyclic_witness = {
            "available": False,
            "reason": "T acts as identity on the two Catalan connectivity states",
            "constant_coefficient": 0,
            "Q_coefficient": 0,
        }
    else:
        right_word = ("T",) * (2 if width == 3 else 1) + ("e0",)
        left = generators["e0"]
        right = _word(generators, right_word)
        commutator = matrix_add(matrix_multiply(left, right), matrix_multiply(right, left), -1)
        constant, q_coefficient = _trace_difference_coefficients(translation, commutator)
        cyclic_witness = {
            "available": True,
            "A": "e0",
            "B": "*".join(right_word),
            "constant_coefficient": constant,
            "Q_coefficient": q_coefficient,
        }

    # The crossed trace obeys a sigma-twisted trace law, but subtracting Q
    # times the ordinary trace destroys it.  A=B=e0 is the smallest witness.
    left = joins[0]
    right = joins[0]
    sigma_inverse_right = matrix_multiply(
        matrix_multiply(inverse_translation, right), translation
    )
    twisted_difference = matrix_add(
        matrix_multiply(left, right), matrix_multiply(sigma_inverse_right, left), -1
    )
    twisted_constant, twisted_q = _trace_difference_coefficients(translation, twisted_difference)

    return {
        "width": width,
        "dimension": dimension,
        "translation_order_residual": matrix_stats(matrix_add(matrix_power(translation, width), unit, -1)),
        "local_generator_residuals": {
            "sites_checked": width,
            "all_sites_equal": True,
            **local_rows[0],
        },
        "trace_laws": {
            "crossed_trace_sigma_twisted_max_abs_residual_on_generators": max(
                map(abs, crossed_twisted_residuals), default=0
            ),
            "D_hom_ordinary_cyclicity_witness": cyclic_witness,
            "D_hom_sigma_twisted_witness_A_equals_B_equals_e0": {
                "constant_coefficient": twisted_constant,
                "Q_coefficient": twisted_q,
            },
        },
    }


def build_certificate(widths: Sequence[int] = (2, 3, 4)) -> dict[str, Any]:
    rows = [width_certificate(width) for width in widths]
    return {
        "schema": SCHEMA,
        "data_class": "exact integer matrices and polynomial-in-Q residuals; no Monte Carlo",
        "representation": {
            "basis": "circular noncrossing FK connectivity partitions",
            "e_i": "idempotent join of cyclic neighbours i and i+1",
            "T": "one-site translation",
            "closure_functionals": {
                "Z_trivial(A)": "Tr(A)",
                "Z_cross(A)": "Tr(T A)",
                "Phi_Q(A)": "Z_cross(A)-Q*Z_trivial(A)=Tr((T-QI)A)",
            },
        },
        "exact_identity": {
            "translation_covariance": "T e_i T^-1=e_(i+1)",
            "crossed_pull_through": "T e_i=e_(i+1) T",
            "candidate": "D_hom=T-QI",
            "candidate_pull_through": "D_hom e_i-e_(i+1)D_hom=Q(e_(i+1)-e_i)",
        },
        "widths": rows,
        "decision": {
            "ordinary_center": "no from the first nondegenerate width 3",
            "single_crossed_defect": "no; T is crossed, I is untwisted, and their Q-weighted difference obeys neither one law",
            "valid_realization": (
                "an exact difference of two closure traces, equivalently a central label only after adjoining "
                "the external direct sum M_cross plus M_trivial"
            ),
            "width_two": "degenerate and therefore not evidence for centrality",
            "continuum_consequence": (
                "the certificate names an additional homology-resolved module/trace; it does not determine "
                "F_t(tau) or the spinful identity-dressing shape"
            ),
        },
        "claim_boundary": (
            "The finite connectivity representation proves the algebraic obstruction and the crossed/ordinary "
            "trace decomposition.  Identifying these two formal closures with the physical Z_2D and Z_0D "
            "requires a calibrated affine-TL transfer representation with loop weights and homology bookkeeping."
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--print", action="store_true", dest="print_certificate")
    args = parser.parse_args(argv)
    certificate = build_certificate()
    if args.check:
        expected = json.loads(args.check.read_text(encoding="utf-8"))
        if expected != certificate:
            raise SystemExit(f"certificate drifted from {args.check}")
    if args.print_certificate:
        print(json.dumps(certificate, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
