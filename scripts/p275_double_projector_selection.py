#!/usr/bin/env python3
"""Exact rank-plane and Jordan-parity certificates for Issue #275.

The calculation is deliberately finite dimensional.  It does not try to
identify a lattice amplitude.  It proves which directions are allowed after
normalization and Alexander duality, and records the imported Potts projector
zero needed by the double-projector staircase.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Iterable, Sequence


Scalar = Fraction
Vector = tuple[Scalar, ...]
Matrix = tuple[Vector, ...]


ZERO = Fraction(0)
ONE = Fraction(1)

I3: Matrix = (
    (ONE, ZERO, ZERO),
    (ZERO, ONE, ZERO),
    (ZERO, ZERO, ONE),
)

# Column-vector convention in the ordered rank basis (P0, P1, P2).
ALEXANDER: Matrix = (
    (ZERO, ZERO, ONE),
    (ZERO, ONE, ZERO),
    (ONE, ZERO, ZERO),
)

TOTAL: Vector = (ONE, ONE, ONE)
A_TOP: Vector = (-ONE, ZERO, ONE)
E_TOP: Vector = (ONE, ZERO, ONE)

EVEN_TANGENT: Vector = (ONE, Fraction(-2), ONE)
ODD_TANGENT: Vector = (-ONE, ZERO, ONE)


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(len(right))) for j in range(len(right[0])))
        for i in range(len(left))
    )


def matvec(matrix: Matrix, vector: Vector) -> Vector:
    return tuple(sum(row[j] * vector[j] for j in range(len(vector))) for row in matrix)


def rowmat(row: Vector, matrix: Matrix) -> Vector:
    return tuple(sum(row[i] * matrix[i][j] for i in range(len(row))) for j in range(len(matrix[0])))


def dot(left: Vector, right: Vector) -> Scalar:
    return sum(a * b for a, b in zip(left, right, strict=True))


def scale(value: Scalar, vector: Vector) -> Vector:
    return tuple(value * entry for entry in vector)


def add(left: Matrix, right: Matrix) -> Matrix:
    return tuple(tuple(a + b for a, b in zip(x, y, strict=True)) for x, y in zip(left, right, strict=True))


def matrix_scale(value: Scalar, matrix: Matrix) -> Matrix:
    return tuple(tuple(value * entry for entry in row) for row in matrix)


def fraction_text(value: Scalar) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def vector_text(vector: Iterable[Scalar]) -> list[str]:
    return [fraction_text(value) for value in vector]


def matrix_text(matrix: Matrix) -> list[list[str]]:
    return [vector_text(row) for row in matrix]


def rank_plane_certificate() -> dict[str, object]:
    even_projector = matrix_scale(Fraction(1, 2), add(I3, ALEXANDER))
    odd_projector = matrix_scale(Fraction(1, 2), add(I3, matrix_scale(-ONE, ALEXANDER)))

    assert matmul(ALEXANDER, ALEXANDER) == I3
    assert matmul(even_projector, even_projector) == even_projector
    assert matmul(odd_projector, odd_projector) == odd_projector
    assert matmul(even_projector, odd_projector) == matrix_scale(ZERO, I3)
    assert rowmat(A_TOP, ALEXANDER) == scale(-ONE, A_TOP)
    assert rowmat(E_TOP, ALEXANDER) == E_TOP

    assert dot(TOTAL, EVEN_TANGENT) == ZERO
    assert dot(TOTAL, ODD_TANGENT) == ZERO
    assert matvec(ALEXANDER, EVEN_TANGENT) == EVEN_TANGENT
    assert matvec(ALEXANDER, ODD_TANGENT) == scale(-ONE, ODD_TANGENT)
    assert dot(A_TOP, EVEN_TANGENT) == ZERO
    assert dot(A_TOP, ODD_TANGENT) == Fraction(2)
    assert dot(E_TOP, ODD_TANGENT) == ZERO

    return {
        "basis": ["P0", "P1", "P2"],
        "alexander_involution": matrix_text(ALEXANDER),
        "even_projector": matrix_text(even_projector),
        "odd_projector": matrix_text(odd_projector),
        "normalized_even_line": vector_text(EVEN_TANGENT),
        "normalized_odd_line": vector_text(ODD_TANGENT),
        "A_top_covector": vector_text(A_TOP),
        "E_top_covector": vector_text(E_TOP),
        "contractions": {
            "A_top_on_even": fraction_text(dot(A_TOP, EVEN_TANGENT)),
            "A_top_on_odd": fraction_text(dot(A_TOP, ODD_TANGENT)),
            "E_top_on_odd": fraction_text(dot(E_TOP, ODD_TANGENT)),
        },
        "theorem": (
            "Every normalized Alexander-even tangent is a multiple of (1,-2,1); "
            "every normalized Alexander-odd tangent is a multiple of (-1,0,1)."
        ),
    }


def vacuum_kdv_certificate() -> dict[str, object]:
    # Arguin at Q=1 gives P2(tau)=P0(tau) for every modulus.  Applying a
    # componentwise linear tau-Ward operator W which kills the normalized
    # total leaves an even, zero-sum vector, hence the unique line below.
    response = EVEN_TANGENT
    assert response[2] == response[0]
    assert sum(response) == ZERO
    assert dot(A_TOP, response) == ZERO
    return {
        "exact_inputs": [
            "Arguin restricted-sector identity Z_cross(Q,tau)=Q Z_trivial(Q,tau)",
            "Q=1 normalization P0+P1+P2=1",
            "componentwise tau-Ward operator W commutes with rank exchange",
            "W[1]=0 (in particular K4=(delta-E2/6)delta)",
        ],
        "response_direction": vector_text(response),
        "identities": {
            "delta_P2_minus_delta_P0": fraction_text(response[2] - response[0]),
            "delta_P0_plus_delta_P1_plus_delta_P2": fraction_text(sum(response)),
            "A_top_response": fraction_text(dot(A_TOP, response)),
        },
        "status": "EXACT_ZERO",
        "conclusion": "The vacuum/KdV rank response is Alexander-even and A_top annihilates it at every tau.",
    }


def thermal_q4_certificate() -> dict[str, object]:
    # For the doubled primal/matching continuum family, eta is chosen odd:
    # P(eta,tau)=C P(-eta,tau).  One eta derivative at zero is C-odd.
    # A componentwise Q4/Ward operation commutes with the finite rank action.
    response = ODD_TANGENT
    assert matvec(ALEXANDER, response) == scale(-ONE, response)
    assert response[1] == ZERO
    assert response[2] == -response[0]
    return {
        "declared_continuum_grading": "P(eta,tau)=C P(-eta,tau), with eta the primal/matching thermal coordinate",
        "differentiated_identity": "partial_eta P|0 = -C partial_eta P|0",
        "ward_compatibility": "A componentwise tau Ward descendant Q4 commutes with C and preserves odd parity.",
        "response_direction": vector_text(response),
        "identities": {
            "delta_P1": fraction_text(response[1]),
            "delta_P2_plus_delta_P0": fraction_text(response[2] + response[0]),
            "A_top_response_per_unit_amplitude": fraction_text(dot(A_TOP, response)),
        },
        "status": "CONDITIONAL_PARITY_THEOREM",
        "condition": (
            "The finite square-site/matching pair is represented by the doubled continuum thermal coordinate eta; "
            "this is not an internal self-duality assertion for the square-site graph alone."
        ),
        "missing_matrix_element": "g[A_top,Q4 epsilon] (parity allows it but does not prove it is nonzero)",
    }


def jordan_parity_certificate() -> dict[str, object]:
    # Basis (q, q_tilde), D q=xq and D q_tilde=x q_tilde+q.
    x = Fraction(21, 4)
    dilation: Matrix = ((x, ONE), (ZERO, x))
    odd_involution: Matrix = ((-ONE, ZERO), (ZERO, -ONE))
    identity2: Matrix = ((ONE, ZERO), (ZERO, ONE))
    assert matmul(odd_involution, dilation) == matmul(dilation, odd_involution)
    assert matmul(odd_involution, odd_involution) == identity2

    return {
        "basis": ["q=Q4 epsilon", "q_tilde=Q4 epsilon_tilde"],
        "dilation": matrix_text(dilation),
        "derived_involution": matrix_text(odd_involution),
        "derivation": [
            "Write J q=-q and J q_tilde=b q+c q_tilde.",
            "[J,D]=0 forces c=-1.",
            "J^2=1 then forces b=0.",
        ],
        "status": "EXACT_GIVEN_GRADING",
        "conclusion": "If the bottom Q4 thermal state is odd, its rank-2 Jordan partner is necessarily odd as well.",
    }


def artifact_payload() -> dict[str, object]:
    return {
        "schema": "matching-one/p275-double-projector-selection/v1",
        "issue": 275,
        "rank_plane": rank_plane_certificate(),
        "vacuum_kdv": vacuum_kdv_certificate(),
        "thermal_q4": thermal_q4_certificate(),
        "thermal_q4_jordan_pair": jordan_parity_certificate(),
        "imported_potts_gate": {
            "status": "EXACT_FOR_REGULAR_UNLABELLED_ONE_INSERTION",
            "dependency_commits": {
                "issue_257_global_matching_spin4_selection": "9320649",
                "issue_262_potts_projector_tomography": "d006f9c",
            },
            "conclusion": "The regular invariant covector annihilates the [2] four-leg sector.",
            "loopholes": ["charged/twisted insertion", "Q-derivative residue", "singular Q->1 normalization"],
        },
        "selection_staircase": [
            {
                "candidate": "vacuum/KdV spin-4",
                "dimension": "4",
                "L_power": "-2",
                "gate": "Alexander-even",
                "A_top": "EXACT_ZERO",
            },
            {
                "candidate": "V_(2,+/-2) four-leg",
                "dimension": "17/4",
                "L_power": "-9/4",
                "gate": "Potts [2]",
                "A_top": "ZERO for regular unlabelled one-insertion",
            },
            {
                "candidate": "thermal Q4 epsilon/Jordan pair",
                "dimension": "21/4",
                "L_power": "-13/4",
                "gate": "Potts singlet and Alexander-odd",
                "A_top": "ALLOWED; nonzero overlap not proved",
            },
        ],
        "strongest_claim": (
            "The first two lower-dimensional spin-4 candidates are removed by independent exact projectors; "
            "thermal Q4 is the first listed candidate allowed by both, conditional on the doubled thermal grading."
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write the exact certificate as JSON")
    args = parser.parse_args(argv)
    rendered = json.dumps(artifact_payload(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
