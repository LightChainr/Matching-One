#!/usr/bin/env python3
"""Exact and numerical gates for the two-clock dilation hypothesis.

The Gaussian transfer is first de-rotated.  At commensurate clock times a
common generator implies an exact matrix-power identity.  A mixed-context
commutator is the basis-independent obstruction that separate spectra miss.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple


Matrix = Tuple[Tuple[complex, complex], Tuple[complex, complex]]


def matrix2(rows: Sequence[Sequence[complex]]) -> Matrix:
    if len(rows) != 2 or any(len(row) != 2 for row in rows):
        raise ValueError("two-clock v1 accepts exactly 2x2 matrices")
    return (
        (complex(rows[0][0]), complex(rows[0][1])),
        (complex(rows[1][0]), complex(rows[1][1])),
    )


def identity() -> Matrix:
    return ((1.0 + 0j, 0j), (0j, 1.0 + 0j))


def multiply(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(2)) for j in range(2))
        for i in range(2)
    )  # type: ignore[return-value]


def subtract(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(left[i][j] - right[i][j] for j in range(2)) for i in range(2)
    )  # type: ignore[return-value]


def power(value: Matrix, exponent: int) -> Matrix:
    if exponent < 0:
        raise ValueError("power exponents must be nonnegative")
    answer = identity()
    base = value
    n = exponent
    while n:
        if n & 1:
            answer = multiply(answer, base)
        base = multiply(base, base)
        n >>= 1
    return answer


def determinant(value: Matrix) -> complex:
    return value[0][0] * value[1][1] - value[0][1] * value[1][0]


def trace(value: Matrix) -> complex:
    return value[0][0] + value[1][1]


def inverse(value: Matrix) -> Matrix:
    det = determinant(value)
    if abs(det) == 0:
        raise ValueError("spin rotation is singular")
    return (
        (value[1][1] / det, -value[0][1] / det),
        (-value[1][0] / det, value[0][0] / det),
    )


def frobenius(value: Matrix) -> float:
    return math.sqrt(sum(abs(value[i][j]) ** 2 for i in range(2) for j in range(2)))


def complex_payload(value: complex) -> Dict[str, float]:
    return {"real": float(value.real), "imag": float(value.imag)}


def matrix_payload(value: Matrix) -> List[List[Dict[str, float]]]:
    return [[complex_payload(value[i][j]) for j in range(2)] for i in range(2)]


def classify(value: Matrix, tolerance: float) -> str:
    disc = trace(value) ** 2 - 4 * determinant(value)
    scale = max(1.0, abs(trace(value)) ** 2, abs(determinant(value)))
    if abs(disc) > tolerance * scale:
        return "distinct_real_or_complex_pair"
    scalar = ((trace(value) / 2, 0j), (0j, trace(value) / 2))
    if frobenius(subtract(value, scalar)) <= tolerance * max(1.0, frobenius(value)):
        return "scalar"
    return "rank2_Jordan_candidate"


def score(payload: Dict[str, Any]) -> Dict[str, Any]:
    gaussian = matrix2(payload["gaussian_transfer"])
    annulus = matrix2(payload["annulus_transfer"])
    rotation = matrix2(payload.get("spin_rotation", [[1, 0], [0, 1]]))
    derotated = multiply(inverse(rotation), gaussian)
    relation = payload["commensurate_relation"]
    gp = int(relation["gaussian_power"])
    ap = int(relation["annulus_power"])
    tolerance = float(payload.get("tolerance", 1e-10))

    commutator = subtract(multiply(derotated, annulus), multiply(annulus, derotated))
    lhs = power(derotated, gp)
    rhs = power(annulus, ap)
    power_gap = subtract(lhs, rhs)
    comm_rel = frobenius(commutator) / max(1.0, frobenius(derotated) * frobenius(annulus))
    power_rel = frobenius(power_gap) / max(1.0, frobenius(lhs), frobenius(rhs))

    return {
        "schema": "matching-one.two-clock-generator-score.v1",
        "derotated_gaussian_transfer": matrix_payload(derotated),
        "annulus_transfer": matrix_payload(annulus),
        "similarity_invariants": {
            "gaussian_trace": complex_payload(trace(derotated)),
            "gaussian_determinant": complex_payload(determinant(derotated)),
            "gaussian_class": classify(derotated, tolerance),
            "annulus_trace": complex_payload(trace(annulus)),
            "annulus_determinant": complex_payload(determinant(annulus)),
            "annulus_class": classify(annulus, tolerance),
        },
        "mixed_context_gate": {
            "commutator": matrix_payload(commutator),
            "relative_frobenius": comm_rel,
            "passes": comm_rel <= tolerance,
        },
        "commensurate_time_gate": {
            "relation": {"gaussian_power": gp, "annulus_power": ap},
            "relative_frobenius": power_rel,
            "passes": power_rel <= tolerance,
        },
        "verdict": (
            "compatible_with_one_generator_at_declared_times"
            if comm_rel <= tolerance and power_rel <= tolerance
            else "one_generator_obstructed_in_declared_realization"
        ),
        "scope": (
            "Exact equality of fitted transfer matrices is a structural gate, not a "
            "field identification. Separate spectra cannot replace mixed contexts."
        ),
    }


def _fraction_matrix(rows: Iterable[Iterable[int]]) -> Tuple[Tuple[Fraction, ...], ...]:
    return tuple(tuple(Fraction(value) for value in row) for row in rows)


def _fraction_multiply(left: Tuple[Tuple[Fraction, ...], ...], right: Tuple[Tuple[Fraction, ...], ...]):
    return tuple(
        tuple(sum((left[i][k] * right[k][j] for k in range(2)), Fraction(0)) for j in range(2))
        for i in range(2)
    )


def _fraction_subtract(left, right):
    return tuple(tuple(left[i][j] - right[i][j] for j in range(2)) for i in range(2))


def exact_oracles() -> Dict[str, Any]:
    cases = {
        "ordinary": (_fraction_matrix([[2, 0], [0, 3]]), _fraction_matrix([[4, 0], [0, 9]])),
        "Jordan": (_fraction_matrix([[2, 1], [0, 2]]), _fraction_matrix([[4, 4], [0, 4]])),
        "complex_pair": (_fraction_matrix([[3, -4], [4, 3]]), _fraction_matrix([[-7, -24], [24, -7]])),
        "same_Jordan_spectrum_wrong_alignment": (
            _fraction_matrix([[2, 1], [0, 2]]),
            _fraction_matrix([[4, 0], [1, 4]]),
        ),
    }
    rendered: Dict[str, Any] = {}
    for name, (a, u) in cases.items():
        a2 = _fraction_multiply(a, a)
        comm = _fraction_subtract(_fraction_multiply(a, u), _fraction_multiply(u, a))
        rendered[name] = {
            "A_squared_equals_U": a2 == u,
            "commutator_zero": comm == _fraction_matrix([[0, 0], [0, 0]]),
            "A_squared_minus_U": [[str(x) for x in row] for row in _fraction_subtract(a2, u)],
            "commutator": [[str(x) for x in row] for row in comm],
        }
    return {
        "schema": "matching-one.two-clock-exact-oracles.v1",
        "theorem": {
            "necessary": "A=exp(-t_A G), U=exp(-t_U G) implies [A,U]=0.",
            "commensurate": "If p*t_A=q*t_U then A^p=U^q after spin removal.",
            "observable_rectangle": "Reachable/observable mixed contexts turn c^T(AU-UA)b=0 for spanning b,c into [A,U]=0.",
            "nonidentifiability": "Matching separate spectra or Jordan class does not align the two latent bases.",
        },
        "cases": rendered,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--exact-oracles", action="store_true")
    args = parser.parse_args()
    if args.exact_oracles == (args.input is not None):
        parser.error("choose exactly one of --exact-oracles or --input")
    rendered = exact_oracles() if args.exact_oracles else score(json.loads(args.input.read_text()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(rendered, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
