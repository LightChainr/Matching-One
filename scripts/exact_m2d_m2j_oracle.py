#!/usr/bin/env python3
"""Exact synthetic witness separating diagonal rank-two M2d from Jordan M2j."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/model-certificates/synthetic/m2d-vs-m2j/latest.json"
SCHEMA = "matching-one.synthetic-model-certificate.m2d-vs-m2j.v1"


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _matmul(left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    return [
        [sum((left[i][k] * right[k][j] for k in range(len(right))), Fraction()) for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def _matvec(matrix: Sequence[Sequence[Fraction]], vector: Sequence[Fraction]) -> list[Fraction]:
    return [sum((entry * value for entry, value in zip(row, vector)), Fraction()) for row in matrix]


def _moments(
    matrix: Sequence[Sequence[Fraction]],
    source: Sequence[Fraction],
    readout: Sequence[Fraction],
    count: int,
) -> list[Fraction]:
    state = list(source)
    values = []
    for _ in range(count):
        values.append(sum((a * b for a, b in zip(readout, state)), Fraction()))
        state = _matvec(matrix, state)
    return values


def build_result() -> dict[str, Any]:
    moments = [Fraction(n + 1) for n in range(5)]
    hankel_minor = moments[0] * moments[2] - moments[1] ** 2
    recurrence_sum = Fraction(2)
    recurrence_product = Fraction(-1)
    recurrence_residuals = [
        moments[n + 2] - recurrence_sum * moments[n + 1] - recurrence_product * moments[n]
        for n in range(3)
    ]
    discriminant = recurrence_sum**2 + 4 * recurrence_product

    matrix = [[Fraction(1), Fraction(1)], [Fraction(0), Fraction(1)]]
    identity = [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(1)]]
    nilpotent = [[matrix[i][j] - identity[i][j] for j in range(2)] for i in range(2)]
    nilpotent_square = _matmul(nilpotent, nilpotent)
    source = [Fraction(0), Fraction(1)]
    readout = [Fraction(1), Fraction(1)]
    reproduced = _moments(matrix, source, readout, len(moments))
    if hankel_minor == 0 or any(recurrence_residuals) or discriminant != 0:
        raise AssertionError("minimal repeated-root recurrence construction failed")
    if nilpotent == [[Fraction(), Fraction()], [Fraction(), Fraction()]] or any(map(any, nilpotent_square)):
        raise AssertionError("common-nilpotent Jordan witness failed")
    if reproduced != moments:
        raise AssertionError("Jordan realization did not reproduce the synthetic moments")

    synthetic_input = {
        "word_family": [f"a^{n}" for n in range(len(moments))],
        "moments": [str(value) for value in moments],
        "dependency_group": "synthetic-exact-m2d-m2j",
    }
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact_synthetic_control",
        "synthetic_input": synthetic_input,
        "synthetic_input_sha256": _sha256(synthetic_input),
        "declared_m2d_class": "minimal reachable/observable one-generator rank-two realization with a diagonalizable generator",
        "m2d_infeasibility_certificate": {
            "nonzero_rank_two_hankel_minor": str(hankel_minor),
            "unique_recurrence": "f_(n+2)=2*f_(n+1)-f_n",
            "recurrence_residuals": [str(value) for value in recurrence_residuals],
            "characteristic_polynomial": "(x-1)^2",
            "discriminant": str(discriminant),
            "reason": "minimal rank two forces the repeated characteristic root; a diagonalizable two-by-two matrix with only that root is scalar and would have Hankel rank one",
            "status": "exactly_infeasible",
        },
        "m2j_extracted_realization": {
            "generator": [[str(value) for value in row] for row in matrix],
            "nilpotent": [[str(value) for value in row] for row in nilpotent],
            "nilpotent_square": [[str(value) for value in row] for row in nilpotent_square],
            "nilpotent_rank": 1,
            "source": [str(value) for value in source],
            "readout": [str(value) for value in readout],
            "reproduced_moments": [str(value) for value in reproduced],
            "status": "exactly_feasible",
        },
        "claim_boundary": {
            "included": "one exact minimal one-generator rational control with a unique order-two recurrence",
            "excluded": "all product-only M2d variants, multi-generator constraints, noisy data, SOS completeness, or physical Jordan-field identification",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = build_result()
    if result != expected:
        raise ValueError("M2d/M2j certificate does not exactly reproduce")
    return {
        "schema": result["schema"],
        "status": "valid_exact_synthetic_control",
        "hankel_minor": result["m2d_infeasibility_certificate"]["nonzero_rank_two_hankel_minor"],
        "discriminant": result["m2d_infeasibility_certificate"]["discriminant"],
        "m2j_status": result["m2j_extracted_realization"]["status"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_result(json.loads(args.validate.read_text())), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_result(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
