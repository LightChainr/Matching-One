#!/usr/bin/env python3
"""Exact synthetic witness separating scalar M1 from diagonal rank-two M2d."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/model-certificates/synthetic/m1-vs-m2d/latest.json"
SCHEMA = "matching-one.synthetic-model-certificate.m1-vs-m2d.v1"


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _matvec(matrix: Sequence[Sequence[Fraction]], vector: Sequence[Fraction]) -> list[Fraction]:
    return [sum((entry * value for entry, value in zip(row, vector)), Fraction()) for row in matrix]


def _response_moments(
    matrix: Sequence[Sequence[Fraction]],
    source: Sequence[Fraction],
    readout: Sequence[Fraction],
    count: int,
) -> list[Fraction]:
    state = list(source)
    moments = []
    for _ in range(count):
        moments.append(sum((left * right for left, right in zip(readout, state)), Fraction()))
        state = _matvec(matrix, state)
    return moments


def build_result() -> dict[str, Any]:
    moments = [Fraction(2), Fraction(3), Fraction(5)]
    rank_one_minor = moments[0] * moments[2] - moments[1] ** 2
    matrix = [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(2)]]
    source = [Fraction(1), Fraction(1)]
    readout = [Fraction(1), Fraction(1)]
    reproduced = _response_moments(matrix, source, readout, len(moments))
    if rank_one_minor == 0 or reproduced != moments:
        raise AssertionError("synthetic M1/M2d witness construction failed")

    synthetic_input = {
        "word_family": ["identity", "a", "a^2"],
        "moments": [str(value) for value in moments],
        "dependency_group": "synthetic-exact-m1-m2d",
    }
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact_synthetic_control",
        "synthetic_input": synthetic_input,
        "synthetic_input_sha256": _sha256(synthetic_input),
        "model_classes": {
            "excluded": "M1: one scalar character f_n=alpha*lambda^n",
            "feasible": "M2d: one diagonalizable two-state generator",
        },
        "m1_infeasibility_certificate": {
            "identity": "f_0*f_2-f_1^2=0 for every scalar character",
            "evaluated_minor": str(rank_one_minor),
            "primitive_witness": [1, -1],
            "status": "exactly_infeasible",
        },
        "m2d_extracted_realization": {
            "generator": [[str(value) for value in row] for row in matrix],
            "source": [str(value) for value in source],
            "readout": [str(value) for value in readout],
            "reproduced_moments": [str(value) for value in reproduced],
            "eigenvalues": ["1", "2"],
            "status": "exactly_feasible",
        },
        "claim_boundary": {
            "included": "one exact one-generator rational positive/negative control",
            "excluded": "SOS completeness, noisy data, multi-generator relations, physical field identification, or exclusion of M1 outside the declared three moments",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = build_result()
    if result != expected:
        raise ValueError("M1/M2d certificate does not exactly reproduce")
    if result["claim_boundary"]["parent_issue"] != "remain open":
        raise ValueError("parent issue boundary drift")
    return {
        "schema": result["schema"],
        "status": "valid_exact_synthetic_control",
        "m1_minor": result["m1_infeasibility_certificate"]["evaluated_minor"],
        "m2d_status": result["m2d_extracted_realization"]["status"],
        "synthetic_input_sha256": result["synthetic_input_sha256"],
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
