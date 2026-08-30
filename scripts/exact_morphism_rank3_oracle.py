#!/usr/bin/env python3
"""Exact synthetic endpoint/morphism witness forcing three predictive states."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/model-certificates/synthetic/morphism-forces-rank3/latest.json"
SCHEMA = "matching-one.synthetic-model-certificate.morphism-rank3.v1"


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


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


def _det3(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


def build_result() -> dict[str, Any]:
    endpoint = [Fraction(n + 1) for n in range(4)]
    morphism = [Fraction(2**n) for n in range(3)]
    stacked_hankel = [endpoint[:3], endpoint[1:4], morphism]
    determinant = _det3(stacked_hankel)

    endpoint_matrix = [[Fraction(1), Fraction(1)], [Fraction(0), Fraction(1)]]
    endpoint_source = [Fraction(0), Fraction(1)]
    endpoint_readout = [Fraction(1), Fraction(1)]
    endpoint_reproduced = _moments(endpoint_matrix, endpoint_source, endpoint_readout, 4)

    rank3_matrix = [
        [Fraction(1), Fraction(1), Fraction(0)],
        [Fraction(0), Fraction(1), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(2)],
    ]
    rank3_endpoint_source = [Fraction(0), Fraction(1), Fraction(0)]
    rank3_endpoint_readout = [Fraction(1), Fraction(1), Fraction(0)]
    rank3_morphism_source = [Fraction(0), Fraction(0), Fraction(1)]
    rank3_morphism_readout = [Fraction(0), Fraction(0), Fraction(1)]
    rank3_endpoint = _moments(rank3_matrix, rank3_endpoint_source, rank3_endpoint_readout, 4)
    rank3_morphism = _moments(rank3_matrix, rank3_morphism_source, rank3_morphism_readout, 3)
    if endpoint_reproduced != endpoint or determinant == 0:
        raise AssertionError("endpoint/rank-lower-bound construction failed")
    if rank3_endpoint != endpoint or rank3_morphism != morphism:
        raise AssertionError("three-state realization did not reproduce both typed rows")

    synthetic_input = {
        "endpoint_moments": [str(value) for value in endpoint],
        "morphism_sensitive_moments": [str(value) for value in morphism],
        "row_types": ["endpoint", "endpoint_shift", "morphism"],
        "dependency_group": "synthetic-exact-morphism-rank3",
    }
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact_synthetic_control",
        "synthetic_input": synthetic_input,
        "synthetic_input_sha256": _sha256(synthetic_input),
        "endpoint_m2j_positive_control": {
            "generator": [[str(value) for value in row] for row in endpoint_matrix],
            "reproduced_moments": [str(value) for value in endpoint_reproduced],
            "status": "exactly_feasible",
        },
        "common_rank_two_infeasibility_certificate": {
            "stacked_typed_hankel": [[str(value) for value in row] for row in stacked_hankel],
            "determinant": str(determinant),
            "reason": "every common r-state realization factors the stacked matrix through an r-dimensional state space",
            "certified_minimum_predictive_rank": 3,
            "status": "exactly_infeasible",
        },
        "rank3_extracted_realization": {
            "generator": [[str(value) for value in row] for row in rank3_matrix],
            "endpoint_source": [str(value) for value in rank3_endpoint_source],
            "endpoint_readout": [str(value) for value in rank3_endpoint_readout],
            "morphism_source": [str(value) for value in rank3_morphism_source],
            "morphism_readout": [str(value) for value in rank3_morphism_readout],
            "reproduced_endpoint_moments": [str(value) for value in rank3_endpoint],
            "reproduced_morphism_moments": [str(value) for value in rank3_morphism],
            "status": "exactly_feasible",
        },
        "claim_boundary": {
            "included": "one synthetic typed-row rank obstruction and one explicit three-state realization",
            "excluded": "a physical cover/deck/Smith realization, general morphism completeness, noisy confidence regions, SOS completeness, or field identification",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = build_result()
    if result != expected:
        raise ValueError("morphism rank-three certificate does not exactly reproduce")
    return {
        "schema": result["schema"],
        "status": "valid_exact_synthetic_control",
        "stacked_determinant": result["common_rank_two_infeasibility_certificate"]["determinant"],
        "minimum_rank": result["common_rank_two_infeasibility_certificate"]["certified_minimum_predictive_rank"],
        "rank3_status": result["rank3_extracted_realization"]["status"],
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
