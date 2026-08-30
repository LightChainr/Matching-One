#!/usr/bin/env python3
"""Exact reachable-source gauge fixture close to a chart boundary."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/model-certificates/synthetic/gauge-chart-boundary/latest.json"
SCHEMA = "matching-one.synthetic-model-certificate.gauge-boundary.v1"


Matrix = list[list[Fraction]]


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _matmul(left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]) -> Matrix:
    return [
        [sum((left[i][k] * right[k][j] for k in range(len(right))), Fraction()) for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def _identity(size: int) -> Matrix:
    return [[Fraction(int(i == j)) for j in range(size)] for i in range(size)]


def _powers(matrix: Matrix, count: int) -> list[Matrix]:
    values = [_identity(len(matrix))]
    for _ in range(1, count):
        values.append(_matmul(values[-1], matrix))
    return values


def _response_rows(readout: Matrix, matrix: Matrix, sources: Matrix, count: int) -> list[list[Fraction]]:
    return [_matmul(_matmul(readout, power), sources)[0] for power in _powers(matrix, count)]


def _render_matrix(matrix: Sequence[Sequence[Fraction]]) -> list[list[str]]:
    return [[str(value) for value in row] for row in matrix]


def build_result() -> dict[str, Any]:
    epsilon = Fraction(1, 1024)
    generator = [[Fraction(1), Fraction(1)], [Fraction(0), Fraction(2)]]
    sources = [[epsilon, Fraction(0)], [Fraction(0), Fraction(1)]]
    source_inverse = [[1 / epsilon, Fraction(0)], [Fraction(0), Fraction(1)]]
    readout = [[Fraction(1), Fraction(1)]]

    normalized_generator = _matmul(_matmul(source_inverse, generator), sources)
    normalized_sources = _matmul(source_inverse, sources)
    normalized_readout = _matmul(readout, sources)
    original_rows = _response_rows(readout, generator, sources, 4)
    normalized_rows = _response_rows(normalized_readout, normalized_generator, normalized_sources, 4)
    source_minor = sources[0][0] * sources[1][1] - sources[0][1] * sources[1][0]
    trace_original = generator[0][0] + generator[1][1]
    trace_normalized = normalized_generator[0][0] + normalized_generator[1][1]
    determinant_original = generator[0][0] * generator[1][1] - generator[0][1] * generator[1][0]
    determinant_normalized = (
        normalized_generator[0][0] * normalized_generator[1][1]
        - normalized_generator[0][1] * normalized_generator[1][0]
    )
    if normalized_sources != _identity(2) or original_rows != normalized_rows:
        raise AssertionError("reachable-source gauge did not preserve observable contractions")
    if source_minor != epsilon or normalized_generator[0][1] != 1 / epsilon:
        raise AssertionError("near-boundary amplification construction failed")
    if (trace_original, determinant_original) != (trace_normalized, determinant_normalized):
        raise AssertionError("similarity invariants drifted")

    synthetic_input = {
        "epsilon": str(epsilon),
        "generator": _render_matrix(generator),
        "source_matrix": _render_matrix(sources),
        "readout": _render_matrix(readout),
        "dependency_group": "synthetic-exact-gauge-boundary",
    }
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact_synthetic_control",
        "synthetic_input": synthetic_input,
        "synthetic_input_sha256": _sha256(synthetic_input),
        "reachable_source_chart": {
            "source_minor": str(source_minor),
            "chart_open": True,
            "boundary_value": "0",
            "distance_parameter": str(epsilon),
            "normalizing_similarity": _render_matrix(source_inverse),
            "normalized_source_matrix": _render_matrix(normalized_sources),
            "normalized_generator": _render_matrix(normalized_generator),
            "largest_normalized_entry": str(1 / epsilon),
        },
        "exact_invariant_checks": {
            "original_response_rows": _render_matrix(original_rows),
            "normalized_response_rows": _render_matrix(normalized_rows),
            "responses_identical": True,
            "trace_before_after": [str(trace_original), str(trace_normalized)],
            "determinant_before_after": [str(determinant_original), str(determinant_normalized)],
        },
        "boundary_diagnostic": {
            "amplification_factor": str(1 / epsilon),
            "at_epsilon_zero": "reachable-source chart undefined",
            "status": "exact_near_boundary_positive_control",
        },
        "claim_boundary": {
            "included": "one exact nonzero reachable-source chart with explicit amplification and invariant observable checks",
            "excluded": "coverage of the epsilon=0 boundary, a numerically stable chart, all similarity orbits, SOS completeness, or physical model selection",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = build_result()
    if result != expected:
        raise ValueError("gauge-boundary certificate does not exactly reproduce")
    return {
        "schema": result["schema"],
        "status": "valid_exact_synthetic_control",
        "source_minor": result["reachable_source_chart"]["source_minor"],
        "amplification_factor": result["boundary_diagnostic"]["amplification_factor"],
        "responses_identical": result["exact_invariant_checks"]["responses_identical"],
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
