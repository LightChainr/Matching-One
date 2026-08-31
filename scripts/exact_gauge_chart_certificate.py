#!/usr/bin/env python3
"""Verify a supplied reachable-source gauge chart and its boundary exactly."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "results/model-certificates/synthetic/gauge-chart-boundary/latest.json"
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/gauge-chart/latest.json"
SCHEMA = "matching-one/exact-gauge-chart-certificate/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix(values: Sequence[Sequence[Any]]) -> list[list[Fraction]]:
    matrix = [[Fraction(value) for value in row] for row in values]
    _require(matrix and matrix[0] and all(len(row) == len(matrix[0]) for row in matrix), "matrix must be rectangular")
    return matrix


def identity(dimension: int) -> list[list[Fraction]]:
    return [[Fraction(row == column) for column in range(dimension)] for row in range(dimension)]


def multiply(left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    _require(left and right and len(left[0]) == len(right), "matrix product shape mismatch")
    _require(all(len(row) == len(left[0]) for row in left) and all(len(row) == len(right[0]) for row in right), "ragged matrix")
    return [
        [sum((left[row][inner] * right[inner][column] for inner in range(len(right))), Fraction()) for column in range(len(right[0]))]
        for row in range(len(left))
    ]


def inverse(matrix: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    dimension = len(matrix)
    _require(dimension and all(len(row) == dimension for row in matrix), "inverse requires a square matrix")
    work = [list(row) + ident for row, ident in zip(matrix, identity(dimension))]
    for column in range(dimension):
        pivot = next((row for row in range(column, dimension) if work[row][column]), None)
        _require(pivot is not None, "source chart minor is singular")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [value / pivot_value for value in work[column]]
        for row in range(dimension):
            if row == column or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [left - factor * right for left, right in zip(work[row], work[column])]
    return [row[dimension:] for row in work]


def determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    dimension = len(matrix)
    _require(dimension and all(len(row) == dimension for row in matrix), "determinant requires a square matrix")
    work = [list(row) for row in matrix]
    result = Fraction(1)
    for column in range(dimension):
        pivot = next((row for row in range(column, dimension) if work[row][column]), None)
        if pivot is None:
            return Fraction()
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        pivot_value = work[column][column]
        result *= pivot_value
        for row in range(column + 1, dimension):
            factor = work[row][column] / pivot_value
            for inner in range(column, dimension):
                work[row][inner] -= factor * work[column][inner]
    return result


def power(matrix: Sequence[Sequence[Fraction]], exponent: int) -> list[list[Fraction]]:
    _require(exponent >= 0 and len(matrix) == len(matrix[0]), "invalid matrix power")
    result = identity(len(matrix))
    base = [list(row) for row in matrix]
    while exponent:
        if exponent & 1:
            result = multiply(result, base)
        base = multiply(base, base)
        exponent //= 2
    return result


def response_rows(readout: Sequence[Sequence[Fraction]], generator: Sequence[Sequence[Fraction]], source: Sequence[Sequence[Fraction]], count: int) -> list[list[Fraction]]:
    return [multiply(multiply(readout, power(generator, exponent)), source)[0] for exponent in range(count)]


def verify_chart(descriptor: Mapping[str, Any]) -> Mapping[str, Any]:
    required = {"generator", "source_matrix", "readout", "normalizing_similarity", "response_count", "boundary_minor_value"}
    _require(set(descriptor) == required, "descriptor fields drift")
    generator = _matrix(descriptor["generator"])
    source = _matrix(descriptor["source_matrix"])
    readout = _matrix(descriptor["readout"])
    similarity = _matrix(descriptor["normalizing_similarity"])
    dimension = len(generator)
    _require(all(len(row) == dimension for row in generator), "generator must be square")
    _require(len(source) == dimension and all(len(row) == dimension for row in source), "source matrix dimension mismatch")
    _require(len(readout) == 1 and len(readout[0]) == dimension, "readout dimension mismatch")
    _require(len(similarity) == dimension and all(len(row) == dimension for row in similarity), "similarity dimension mismatch")
    source_minor = determinant(source)
    _require(source_minor != 0, "source chart minor is singular")
    _require(similarity == inverse(source), "normalizing similarity is not the exact source inverse")
    similarity_inverse = inverse(similarity)
    normalized_source = multiply(similarity, source)
    normalized_generator = multiply(multiply(similarity, generator), similarity_inverse)
    normalized_readout = multiply(readout, similarity_inverse)
    count = descriptor["response_count"]
    _require(isinstance(count, int) and count >= 1, "response count must be positive")
    original_responses = response_rows(readout, generator, source, count)
    normalized_responses = response_rows(normalized_readout, normalized_generator, normalized_source, count)
    _require(normalized_source == identity(dimension), "source matrix did not normalize to identity")
    _require(original_responses == normalized_responses, "gauge transformation changed responses")
    boundary = Fraction(descriptor["boundary_minor_value"])
    _require(boundary == 0, "reachable-source boundary must be the zero-minor locus")
    amplification = max(abs(value) for row in similarity for value in row)
    trace_before = sum((generator[index][index] for index in range(dimension)), Fraction())
    trace_after = sum((normalized_generator[index][index] for index in range(dimension)), Fraction())
    determinant_before = determinant(generator)
    determinant_after = determinant(normalized_generator)
    _require(trace_before == trace_after and determinant_before == determinant_after, "similarity invariants changed")
    return {
        "dimension": dimension,
        "source_minor": str(source_minor),
        "normalized_source": [[str(value) for value in row] for row in normalized_source],
        "normalized_generator": [[str(value) for value in row] for row in normalized_generator],
        "original_response_rows": [[str(value) for value in row] for row in original_responses],
        "normalized_response_rows": [[str(value) for value in row] for row in normalized_responses],
        "responses_identical": True,
        "trace_before_after": [str(trace_before), str(trace_after)],
        "determinant_before_after": [str(determinant_before), str(determinant_after)],
        "amplification_factor": str(amplification),
        "chart_coverage": "nonzero_source_minor_only",
        "uncovered_boundary": "source_minor=0",
        "status": "exact_reachable_source_chart_verified",
    }


def build_result(source_path: Path = DEFAULT_SOURCE) -> dict[str, Any]:
    source = json.loads(source_path.read_text(encoding="utf-8"))
    descriptor = {
        "generator": source["synthetic_input"]["generator"],
        "source_matrix": source["synthetic_input"]["source_matrix"],
        "readout": source["synthetic_input"]["readout"],
        "normalizing_similarity": source["reachable_source_chart"]["normalizing_similarity"],
        "response_count": len(source["exact_invariant_checks"]["original_response_rows"]),
        "boundary_minor_value": source["reachable_source_chart"]["boundary_value"],
    }
    verification = verify_chart(descriptor)
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact",
        "source": {
            "path": str(source_path.relative_to(ROOT)),
            "sha256": _sha256_file(source_path),
            "dependency_group": source["synthetic_input"]["dependency_group"],
        },
        "descriptor": descriptor,
        "verification": verification,
        "claim_boundary": {
            "included": "exact verification of one supplied reachable-source chart, response invariance, amplification, and its zero-minor boundary",
            "excluded": "coverage of the zero-minor boundary, a complete gauge atlas, numerical conditioning guarantees, model search, noisy data, or physical interpretation",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], source_path: Path = DEFAULT_SOURCE) -> Mapping[str, Any]:
    expected = build_result(source_path)
    _require(result == expected, "gauge-chart certificate does not exactly reproduce")
    verification = result["verification"]
    return {
        "schema": result["schema"],
        "status": "valid_exact_gauge_chart_certificate",
        "source_minor": verification["source_minor"],
        "amplification_factor": verification["amplification_factor"],
        "coverage": verification["chart_coverage"],
        "source_sha256": result["source"]["sha256"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        value = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_result(value, args.source), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_result(args.source), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
