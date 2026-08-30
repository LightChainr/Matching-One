#!/usr/bin/env python3
"""Compile and verify gauge-free exact Hankel-minor certificates."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "results/model-certificates/synthetic/m2d-vs-m2j/latest.json"
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/hankel-minor/latest.json"
SCHEMA = "matching-one/exact-hankel-minor-certificate/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    size = len(matrix)
    _require(size > 0 and all(len(row) == size for row in matrix), "determinant requires a nonempty square matrix")
    work = [list(row) for row in matrix]
    value = Fraction(1)
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return Fraction()
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            value = -value
        pivot_value = work[column][column]
        value *= pivot_value
        for row in range(column + 1, size):
            scale = work[row][column] / pivot_value
            for entry in range(column + 1, size):
                work[row][entry] -= scale * work[column][entry]
    return value


def hankel_matrix(moments: Sequence[Fraction], rows: int, columns: int) -> list[list[Fraction]]:
    _require(rows >= 1 and columns >= 1, "Hankel shape must be positive")
    _require(len(moments) >= rows + columns - 1, "insufficient moments for Hankel shape")
    return [[moments[row + column] for column in range(columns)] for row in range(rows)]


def compile_minors(matrix: Sequence[Sequence[Fraction]], order: int) -> list[dict[str, Any]]:
    rows = len(matrix)
    columns = len(matrix[0])
    _require(1 <= order <= min(rows, columns), "minor order outside matrix shape")
    records = []
    for row_indices in combinations(range(rows), order):
        for column_indices in combinations(range(columns), order):
            minor = [[matrix[row][column] for column in column_indices] for row in row_indices]
            records.append(
                {
                    "rows": list(row_indices),
                    "columns": list(column_indices),
                    "matrix": [[str(value) for value in row] for row in minor],
                    "determinant": str(determinant(minor)),
                }
            )
    return records


def build_result(source_path: Path = DEFAULT_SOURCE) -> dict[str, Any]:
    source = json.loads(source_path.read_text(encoding="utf-8"))
    moments = [Fraction(value) for value in source["synthetic_input"]["moments"]]
    matrix = hankel_matrix(moments, rows=2, columns=4)
    minors = compile_minors(matrix, order=2)
    nonzero = [record for record in minors if Fraction(record["determinant"]) != 0]
    _require(len(minors) == 6 and len(nonzero) == 6, "frozen Hankel witness count drift")
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact",
        "source": {
            "path": str(source_path.relative_to(ROOT)),
            "sha256": _sha256_file(source_path),
            "dependency_group": source["synthetic_input"]["dependency_group"],
        },
        "moment_sequence": [str(value) for value in moments],
        "hankel_shape": [2, 4],
        "hankel_matrix": [[str(value) for value in row] for row in matrix],
        "minor_order": 2,
        "minors": minors,
        "summary": {
            "minor_count": len(minors),
            "nonzero_minor_count": len(nonzero),
            "certified_rank_lower_bound": 2,
            "status": "rank_at_most_one_exactly_excluded",
        },
        "claim_boundary": {
            "included": "all order-two minors of the declared 2x4 exact Hankel matrix",
            "excluded": "a rank upper bound, flat extension, latent realization, noisy confidence region, or physical state dimension",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], source_path: Path = DEFAULT_SOURCE) -> Mapping[str, Any]:
    expected = build_result(source_path)
    _require(result == expected, "Hankel-minor certificate does not exactly reproduce")
    return {
        "schema": result["schema"],
        "status": "valid_exact_hankel_minor_certificate",
        "minor_count": result["summary"]["minor_count"],
        "nonzero_minor_count": result["summary"]["nonzero_minor_count"],
        "rank_lower_bound": result["summary"]["certified_rank_lower_bound"],
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
