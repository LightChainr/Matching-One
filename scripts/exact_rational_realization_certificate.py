#!/usr/bin/env python3
"""Verify explicit typed rational state-space realizations exactly."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "results/model-certificates/synthetic/morphism-forces-rank3/latest.json"
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/rational-realization/latest.json"
SCHEMA = "matching-one/exact-rational-realization-certificate/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fraction_vector(values: Sequence[Any]) -> list[Fraction]:
    return [Fraction(value) for value in values]


def _fraction_matrix(values: Sequence[Sequence[Any]]) -> list[list[Fraction]]:
    return [_fraction_vector(row) for row in values]


def matvec(matrix: Sequence[Sequence[Fraction]], vector: Sequence[Fraction]) -> list[Fraction]:
    _require(matrix and all(len(row) == len(vector) for row in matrix), "matrix/vector shape mismatch")
    return [sum((entry * value for entry, value in zip(row, vector)), Fraction()) for row in matrix]


def row_times_matrix(row: Sequence[Fraction], matrix: Sequence[Sequence[Fraction]]) -> list[Fraction]:
    _require(matrix and len(row) == len(matrix), "row/matrix shape mismatch")
    return [sum((row[index] * matrix[index][column] for index in range(len(row))), Fraction()) for column in range(len(matrix[0]))]


def matrix_rank(matrix: Sequence[Sequence[Fraction]]) -> int:
    _require(matrix and matrix[0] and all(len(row) == len(matrix[0]) for row in matrix), "rank requires a rectangular matrix")
    work = [[Fraction(value) for value in row] for row in matrix]
    rank = 0
    for column in range(len(work[0])):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        work[rank] = [value / pivot_value for value in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [left - scale * right for left, right in zip(work[row], work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


def verify_realization(descriptor: Mapping[str, Any]) -> Mapping[str, Any]:
    generator = _fraction_matrix(descriptor["generator"])
    dimension = len(generator)
    _require(dimension >= 1 and all(len(row) == dimension for row in generator), "generator must be square")
    channels = descriptor["channels"]
    _require(isinstance(channels, list) and channels, "typed channels are required")
    reproduced = []
    reachability_columns: list[list[Fraction]] = []
    observability_rows: list[list[Fraction]] = []
    for channel in channels:
        _require(set(channel) == {"id", "source", "readout", "moments"}, "channel fields drift")
        source = _fraction_vector(channel["source"])
        readout = _fraction_vector(channel["readout"])
        expected = _fraction_vector(channel["moments"])
        _require(len(source) == dimension and len(readout) == dimension, "channel dimension mismatch")
        state = source
        row = readout
        actual = []
        for _ in range(len(expected)):
            actual.append(sum((left * right for left, right in zip(readout, state)), Fraction()))
            reachability_columns.append(list(state))
            observability_rows.append(list(row))
            state = matvec(generator, state)
            row = row_times_matrix(row, generator)
        _require(actual == expected, f"response mismatch in channel {channel['id']}")
        reproduced.append({"id": channel["id"], "moments": [str(value) for value in actual]})
    reachability_matrix = [
        [column[row] for column in reachability_columns]
        for row in range(dimension)
    ]
    reachability_rank = matrix_rank(reachability_matrix)
    observability_rank = matrix_rank(observability_rows)
    return {
        "dimension": dimension,
        "channels": reproduced,
        "reachability_rank_across_typed_channels": reachability_rank,
        "observability_rank_across_typed_channels": observability_rank,
        "minimal_on_declared_typed_rows": reachability_rank == dimension and observability_rank == dimension,
        "status": "exact_rational_realization_verified",
    }


def build_result(source_path: Path = DEFAULT_SOURCE) -> dict[str, Any]:
    source = json.loads(source_path.read_text(encoding="utf-8"))
    realization = source["rank3_extracted_realization"]
    descriptor = {
        "generator": realization["generator"],
        "channels": [
            {
                "id": "endpoint",
                "source": realization["endpoint_source"],
                "readout": realization["endpoint_readout"],
                "moments": realization["reproduced_endpoint_moments"],
            },
            {
                "id": "morphism",
                "source": realization["morphism_source"],
                "readout": realization["morphism_readout"],
                "moments": realization["reproduced_morphism_moments"],
            },
        ],
    }
    verification = verify_realization(descriptor)
    _require(verification["minimal_on_declared_typed_rows"], "frozen realization lost typed minimality")
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
            "included": "exact response, reachability, and observability verification of the supplied typed rational realization",
            "excluded": "realization search or extraction, uniqueness up to similarity, other contexts, noisy data, cover/deck/Smith semantics, or physical state dimension",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], source_path: Path = DEFAULT_SOURCE) -> Mapping[str, Any]:
    expected = build_result(source_path)
    _require(result == expected, "rational-realization certificate does not exactly reproduce")
    verification = result["verification"]
    return {
        "schema": result["schema"],
        "status": "valid_exact_rational_realization_certificate",
        "dimension": verification["dimension"],
        "channel_count": len(verification["channels"]),
        "reachability_rank": verification["reachability_rank_across_typed_channels"],
        "observability_rank": verification["observability_rank_across_typed_channels"],
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
