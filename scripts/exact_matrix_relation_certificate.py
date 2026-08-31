#!/usr/bin/env python3
"""Verify supplied rational matrix relations without numerical tolerances."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "results/model-certificates/synthetic/m2d-vs-m2j/latest.json"
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/matrix-relations/latest.json"
SCHEMA = "matching-one/exact-matrix-relation-certificate/v1"


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


def add(left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    _require(len(left) == len(right) and all(len(a) == len(b) for a, b in zip(left, right)), "matrix shape mismatch")
    return [[a + b for a, b in zip(left_row, right_row)] for left_row, right_row in zip(left, right)]


def scale(value: Fraction, matrix: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    return [[value * entry for entry in row] for row in matrix]


def multiply(left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    _require(left and right and len(left[0]) == len(right), "matrix product shape mismatch")
    _require(all(len(row) == len(left[0]) for row in left) and all(len(row) == len(right[0]) for row in right), "ragged matrix")
    return [
        [sum((left[row][inner] * right[inner][column] for inner in range(len(right))), Fraction()) for column in range(len(right[0]))]
        for row in range(len(left))
    ]


def matrix_rank(values: Sequence[Sequence[Fraction]]) -> int:
    work = [list(row) for row in values]
    _require(work and work[0] and all(len(row) == len(work[0]) for row in work), "rank requires a rectangular matrix")
    rank = 0
    for column in range(len(work[0])):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        work[rank] = [entry / pivot_value for entry in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [left - factor * right for left, right in zip(work[row], work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


def _zero(dimension: int) -> list[list[Fraction]]:
    return [[Fraction() for _ in range(dimension)] for _ in range(dimension)]


def verify_relations(descriptor: Mapping[str, Any]) -> Mapping[str, Any]:
    _require(set(descriptor) == {"generator", "identity_coefficient", "nilpotent", "declared_nilpotent_rank"}, "descriptor fields drift")
    generator = _matrix(descriptor["generator"])
    nilpotent = _matrix(descriptor["nilpotent"])
    dimension = len(generator)
    _require(all(len(row) == dimension for row in generator), "generator must be square")
    _require(len(nilpotent) == dimension and all(len(row) == dimension for row in nilpotent), "nilpotent dimension mismatch")
    coefficient = Fraction(descriptor["identity_coefficient"])
    zero = _zero(dimension)
    affine_residual = add(generator, scale(Fraction(-1), add(scale(coefficient, identity(dimension)), nilpotent)))
    square = multiply(nilpotent, nilpotent)
    commutator = add(multiply(generator, nilpotent), scale(Fraction(-1), multiply(nilpotent, generator)))
    rank = matrix_rank(nilpotent)
    declared_rank = descriptor["declared_nilpotent_rank"]
    _require(affine_residual == zero, "generator is not the declared identity-plus-nilpotent form")
    _require(square == zero, "declared nilpotent does not square to zero")
    _require(commutator == zero, "generator and nilpotent do not commute")
    _require(isinstance(declared_rank, int) and rank == declared_rank, "nilpotent rank mismatch")
    return {
        "dimension": dimension,
        "affine_relation_residual": [[str(value) for value in row] for row in affine_residual],
        "nilpotent_square": [[str(value) for value in row] for row in square],
        "commutator": [[str(value) for value in row] for row in commutator],
        "nilpotent_rank": rank,
        "status": "exact_matrix_relations_verified",
    }


def build_result(source_path: Path = DEFAULT_SOURCE) -> dict[str, Any]:
    source = json.loads(source_path.read_text(encoding="utf-8"))
    realization = source["m2j_extracted_realization"]
    descriptor = {
        "generator": realization["generator"],
        "identity_coefficient": "1",
        "nilpotent": realization["nilpotent"],
        "declared_nilpotent_rank": realization["nilpotent_rank"],
    }
    verification = verify_relations(descriptor)
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
        "solver_invoked": False,
        "claim_boundary": {
            "included": "exact verification of the supplied rational identity-plus-nilpotent matrix relations, commutator, and rank",
            "excluded": "matrix search, diagonalizability of other generators, multi-generator relations, noisy data, SOS completeness, or physical Jordan identification",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], source_path: Path = DEFAULT_SOURCE) -> Mapping[str, Any]:
    expected = build_result(source_path)
    _require(result == expected, "matrix-relation certificate does not exactly reproduce")
    return {
        "schema": result["schema"],
        "status": "valid_exact_matrix_relation_certificate",
        "dimension": result["verification"]["dimension"],
        "nilpotent_rank": result["verification"]["nilpotent_rank"],
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
