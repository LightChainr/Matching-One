#!/usr/bin/env python3
"""Exact projection-leakage commutator control for Issue 400."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Union


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "projection_leakage_commutator_certificate.json"
SCHEMA = "matching-one/projection-leakage-commutator/v1"
ExactInput = Union[int, str, Fraction]
Matrix = tuple[tuple[Fraction, ...], ...]


def exact_fraction(value: ExactInput, *, field: str) -> Fraction:
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError(f"{field} must be exact; floats and booleans are forbidden")
    try:
        return value if isinstance(value, Fraction) else Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"invalid exact value for {field}") from exc


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}" if value.denominator != 1 else str(value.numerator)


def matrix(source: Sequence[Sequence[ExactInput]], *, field: str = "matrix") -> Matrix:
    if not isinstance(source, Sequence) or isinstance(source, (str, bytes)) or not source:
        raise ValueError(f"{field} must be a nonempty row sequence")
    rows = []
    for row_index, row in enumerate(source):
        if not isinstance(row, Sequence) or isinstance(row, (str, bytes)) or not row:
            raise ValueError(f"{field}[{row_index}] must be a nonempty sequence")
        rows.append(
            tuple(
                exact_fraction(value, field=f"{field}[{row_index}][{column_index}]")
                for column_index, value in enumerate(row)
            )
        )
    if any(len(row) != len(rows[0]) for row in rows):
        raise ValueError(f"{field} must be rectangular")
    return tuple(rows)


def square(source: Sequence[Sequence[ExactInput]], *, field: str = "matrix") -> Matrix:
    result = matrix(source, field=field)
    if len(result) != len(result[0]):
        raise ValueError(f"{field} must be square")
    return result


def identity(size: int) -> Matrix:
    if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
        raise ValueError("identity size must be a positive integer")
    return tuple(
        tuple(Fraction(1 if row == column else 0) for column in range(size))
        for row in range(size)
    )


def transpose(left: Matrix) -> Matrix:
    return tuple(tuple(left[row][column] for row in range(len(left))) for column in range(len(left[0])))


def add(left: Matrix, right: Matrix) -> Matrix:
    if len(left) != len(right) or len(left[0]) != len(right[0]):
        raise ValueError("matrix addition requires equal shapes")
    return tuple(
        tuple(left[row][column] + right[row][column] for column in range(len(left[0])))
        for row in range(len(left))
    )


def subtract(left: Matrix, right: Matrix) -> Matrix:
    if len(left) != len(right) or len(left[0]) != len(right[0]):
        raise ValueError("matrix subtraction requires equal shapes")
    return tuple(
        tuple(left[row][column] - right[row][column] for column in range(len(left[0])))
        for row in range(len(left))
    )


def multiply(left: Matrix, right: Matrix) -> Matrix:
    if len(left[0]) != len(right):
        raise ValueError("matrix multiplication dimensions do not close")
    return tuple(
        tuple(
            sum(left[row][inner] * right[inner][column] for inner in range(len(right)))
            for column in range(len(right[0]))
        )
        for row in range(len(left))
    )


def product(*factors: Matrix) -> Matrix:
    if not factors:
        raise ValueError("matrix product requires at least one factor")
    result = factors[0]
    for factor in factors[1:]:
        result = multiply(result, factor)
    return result


def commutator(left: Matrix, right: Matrix) -> Matrix:
    return subtract(multiply(left, right), multiply(right, left))


def frobenius_norm_squared(value: Matrix) -> Fraction:
    return sum(entry * entry for row in value for entry in row)


def matrix_text(value: Matrix) -> list[list[str]]:
    return [[fraction_text(entry) for entry in row] for row in value]


def witness_matrices() -> tuple[Matrix, Matrix, Matrix]:
    u = square(((1, 0, 0), (0, -1, 0), (0, 0, 1)), field="U")
    v = square(((1, 0, 0), (0, 1, 0), (0, 0, -1)), field="V")
    p = square(
        (
            ("2/3", "1/3", "-1/3"),
            ("1/3", "2/3", "1/3"),
            ("-1/3", "1/3", "2/3"),
        ),
        field="P",
    )
    return u, v, p


def analyze(u_source: Sequence[Sequence[ExactInput]], v_source: Sequence[Sequence[ExactInput]], p_source: Sequence[Sequence[ExactInput]]) -> dict[str, Any]:
    u, v, p = square(u_source, field="U"), square(v_source, field="V"), square(p_source, field="P")
    if not (len(u) == len(v) == len(p)):
        raise ValueError("U, V, and P must have the same dimension")
    size = len(u)
    i = identity(size)
    q = subtract(i, p)
    if multiply(p, p) != p or transpose(p) != p:
        raise ValueError("P must be an orthogonal projection")

    full_commutator = commutator(u, v)
    compressed_u = product(p, u, p)
    compressed_v = product(p, v, p)
    compressed_commutator = commutator(compressed_u, compressed_v)
    intrinsic = product(p, full_commutator, p)
    leakage_vu = product(p, v, q, u, p)
    leakage_uv = product(p, u, q, v, p)
    leakage = subtract(leakage_vu, leakage_uv)
    reconstructed = add(intrinsic, leakage)
    if compressed_commutator != reconstructed:
        raise ArithmeticError("projection-leakage commutator identity failed")

    return {
        "dimension": size,
        "U": matrix_text(u),
        "V": matrix_text(v),
        "P": matrix_text(p),
        "Q": matrix_text(q),
        "full_commutator": matrix_text(full_commutator),
        "compressed_commutator": matrix_text(compressed_commutator),
        "intrinsic_projected_commutator": matrix_text(intrinsic),
        "leakage_PVQUP_minus_PUQVP": matrix_text(leakage),
        "compressed_commutator_frobenius_norm_squared": fraction_text(
            frobenius_norm_squared(compressed_commutator)
        ),
        "exact_checks": {
            "P_idempotent": True,
            "P_symmetric": True,
            "U_involution": multiply(u, u) == i,
            "V_involution": multiply(v, v) == i,
            "U_orthogonal": product(transpose(u), u) == i,
            "V_orthogonal": product(transpose(v), v) == i,
            "full_UV_commutator_zero": full_commutator == subtract(i, i),
            "compressed_equals_intrinsic_plus_leakage": True,
        },
    }


def build_artifact() -> dict[str, Any]:
    u, v, p = witness_matrices()
    witness = analyze(u, v, p)
    expected = [
        ["0", "4/9", "4/9"],
        ["-4/9", "0", "4/9"],
        ["-4/9", "-4/9", "0"],
    ]
    if witness["compressed_commutator"] != expected:
        raise AssertionError("frozen compressed commutator changed")
    if witness["compressed_commutator_frobenius_norm_squared"] != "32/27":
        raise AssertionError("frozen commutator norm changed")
    if not all(witness["exact_checks"].values()):
        raise AssertionError("frozen witness exact check failed")

    commuting_projection = square(((1, 0, 0), (0, 1, 0), (0, 0, 0)), field="P0")
    zero_control = analyze(u, v, commuting_projection)
    if zero_control["compressed_commutator_frobenius_norm_squared"] != "0":
        raise AssertionError("commuting-projection control produced curvature")
    return {
        "schema": SCHEMA,
        "issue": 400,
        "status": "exact_projection_leakage_commutator_control",
        "identity": "[PUP,PVP]=P[U,V]P+PVQUP-PUQVP",
        "witness": witness,
        "commuting_projection_zero_control": zero_control,
        "claim_boundary": {
            "included": "finite exact projection-leakage algebra",
            "excluded": "P250 production data, inferred Q-sector excursions, microscopic noncommutation, or memory claims",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("projection-leakage certificate does not exactly reproduce")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_projection_leakage_commutator_control",
        "dimension": expected["witness"]["dimension"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        artifact = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_artifact(artifact), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
