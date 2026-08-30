#!/usr/bin/env python3
"""Exact two-sided Pell defect arithmetic around the hexagonal modulus."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Tuple


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "two_sided_hexagonal_pell_contract.json"


@dataclass(frozen=True)
class Quadratic:
    """An exact element a + b*sqrt(3)."""

    rational: Fraction
    sqrt3: Fraction = Fraction(0)

    def __add__(self, other: "Quadratic") -> "Quadratic":
        return Quadratic(self.rational + other.rational, self.sqrt3 + other.sqrt3)

    def __sub__(self, other: "Quadratic") -> "Quadratic":
        return Quadratic(self.rational - other.rational, self.sqrt3 - other.sqrt3)

    def __mul__(self, other: "Quadratic") -> "Quadratic":
        return Quadratic(
            self.rational * other.rational + 3 * self.sqrt3 * other.sqrt3,
            self.rational * other.sqrt3 + self.sqrt3 * other.rational,
        )


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def pell_family(residual: int, count: int) -> List[Tuple[int, int]]:
    if count < 0:
        raise ValueError("count must be nonnegative")
    if residual == 1:
        p, q = 2, 1
    elif residual == -2:
        p, q = 1, 1
    else:
        raise ValueError("supported Pell residuals are +1 and -2")
    rows: List[Tuple[int, int]] = []
    for _ in range(count):
        rows.append((p, q))
        p, q = 2 * p + 3 * q, p + 2 * q
    return rows


def pell_residual(p: int, q: int) -> int:
    return p * p - 3 * q * q


def site_count(p: int, q: int) -> int:
    return 2 * p * q


def scaled_shape_defect(p: int, q: int) -> Quadratic:
    """Return N*(p/(2q)-sqrt(3)/2), where N=2pq."""

    return Quadratic(Fraction(p * p), Fraction(-p * q))


def exact_limit_error_identity(p: int, q: int) -> bool:
    """Verify (V-eta/2)*(p+q*sqrt(3))^2 = eta^2/2 exactly."""

    eta = pell_residual(p, q)
    if eta not in (1, -2):
        raise ValueError("pair is not in either supported Pell family")
    value = scaled_shape_defect(p, q)
    error = value - Quadratic(Fraction(eta, 2))
    denominator = Quadratic(Fraction(p), Fraction(q))
    return error * denominator * denominator == Quadratic(Fraction(eta * eta, 2))


def row(p: int, q: int) -> Dict[str, object]:
    eta = pell_residual(p, q)
    if eta not in (1, -2):
        raise ValueError("pair is not in either supported Pell family")
    value = scaled_shape_defect(p, q)
    if not exact_limit_error_identity(p, q):
        raise AssertionError("exact Pell limit-error identity failed")
    return {
        "p": p,
        "q": q,
        "pell_residual": eta,
        "site_count": site_count(p, q),
        "scaled_shape_defect": {
            "rational": fraction_text(value.rational),
            "sqrt3": fraction_text(value.sqrt3),
        },
        "shape_side": "above" if eta > 0 else "below",
        "limit": fraction_text(Fraction(eta, 2)),
        "approaches_limit_from": "above",
    }


def build_contract(count: int = 4) -> Dict[str, object]:
    positive = [row(p, q) for p, q in pell_family(1, count)]
    negative = [row(p, q) for p, q in pell_family(-2, count)]
    return {
        "schema": "matching-one/two-sided-hexagonal-pell/v1",
        "status": "valid_exact_two_sided_pell_defect",
        "parent_issue": "remain open",
        "recurrence": "(p,q)->(2p+3q,p+2q)",
        "site_count_identity": "N=2pq",
        "scaled_defect_identity": "N*delta=p^2-pq*sqrt(3)",
        "positive_family": positive,
        "negative_family": negative,
        "positive_limit": "1/2",
        "negative_limit": "-1",
        "negative_over_positive_limit": "-2",
        "makes_e4_or_production_claim": False,
    }


def validate_contract(path: Path = CONTRACT) -> Dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        frozen = json.load(handle)
    actual = build_contract()
    if frozen != actual:
        raise AssertionError("checked-in two-sided Pell contract drifted")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    print(json.dumps(validate_contract(args.contract), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
