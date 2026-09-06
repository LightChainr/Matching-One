#!/usr/bin/env python3
"""Exact angular/Smith preflight for Issues #589/#590.

No percolation target is read.  The script only uses Gaussian-integer arithmetic
to answer three design questions:

1. what Smith class does each N=650 orientation realize?
2. how well conditioned is the exact H0/H4/H8 interpolation?
3. what is the first shared site count for which both the square and 2:1
   rectangular families have at least three primitive/cyclic directions?
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path

SCHEMA = "matching-one/n650-angular-smith-preflight/v1"


def d4_representations(norm: int, primitive: bool | None = None) -> list[tuple[int, int]]:
    """Gaussian representations a^2+b^2=norm, one representative modulo D4."""
    out: list[tuple[int, int]] = []
    for b in range(isqrt(norm) + 1):
        a2 = norm - b * b
        a = isqrt(a2)
        if a * a != a2 or a < b:
            continue
        is_primitive = gcd(a, b) == 1
        if primitive is None or is_primitive == primitive:
            out.append((a, b))
    return out


def smith_square(a: int, b: int) -> tuple[int, int]:
    n = a * a + b * b
    d1 = gcd(a, b)
    return d1, n // d1


def smith_rectangle_2to1(a: int, b: int) -> tuple[int, int]:
    # Period matrix columns (a,b) and (-2b,2a).
    n = 2 * (a * a + b * b)
    d1 = gcd(a, b)
    return d1, n // d1


def cos4(a: int, b: int) -> Fraction:
    n = a * a + b * b
    return Fraction(a**4 - 6 * a * a * b * b + b**4, n * n)


def cos8(a: int, b: int) -> Fraction:
    c4 = cos4(a, b)
    return 2 * c4 * c4 - 1


def invert3(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    a = [
        list(map(Fraction, row))
        + [Fraction(int(i == j)) for j in range(3)]
        for i, row in enumerate(matrix)
    ]
    for col in range(3):
        pivot = next(row for row in range(col, 3) if a[row][col])
        a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        a[col] = [value / scale for value in a[col]]
        for row in range(3):
            if row == col:
                continue
            scale = a[row][col]
            if scale:
                a[row] = [x - scale * y for x, y in zip(a[row], a[col])]
    return [row[3:] for row in a]


def harmonic_design(reps: tuple[tuple[int, int], ...]) -> dict[str, object]:
    if len(reps) != 3:
        raise ValueError("H0/H4/H8 design requires exactly three rows")
    matrix = [[Fraction(1), cos4(a, b), cos8(a, b)] for a, b in reps]
    inverse = invert3(matrix)
    variance_amplification = [sum(weight * weight for weight in row) for row in inverse]
    return {
        "representations": [list(rep) for rep in reps],
        "matrix": [[str(value) for value in row] for row in matrix],
        "inverse_weights_for_C_A4_A8": [
            [str(value) for value in row] for row in inverse
        ],
        "equal_unit_variance_amplification_C_A4_A8": [
            float(value) for value in variance_amplification
        ],
        "worst_amplification": float(max(variance_amplification)),
    }


def best_triplet(reps: list[tuple[int, int]]) -> dict[str, object]:
    best = None
    for triple in combinations(reps, 3):
        design = harmonic_design(tuple(triple))
        score = design["worst_amplification"]
        if best is None or score < best[0]:
            best = (score, triple, design)
    if best is None:
        raise ValueError("fewer than three representations")
    return best[2]


def first_common_site_count(limit: int, primitive: bool | None) -> dict[str, object]:
    for site_count in range(2, limit + 1, 2):
        square = d4_representations(site_count, primitive)
        rectangle = d4_representations(site_count // 2, primitive)
        if len(square) >= 3 and len(rectangle) >= 3:
            return {
                "site_count": site_count,
                "square_representations": [list(rep) for rep in square],
                "rectangle_underlying_representations": [list(rep) for rep in rectangle],
            }
    raise ValueError(f"no common design through {limit}")


def row_payload(rep: tuple[int, int], family: str) -> dict[str, object]:
    a, b = rep
    smith = smith_square(a, b) if family == "square" else smith_rectangle_2to1(a, b)
    return {
        "representation": [a, b],
        "gcd": gcd(a, b),
        "smith": list(smith),
        "cyclic": smith[0] == 1,
        "cos4": str(cos4(a, b)),
        "cos8": str(cos8(a, b)),
    }


def build() -> dict[str, object]:
    square_650 = ((25, 5), (23, 11), (19, 17))
    rectangle_650 = ((18, 1), (17, 6), (15, 10))

    first_any = first_common_site_count(650, primitive=None)
    first_cyclic = first_common_site_count(2210, primitive=True)

    assert first_any["site_count"] == 650
    assert first_cyclic["site_count"] == 2210
    assert [smith_square(*rep) for rep in square_650] == [(5, 130), (1, 650), (1, 650)]
    assert [smith_rectangle_2to1(*rep) for rep in rectangle_650] == [
        (1, 650), (1, 650), (5, 130)
    ]

    square_2210 = [tuple(rep) for rep in first_cyclic["square_representations"]]
    rectangle_2210 = [
        tuple(rep) for rep in first_cyclic["rectangle_underlying_representations"]
    ]

    return {
        "schema": SCHEMA,
        "issues": [589, 590],
        "evidence_type": "exact geometry/arithmetic only; no percolation target read",
        "n650": {
            "square": [row_payload(rep, "square") for rep in square_650],
            "rectangle_2to1": [
                row_payload(rep, "rectangle") for rep in rectangle_650
            ],
            "square_h0_h4_h8_design": harmonic_design(square_650),
            "rectangle_h0_h4_h8_design": harmonic_design(rectangle_650),
            "interpretation": (
                "Each family contains two cyclic rows and one Smith=(5,130) "
                "noncyclic row. A saturated three-row H0/H4/H8 fit therefore "
                "cannot by itself distinguish angular H8 from Smith/quotient loading."
            ),
        },
        "first_common_three_direction_design": first_any,
        "first_common_three_primitive_cyclic_design": first_cyclic,
        "n2210_clean_cyclic_design": {
            "square_best_three": best_triplet(square_2210),
            "rectangle_best_three": best_triplet(rectangle_2210),
            "square_all_four": [list(rep) for rep in square_2210],
            "rectangle_all_four": [list(rep) for rep in rectangle_2210],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/n650-angular-smith-preflight/latest.json"),
    )
    args = parser.parse_args()
    payload = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
