#!/usr/bin/env python3
"""Construct finite-size annihilator weights at arbitrary precision.

For values with an expansion

    q(L) = q_inf + sum_r a_r L**(-alpha_r) + ...,

choose weights w_j satisfying

    sum_j w_j = 1,
    sum_j w_j L_j**(-alpha_r) = 0

for every requested correction exponent alpha_r. The weighted combination then
removes those modeled corrections. The same weights can be used in a root
condition sum_j w_j M_{L_j}(p) = 0 for finite-size matching functions.

When more sizes than constraints are supplied, the script returns the minimum-L2
norm solution. Large weight norms or a large constraint condition number warn
that statistical noise and exponent misspecification will be amplified.

Examples:
    # Equivalent, up to normalization, to the known two-size matching equation.
    python scripts/matching_annihilator.py --sizes 15 16 --cancel 13/4

    # Cancel the L^(-13/4) and L^(-25/4) terms with three sizes.
    python scripts/matching_annihilator.py --sizes 14 15 16 --cancel 13/4 25/4
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import mpmath as mp


@dataclass(frozen=True)
class Result:
    sizes: list[int]
    cancel_exponents: list[str]
    weights: list[str]
    residuals: list[str]
    l1_norm: str
    l2_norm: str
    max_abs_weight: str
    constraint_condition: str
    weighted_value: str | None


def parse_exponent(text: str) -> mp.mpf:
    """Parse a decimal/scientific number or an exact fraction such as 13/4."""
    try:
        if "/" in text:
            fraction = Fraction(text)
            return mp.mpf(fraction.numerator) / fraction.denominator
        return mp.mpf(text)
    except (ValueError, ZeroDivisionError) as exc:
        raise argparse.ArgumentTypeError(f"invalid exponent {text!r}") from exc


def format_mpf(value: mp.mpf, digits: int = 40) -> str:
    return mp.nstr(value, n=digits, strip_zeros=False)


def constraint_matrix(sizes: Sequence[int], exponents: Sequence[mp.mpf]) -> mp.matrix:
    rows: list[list[mp.mpf]] = [[mp.mpf(1) for _ in sizes]]
    for exponent in exponents:
        rows.append([mp.power(size, -exponent) for size in sizes])
    return mp.matrix(rows)


def solve_weights(matrix: mp.matrix) -> mp.matrix:
    constraints, sizes_count = matrix.rows, matrix.cols
    if sizes_count < constraints:
        raise ValueError(
            f"need at least {constraints} sizes for {constraints - 1} cancellation exponents"
        )
    target = mp.matrix([mp.mpf(1), *([mp.mpf(0)] * (constraints - 1))])
    if sizes_count == constraints:
        return mp.lu_solve(matrix, target)

    # Minimum Euclidean norm solution to A w = b:
    # w = A^T (A A^T)^(-1) b.
    gram = matrix * matrix.T
    multipliers = mp.lu_solve(gram, target)
    return matrix.T * multipliers


def matrix_condition(matrix: mp.matrix) -> mp.mpf:
    """2-norm condition estimate of the row-constraint map via its Gram matrix."""
    gram = matrix * matrix.T
    eigenvalues = mp.eigsy(gram, eigvals_only=True)
    positive = [value for value in eigenvalues if value > 0]
    if len(positive) != gram.rows:
        return mp.inf
    return mp.sqrt(max(positive) / min(positive))


def compute(
    sizes: Sequence[int],
    exponents: Sequence[mp.mpf],
    values: Sequence[mp.mpf] | None,
) -> Result:
    if not sizes:
        raise ValueError("at least one size is required")
    if any(size <= 0 for size in sizes):
        raise ValueError("sizes must be positive")
    if len(set(sizes)) != len(sizes):
        raise ValueError("sizes must be distinct")
    if any(exponent <= 0 for exponent in exponents):
        raise ValueError("cancellation exponents must be positive")
    if values is not None and len(values) != len(sizes):
        raise ValueError("--values must contain one value per size")

    matrix = constraint_matrix(sizes, exponents)
    weights = solve_weights(matrix)
    target = mp.matrix([mp.mpf(1), *([mp.mpf(0)] * len(exponents))])
    residual_vector = matrix * weights - target

    l1_norm = mp.fsum(abs(weight) for weight in weights)
    l2_norm = mp.sqrt(mp.fsum(weight * weight for weight in weights))
    max_abs = max(abs(weight) for weight in weights)
    weighted = None
    if values is not None:
        if len(weights) != len(values):
            raise RuntimeError("weight and value counts differ")
        weighted = mp.fsum(weight * value for weight, value in zip(weights, values))

    return Result(
        sizes=list(sizes),
        cancel_exponents=[format_mpf(exponent, 20) for exponent in exponents],
        weights=[format_mpf(weight) for weight in weights],
        residuals=[format_mpf(value, 15) for value in residual_vector],
        l1_norm=format_mpf(l1_norm, 20),
        l2_norm=format_mpf(l2_norm, 20),
        max_abs_weight=format_mpf(max_abs, 20),
        constraint_condition=format_mpf(matrix_condition(matrix), 20),
        weighted_value=None if weighted is None else format_mpf(weighted),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", nargs="+", type=int, required=True)
    parser.add_argument(
        "--cancel",
        nargs="*",
        type=parse_exponent,
        default=[],
        metavar="ALPHA",
        help="correction exponents to annihilate; accepts decimals or fractions",
    )
    parser.add_argument(
        "--values",
        nargs="+",
        type=mp.mpf,
        default=None,
        help="optional q(L) values to combine with the computed weights",
    )
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument("--json", type=Path, default=None)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.dps < 50:
        raise SystemExit("use at least 50 decimal digits")
    mp.mp.dps = args.dps
    try:
        result = compute(args.sizes, args.cancel, args.values)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    print("sizes:", " ".join(map(str, result.sizes)))
    print("cancel exponents:", ", ".join(result.cancel_exponents) or "none")
    print("weights:")
    if len(result.sizes) != len(result.weights):
        raise RuntimeError("size and weight counts differ")
    for size, weight in zip(result.sizes, result.weights):
        print(f"  L={size}: {weight}")
    print("constraint residuals:", ", ".join(result.residuals))
    print("L1 noise factor:", result.l1_norm)
    print("L2 independent-noise factor:", result.l2_norm)
    print("max |weight|:", result.max_abs_weight)
    print("constraint condition:", result.constraint_condition)
    if result.weighted_value is not None:
        print("weighted value:", result.weighted_value)

    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(asdict(result), indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
