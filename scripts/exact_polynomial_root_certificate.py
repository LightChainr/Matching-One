#!/usr/bin/env python3
"""Exact rational Sturm isolation and stationary-point classification."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Sequence

from threshold_histogram_profile import density_coefficients, mixture_weights, parse_histogram


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "threshold_histogram_profile_contract.json"


def trim(polynomial: Sequence[Fraction]) -> list[Fraction]:
    result = [Fraction(value) for value in polynomial]
    while result and result[-1] == 0:
        result.pop()
    return result


def derivative(polynomial: Sequence[Fraction]) -> list[Fraction]:
    return trim([Fraction(index) * Fraction(value) for index, value in enumerate(polynomial) if index])


def evaluate(polynomial: Sequence[Fraction], value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(polynomial):
        result = result * value + Fraction(coefficient)
    return result


def polynomial_divmod(
    numerator: Sequence[Fraction], denominator: Sequence[Fraction]
) -> tuple[list[Fraction], list[Fraction]]:
    top, bottom = trim(numerator), trim(denominator)
    if not bottom:
        raise ZeroDivisionError("polynomial division by zero")
    if len(top) < len(bottom):
        return [], top
    quotient = [Fraction(0) for _ in range(len(top) - len(bottom) + 1)]
    while top and len(top) >= len(bottom):
        degree = len(top) - len(bottom)
        scale = top[-1] / bottom[-1]
        quotient[degree] = scale
        for index, coefficient in enumerate(bottom):
            top[degree + index] -= scale * coefficient
        top = trim(top)
    return trim(quotient), top


def monic(polynomial: Sequence[Fraction]) -> list[Fraction]:
    values = trim(polynomial)
    if not values:
        return []
    return [value / values[-1] for value in values]


def polynomial_gcd(first: Sequence[Fraction], second: Sequence[Fraction]) -> list[Fraction]:
    left, right = trim(first), trim(second)
    while right:
        _, remainder = polynomial_divmod(left, right)
        left, right = right, remainder
    return monic(left)


def square_free_part(polynomial: Sequence[Fraction]) -> list[Fraction]:
    values = trim(polynomial)
    if len(values) <= 1:
        return values
    common = polynomial_gcd(values, derivative(values))
    quotient, remainder = polynomial_divmod(values, common)
    if remainder:
        raise ArithmeticError("square-free division left a remainder")
    return monic(quotient)


def sturm_sequence(polynomial: Sequence[Fraction]) -> list[list[Fraction]]:
    first = monic(polynomial)
    if len(first) <= 1:
        return [first]
    sequence = [first, derivative(first)]
    while sequence[-1]:
        _, remainder = polynomial_divmod(sequence[-2], sequence[-1])
        if not remainder:
            break
        sequence.append([-value for value in remainder])
    return sequence


def sign_variations(sequence: Sequence[Sequence[Fraction]], point: Fraction) -> int:
    signs = []
    for polynomial in sequence:
        value = evaluate(polynomial, point)
        if value:
            signs.append(1 if value > 0 else -1)
    return sum(left != right for left, right in zip(signs, signs[1:]))


def open_root_count(
    sequence: Sequence[Sequence[Fraction]], left: Fraction, right: Fraction
) -> int:
    if not left < right:
        raise ValueError("root-count interval must be ordered")
    if evaluate(sequence[0], left) == 0 or evaluate(sequence[0], right) == 0:
        raise ValueError("root-count endpoints must not be roots")
    return sign_variations(sequence, left) - sign_variations(sequence, right)


def _divide_linear(polynomial: Sequence[Fraction], root: Fraction) -> list[Fraction]:
    quotient, remainder = polynomial_divmod(polynomial, [-root, Fraction(1)])
    if remainder:
        raise ArithmeticError("declared exact root did not divide polynomial")
    return quotient


def isolate_roots(
    polynomial: Sequence[Fraction],
    left: Fraction = Fraction(0),
    right: Fraction = Fraction(1),
    bits: int = 20,
) -> list[tuple[Fraction, Fraction]]:
    """Isolate every distinct root in [left,right] by exact dyadic bisection."""
    if bits < 1:
        raise ValueError("bits must be positive")
    work = square_free_part(polynomial)
    if len(work) <= 1:
        return []
    exact: list[tuple[Fraction, Fraction]] = []
    for endpoint in (left, right):
        if evaluate(work, endpoint) == 0:
            exact.append((endpoint, endpoint))
            work = _divide_linear(work, endpoint)
    if len(work) <= 1:
        return sorted(exact)

    target_width = Fraction(1, 1 << bits)
    while True:
        sequence = sturm_sequence(work)
        found_exact: Fraction | None = None
        isolated: list[tuple[Fraction, Fraction]] = []

        def visit(lo: Fraction, hi: Fraction, count: int) -> None:
            nonlocal found_exact
            if count == 0 or found_exact is not None:
                return
            if count == 1 and hi - lo <= target_width:
                isolated.append((lo, hi))
                return
            midpoint = (lo + hi) / 2
            if evaluate(work, midpoint) == 0:
                found_exact = midpoint
                return
            left_count = open_root_count(sequence, lo, midpoint)
            right_count = open_root_count(sequence, midpoint, hi)
            if left_count + right_count != count:
                raise ArithmeticError("Sturm root count did not split additively")
            visit(lo, midpoint, left_count)
            visit(midpoint, hi, right_count)

        total = open_root_count(sequence, left, right)
        visit(left, right, total)
        if found_exact is None:
            return sorted(exact + isolated)
        exact.append((found_exact, found_exact))
        work = _divide_linear(work, found_exact)
        if len(work) <= 1:
            return sorted(exact)


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def classify_stationary_points(
    density: Sequence[Fraction], bits: int = 20
) -> list[dict[str, str]]:
    slope = derivative(density)
    roots = isolate_roots(slope, bits=bits)
    output = []
    for index, (lo, hi) in enumerate(roots):
        previous = roots[index - 1][1] if index else Fraction(0)
        following = roots[index + 1][0] if index + 1 < len(roots) else Fraction(1)
        left_probe = (previous + lo) / 2
        right_probe = (hi + following) / 2
        left_sign = evaluate(slope, left_probe)
        right_sign = evaluate(slope, right_probe)
        kind = (
            "strict_maximum"
            if left_sign > 0 > right_sign
            else "strict_minimum"
            if left_sign < 0 < right_sign
            else "stationary_no_sign_change"
        )
        output.append(
            {
                "left": fraction_text(lo),
                "right": fraction_text(hi),
                "classification": kind,
            }
        )
    return output


def build_artifact(contract_path: Path = DEFAULT_CONTRACT) -> dict:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    n = contract["N"]
    minus = parse_histogram(contract["K_minus_counts"], n, "K_minus")
    plus = parse_histogram(contract["K_plus_counts"], n, "K_plus")
    density = density_coefficients(mixture_weights(minus, plus, n))
    stationary = classify_stationary_points(density)
    if stationary != [{"left": "1/2", "right": "1/2", "classification": "strict_maximum"}]:
        raise ValueError("frozen profile stationary-point certificate drifted")
    return {
        "schema": "matching-one/exact-polynomial-root-certificate/v1",
        "issue": 28,
        "data_class": "exact synthetic rational polynomial",
        "density_power_coefficients": [fraction_text(value) for value in density],
        "derivative_power_coefficients": [fraction_text(value) for value in derivative(density)],
        "isolation_bits": 20,
        "distinct_stationary_points": stationary,
        "unique_internal_mode": "1/2",
        "boundary": (
            "Exact rational-polynomial certificate only: no production histogram, empirical mode, "
            "tail fit, or universality claim."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(build_artifact(args.contract), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

