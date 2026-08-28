#!/usr/bin/env python3
"""Exact tangent algebra for a staggered self-matching checkerboard family.

On the degree-8/even and degree-4/odd sublattices set

    p_even = 1/2 + t + lambda,
    p_odd  = 1/2 + t - lambda.

Occupation complement maps ``(t, lambda)`` exactly to ``(-t, -lambda)``.
This script exhausts the N=10 Gaussian quotient, constructs the exact
two-variable wrapping polynomials, and reports the exchange-even/odd response.
It is a local/RG-tangent gate, not a large-N scaling-field identification.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Dict, Sequence, Tuple

from c4_self_matching_exact import (
    CHANNELS,
    _active,
    c4_self_matching_torus,
)
from integer_period_torus import classify_configuration


Exponent = Tuple[int, int]  # powers of (t, lambda)
Polynomial = Dict[Exponent, Fraction]


def _add(left: Polynomial, right: Polynomial) -> Polynomial:
    output = dict(left)
    for exponent, coefficient in right.items():
        output[exponent] = output.get(exponent, Fraction()) + coefficient
        if output[exponent] == 0:
            del output[exponent]
    return output


def _multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    output: Polynomial = {}
    for (i, j), x in left.items():
        for (k, ell), y in right.items():
            exponent = (i + k, j + ell)
            output[exponent] = output.get(exponent, Fraction()) + x * y
    return {exponent: value for exponent, value in output.items() if value}


def _power(base: Polynomial, exponent: int) -> Polynomial:
    output: Polynomial = {(0, 0): Fraction(1)}
    for _ in range(exponent):
        output = _multiply(output, base)
    return output


def _scale(polynomial: Polynomial, value: int) -> Polynomial:
    return {exponent: value * coefficient for exponent, coefficient in polynomial.items()}


def _fraction(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def sublattice_bernstein_counts(a: int = 3, b: int = 1) -> dict[str, list[list[int]]]:
    """Count wrapping masks by occupied even/odd sublattice populations."""

    geometry = c4_self_matching_torus(a, b)
    if geometry.n > 22:
        raise ValueError("exact staggered enumeration is restricted to N<=22")
    parity = tuple((x + y) % 2 for x, y in geometry.coordinates)
    n_even = parity.count(0)
    n_odd = parity.count(1)
    counts = {
        channel: [[0] * (n_odd + 1) for _ in range(n_even + 1)]
        for channel in CHANNELS
    }
    for mask in range(1 << geometry.n):
        wrapping, _ = classify_configuration(geometry, _active(mask, geometry.n))
        k_even = sum(bool(mask & (1 << i)) for i, value in enumerate(parity) if value == 0)
        k_odd = sum(bool(mask & (1 << i)) for i, value in enumerate(parity) if value == 1)
        for channel in CHANNELS:
            counts[channel][k_even][k_odd] += int(getattr(wrapping, channel))
    return counts


def wrapping_polynomial(counts: list[list[int]]) -> Polynomial:
    """Return exact R(t,lambda) from sublattice Bernstein counts."""

    n_even = len(counts) - 1
    n_odd = len(counts[0]) - 1
    half = Fraction(1, 2)
    p_even = {(0, 0): half, (1, 0): Fraction(1), (0, 1): Fraction(1)}
    q_even = {(0, 0): half, (1, 0): Fraction(-1), (0, 1): Fraction(-1)}
    p_odd = {(0, 0): half, (1, 0): Fraction(1), (0, 1): Fraction(-1)}
    q_odd = {(0, 0): half, (1, 0): Fraction(-1), (0, 1): Fraction(1)}
    output: Polynomial = {}
    for k_even, row in enumerate(counts):
        for k_odd, count in enumerate(row):
            term = _multiply(
                _multiply(_power(p_even, k_even), _power(q_even, n_even - k_even)),
                _multiply(_power(p_odd, k_odd), _power(q_odd, n_odd - k_odd)),
            )
            output = _add(output, _scale(term, count))
    return output


def _terms(polynomial: Polynomial) -> list[dict[str, object]]:
    return [
        {"t_power": i, "lambda_power": j, "coefficient": _fraction(value)}
        for (i, j), value in sorted(polynomial.items(), key=lambda item: (sum(item[0]), item[0]))
    ]


def _lambda_coefficients(polynomial: Polynomial, exchange_parity: int) -> list[str]:
    """Restrict t=0 and select total exchange parity (+1 or -1)."""

    degree = max((j for (i, j) in polynomial if i == 0), default=0)
    coefficients = []
    for j in range(degree + 1):
        value = polynomial.get((0, j), Fraction())
        if (-1) ** j != exchange_parity:
            value = Fraction()
        coefficients.append(_fraction(value))
    while len(coefficients) > 1 and coefficients[-1] == "0":
        coefficients.pop()
    return coefficients


def exact_tangent_report(a: int = 3, b: int = 1) -> dict[str, object]:
    counts = sublattice_bernstein_counts(a, b)
    geometry = c4_self_matching_torus(a, b)
    channels: dict[str, object] = {}
    for channel in CHANNELS:
        polynomial = wrapping_polynomial(counts[channel])
        r_t = polynomial.get((1, 0), Fraction())
        r_lambda = polynomial.get((0, 1), Fraction())
        channels[channel] = {
            "R_polynomial_terms": _terms(polynomial),
            "R_plus_at_t0_lambda_coefficients_ascending": _lambda_coefficients(polynomial, +1),
            "R_minus_at_t0_lambda_coefficients_ascending": _lambda_coefficients(polynomial, -1),
            "response_matrix_rows_Rplus_Rminus_columns_t_lambda": [
                ["0", "0"],
                [_fraction(r_t), _fraction(r_lambda)],
            ],
        }

    expected_odd = ["0", "5/4", "0", "0", "0", "-4"]
    if any(
        row["R_minus_at_t0_lambda_coefficients_ascending"] != expected_odd
        for row in channels.values()
    ):
        raise AssertionError("N=10 channel-independent odd tangent polynomial changed")

    nonphysical_root = (5.0 / 16.0) ** 0.25
    return {
        "schema": "matching-one/c4-self-matching-tangent/v1",
        "family": {
            "p_even": "1/2 + t + lambda",
            "p_odd": "1/2 + t - lambda",
            "legal_domain": "abs(t+lambda)<=1/2 and abs(t-lambda)<=1/2",
            "matching_exchange": "(t,lambda)->(-t,-lambda)",
            "pair_convention": "Rplus=(R(t,lambda)+R(-t,-lambda))/2; Rminus=(R(t,lambda)-R(-t,-lambda))/2",
        },
        "geometry": {"a": a, "b": b, "N": geometry.n, "even_sites": geometry.n // 2, "odd_sites": geometry.n // 2},
        "selection_rule": "A Taylor monomial t^m lambda^n is exchange-even iff m+n is even, exchange-odd iff m+n is odd.",
        "channels": channels,
        "exact_odd_root_gate": {
            "Rminus_all_channels": "(5/4)*lambda - 4*lambda^5",
            "factorization": "lambda*(5/4 - 4*lambda^4)",
            "legal_interval_at_t0": "-1/2 <= lambda <= 1/2",
            "only_legal_root": "0",
            "other_real_root_magnitudes": nonphysical_root,
            "other_roots_are_legal": False,
        },
        "nontrivial_improved_action_protocol": {
            "target": "matching-even H4 amplitude",
            "definition": "A_T4_plus(N,lambda)=N*(Rplus_(a1,b1)-Rplus_(a2,b2))/(cos(4theta1)-cos(4theta2))",
            "minimum_orientation_design": {"N": 130, "first": [11, 3], "second": [9, 7]},
            "frozen_nonnegative_scan": ["0", "1/8", "1/4", "3/8"],
            "fit": "Fit the first three points to a0+a2*z+a4*z^2 with z=lambda^2; keep lambda=3/8 as a no-refit lack-of-fit test.",
            "root_acceptance": "Accept lambda_star=sqrt(z_star) only if 0<z_star<1/4, the held-out 3/8 residual passes, and an antithetic run at +/-lambda_star brackets zero.",
            "replication": "Repeat the frozen root without refitting at N=170, (13,1) versus (11,7).",
        },
        "passed": True,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a", type=int, default=3)
    parser.add_argument("--b", type=int, default=1)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args(argv)
    report = exact_tangent_report(args.a, args.b)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
