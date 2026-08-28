#!/usr/bin/env python3
"""Exact combinatorial baselines for Matching One threshold signatures.

Stdlib only.  This module deliberately avoids reading target result files: it
checks structural identities that can be compared with exact/Monte-Carlo data.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from math import comb, factorial


def beta_majority_layer_count(N: int, s: int, k: int) -> int:
    """Successful k-subsets for Beta(s,s) after degree elevation to N.

    Equivalent to majority on m=2s-1 relevant variables and N-m dummies.
    """
    m = 2 * s - 1
    if N < m:
        raise ValueError("N must be at least 2s-1")
    total = 0
    for j in range(s, m + 1):
        if 0 <= k - j <= N - m:
            total += comb(m, j) * comb(N - m, k - j)
    return total


def beta_majority_domination(N: int, s: int) -> list[Fraction]:
    return [
        Fraction(beta_majority_layer_count(N, s, k), comb(N, k))
        for k in range(N + 1)
    ]


def activation_signature_from_domination(d: list[Fraction]) -> list[Fraction]:
    out = []
    previous = Fraction(0)
    for value in d:
        out.append(value - previous)
        previous = value
    # d_0 normally contributes q_0=0 for nontrivial wrapping events.  Return
    # the full 0..N indexing because it makes exact layer conversions explicit.
    return out


def beta_center_density(s: int) -> Fraction:
    """Beta(s,s) density at 1/2 as an exact rational number."""
    beta = Fraction(factorial(s - 1) ** 2, factorial(2 * s - 1))
    return Fraction(1, 2 ** (2 * s - 2)) / beta


def beta_kappas(s: int) -> tuple[Fraction, Fraction]:
    """Return exact kappa3,kappa5 for M=2 I_p(s,s)-1 at p=1/2."""
    if s < 2:
        raise ValueError("s must be >=2")
    d = beta_center_density(s)
    k3 = Fraction(-2 * (s - 1), 1) / (d * d)
    k5 = Fraction(12 * (s - 1) * (s - 2), 1) / (d ** 4)
    return k3, k5


def signed_fourier_level_sums_from_lambda_coefficients(
    coefficients: dict[int, Fraction],
) -> dict[int, Fraction]:
    """Convert product-measure lambda coefficients to signed Walsh level sums.

    If E_lambda[f] = sum_r c_r lambda^r for p_i=1/2+sigma_i lambda,
    then c_r = 2^r sum_{|S|=r} sigma_S fhat(S).
    """
    return {r: c / (2 ** r) for r, c in sorted(coefficients.items())}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    n10_s3 = beta_majority_layer_count(10, 3, 3)
    n26_s5 = beta_majority_layer_count(26, 5, 5)
    n26_s7_k5 = beta_majority_layer_count(26, 7, 5)

    kappa_rows = {}
    for s in (3, 5, 7):
        k3, k5 = beta_kappas(s)
        kappa_rows[str(s)] = {
            "kappa3": str(k3),
            "kappa3_decimal": float(k3),
            "kappa5": str(k5),
            "kappa5_decimal": float(k5),
        }

    n10_interactions = signed_fourier_level_sums_from_lambda_coefficients(
        {1: Fraction(5, 4), 3: Fraction(0), 5: Fraction(-4)}
    )

    out = {
        "beta_majority_minimum_layer": {
            "N10_s3_a3": n10_s3,
            "N26_s5_a5": n26_s5,
            "N26_s7_a5": n26_s7_k5,
        },
        "beta_kappas": kappa_rows,
        "gaussian_majority_limit": {
            "kappa3": -math.pi / 2,
            "kappa5": 3 * math.pi * math.pi / 4,
        },
        "n10_tangent_signed_fourier_level_sums": {
            str(k): str(v) for k, v in n10_interactions.items()
        },
    }

    # Internal exact guards.
    assert n10_s3 == 10
    assert n26_s5 == 126
    assert n26_s7_k5 == 0
    assert beta_kappas(3)[0] == Fraction(-256, 225)
    assert beta_kappas(3)[1] == Fraction(32768, 16875)
    assert n10_interactions == {
        1: Fraction(5, 8),
        3: Fraction(0),
        5: Fraction(-1, 8),
    }

    if args.json:
        print(json.dumps(out, indent=2, sort_keys=True))
    else:
        print("PASS reliability-signature baselines")
        print("N10 Beta(3,3) minimum-layer count:", n10_s3)
        print("N26 Beta(5,5) minimum-layer count:", n26_s5)
        print("Gaussian kappa3 baseline:", -math.pi / 2)
        print("Gaussian kappa5 baseline:", 3 * math.pi * math.pi / 4)


if __name__ == "__main__":
    main()
