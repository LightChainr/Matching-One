#!/usr/bin/env python3
"""Exact all-order Q-jet implied by the Arguin homology identity.

If H(Q,tau)=(Q-1) P(Q,tau), ordinary Q derivatives and logarithmic
derivatives D=Q*d/dQ obey different coefficient tables.  This oracle builds
both tables exactly and checks their equivalence through Stirling transforms.
"""

from __future__ import annotations

import argparse
from math import comb, factorial
import json
from pathlib import Path


def stirling_second(n: int, k: int) -> int:
    if n == k == 0:
        return 1
    if n == 0 or k == 0:
        return 0
    table = [[0] * (k + 1) for _ in range(n + 1)]
    table[0][0] = 1
    for row in range(1, n + 1):
        for column in range(1, min(row, k) + 1):
            table[row][column] = (
                table[row - 1][column - 1]
                + column * table[row - 1][column]
            )
    return table[n][k]


def stirling_first_signed(n: int, k: int) -> int:
    if n == k == 0:
        return 1
    if n == 0 or k == 0:
        return 0
    table = [[0] * (k + 1) for _ in range(n + 1)]
    table[0][0] = 1
    for row in range(1, n + 1):
        for column in range(1, min(row, k) + 1):
            table[row][column] = (
                table[row - 1][column - 1]
                - (row - 1) * table[row - 1][column]
            )
    return table[n][k]


def ordinary_to_log_coefficients(order: int) -> list[int]:
    """Convert d_Q^n H=n d_Q^(n-1)P to D derivatives of P at Q=1."""

    coefficients = [0] * order
    for derivative_h in range(1, order + 1):
        weight = stirling_second(order, derivative_h) * derivative_h
        ordinary_p = derivative_h - 1
        if ordinary_p == 0:
            coefficients[0] += weight
            continue
        for log_p in range(1, ordinary_p + 1):
            coefficients[log_p] += weight * stirling_first_signed(
                ordinary_p, log_p
            )
    return coefficients


def pascal_coefficients(order: int) -> list[int]:
    """D^n[(Q-1)P]|_1 in the basis P,DP,...,D^(n-1)P."""

    return [comb(order, derivative) for derivative in range(order)]


def polynomial_basis_check(order: int, degree: int) -> dict:
    """Check P(Q)=(Q-1)^degree with exact exponential generating functions."""

    left = factorial(degree + 1) * stirling_second(order, degree + 1)
    right = sum(
        comb(order, derivative)
        * factorial(degree)
        * stirling_second(derivative, degree)
        for derivative in range(order)
    )
    return {
        "P": f"(Q-1)^{degree}",
        "D_order_H_at_Q1": left,
        "pascal_reconstruction": right,
        "equal": left == right,
    }


def build_oracle(max_order: int = 10) -> dict:
    if max_order < 1:
        raise ValueError("max_order must be positive")

    rows = []
    for order in range(1, max_order + 1):
        pascal = pascal_coefficients(order)
        stirling = ordinary_to_log_coefficients(order)
        rows.append(
            {
                "order": order,
                "ordinary_Q_identity": {
                    "left": f"partial_Q^{order} H|_1",
                    "right": f"{order} partial_Q^{order - 1} P|_1",
                },
                "log_Q_identity": {
                    "left": f"D^{order} H|_1",
                    "basis": [
                        "P" if derivative == 0 else f"D^{derivative}P"
                        for derivative in range(order)
                    ],
                    "coefficients": pascal,
                },
                "coefficients_from_ordinary_identity_via_stirling": stirling,
                "stirling_conversion_equals_pascal": stirling == pascal,
                "score_mode_combination": "+".join(
                    f"{coefficient}*E[I_trivial*H{derivative}]"
                    for derivative, coefficient in enumerate(pascal)
                ),
                "polynomial_basis_checks": [
                    polynomial_basis_check(order, degree)
                    for degree in range(order)
                ],
            }
        )

    return {
        "schema": "matching-one.p275-arguin-qjet.v1",
        "issues": [252, 258, 263, 275],
        "identity": {
            "definition": "H(Q,tau)=pi_Q(cross;tau)-pi_Q(trivial;tau)",
            "arguin": "H(Q,tau)=(Q-1)P(Q,tau), P=pi_Q(trivial;tau)",
            "ordinary_all_order": "partial_Q^n H|_1=n partial_Q^(n-1)P|_1",
            "logarithmic_all_order": "D^n H|_1=sum_{k=0}^{n-1} binom(n,k) D^k P|_1, D=Q partial_Q",
            "modulus_scope": "holds pointwise for every tau where the Arguin identity holds",
        },
        "basis_warning": {
            "ordinary_third_order": "D^3 H|_1=P+6 partial_Q P+3 partial_Q^2 P",
            "logarithmic_third_order": "D^3 H|_1=P+3 D P+3 D^2 P",
            "do_not_mix": "the coefficient 6 belongs to the ordinary-Q derivative basis, not the D-derivative basis",
        },
        "orders": rows,
        "all_stirling_checks_pass": all(
            row["stirling_conversion_equals_pascal"] for row in rows
        ),
        "all_polynomial_basis_checks_pass": all(
            check["equal"]
            for row in rows
            for check in row["polynomial_basis_checks"]
        ),
        "research_consequence": {
            "forced_jet": "Higher Q-score modes of the homology contrast contain a topology-forced Pascal jet even without a new logarithmic field.",
            "residual_definition": "R_n(tau)=D^n H|_1-sum_{k<n}binom(n,k)D^k P|_1",
            "interpretation": "Only a typed nonzero residual after the forced jet is removed can diagnose an explicit field/projector derivative, a singular collision residue, or finite-lattice violation of the continuum identity.",
            "rank_warning": "Do not count the forced higher jet as evidence for a higher-rank LCFT module.",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-order", type=int, default=10)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(build_oracle(args.max_order), indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
