#!/usr/bin/env python3
"""Exact neutral-area Krawtchouk covector and bounded mode-front diagnostic.

For the neutral indicator U_k=1-m_k^2, marginal threshold histograms give

    q_k = P(K_minus <= k) - P(K_plus <= k).

The canonical area is both E[K_plus-K_minus]/(N+1) and a known covector of
the complete Krawtchouk coefficient vector.  This script verifies the first
identity exactly from histogram totals, evaluates the covector through a
frozen maximum order, and reports whether the remaining exact area tail has
fallen below predeclared relative tolerances.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import mpmath as mp

from analyze_matching_parity_derivatives_fast import H, combine, read
from c4_self_matching_exact import enumerate_exact
from threshold_score_modes import binomial_weights, center, krawtchouk_mode


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def neutral_curve(row: H) -> list[mp.mpf]:
    """Return q_k=P(Kminus<=k<Kplus) from threshold marginals."""

    minus_cdf = 0
    plus_cdf = 0
    curve: list[mp.mpf] = []
    for occupation in range(row.n + 1):
        minus_cdf += row.minus[occupation]
        plus_cdf += row.plus[occupation]
        value = mp.mpf(minus_cdf - plus_cdf) / row.samples
        if value < 0 or value > 1:
            raise ValueError("threshold ordering produced an invalid neutral curve")
        curve.append(value)
    return curve


def mean_rank_gap(row: H) -> mp.mpf:
    return mp.fsum(
        occupation * (row.plus[occupation] - row.minus[occupation])
        for occupation in range(row.n + 1)
    ) / row.samples


def area_covector(n: int, p0: mp.mpf, order: int) -> mp.mpf:
    """Return the exact integral covector multiplying Krawtchouk c_order."""

    if not 0 < p0 < 1 or not 0 <= order <= n:
        raise ValueError("invalid center or mode order")
    moment = (
        mp.power(1 - p0, order + 1) - mp.power(-p0, order + 1)
    ) / (order + 1)
    return (
        mp.sqrt(math.comb(n, order))
        * moment
        / mp.power(p0 * (1 - p0), mp.mpf(order) / 2)
    )


def score_coefficients(
    curve: Sequence[mp.mpf], n: int, p0: mp.mpf, max_order: int
) -> list[mp.mpf]:
    if len(curve) != n + 1 or not 0 <= max_order <= n:
        raise ValueError("curve length or maximum order is invalid")
    weights = binomial_weights(n, p0)
    return [
        mp.fsum(
            weights[k]
            * curve[k]
            * krawtchouk_mode(n, k, order, p0)
            for k in range(n + 1)
        )
        for order in range(max_order + 1)
    ]


def first_mode_within_tolerance(
    relative_tails: Sequence[mp.mpf], epsilon: mp.mpf
) -> int | None:
    """Return the first truncation order whose exact aggregate tail is small."""

    if not 0 < epsilon < 1:
        raise ValueError("epsilon must lie strictly between zero and one")
    return next(
        (order for order, value in enumerate(relative_tails) if value <= epsilon),
        None,
    )


def mode_front(
    curve: Sequence[mp.mpf],
    n: int,
    p0: mp.mpf,
    exact_area: mp.mpf,
    max_order: int,
    epsilons: Sequence[mp.mpf],
) -> dict[str, object]:
    coefficients = score_coefficients(curve, n, p0, max_order)
    contributions = [
        area_covector(n, p0, order) * coefficient
        for order, coefficient in enumerate(coefficients)
    ]
    partial = mp.mpf(0)
    tails: list[mp.mpf] = []
    relative_tails: list[mp.mpf] = []
    for contribution in contributions:
        partial += contribution
        tail = exact_area - partial
        tails.append(tail)
        relative_tails.append(abs(tail) / abs(exact_area))
    thresholds: dict[str, object] = {}
    for epsilon in epsilons:
        order = first_mode_within_tolerance(relative_tails, epsilon)
        thresholds[mp.nstr(epsilon, 8)] = {
            "R_epsilon": order,
            "R_epsilon_over_N_one_quarter": (
                None if order is None else mp.nstr(order / mp.power(n, mp.mpf(1) / 4), 15)
            ),
            "resolved_within_max_order": order is not None,
        }
    return {
        "coefficients": [mp.nstr(value, 30) for value in coefficients],
        "area_contributions": [mp.nstr(value, 30) for value in contributions],
        "tail_after_order": [mp.nstr(value, 30) for value in tails],
        "relative_abs_tail_after_order": [
            mp.nstr(value, 15) for value in relative_tails
        ],
        "epsilon_diagnostics": thresholds,
    }


def exact_n10_oracle() -> dict[str, object]:
    """Return an exact p0=1/2 oracle from the self-matching N=10 quotient."""

    result = enumerate_exact(3, 1)
    n = int(result["geometry"]["N"])
    wrapping = result["channels"]["cross"]["R_bernstein_integer_coefficients"]
    curve = [
        Fraction(
            math.comb(n, k) - wrapping[k] - wrapping[n - k],
            math.comb(n, k),
        )
        for k in range(n + 1)
    ]
    scaled_coefficients: list[Fraction] = []
    contributions: list[Fraction] = []
    for order in range(n + 1):
        value = Fraction()
        for k in range(n + 1):
            lower = max(0, order - (n - k))
            upper = min(order, k)
            polynomial = sum(
                math.comb(k, j)
                * math.comb(n - k, order - j)
                * (-1) ** j
                for j in range(lower, upper + 1)
            )
            value += (
                Fraction(math.comb(n, k), 2**n)
                * curve[k]
                * (-1) ** order
                * polynomial
            )
        scaled_coefficients.append(value)
        contributions.append(value / (order + 1) if order % 2 == 0 else Fraction())
    area = sum(curve, Fraction()) / (n + 1)
    if any(scaled_coefficients[order] for order in range(1, n + 1, 2)):
        raise AssertionError("N=10 half-filled oracle has a nonzero odd mode")
    if sum(contributions, Fraction()) != area:
        raise AssertionError("N=10 exact area covector failed")
    partial = Fraction()
    relative_tails: list[Fraction] = []
    for contribution in contributions:
        partial += contribution
        relative_tails.append(abs(area - partial) / area)
    exact_epsilon_orders = {
        "0.05": next(
            order for order, value in enumerate(relative_tails) if value <= Fraction(1, 20)
        ),
        "0.10": next(
            order for order, value in enumerate(relative_tails) if value <= Fraction(1, 10)
        ),
    }
    return {
        "N": n,
        "p0": "1/2",
        "channel": "cross",
        "neutral_curve": [str(value) for value in curve],
        "c_r_times_sqrt_binomial": [str(value) for value in scaled_coefficients],
        "area_contributions": [str(value) for value in contributions],
        "canonical_neutral_area": str(area),
        "mean_rank_gap": str((n + 1) * area),
        "odd_modes_zero": True,
        "R_epsilon": exact_epsilon_orders,
    }


def analyze_histogram(path: Path, max_order: int, epsilons: Sequence[mp.mpf]) -> dict:
    data = read(path)
    sizes = {key[0] for key in data}
    if len(sizes) != 1:
        raise ValueError(f"{path}: expected one size per input")
    n = sizes.pop()
    grouped = {
        orientation: [
            data[key]
            for key in sorted(data)
            if key[0] == n and key[1] == orientation
        ]
        for orientation in ("first", "second")
    }
    if len(grouped["first"]) != len(grouped["second"]) or not grouped["first"]:
        raise ValueError(f"N={n}: orientations are not aligned")
    first = combine(grouped["first"])
    second = combine(grouped["second"])
    p0 = center(first, second)
    first_curve = neutral_curve(first)
    second_curve = neutral_curve(second)
    pooled_curve = [
        (left + right) / 2 for left, right in zip(first_curve, second_curve)
    ]
    gap = (mean_rank_gap(first) + mean_rank_gap(second)) / 2
    area_from_gap = gap / (n + 1)
    area_from_layers = mp.fsum(pooled_curve) / (n + 1)
    if abs(area_from_gap - area_from_layers) > mp.mpf("1e-60"):
        raise AssertionError(f"N={n}: exact histogram area bridge failed")
    return {
        "N": n,
        "representations": {
            "first": [first.a, first.b],
            "second": [second.a, second.b],
        },
        "samples_per_orientation": first.samples,
        "batches": len(grouped["first"]),
        "p0": mp.nstr(p0, 30),
        "mean_rank_gap": mp.nstr(gap, 30),
        "canonical_neutral_area": mp.nstr(area_from_gap, 30),
        "histogram_area_bridge_abs_error": mp.nstr(
            abs(area_from_gap - area_from_layers), 10
        ),
        "mode_front": mode_front(
            pooled_curve, n, p0, area_from_gap, max_order, epsilons
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("histograms", nargs="+", type=Path)
    parser.add_argument("--max-order", type=int, default=16)
    parser.add_argument("--epsilon", nargs="+", default=("0.05", "0.10"))
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    if not 1 <= args.max_order <= 32:
        raise SystemExit("--max-order must lie in 1..32")
    epsilons = tuple(mp.mpf(value) for value in args.epsilon)
    mp.mp.dps = args.dps
    rows = [analyze_histogram(path, args.max_order, epsilons) for path in args.histograms]
    if len({row["N"] for row in rows}) != len(rows):
        raise SystemExit("input histograms contain duplicate sizes")
    payload = {
        "schema": "matching-one/neutral-area-mode-front/v1",
        "status": "retrospective_existing_data_diagnostic",
        "exact_covector": (
            "E[G]/(N+1)=sum_r c_r*sqrt(C(N,r))*"
            "((1-p0)^(r+1)-(-p0)^(r+1))/"
            "((r+1)*(p0*(1-p0))^(r/2))"
        ),
        "half_filling_reduction": (
            "E[G]/(N+1)=sum_even_r sqrt(C(N,r))*c_r/(r+1)"
        ),
        "max_order": args.max_order,
        "epsilons": [mp.nstr(value, 8) for value in epsilons],
        "inputs": [
            {"path": str(path), "sha256": sha256(path)} for path in args.histograms
        ],
        "exact_n10_oracle": exact_n10_oracle(),
        "by_N": {str(row["N"]): row for row in sorted(rows, key=lambda row: row["N"])},
        "interpretation_boundary": (
            "The covector identity is exact. R_epsilon is a bounded-order, "
            "retrospective diagnostic; unresolved means only that r<=max_order "
            "does not reconstruct the global area in this local basis."
        ),
    }
    rendered = json.dumps(payload, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(args.output)
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
