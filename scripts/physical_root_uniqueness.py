#!/usr/bin/env python3
"""Exact positive-density and activation-median bracket oracle for Issue #113."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Mapping

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from integer_period_torus import axis_integer_torus  # noqa: E402
from threshold_rank_nz import enumerate_exact  # noqa: E402


def binomial_tail(N: int, k: int, p: float) -> float:
    if k == 0:
        return 1.0
    if k == N + 1:
        return 0.0
    if not 1 <= k <= N:
        raise ValueError("activation rank must lie in 0..N+1")
    return sum(
        math.comb(N, n) * p**n * (1.0 - p) ** (N - n)
        for n in range(k, N + 1)
    )


def beta_density(N: int, k: int, p: float) -> float:
    if k in (0, N + 1):
        return 0.0
    if not 1 <= k <= N:
        raise ValueError("activation rank must lie in 0..N+1")
    return (
        N
        * math.comb(N - 1, k - 1)
        * p ** (k - 1)
        * (1.0 - p) ** (N - k)
    )


def normalized_histogram(histogram: Mapping[int, int]) -> dict[int, float]:
    total = sum(histogram.values())
    if total <= 0 or any(count < 0 for count in histogram.values()):
        raise ValueError("histogram counts must be nonnegative with positive total")
    return {int(rank): count / total for rank, count in histogram.items() if count}


def activation_cdf(N: int, histogram: Mapping[int, int], p: float) -> float:
    return sum(
        weight * binomial_tail(N, rank, p)
        for rank, weight in normalized_histogram(histogram).items()
    )


def activation_density(N: int, histogram: Mapping[int, int], p: float) -> float:
    return sum(
        weight * beta_density(N, rank, p)
        for rank, weight in normalized_histogram(histogram).items()
    )


def monotone_root(function, target: float = 0.0) -> float:
    left, right = 0.0, 1.0
    if function(left) > target or function(right) < target:
        raise ValueError("target is not bracketed on [0,1]")
    for _ in range(120):
        middle = (left + right) / 2.0
        if function(middle) < target:
            left = middle
        else:
            right = middle
    return (left + right) / 2.0


def activation_bracket(
    N: int, onset: Mapping[int, int], completion: Mapping[int, int]
) -> dict[str, float]:
    onset_cdf = lambda p: activation_cdf(N, onset, p)
    completion_cdf = lambda p: activation_cdf(N, completion, p)
    matching = lambda p: onset_cdf(p) + completion_cdf(p) - 1.0
    m1 = monotone_root(onset_cdf, 0.5)
    root = monotone_root(matching, 0.0)
    m2 = monotone_root(completion_cdf, 0.5)
    if not m1 <= root <= m2:
        raise AssertionError("stochastic-order root bracket failed")
    derivative = activation_density(N, onset, root) + activation_density(
        N, completion, root
    )
    if derivative <= 0.0:
        raise AssertionError("physical root is not simple")
    return {
        "onset_mixture_median": m1,
        "physical_matching_root": root,
        "completion_mixture_median": m2,
        "M_derivative_at_root": derivative,
        "root_balance_residual": matching(root),
    }


def exact_axis_L2_oracle() -> dict[str, object]:
    geometry = axis_integer_torus(2)
    counts = enumerate_exact(geometry)
    onset = {rank: count for rank, count in enumerate(counts.kminus) if count}
    completion = {rank: count for rank, count in enumerate(counts.kplus) if count}
    bracket = activation_bracket(geometry.n, onset, completion)
    strict_grid = all(
        activation_density(geometry.n, onset, index / 20.0)
        + activation_density(geometry.n, completion, index / 20.0)
        > 0.0
        for index in range(1, 20)
    )
    return {
        "schema": "matching-one.physical-root-uniqueness.v1",
        "issue": 113,
        "dependency": "Issue 28 two-activation Beta mixture and Issue 269 rank theorem",
        "theorem": {
            "honest_torus_endpoints": {"M(0)": -1, "M(1)": 1},
            "open_interval_derivative": "E[Beta(K1,N-K1+1)+Beta(K2,N-K2+1)]>0",
            "conclusion": "exactly one simple real root in (0,1)",
        },
        "degenerate_extension": {
            "K=0": "preactivated endpoint atom; no derivative contribution",
            "K=N+1": "never-activated endpoint atom; no derivative contribution",
            "M(0)": "-1+P(K1=0)+P(K2=0)",
            "M(1)": "1-P(K1=N+1)-P(K2=N+1)",
        },
        "axis_L2_exact_histogram": {
            "N": geometry.n,
            "samples": counts.sample_count,
            "onset_histogram": onset,
            "completion_histogram": completion,
            "bracket": bracket,
            "strict_positive_derivative_grid_check": strict_grid,
        },
        "complex_root_boundary": "no constraint on nonreal or out-of-interval roots; matching zeros are not automatically Fisher or Lee-Yang zeros",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/physical_root_uniqueness_exact.json"),
    )
    args = parser.parse_args()
    payload = exact_axis_L2_oracle()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
