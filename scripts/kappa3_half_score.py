#!/usr/bin/env python3
"""Utilities for kappa_3 at an exact Bernoulli threshold p=1/2.

For an observable D(C) with E_p[D] = M(p), N independent Bernoulli
variables, K occupied variables, and x=2K-N, the p=1/2 score identities are

    M'(1/2)   = E[ 2 x D ],
    M'''(1/2) = E[ 8 (x^3 - (3N-2)x) D ].

These formulas are useful for exact-threshold controls such as square bond
percolation.  They avoid numerical differentiation.  This file intentionally
does not define a topology observable; it provides the score polynomials and
an aggregation helper to be reused by validated geometry-specific code.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import fsum
from typing import Iterable


@dataclass(frozen=True)
class ScoreMoments:
    samples: int
    mean_d: float
    first_derivative: float
    third_derivative: float
    kappa3: float


def first_score(n: int, k: int) -> int:
    """Bernoulli likelihood first derivative score at p=1/2."""
    x = 2 * k - n
    return 2 * x


def third_weight(n: int, k: int) -> int:
    """Weight whose expectation against D equals d^3 E_p[D]/dp^3 at p=1/2."""
    x = 2 * k - n
    return 8 * (x**3 - (3 * n - 2) * x)


def aggregate(n: int, samples: Iterable[tuple[int, float]]) -> ScoreMoments:
    """Aggregate `(occupied_count, observable)` samples.

    Samples must be drawn from the p=1/2 Bernoulli ensemble.  The function
    treats them as iid; callers using correlated chains must handle effective
    sample sizes separately.
    """

    rows = list(samples)
    if not rows:
        raise ValueError("no samples")

    d_values = [float(d) for _k, d in rows]
    d1_values = [first_score(n, k) * float(d) for k, d in rows]
    d3_values = [third_weight(n, k) * float(d) for k, d in rows]
    count = len(rows)

    mean_d = fsum(d_values) / count
    d1 = fsum(d1_values) / count
    d3 = fsum(d3_values) / count
    if d1 == 0.0:
        raise ZeroDivisionError("estimated first derivative is zero")
    return ScoreMoments(
        samples=count,
        mean_d=mean_d,
        first_derivative=d1,
        third_derivative=d3,
        kappa3=d3 / d1**3,
    )


def _self_test() -> None:
    # For D(K)=K/N, E_p[D]=p, so M'(1/2)=1 and M'''(1/2)=0 exactly.
    # Exhaustively average all 2^N Bernoulli configurations through the
    # binomial multiplicities rather than enumerate bit strings.
    from math import comb

    for n in range(1, 10):
        weighted: list[tuple[int, float]] = []
        # Repeat according to binomial multiplicity only in this tiny self-test.
        for k in range(n + 1):
            weighted.extend([(k, k / n)] * comb(n, k))
        result = aggregate(n, weighted)
        assert abs(result.first_derivative - 1.0) < 1e-12
        assert abs(result.third_derivative) < 1e-10


if __name__ == "__main__":
    _self_test()
    print("score identities: self-test passed")
