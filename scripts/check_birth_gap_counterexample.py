#!/usr/bin/env python3
"""One cheap exact check for the birth-gap note; no percolation simulation."""

from fractions import Fraction as F
from math import comb
import json


def average(pairs, fn):
    return sum((fn(*pair) for pair in pairs), F(0)) / len(pairs)


def beta_cdf(alpha, beta, x):
    n = alpha + beta - 1
    return sum((comb(n, k) * x**k * (1 - x)**(n - k)
                for k in range(alpha, n + 1)), F(0))


def beta_soft_direct(n, d, a, delta):
    """Integrate (1-a*t/delta) times beta density, expanded as a polynomial."""
    if d == 0:
        return F(1)
    x = min(delta / a, F(1))
    return d * comb(n, d) * sum((
        (-1)**k * comb(n - d, k)
        * (x**(d + k) / (d + k)
           - a / delta * x**(d + k + 1) / (d + k + 1))
        for k in range(n - d + 1)), F(0))


def main():
    counterexamples = []
    for n in (2, 4, 16, 64):
        a = [(F(0), 1 + F(1, n)), (F(1), F(2))]
        b = [(F(0), F(2)), (F(1), 1 + F(1, n))]
        assert sorted(s for s, _ in a) == sorted(s for s, _ in b)
        assert sorted(t for _, t in a) == sorted(t for _, t in b)
        assert all(s < t for s, t in a + b)
        endpoints = sorted({x for pair in a + b for x in pair})
        for t in endpoints:
            for rank in range(3):
                assert average(a, lambda s, u: F(int((s <= t) + (u <= t) == rank))) == average(
                    b, lambda s, u: F(int((s <= t) + (u <= t) == rank)))
        assert average(a, lambda s, t: t - s) == 1 + F(1, 2*n)
        assert average(a, lambda s, t: t - s) == average(b, lambda s, t: t - s)
        delta = F(1, 4)
        z_a = average(a, lambda s, t: max(F(0), 1 - (t-s)/delta))
        z_b = average(b, lambda s, t: max(F(0), 1 - (t-s)/delta))
        assert z_a == 0
        assert z_b == max(F(0), F(1, 2) * (1 - F(4, n)))
        counterexamples.append({"n": n, "Z_A_quarter": str(z_a), "Z_B_quarter": str(z_b)})

    n = 9
    hist = {0: F(1, 5), 1: F(1, 5), 3: F(2, 5), 7: F(1, 5)}
    overlap = lambda m: sum((h * max(d - m, 0) for d, h in hist.items()), F(0))
    for d in range(1, n + 1):
        assert overlap(d-1) - 2*overlap(d) + overlap(d+1) == hist.get(d, 0)
    for m in range(1, n + 1):
        assert 1 - (overlap(0)-overlap(m))/m == sum(
            (h * max(F(0), 1-F(d, m)) for d, h in hist.items()), F(0))
    assert 1 - (overlap(0)-overlap(1)) == hist[0]

    scale, delta = F(2), F(3, 4)
    x = min(delta/scale, F(1))
    soft = F(0)
    for d, weight in hist.items():
        closed = F(1) if d == 0 else (
            beta_cdf(d, n+1-d, x)
            - scale/delta * F(d, n+1) * beta_cdf(d+1, n+1-d, x))
        assert closed == beta_soft_direct(n, d, scale, delta)
        soft += weight * closed
    print(json.dumps({"status": "passed_exact_fraction_checks", "counterexamples": counterexamples,
                      "beta_mixture_soft_mass": str(soft), "direct_atom": str(hist[0])}, indent=2))


if __name__ == "__main__":
    main()
