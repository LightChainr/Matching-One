#!/usr/bin/env python3
"""Exact verifier for the P429 parallel-gadget predictive-state lower bound.

The only scientific inputs are the exact PR #435 base-gadget counts:

    safe-subset counts b=(1,7,18,20,8,0,0,0,0)
    A safe-successor exit counts: 1x3, 2x2, 3x2
    B safe-successor exit counts: 1x1, 2x6

No random sampling is used.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import comb

BASE_SAFE = (1, 7, 18, 20, 8, 0, 0, 0, 0)
EXIT_A = (1, 1, 1, 2, 2, 3, 3)
EXIT_B = (1, 2, 2, 2, 2, 2, 2)


def convolve(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return tuple(out)


def safe_counts(k: int) -> tuple[int, ...]:
    if k < 1:
        raise ValueError("k must be >=1")
    out = (1,)
    for _ in range(k):
        out = convolve(out, BASE_SAFE)
    return out


def survival_law(k: int) -> tuple[Fraction, ...]:
    counts = safe_counts(k)
    n = 8 * k
    return tuple(Fraction(c, comb(n, m)) for m, c in enumerate(counts))


def fork_probability(k: int, a_count: int) -> Fraction:
    if k < 1 or not 0 <= a_count <= k:
        raise ValueError("require k>=1 and 0<=a_count<=k")
    numerator = 343 * k**3 - 182 * k**2 + 25 * k + 4 * a_count
    denominator = 8 * k * (8 * k - 1) ** 2
    return Fraction(numerator, denominator)


def fork_probability_by_successor_sum(k: int, a_count: int) -> Fraction:
    exits = [EXIT_A] * a_count + [EXIT_B] * (k - a_count)
    total = 0
    for local in exits:
        for x in local:
            safe_second = 7 * k - x
            total += safe_second * safe_second
    return Fraction(total, 8 * k * (8 * k - 1) ** 2)


def verify(k_max: int = 12) -> None:
    assert sum(EXIT_A) == sum(EXIT_B) == 13
    assert sum(x * x for x in EXIT_A) == 29
    assert sum(x * x for x in EXIT_B) == 25

    assert fork_probability(1, 0) == Fraction(93, 196)
    assert fork_probability(1, 1) == Fraction(95, 196)
    assert fork_probability(1, 1) - fork_probability(1, 0) == Fraction(1, 98)

    for k in range(1, k_max + 1):
        law = survival_law(k)
        assert law[0] == 1
        for a in range(k + 1):
            direct = fork_probability_by_successor_sum(k, a)
            closed = fork_probability(k, a)
            assert direct == closed
            if a < k:
                gap = fork_probability(k, a + 1) - closed
                assert gap == Fraction(1, 2 * k * (8 * k - 1) ** 2)
                assert gap > 0

        # The complete unbranched survival law is composition-independent:
        # the only input is BASE_SAFE, so there is one law for all a=0..k.
        assert len(law) == 8 * k + 1


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k-max", type=int, default=12)
    parser.add_argument("--show-k", type=int)
    args = parser.parse_args()
    verify(args.k_max)
    print(f"verified exact lower-bound identities for k=1..{args.k_max}")
    if args.show_k is not None:
        k = args.show_k
        print(f"k={k}; future_vertices={8*k}; predictive_classes>={k+1}")
        print("fork probabilities:")
        for a in range(k + 1):
            print(f"  a={a}: {fork_probability(k, a)}")


if __name__ == "__main__":
    main()
