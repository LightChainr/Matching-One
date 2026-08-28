#!/usr/bin/env python3
"""Combinatorial check of the c=0, h=5/8 thermal spin-4 tower.

For the h=5/8 Verma module the first singular vector is at level 2.
The next Kac singular level is 10, so through level 9 the quotient
character is p(n)-p(n-2).  For h>0, quasiprimary counts below that
second singular level are obtained by subtracting the L_-1 image:
q(n)=d(n)-d(n-1).

A non-redundant bulk spin-4 descendant built from chiral quasiprimaries
at levels (n,m) requires n-m=4 and q(n),q(m)>0.
"""

from functools import lru_cache
from fractions import Fraction


@lru_cache(None)
def partitions(n: int, largest: int | None = None) -> int:
    if n < 0:
        return 0
    if n == 0:
        return 1
    if largest is None or largest > n:
        largest = n
    if largest <= 0:
        return 0
    return partitions(n, largest - 1) + partitions(n - largest, largest)


def quotient_dim(level: int) -> int:
    return partitions(level) - partitions(level - 2)


def quasiprimary_count(level: int) -> int:
    if level == 0:
        return 1
    return quotient_dim(level) - quotient_dim(level - 1)


def kac_singular_levels(limit_r: int = 20, limit_s: int = 30) -> list[int]:
    # c=0: h_{r,s}=((3r-2s)^2-1)/24.  h=5/8 iff |3r-2s|=4.
    levels = {
        r * s
        for r in range(1, limit_r + 1)
        for s in range(1, limit_s + 1)
        if abs(3 * r - 2 * s) == 4
    }
    return sorted(levels)


def main() -> int:
    singular = kac_singular_levels()
    assert singular[:3] == [2, 10, 16]

    counts = {level: quasiprimary_count(level) for level in range(10)}
    assert counts[0] == 1
    assert counts[1] == 0
    assert counts[2] == 0
    assert counts[3] > 0
    assert counts[4] > 0
    assert counts[7] > 0

    pairs = []
    for m in range(10):
        n = m + 4
        if n >= 10:
            break
        if counts[m] > 0 and counts[n] > 0:
            total_level = n + m
            x = Fraction(5, 4) + total_level
            residual_power_L = x - 2
            root_power_L = x - Fraction(5, 4)
            pairs.append((n, m, total_level, x, residual_power_L, root_power_L))

    # Below the second chiral singular level, the first two non-redundant
    # spin-4 bulk combinations are (4,0) and (7,3).
    assert pairs[0][:3] == (4, 0, 4)
    assert pairs[1][:3] == (7, 3, 10)
    assert pairs[0][3] == Fraction(21, 4)
    assert pairs[1][3] == Fraction(45, 4)

    print("c=0 thermal h=hbar=5/8")
    print("Kac singular levels:", singular[:6])
    print("level : quotient_dim : quasiprimary_count")
    for level in range(10):
        print(level, quotient_dim(level), counts[level])
    print("spin-4 quasiprimary pairs (n,m):")
    for row in pairs:
        n, m, total, x, residual, root = row
        print(
            f"({n},{m}) total={total} x={x} "
            f"M(pc)~L^-{residual} root~L^-{root}"
        )
    print("leading-to-next same-family relative correction: L^-6 = N^-3")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
