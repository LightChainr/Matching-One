#!/usr/bin/env python3
"""Exact arithmetic behind the minimal N=1105 four-orientation torus.

For a primitive Gaussian square torus represented by a+ib with gcd(a,b)=1,
N=a^2+b^2.  Modulo square D4 symmetry, the number of primitive orientations is
2^(k-1), where k is the number of distinct prime divisors p == 1 (mod 4).
A primitive representation forbids p == 3 (mod 4) divisors and 4 | N.

Hence four D4-inequivalent primitive orientations require k>=3.  The smallest
possible N is 5*13*17=1105.
"""

from __future__ import annotations

import math
from typing import List, Tuple


def factorize(n: int) -> List[Tuple[int, int]]:
    if n <= 0:
        raise ValueError("n must be positive")
    out: List[Tuple[int, int]] = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            out.append((p, e))
        p = 3 if p == 2 else p + 2
    if n > 1:
        out.append((n, 1))
    return out


def primitive_first_octant_representations(n: int) -> List[Tuple[int, int]]:
    """Return primitive a>=b>=0 representatives of a^2+b^2=n."""
    out: List[Tuple[int, int]] = []
    for b in range(math.isqrt(n) + 1):
        a2 = n - b * b
        a = math.isqrt(a2)
        if a * a == a2 and a >= b and a > 0 and math.gcd(a, b) == 1:
            out.append((a, b))
    return out


def primitive_orientation_count_formula(n: int) -> int:
    """Count primitive D4 orientation orbits using Gaussian factorization."""
    factors = factorize(n)
    for p, e in factors:
        if p % 4 == 3:
            return 0
        if p == 2 and e >= 2:
            return 0
    k = sum(1 for p, _e in factors if p % 4 == 1)
    if k == 0:
        # N=1 or 2 has the axis/diagonal exceptional primitive orbit.
        return 1 if n in (1, 2) else 0
    return 1 << (k - 1)


def minimal_n_for_at_least_four_orientations(limit: int = 1105) -> int:
    for n in range(1, limit + 1):
        if primitive_orientation_count_formula(n) >= 4:
            return n
    raise ValueError("no four-orientation primitive torus within limit")


def _self_test() -> None:
    for n in range(1, 1106):
        brute = len(primitive_first_octant_representations(n))
        formula = primitive_orientation_count_formula(n)
        assert brute == formula, (n, brute, formula)
    assert minimal_n_for_at_least_four_orientations() == 1105
    assert factorize(1105) == [(5, 1), (13, 1), (17, 1)]
    assert primitive_first_octant_representations(1105) == [
        (33, 4),
        (32, 9),
        (31, 12),
        (24, 23),
    ]


if __name__ == "__main__":
    _self_test()
    print("minimal N with >=4 primitive D4 orientation orbits: 1105")
    print("factorization: 5 * 13 * 17")
    print("orientations:", primitive_first_octant_representations(1105))
