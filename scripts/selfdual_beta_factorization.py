#!/usr/bin/env python3
"""Exact self-dual reliability factorization checks.

For complement-odd M, write M=(2p-1) H(p(1-p)).  This script verifies the
N=10 and post-target PR #152 N=26 decompositions using integer arithmetic only.
"""

from __future__ import annotations

import json
from math import comb


def central_binomial_prefix(s: int) -> list[int]:
    return [comb(2 * j, j) for j in range(s)]


def expand_m_from_h(h: list[int]) -> list[int]:
    """Expand (2p-1) sum_j h_j [p(1-p)]^j in the power basis."""
    out = [0] * (2 * (len(h) - 1) + 2)
    for j, c in enumerate(h):
        # u^j = p^j (1-p)^j = sum_r (-1)^r C(j,r) p^(j+r)
        for r in range(j + 1):
            base = c * ((-1) ** r) * comb(j, r)
            degree = j + r
            out[degree] -= base          # multiply by -1
            out[degree + 1] += 2 * base  # multiply by 2p
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def minimum_success_count_from_h(s: int, c_s: int) -> int:
    numerator = comb(2 * s, s) - c_s
    assert numerator % 2 == 0
    return numerator // 2


def main() -> None:
    h10 = [1, 2, 6]
    m10 = expand_m_from_h(h10)
    expected10 = [-1, 0, 0, 20, -30, 12]
    assert m10 == expected10
    assert h10 == central_binomial_prefix(3)

    h26 = [1, 2, 6, 20, 70, 96, 170, 260, 260, 78]
    m26 = expand_m_from_h(h26)
    expected26 = [
        -1, 0, 0, 0, 0, 156, -338, 260, -260, -338,
        1144, 3536, -13702, 15628, -3016, -10088,
        11492, -5798, 1482, -156,
    ]
    assert m26 == expected26
    assert h26[:5] == central_binomial_prefix(5)
    assert minimum_success_count_from_h(5, h26[5]) == 78
    assert all(c > 0 for c in h26)

    beta5_h = central_binomial_prefix(5)
    # The deformation begins at u^5 with the post-target N=26 coefficients.
    deformation = h26[5:]
    assert deformation == [96, 170, 260, 260, 78]

    print(json.dumps({
        "N10_H": h10,
        "N26_H": h26,
        "N26_universal_prefix": central_binomial_prefix(5),
        "N26_first_free_coefficient": h26[5],
        "N26_minimum_success_count": minimum_success_count_from_h(5, h26[5]),
        "N26_deformation_after_u5": deformation,
        "N26_H_coefficients_positive": True,
    }, indent=2))


if __name__ == "__main__":
    main()
