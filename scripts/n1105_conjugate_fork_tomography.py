#!/usr/bin/env python3
"""Exact Gaussian arithmetic for the N=221 -> N=1105 conjugate norm-5 fork.

No percolation outcome is read.  The purpose is to expose a hidden genealogy in
the existing N=1105 four-orientation design:

    parent N=221: (14,5), (11,10)
    multiply by 2+i -> (24,23), (31,12)
    multiply by 2-i -> (33,4),  (32,9)

The two child arms have the same area.  For a pure H4 orientation response the
radial amplitude therefore cancels between arms and the exact child-contrast
ratio is -33/13, equivalently

    13 D_plus + 33 D_minus = 0.

This is independent of the H4 radial exponent and of the nonuniversal amplitude.
"""

from __future__ import annotations

import json
from fractions import Fraction

PARENT_N = 221
CHILD_N = 1105
PARENT = ((14, 5), (11, 10))
M_PLUS = (2, 1)
M_MINUS = (2, -1)
EXPECTED_PLUS = ((24, 23), (31, 12))
EXPECTED_MINUS = ((33, 4), (32, 9))
ISSUE74_ORIENTATIONS = ((33, 4), (32, 9), (31, 12), (24, 23))


def gaussian_mul(z: tuple[int, int], w: tuple[int, int]) -> tuple[int, int]:
    a, b = z
    c, d = w
    return a * c - b * d, a * d + b * c


def gaussian_pow(z: tuple[int, int], power: int) -> tuple[int, int]:
    result = (1, 0)
    base = z
    while power:
        if power & 1:
            result = gaussian_mul(result, base)
        base = gaussian_mul(base, base)
        power >>= 1
    return result


def norm(z: tuple[int, int]) -> int:
    return z[0] * z[0] + z[1] * z[1]


def canonical_orientation(z: tuple[int, int]) -> tuple[int, int]:
    a, b = abs(z[0]), abs(z[1])
    return (a, b) if a >= b else (b, a)


def chi4(z: tuple[int, int]) -> tuple[Fraction, Fraction]:
    """Return ((z/|z|)^4).real/.imag exactly."""
    a, b = z
    n = norm(z)
    return (
        Fraction(a**4 - 6 * a * a * b * b + b**4, n * n),
        Fraction(4 * a * b * (a * a - b * b), n * n),
    )


def cos_harmonic(z: tuple[int, int], spin: int) -> Fraction:
    if spin <= 0 or spin % 4:
        raise ValueError("spin must be a positive multiple of 4")
    real, _imag = gaussian_pow(z, spin)
    return Fraction(real, norm(z) ** (spin // 2))


def contrast(pair: tuple[tuple[int, int], tuple[int, int]], spin: int) -> Fraction:
    return cos_harmonic(pair[0], spin) - cos_harmonic(pair[1], spin)


def child_pair(multiplier: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    return tuple(gaussian_mul(z, multiplier) for z in PARENT)  # type: ignore[return-value]


def canonical_child_pair(multiplier: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    return tuple(canonical_orientation(z) for z in child_pair(multiplier))  # type: ignore[return-value]


def fork_ratio(spin: int) -> Fraction:
    plus = contrast(child_pair(M_PLUS), spin)
    minus = contrast(child_pair(M_MINUS), spin)
    if not minus:
        raise ZeroDivisionError("minus-arm contrast vanishes")
    return plus / minus


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def render() -> dict[str, object]:
    if norm(PARENT[0]) != PARENT_N or norm(PARENT[1]) != PARENT_N:
        raise AssertionError("parent norms changed")
    plus = canonical_child_pair(M_PLUS)
    minus = canonical_child_pair(M_MINUS)
    if plus != EXPECTED_PLUS or minus != EXPECTED_MINUS:
        raise AssertionError("conjugate fork genealogy changed")
    if any(norm(z) != CHILD_N for z in (*plus, *minus)):
        raise AssertionError("child norm changed")
    if set((*plus, *minus)) != set(ISSUE74_ORIENTATIONS):
        raise AssertionError("fork does not reproduce the Issue #74 four-angle set")

    parent_re, parent_im = (
        chi4(PARENT[0])[0] - chi4(PARENT[1])[0],
        chi4(PARENT[0])[1] - chi4(PARENT[1])[1],
    )
    if parent_re != Fraction(57600, 48841):
        raise AssertionError("unexpected parent DeltaCos4")
    if parent_im != Fraction(38640, 48841):
        raise AssertionError("unexpected parent DeltaSin4")
    if parent_im / parent_re != Fraction(161, 240):
        raise AssertionError("unexpected parent spin-4 quadrature ratio")

    plus_re = contrast(child_pair(M_PLUS), 4)
    minus_re = contrast(child_pair(M_MINUS), 4)
    plus_transfer = plus_re / parent_re
    minus_transfer = minus_re / parent_re
    if plus_transfer != Fraction(-231, 250):
        raise AssertionError("unexpected 2+i H4 angular transfer")
    if minus_transfer != Fraction(91, 250):
        raise AssertionError("unexpected 2-i H4 angular transfer")
    if plus_re / minus_re != Fraction(-33, 13):
        raise AssertionError("unexpected H4 fork ratio")
    if 13 * plus_re + 33 * minus_re:
        raise AssertionError("H4 fork null failed")

    harmonics = {}
    for spin in (4, 8, 12, 16, 20):
        p = contrast(child_pair(M_PLUS), spin)
        m = contrast(child_pair(M_MINUS), spin)
        harmonics[f"H{spin}"] = {
            "plus_contrast": fraction_text(p),
            "minus_contrast": fraction_text(m),
            "plus_over_minus": fraction_text(p / m),
        }

    return {
        "schema": "matching-one/N1105-conjugate-norm5-fork/v1",
        "classification": "exact geometry/design; no target data",
        "parent": {
            "N": PARENT_N,
            "orientations": [list(z) for z in PARENT],
            "DeltaCos4": fraction_text(parent_re),
            "DeltaSin4": fraction_text(parent_im),
            "DeltaSin4_over_DeltaCos4": "161/240",
        },
        "multipliers": {
            "plus": {"m": "2+i", "chi4": ["-7/25", "24/25"]},
            "minus": {"m": "2-i", "chi4": ["-7/25", "-24/25"]},
        },
        "children": {
            "N": CHILD_N,
            "plus_genealogy": [list(z) for z in EXPECTED_PLUS],
            "minus_genealogy": [list(z) for z in EXPECTED_MINUS],
            "union_equals_issue74_four_angles": True,
        },
        "H4": {
            "plus_over_parent_DeltaCos4": "-231/250",
            "minus_over_parent_DeltaCos4": "91/250",
            "plus_over_minus_same_N": "-33/13",
            "amplitude_free_null": "13*D_plus + 33*D_minus = 0",
            "root_plus_over_parent_if_quartic_character": "-231/6250",
            "root_minus_over_parent_if_quartic_character": "91/6250",
        },
        "harmonic_fork_directions": harmonics,
        "interpretation": (
            "The same N=1105 four-angle data can be viewed both as an H0/H4/H8/H12 "
            "projector and as two conjugate degree-5 descendants of one N=221 parent."
        ),
    }


def main() -> int:
    print(json.dumps(render(), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
