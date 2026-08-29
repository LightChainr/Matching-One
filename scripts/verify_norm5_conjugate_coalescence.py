#!/usr/bin/env python3
"""Verify exact Gaussian coalescence and same-N harmonic interpolation targets."""
from __future__ import annotations

from fractions import Fraction
from math import gcd


def mul(z: tuple[int, int], w: tuple[int, int]) -> tuple[int, int]:
    a, b = z
    c, d = w
    return a * c - b * d, a * d + b * c


def norm(z: tuple[int, int]) -> int:
    a, b = z
    return a * a + b * b


def d4_canonical(z: tuple[int, int]) -> tuple[int, int]:
    a, b = map(abs, z)
    return max(a, b), min(a, b)


def smith_for_gaussian(z: tuple[int, int]) -> tuple[int, int]:
    a, b = map(abs, z)
    d1 = gcd(a, b)
    return d1, norm(z) // d1


def harmonic_cos(z: tuple[int, int], spin: int) -> Fraction:
    if spin % 2:
        raise ValueError("spin must be even")
    a, b = z
    x, y = 1, 0
    for _ in range(spin):
        x, y = x * a - y * b, x * b + y * a
    return Fraction(x, norm(z) ** (spin // 2))


def harmonic_sin(z: tuple[int, int], spin: int = 4) -> Fraction:
    if spin % 2:
        raise ValueError("spin must be even")
    a, b = z
    x, y = 1, 0
    for _ in range(spin):
        x, y = x * a - y * b, x * b + y * a
    return Fraction(y, norm(z) ** (spin // 2))


def interpolation_weights(
    A: tuple[int, int], B: tuple[int, int], C: tuple[int, int], spin: int
) -> tuple[Fraction, Fraction]:
    cA, cB, cC = (harmonic_cos(z, spin) for z in (A, B, C))
    if cA == cB:
        raise ValueError("source angles do not identify this harmonic")
    wA = (cC - cB) / (cA - cB)
    return wA, 1 - wA


def verified_payload() -> dict:
    h_minus = (2, -1)
    h_plus = (2, 1)

    cases = {
        325: {
            "parents": ((8, 1), (7, 4)),
            "observed_h": h_minus,
            "conjugate_h": h_plus,
            "A": (17, 6),
            "B": (18, 1),
            "C": (15, 10),
            "smith": (5, 65),
            "sine_ratio": Fraction(6, 11),
            "weights": {
                4: (Fraction(11, 5), Fraction(-6, 5)),
                8: (Fraction(22517, 44795), Fraction(22278, 44795)),
                12: (Fraction(363263, 7144145), Fraction(6780882, 7144145)),
            },
        },
        425: {
            "parents": ((9, 2), (7, 6)),
            "observed_h": h_plus,
            "conjugate_h": h_minus,
            "A": (16, 13),
            "B": (19, 8),
            "C": (20, 5),
            "smith": (5, 85),
            "sine_ratio": Fraction(33, 13),
            "weights": {
                4: (Fraction(-13, 20), Fraction(33, 20)),
                8: (Fraction(89531, 242420), Fraction(152889, 242420)),
                12: (Fraction(181189, 68620), Fraction(-112569, 68620)),
            },
        },
    }

    output = {}
    for N, row in cases.items():
        p1, p2 = row["parents"]
        obs1 = mul(p1, row["observed_h"])
        obs2 = mul(p2, row["observed_h"])
        con1 = mul(p1, row["conjugate_h"])
        con2 = mul(p2, row["conjugate_h"])
        assert d4_canonical(obs1) == row["A"]
        assert d4_canonical(obs2) == row["B"]
        assert d4_canonical(con1) == row["C"]
        assert d4_canonical(con2) == row["C"]
        assert norm(row["A"]) == norm(row["B"]) == norm(row["C"]) == N
        assert smith_for_gaussian(row["C"]) == row["smith"]
        assert harmonic_sin(p1) / harmonic_sin(p2) == row["sine_ratio"]
        for spin, expected in row["weights"].items():
            assert interpolation_weights(row["A"], row["B"], row["C"], spin) == expected
        output[N] = {
            "products": [obs1, obs2, con1, con2],
            "canonical_children": [
                d4_canonical(obs1), d4_canonical(obs2),
                d4_canonical(con1), d4_canonical(con2),
            ],
            "C_smith": row["smith"],
            "parent_sine_ratio": str(row["sine_ratio"]),
            "weights": {
                f"H{spin}": [str(weight) for weight in row["weights"][spin]]
                for spin in (4, 8, 12)
            },
        }
    return output


if __name__ == "__main__":
    for N, row in verified_payload().items():
        print(N, row)
