#!/usr/bin/env python3
"""Verify exact same-parent norm-5 conjugate coalescence and harmonic weights."""
from __future__ import annotations

import json
from fractions import Fraction
from math import gcd

Gaussian = tuple[int, int]


def mul(z: Gaussian, w: Gaussian) -> Gaussian:
    a, b = z
    c, d = w
    return a * c - b * d, a * d + b * c


def norm(z: Gaussian) -> int:
    a, b = z
    return a * a + b * b


def power(z: Gaussian, exponent: int) -> Gaussian:
    result = (1, 0)
    for _ in range(exponent):
        result = mul(result, z)
    return result


def d4_canonical(z: Gaussian) -> Gaussian:
    a, b = sorted((abs(z[0]), abs(z[1])), reverse=True)
    return a, b


def smith_invariants(z: Gaussian) -> tuple[int, int]:
    divisor = gcd(abs(z[0]), abs(z[1]))
    return divisor, norm(z) // divisor


def cosine_harmonic(z: Gaussian, spin: int) -> Fraction:
    if spin <= 0 or spin % 4:
        raise ValueError("spin must be a positive multiple of four")
    real, _ = power(z, spin)
    return Fraction(real, norm(z) ** (spin // 2))


def affine_weights(a: Gaussian, b: Gaussian, c: Gaussian, spin: int) -> tuple[Fraction, Fraction]:
    ca = cosine_harmonic(a, spin)
    cb = cosine_harmonic(b, spin)
    cc = cosine_harmonic(c, spin)
    if ca == cb:
        raise ValueError("A and B do not separate this harmonic")
    weight_a = Fraction(cc - cb, ca - cb)
    return weight_a, 1 - weight_a


def ftext(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


DESIGNS = {
    325: {
        "parents": ((8, 1), (7, 4)),
        "observed_multiplier": (2, -1),
        "conjugate_multiplier": (2, 1),
        "A": (17, 6),
        "B": (18, 1),
        "C": (15, 10),
        "expected_weights": {
            4: (Fraction(11, 5), Fraction(-6, 5)),
            8: (Fraction(22517, 44795), Fraction(22278, 44795)),
            12: (Fraction(363263, 7144145), Fraction(6780882, 7144145)),
        },
        "h4_integer_residual_C_A_B": (5, -11, 6),
        "h4_conjugate_ratio": Fraction(6, 11),
    },
    425: {
        "parents": ((9, 2), (7, 6)),
        "observed_multiplier": (2, 1),
        "conjugate_multiplier": (2, -1),
        "A": (16, 13),
        "B": (19, 8),
        "C": (20, 5),
        "expected_weights": {
            4: (Fraction(-13, 20), Fraction(33, 20)),
            8: (Fraction(89531, 242420), Fraction(152889, 242420)),
            12: (Fraction(181189, 68620), Fraction(-112569, 68620)),
        },
        "h4_integer_residual_C_A_B": (20, 13, -33),
        "h4_conjugate_ratio": Fraction(33, 13),
    },
}


def verify_design(n: int, design: dict[str, object]) -> dict[str, object]:
    parents = design["parents"]
    observed = design["observed_multiplier"]
    conjugate = design["conjugate_multiplier"]
    a = design["A"]
    b = design["B"]
    c = design["C"]
    assert isinstance(parents, tuple) and isinstance(observed, tuple) and isinstance(conjugate, tuple)
    assert isinstance(a, tuple) and isinstance(b, tuple) and isinstance(c, tuple)

    observed_children = tuple(d4_canonical(mul(parent, observed)) for parent in parents)
    conjugate_children = tuple(d4_canonical(mul(parent, conjugate)) for parent in parents)
    assert observed_children == (a, b)
    assert conjugate_children == (c, c)
    assert norm(a) == norm(b) == norm(c) == n
    assert smith_invariants(a) == (1, n)
    assert smith_invariants(b) == (1, n)
    assert smith_invariants(c) == (5, n // 5)

    harmonics: dict[str, object] = {}
    expected = design["expected_weights"]
    assert isinstance(expected, dict)
    for spin in (4, 8, 12):
        weights = affine_weights(a, b, c, spin)
        assert weights == expected[spin]
        harmonics[f"H{spin}"] = {
            "cos_A": ftext(cosine_harmonic(a, spin)),
            "cos_B": ftext(cosine_harmonic(b, spin)),
            "cos_C": ftext(cosine_harmonic(c, spin)),
            "weight_A": ftext(weights[0]),
            "weight_B": ftext(weights[1]),
        }

    ratio = Fraction(cosine_harmonic(c, 4) - cosine_harmonic(a, 4),
                     cosine_harmonic(c, 4) - cosine_harmonic(b, 4))
    assert ratio == design["h4_conjugate_ratio"]
    residual = design["h4_integer_residual_C_A_B"]
    assert isinstance(residual, tuple) and sum(residual) == 0

    return {
        "N": n,
        "parents": [list(parent) for parent in parents],
        "observed_multiplier": list(observed),
        "conjugate_multiplier": list(conjugate),
        "observed_children_D4": [list(child) for child in observed_children],
        "conjugate_children_D4": [list(child) for child in conjugate_children],
        "A": list(a),
        "B": list(b),
        "C": list(c),
        "smith_A": list(smith_invariants(a)),
        "smith_B": list(smith_invariants(b)),
        "smith_C": list(smith_invariants(c)),
        "harmonics": harmonics,
        "H4_integer_residual_C_A_B": list(residual),
        "H4_conjugate_ratio": ftext(ratio),
    }


def main() -> int:
    payload = {
        "schema": "matching-one/norm5-conjugate-coalescence/v1",
        "interpretation": (
            "same-parent conjugate norm-5 branches coalesce to one noncyclic D4 geometry; "
            "the H4 test is same-N and uses no radial exponent or fitted amplitude"
        ),
        "designs": [verify_design(n, design) for n, design in sorted(DESIGNS.items())],
    }
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
