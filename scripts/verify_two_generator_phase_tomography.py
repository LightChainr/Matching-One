#!/usr/bin/env python3
"""Exact Gaussian spin-4 phase arithmetic for the two-generator protocol.

This script is deterministic and outcome-free.  It verifies the phase nodes,
Gaussian-character composition, and the norm-5/norm-10 two-row quadrature
reconstruction used by predictions/norm4_two_generator_transfer_20260829.yaml.
"""
from __future__ import annotations

import json
from fractions import Fraction


def mul(z: tuple[int, int], w: tuple[int, int]) -> tuple[int, int]:
    a, b = z
    c, d = w
    return a * c - b * d, a * d + b * c


def fourth(z: tuple[int, int]) -> tuple[int, int]:
    square = mul(z, z)
    return mul(square, square)


def chi4(z: tuple[int, int]) -> tuple[Fraction, Fraction]:
    a, b = z
    real, imag = fourth(z)
    norm = a * a + b * b
    denom = norm * norm
    return Fraction(real, denom), Fraction(imag, denom)


def ftext(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main() -> int:
    names = {
        "1+i": (1, 1),
        "1-i": (1, -1),
        "2i": (0, 2),
        "2+i": (2, 1),
        "2-i": (2, -1),
        "3+i": (3, 1),
        "3-i": (3, -1),
    }
    table = {
        name: {
            "norm": z[0] * z[0] + z[1] * z[1],
            "real": ftext(chi4(z)[0]),
            "imag": ftext(chi4(z)[1]),
        }
        for name, z in names.items()
    }

    assert table["1+i"]["imag"] == "0"
    assert table["1-i"]["imag"] == "0"
    assert table["2i"] == {"norm": 4, "real": "1", "imag": "0"}
    assert chi4(mul(names["1+i"], names["2-i"])) == chi4(names["3+i"])

    # Phase tomography rows for m5=2+i and m10=3+i, after radial normalization.
    a = Fraction(-7, 25)
    b = Fraction(24, 25)
    c = Fraction(7, 25)
    d = Fraction(24, 25)
    determinant = a * d - b * c
    assert determinant == Fraction(-336, 625)

    # A^T A is diagonal because the cosine columns have opposite signs while
    # the sine columns agree.  The exact 2-norm condition number is 24/7.
    ata_cos = a * a + c * c
    ata_sin = b * b + d * d
    assert ata_cos == Fraction(98, 625)
    assert ata_sin == Fraction(1152, 625)

    # Symbolic reconstruction check using arbitrary rational quadratures.
    cosine = Fraction(7, 13)
    sine = Fraction(-5, 11)
    y5 = a * cosine + b * sine
    y10 = c * cosine + d * sine
    cosine_back = Fraction(25, 14) * (y10 - y5)
    sine_back = Fraction(25, 48) * (y10 + y5)
    assert cosine_back == cosine
    assert sine_back == sine

    payload = {
        "schema": "matching-one/two-generator-phase-tomography/v1",
        "chi4": table,
        "semigroup": {
            "identity": "(1+i)(2-i)=3+i",
            "character_composes": True,
        },
        "norm5_norm10_phase_matrix": [
            ["-7/25", "24/25"],
            ["7/25", "24/25"],
        ],
        "determinant": "-336/625",
        "condition_number_2norm": "24/7",
        "reconstruction": {
            "cosine": "(25/14)*(Y10-Y5)",
            "sine": "(25/48)*(Y10+Y5)",
        },
        "phase_nodes": {
            "norm2_1_plusminus_i_sine": "0",
            "norm4_2i_sine": "0",
        },
    }
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
