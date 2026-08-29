#!/usr/bin/env python3
"""Derive exact weight-4 Hecke ratios at the square CM point.

The calculation uses only the modular transformation law of a weight-4
form and the normalized Hecke eigenvalue of E4 at p=2.  It deliberately
does not assert that a lattice percolation observable equals E4.
"""

from __future__ import annotations

import json
from fractions import Fraction


def exact_ratios() -> dict[str, Fraction]:
    """Return the exact holomorphic and area-normalized E4 ratios.

    Put x=E4(2i)/E4(i).  The weight-4 modular transformations give

      E4(i/2)/E4(i)       = 16 x,
      E4((1+i)/2)/E4(i)   = -4.

    Since T_2 E4=(1+2^3)E4, evaluation at i gives

      8x + (16x-4)/2 = 9,

    hence x=11/16.  The area-normalized spin-4 shape is
    E4hat(tau)=(Im tau)^2 E4(tau).
    """

    weight = 4
    hecke_eigenvalue = Fraction(1 + 2 ** (weight - 1), 1)
    e4_i_over_2_over_e4_2i = Fraction(2**weight, 1)
    e4_diagonal_over_e4_i = Fraction(-4, 1)

    # T_2 f(i) / f(i) = 2^(k-1) x
    #                    + 1/2 [2^k x + f((i+1)/2)/f(i)].
    coefficient_of_x = Fraction(2 ** (weight - 1), 1) + Fraction(
        2**weight, 2
    )
    constant_term = e4_diagonal_over_e4_i / 2
    e4_2i_over_e4_i = (
        hecke_eigenvalue - constant_term
    ) / coefficient_of_x
    e4_i_over_2_over_e4_i = (
        e4_i_over_2_over_e4_2i * e4_2i_over_e4_i
    )

    # E4hat=(Im tau)^2 E4.  Im(2i)=2 and both half-moduli have Im=1/2.
    e4hat_2i_over_e4hat_i = Fraction(4, 1) * e4_2i_over_e4_i
    e4hat_i_over_2_over_e4hat_i = (
        Fraction(1, 4) * e4_i_over_2_over_e4_i
    )
    e4hat_diagonal_over_e4hat_i = (
        Fraction(1, 4) * e4_diagonal_over_e4_i
    )

    return {
        "weight": Fraction(weight, 1),
        "T2_eigenvalue": hecke_eigenvalue,
        "E4_2i_over_E4_i": e4_2i_over_e4_i,
        "E4_i_over_2_over_E4_i": e4_i_over_2_over_e4_i,
        "E4_diagonal_over_E4_i": e4_diagonal_over_e4_i,
        "E4hat_2i_over_E4hat_i": e4hat_2i_over_e4hat_i,
        "E4hat_i_over_2_over_E4hat_i": e4hat_i_over_2_over_e4hat_i,
        "E4hat_diagonal_over_E4hat_i": e4hat_diagonal_over_e4hat_i,
    }


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def main() -> int:
    payload = {name: _fraction_text(value) for name, value in exact_ratios().items()}
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
