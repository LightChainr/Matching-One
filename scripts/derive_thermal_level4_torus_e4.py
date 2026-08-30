#!/usr/bin/env python3
"""Exact coefficient audit for the c=0,h=5/8 thermal level-4 torus one-point.

Inputs encoded here are the exact algebraic consequences documented in
notes/thermal-level4-torus-onepoint-E4.md:

  <L_-3 L_-1> / <L_-4> = -2
  <L_-2^2> / <L_-4> = 4/3
  Q4 = 40 L_-2^2 - 60 L_-3 L_-1 - 9 L_-4
  <L_-4 phi>/<phi> = h g2 / 20, h=5/8.
"""

from fractions import Fraction


def coefficients() -> dict[str, Fraction]:
    l31 = Fraction(-2, 1)
    l22 = Fraction(4, 3)
    q4_over_l4 = 40 * l22 - 60 * l31 - 9
    h = Fraction(5, 8)
    q4_over_primary_g2 = q4_over_l4 * h / 20
    # g2=(4*pi^4/3)E4 for periods (1,tau): this is the rational
    # coefficient multiplying pi^4 E4.
    q4_over_primary_pi4_e4 = q4_over_primary_g2 * Fraction(4, 3)
    return {
        "L31_over_L4": l31,
        "L22_over_L4": l22,
        "Q4_over_L4": q4_over_l4,
        "Q4_over_primary_g2": q4_over_primary_g2,
        "Q4_over_primary_pi4_E4": q4_over_primary_pi4_e4,
    }


def main() -> int:
    for name, value in coefficients().items():
        print(f"{name}={value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
