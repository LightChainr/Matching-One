#!/usr/bin/env python3
"""Exact arithmetic for the V_<1,4> scalar post-L^-7 mechanism.

The Potts labels used by this repository put the percolation thermal field at
V_<1,2>, so at Q=1

    h_{r,s}=((2r-3s)^2-1)/24.

Then V_<1,4> is a diagonal singlet scalar with x=33/4.  If the matching-pair
RG/OPE action assigns the interchiral recursion eta_s=(-1)^(s-1), V_<1,4>
is matching odd.  Its central torus contribution is L^-25/4, exactly three
length powers below the leading L^-13/4 H4 term; the Mertens-Ziff leading-term
annihilator then generates an L^-7 root correction.
"""
from __future__ import annotations

from fractions import Fraction as F


def critical_potts_h(r: int, s: int) -> F:
    return F((2 * r - 3 * s) ** 2 - 1, 24)


def diagonal_x(s: int) -> F:
    return 2 * critical_potts_h(1, s)


def conditional_interchiral_parity(s: int) -> int:
    """Return eta_s=(-1)^(s-1), conditional on the matching automorphism."""
    return 1 if s % 2 else -1


def root_bias_exponent_in_L(x: F, thermal_y: F = F(3, 4)) -> F:
    """If M(pc)~L^(2-x) and M'~L^y, root bias is L^-[(x-2)+y]."""
    return (x - 2) + thermal_y


def main() -> int:
    thermal_h = critical_potts_h(1, 2)
    leading_x = F(21, 4)
    v14_x = diagonal_x(4)
    central_power = v14_x - 2
    relative_q = v14_x - leading_x
    ordinary_root_power = root_bias_exponent_in_L(v14_x)
    annihilated_root_power = F(4) + relative_q

    assert thermal_h == F(5, 8)
    assert v14_x == F(33, 4)
    assert conditional_interchiral_parity(4) == -1
    assert central_power == F(25, 4)
    assert relative_q == 3
    assert ordinary_root_power == 7
    assert annihilated_root_power == 7

    print("thermal V_<1,2>: h=5/8")
    print("V_<1,4>: h=hbar=33/8, x=33/4, spin=0")
    print("conditional matching parity: odd")
    print("M(pc) scalar contribution: L^-25/4 = N^-25/8")
    print("ordinary scalar-induced root bias: L^-7 = N^-7/2")
    print("relative to leading L^-13/4: q=3")
    print("leading-H4-annihilated root: L^-7")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
