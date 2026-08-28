#!/usr/bin/env python3
"""Exact arithmetic for the V_<1,4> scalar post-L^-7 mechanism.

On the critical Potts branch at Q=1,

    h_{r,s}=((2r-3s)^2-1)/24.

The diagonal singlet V_<1,4> has x=33/4.  Under the conditional
interchiral parity rule eta_s=(-1)^(s-1), it is matching odd.  Its linear
contribution to a dimensionless torus observable is L^(2-x)=L^-25/4,
three length powers below the leading thermal spin-4 L^-13/4 term.  A
leading-annihilated root therefore scales as L^-(4+3)=L^-7.
"""
from __future__ import annotations

from fractions import Fraction as F


def critical_potts_h(r: int, s: int) -> F:
    return F((2 * r - 3 * s) ** 2 - 1, 24)


def diagonal_x(s: int) -> F:
    return 2 * critical_potts_h(1, s)


def interchiral_parity(s: int) -> int:
    """Conditional parity from eta_s=(-1)^(s-1)."""
    return 1 if s % 2 else -1


def main() -> int:
    leading_x = F(21, 4)
    v14_x = diagonal_x(4)
    q = v14_x - leading_x
    residual_power = v14_x - 2
    root_power = q + 4

    assert v14_x == F(33, 4)
    assert interchiral_parity(4) == -1
    assert residual_power == F(25, 4)
    assert q == 3
    assert root_power == 7

    print("V_<1,4>: h=hbar=33/8, x=33/4, spin=0")
    print("conditional matching parity: odd")
    print("M(pc) contribution: L^-25/4")
    print("relative to L^-13/4: q=3")
    print("leading-annihilated root: L^-7")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
