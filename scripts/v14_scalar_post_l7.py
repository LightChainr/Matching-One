#!/usr/bin/env python3
"""Exact arithmetic for the V_<1,4> scalar post-L^-7 mechanism.

Critical Potts branch at Q=1:

    h_{r,s}=((2r-3s)^2-1)/24.

The diagonal singlet V_<1,4> has x=33/4 and spin 0.  If the stronger
interchiral/OPE matching automorphism gives eta_s=(-1)^(s-1), it is
matching-odd.  Its dimensionless torus contribution is L^-25/4, exactly three
length powers below the leading candidate thermal spin-4 L^-13/4 term.  After
annihilating the leading term, the induced root correction is L^-7.

The parity/coupling statement is conditional; the spectrum and exponent
arithmetic are exact within the stated critical-Potts convention.
"""
from __future__ import annotations

from fractions import Fraction as F


def critical_potts_h(r: int, s: int) -> F:
    return F((2 * r - 3 * s) ** 2 - 1, 24)


def diagonal_x(s: int) -> F:
    return 2 * critical_potts_h(1, s)


def conditional_interchiral_parity(s: int) -> int:
    """Conditional OPE-level parity eta_s=(-1)^(s-1)."""
    return 1 if s % 2 else -1


def torus_residual_length_power(x: F) -> F:
    """Return positive rho with dimensionless correction L^-rho."""
    return x - 2


def annihilated_root_power(relative_q: F, leading_root_power: F = F(4)) -> F:
    return leading_root_power + relative_q


def mechanism_table() -> dict[str, dict[str, object]]:
    leading_x = F(21, 4)
    v14_x = diagonal_x(4)
    thermal_next_x = F(45, 4)

    return {
        "leading_thermal_H4": {
            "x": leading_x,
            "spin": 4,
            "harmonic": "H4",
            "residual_power_L": torus_residual_length_power(leading_x),
            "residual_power_N": torus_residual_length_power(leading_x) / 2,
            "relative_q_L": F(0),
            "root_power_L": F(4),
        },
        "V14_scalar_H0": {
            "x": v14_x,
            "spin": 0,
            "harmonic": "H0",
            "residual_power_L": torus_residual_length_power(v14_x),
            "residual_power_N": torus_residual_length_power(v14_x) / 2,
            "relative_q_L": v14_x - leading_x,
            "root_power_L": annihilated_root_power(v14_x - leading_x),
            "conditional_matching_parity": conditional_interchiral_parity(4),
        },
        "thermal_next_H4": {
            "x": thermal_next_x,
            "spin": 4,
            "harmonic": "H4",
            "residual_power_L": torus_residual_length_power(thermal_next_x),
            "residual_power_N": torus_residual_length_power(thermal_next_x) / 2,
            "relative_q_L": thermal_next_x - leading_x,
            "root_power_L": annihilated_root_power(thermal_next_x - leading_x),
        },
    }


def main() -> int:
    table = mechanism_table()
    v14 = table["V14_scalar_H0"]
    nxt = table["thermal_next_H4"]

    assert diagonal_x(4) == F(33, 4)
    assert critical_potts_h(1, 4) == F(33, 8)
    assert v14["conditional_matching_parity"] == -1
    assert v14["residual_power_L"] == F(25, 4)
    assert v14["residual_power_N"] == F(25, 8)
    assert v14["relative_q_L"] == 3
    assert v14["root_power_L"] == 7
    assert nxt["residual_power_L"] == F(37, 4)
    assert nxt["relative_q_L"] == 6
    assert nxt["root_power_L"] == 10

    print("V_<1,4>: h=hbar=33/8, x=33/4, spin=0, H0")
    print("conditional OPE-level matching parity: odd")
    print("M(pc) contribution: L^-25/4 = N^-25/8")
    print("relative to leading L^-13/4: q=3")
    print("leading-annihilated root: L^-7")
    print("next ordinary thermal H4: L^-37/4, q=6, root L^-10")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
