#!/usr/bin/env python3
"""Exact modular fixed-point selection rules for the c=0 vacuum KdV I3 response."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import mpmath as mp

from pinson_arguin_kdv import SECTORS, primitive_k4_holomorphic_series


@dataclass(frozen=True)
class C3:
    """Element a+b*omega of Q(omega), omega^2+omega+1=0."""

    a: Fraction
    b: Fraction

    def __add__(self, other: "C3") -> "C3":
        return C3(self.a + other.a, self.b + other.b)

    def __sub__(self, other: "C3") -> "C3":
        return C3(self.a - other.a, self.b - other.b)

    def __mul__(self, other: "C3") -> "C3":
        # (a+bw)(c+dw)=(ac-bd)+(ad+bc-bd)w, since w^2=-1-w.
        return C3(
            self.a * other.a - self.b * other.b,
            self.a * other.b + self.b * other.a - self.b * other.b,
        )

    def scale(self, scalar: Fraction) -> "C3":
        return C3(self.a * scalar, self.b * scalar)

    def conjugate(self) -> "C3":
        # conjugate(w)=w^2=-1-w.
        return C3(self.a - self.b, -self.b)

    def text(self) -> str:
        if self.b == 0:
            return str(self.a)
        return f"({self.a})+({self.b})*omega"


ZERO = C3(Fraction(0), Fraction(0))
ONE = C3(Fraction(1), Fraction(0))
OMEGA = C3(Fraction(0), Fraction(1))
OMEGA2 = OMEGA * OMEGA


def cycle(vector: list[C3]) -> list[C3]:
    """Registered forward C3 transport: (v0,v1,v2)->(v1,v2,v0)."""
    return [vector[1], vector[2], vector[0]]


def vector_scale(scalar: C3, vector: list[C3]) -> list[C3]:
    return [scalar * value for value in vector]


def vector_add(*vectors: list[C3]) -> list[C3]:
    return [sum((vector[i] for vector in vectors), ZERO) for i in range(len(vectors[0]))]


def projector_numerator(character: C3, vector: list[C3]) -> list[C3]:
    """Three times the projector onto the stated cycle eigenvalue."""
    if character == ONE:
        return vector_add(vector, cycle(vector), cycle(cycle(vector)))
    if character == OMEGA:
        return vector_add(vector, vector_scale(OMEGA2, cycle(vector)), vector_scale(OMEGA, cycle(cycle(vector))))
    if character == OMEGA2:
        return vector_add(vector, vector_scale(OMEGA, cycle(vector)), vector_scale(OMEGA2, cycle(cycle(vector))))
    raise ValueError("character must be 1, omega, or omega^2")


def reflection_even(vector: list[C3]) -> list[C3]:
    return [value + value.conjugate() for value in vector]


def exact_selection() -> dict:
    chiral = [ONE, OMEGA, OMEGA2]
    even = reflection_even(chiral)
    return {
        "chiral_vector": [value.text() for value in chiral],
        "cycle_eigenvalue": "omega",
        "cycle_check": [value.text() for value in cycle(chiral)],
        "omega_times_vector": [value.text() for value in vector_scale(OMEGA, chiral)],
        "projector_numerators": {
            "trivial": [value.text() for value in projector_numerator(ONE, chiral)],
            "omega": [value.text() for value in projector_numerator(OMEGA, chiral)],
            "omega2": [value.text() for value in projector_numerator(OMEGA2, chiral)],
        },
        "reflection_even_vector": [value.text() for value in even],
        "reflection_even_normalized": ["1", "-1/2", "-1/2"],
        "C_Q_S": ["3/2", "0", "0"],
    }


def numerical_rho_oracle(dps: int = 90) -> dict:
    with mp.workdps(dps + 20):
        rho = (1 + mp.sqrt(3) * mp.j) / 2
        values = [primitive_k4_holomorphic_series(a, b, rho, dps=dps) for a, b in SECTORS]
        omega = mp.exp(2 * mp.pi * mp.j / 3)
        scale = values[0]
        even = [2 * mp.re(value) for value in values]
        c = even[0] - (even[1] + even[2]) / 2
        q = mp.sqrt(3) * (even[2] - even[1]) / 2
        s = sum(even)
        return {
            "dps": dps,
            "sector_order": [list(sector) for sector in SECTORS],
            "holomorphic_values": [mp.nstr(value, 70) for value in values],
            "amplitude": mp.nstr(scale, 70),
            "ratio_to_first": [mp.nstr(value / scale, 70) for value in values],
            "errors": {
                "sector_01_minus_omega": mp.nstr(abs(values[1] / scale - omega), 12),
                "sector_11_minus_omega2": mp.nstr(abs(values[2] / scale - omega**2), 12),
                "reflection_Q_abs": mp.nstr(abs(q), 12),
                "reflection_S_abs": mp.nstr(abs(s), 12),
                "C_minus_3A_abs": mp.nstr(abs(c - 3 * scale), 12),
            },
        }


def analyze(dps: int = 90) -> dict:
    return {
        "schema": "matching-one/vacuum-kdv-fixedpoint-selection/v1",
        "issue": 231,
        "status": "exact_modular_selection_and_frozen_sector_character",
        "vacuum_degeneracy": {
            "Q1_total_partition": "Z_total=1",
            "I3_mean": "D2 D0[1]=0",
            "I3_second_moment": "[D^4+(1/18)E4 D^2-(11/1080)E6 D][1]=0",
            "selection_rule": "At c=0 the normalized unconditioned vacuum KdV response vanishes; a nonzero I3 signal must be sector-resolved, defect-twisted, or a lattice contact/nonuniversal term.",
        },
        "fixed_point_covariance": {
            "rho": "exp(i*pi/3)",
            "stabilizer": "gamma(tau)=(tau-1)/tau with (c tau+d)=tau",
            "weight_w_rule": "R_gamma v(rho)=rho^(-w) v(rho)",
            "weight4_allowed_character": "omega=exp(2*pi*i/3)",
            "weight8_allowed_character": "omega^2",
            "square_rule": "at tau=i, weight 4 lies in the +1 eigenspace of the sector swap",
        },
        "exact_rho_sector_oracle": exact_selection(),
        "frozen_predictions": {
            "chiral_ratios": "K4_01/K4_10=omega and K4_11/K4_10=omega^2",
            "reflection_even_ratios": "R_10:R_01:R_11=2:-1:-1",
            "contrast_nulls": "Q_rho=0 and S_rho=0 exactly; C_rho is allowed and nonzero",
            "character_purity": "P_omega K4=K4, P_1 K4=P_omega2 K4=0",
            "thermal_Q4_discriminator": "a scalar/trivial-sector weight-4 one-point is forced to zero at rho (the E4 node), while the sector-valued vacuum KdV omega character is allowed",
            "four_leg_discriminator": "an unmarked singlet readout cannot linearly see the [2] four-leg primary; the nonzero rho omega-sector response requires topological-sector resolution rather than a hidden scalar V_(2,2) amplitude",
        },
        "minimal_score": {
            "measure": "the complex primitive-sector response in registered order [(1,0),(0,1),(1,1)] at rho",
            "zero_parameter_residuals": [
                "K4_01-omega*K4_10",
                "K4_11-omega^2*K4_10",
                "R_10+R_01+R_11",
                "R_01-R_11",
            ],
            "weight8_adversary": "the chiral character is conjugated (omega^2), even though the reflection-even 2:-1:-1 projection alone cannot distinguish it",
        },
        "numerical_oracle": numerical_rho_oracle(dps),
        "claim_boundary": {
            "exact": "c=0 total-partition degeneracy, vector-valued modular fixed-point selection, and Pinson-sector character ratios",
            "conditional_bridge": "matching-even H4 must couple to the same sector-resolved I3 insertion and registered modular frame",
            "not_claimed": "that an unconditioned matching observable is I3, or that the real projection alone distinguishes weight4 from its weight8 conjugate",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.dumps(analyze(args.dps), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
