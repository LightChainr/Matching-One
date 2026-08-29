#!/usr/bin/env python3
"""Exact triangular-site E6 control using the nondegenerate center score."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import mpmath as mp

from derive_hexagonal_degree2_hecke import normalized_eisenstein


@dataclass(frozen=True)
class Eisenstein:
    """a+b*omega with omega=exp(i*pi/3), omega^2=omega-1."""

    a: int
    b: int

    def __mul__(self, other: "Eisenstein") -> "Eisenstein":
        return Eisenstein(
            self.a * other.a - self.b * other.b,
            self.a * other.b + self.b * other.a + self.b * other.b,
        )

    def __pow__(self, exponent: int) -> "Eisenstein":
        if exponent < 0:
            raise ValueError("use inverse_power_record for negative powers")
        result = Eisenstein(1, 0)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent //= 2
        return result

    def norm(self) -> int:
        return self.a * self.a + self.a * self.b + self.b * self.b

    def conjugate(self) -> "Eisenstein":
        return Eisenstein(self.a + self.b, -self.b)

    def multiplication_matrix(self) -> list[list[int]]:
        return [[self.a, -self.b], [self.b, self.a + self.b]]


def inverse_sixth_record(value: Eisenstein) -> dict:
    numerator = value.conjugate() ** 6
    denominator = value.norm() ** 6
    phase_denominator = value.norm() ** 3
    return {
        "numerator_basis_1_omega": [numerator.a, numerator.b],
        "denominator": denominator,
        "text": f"({numerator.a}+{numerator.b}*omega)/{denominator}",
        "N_cubed_times_inverse_sixth": f"({numerator.a}+{numerator.b}*omega)/{phase_denominator}",
    }


def normalized_sixth_phase(value: Eisenstein) -> dict:
    numerator = value.conjugate() ** 6
    denominator = value.norm() ** 3
    # a+b omega has real part a+b/2 and imaginary part b*sqrt(3)/2.
    real = Fraction(2 * numerator.a + numerator.b, 2 * denominator)
    imag_sqrt3 = Fraction(numerator.b, 2 * denominator)
    return {
        "basis_1_omega": [str(Fraction(numerator.a, denominator)), str(Fraction(numerator.b, denominator))],
        "real": str(real),
        "imaginary_coefficient_of_sqrt3": str(imag_sqrt3),
    }


def q6_coefficients() -> dict:
    q6_over_primary_g3 = Fraction(-3975, 224)
    # g3=(8*pi^6/27)E6 in the repository period convention.
    q6_over_primary_pi6_e6 = q6_over_primary_g3 * Fraction(8, 27)
    return {
        "representative": "-25 L_-6 + 28 L_-5 L_-1 - 56 L_-4 L_-2 + 35 L_-3^2",
        "quotient_dimension": 2,
        "Ward_map_rank": 1,
        "Q6_over_primary_g3": str(q6_over_primary_g3),
        "Q6_over_primary_pi6_E6": str(q6_over_primary_pi6_e6),
    }


def exact_center_derivatives(root: Path) -> list[dict]:
    source = json.loads((root / "results/server-20260828/C04/exact_regression/triangular.json").read_text())
    return [
        {
            "L": row["L"],
            "M_at_half": row["exact"]["mean_observable"],
            "Mprime_at_half": row["exact"]["first_derivative"],
        }
        for row in source["results"]
    ]


def numerical_e6_oracle(dps: int = 90) -> dict:
    with mp.workdps(dps):
        rho = (1 + mp.sqrt(3) * mp.j) / 2
        square = mp.j
        children = [2 * rho, rho / 2, (rho + 1) / 2]
        parent = normalized_eisenstein(6, rho, dps=dps)[0]
        child_values = [normalized_eisenstein(6, tau, dps=dps)[0] for tau in children]
        square_value = normalized_eisenstein(6, square, dps=dps)[0]
        return {
            "dps": dps,
            "parent": mp.nstr(parent, 70),
            "children": [mp.nstr(value, 70) for value in child_values],
            "square": mp.nstr(square_value, 70),
            "errors": {
                "equal_children": mp.nstr(max(abs(value / child_values[0] - 1) for value in child_values), 12),
                "child_parent_11_over_4": mp.nstr(max(abs(value / parent - mp.mpf(11) / 4) for value in child_values), 12),
                "square_zero": mp.nstr(abs(square_value), 12),
            },
        }


def analyze(root: Path, dps: int = 90) -> dict:
    multiplier = Eisenstein(1, 1)
    g1 = Eisenstein(1, 9)
    g2 = Eisenstein(5, 6)
    p1 = normalized_sixth_phase(g1)
    p2 = normalized_sixth_phase(g2)
    delta_cos = Fraction(p1["real"]) - Fraction(p2["real"])
    return {
        "schema": "matching-one/triangular-e6-center-score/v1",
        "issue": 165,
        "status": "exact_non_degenerate_observable_and_frozen_E6_controls",
        "self_matching_correction": {
            "identity": "M_L(1-p)=-M_L(p)",
            "degenerate": "M_L(1/2)=0 and the matching root is pinned",
            "nondegenerate": "K_L=M'_L(1/2) is complement-even and need not vanish",
            "exact_repository_values": exact_center_derivatives(root),
            "replacement_observable": "the same-N sextic orientation projector of K_L, divided by the mean center slope",
        },
        "level6_Ward_anchor": q6_coefficients(),
        "Eisenstein_convention": {
            "omega": "exp(i*pi/3), omega^2=omega-1",
            "norm": "N(a+b omega)=a^2+ab+b^2",
            "multiplication_matrix": "[[a,-b],[b,a+b]]",
            "determinant_equals_norm": True,
        },
        "same_N91_projector": {
            "orientations": [
                {"g": [g1.a, g1.b], "norm": g1.norm(), "normalized_g^-6": p1},
                {"g": [g2.a, g2.b], "norm": g2.norm(), "normalized_g^-6": p2},
            ],
            "delta_cos6": str(delta_cos),
            "statistic": "D6=[K(1+9omega)-K(5+6omega)]/delta_cos6",
            "normalized_statistic": "R6=D6/Kbar, Kbar=(K1+K2)/2",
            "scaled_amplitude": "A6=N^3 R6",
            "reason": "same N cancels scalar size dependence; real C6 coupling removes the sine quadrature",
        },
        "multiplier_control": {
            "m": [1, 1],
            "norm": multiplier.norm(),
            "matrix": multiplier.multiplication_matrix(),
            "m_inverse_6": inverse_sixth_record(multiplier),
            "raw_relative_response": "R6(mg)/R6(g)=-1/27",
            "scaled_amplitude": "A6(mg)/A6(g)=-1",
            "alias_discriminator": "spin 12 has scaled factor +1 under the same multiplier",
            "source_target": "N=91 orientation pair -> N=273 child pair",
        },
        "modulus_control": {
            "natural_hex": "E6hat(rho) is nonzero",
            "square_zero": "E6hat(i)=0",
            "degree2_children": "E6hat(child_j)/E6hat(rho)=11/4 for all three children",
            "frozen_nulls": [
                "W6(child_1)-W6(child_0)=0",
                "W6(child_2)-W6(child_0)=0",
                "4 W6(child_0)-11 W6(parent)=0",
                "W6(i)=0",
            ],
        },
        "minimal_execution": {
            "phase_A": "use exact p=1/2 score K=E[M*4(k-N/2)] on the two N=91 orientations",
            "phase_B": "repeat with common random fields on the m=1+omega N=273 children",
            "primary_scores": ["A6_child/A6_parent+1", "R6_child/R6_parent+1/27"],
            "controls": ["M(1/2)=0 exactly", "Kbar nonzero", "orientation swap/conjugation", "spin12 scaled alternative +1"],
            "compute_class": "small paired Monte Carlo; no threshold search and no off-critical scan",
        },
        "interpretation_boundary": {
            "exact": "Q6 Ward image, Eisenstein phases, multiplier factors, E6 Hecke ratios, and existing tiny center derivatives",
            "conditional": "K-slope sextic projector must overlap the Ward-active level-6 line rather than the Ward-null quotient line",
            "not_claimed": "that raw K itself is E6, that M is identically zero away from p=1/2, or that one real orientation pair resolves arbitrary chiral mixing",
        },
        "numerical_E6_oracle": numerical_e6_oracle(dps),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    payload = json.dumps(analyze(root, args.dps), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
