#!/usr/bin/env python3
"""Exact null-residue and Q4 Ward fingerprints at c=0, h=5/8."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path


Q = Fraction
H = Q(5, 8)
C = Q(0)


@dataclass(frozen=True)
class Dual:
    value: Q
    tangent: Q = Q(0)

    def __add__(self, other):
        other = as_dual(other)
        return Dual(self.value + other.value, self.tangent + other.tangent)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.value, -self.tangent)

    def __sub__(self, other):
        return self + (-as_dual(other))

    def __rsub__(self, other):
        return as_dual(other) - self

    def __mul__(self, other):
        other = as_dual(other)
        return Dual(
            self.value * other.value,
            self.tangent * other.value + self.value * other.tangent,
        )

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = as_dual(other)
        return Dual(
            self.value / other.value,
            (self.tangent * other.value - self.value * other.tangent)
            / (other.value * other.value),
        )


def as_dual(value) -> Dual:
    return value if isinstance(value, Dual) else Dual(Q(value))


def fraction_text(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def null_coefficient(h):
    """Coefficient a(h) in N2(h)=L_-2-a(h)L_-1^2."""

    return Q(3, 2) / (2 * h + 1)


def level2_kac_central_charge(h):
    """Level-two Kac curve on which N2(h) is singular."""

    return 2 * h * (5 - 8 * h) / (2 * h + 1)


def ordinary_null_jet(h: Q) -> tuple[Q, Q, Q]:
    """Coefficients of L1 R, L1^2 R, L2 R for R=N2(h0)|h>."""

    first = (5 - 8 * h) / 3
    return first, 2 * h * first, Q(0)


def fixed_c_jordan_null_jet() -> tuple[Q, Q, Q]:
    # L0|psi>=h|psi>+|phi>, Ln|psi>=0 for n>0.
    return Q(-8, 3), Q(-10, 3), Q(0)


def kac_tangent_null_jet() -> tuple[Q, Q, Q]:
    # Differentiate N2(h)|phi_h>=0 along the level-two Kac curve c(h).
    return Q(-8, 3), Q(-10, 3), Q(-20, 9)


def q4_scalar_ward_vector(h, c):
    """Return (L2^2, L1 L3, L4) Q4|h> as scalar polynomials.

    Q4=40 L_-2^2-60 L_-3 L_-1-9 L_-4 is kept in the repository
    normalization.  Inputs may be Fractions or first-order Dual numbers.
    """

    h = as_dual(h)
    c = as_dual(c)
    l2_coefficient = 320 * h + 266 + 40 * c
    w22 = l2_coefficient * (4 * h + c / 2) - 1800 * h
    w13 = 2 * h * (177 - 360 * h - 120 * c)
    w4 = 48 * h - 45 * c
    return w22, w13, w4


def covectors(vector: tuple[Q, Q, Q]) -> tuple[Q, Q]:
    """Two covectors annihilating the ordinary Q4 Ward line (40,-60,30)."""

    return 3 * vector[0] + 2 * vector[1], vector[1] + 2 * vector[2]


def singular_levels(limit_r: int = 20, limit_s: int = 30) -> list[int]:
    return sorted(
        {
            r * s
            for r in range(1, limit_r + 1)
            for s in range(1, limit_s + 1)
            if abs(3 * r - 2 * s) == 4
        }
    )


def jet_record(values: tuple[Q, Q, Q]) -> dict[str, object]:
    first, second, l2 = values
    return {
        "a_L1R_over_Lminus1phi": fraction_text(first),
        "b_L1squaredR_over_phi": fraction_text(second),
        "d_L2R_over_phi": fraction_text(l2),
        "b_over_a_when_defined": None if first == 0 else fraction_text(second / first),
        "d_over_b_when_defined": None if second == 0 else fraction_text(l2 / second),
    }


def build_artifact() -> dict[str, object]:
    h_dual = Dual(H, Q(1))
    c_fixed = Dual(C, Q(0))
    c_kac = level2_kac_central_charge(h_dual)
    if c_kac.value != 0 or c_kac.tangent != Q(-40, 9):
        raise AssertionError("unexpected level-two Kac tangent")

    base_dual = q4_scalar_ward_vector(Dual(H), Dual(C))
    base = tuple(value.value for value in base_dual)
    fixed_dual = q4_scalar_ward_vector(h_dual, c_fixed)
    fixed_tangent = tuple(value.tangent for value in fixed_dual)
    kac_dual = q4_scalar_ward_vector(h_dual, c_kac)
    kac_tangent = tuple(value.tangent for value in kac_dual)
    fixed_residue = covectors(fixed_tangent)
    kac_residue = covectors(kac_tangent)

    negative_control_h = Q(-1, 24)
    four_leg_spin_plus_h = Q(33, 8)
    four_leg_spin_minus_h = Q(1, 8)

    return {
        "schema": "matching-one/null-field-residue-ward/v1",
        "issue": 252,
        "status": "exact Virasoro algebra and frozen conditional lattice prediction",
        "thermal_level2": {
            "c": "0",
            "h": "5/8",
            "null_operator": "N2=L_-2-(2/3)L_-1^2",
            "coefficient_a_h": "3/[2(2h+1)]",
            "level2_kac_curve": "c(h)=2h(5-8h)/(2h+1)",
            "kac_curve_derivative_at_5_over_8": fraction_text(c_kac.tangent),
        },
        "inhomogeneous_null_jets": {
            "definition": (
                "R=N2|psi>; L1 R=a L_-1|phi>, L1^2 R=b|phi>, L2 R=d|phi>, "
                "with (L0-h)|psi>=|phi> for logarithmic rows"
            ),
            "ordinary_thermal_primary": jet_record(ordinary_null_jet(H)),
            "fixed_c_jordan_top": jet_record(fixed_c_jordan_null_jet()),
            "kac_curve_parameter_tangent": jet_record(kac_tangent_null_jet()),
            "negative_dimension_control_h_minus_1_over_24": jet_record(
                ordinary_null_jet(negative_control_h)
            ),
            "four_leg_V22_spin_plus_chiral_h_33_over_8": jet_record(
                ordinary_null_jet(four_leg_spin_plus_h)
            ),
            "four_leg_V22_spin_minus_chiral_h_1_over_8": jet_record(
                ordinary_null_jet(four_leg_spin_minus_h)
            ),
            "exact_discriminator": {
                "ordinary_thermal": "(a,b,d)=(0,0,0)",
                "fixed_c_jordan": "b/a=5/4 and d/b=0",
                "kac_curve_log_tangent": "b/a=5/4 and d/b=2/3",
                "negative_h_primary": "a>0, b<0",
                "V22_spin_plus": "b/a=33/4 with (a,b)=(-28/3,-77)",
            },
        },
        "direct_spin4_ward_gate": {
            "q4_state": "Q4|epsilon>=(40L_-2^2-60L_-3L_-1-9L_-4)|epsilon>",
            "scalar_mode_order": ["L2^2", "L1 L3", "L4"],
            "ordinary_q4_vector": [fraction_text(value) for value in base],
            "primitive_ratio": "4:-6:3",
            "four_leg_primary_vector": ["0", "0", "0"],
            "selection_rule": (
                "A spin-4 primary such as V_(2,-2) is killed by every positive mode, whereas "
                "the thermal Q4 descendant has the nonzero exact vector (40,-60,30)."
            ),
        },
        "q4_jordan_residue": {
            "ordinary_line_annihilating_covectors": [
                "J_A=3 W_(L2^2)+2 W_(L1L3)",
                "J_B=W_(L1L3)+2 W_(L4)",
            ],
            "fixed_c_tangent_vector": [fraction_text(value) for value in fixed_tangent],
            "fixed_c_covector_residue": [fraction_text(value) for value in fixed_residue],
            "fixed_c_ratio_JA_over_JB": fraction_text(fixed_residue[0] / fixed_residue[1]),
            "kac_curve_tangent_vector": [fraction_text(value) for value in kac_tangent],
            "kac_curve_covector_residue": [fraction_text(value) for value in kac_residue],
            "kac_curve_ratio_JA_over_JB": fraction_text(kac_residue[0] / kac_residue[1]),
            "basis_invariance": (
                "psi->psi+lambda phi adds lambda*(40,-60,30); both covectors annihilate that shift"
            ),
        },
        "thermal_tower_boundary": {
            "first_singular_levels": singular_levels()[:3],
            "leading_spin4_pair": [4, 0],
            "next_ordinary_spin4_pair": [7, 3],
            "next_total_level": 10,
            "next_dimension": "45/4",
            "relative_to_leading": "L^-6=N^-3",
            "consequence": (
                "The ordinary level-10 correction stays in the N2 quotient and cannot create a "
                "nonzero inhomogeneous null residue."
            ),
        },
        "frozen_high_risk_prediction": {
            "target": "source-frozen local/torus spin-4 Ward projection of the leading matching-odd row",
            "stage_1_operator_gate": (
                "after common normalization, the three positive-mode contractions lie on 4:-6:3 "
                "rather than the zero vector of the x=17/4 four-leg primary"
            ),
            "stage_2_module_gate": (
                "after subtracting the ordinary Q4 line, (J_A/J_B)=-10/3 and the parent null jet "
                "has d/b=0, selecting a fixed-c Jordan lift"
            ),
            "frozen_adversaries": {
                "ordinary_Q4": "J_A=J_B=0",
                "Kac_curve_parameter_log": "J_A/J_B=-482/185 and parent d/b=2/3",
                "V22_four_leg_primary": "positive-mode vector=0; parent b/a=33/4 for spin +4",
                "negative_dimension_primary": "parent a>0,b<0",
            },
            "claim_status": "conditional until the lattice-to-Virasoro Ward projection is calibrated",
        },
        "minimal_execution": {
            "first_control": (
                "calibrate the same positive-mode/Ward projection on an exact-critical model with a "
                "known thermal degenerate insertion; do not infer modes from the target"
            ),
            "saved_joint_vector": [
                "W_L2squared",
                "W_L1L3",
                "W_L4",
                "parent_L1R",
                "parent_L1squaredR",
                "parent_L2R",
                "bottom normalization rows",
            ],
            "covariance": "one same-stream full covariance; score linear constraints before ratios",
            "cover_test": (
                "transport the source-frozen Ward covectors over one norm-2 and the existing norm-5 "
                "handed pair; the coefficients are representation algebra and must not be refit by cover"
            ),
            "no_new_large_compute_before": "the control reproduces the exact Q4 vector (40,-60,30)",
        },
        "claim_boundary": {
            "exact": "all Virasoro mode coefficients, residue covectors, ratios, and spectrum arithmetic",
            "prediction": "the observed lattice row selects the fixed-c Q4 Jordan fingerprint",
            "not_claimed": "that the existing torus or pivotal rows already implement the required Ward projection",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = build_artifact()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
