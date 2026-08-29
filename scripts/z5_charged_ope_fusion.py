#!/usr/bin/env python3
"""Exact Z5 charged-fusion oracle and norm-five cubic phase predictions."""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations_with_replacement
import json
from pathlib import Path

from norm5_chiral_hecke_phase import gaussian_ratio_power


CHARGES = tuple(range(5))
NONTRIVIAL = tuple(range(1, 5))
SPIN_HYPOTHESES = (4, 8, 12)


def conjugate_charge(charge: int) -> int:
    return (-charge) % 5


def conjugate_triple(triple: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(sorted(conjugate_charge(charge) for charge in triple))


def neutral(triple: tuple[int, ...]) -> bool:
    return sum(triple) % 5 == 0


def primitive_neutral_triples() -> list[tuple[int, int, int]]:
    """Neutral unordered triples with no neutral one- or two-point subchannel."""

    rows = []
    for triple in combinations_with_replacement(NONTRIVIAL, 3):
        if not neutral(triple):
            continue
        if any((triple[i] + triple[j]) % 5 == 0 for i in range(3) for j in range(i + 1, 3)):
            continue
        rows.append(triple)
    return rows


def determinant(matrix: list[list[int]]) -> int:
    values = [[Fraction(value) for value in row] for row in matrix]
    result = Fraction(1)
    for column in range(len(values)):
        pivot = next(row for row in range(column, len(values)) if values[row][column])
        if pivot != column:
            values[column], values[pivot] = values[pivot], values[column]
            result = -result
        scale = values[column][column]
        result *= scale
        values[column] = [value / scale for value in values[column]]
        for row in range(column + 1, len(values)):
            factor = values[row][column]
            values[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(values[row], values[column])
            ]
    if result.denominator != 1:
        raise AssertionError("integer matrix has noninteger determinant")
    return result.numerator


def phase_record(spin: int) -> dict[str, object]:
    real, imag, denominator = gaussian_ratio_power(3 * spin)
    return {
        "external_spin_each": spin,
        "cubic_total_spin": 3 * spin,
        "exact_handed_ratio": {
            "real": f"{real}/{denominator}",
            "imag": f"{imag}/{denominator}",
        },
        "unit_modulus_check": real * real + imag * imag == denominator * denominator,
    }


def build_artifact() -> dict[str, object]:
    triples = primitive_neutral_triples()
    gauge_matrix = [
        # Columns are C_113, C_122, C_442, C_433.  Rows are field
        # rescalings lambda_1,...,lambda_4.
        [2, 1, 0, 0],
        [0, 2, 1, 0],
        [1, 0, 0, 2],
        [0, 0, 2, 1],
    ]
    return {
        "schema": "matching-one/z5-charged-ope-fusion/v1",
        "issue": 250,
        "status": "exact_fusion_and_frozen_conditional_prediction",
        "dependencies": {
            "chiral_gls_branch": "analysis/norm5-chiral-hecke-phase-20260829",
            "selection_rule_issue": 244,
        },
        "exact_character_algebra": {
            "group": "Z/5Z",
            "fusion_rule": "chi_r tensor chi_s = chi_(r+s mod 5)",
            "fusion_target_table": [
                [(left + right) % 5 for right in CHARGES] for left in CHARGES
            ],
            "two_point_support": [[charge, conjugate_charge(charge)] for charge in NONTRIVIAL],
            "primitive_connected_neutral_triples": [list(row) for row in triples],
            "conjugation_pairs": [
                [list(row), list(conjugate_triple(row))]
                for row in triples
                if row <= conjugate_triple(row)
            ],
            "fixed_generator_cubic_channels": {
                "A": [1, 1, 3],
                "B": [1, 2, 2],
                "A_conjugate": [2, 4, 4],
                "B_conjugate": [3, 3, 4],
            },
            "connectedness_reason": (
                "Every proper nonempty subset of A or B has nonzero Z5 charge, "
                "so all one-point and disconnected two-point factors vanish exactly."
            ),
        },
        "normalization_and_phase": {
            "field_rescaling_exponent_matrix_rows_lambda1_to_lambda4_columns_A_B_Abar_Bbar": gauge_matrix,
            "determinant": determinant(gauge_matrix),
            "phase_no_go": (
                "det=-15 != 0: no nonconstant monomial of only the four primitive cubic "
                "coefficients is invariant under independent charged-field rescalings."
            ),
            "normalization_free_squared_coupling": (
                "I_rst=C_rst*C_(-r,-s,-t)*C_000/"
                "(C_(0,r,-r)*C_(0,s,-s)*C_(0,t,-t))"
            ),
            "interpretation": (
                "I_rst is universal without choosing charged-field phases; a complex OPE "
                "phase requires an exact transported deck basis or a larger closed fusion loop."
            ),
        },
        "frozen_next_score": {
            "primary_channels": ["C_113_plus", "C_113_minus", "C_122_plus", "C_122_minus"],
            "zero_parameter_complex_null": (
                "C_113_plus*C_122_minus-C_113_minus*C_122_plus=0"
            ),
            "null_degrees_of_freedom": 2,
            "joint_gls_model": (
                "C_A_plus=q_(3s)*a, C_A_minus=a, "
                "C_B_plus=q_(3s)*b, C_B_minus=b"
            ),
            "joint_gls_nuisance_parameters": "two complex amplitudes a,b (4 real)",
            "joint_gls_degrees_of_freedom": 4,
            "candidate_cubic_phases": {
                f"H{spin}": phase_record(spin) for spin in SPIN_HYPOTHESES
            },
            "score_order": [
                "score the two-real-component cross-product null first",
                "then rank the frozen H4/H8/H12 q_(3s) joint GLS models",
                "report I_113 and I_122 with joint covariance as phase-gauge-free OPE magnitudes",
            ],
            "mechanism_assumptions": [
                "the same deck character basis is transported exactly between the two hands",
                "the same local angular spin-s insertion is used in all three legs",
                "one Hecke eigenfield dominates each of the A and B cubic channels",
                "the marked three-point geometry is matched between the two hands",
            ],
        },
        "minimal_experiment": {
            "geometry": "the existing N325 same-parent children (8+i)(2+i) and (8+i)(2-i)",
            "insertions": (
                "one randomly translated fixed-shape three-anchor pattern per replica; "
                "evaluate the existing charged landing-marked local row at all three anchors"
            ),
            "charge_patterns": [[1, 1, 3], [1, 2, 2]],
            "required_covariance": "full 8x8 covariance of real and imaginary parts of four complex channels",
            "controls": [
                "configurationwise global deck-translation charge conservation",
                "reflected conjugates C_244=conj(C_113) and C_334=conj(C_122)",
                "forbidden nonneutral triples remain exact zero or labeling-error diagnostics",
            ],
        },
        "claim_boundary": {
            "exact": "Z5 fusion/support, connectedness, determinant phase no-go, and invariant I_rst formula",
            "conditional_prediction": "the q_(3s) handed law and cross-channel closure null",
            "not_claimed": "that the lattice charged row is already a particular loop-CFT primary",
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
