#!/usr/bin/env python3
"""Minimal exact algebra oracle for a two-channel Q=1 crossing confluence."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


Matrix = list[list[Fraction]]


def rank(matrix: Matrix) -> int:
    rows = [row[:] for row in matrix]
    out = 0
    for col in range(len(rows[0]) if rows else 0):
        pivot = next((i for i in range(out, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[out], rows[pivot] = rows[pivot], rows[out]
        value = rows[out][col]
        rows[out] = [entry / value for entry in rows[out]]
        for i in range(len(rows)):
            if i == out or not rows[i][col]:
                continue
            value = rows[i][col]
            rows[i] = [a - value * b for a, b in zip(rows[i], rows[out])]
        out += 1
    return out


def matvec(matrix: Matrix, vector: list[Fraction]) -> list[Fraction]:
    return [sum(a * b for a, b in zip(row, vector)) for row in matrix]


def dot(left: list[Fraction], right: list[Fraction]) -> Fraction:
    return sum(a * b for a, b in zip(left, right))


def toy_jet_columns(t_values: list[int]) -> dict[str, list[Fraction]]:
    return {
        "block": [Fraction(1) for _ in t_values],
        "d_delta_block": [Fraction(t) for t in t_values],
        "half_d2_delta_block": [Fraction(t * t, 2) for t in t_values],
    }


def analyze() -> dict:
    t_values = [-2, -1, 1, 2]
    columns = toy_jet_columns(t_values)
    rank2_matrix = [[columns["block"][i], columns["d_delta_block"][i]] for i in range(4)]
    rank3_matrix = [row + [columns["half_d2_delta_block"][i]] for i, row in enumerate(rank2_matrix)]
    nulls_rank2 = [
        [Fraction(-2), Fraction(3), Fraction(-1), Fraction(0)],
        [Fraction(0), Fraction(-1), Fraction(3), Fraction(-2)],
    ]
    rank3_null = [a - b for a, b in zip(nulls_rank2[0], nulls_rank2[1])]
    quadratic_responses = [dot(row, columns["half_d2_delta_block"]) for row in nulls_rank2]
    logistic_cross_ratios = [math.exp(t) / (1 + math.exp(t)) for t in t_values]
    return {
        "schema": "matching-one/linearized-potts-crossing-toy/v1",
        "issue": 263,
        "status": "exact_two_channel_bookkeeping_and_frozen_jet_rank_score",
        "two_channel_confluence": {
            "laurent_input": [
                "A_+(e)=R/e+a_+ + b_+ e+O(e^2)",
                "A_-(e)=-R/e+a_- + b_- e+O(e^2)",
                "Delta_+(e)=Delta_*+v_+ e+w_+ e^2/2+O(e^3)",
                "Delta_-(e)=Delta_*+v_- e+w_- e^2/2+O(e^3)",
            ],
            "finite_Q1_coefficients": {
                "block": "a_+ + a_-",
                "d_delta_block": "R(v_+-v_-)",
            },
            "Q_tangent_coefficients": {
                "block": "b_+ + b_-",
                "d_delta_block": "a_+v_+ + a_-v_- + R(w_+-w_-)/2",
                "d2_delta_block": "R(v_+^2-v_-^2)/2",
            },
            "regularity_gate": "If Delta_+(1)!=Delta_-(1), the uncancelled R[B(Delta_+)-B(Delta_-)]/e divergence forbids this simple opposite-residue confluence.",
        },
        "exact_positive_control": {
            "channel_plus": "energy singlet, residue +R",
            "channel_minus": "two-cluster [2], residue -R",
            "dimensions_at_Q1": "x_energy=x_2cluster=5/4",
            "v_energy": "-9*sqrt(3)/(16*pi)",
            "v_2cluster": "7*sqrt(3)/(16*pi)",
            "v_energy_minus_v_2cluster": "-sqrt(3)/pi",
            "finite_d_delta_block_coefficient": "-R*sqrt(3)/pi",
            "Q_tangent_d2_delta_block_coefficient": "3R/(16*pi^2)",
            "derivation_boundary": "uses the VJS gap derivative d_Q(x_2cluster-x_energy)=sqrt(3)/pi and the #261 energy-family velocity",
        },
        "spin4_adversary": {
            "thermal_Q4": {"x": "21/4", "dx_dQ": "-9*sqrt(3)/(16*pi)"},
            "V_2_2_[2]": {"x": "17/4", "dx_dQ": "-5*sqrt(3)/(16*pi)"},
            "velocity_gap": "sqrt(3)/(4*pi)",
            "dimension_gap_at_Q1": "-1",
            "conclusion": "These fields are not a confluent pair. Their opposite projector residues cannot cancel between their distinct blocks; they must be kept as separate crossing columns.",
        },
        "coefficient_sources": {
            "projector_issue_262": "fixes opposite residue R=2J and its sign",
            "velocity_issue_261": "fixes v for declared generic-Q families",
            "normalized_OPE_issue_250": "kappa=d_Q log(C12p*C34p) shifts finite a and tangent d_delta coefficients but cannot create the d2_delta coefficient",
            "measure_score_issue_258": "supplies the lattice measure component Cov(G,T); it must be added to projector and explicit-field terms before the crossing vector is scored",
        },
        "finite_cross_ratio_toy": {
            "meaning": "a stripped block jet with t=(d_Delta B)/B; this is a rank oracle, not a numerical conformal block approximation",
            "t_values": t_values,
            "logistic_cross_ratios_for_display": logistic_cross_ratios,
            "rank2_columns": ["1", "t"],
            "rank3_columns": ["1", "t", "t^2/2"],
            "matrix_ranks": {"rank2": rank(rank2_matrix), "rank3": rank(rank3_matrix)},
            "rank2_left_null_covectors": [[str(value) for value in row] for row in nulls_rank2],
            "rank3_left_null_covector": [str(value) for value in rank3_null],
            "pure_half_d2_responses_under_rank2_nulls": [str(value) for value in quadratic_responses],
        },
        "frozen_score": {
            "input": "four reconstructed, block-normalized tangent values y_i at ordered t=(-2,-1,1,2)",
            "residuals": {
                "r_left": "-2 y_1 + 3 y_2 - y_3",
                "r_right": "-y_2 + 3 y_3 - 2 y_4",
            },
            "rank2_prediction": "r_left=r_right=0",
            "rank3_second_derivative_prediction": "r_left=r_right=-3 C for y=A+Bt+C t^2/2; hence r_left-r_right=0",
            "higher_rank_falsifier": "r_left-r_right != 0 after covariance propagation and declared ordinary-column subtraction",
            "general_block_form": "At arbitrary cross ratios, rank2 requires det[B,d_Delta B,Y]=0 on every triple; rank3 requires det[B,d_Delta B,d_Delta^2 B,Y]=0 on every quadruple.",
        },
        "claim_boundary": {
            "exact": "Laurent cancellation algebra, VJS jet coefficients conditional on its published velocity gap, and finite matrix/null relations",
            "not_claimed": "that Q4 and V_(2,2) collide, that the toy t coordinate is an exact Potts block, or that the measure score alone is the tangent field",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.dumps(analyze(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
