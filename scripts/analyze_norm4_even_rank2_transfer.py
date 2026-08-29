#!/usr/bin/env python3
"""Fit the minimal post-reveal norm-4 Jordan-plus-one-even-mode model."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Sequence

import mpmath as mp


ORDERS = (2, 3, 4, 5, 6)
LINEAGES = ((65, 130, 260, 520), (85, 170, 340, 680))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def chi_square_survival(value: mp.mpf, degrees: int) -> float:
    return float(
        mp.gammainc(mp.mpf(degrees) / 2, value / 2, mp.inf)
        / mp.gamma(mp.mpf(degrees) / 2)
    )


def rank_one_gls(
    residual: Sequence[float], covariance: Sequence[Sequence[float]], ratio: mp.mpf
):
    values = mp.matrix([mp.mpf(str(value)) for value in residual])
    precision = mp.matrix(
        [[mp.mpf(str(value)) for value in row] for row in covariance]
    ) ** -1
    design = mp.matrix(10, 5)
    for index in range(5):
        design[index, index] = 1
        design[5 + index, index] = ratio
    direction = mp.lu_solve(
        design.T * precision * design, design.T * precision * values
    )
    error = values - design * direction
    chi_square = (error.T * precision * error)[0]
    return chi_square, direction, error


def fit_ratio(residual, covariance):
    def objective(value: mp.mpf) -> mp.mpf:
        return rank_one_gls(residual, covariance, value)[0]

    grid = [mp.mpf(-8) + mp.mpf(index) / 100 for index in range(1601)]
    best_index = min(range(len(grid)), key=lambda index: objective(grid[index]))
    lower = grid[max(0, best_index - 2)]
    upper = grid[min(len(grid) - 1, best_index + 2)]
    golden = (mp.sqrt(5) - 1) / 2
    left = upper - golden * (upper - lower)
    right = lower + golden * (upper - lower)
    for _ in range(100):
        if objective(left) <= objective(right):
            upper, right = right, left
            left = upper - golden * (upper - lower)
        else:
            lower, left = left, right
            right = lower + golden * (upper - lower)
    ratio = (lower + upper) / 2
    return ratio, rank_one_gls(residual, covariance, ratio)


def profile_interval(residual, covariance, best_ratio, best_chi, delta):
    def equation(value: mp.mpf) -> mp.mpf:
        return rank_one_gls(residual, covariance, value)[0] - best_chi - delta

    def bisect(lower: mp.mpf, upper: mp.mpf) -> mp.mpf:
        f_lower = equation(lower)
        for _ in range(120):
            midpoint = (lower + upper) / 2
            f_midpoint = equation(midpoint)
            if f_lower * f_midpoint <= 0:
                upper = midpoint
            else:
                lower, f_lower = midpoint, f_midpoint
        return (lower + upper) / 2

    def outward(sign: int):
        inner = best_ratio
        step = mp.mpf("0.05")
        for _ in range(80):
            outer = best_ratio + sign * step
            if equation(outer) >= 0:
                return bisect(outer, inner) if sign < 0 else bisect(inner, outer)
            inner = outer
            step *= mp.mpf("1.25")
        return None

    lower = outward(-1)
    upper = outward(+1)
    return [float(lower) if lower is not None else None, float(upper) if upper is not None else None]


def next_point(previous, current, second_difference, eigenvalue):
    return [
        2 * float(current[index])
        - float(previous[index])
        + eigenvalue * float(second_difference[index])
        for index in range(len(current))
    ]


def render(scalar_path: Path, thermal_path: Path, eigenvalue: float) -> dict:
    scalar = json.loads(scalar_path.read_text(encoding="utf-8"))
    thermal = json.loads(thermal_path.read_text(encoding="utf-8"))
    if scalar.get("schema") != "matching-one/norm4-production-score/v1":
        raise ValueError("unexpected scalar score schema")
    if thermal.get("schema") != "matching-one/norm4-one-generator-thermal-jet/v1":
        raise ValueError("unexpected thermal-jet score schema")
    jordan = thermal["models_in_frozen_order"][1]
    if jordan["label"] != "Jordan_one_even_generator":
        raise ValueError("Jordan thermal score is not in the frozen second slot")
    if jordan["residual_order"] != [
        f"N{start}_r{order}" for start in (65, 85) for order in ORDERS
    ]:
        raise ValueError("unexpected Jordan residual order")

    mp.mp.dps = 80
    ratio, (chi_square, direction, error) = fit_ratio(
        jordan["residual"], jordan["covariance"]
    )
    original_chi = mp.mpf(str(jordan["score"]["chi_square"]))
    fitted = [direction[index] for index in range(5)]
    fitted_by_lineage = [fitted, [ratio * value for value in fitted]]

    predictions = []
    scalar_jordan = scalar["primary_models_in_frozen_order"][1]
    for lineage_index, (n0, n1, n2, n3) in enumerate(LINEAGES):
        thermal_points = thermal["point"]
        scalar_points = scalar["point"]
        scalar_second_difference = (
            float(scalar_points[str(n2)]["U"])
            - 2 * float(scalar_points[str(n1)]["U"])
            + float(scalar_points[str(n0)]["U"])
        )
        scalar_prediction = (
            2 * float(scalar_points[str(n2)]["U"])
            - float(scalar_points[str(n1)]["U"])
            + eigenvalue * scalar_second_difference
        )
        predictions.append(
            {
                "lineage": [n0, n1, n2, n3],
                "next_N": n3,
                "thermal_jet_orders": list(ORDERS),
                "fitted_current_second_difference": [
                    float(value) for value in fitted_by_lineage[lineage_index]
                ],
                "next_thermal_jet_point_if_lambda_half": next_point(
                    thermal_points[str(n1)],
                    thermal_points[str(n2)],
                    fitted_by_lineage[lineage_index],
                    eigenvalue,
                ),
                "observed_scalar_U_second_difference": scalar_second_difference,
                "next_scalar_U_point_if_lambda_half": scalar_prediction,
                "future_null": (
                    "x3 - 2*x2 + x1 - lambda_secondary*"
                    "(x2 - 2*x1 + x0) = 0"
                ),
            }
        )

    return {
        "schema": "matching-one/norm4-even-rank2-transfer/v1",
        "status": "post_reveal_exploratory_fit_and_frozen_future_fork",
        "inputs": {
            "scalar_score": str(scalar_path),
            "scalar_sha256": sha256(scalar_path),
            "thermal_jet_score": str(thermal_path),
            "thermal_jet_sha256": sha256(thermal_path),
        },
        "model": {
            "fixed_primary": "rank2_Jordan; affine in generation k=log2(N/N0)",
            "secondary": (
                "one conjugation-even mode; its Jordan second-difference vector "
                "is rank one across the two lineages"
            ),
            "fitted_parameters": "five mode-direction entries plus one lineage amplitude ratio",
            "degrees_of_freedom": 4,
            "identifiability": (
                "three generations identify the shared residual direction but not "
                "the secondary eigenvalue"
            ),
        },
        "covariance_aware_rank_one_fit": {
            "lineage_85_over_65_amplitude_ratio": float(ratio),
            "profile_68_percent_interval_delta_chi2_1": profile_interval(
                jordan["residual"], jordan["covariance"], ratio, chi_square, mp.mpf(1)
            ),
            "profile_95_percent_interval_delta_chi2_3p841": profile_interval(
                jordan["residual"],
                jordan["covariance"],
                ratio,
                chi_square,
                mp.mpf("3.841458820694124"),
            ),
            "thermal_mode_direction_orders_2_to_6": [float(value) for value in fitted],
            "chi_square": float(chi_square),
            "degrees_of_freedom": 4,
            "chi_square_survival": chi_square_survival(chi_square, 4),
            "improvement_over_zero_secondary_Jordan": {
                "delta_chi_square": float(original_chi - chi_square),
                "added_parameters": 6,
                "chi_square_survival_if_used_as_nested_diagnostic": chi_square_survival(
                    original_chi - chi_square, 6
                ),
            },
            "post_fit_residual": [float(value) for value in error],
        },
        "future_branch": {
            "frozen_secondary_eigenvalue": eigenvalue,
            "basis": (
                "lambda=1/2 is inherited from the ordinary analytic q=2 "
                "irrelevant mode; it is not fitted from the three generations"
            ),
            "predictions": predictions,
            "decision_table": {
                "next_second_difference_over_current_near_0": "pure rank2 Jordan reached",
                "near_1_over_2": "Jordan plus analytic even irrelevant mode",
                "near_1": "persistent curvature; rank3/log-squared candidate",
                "other_common_ratio": "different coherent even eigenmode",
                "no_common_ratio_across_orders_or_lineages": "one-secondary-mode transfer falsified",
            },
        },
        "scalar_side_view": {
            "Jordan_residuals": scalar_jordan["residual"],
            "guard": (
                "U and the thermal jet share curves but their cross-covariance is "
                "not present in the submitted score JSONs, so scalar U is predicted "
                "separately and is not added to the rank-one chi-square"
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scalar-score", type=Path, required=True)
    parser.add_argument("--thermal-score", type=Path, required=True)
    parser.add_argument("--secondary-eigenvalue", type=float, default=0.5)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = render(args.scalar_score, args.thermal_score, args.secondary_eigenvalue)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
