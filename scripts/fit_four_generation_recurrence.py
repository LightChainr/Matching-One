#!/usr/bin/env python3
"""Fit the fixed-leading two-mode recurrence to N85 through N680."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Mapping, Optional, Sequence

from analyze_three_generation_recurrence import amplitude_from_pair


LAMBDA0 = 2.0 ** (-13.0 / 8.0)
PROFILE_95_DELTA = 3.841458820694124


def solve_two_by_two(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> list[float]:
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    if abs(determinant) < 1e-24:
        raise ValueError("singular two-column design")
    return [
        (rhs[0] * matrix[1][1] - matrix[0][1] * rhs[1]) / determinant,
        (matrix[0][0] * rhs[1] - rhs[0] * matrix[1][0]) / determinant,
    ]


def inverse(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    size = len(matrix)
    augmented = [list(row) + [1.0 if i == j else 0.0 for j in range(size)]
                 for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1e-20:
            raise ValueError("singular information matrix")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [value - factor * fixed
                              for value, fixed in zip(augmented[row], augmented[column])]
    return [row[size:] for row in augmented]


def quadratic(residual: Sequence[float], variances: Sequence[float]) -> float:
    return math.fsum(value * value / variance for value, variance in zip(residual, variances))


def chi_square_survival(q: float, df: int) -> float:
    if df == 1:
        return math.erfc(math.sqrt(q / 2.0))
    if df == 2:
        return math.exp(-q / 2.0)
    if df == 3:
        t = math.sqrt(q / 2.0)
        return math.erfc(t) + 2.0 / math.sqrt(math.pi) * t * math.exp(-q / 2.0)
    raise ValueError("only df 1, 2, 3 are used")


def two_mode_at(values: Sequence[float], variances: Sequence[float], lambda1: float) -> tuple[float, list[float]]:
    rows = [(LAMBDA0**n, lambda1**n) for n in range(len(values))]
    information = [[math.fsum(row[i] * row[j] / variance for row, variance in zip(rows, variances))
                    for j in range(2)] for i in range(2)]
    rhs = [math.fsum(row[i] * value / variance
                     for row, value, variance in zip(rows, values, variances)) for i in range(2)]
    coefficients = solve_two_by_two(information, rhs)
    residual = [value - coefficients[0] * row[0] - coefficients[1] * row[1]
                for value, row in zip(values, rows)]
    return quadratic(residual, variances), coefficients


def golden_minimize(function, left: float, right: float) -> float:
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    c = right - ratio * (right - left)
    d = left + ratio * (right - left)
    fc, fd = function(c), function(d)
    for _ in range(180):
        if fc < fd:
            right, d, fd = d, c, fc
            c = right - ratio * (right - left)
            fc = function(c)
        else:
            left, c, fc = c, d, fd
            d = left + ratio * (right - left)
            fd = function(d)
    return 0.5 * (left + right)


def bisect(function, left: float, right: float) -> float:
    f_left, f_right = function(left), function(right)
    if f_left * f_right > 0.0:
        raise ValueError("profile endpoint does not bracket a root")
    for _ in range(160):
        middle = 0.5 * (left + right)
        f_middle = function(middle)
        if f_left * f_middle <= 0.0:
            right, f_right = middle, f_middle
        else:
            left, f_left = middle, f_middle
    return 0.5 * (left + right)


def fixed_single(values: Sequence[float], variances: Sequence[float], lam: float) -> dict[str, object]:
    design = [lam**n for n in range(len(values))]
    precision = math.fsum(x * x / variance for x, variance in zip(design, variances))
    coefficient = math.fsum(x * value / variance for x, value, variance in zip(design, values, variances)) / precision
    coefficient_variance = 1.0 / precision
    residual = [value - coefficient * x for value, x in zip(values, design)]
    q = quadratic(residual, variances)
    prediction = coefficient * lam**4
    return {
        "lambda": lam, "coefficient": coefficient,
        "coefficient_standard_error": math.sqrt(coefficient_variance),
        "quadratic": q, "df": 3, "gof_p": chi_square_survival(q, 3),
        "AIC_descriptive": q + 2.0,
        "N1360_H4_amplitude_prediction": prediction,
        "N1360_prediction_standard_error": abs(lam**4) * math.sqrt(coefficient_variance),
    }


def free_single(values: Sequence[float], variances: Sequence[float]) -> dict[str, object]:
    objective = lambda lam: fixed_single(values, variances, lam)["quadratic"]
    lam = golden_minimize(objective, 1e-9, 1.0 - 1e-9)
    base = fixed_single(values, variances, lam)
    coefficient = base["coefficient"]
    rows = [(lam**n, coefficient * n * lam ** (n - 1) if n else 0.0)
            for n in range(len(values))]
    information = [[math.fsum(row[i] * row[j] / variance for row, variance in zip(rows, variances))
                    for j in range(2)] for i in range(2)]
    covariance = inverse(information)
    gradient = [lam**4, 4.0 * coefficient * lam**3]
    prediction_variance = math.fsum(gradient[i] * covariance[i][j] * gradient[j]
                                    for i in range(2) for j in range(2))
    q = base["quadratic"]
    return {
        "lambda": lam, "lambda_standard_error": math.sqrt(covariance[1][1]),
        "coefficient": coefficient, "parameter_covariance": covariance,
        "quadratic": q, "df": 2, "gof_p": chi_square_survival(q, 2),
        "AIC_descriptive": q + 4.0,
        "N1360_H4_amplitude_prediction": coefficient * lam**4,
        "N1360_prediction_standard_error": math.sqrt(prediction_variance),
        "interpretation": "same-lineage transfer only; not a cross-geometry exponent",
    }


def analyze(n85_path: Path, n170_path: Path, n340_path: Path, n680_path: Path) -> dict[str, object]:
    scores = [json.loads(path.read_text(encoding="utf-8"))
              for path in (n85_path, n170_path, n340_path, n680_path)]
    a0, v0 = amplitude_from_pair(scores[0], [Fraction(4633, 7225), Fraction(-6887, 7225)])
    a1 = scores[1]["curvature_projective_decomposition"]["observed"][0]
    v1 = scores[1]["curvature_projective_decomposition"]["measurement_covariance"][0][0]
    a2 = scores[2]["decomposition"]["observed"][0]
    v2 = scores[2]["decomposition"]["measurement_covariance"][0][0]
    a3 = scores[3]["decomposition"]["observed"][0]
    v3 = scores[3]["decomposition"]["measurement_covariance"][0][0]
    values, variances = [a0, a1, a2, a3], [v0, v1, v2, v3]

    objective = lambda lam: two_mode_at(values, variances, lam)[0]
    lambda1 = golden_minimize(objective, 1e-9, 1.0 - 1e-9)
    q, (c0, c1) = two_mode_at(values, variances, lambda1)
    rows = [
        (LAMBDA0**n, lambda1**n, c1 * n * lambda1 ** (n - 1) if n else 0.0)
        for n in range(4)
    ]
    information = [[math.fsum(row[i] * row[j] / variance for row, variance in zip(rows, variances))
                    for j in range(3)] for i in range(3)]
    covariance = inverse(information)
    standard_error = [math.sqrt(covariance[i][i]) for i in range(3)]
    profile_target = q + PROFILE_95_DELTA
    profile_function = lambda lam: objective(lam) - profile_target
    profile_interval = [
        bisect(profile_function, 1e-12, lambda1),
        bisect(profile_function, lambda1, 1.0 - 1e-9),
    ]
    leading = [c0 * LAMBDA0**n for n in range(5)]
    correction = [c1 * lambda1**n for n in range(5)]
    ratios = [abs(right) / abs(left) for left, right in zip(leading, correction)]
    prediction = leading[4] + correction[4]
    prediction_gradient = [LAMBDA0**4, lambda1**4, 4.0 * c1 * lambda1**3]
    prediction_variance = math.fsum(
        prediction_gradient[i] * covariance[i][j] * prediction_gradient[j]
        for i in range(3) for j in range(3)
    )

    frozen = fixed_single(values, variances, LAMBDA0)
    free = free_single(values, variances)
    neutral = fixed_single(values, variances, 1.0)
    recurrence = {
        "definition": "A_n=c0*lambda0^n+c1*lambda1^n",
        "lambda0_frozen": LAMBDA0,
        "lambda1": lambda1, "lambda1_standard_error_wald": standard_error[2],
        "lambda1_95pct_wald_interval": [lambda1 - 1.96 * standard_error[2], lambda1 + 1.96 * standard_error[2]],
        "lambda1_95pct_profile_interval": profile_interval,
        "profile_boundary_note": "lower endpoint is close to zero and should be read as marginal exclusion, not a robust sign discovery",
        "c0": c0, "c0_standard_error": standard_error[0],
        "c1": c1, "c1_standard_error": standard_error[1],
        "parameter_order": ["c0", "c1", "lambda1"],
        "parameter_covariance": covariance,
        "quadratic": q, "df": 1, "gof_p": chi_square_survival(q, 1),
        "AIC_descriptive": q + 6.0,
        "leading_term_n0_to_n4": leading,
        "correction_term_n0_to_n4": correction,
        "absolute_correction_to_leading_ratio": ratios,
        "N1360_H4_amplitude_prediction": prediction,
        "N1360_prediction_standard_error": math.sqrt(prediction_variance),
    }
    models = {
        "fixed_lambda0_single": frozen,
        "free_single_lambda": free,
        "scale_neutral": neutral,
        "fixed_lambda0_plus_correction": recurrence,
    }
    aic_min = min(model["AIC_descriptive"] for model in models.values())
    for model in models.values():
        model["delta_AIC_descriptive"] = model["AIC_descriptive"] - aic_min

    covectors = [Fraction(4633, 7225), Fraction(-6887, 7225)]
    return {
        "schema": "matching-one/P337-four-generation-H4-recurrence-fit/v1",
        "status": "post-heldout covariance-aware same-lineage fit; no new simulation",
        "source": {
            "N85": str(n85_path), "N170": str(n170_path),
            "N340": str(n340_path), "N680": str(n680_path),
            "N680_reveal_commit": "02080a4",
            "independence": "four independent seed/counter blocks; full within-pair covariance retained before H4 projection",
        },
        "data": {
            "generation_order": ["N85", "N170", "N340", "N680"],
            "H4_amplitude": values,
            "standard_error": [math.sqrt(value) for value in variances],
            "covariance": [[variance if i == j else 0.0 for j in range(4)] for i, variance in enumerate(variances)],
        },
        "models": models,
        "N1360_forecast_freeze": {
            "status": "frozen for a future heldout; production not authorized or started",
            "geometry": [[36, 8], [28, 24]],
            "period_matrices": [[[36, -8], [8, 36]], [[28, -24], [24, 28]]],
            "Smith_classes": [[4, 340], [4, 340]],
            "H4_covectors_exact": [str(value) for value in covectors],
            "pair_expected_sign": "positive",
            "model_order": list(models),
            "predictions": {
                name: {
                    "H4_amplitude": model["N1360_H4_amplitude_prediction"],
                    "standard_error": model["N1360_prediction_standard_error"],
                    "pair_second_minus_first": model["N1360_H4_amplitude_prediction"] * float(covectors[1] - covectors[0]),
                }
                for name, model in models.items()
            },
            "projective_scalar_target": 0.0,
        },
        "reading": {
            "heldout_GOF": "the recurrence passes the fourth-generation one-df GOF; fixed single H4 and scale-neutral fail, free-single remains acceptable",
            "AIC_boundary": "AIC is descriptive only and leaves recurrence/free-single tied; the heldout residual and GOF carry the scientific weight",
            "lambda1": "positive decaying correction is now marginally separated from zero by profile likelihood, but the boundary is too close for a strong sign claim",
        },
        "claim_boundary": "four-generation same-lineage recurrence GOF and frozen N1360 forecasts; no cross-geometry exponent or unique correction-field identification",
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    recurrence = payload["models"]["fixed_lambda0_plus_correction"]
    lines = [
        "# Four-generation covariance-aware H4 recurrence", "",
        "No new simulation is used. Exact H4 geometry signs are removed before fitting the one Gaussian lineage.", "",
        f"Fixed `lambda0=2^-13/8`; fit `lambda1={recurrence['lambda1']:.6f} +/- {recurrence['lambda1_standard_error_wald']:.6f}` with 95% profile interval `[{recurrence['lambda1_95pct_profile_interval'][0]:.6g},{recurrence['lambda1_95pct_profile_interval'][1]:.6g}]`. The lower endpoint is boundary-close and only marginal evidence for positive lambda1.", "",
        f"Recurrence GOF is `{recurrence['quadratic']:.3f}/1` (`p={recurrence['gof_p']:.3f}`). Correction/leading magnitude falls `" + " -> ".join(f"{value:.3f}" for value in recurrence["absolute_correction_to_leading_ratio"]) + "` through predicted N1360.", "",
        "| model | q/df | GOF p | descriptive AIC | delta AIC | N1360 A_H | SE |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for name, model in payload["models"].items():
        lines.append(
            f"| {name} | {model['quadratic']:.3f}/{model['df']} | {model['gof_p']:.3g} | "
            f"{model['AIC_descriptive']:.3f} | {model['delta_AIC_descriptive']:.3f} | "
            f"{model['N1360_H4_amplitude_prediction']:+.7f} | {model['N1360_prediction_standard_error']:.3g} |"
        )
    lines.extend([
        "",
        "The recurrence passes the genuinely heldout fourth-generation shape. Free-single also has acceptable GOF and is essentially AIC-tied; fixed single H4 and scale-neutral fail. AIC is descriptive, not the core claim.", "",
        "N1360 forecasts and exact `(4,340)` child geometry are frozen in the machine-readable result. No N1360 production is authorized or running.", "",
    ])
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n85", type=Path, required=True)
    parser.add_argument("--n170", type=Path, required=True)
    parser.add_argument("--n340", type=Path, required=True)
    parser.add_argument("--n680", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    payload = analyze(args.n85, args.n170, args.n340, args.n680)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
