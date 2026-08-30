#!/usr/bin/env python3
"""Compare fixed identity dressing and minimal rank-3 Jordan adversaries."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Mapping, Optional, Sequence

from analyze_three_generation_recurrence import amplitude_from_pair
from fit_four_generation_recurrence import (
    LAMBDA0,
    chi_square_survival,
    fixed_single,
    free_single,
    golden_minimize,
    inverse,
    quadratic,
    two_mode_at,
)


LAMBDA_ID4 = 2.0 ** (-21.0 / 8.0)


def matvec(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> list[float]:
    return [math.fsum(a * b for a, b in zip(row, vector)) for row in matrix]


def linear_fit(
    values: Sequence[float], variances: Sequence[float],
    design: Sequence[Sequence[float]], forecast_row: Sequence[float],
) -> tuple[dict[str, object], list[float]]:
    width = len(design[0])
    information = [[math.fsum(row[i] * row[j] / variance
                              for row, variance in zip(design, variances))
                    for j in range(width)] for i in range(width)]
    covariance = inverse(information)
    rhs = [math.fsum(row[i] * value / variance
                     for row, value, variance in zip(design, values, variances))
           for i in range(width)]
    coefficients = matvec(covariance, rhs)
    fitted = [math.fsum(a * b for a, b in zip(row, coefficients)) for row in design]
    residual = [value - fixed for value, fixed in zip(values, fitted)]
    q = quadratic(residual, variances)
    influence = [
        math.fsum(forecast_row[i] * covariance[i][j] for i in range(width))
        for j in range(width)
    ]
    data_gradient = [
        math.fsum(influence[j] * row[j] for j in range(width)) / variance
        for row, variance in zip(design, variances)
    ]
    forecast = math.fsum(a * b for a, b in zip(forecast_row, coefficients))
    forecast_variance = math.fsum(weight * weight * variance
                                  for weight, variance in zip(data_gradient, variances))
    df = len(values) - width
    return ({
        "coefficients": coefficients,
        "parameter_covariance": covariance,
        "quadratic": q, "df": df, "gof_p": chi_square_survival(q, df),
        "AIC_descriptive": q + 2.0 * width,
        "N1360_H4_amplitude_prediction": forecast,
        "N1360_prediction_standard_error": math.sqrt(forecast_variance),
        "N1360_data_influence": data_gradient,
    }, data_gradient)


def recurrence_fit(values: Sequence[float], variances: Sequence[float]) -> tuple[dict[str, object], list[float]]:
    lambda1 = golden_minimize(lambda lam: two_mode_at(values, variances, lam)[0], 1e-9, 1.0 - 1e-9)
    q, coefficients = two_mode_at(values, variances, lambda1)
    c0, c1 = coefficients

    def forecast(data: Sequence[float]) -> float:
        fitted_lambda = golden_minimize(
            lambda lam: two_mode_at(data, variances, lam)[0], 1e-9, 1.0 - 1e-9
        )
        _, (left, right) = two_mode_at(data, variances, fitted_lambda)
        return left * LAMBDA0**4 + right * fitted_lambda**4

    data_gradient = []
    for index, value in enumerate(values):
        step = max(abs(value), 1e-3) * 1e-4
        plus, minus = list(values), list(values)
        plus[index] += step
        minus[index] -= step
        data_gradient.append((forecast(plus) - forecast(minus)) / (2.0 * step))
    prediction_variance = math.fsum(weight * weight * variance
                                    for weight, variance in zip(data_gradient, variances))
    model = {
        "lambda1": lambda1, "coefficients": coefficients,
        "quadratic": q, "df": 1, "gof_p": chi_square_survival(q, 1),
        "AIC_descriptive": q + 6.0,
        "N1360_H4_amplitude_prediction": forecast(values),
        "N1360_prediction_standard_error": math.sqrt(prediction_variance),
        "N1360_data_influence": data_gradient,
        "influence_method": "central refit derivative with full four-generation covariance",
    }
    return model, data_gradient


def free_single_fit(values: Sequence[float], variances: Sequence[float]) -> tuple[dict[str, object], list[float]]:
    base = free_single(values, variances)

    def forecast(data: Sequence[float]) -> float:
        return free_single(data, variances)["N1360_H4_amplitude_prediction"]

    data_gradient = []
    for index, value in enumerate(values):
        step = max(abs(value), 1e-3) * 1e-4
        plus, minus = list(values), list(values)
        plus[index] += step
        minus[index] -= step
        data_gradient.append((forecast(plus) - forecast(minus)) / (2.0 * step))
    prediction_variance = math.fsum(weight * weight * variance
                                    for weight, variance in zip(data_gradient, variances))
    base["N1360_H4_amplitude_prediction"] = forecast(values)
    base["N1360_prediction_standard_error"] = math.sqrt(prediction_variance)
    base["N1360_data_influence"] = data_gradient
    base["influence_method"] = "central refit derivative with full four-generation covariance"
    return base, data_gradient


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
    generations = list(range(4))

    models: dict[str, dict[str, object]] = {}
    gradients: dict[str, list[float]] = {}
    identity_design = [[LAMBDA0**n, LAMBDA_ID4**n] for n in generations]
    identity, gradients["fixed_identity_dressing"] = linear_fit(
        values, variances, identity_design, [LAMBDA0**4, LAMBDA_ID4**4]
    )
    identity.update({
        "definition": "c0*lambda0^n+c_id4*lambda_id4^n",
        "lambda0": LAMBDA0, "lambda_id4": LAMBDA_ID4,
        "coefficient_order": ["c0", "c_id4"],
        "absolute_dressing_to_leading_ratio_n0_to_n4": [
            abs(identity["coefficients"][1] * LAMBDA_ID4**n)
            / abs(identity["coefficients"][0] * LAMBDA0**n)
            for n in range(5)
        ],
    })
    models["fixed_identity_dressing"] = identity

    jordan_design = [[LAMBDA0**n, n * LAMBDA0**n,
                       n * (n - 1) / 2.0 * LAMBDA0**n] for n in generations]
    jordan, gradients["rank3_same_base_jordan"] = linear_fit(
        values, variances, jordan_design,
        [LAMBDA0**4, 4.0 * LAMBDA0**4, 6.0 * LAMBDA0**4],
    )
    jordan.update({
        "definition": "lambda0^n*(c0+c1*n+c2*n*(n-1)/2)",
        "lambda0": LAMBDA0,
        "coefficient_order": ["c0", "c1", "c2"],
        "rank": 3,
    })
    models["rank3_same_base_jordan"] = jordan

    recurrence, gradients["free_lambda_recurrence"] = recurrence_fit(values, variances)
    recurrence["definition"] = "c0*lambda0^n+c1*lambda1^n; lambda1 fitted"
    models["free_lambda_recurrence"] = recurrence
    free, gradients["free_single_lambda"] = free_single_fit(values, variances)
    models["free_single_lambda"] = free

    fixed = fixed_single(values, variances, LAMBDA0)
    fixed_design = [[LAMBDA0**n] for n in generations]
    _, gradients["fixed_single_lambda0"] = linear_fit(
        values, variances, fixed_design, [LAMBDA0**4]
    )
    fixed["N1360_data_influence"] = gradients["fixed_single_lambda0"]
    models["fixed_single_lambda0"] = fixed

    neutral = fixed_single(values, variances, 1.0)
    neutral_design = [[1.0] for _ in generations]
    _, gradients["scale_neutral"] = linear_fit(values, variances, neutral_design, [1.0])
    neutral["N1360_data_influence"] = gradients["scale_neutral"]
    models["scale_neutral"] = neutral

    order = list(models)
    forecast_covariance = [
        [math.fsum(gradients[left][k] * variances[k] * gradients[right][k]
                   for k in range(4)) for right in order]
        for left in order
    ]
    forecasts = [models[name]["N1360_H4_amplitude_prediction"] for name in order]
    pairwise = []
    identity_index = order.index("fixed_identity_dressing")
    for right_index, name in enumerate(order):
        if right_index == identity_index:
            continue
        difference = forecasts[identity_index] - forecasts[right_index]
        difference_variance = (
            forecast_covariance[identity_index][identity_index]
            + forecast_covariance[right_index][right_index]
            - 2.0 * forecast_covariance[identity_index][right_index]
        )
        source_se = math.sqrt(max(difference_variance, 0.0))
        maximum_z = abs(difference) / source_se
        measurement_variance_ceiling = (abs(difference) / 3.0) ** 2 - difference_variance
        pairwise.append({
            "left": "fixed_identity_dressing", "right": name,
            "forecast_difference": difference,
            "source_difference_standard_error": source_se,
            "maximum_source_limited_z": maximum_z,
            "future_measurement_standard_error_ceiling_for_3sigma": (
                math.sqrt(measurement_variance_ceiling)
                if measurement_variance_ceiling > 0.0 else None
            ),
            "three_sigma_possible_without_refitting_sources": measurement_variance_ceiling > 0.0,
        })

    aic_min = min(model["AIC_descriptive"] for model in models.values())
    for model in models.values():
        model["delta_AIC_descriptive"] = model["AIC_descriptive"] - aic_min
    return {
        "schema": "matching-one/P337-theory-fixed-identity-dressing-adversary/v1",
        "status": "existing-data theory adversary; no new simulation",
        "source": {
            "N85": str(n85_path), "N170": str(n170_path),
            "N340": str(n340_path), "N680": str(n680_path),
            "independence": "four independent random blocks; full within-pair covariance projected before diagonal inter-generation fit",
        },
        "data": {
            "generation_order": ["N85", "N170", "N340", "N680"],
            "H4_amplitude": values,
            "covariance": [[variance if i == j else 0.0 for j in range(4)]
                           for i, variance in enumerate(variances)],
        },
        "models_in_fixed_order": order,
        "models": models,
        "N1360_forecast_joint": {
            "order": order,
            "value": forecasts,
            "covariance": forecast_covariance,
            "pairwise_identity_dressing_separation": pairwise,
            "geometry": [[36, 8], [28, 24]],
            "Smith_classes": [[4, 340], [4, 340]],
            "production_status": "not authorized or started",
        },
        "reading": {
            "identity_dressing": "theory-fixed identity dressing passes two-df GOF and has the lowest descriptive AIC",
            "rank3": "minimal same-base rank-3 Jordan also passes but is forecast-degenerate with free-lambda recurrence",
            "N1360": "one N1360 measurement cannot broadly separate identity dressing from recurrence, rank-3, or free-single under current source covariance; it can target fixed-single and neutral",
            "AIC_boundary": "AIC is descriptive; fixed-theory GOF and preregistered lineage behavior carry the mechanism claim",
        },
        "claim_boundary": "same-lineage theory adversary and N1360 identifiability audit; no field identity, new exponent, or production authorization",
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    lines = [
        "# Theory-fixed identity-dressing adversary", "",
        "No new simulation is used. The identity-dressing eigenvalues are frozen to `2^-13/8` and `2^-21/8`; only their amplitudes are fitted.", "",
        "| model | q/df | GOF p | descriptive AIC | delta AIC | N1360 A_H | target SE |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for name in payload["models_in_fixed_order"]:
        model = payload["models"][name]
        lines.append(
            f"| {name} | {model['quadratic']:.3f}/{model['df']} | {model['gof_p']:.3g} | "
            f"{model['AIC_descriptive']:.3f} | {model['delta_AIC_descriptive']:.3f} | "
            f"{model['N1360_H4_amplitude_prediction']:+.7f} | {model['N1360_prediction_standard_error']:.3g} |"
        )
    identity = payload["models"]["fixed_identity_dressing"]
    lines.extend([
        "",
        "The fixed identity dressing passes GOF and is the descriptive AIC leader. Its dressing/leading magnitude halves each generation by construction and is `" + " -> ".join(f"{value:.3f}" for value in identity["absolute_dressing_to_leading_ratio_n0_to_n4"]) + "`.", "",
        "N1360 source-covariance ceiling against identity dressing:", "",
    ])
    for row in payload["N1360_forecast_joint"]["pairwise_identity_dressing_separation"]:
        ceiling = row["future_measurement_standard_error_ceiling_for_3sigma"]
        lines.append(
            f"- `{row['right']}`: maximum `{row['maximum_source_limited_z']:.3f}` sigma; "
            + (f"3-sigma requires measurement SE below `{ceiling:.3g}`."
               if ceiling is not None else "3-sigma is impossible without reducing source uncertainty.")
        )
    lines.extend([
        "",
        "Therefore N1360 is not yet a universal discriminator: recurrence and rank-3 Jordan are nearly forecast-identical, and free-single remains source-limited. It can efficiently reject neutral, while separating fixed-single would require unusually small measurement error. No N1360 production is started.", "",
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
