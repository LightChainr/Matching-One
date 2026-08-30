#!/usr/bin/env python3
"""Resolve the N85/N170/N340 same-lineage H4 amplitude recurrence."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Mapping, Optional, Sequence


LAMBDA0 = 2.0 ** (-13.0 / 8.0)


def quadratic_diagonal(residual: Sequence[float], variances: Sequence[float]) -> float:
    return math.fsum(value * value / variance for value, variance in zip(residual, variances))


def amplitude_from_pair(score: Mapping[str, object], covectors: Sequence[Fraction]) -> tuple[float, float]:
    values = score["natural_coordinate"]["value"][:2]
    covariance = [row[:2] for row in score["natural_coordinate"]["covariance"][:2]]
    dc = float(covectors[1] - covectors[0])
    weights = [-1.0 / dc, 1.0 / dc]
    amplitude = math.fsum(weight * value for weight, value in zip(weights, values))
    variance = math.fsum(
        weights[i] * covariance[i][j] * weights[j]
        for i in range(2) for j in range(2)
    )
    return amplitude, variance


def recurrence(values: Sequence[float]) -> list[float]:
    a0, a1, a2 = values
    denominator = a1 - LAMBDA0 * a0
    lambda1 = (a2 - LAMBDA0 * a1) / denominator
    c0 = (a1 - lambda1 * a0) / (LAMBDA0 - lambda1)
    c1 = a0 - c0
    a3 = c0 * LAMBDA0**3 + c1 * lambda1**3
    return [lambda1, c0, c1, a3]


def numerical_jacobian(function, values: Sequence[float]) -> list[list[float]]:
    base = list(values)
    output_size = len(function(base))
    jacobian = [[0.0 for _ in base] for _ in range(output_size)]
    for column, value in enumerate(base):
        step = max(abs(value), 1e-3) * 1e-5
        plus = list(base)
        minus = list(base)
        plus[column] += step
        minus[column] -= step
        high = function(plus)
        low = function(minus)
        for row in range(output_size):
            jacobian[row][column] = (high[row] - low[row]) / (2.0 * step)
    return jacobian


def propagate_diagonal(jacobian: Sequence[Sequence[float]], variances: Sequence[float]) -> list[list[float]]:
    return [
        [math.fsum(left[k] * variances[k] * right[k] for k in range(len(variances)))
         for right in jacobian]
        for left in jacobian
    ]


def fixed_lambda_fit(values: Sequence[float], variances: Sequence[float], lam: float) -> dict[str, object]:
    design = [1.0, lam, lam * lam]
    precision = math.fsum(x * x / variance for x, variance in zip(design, variances))
    coefficient = math.fsum(x * y / variance for x, y, variance in zip(design, values, variances)) / precision
    coefficient_variance = 1.0 / precision
    fitted = [coefficient * x for x in design]
    residual = [value - fixed for value, fixed in zip(values, fitted)]
    prediction = coefficient * lam**3
    return {
        "lambda": lam,
        "coefficient": coefficient,
        "coefficient_standard_error": math.sqrt(coefficient_variance),
        "fitted": fitted,
        "quadratic": quadratic_diagonal(residual, variances),
        "df": 2,
        "N680_H4_amplitude_prediction": prediction,
        "N680_prediction_standard_error": abs(lam**3) * math.sqrt(coefficient_variance),
    }


def free_lambda_fit(values: Sequence[float], variances: Sequence[float]) -> dict[str, object]:
    def at(lam: float) -> tuple[float, float]:
        fit = fixed_lambda_fit(values, variances, lam)
        return fit["quadratic"], fit["coefficient"]

    left, right = 1e-9, 1.0 - 1e-9
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    c = right - ratio * (right - left)
    d = left + ratio * (right - left)
    fc, _ = at(c)
    fd, _ = at(d)
    for _ in range(160):
        if fc < fd:
            right, d, fd = d, c, fc
            c = right - ratio * (right - left)
            fc, _ = at(c)
        else:
            left, c, fc = c, d, fd
            d = left + ratio * (right - left)
            fd, _ = at(d)
    lam = 0.5 * (left + right)
    quadratic, coefficient = at(lam)

    # Local nonlinear GLS covariance for parameters (coefficient, lambda).
    derivatives = [(1.0, 0.0), (lam, coefficient), (lam * lam, 2.0 * coefficient * lam)]
    information = [[math.fsum(derivatives[n][i] * derivatives[n][j] / variances[n]
                              for n in range(3)) for j in range(2)] for i in range(2)]
    determinant = information[0][0] * information[1][1] - information[0][1] ** 2
    covariance = [
        [information[1][1] / determinant, -information[0][1] / determinant],
        [-information[1][0] / determinant, information[0][0] / determinant],
    ]
    prediction_gradient = [lam**3, 3.0 * coefficient * lam * lam]
    prediction_variance = math.fsum(
        prediction_gradient[i] * covariance[i][j] * prediction_gradient[j]
        for i in range(2) for j in range(2)
    )
    return {
        "lambda": lam,
        "lambda_standard_error": math.sqrt(covariance[1][1]),
        "coefficient": coefficient,
        "parameter_covariance": covariance,
        "quadratic": quadratic,
        "df": 1,
        "N680_H4_amplitude_prediction": coefficient * lam**3,
        "N680_prediction_standard_error": math.sqrt(prediction_variance),
        "interpretation": "same-lineage descriptive transfer, not a cross-geometry exponent fit",
    }


def analyze(n85_path: Path, n170_path: Path, n340_path: Path) -> dict[str, object]:
    n85 = json.loads(n85_path.read_text(encoding="utf-8"))
    n170 = json.loads(n170_path.read_text(encoding="utf-8"))
    n340 = json.loads(n340_path.read_text(encoding="utf-8"))
    a0, v0 = amplitude_from_pair(n85, [Fraction(4633, 7225), Fraction(-6887, 7225)])
    a1 = n170["curvature_projective_decomposition"]["observed"][0]
    v1 = n170["curvature_projective_decomposition"]["measurement_covariance"][0][0]
    a2 = n340["decomposition"]["observed"][0]
    v2 = n340["decomposition"]["measurement_covariance"][0][0]
    values = [a0, a1, a2]
    variances = [v0, v1, v2]
    covariance = [[variance if i == j else 0.0 for j in range(3)] for i, variance in enumerate(variances)]

    resolved = recurrence(values)
    jacobian = numerical_jacobian(recurrence, values)
    resolved_covariance = propagate_diagonal(jacobian, variances)
    resolved_se = [math.sqrt(max(row[i], 0.0)) for i, row in enumerate(resolved_covariance)]
    lambda1, c0, c1, a3 = resolved
    denominator = a1 - LAMBDA0 * a0
    denominator_variance = v1 + LAMBDA0 * LAMBDA0 * v0
    leading = [c0 * LAMBDA0**n for n in range(4)]
    correction = [c1 * lambda1**n for n in range(4)]
    correction_ratio = [abs(b) / abs(a) for a, b in zip(leading, correction)]
    lambda_gradient = jacobian[0]
    lambda_variance_contributions = [lambda_gradient[i] ** 2 * variances[i] for i in range(3)]
    lambda_variance = resolved_covariance[0][0]

    return {
        "schema": "matching-one/P337-three-generation-H4-recurrence/v1",
        "status": "post-reveal same-lineage algebra; no new simulation",
        "source": {
            "N85": str(n85_path), "N170": str(n170_path), "N340": str(n340_path),
            "N340_reveal_commit": "e819f5e",
            "independence": "three independent seed/counter blocks; each amplitude retains its full within-pair covariance",
        },
        "sign_alignment": {
            "method": "divide each pair contrast by its exact alternating H4 covector difference",
            "result": "geometry sign is removed; A0,A1,A2 are one signed H4-amplitude lineage",
        },
        "data": {
            "generation_order": ["N85", "N170", "N340"],
            "H4_amplitude": values,
            "standard_error": [math.sqrt(value) for value in variances],
            "covariance": covariance,
        },
        "two_mode_recurrence": {
            "definition": "A_n=c0*lambda0^n+c1*lambda1^n",
            "lambda0_frozen": LAMBDA0,
            "lambda1_exact_formula": "(A2-lambda0*A1)/(A1-lambda0*A0)",
            "lambda1": lambda1,
            "lambda1_standard_error_delta": resolved_se[0],
            "lambda1_95pct_delta_interval": [lambda1 - 1.96 * resolved_se[0], lambda1 + 1.96 * resolved_se[0]],
            "point_in_open_unit_interval": 0.0 < lambda1 < 1.0,
            "z_above_zero": lambda1 / resolved_se[0],
            "z_below_one": (1.0 - lambda1) / resolved_se[0],
            "ratio_denominator_A1_minus_lambda0_A0": denominator,
            "ratio_denominator_standard_error": math.sqrt(denominator_variance),
            "ratio_denominator_z": denominator / math.sqrt(denominator_variance),
            "c0": c0, "c0_standard_error_delta": resolved_se[1],
            "c1": c1, "c1_standard_error_delta": resolved_se[2],
            "parameter_and_prediction_order": ["lambda1", "c0", "c1", "A3_N680"],
            "delta_covariance": resolved_covariance,
            "lambda1_variance_fraction_by_generation": [part / lambda_variance for part in lambda_variance_contributions],
            "leading_term_n0_to_n3": leading,
            "correction_term_n0_to_n3": correction,
            "absolute_correction_to_leading_ratio": correction_ratio,
            "N680_H4_amplitude_prediction": a3,
            "N680_prediction_standard_error_delta": resolved_se[3],
            "N680_geometry": [[22, 14], [26, 2]],
            "N680_H4_covectors_exact": ["-4633/7225", "6887/7225"],
            "N680_pair_prediction": a3 * float(Fraction(6887 + 4633, 7225)),
        },
        "comparators": {
            "single_frozen_lambda0": fixed_lambda_fit(values, variances, LAMBDA0),
            "single_free_lambda": free_lambda_fit(values, variances),
            "scale_neutral": fixed_lambda_fit(values, variances, 1.0),
            "two_mode": {
                "quadratic": 0.0, "df": 0,
                "N680_H4_amplitude_prediction": a3,
                "N680_prediction_standard_error": resolved_se[3],
                "warning": "three parameters interpolate three generations; use the N680 prediction, not in-sample q, as the discriminator",
            },
        },
        "reading": {
            "point_estimate": "lambda1 is a decaying opposite-sign correction faster than frozen lambda0",
            "uncertainty": "lambda1 lies in (0,1) at the point estimate but is not resolved above zero by delta propagation",
            "cancellation": "the correction cancels the leading H4 term and its relative magnitude decreases generation by generation",
            "dominant_lambda1_uncertainty_source": "N340",
        },
        "claim_boundary": "same-lineage three-generation mechanism coordinate and N680 prediction; no cross-geometry exponent fit and no unique asymptotic correction claim",
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    model = payload["two_mode_recurrence"]
    comparators = payload["comparators"]
    lines = [
        "# Three-generation same-lineage H4 recurrence", "",
        "No new simulation is used. Exact alternating H4 geometry signs are divided out before analysis.", "",
        f"With `lambda0=2^(-13/8)={model['lambda0_frozen']:.9g}`, the exact three-point identity gives `lambda1={model['lambda1']:.6f} +/- {model['lambda1_standard_error_delta']:.6f}` (delta method). The point lies in `(0,1)`, but is only `{model['z_above_zero']:.3f}` SE above zero; the 95% delta interval is `[{model['lambda1_95pct_delta_interval'][0]:.3f},{model['lambda1_95pct_delta_interval'][1]:.3f}]`.", "",
        f"The amplitudes are `c0={model['c0']:+.6f}` and `c1={model['c1']:+.6f}`. They have opposite signs. The correction/leading magnitude ratio falls `{model['absolute_correction_to_leading_ratio'][0]:.3f} -> {model['absolute_correction_to_leading_ratio'][1]:.3f} -> {model['absolute_correction_to_leading_ratio'][2]:.3f} -> {model['absolute_correction_to_leading_ratio'][3]:.3f}` from N85 through predicted N680.", "",
        f"Two-mode N680 prediction: `A_H={model['N680_H4_amplitude_prediction']:+.7f} +/- {model['N680_prediction_standard_error_delta']:.3g}`; its exact child geometry makes the pair negative, `{model['N680_pair_prediction']:+.7f}`.", "",
        "| model | in-sample q/df | N680 A_H | prediction SE |", "|---|---:|---:|---:|",
        f"| frozen single lambda0 | {comparators['single_frozen_lambda0']['quadratic']:.3f}/2 | {comparators['single_frozen_lambda0']['N680_H4_amplitude_prediction']:+.7f} | {comparators['single_frozen_lambda0']['N680_prediction_standard_error']:.3g} |",
        f"| free single lambda | {comparators['single_free_lambda']['quadratic']:.3f}/1 | {comparators['single_free_lambda']['N680_H4_amplitude_prediction']:+.7f} | {comparators['single_free_lambda']['N680_prediction_standard_error']:.3g} |",
        f"| scale-neutral | {comparators['scale_neutral']['quadratic']:.3f}/2 | {comparators['scale_neutral']['N680_H4_amplitude_prediction']:+.7f} | {comparators['scale_neutral']['N680_prediction_standard_error']:.3g} |",
        f"| fixed-lambda0 plus correction | interpolation/0 | {model['N680_H4_amplitude_prediction']:+.7f} | {model['N680_prediction_standard_error_delta']:.3g} |", "",
        "The useful new coordinate is a point-estimate fast, opposite-sign correction that naturally produces the N170 overshoot and N340 return. Its sign/rate is not yet resolved: about 80% of lambda1 variance comes from N340. N680, not the zero-df interpolation, is the clean discriminator. No exponent is fitted across unrelated geometries.", "",
    ]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n85", type=Path, required=True)
    parser.add_argument("--n170", type=Path, required=True)
    parser.add_argument("--n340", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    payload = analyze(args.n85, args.n170, args.n340)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
