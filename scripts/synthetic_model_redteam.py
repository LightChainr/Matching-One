#!/usr/bin/env python3
"""Synthetic red-team for the frozen same-N orientation model score.

The benchmark deliberately uses the five P31/P32 sizes, their reported
pooled standard errors and the frozen 3+2 train/held-out split.  It asks a
narrow question: if one of six declared mean laws generated fresh data at the
current noise level, which law would the present scoring design select?

This is a design-power calculation, not evidence about which law is true.
Alternative components are normalized to a declared fraction of the H4
signal in the covariance metric so that the stress test is reproducible and
does not hide arbitrary unit choices.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
import math
from pathlib import Path
import random
from typing import Dict, Iterable, List, Sequence, Tuple


MODEL_ORDER = (
    "H4",
    "H12",
    "H4+H12",
    "two_radial_powers",
    "log_Jordan",
    "ordinary_correction",
)

Vector = List[float]
Matrix = List[List[float]]


@dataclass(frozen=True)
class Observation:
    n: int
    first: Tuple[int, int]
    second: Tuple[int, int]
    se: float
    split: str
    delta_cos4: float
    delta_cos12: float


def gaussian_multiply(left: Tuple[int, int], right: Tuple[int, int]) -> Tuple[int, int]:
    a, b = left
    c, d = right
    return a * c - b * d, a * d + b * c


def gaussian_power(value: Tuple[int, int], exponent: int) -> Tuple[int, int]:
    output = (1, 0)
    base = value
    while exponent:
        if exponent & 1:
            output = gaussian_multiply(output, base)
        base = gaussian_multiply(base, base)
        exponent //= 2
    return output


def cosine_harmonic(value: Tuple[int, int], spin: int) -> Fraction:
    if spin <= 0 or spin % 4:
        raise ValueError("spin must be a positive multiple of four")
    norm = value[0] * value[0] + value[1] * value[1]
    if norm == 0:
        raise ValueError("zero Gaussian integer")
    real, _imaginary = gaussian_power(value, spin)
    return Fraction(real, norm ** (spin // 2))


def load_config(path: Path) -> Tuple[dict, List[Observation]]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    observations = []
    for row in raw["observations"]:
        first = tuple(map(int, row["first"]))
        second = tuple(map(int, row["second"]))
        n = int(row["N"])
        if sum(value * value for value in first) != n or sum(value * value for value in second) != n:
            raise ValueError("Gaussian representations must have the declared common norm")
        observations.append(Observation(
            n=n,
            first=first,
            second=second,
            se=float(row["standard_error"]),
            split=str(row["split"]),
            delta_cos4=float(cosine_harmonic(first, 4) - cosine_harmonic(second, 4)),
            delta_cos12=float(cosine_harmonic(first, 12) - cosine_harmonic(second, 12)),
        ))
    if {row.split for row in observations} != {"train", "heldout"}:
        raise ValueError("config must contain train and heldout rows")
    if any(row.se <= 0 for row in observations):
        raise ValueError("standard errors must be positive")
    return raw, observations


def transpose(matrix: Matrix) -> Matrix:
    return [list(column) for column in zip(*matrix)]


def matmul(first: Matrix, second: Matrix) -> Matrix:
    second_t = transpose(second)
    return [[math.fsum(a * b for a, b in zip(row, column)) for column in second_t]
            for row in first]


def matvec(matrix: Matrix, vector: Vector) -> Vector:
    return [math.fsum(a * b for a, b in zip(row, vector)) for row in matrix]


def solve(matrix: Matrix, vector: Vector) -> Vector:
    size = len(vector)
    work = [list(map(float, matrix[index])) + [float(vector[index])]
            for index in range(size)]
    scale = max(abs(value) for row in matrix for value in row)
    tolerance = max(scale * 1e-13, 1e-300)
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(work[row][column]))
        if abs(work[pivot][column]) <= tolerance:
            raise ArithmeticError("singular model design")
        work[column], work[pivot] = work[pivot], work[column]
        divisor = work[column][column]
        for entry in range(column, size + 1):
            work[column][entry] /= divisor
        for row in range(size):
            if row == column:
                continue
            factor = work[row][column]
            for entry in range(column, size + 1):
                work[row][entry] -= factor * work[column][entry]
    return [work[row][-1] for row in range(size)]


def inverse(matrix: Matrix) -> Matrix:
    size = len(matrix)
    columns = [solve(matrix, [float(row == column) for row in range(size)])
               for column in range(size)]
    return [[columns[column][row] for column in range(size)] for row in range(size)]


def determinant(matrix: Matrix) -> float:
    size = len(matrix)
    work = [row[:] for row in matrix]
    output = 1.0
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(work[row][column]))
        if work[pivot][column] == 0:
            return 0.0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            output *= -1.0
        diagonal = work[column][column]
        output *= diagonal
        for row in range(column + 1, size):
            factor = work[row][column] / diagonal
            for entry in range(column + 1, size):
                work[row][entry] -= factor * work[column][entry]
    return output


def quadratic(vector: Vector, precision: Matrix) -> float:
    return math.fsum(a * b for a, b in zip(vector, matvec(precision, vector)))


def basis_rows(model: str, observations: Sequence[Observation], reference_n: float) -> Matrix:
    rows = []
    for observation in observations:
        n = float(observation.n)
        h4 = observation.delta_cos4 * n ** (-13.0 / 8.0)
        h12 = observation.delta_cos12 * n ** (-13.0 / 8.0)
        if model == "H4":
            row = [h4]
        elif model == "H12":
            row = [h12]
        elif model == "H4+H12":
            row = [h4, h12]
        elif model == "two_radial_powers":
            row = [h4, observation.delta_cos4 * n ** (-9.0 / 8.0)]
        elif model == "log_Jordan":
            row = [h4, h4 * math.log(n / reference_n)]
        elif model == "ordinary_correction":
            row = [h4, observation.delta_cos4 * n ** (-21.0 / 8.0)]
        else:
            raise ValueError("unknown model " + model)
        rows.append(row)
    return rows


def weighted_norm(vector: Sequence[float], observations: Sequence[Observation]) -> float:
    return math.sqrt(math.fsum((value / row.se) ** 2 for value, row in zip(vector, observations)))


def column(matrix: Matrix, index: int) -> Vector:
    return [row[index] for row in matrix]


def scaled_second_coefficient(
    base_mean: Vector,
    second_basis: Vector,
    observations: Sequence[Observation],
    fraction: float,
) -> float:
    denominator = weighted_norm(second_basis, observations)
    if denominator == 0:
        raise ArithmeticError("zero alternative component")
    return fraction * weighted_norm(base_mean, observations) / denominator


def truth_means(
    observations: Sequence[Observation], reference_n: float, amplitude: float, fraction: float
) -> Tuple[Dict[str, Vector], Dict[str, Vector]]:
    designs = {model: basis_rows(model, observations, reference_n) for model in MODEL_ORDER}
    h4_basis = column(designs["H4"], 0)
    h4_mean = [amplitude * value for value in h4_basis]
    output: Dict[str, Vector] = {"H4": h4_mean}
    coefficients: Dict[str, Vector] = {"H4": [amplitude]}

    h12_basis = column(designs["H12"], 0)
    a12 = scaled_second_coefficient(h4_mean, h12_basis, observations, 1.0)
    output["H12"] = [a12 * value for value in h12_basis]
    coefficients["H12"] = [a12]

    for model in MODEL_ORDER[2:]:
        first = column(designs[model], 0)
        second = column(designs[model], 1)
        base = [amplitude * value for value in first]
        coefficient = scaled_second_coefficient(base, second, observations, fraction)
        output[model] = [left + coefficient * right for left, right in zip(base, second)]
        coefficients[model] = [amplitude, coefficient]
    return output, coefficients


def fit_and_score(
    values: Vector,
    observations: Sequence[Observation],
    design: Matrix,
) -> Dict[str, float]:
    train = [index for index, row in enumerate(observations) if row.split == "train"]
    heldout = [index for index, row in enumerate(observations) if row.split == "heldout"]
    train_x = [design[index] for index in train]
    heldout_x = [design[index] for index in heldout]
    weights = [1.0 / observations[index].se**2 for index in train]
    normal = [[math.fsum(weights[i] * train_x[i][j] * train_x[i][k]
                         for i in range(len(train)))
               for k in range(len(train_x[0]))]
              for j in range(len(train_x[0]))]
    rhs = [math.fsum(weights[i] * train_x[i][j] * values[train[i]]
                     for i in range(len(train)))
           for j in range(len(train_x[0]))]
    parameter_covariance = inverse(normal)
    parameters = matvec(parameter_covariance, rhs)
    predictions = matvec(heldout_x, parameters)
    residual = [values[index] - prediction for index, prediction in zip(heldout, predictions)]
    prediction_covariance = matmul(matmul(heldout_x, parameter_covariance), transpose(heldout_x))
    residual_covariance = [row[:] for row in prediction_covariance]
    for i, index in enumerate(heldout):
        residual_covariance[i][i] += observations[index].se**2
    precision = inverse(residual_covariance)
    chi_square = quadratic(residual, precision)
    det = determinant(residual_covariance)
    # The determinant term makes uncertain/ill-conditioned predictive models
    # pay for their broad forecast.  A small BIC-style parameter term breaks
    # residual ties in favor of the declared simpler model.
    predictive_deviance = chi_square + math.log(det) + len(parameters) * math.log(len(train))
    return {
        "heldout_chi_square": chi_square,
        "predictive_deviance": predictive_deviance,
        "parameter_count": float(len(parameters)),
    }


def select(scores: Dict[str, Dict[str, float]], criterion: str) -> str:
    return min(MODEL_ORDER, key=lambda model: (scores[model][criterion], MODEL_ORDER.index(model)))


def simulate(
    config: dict,
    observations: Sequence[Observation],
    replicates: int,
    seed: int,
    noise_scale: float,
    fraction: float,
) -> dict:
    reference_n = float(config["reference_N"])
    amplitude = float(config["reference_amplitude_h4"])
    truths, coefficients = truth_means(observations, reference_n, amplitude, fraction)
    designs = {model: basis_rows(model, observations, reference_n) for model in MODEL_ORDER}
    rng = random.Random(seed)
    criteria = ("heldout_chi_square", "predictive_deviance")
    confusion = {
        criterion: {truth: {candidate: 0 for candidate in MODEL_ORDER} for truth in MODEL_ORDER}
        for criterion in criteria
    }
    h4_acceptance = {truth: 0 for truth in MODEL_ORDER}
    for truth in MODEL_ORDER:
        mean = truths[truth]
        for _replica in range(replicates):
            values = [value + rng.gauss(0.0, row.se * noise_scale)
                      for value, row in zip(mean, observations)]
            # Score with the covariance belonging to the simulated precision.
            scaled_rows = [Observation(
                row.n, row.first, row.second, row.se * noise_scale, row.split,
                row.delta_cos4, row.delta_cos12,
            ) for row in observations]
            scores = {model: fit_and_score(values, scaled_rows, designs[model])
                      for model in MODEL_ORDER}
            for criterion in criteria:
                confusion[criterion][truth][select(scores, criterion)] += 1
            if scores["H4"]["heldout_chi_square"] <= 5.991464547107979:
                h4_acceptance[truth] += 1

    normalized = {
        criterion: {
            truth: {candidate: count / replicates for candidate, count in row.items()}
            for truth, row in matrix.items()
        }
        for criterion, matrix in confusion.items()
    }
    return {
        "replicates_per_truth": replicates,
        "seed": seed,
        "noise_scale_relative_to_current_se": noise_scale,
        "admixture_fraction_in_covariance_norm": fraction,
        "truth_coefficients": coefficients,
        "confusion_counts": confusion,
        "confusion_rates": normalized,
        "correct_selection_rate": {
            criterion: {truth: normalized[criterion][truth][truth] for truth in MODEL_ORDER}
            for criterion in criteria
        },
        "H4_nonrejection_rate_at_5pct": {
            truth: count / replicates for truth, count in h4_acceptance.items()
        },
    }


def markdown_report(config: dict, observations: Sequence[Observation], runs: Sequence[dict]) -> str:
    lines = [
        "# Synthetic model red-team at the current Gaussian-orientation precision",
        "",
        "This is a design-power calculation, not a fit to target data. It mirrors the frozen",
        "P32 split (`65,85,130` train; `145,170` held out), uses the pooled P31 standard",
        "errors, and preserves the exact Gaussian-pair H4/H12 angular columns.",
        "",
        "Two selectors are reported: raw held-out chi-square (closest to the current fixed-model",
        "score) and predictive deviance (chi-square plus forecast-volume and parameter penalties).",
        "A model can remain statistically acceptable without being selected.",
        "",
        "## Design",
        "",
        "| N | split | SE | delta cos4 | delta cos12 |",
        "|---:|:---|---:|---:|---:|",
    ]
    for row in observations:
        lines.append(f"| {row.n} | {row.split} | {row.se:.3g} | {row.delta_cos4:.6f} | {row.delta_cos12:.6f} |")
    for run in runs:
        lines.extend([
            "",
            "## Noise scale {:.3g}x current SE; admixture {:.3g}".format(
                run["noise_scale_relative_to_current_se"],
                run["admixture_fraction_in_covariance_norm"],
            ),
            "",
            "| hidden truth | correct: chi-square | correct: predictive deviance | H4 not rejected |",
            "|:---|---:|---:|---:|",
        ])
        for model in MODEL_ORDER:
            lines.append("| {} | {:.1%} | {:.1%} | {:.1%} |".format(
                model,
                run["correct_selection_rate"]["heldout_chi_square"][model],
                run["correct_selection_rate"]["predictive_deviance"][model],
                run["H4_nonrejection_rate_at_5pct"][model],
            ))
        lines.extend(["", "Predictive-deviance confusion matrix (row truth, column selected):", ""])
        lines.append("| truth \\ selected | " + " | ".join(MODEL_ORDER) + " |")
        lines.append("|:---|" + "---:|" * len(MODEL_ORDER))
        matrix = run["confusion_rates"]["predictive_deviance"]
        for truth in MODEL_ORDER:
            lines.append("| {} | {} |".format(
                truth, " | ".join("{:.1%}".format(matrix[truth][candidate]) for candidate in MODEL_ORDER)
            ))
    lines.extend([
        "",
        "## Interpretation guardrails",
        "",
        "- Pure H12 is amplitude-normalized to the H4 signal norm; mixed alternatives add a",
        "  component whose covariance-weighted norm is the declared admixture fraction of H4.",
        "- `two_radial_powers` means H4 angular structure with N^-13/8 and N^-9/8 terms.",
        "  `ordinary_correction` uses N^-13/8 and N^-21/8; `log_Jordan` uses",
        "  N^-13/8 times 1 and log(N/Nref).",
        "- Only two held-out points exist. Low correct-selection power or high H4 non-rejection",
        "  is therefore a property of the current design, not evidence that the mechanisms are equal.",
        "- Re-run after target covariance is available; do not tune admixture after seeing targets.",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--replicates", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=2026082801)
    parser.add_argument("--admixture", type=float)
    parser.add_argument("--noise-scales", type=float, nargs="+", default=[1.0, 0.5])
    parser.add_argument("--json", type=Path)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args()
    if args.replicates <= 0 or any(scale <= 0 for scale in args.noise_scales):
        raise SystemExit("replicates and noise scales must be positive")
    config, observations = load_config(args.config)
    fraction = (float(config["default_admixture_fraction"])
                if args.admixture is None else args.admixture)
    if fraction <= 0:
        raise SystemExit("admixture must be positive")
    runs = [simulate(config, observations, args.replicates, args.seed + index,
                     scale, fraction)
            for index, scale in enumerate(args.noise_scales)]
    payload = {
        "schema": "synthetic orientation model red-team result v1",
        "config": str(args.config),
        "provenance": config["provenance"],
        "models": list(MODEL_ORDER),
        "runs": runs,
    }
    text = json.dumps(payload, indent=2)
    print(text)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(markdown_report(config, observations, runs), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
