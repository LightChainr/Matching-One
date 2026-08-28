#!/usr/bin/env python3
"""Pre-register predictions for Jacobsen-style widths 22--24 using only n<=21.

The model ensemble and bias calibration are intentionally hard-coded. Changing them
creates a different prediction protocol and must occur before importing the target
values. Output is deterministic YAML and contains the SHA-256 of the input CSV.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import mpmath as mp


SELECTED_MODELS: tuple[tuple[int, tuple[int, ...]], ...] = (
    (8, (4, 6, 8, 10, 12, 14)),
    (9, (4, 6, 8, 10, 12)),
    (10, (4, 6, 8, 10)),
    (8, (4, 6, 8, 10, 12)),
)
CALIBRATION_CUTOFFS = (15, 16, 17, 18)
TARGET_CUTOFF = 21
TARGET_WIDTHS = (22, 23, 24)


@dataclass(frozen=True)
class Observation:
    n: int
    value: mp.mpf


def load(path: Path) -> list[Observation]:
    rows: list[Observation] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(Observation(int(row["n"]), mp.mpf(row["value"])))
    rows.sort(key=lambda row: row.n)
    if [row.n for row in rows] != list(range(1, 22)):
        raise ValueError("pre-registration input must contain exactly n=1..21")
    return rows


def fit(rows: Sequence[Observation], powers: Sequence[int]) -> tuple[mp.matrix, mp.mpf]:
    reference = mp.mpf(max(row.n for row in rows))
    design = mp.matrix(
        [
            [mp.mpf(1), *((reference / row.n) ** power for power in powers)]
            for row in rows
        ]
    )
    target = mp.matrix([row.value for row in rows])
    coefficients, _residual = mp.qr_solve(design, target)
    return coefficients, reference


def predict(
    n: int, coefficients: mp.matrix, reference: mp.mpf, powers: Sequence[int]
) -> mp.mpf:
    return coefficients[0] + mp.fsum(
        coefficients[index + 1] * (reference / n) ** power
        for index, power in enumerate(powers)
    )


def ensemble(rows: Sequence[Observation], cutoff: int, targets: Sequence[int]) -> list[list[mp.mpf]]:
    model_predictions: list[list[mp.mpf]] = []
    for n_min, powers in SELECTED_MODELS:
        training = [row for row in rows if n_min <= row.n <= cutoff]
        if len(training) < len(powers) + 1:
            raise ValueError(f"insufficient data for model n_min={n_min}, powers={powers}")
        coefficients, reference = fit(training, powers)
        model_predictions.append(
            [predict(n, coefficients, reference, powers) for n in targets]
        )
    return model_predictions


def mean(values: Sequence[mp.mpf]) -> mp.mpf:
    return mp.fsum(values) / len(values)


def linear_extrapolation(xs: Sequence[int], ys: Sequence[mp.mpf], target: int) -> mp.mpf:
    x_mean = mean([mp.mpf(x) for x in xs])
    y_mean = mean(list(ys))
    denominator = mp.fsum((mp.mpf(x) - x_mean) ** 2 for x in xs)
    slope = mp.fsum(
        (mp.mpf(x) - x_mean) * (y - y_mean) for x, y in zip(xs, ys, strict=True)
    ) / denominator
    return y_mean + slope * (mp.mpf(target) - x_mean)


def q(value: object) -> str:
    text = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{text}"'


def number(value: mp.mpf, digits: int = 35) -> str:
    return mp.nstr(value, n=digits, strip_zeros=False)


def render_yaml(path: Path, rows: Sequence[Observation]) -> str:
    actual = {row.n: row.value for row in rows}
    historical_errors: dict[int, list[mp.mpf]] = {1: [], 2: [], 3: []}

    for cutoff in CALIBRATION_CUTOFFS:
        targets = (cutoff + 1, cutoff + 2, cutoff + 3)
        predictions = ensemble(rows, cutoff, targets)
        for horizon, target in enumerate(targets, start=1):
            ensemble_mean = mean([model[horizon - 1] for model in predictions])
            historical_errors[horizon].append(ensemble_mean - actual[target])

    final_models = ensemble(rows, TARGET_CUTOFF, TARGET_WIDTHS)
    input_sha = hashlib.sha256(path.read_bytes()).hexdigest()

    lines = [
        "version: 1",
        "status: preregistered",
        "target_sequence: jacobsen_periodic_cylinder_site_percolation",
        f"input_csv: {q(str(path))}",
        f"input_sha256: {q(input_sha)}",
        "knowledge_barrier: >-",
        "  Generated using only the n=1..21 CSV. The n=22..24 values must not be",
        "  imported or used to alter this artifact before scoring.",
        "model_family: >-",
        "  Equal-weight ensemble of fixed correction expansions",
        "  p(n)=p_inf+sum_k a_k n^(-Delta_k).",
        "selected_models:",
    ]
    for n_min, powers in SELECTED_MODELS:
        lines.extend(
            [
                f"  - n_min: {n_min}",
                "    powers: [" + ", ".join(str(power) for power in powers) + "]",
            ]
        )

    lines.extend(
        [
            "bias_calibration:",
            "  cutoffs: [15, 16, 17, 18]",
            "  method: >-",
            "    For each forecast horizon h=1,2,3, fit a line to the signed ensemble",
            "    forecast errors versus training cutoff and extrapolate that bias to cutoff 21.",
            "interval_rule: >-",
            "  half_width = 2 * max historical absolute error at that horizon",
            "  + spread of the four final model predictions. This is a calibration band,",
            "  not a probabilistic confidence interval.",
            "predictions:",
        ]
    )

    for index, width in enumerate(TARGET_WIDTHS):
        horizon = index + 1
        model_values = [model[index] for model in final_models]
        raw_mean = mean(model_values)
        bias = linear_extrapolation(
            CALIBRATION_CUTOFFS, historical_errors[horizon], TARGET_CUTOFF
        )
        primary = raw_mean - bias
        spread = max(model_values) - min(model_values)
        half_width = 2 * max(abs(error) for error in historical_errors[horizon]) + spread
        lower = primary - half_width
        upper = primary + half_width
        lines.extend(
            [
                f"  - n: {width}",
                f"    primary: {q(number(primary))}",
                f"    lower: {q(number(lower))}",
                f"    upper: {q(number(upper))}",
                f"    raw_ensemble_mean: {q(number(raw_mean))}",
                f"    projected_bias_subtracted: {q(number(bias))}",
                f"    half_width: {q(number(half_width))}",
                "    model_predictions:",
            ]
        )
        for (n_min, powers), value in zip(SELECTED_MODELS, model_values, strict=True):
            lines.extend(
                [
                    f"      - model: {q(f'nmin={n_min};powers={','.join(map(str, powers))}')}",
                    f"        value: {q(number(value))}",
                ]
            )

    lines.extend(
        [
            "scoring_after_reveal:",
            "  per_width: absolute_error",
            "  aggregate: root_mean_square_error",
            "  rule: >-",
            "    Score the primary values exactly as committed. Also report raw-ensemble",
            "    and every component-model error. Do not retune before scoring.",
            "warning: >-",
            "  These are predictions of finite-width p_c(n), not an estimate of the",
            "  infinite-lattice threshold.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "csv",
        type=Path,
        nargs="?",
        default=Path("data/jacobsen_2015_square_site_cylinder.csv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("predictions/polynomial_widths_22_24.yaml"),
    )
    parser.add_argument("--dps", type=int, default=100)
    args = parser.parse_args()
    if args.dps < 60:
        raise SystemExit("use at least 60 decimal digits")
    mp.mp.dps = args.dps
    rows = load(args.csv)
    output = render_yaml(args.csv, rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output, encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
