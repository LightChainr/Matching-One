#!/usr/bin/env python3
"""Pilot-frozen minimum-variance combination of equal-mean channels.

Rows are simultaneous measurements of estimators that have the same expected
value (for example the torus matching differences ``c``, ``b``, ``e``, and a
directional channel).  Learn covariance weights on an independent pilot set,
freeze them, and apply them to an evaluation set to avoid adaptive bias.

The implementation uses only the Python standard library so it can serve as a
small CPU/server reference.  A diagonal ridge is introduced only when asked
for or when the sample covariance is singular at machine precision.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


Matrix = list[list[float]]


def _validate_rows(rows: Sequence[Sequence[float]]) -> int:
    if len(rows) < 2:
        raise ValueError("at least two rows are required")
    width = len(rows[0])
    if width == 0 or any(len(row) != width for row in rows):
        raise ValueError("rows must form a nonempty rectangular matrix")
    if any(not math.isfinite(float(value)) for row in rows for value in row):
        raise ValueError("all channel values must be finite")
    return width


def sample_covariance(rows: Sequence[Sequence[float]]) -> Matrix:
    """Unbiased sample covariance of row-wise simultaneous observations."""

    width = _validate_rows(rows)
    means = [
        math.fsum(float(row[j]) for row in rows) / len(rows) for j in range(width)
    ]
    covariance = [[0.0] * width for _ in range(width)]
    for i in range(width):
        for j in range(i, width):
            value = math.fsum(
                (float(row[i]) - means[i]) * (float(row[j]) - means[j])
                for row in rows
            ) / (len(rows) - 1)
            covariance[i][j] = covariance[j][i] = value
    return covariance


def _solve(matrix: Matrix, vector: Sequence[float]) -> list[float]:
    """Solve a dense system by scaled partial-pivot Gaussian elimination."""

    n = len(vector)
    augmented = [list(map(float, matrix[i])) + [float(vector[i])] for i in range(n)]
    scale = max((abs(value) for row in matrix for value in row), default=0.0)
    tolerance = max(scale * 1e-14, 1e-300)
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= tolerance:
            raise ArithmeticError("singular covariance matrix")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        pivot_value = augmented[column][column]
        for entry in range(column, n + 1):
            augmented[column][entry] /= pivot_value
        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor == 0.0:
                continue
            for entry in range(column, n + 1):
                augmented[row][entry] -= factor * augmented[column][entry]
    return [augmented[i][-1] for i in range(n)]


def minimum_variance_weights(
    covariance: Sequence[Sequence[float]], ridge: float = 0.0
) -> tuple[list[float], float]:
    """Return sum-one GLS weights and the diagonal ridge actually applied."""

    n = len(covariance)
    if n == 0 or any(len(row) != n for row in covariance):
        raise ValueError("covariance must be a nonempty square matrix")
    if ridge < 0 or not math.isfinite(ridge):
        raise ValueError("ridge must be a finite nonnegative number")
    matrix = [list(map(float, row)) for row in covariance]
    if any(not math.isfinite(value) for row in matrix for value in row):
        raise ValueError("covariance entries must be finite")
    for i in range(n):
        for j in range(i):
            if not math.isclose(
                matrix[i][j], matrix[j][i], rel_tol=1e-10, abs_tol=1e-300
            ):
                raise ValueError("covariance must be symmetric")

    scale = math.fsum(abs(matrix[i][i]) for i in range(n)) / n
    if scale == 0.0:
        scale = max((abs(value) for row in matrix for value in row), default=0.0)
    if scale == 0.0:
        scale = 1.0
    applied = ridge
    for attempt in range(9):
        trial = [row[:] for row in matrix]
        for i in range(n):
            trial[i][i] += applied
        try:
            inverse_ones = _solve(trial, [1.0] * n)
            normalizer = math.fsum(inverse_ones)
            if not math.isfinite(normalizer) or normalizer <= 0.0:
                raise ArithmeticError("unstable GLS normalization")
            return [value / normalizer for value in inverse_ones], applied
        except ArithmeticError:
            if ridge > 0:
                raise
            applied = scale * 1e-12 * (10**attempt)
    raise ArithmeticError("could not regularize covariance matrix")


def combine_rows(
    rows: Sequence[Sequence[float]], weights: Sequence[float]
) -> list[float]:
    width = _validate_rows(rows)
    if len(weights) != width:
        raise ValueError("weight count does not match channel count")
    if not math.isclose(math.fsum(weights), 1.0, rel_tol=1e-10, abs_tol=1e-10):
        raise ValueError("weights must sum to one")
    return [
        math.fsum(float(value) * weight for value, weight in zip(row, weights))
        for row in rows
    ]


@dataclass(frozen=True)
class FrozenEstimator:
    channel_names: tuple[str, ...]
    weights: tuple[float, ...]
    pilot_covariance: tuple[tuple[float, ...], ...]
    applied_ridge: float

    @classmethod
    def fit(
        cls,
        channel_names: Sequence[str],
        pilot_rows: Sequence[Sequence[float]],
        ridge: float = 0.0,
    ) -> "FrozenEstimator":
        if len(set(channel_names)) != len(channel_names):
            raise ValueError("channel names must be unique")
        width = _validate_rows(pilot_rows)
        if len(channel_names) != width:
            raise ValueError("channel name count does not match pilot data")
        covariance = sample_covariance(pilot_rows)
        weights, applied_ridge = minimum_variance_weights(covariance, ridge)
        return cls(
            tuple(channel_names),
            tuple(weights),
            tuple(tuple(row) for row in covariance),
            applied_ridge,
        )

    def evaluate(self, rows: Sequence[Sequence[float]]) -> dict[str, object]:
        values = combine_rows(rows, self.weights)
        mean = math.fsum(values) / len(values)
        variance = math.fsum((value - mean) ** 2 for value in values) / (len(values) - 1)
        return {
            "samples": len(values),
            "channel_names": list(self.channel_names),
            "weights": list(self.weights),
            "applied_ridge": self.applied_ridge,
            "mean": mean,
            "sample_variance": variance,
            "standard_error": math.sqrt(variance / len(values)),
        }


@dataclass(frozen=True)
class FrozenZeroMeanControls:
    """Pilot-frozen coefficients for a target plus exact zero-mean controls.

    Unlike :class:`FrozenEstimator`, the columns do not have equal means and
    their coefficients need not sum to one.  If ``Z`` has known expectation
    zero, ``target + beta^T Z`` remains unbiased for every fixed ``beta``.
    """

    control_names: tuple[str, ...]
    coefficients: tuple[float, ...]
    pilot_covariance: tuple[tuple[float, ...], ...]
    applied_ridge: float

    @classmethod
    def fit(
        cls,
        control_names: Sequence[str],
        pilot_targets: Sequence[float],
        pilot_controls: Sequence[Sequence[float]],
        ridge: float = 0.0,
    ) -> "FrozenZeroMeanControls":
        if len(set(control_names)) != len(control_names) or not control_names:
            raise ValueError("control names must be nonempty and unique")
        width = _validate_rows(pilot_controls)
        if len(control_names) != width:
            raise ValueError("control name count does not match pilot data")
        if len(pilot_targets) != len(pilot_controls):
            raise ValueError("pilot target/control row counts differ")
        if any(not math.isfinite(float(value)) for value in pilot_targets):
            raise ValueError("pilot targets must be finite")
        if ridge < 0 or not math.isfinite(ridge):
            raise ValueError("ridge must be a finite nonnegative number")

        joint_rows = [
            [float(target)] + [float(value) for value in controls]
            for target, controls in zip(pilot_targets, pilot_controls)
        ]
        covariance = sample_covariance(joint_rows)
        control_covariance = [row[1:] for row in covariance[1:]]
        target_covariance = [covariance[index][0] for index in range(1, width + 1)]
        scale = math.fsum(
            abs(control_covariance[index][index]) for index in range(width)
        ) / width
        if scale == 0.0:
            scale = 1.0
        applied = ridge
        for attempt in range(9):
            trial = [row[:] for row in control_covariance]
            for index in range(width):
                trial[index][index] += applied
            try:
                solution = _solve(trial, target_covariance)
                coefficients = [-value for value in solution]
                if any(not math.isfinite(value) for value in coefficients):
                    raise ArithmeticError("nonfinite control coefficient")
                return cls(
                    tuple(control_names),
                    tuple(coefficients),
                    tuple(tuple(row) for row in covariance),
                    applied,
                )
            except ArithmeticError:
                if ridge > 0:
                    raise
                applied = scale * 1e-12 * (10**attempt)
        raise ArithmeticError("could not regularize control covariance matrix")

    def adjusted_values(
        self,
        targets: Sequence[float],
        controls: Sequence[Sequence[float]],
    ) -> list[float]:
        width = _validate_rows(controls)
        if width != len(self.coefficients):
            raise ValueError("control count does not match frozen coefficients")
        if len(targets) != len(controls):
            raise ValueError("target/control row counts differ")
        return [
            float(target)
            + math.fsum(
                coefficient * float(value)
                for coefficient, value in zip(self.coefficients, row)
            )
            for target, row in zip(targets, controls)
        ]

    def evaluate(
        self,
        targets: Sequence[float],
        controls: Sequence[Sequence[float]],
    ) -> dict[str, object]:
        adjusted = self.adjusted_values(targets, controls)
        mean = math.fsum(adjusted) / len(adjusted)
        variance = math.fsum((value - mean) ** 2 for value in adjusted) / (
            len(adjusted) - 1
        )
        return {
            "samples": len(adjusted),
            "control_names": list(self.control_names),
            "coefficients": list(self.coefficients),
            "applied_ridge": self.applied_ridge,
            "mean": mean,
            "sample_variance": variance,
            "standard_error": math.sqrt(variance / len(adjusted)),
        }


def _read_csv(path: Path, channels: Sequence[str]) -> list[list[float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = [channel for channel in channels if channel not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        return [[float(row[channel]) for channel in channels] for row in reader]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--evaluation", type=Path, required=True)
    parser.add_argument("--channels", nargs="+", required=True)
    parser.add_argument("--ridge", type=float, default=0.0)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    try:
        pilot = _read_csv(args.pilot, args.channels)
        evaluation = _read_csv(args.evaluation, args.channels)
        estimator = FrozenEstimator.fit(args.channels, pilot, args.ridge)
        result = estimator.evaluate(evaluation)
        result["pilot_samples"] = len(pilot)
        result["pilot_covariance"] = [list(row) for row in estimator.pilot_covariance]
    except (ValueError, ArithmeticError) as exc:
        raise SystemExit(str(exc)) from exc
    output = json.dumps(result, indent=2)
    print(output)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(output + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
