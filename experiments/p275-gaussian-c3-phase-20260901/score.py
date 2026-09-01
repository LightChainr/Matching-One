#!/usr/bin/env python3
"""GLS signed-real-gain scorer for the P275 paired Gaussian C3 phase gate."""

from __future__ import annotations

import argparse
import csv
import json
from math import atan2, cos, erfc, inf, pi, sin, sqrt, tan
from pathlib import Path
from typing import Callable, Sequence


ALPHA = 0.01
DELTA = atan2(5.0, 12.0)
COORDINATES = ("z1_re", "z1_im", "z2_re", "z2_im")


def covariance_of_mean(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    count = len(rows)
    if count <= 4:
        raise ValueError("at least five paired batches are required")
    means = [sum(row[j] for row in rows) / count for j in range(4)]
    return [
        [
            sum((row[i] - means[i]) * (row[j] - means[j]) for row in rows)
            / (count * (count - 1))
            for j in range(4)
        ]
        for i in range(4)
    ]


def invert(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    size = len(matrix)
    work = [
        [float(value) for value in row]
        + [1.0 if index == column else 0.0 for column in range(size)]
        for index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(work[row][column]))
        if abs(work[pivot][column]) < 1e-30:
            raise ValueError("paired covariance is singular")
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [value / scale for value in work[column]]
        for row in range(size):
            if row == column:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][j] - factor * work[column][j] for j in range(2 * size)
            ]
    return [row[size:] for row in work]


def mat_vec(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> list[float]:
    return [sum(a * b for a, b in zip(row, vector)) for row in matrix]


def dot(first: Sequence[float], second: Sequence[float]) -> float:
    return sum(a * b for a, b in zip(first, second))


def rotation(angle: float) -> list[list[float]]:
    return [[cos(angle), -sin(angle)], [sin(angle), cos(angle)]]


def profile_at_t(
    y: Sequence[float], precision: Sequence[Sequence[float]], rot: Sequence[Sequence[float]], t: float
) -> tuple[float, list[float], list[float]]:
    c, s = cos(t), sin(t)
    design = [
        [c, 0.0],
        [0.0, c],
        [s * rot[0][0], s * rot[0][1]],
        [s * rot[1][0], s * rot[1][1]],
    ]
    wx = [[sum(precision[i][k] * design[k][j] for k in range(4)) for j in range(2)] for i in range(4)]
    normal = [[sum(design[i][a] * wx[i][b] for i in range(4)) for b in range(2)] for a in range(2)]
    wy = mat_vec(precision, y)
    rhs = [sum(design[i][a] * wy[i] for i in range(4)) for a in range(2)]
    determinant = normal[0][0] * normal[1][1] - normal[0][1] * normal[1][0]
    if abs(determinant) < 1e-30:
        return inf, [0.0, 0.0], [0.0] * 4
    amplitude = [
        (rhs[0] * normal[1][1] - rhs[1] * normal[0][1]) / determinant,
        (normal[0][0] * rhs[1] - normal[1][0] * rhs[0]) / determinant,
    ]
    fitted = [sum(design[i][j] * amplitude[j] for j in range(2)) for i in range(4)]
    residual = [y[i] - fitted[i] for i in range(4)]
    chi2 = dot(residual, mat_vec(precision, residual))
    return max(0.0, chi2), amplitude, fitted


def golden_minimize(function: Callable[[float], float], left: float, right: float) -> float:
    ratio = (sqrt(5.0) - 1.0) / 2.0
    x1 = right - ratio * (right - left)
    x2 = left + ratio * (right - left)
    f1, f2 = function(x1), function(x2)
    for _ in range(100):
        if right - left < 1e-13:
            break
        if f1 <= f2:
            right, x2, f2 = x2, x1, f1
            x1 = right - ratio * (right - left)
            f1 = function(x1)
        else:
            left, x1, f1 = x1, x2, f2
            x2 = left + ratio * (right - left)
            f2 = function(x2)
    return (left + right) / 2.0


def profile_model(
    y: Sequence[float], covariance: Sequence[Sequence[float]], phase: float
) -> dict[str, object]:
    precision = invert(covariance)
    rot = rotation(phase)
    grid_count = 4096
    step = pi / grid_count
    grid = [-pi / 2 + index * step for index in range(grid_count + 1)]
    values = [profile_at_t(y, precision, rot, t)[0] for t in grid]
    best_index = min(range(len(grid)), key=values.__getitem__)
    left = grid[max(0, best_index - 1)]
    right = grid[min(grid_count, best_index + 1)]
    objective = lambda t: profile_at_t(y, precision, rot, t)[0]
    optimum = golden_minimize(objective, left, right)
    chi2, amplitude, fitted = profile_at_t(y, precision, rot, optimum)
    cosine = cos(optimum)
    gain = tan(optimum) if abs(cosine) > 1e-10 else (inf if sin(optimum) >= 0 else -inf)
    return {
        "phase_radians": phase,
        "profile_t_radians": optimum,
        "signed_real_gain": gain,
        "base_amplitude_re_im": amplitude,
        "fitted_coordinates": fitted,
        "chi_square": chi2,
        "reference_df": 1,
        "p_value": erfc(sqrt(chi2 / 2.0)),
    }


def score(rows: Sequence[Sequence[float]]) -> dict[str, object]:
    means = [sum(row[j] for row in rows) / len(rows) for j in range(4)]
    covariance = covariance_of_mean(rows)
    models = {
        "H4": profile_model(means, covariance, +4.0 * DELTA),
        "H8": profile_model(means, covariance, -8.0 * DELTA),
    }
    passes = {name: model["p_value"] >= ALPHA for name, model in models.items()}
    if passes == {"H4": True, "H8": False}:
        decision = "H4_SELECTED_H8_STOP"
    elif passes == {"H4": False, "H8": True}:
        decision = "H8_SELECTED_H4_STOP"
    elif passes == {"H4": True, "H8": True}:
        decision = "UNRESOLVED_BOTH_SURVIVE_STOP"
    else:
        decision = "SIGNED_REAL_TRANSPORT_FAILED_STOP"
    return {
        "schema": "matching-one/p275-gaussian-c3-phase-score/v1",
        "batch_count": len(rows),
        "coordinate_order": list(COORDINATES),
        "mean": means,
        "covariance_of_mean": covariance,
        "delta_radians": DELTA,
        "delta_degrees": DELTA * 180.0 / pi,
        "alpha": ALPHA,
        "models": models,
        "decision": decision,
        "top_up_authorized": False,
    }


def read_rows(path: Path) -> list[list[float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or any(name not in reader.fieldnames for name in COORDINATES):
            raise ValueError("batch CSV lacks the four frozen C3 coordinates")
        return [[float(row[name]) for name in COORDINATES] for row in reader]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("batches", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    payload = score(read_rows(args.batches))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"decision": payload["decision"], "output": str(args.out)}))


if __name__ == "__main__":
    main()
