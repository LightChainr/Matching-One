#!/usr/bin/env python3
"""Score the frozen N365 third-readout annulus acquisition for Issue 253."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Sequence

try:
    from scripts.analyze_annulus_channel_recurrence import (
        chi_square_survival,
        finite_jacobian,
        inverse,
        propagated_covariance,
        quadratic,
        recurrence,
        submatrix,
        subvector,
    )
except ModuleNotFoundError:  # Direct `python scripts/...` execution.
    from analyze_annulus_channel_recurrence import (
        chi_square_survival,
        finite_jacobian,
        inverse,
        propagated_covariance,
        quadratic,
        recurrence,
        submatrix,
        subvector,
    )


CHANNELS = ("A_plus", "A_minus")
CSV_COLUMNS = {"A_plus": "h4_plus", "A_minus": "h4_minus"}
DYADIC_RADII = (2, 4, 8)
ALL_RADII = (2, 4, 7, 8)
OLD_GEOMETRIES = (325, 425)


def covariance_of_mean(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    means = [sum(row[j] for row in rows) / len(rows) for j in range(len(rows[0]))]
    return [
        [
            sum((row[i] - means[i]) * (row[j] - means[j]) for row in rows)
            / (len(rows) * (len(rows) - 1))
            for j in range(len(means))
        ]
        for i in range(len(means))
    ]


def block_diagonal(first: Sequence[Sequence[float]], second: Sequence[Sequence[float]]) -> list[list[float]]:
    n, m = len(first), len(second)
    result = [[0.0] * (n + m) for _ in range(n + m)]
    for i in range(n):
        for j in range(n):
            result[i][j] = float(first[i][j])
    for i in range(m):
        for j in range(m):
            result[n + i][n + j] = float(second[i][j])
    return result


def read_n365(metadata_path: Path) -> dict:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata["radii"] != list(ALL_RADII):
        raise ValueError("N365 radii differ from the frozen acquisition")
    designs = [row["label"] for row in metadata["designs"]]
    if designs != ["n365_first", "n365_second"]:
        raise ValueError("N365 design order differs from the frozen acquisition")
    batch_path = metadata_path.parent / Path(metadata["batch_csv"]).name
    with batch_path.open(newline="", encoding="utf-8") as handle:
        raw = list(csv.DictReader(handle))
    indexed = {(row["label"], int(row["radius"]), int(row["batch"])): row for row in raw}
    batch_ids = sorted({int(row["batch"]) for row in raw})
    order = [f"N365_R{radius}_Delta_{channel}" for channel in CHANNELS for radius in ALL_RADII]
    batch_vectors = []
    for batch in batch_ids:
        vector = []
        for channel in CHANNELS:
            column = CSV_COLUMNS[channel]
            for radius in ALL_RADII:
                first = indexed[("n365_first", radius, batch)]
                second = indexed[("n365_second", radius, batch)]
                if first["samples"] != second["samples"]:
                    raise ValueError("paired N365 design rows have different sample counts")
                vector.append((int(first[column]) - int(second[column])) / int(first["samples"]))
        batch_vectors.append(vector)
    point = [sum(row[j] for row in batch_vectors) / len(batch_vectors) for j in range(len(order))]
    return {
        "order": order,
        "point": point,
        "covariance": covariance_of_mean(batch_vectors),
        "batches": len(batch_vectors),
        "metadata": metadata,
        "batch_csv": str(batch_path),
    }


def old_labels(channel: str) -> list[str]:
    return [
        f"N{geometry}_R{radius}_Delta_{channel}"
        for geometry in OLD_GEOMETRIES
        for radius in DYADIC_RADII
    ]


def new_labels(channel: str, radii: Sequence[int]) -> list[str]:
    return [f"N365_R{radius}_Delta_{channel}" for radius in radii]


def dyadic_residual(values: Sequence[float]) -> list[float]:
    rec = recurrence(values[:6])
    y0, y1, y2 = values[6:9]
    return [y2 - rec["T"] * y1 + rec["D"] * y0]


def fractional_prediction(trace: float, determinant: float, y0: float, y1: float, radial_class: str) -> float:
    step = math.log(7.0 / 2.0, 2.0)
    discriminant = trace * trace - 4.0 * determinant
    if radial_class == "R2":
        root = math.sqrt(discriminant)
        first, second = (trace + root) / 2.0, (trace - root) / 2.0
        if min(first, second) <= 0.0:
            raise ValueError("R2 roots must be positive")
        a = (y1 - second * y0) / (first - second)
        b = (first * y0 - y1) / (first - second)
        return a * first**step + b * second**step
    if radial_class == "C2":
        modulus = math.sqrt(determinant)
        theta = math.acos(trace / (2.0 * modulus))
        sine_amplitude = (y1 / modulus - y0 * math.cos(theta)) / math.sin(theta)
        return modulus**step * (y0 * math.cos(theta * step) + sine_amplitude * math.sin(theta * step))
    raise ValueError(radial_class)


def radius7_residual(values: Sequence[float], radial_class: str) -> list[float]:
    rec = recurrence(values[:6])
    y0, y1, y7 = values[6:9]
    return [y7 - fractional_prediction(rec["T"], rec["D"], y0, y1, radial_class)]


def scalar_score(function, point: Sequence[float], covariance: Sequence[Sequence[float]]) -> dict:
    residual = list(function(point))
    jacobian = finite_jacobian(function, point, covariance)
    residual_covariance = propagated_covariance(jacobian, covariance)
    chi_square = quadratic(residual, inverse(residual_covariance))
    return {
        "residual": residual[0],
        "standard_error_delta_method": math.sqrt(residual_covariance[0][0]),
        "z": residual[0] / math.sqrt(residual_covariance[0][0]),
        "chi_square": chi_square,
        "degrees_of_freedom": 1,
        "chi_square_survival": chi_square_survival(chi_square, 1),
    }


def channel_score(old: dict, new: dict, channel: str) -> dict:
    old_source_labels = old_labels(channel)
    old_point = subvector(old["order"], old["point"], old_source_labels)
    old_covariance = submatrix(old["order"], old["covariance"], old_source_labels)
    rec = recurrence(old_point)
    dyadic_labels = new_labels(channel, DYADIC_RADII)
    dyadic_point = subvector(new["order"], new["point"], dyadic_labels)
    dyadic_covariance = submatrix(new["order"], new["covariance"], dyadic_labels)
    dyadic_joined = old_point + dyadic_point
    dyadic_score = scalar_score(
        dyadic_residual, dyadic_joined, block_diagonal(old_covariance, dyadic_covariance)
    )
    if rec["Delta"] > 0.0:
        radial_class = "R2"
    elif rec["D"] > 0.0:
        radial_class = "C2"
    else:
        radial_class = None
    radius7_score = None
    if radial_class is not None:
        target_labels = new_labels(channel, (2, 4, 7))
        target_point = subvector(new["order"], new["point"], target_labels)
        target_covariance = submatrix(new["order"], new["covariance"], target_labels)
        joined = old_point + target_point
        covariance = block_diagonal(old_covariance, target_covariance)
        radius7_score = scalar_score(
            lambda values: radius7_residual(values, radial_class), joined, covariance
        )
        radius7_score["branch"] = radial_class
    return {
        "channel": channel,
        "old_saturated_recurrence": rec,
        "N365_dyadic_point_R2_R4_R8": dyadic_point,
        "heldout_dyadic_recurrence": dyadic_score,
        "heldout_R7_fractional_propagation": radius7_score,
    }


def joint_dyadic_score(old: dict, new: dict) -> dict:
    old_order = [label for channel in CHANNELS for label in old_labels(channel)]
    new_order = [label for channel in CHANNELS for label in new_labels(channel, DYADIC_RADII)]
    old_point = subvector(old["order"], old["point"], old_order)
    new_point = subvector(new["order"], new["point"], new_order)
    covariance = block_diagonal(
        submatrix(old["order"], old["covariance"], old_order),
        submatrix(new["order"], new["covariance"], new_order),
    )
    point = old_point + new_point

    def function(values: Sequence[float]) -> list[float]:
        plus_old, minus_old = values[:6], values[6:12]
        plus_new, minus_new = values[12:15], values[15:18]
        return [
            dyadic_residual(plus_old + plus_new)[0],
            dyadic_residual(minus_old + minus_new)[0],
        ]

    residual = function(point)
    jacobian = finite_jacobian(function, point, covariance)
    residual_covariance = propagated_covariance(jacobian, covariance)
    chi_square = quadratic(residual, inverse(residual_covariance))
    return {
        "order": list(CHANNELS),
        "residual": residual,
        "covariance_delta_method": residual_covariance,
        "chi_square": chi_square,
        "degrees_of_freedom": 2,
        "chi_square_survival": chi_square_survival(chi_square, 2),
        "assumption": "separate T,D inherited from the old two-readout saturation; only the third geometry is held out",
    }


def render(old_path: Path, metadata_path: Path) -> dict:
    old_payload = json.loads(old_path.read_text(encoding="utf-8"))["contrast_vector"]
    new = read_n365(metadata_path)
    return {
        "schema": "matching-one.p253-n365-heldout.v1",
        "issue": 253,
        "status": "scored_frozen_third_geometry_acquisition",
        "old_source": str(old_path),
        "new_source": str(metadata_path),
        "N365_contrast_block": new,
        "channels": [channel_score(old_payload, new, channel) for channel in CHANNELS],
        "joint_heldout_dyadic_recurrence": joint_dyadic_score(old_payload, new),
        "scope": [
            "N325/N425 determine T,D before N365 enters the score.",
            "N365 is an independent replica-counter block, so old/new covariance is block diagonal.",
            "The R7 check uses the same N365 block and is correlated design information, not another evidence row.",
            "Passing or failing this recurrence tests one shared two-state radial generator per channel; it does not by itself identify Jordan, two-real, or complex scaling dimensions.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old-analysis", required=True, type=Path)
    parser.add_argument("--n365-metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = render(args.old_analysis, args.n365_metadata)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
