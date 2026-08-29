#!/usr/bin/env python3
"""Score the frozen #263 boundary Q-score pilot against the ODE tangent."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

import mpmath as mp

from p263_boundary_lattice_qscore import ANCHOR_INDEX, LINK_PATTERNS, render as frozen_render


GEOMETRY_ORDER = ("lambda_1_4", "lambda_1_3", "lambda_2_3", "lambda_3_4")
HIGH_CHANNEL = "14_23"


def read_rows(paths: Iterable[Path]) -> list[dict[str, int | str]]:
    rows: list[dict[str, int | str]] = []
    integer_fields = {
        "lambda_num", "lambda_den", "level", "span_L", "nx", "ny", "vertices",
        "edges", "batch", "seed", "sample_begin", "samples", "sum_J", "sum_J2", "count_1234",
        "sum_J_1234", "sum_J2_1234", "count_12_34", "sum_J_12_34",
        "sum_J2_12_34", "count_14_23", "sum_J_14_23", "sum_J2_14_23",
    }
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            for raw in csv.DictReader(handle):
                row: dict[str, int | str] = {}
                for key, value in raw.items():
                    row[key] = int(value) if key in integer_fields else value
                rows.append(row)
    return rows


def _aggregate(rows: list[dict[str, int | str]], omitted_batch: int | None) -> dict:
    frozen = frozen_render()
    target = frozen["amplitude_gauge_and_score"]["frozen_target"]
    expected_lambda = target["lambda_order"]
    by_geometry: dict[str, list[dict[str, int | str]]] = {key: [] for key in GEOMETRY_ORDER}
    for row in rows:
        if omitted_batch is not None and row["batch"] == omitted_batch:
            continue
        by_geometry[str(row["geometry_id"])].append(row)
    estimates = []
    h_prime = math.sqrt(3) / (3 * math.pi)
    for index, geometry_id in enumerate(GEOMETRY_ORDER):
        selected = by_geometry[geometry_id]
        if not selected:
            raise ValueError(f"missing rows for {geometry_id}")
        lambda_text = f"{selected[0]['lambda_num']}/{selected[0]['lambda_den']}"
        if lambda_text != expected_lambda[index]:
            raise ValueError(f"lambda mismatch for {geometry_id}")
        invariant = (selected[0]["level"], selected[0]["span_L"], selected[0]["vertices"], selected[0]["edges"])
        if any(
            (row["level"], row["span_L"], row["vertices"], row["edges"]) != invariant
            for row in selected
        ):
            raise ValueError("geometry metadata changed across batches")
        samples = sum(int(row["samples"]) for row in selected)
        sum_j = sum(int(row["sum_J"]) for row in selected)
        count = sum(int(row[f"count_{HIGH_CHANNEL}"]) for row in selected)
        sum_j_channel = sum(int(row[f"sum_J_{HIGH_CHANNEL}"]) for row in selected)
        if count == 0:
            raise ValueError(f"zero high-channel events for {geometry_id}")
        probability = count / samples
        measure = sum_j_channel / (2 * samples) - probability * sum_j / (2 * samples)
        d_log_probability = measure / probability
        span = int(selected[0]["span_L"])
        k_prefactor = float(
            frozen["frozen_geometries"][index]["K_decimal"]
        )
        z = d_log_probability + 4 * h_prime * math.log(span) - 2 * h_prime * math.log(k_prefactor)
        estimates.append(
            {
                "geometry_id": geometry_id,
                "lambda": lambda_text,
                "span_L": span,
                "samples": samples,
                "events_14_23": count,
                "probability_14_23": probability,
                "measure_tangent_14_23": measure,
                "d_log_probability": d_log_probability,
                "z_before_amplitude_projection": z,
            }
        )
    anchored = [row["z_before_amplitude_projection"] - estimates[ANCHOR_INDEX]["z_before_amplitude_projection"] for row in estimates]
    residual_full = [value - target_value for value, target_value in zip(anchored, target["anchored_dQ_logU"])]
    active = [index for index in range(4) if index != ANCHOR_INDEX]
    return {
        "estimates": estimates,
        "anchored_lattice_tangent": anchored,
        "frozen_target": target["anchored_dQ_logU"],
        "residual": [residual_full[index] for index in active],
    }


def _jackknife_covariance(full: list[float], deleted: list[list[float]]) -> list[list[float]]:
    batches = len(deleted)
    pseudo = [
        [batches * full[i] - (batches - 1) * row[i] for i in range(len(full))]
        for row in deleted
    ]
    means = [sum(row[i] for row in pseudo) / batches for i in range(len(full))]
    return [
        [
            sum((row[i] - means[i]) * (row[j] - means[j]) for row in pseudo)
            / (batches * (batches - 1))
            for j in range(len(full))
        ]
        for i in range(len(full))
    ]


def _gls(residual: list[float], covariance: list[list[float]]) -> dict:
    matrix = mp.matrix(covariance)
    values, vectors = mp.eigsy(matrix)
    maximum = max(abs(values[index]) for index in range(len(residual)))
    cutoff = maximum * mp.mpf("1e-10")
    inverse = mp.zeros(len(residual))
    rank = 0
    for index in range(len(residual)):
        if values[index] > cutoff:
            column = vectors[:, index]
            inverse += (column * column.T) / values[index]
            rank += 1
    vector = mp.matrix(residual)
    chi_square = (vector.T * inverse * vector)[0] if rank else mp.nan
    return {
        "chi_square": float(chi_square),
        "degrees_of_freedom": rank,
        "covariance_eigenvalues": [float(values[index]) for index in range(len(residual))],
        "relative_pseudoinverse_cutoff": 1e-10,
    }


def score(rows: list[dict[str, int | str]]) -> dict:
    batch_sets = {
        geometry: {int(row["batch"]) for row in rows if row["geometry_id"] == geometry}
        for geometry in GEOMETRY_ORDER
    }
    if not batch_sets or any(value != next(iter(batch_sets.values())) for value in batch_sets.values()):
        raise ValueError("all geometries must have identical synchronized batch ids")
    batch_ids = sorted(next(iter(batch_sets.values())))
    if len(batch_ids) < 3:
        raise ValueError("at least three batches are required")
    full = _aggregate(rows, None)
    deleted = [_aggregate(rows, batch)["residual"] for batch in batch_ids]
    covariance = _jackknife_covariance(full["residual"], deleted)
    return {
        "schema": "matching-one.p263-boundary-qscore-pilot-score.v1",
        "issue": 263,
        "status": "pilot_not_continuum_evidence",
        "batch_count": len(batch_ids),
        "state_order": ["lambda_1_4", "lambda_2_3", "lambda_3_4"],
        **full,
        "residual_covariance": covariance,
        "joint_gls": _gls(full["residual"], covariance),
        "warnings": [
            "The pilot rectangles are finite boxes, not the upper half-plane.",
            "A single resolution level cannot separate scaling correction from ODE-shape error.",
            "The GLS result is a software/sensitivity diagnostic until both frozen levels are acquired.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = score(read_rows(args.batches))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
