#!/usr/bin/env python3
"""Score the exploratory stopped-transcript #263 Q-tangent estimator."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

from p263_boundary_lattice_qscore import ANCHOR_INDEX, render as frozen_render
from score_p263_boundary_qscore_pilot import _gls, _jackknife_covariance


GEOMETRY_ORDER = ("lambda_1_4", "lambda_1_3", "lambda_2_3", "lambda_3_4")
INTEGER_FIELDS = {
    "lambda_num", "lambda_den", "level", "span_L", "nx", "ny", "vertices",
    "edges", "batch", "outer_seed", "completion_seed", "inner_replicates",
    "sample_begin", "samples", "count_14_23", "sum_delta_J_14_23",
    "sum_delta_J_inner_square_14_23", "sum_delta_J2_individual_14_23",
    "sum_revealed", "sum_revealed2", "sum_revealed_14_23", "max_revealed",
}


def read_rows(paths: Iterable[Path]) -> list[dict[str, int | str]]:
    rows: list[dict[str, int | str]] = []
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            for raw in csv.DictReader(handle):
                missing = INTEGER_FIELDS.difference(raw)
                if missing:
                    raise ValueError(f"missing fields: {sorted(missing)}")
                rows.append({
                    key: int(value) if key in INTEGER_FIELDS else value
                    for key, value in raw.items()
                })
    return rows


def _aggregate(rows: list[dict[str, int | str]], omitted_batch: int | None) -> dict:
    frozen = frozen_render()
    target = frozen["amplitude_gauge_and_score"]["frozen_target"]
    by_geometry = {geometry: [] for geometry in GEOMETRY_ORDER}
    for row in rows:
        geometry = str(row["geometry_id"])
        if geometry not in by_geometry:
            raise ValueError(f"unknown geometry {geometry}")
        if omitted_batch is None or int(row["batch"]) != omitted_batch:
            by_geometry[geometry].append(row)

    h_prime = math.sqrt(3) / (3 * math.pi)
    estimates = []
    for index, geometry in enumerate(GEOMETRY_ORDER):
        selected = by_geometry[geometry]
        if not selected:
            raise ValueError(f"missing rows for {geometry}")
        lambda_text = f"{selected[0]['lambda_num']}/{selected[0]['lambda_den']}"
        if lambda_text != target["lambda_order"][index]:
            raise ValueError(f"lambda mismatch for {geometry}")
        invariant_fields = (
            "level", "span_L", "vertices", "edges", "outer_seed",
            "completion_seed", "inner_replicates",
        )
        invariant = tuple(selected[0][field] for field in invariant_fields)
        if any(tuple(row[field] for field in invariant_fields) != invariant for row in selected):
            raise ValueError(f"metadata changed across batches for {geometry}")
        if int(selected[0]["outer_seed"]) == int(selected[0]["completion_seed"]):
            raise ValueError("outer and completion RNG seeds must be distinct")

        samples = sum(int(row["samples"]) for row in selected)
        count = sum(int(row["count_14_23"]) for row in selected)
        inner = int(selected[0]["inner_replicates"])
        sum_delta = sum(int(row["sum_delta_J_14_23"]) for row in selected)
        sum_outer_delta2 = sum(
            int(row["sum_delta_J_inner_square_14_23"]) for row in selected
        )
        sum_individual_delta2 = sum(
            int(row["sum_delta_J2_individual_14_23"]) for row in selected
        )
        if samples <= 1 or count == 0:
            raise ValueError(f"insufficient samples/events for {geometry}")

        measure = sum_delta / (2 * inner * samples)
        d_log_probability = sum_delta / (2 * inner * count)
        span = int(selected[0]["span_L"])
        k_prefactor = float(frozen["frozen_geometries"][index]["K_decimal"])
        z = (
            d_log_probability
            + 4 * h_prime * math.log(span)
            - 2 * h_prime * math.log(k_prefactor)
        )

        sum_y2 = sum_outer_delta2 / (4 * inner * inner)
        outer_variance = (sum_y2 - samples * measure * measure) / (samples - 1)
        conditional_variance_sum = None
        completion_noise = None
        ideal_variance = None
        if inner > 1:
            conditional_variance_sum = (
                sum_individual_delta2 - sum_outer_delta2 / inner
            ) / (inner - 1)
            completion_noise = conditional_variance_sum / (4 * inner * samples)
            ideal_variance = outer_variance - completion_noise

        edges = int(selected[0]["edges"])
        revealed = sum(int(row["sum_revealed"]) for row in selected)
        revealed_high = sum(int(row["sum_revealed_14_23"]) for row in selected)
        estimates.append({
            "geometry_id": geometry,
            "lambda": lambda_text,
            "span_L": span,
            "samples": samples,
            "events_14_23": count,
            "probability_14_23": count / samples,
            "inner_replicates": inner,
            "measure_tangent_14_23": measure,
            "d_log_probability": d_log_probability,
            "z_before_amplitude_projection": z,
            "mean_revealed_edges": revealed / samples,
            "mean_revealed_fraction": revealed / (samples * edges),
            "mean_revealed_edges_on_14_23": revealed_high / count,
            "max_revealed_edges": max(int(row["max_revealed"]) for row in selected),
            "local_outer_estimator_sample_variance": outer_variance,
            "estimated_completion_noise_in_outer_variance": completion_noise,
            "estimated_ideal_stopped_outer_variance": ideal_variance,
            "mean_conditional_delta_J_variance_on_14_23": (
                conditional_variance_sum / count
                if conditional_variance_sum is not None else None
            ),
        })

    anchored = [
        row["z_before_amplitude_projection"]
        - estimates[ANCHOR_INDEX]["z_before_amplitude_projection"]
        for row in estimates
    ]
    residual_full = [
        value - target_value
        for value, target_value in zip(anchored, target["anchored_dQ_logU"])
    ]
    active = [index for index in range(4) if index != ANCHOR_INDEX]
    return {
        "estimates": estimates,
        "anchored_lattice_tangent": anchored,
        "frozen_target": target["anchored_dQ_logU"],
        "residual": [residual_full[index] for index in active],
    }


def score(rows: list[dict[str, int | str]]) -> dict:
    if not rows:
        raise ValueError("no rows")
    for field in ("level", "outer_seed", "completion_seed", "inner_replicates"):
        if len({int(row[field]) for row in rows}) != 1:
            raise ValueError(f"{field} must be shared across all geometries")
    batch_sets = {
        geometry: {
            int(row["batch"]) for row in rows if row["geometry_id"] == geometry
        }
        for geometry in GEOMETRY_ORDER
    }
    first = batch_sets[GEOMETRY_ORDER[0]]
    if len(first) < 3 or any(batch_sets[geometry] != first for geometry in GEOMETRY_ORDER):
        raise ValueError("all geometries need at least three synchronized batches")
    reference_counters = {
        int(row["batch"]): (int(row["sample_begin"]), int(row["samples"]))
        for row in rows if row["geometry_id"] == GEOMETRY_ORDER[0]
    }
    for geometry in GEOMETRY_ORDER:
        selected = [row for row in rows if row["geometry_id"] == geometry]
        if len(selected) != len(first):
            raise ValueError(f"duplicate batch rows for {geometry}")
        counters = {
            int(row["batch"]): (int(row["sample_begin"]), int(row["samples"]))
            for row in selected
        }
        if counters != reference_counters:
            raise ValueError("batch counter domains are not synchronized")
    batch_ids = sorted(first)
    full = _aggregate(rows, None)
    deleted = [_aggregate(rows, batch)["residual"] for batch in batch_ids]
    covariance = _jackknife_covariance(full["residual"], deleted)
    return {
        "schema": "matching-one.p263-local-stopped-qscore-score.v1",
        "issue": 263,
        "status": "exploratory_mechanism_pilot_not_primary_evidence",
        "batch_count": len(batch_ids),
        "state_order": ["lambda_1_4", "lambda_2_3", "lambda_3_4"],
        **full,
        "residual_covariance": covariance,
        "covariance_trace": sum(covariance[index][index] for index in range(3)),
        "joint_gls": _gls(full["residual"], covariance),
        "exact_identity": (
            "For event-determining transcript T and independent completion C, "
            "E[I(T)*(J(C<-T)-J(C))/2]=Cov(I,J/2)."
        ),
        "warnings": [
            "The covariance numerator is unbiased; d_log_probability is a finite-sample ratio estimator.",
            "Finite inner replication adds completion noise even though untouched edges cancel sample by sample.",
            "Arm/event-scale variance is a scaling conjecture; transcript size can remain large in finite rectangles.",
            "This secondary must not replace or be combined with the revealed Phase-E primary score.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = score(read_rows(args.batches))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
