#!/usr/bin/env python3
"""Score the frozen one-generator norm-4 Hermite--Krawtchouk thermal jet."""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Sequence

import mpmath as mp

from analyze_matching_parity_derivatives_fast import H, combine, read, remove
from analyze_norm4_variance_pilot import covariance_of_mean
from analyze_rank_gap_thermal_window import FIELDS
from hermite_krawtchouk_scaling_jet import (
    canonical_dimensionless_width,
    scaling_derivative_jet,
    width_normalized_jet,
)
from score_p50_fullcurve_n290 import generalized_covariance_score
from threshold_score_modes import project


SOURCE_ORDER = (65, 130, 85, 170)
SIZE_ORDER = (65, 130, 260, 85, 170, 340)
ORDERS = (2, 3, 4, 5, 6)


@dataclass(frozen=True)
class SourceRun:
    n: int
    base_histogram: Path
    base_moments: Path
    extension_histogram: Path
    extension_moments: Path


@dataclass(frozen=True)
class TargetRun:
    n: int
    histogram: Path
    moments: Path


def parse_source(text: str) -> SourceRun:
    fields = text.split(":", 4)
    if len(fields) != 5:
        raise argparse.ArgumentTypeError(
            "source must be N:BASE_HIST:BASE_MOMENTS:EXT_HIST:EXT_MOMENTS"
        )
    return SourceRun(int(fields[0]), *(Path(value) for value in fields[1:]))


def parse_target(text: str) -> TargetRun:
    fields = text.split(":", 2)
    if len(fields) != 3:
        raise argparse.ArgumentTypeError("target must be N:HISTOGRAM:MOMENTS")
    return TargetRun(int(fields[0]), Path(fields[1]), Path(fields[2]))


def merge_histogram_blocks(paths: Sequence[Path], n: int):
    blocks = [read(path) for path in paths]
    keys = set(blocks[0])
    if any(set(block) != keys for block in blocks[1:]) or {key[0] for key in keys} != {n}:
        raise ValueError(f"N={n}: histogram blocks have different keys")
    output = {}
    for key in keys:
        rows = [block[key] for block in blocks]
        first = rows[0]
        if any(
            (row.a, row.b, row.orientation, row.batch)
            != (first.a, first.b, first.orientation, first.batch)
            for row in rows[1:]
        ):
            raise ValueError(f"N={n}: histogram descriptors differ")
        output[key] = H(
            n,
            first.a,
            first.b,
            first.orientation,
            first.batch,
            sum(row.samples for row in rows),
            [sum(row.minus[k] for row in rows) for k in range(n + 1)],
            [sum(row.plus[k] for row in rows) for k in range(n + 1)],
        )
    return output


def read_gap_batches(paths: Sequence[Path], n: int):
    output = {(orientation, batch): {"samples": 0, "sum_gap": 0}
              for orientation in ("first", "second") for batch in range(100)}
    for path in paths:
        seen = set()
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            required = {"n", "orientation", "batch", "samples", *FIELDS}
            if required - set(reader.fieldnames or ()):
                raise ValueError(f"{path}: incomplete moments schema")
            for raw in reader:
                if int(raw["n"]) != n:
                    raise ValueError(f"{path}: unexpected size")
                key = (raw["orientation"], int(raw["batch"]))
                if key not in output or key in seen:
                    raise ValueError(f"{path}: invalid orientation/batch")
                seen.add(key)
                output[key]["samples"] += int(raw["samples"])
                output[key]["sum_gap"] += int(raw["sum_gap"])
        if seen != set(output):
            raise ValueError(f"{path}: moments batches are incomplete")
    return output


def grouped_histograms(records, n: int):
    return {
        orientation: [
            records[key]
            for key in sorted(records)
            if key[0] == n and key[1] == orientation
        ]
        for orientation in ("first", "second")
    }


def thermal_coefficients(point):
    return [
        mp.mpf(point["P4_D_modes"][order] if order % 2 == 0 else point["P4_S_modes"][order])
        for order in range(max(ORDERS) + 1)
    ]


def jet_state(groups, gap_batches, n: int, omitted: int = -1):
    totals = {
        orientation: combine(
            [row for row in groups[orientation] if row.batch != omitted]
        )
        for orientation in ("first", "second")
    }
    projected = project(totals["first"], totals["second"], max(ORDERS))
    jet = scaling_derivative_jet(
        thermal_coefficients(projected), n, mp.mpf(projected["p0"]), mp.mpf(13) / 8
    )
    selected_gaps = [
        value
        for (orientation, batch), value in gap_batches.items()
        if batch != omitted
    ]
    gap_mean = mp.mpf(sum(value["sum_gap"] for value in selected_gaps)) / sum(
        value["samples"] for value in selected_gaps
    )
    width = canonical_dimensionless_width(n, gap_mean)
    normalized = width_normalized_jet(jet, width)
    return [float(normalized[order]) for order in ORDERS]


def estimate_aligned(groups, gaps, sizes):
    points = {n: jet_state(groups[n], gaps[n], n) for n in sizes}
    full = [value for n in sizes for value in points[n]]
    pseudovalue_rows = []
    for batch in range(100):
        deleted = [
            value
            for n in sizes
            for value in jet_state(groups[n], gaps[n], n, omitted=batch)
        ]
        pseudovalue_rows.append(
            [100 * a - 99 * b for a, b in zip(full, deleted)]
        )
    return points, covariance_of_mean(pseudovalue_rows)


def assemble_covariance(source_cov, target_covariances):
    dimension = len(SIZE_ORDER) * len(ORDERS)
    output = [[0.0] * dimension for _ in range(dimension)]
    for i, n_i in enumerate(SOURCE_ORDER):
        for j, n_j in enumerate(SOURCE_ORDER):
            p_i, p_j = SIZE_ORDER.index(n_i), SIZE_ORDER.index(n_j)
            for a in range(len(ORDERS)):
                for b in range(len(ORDERS)):
                    output[p_i * 5 + a][p_j * 5 + b] = source_cov[i * 5 + a][j * 5 + b]
    for n, covariance in target_covariances.items():
        position = SIZE_ORDER.index(n)
        for a in range(5):
            for b in range(5):
                output[position * 5 + a][position * 5 + b] = covariance[a][b]
    return output


def cocycle_transform(multiplier: float):
    output = []
    for lineage_start in (0, 3):
        for mode in range(len(ORDERS)):
            row = [0.0] * (len(SIZE_ORDER) * len(ORDERS))
            row[(lineage_start + 0) * 5 + mode] = multiplier - 1.0
            row[(lineage_start + 1) * 5 + mode] = -multiplier
            row[(lineage_start + 2) * 5 + mode] = 1.0
            output.append(row)
    return output


def transform(vector, covariance, matrix):
    residual = [math.fsum(row[i] * vector[i] for i in range(len(vector))) for row in matrix]
    residual_covariance = [
        [
            math.fsum(
                matrix[i][a] * covariance[a][b] * matrix[j][b]
                for a in range(len(vector))
                for b in range(len(vector))
            )
            for j in range(len(matrix))
        ]
        for i in range(len(matrix))
    ]
    return residual, residual_covariance


def model_score(vector, covariance, multiplier, label):
    residual, residual_covariance = transform(
        vector, covariance, cocycle_transform(multiplier)
    )
    return {
        "label": label,
        "cocycle_multiplier": multiplier,
        "residual_order": [
            f"N{start}_r{order}" for start in (65, 85) for order in ORDERS
        ],
        "residual": residual,
        "covariance": residual_covariance,
        "score": generalized_covariance_score(residual, residual_covariance),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run", action="append", type=parse_source, required=True)
    parser.add_argument("--target-run", action="append", type=parse_target, required=True)
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    sources = {run.n: run for run in args.source_run}
    targets = {run.n: run for run in args.target_run}
    if tuple(sources) != SOURCE_ORDER or tuple(targets) != (260, 340):
        raise SystemExit("source order must be 65,130,85,170 and targets 260,340")

    groups = {}
    gaps = {}
    for n, run in sources.items():
        groups[n] = grouped_histograms(
            merge_histogram_blocks((run.base_histogram, run.extension_histogram), n), n
        )
        gaps[n] = read_gap_batches((run.base_moments, run.extension_moments), n)
    source_points, source_covariance = estimate_aligned(groups, gaps, SOURCE_ORDER)
    target_points = {}
    target_covariances = {}
    for n, run in targets.items():
        target_group = grouped_histograms(merge_histogram_blocks((run.histogram,), n), n)
        target_gap = read_gap_batches((run.moments,), n)
        point, covariance = estimate_aligned({n: target_group}, {n: target_gap}, (n,))
        target_points[n] = point[n]
        target_covariances[n] = covariance

    points = {**source_points, **target_points}
    vector = [value for n in SIZE_ORDER for value in points[n]]
    covariance = assemble_covariance(source_covariance, target_covariances)
    payload = {
        "schema": "matching-one/norm4-one-generator-thermal-jet/v1",
        "status": "prospective_correlated_secondary_score",
        "size_order": list(SIZE_ORDER),
        "orders": list(ORDERS),
        "coordinate": "finite-N Hermite-Krawtchouk thermal jet times w_can^r",
        "point": {str(n): points[n] for n in SIZE_ORDER},
        "models_in_frozen_order": [
            model_score(vector, covariance, 1.5, "q2_one_even_generator"),
            model_score(vector, covariance, 2.0, "Jordan_one_even_generator"),
        ],
        "evidence_guard": (
            "One common cocycle multiplier acts on ranks 2..6. No rank-specific "
            "amplitudes are fitted. This vector shares data with the scalar U score."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
