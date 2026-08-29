#!/usr/bin/env python3
"""Project norm-4 pilot variance and CPU cost without scoring pilot means."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from analyze_p48_retrospective import project_size, read_histograms
from score_p50_fullcurve_n290 import grouped
from score_p49_fullcurve_doubling import aggregate, orientation_values


COMPLETION_GAPS = {260: 0.1229742463, 340: 0.2046281595}
PRODUCTION_SAMPLES = (500_000_000, 1_000_000_000, 2_000_000_000)
SOURCE_PRODUCTION_SAMPLES = (100_000_000, 500_000_000, 1_000_000_000, 2_000_000_000)
SOURCE_SIZES = (65, 130, 85, 170)


def u_value(n: int, state: Mapping[str, float]) -> float:
    return math.pow(n, 13.0 / 8.0) * float(state["P4_S_prime"]) / float(
        state["mean_slope"]
    )


def jackknife_se(full: float, deleted: list[float]) -> float:
    batches = len(deleted)
    if batches < 2:
        raise ValueError("at least two delete-one replicates are required")
    pseudo = [batches * full - (batches - 1) * value for value in deleted]
    mean = math.fsum(pseudo) / batches
    variance = math.fsum((value - mean) ** 2 for value in pseudo) / (
        batches * (batches - 1)
    )
    return math.sqrt(variance)


def covariance_of_mean(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    batches = len(rows)
    if batches < 2 or not rows or not rows[0]:
        raise ValueError("jackknife covariance needs at least two nonempty rows")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("jackknife rows have inconsistent widths")
    means = [math.fsum(row[j] for row in rows) / batches for j in range(width)]
    scale = 1.0 / (batches * (batches - 1))
    return [
        [
            scale
            * math.fsum(
                (row[i] - means[i]) * (row[j] - means[j]) for row in rows
            )
            for j in range(width)
        ]
        for i in range(width)
    ]


def quadratic_2(vector: Sequence[float], covariance: Sequence[Sequence[float]]) -> float:
    a, b = covariance[0]
    c, d = covariance[1]
    determinant = a * d - b * c
    if determinant <= 0.0:
        raise ValueError("2x2 covariance is not positive definite")
    x, y = vector
    return (d * x * x - (b + c) * x * y + a * y * y) / determinant


def transform_covariance(
    matrix: Sequence[Sequence[float]], transform: Sequence[Sequence[float]]
) -> list[list[float]]:
    return [
        [
            math.fsum(
                transform[i][a] * matrix[a][b] * transform[j][b]
                for a in range(len(matrix))
                for b in range(len(matrix))
            )
            for j in range(len(transform))
        ]
        for i in range(len(transform))
    ]


def u_state(by_orientation, *, omitted: int = -1) -> dict[str, float]:
    """Return only the center observables needed for U.

    At N=340 the standard-double binomial recurrence underflows at p=0.9,
    outside the physical region. ``project_size`` bisects inward from p=0.5
    without evaluating that endpoint, so its center remains stable.
    """
    projected = project_size(by_orientation, omitted)
    p0 = projected["p0"]
    rows = {
        name: aggregate(by_orientation[name], omitted)
        for name in ("first", "second")
    }
    n = by_orientation["first"][0].n
    slopes = [
        orientation_values(n, rows[name], p0)["M_prime"]
        for name in ("first", "second")
    ]
    return {
        "mean_slope": math.fsum(slopes) / 2.0,
        "P4_S_prime": float(projected["P4_S_prime"]),
    }


def source_covariance(
    histograms: Sequence[Path], metadata_paths: Sequence[Path]
) -> dict[str, object]:
    records = {}
    for path in histograms:
        incoming = read_histograms(path)
        overlap = set(records) & set(incoming)
        if overlap:
            raise ValueError(f"duplicate source histogram keys: {sorted(overlap)[:2]}")
        records.update(incoming)
    found = tuple(sorted({key[0] for key in records}))
    if found != tuple(sorted(SOURCE_SIZES)):
        raise ValueError(f"expected source sizes {SOURCE_SIZES}, got {found}")
    groups = {n: grouped(records, n) for n in SOURCE_SIZES}
    batch_ids = [row.batch for row in groups[SOURCE_SIZES[0]]["first"]]
    for n in SOURCE_SIZES[1:]:
        if [row.batch for row in groups[n]["first"]] != batch_ids:
            raise ValueError("source batches are not aligned across sizes")
    full = [u_value(n, u_state(groups[n])) for n in SOURCE_SIZES]
    pseudovalue_rows = []
    batches = len(batch_ids)
    for batch in batch_ids:
        deleted = [
            u_value(n, u_state(groups[n], omitted=batch)) for n in SOURCE_SIZES
        ]
        pseudovalue_rows.append(
            [batches * x - (batches - 1) * y for x, y in zip(full, deleted)]
        )
    metadata = [json.loads(path.read_text(encoding="utf-8")) for path in metadata_paths]
    signatures = {
        (
            int(row["seed"]),
            int(row["replica_counter_first"]),
            int(row["replica_counter_last_exclusive"]),
            int(row["samples_per_pair"]),
            int(row["batches"]),
        )
        for row in metadata
    }
    if len(signatures) != 1:
        raise ValueError("source metadata do not declare one aligned RNG block")
    cpu_seconds_per_million = math.fsum(
        float(row["elapsed_seconds"])
        * int(row["threads_requested"])
        * 1_000_000
        / int(row["samples_per_pair"])
        for row in metadata
    )
    return {
        "order": list(SOURCE_SIZES),
        "samples": groups[SOURCE_SIZES[0]]["first"][0].samples * batches,
        "batches": batches,
        "U": full,
        "covariance": covariance_of_mean(pseudovalue_rows),
        "cpu_seconds_per_million_all_four_sizes": cpu_seconds_per_million,
        "metadata": [str(path) for path in metadata_paths],
    }


def joint_forecasts(source: Mapping[str, object], targets: Sequence[Mapping[str, object]]):
    source_cov = source["covariance"]
    q2_transform = ((-0.5, 1.5, 0.0, 0.0), (0.0, 0.0, -0.5, 1.5))
    jordan_transform = ((-1.0, 2.0, 0.0, 0.0), (0.0, 0.0, -1.0, 2.0))
    source_q2 = transform_covariance(source_cov, q2_transform)
    source_jordan = transform_covariance(source_cov, jordan_transform)
    gaps = [COMPLETION_GAPS[int(row["N"])] for row in targets]
    output = []
    source_samples_observed = int(source["samples"])
    source_cpu_rate = float(source["cpu_seconds_per_million_all_four_sizes"])
    target_cpu_rate = math.fsum(float(row["cpu_seconds_per_million"]) for row in targets)
    for source_samples in SOURCE_PRODUCTION_SAMPLES:
        source_scale = source_samples_observed / source_samples
        for target_samples in PRODUCTION_SAMPLES:
            target_variances = [
                float(row["jackknife_se_U"]) ** 2
                * float(row["samples"])
                / target_samples
                for row in targets
            ]
            covariance_q2 = [
                [source_scale * value for value in row] for row in source_q2
            ]
            covariance_jordan = [
                [source_scale * value for value in row] for row in source_jordan
            ]
            for i, variance in enumerate(target_variances):
                covariance_q2[i][i] += variance
                covariance_jordan[i][i] += variance
            q2_chi2 = quadratic_2(gaps, covariance_q2)
            jordan_chi2 = quadratic_2(gaps, covariance_jordan)
            output.append(
                {
                    "samples_per_target": target_samples,
                    "samples_per_source_size_total": source_samples,
                    "additional_source_samples_per_size": max(
                        0, source_samples - source_samples_observed
                    ),
                    "projected_incremental_cpu_seconds": (
                        max(0, source_samples - source_samples_observed)
                        / 1_000_000
                        * source_cpu_rate
                        + target_samples / 1_000_000 * target_cpu_rate
                    ),
                    "target_covariance": [
                        [target_variances[0], 0.0],
                        [0.0, target_variances[1]],
                    ],
                    "q2_null_if_jordan_true": {
                        "expected_noncentral_chi_square": q2_chi2,
                        "sqrt_noncentrality": math.sqrt(q2_chi2),
                        "covariance": covariance_q2,
                    },
                    "jordan_null_if_q2_true": {
                        "expected_noncentral_chi_square": jordan_chi2,
                        "sqrt_noncentrality": math.sqrt(jordan_chi2),
                        "covariance": covariance_jordan,
                    },
                }
            )
    return output


def analyze_one(histogram: Path, metadata: Path, n: int) -> dict[str, object]:
    data = read_histograms(histogram)
    by_orientation = grouped(data, n)
    point = u_state(by_orientation)
    point_u = u_value(n, point)
    batch_ids = [row.batch for row in by_orientation["first"]]
    deleted_u = [
        u_value(
            n,
            u_state(by_orientation, omitted=batch),
        )
        for batch in batch_ids
    ]
    se_u = jackknife_se(point_u, deleted_u)
    meta = json.loads(metadata.read_text(encoding="utf-8"))
    samples = int(meta["samples_per_pair"])
    elapsed = float(meta["elapsed_seconds"])
    threads = int(meta["threads_requested"])
    forecasts = []
    for target_samples in PRODUCTION_SAMPLES:
        projected_se = se_u * math.sqrt(samples / target_samples)
        forecasts.append(
            {
                "samples": target_samples,
                "target_only_projected_se_U": projected_se,
                "completion_gap": COMPLETION_GAPS[n],
                "optimistic_gap_sigma": COMPLETION_GAPS[n] / projected_se,
                "projected_wall_seconds_at_same_threads": elapsed
                * target_samples
                / samples,
                "projected_cpu_seconds": elapsed
                * threads
                * target_samples
                / samples,
            }
        )
    return {
        "N": n,
        "samples": samples,
        "batches": len(batch_ids),
        "mean_slope": float(point["mean_slope"]),
        "P4_S_prime": float(point["P4_S_prime"]),
        "U": point_u,
        "jackknife_se_U": se_u,
        "elapsed_seconds": elapsed,
        "threads": threads,
        "cpu_seconds_per_million": elapsed * threads * 1_000_000 / samples,
        "forecasts": forecasts,
        "metadata": str(metadata),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n260-hist", type=Path, required=True)
    parser.add_argument("--n260-metadata", type=Path, required=True)
    parser.add_argument("--n340-hist", type=Path, required=True)
    parser.add_argument("--n340-metadata", type=Path, required=True)
    parser.add_argument("--source-histograms", type=Path, nargs=4, required=True)
    parser.add_argument("--source-metadata", type=Path, nargs=4, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = [
        analyze_one(args.n260_hist, args.n260_metadata, 260),
        analyze_one(args.n340_hist, args.n340_metadata, 340),
    ]
    source = source_covariance(args.source_histograms, args.source_metadata)
    payload = {
        "schema": "matching-one/norm4-variance-pilot/v1",
        "status": "variance and cost planning only; pilot means are not model scores",
        "sizes": {str(row["N"]): row for row in rows},
        "source": source,
        "joint_forecasts": joint_forecasts(source, rows),
        "production_recommendation": {
            "samples_per_source_size_total": 2000000000,
            "additional_source_samples_per_size": 1900000000,
            "samples_per_target_size": 1000000000,
            "reason": (
                "lowest tested CPU cost with approximately three-sigma expected "
                "separation in both directional fixed-model null scores"
            ),
        },
        "forecast_scope": (
            "inverse-sample target variance scaling plus the full frozen common-"
            "random-number covariance of N=65,130,85,170 source U values"
        ),
        "decision": (
            "use the smallest scale that retains useful separation after source "
            "covariance is added; arithmetic-class closure may independently justify promotion"
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
