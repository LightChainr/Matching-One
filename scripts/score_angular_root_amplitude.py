#!/usr/bin/env python3
"""Score the frozen angular-normalized threshold-root amplitude."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, Iterable, Sequence

import mpmath as mp

from analyze_threshold_rank_orientation import (
    Key,
    add_histograms,
    cos4,
    read_histograms,
    validate_moments,
)
from analyze_threshold_ranks import matching_derivative, matching_root, matching_value


METRICS = ("delta_M", "mean_M_prime", "root_gap", "closure_C", "A_M", "B", "A_p")


def merge_inputs(
    histogram_paths: Sequence[Path], moment_paths: Sequence[Path]
) -> Dict[Key, Dict[str, object]]:
    if len(histogram_paths) != len(moment_paths) or not histogram_paths:
        raise ValueError("histogram and moment path counts must agree")
    merged: Dict[Key, Dict[str, object]] = {}
    for histogram_path, moment_path in zip(histogram_paths, moment_paths):
        records = read_histograms(histogram_path)
        validate_moments(moment_path, records)
        overlap = set(merged).intersection(records)
        if overlap:
            raise ValueError(f"duplicate histogram keys: {sorted(overlap)[:3]}")
        merged.update(records)
    return merged


def validate_metadata(
    paths: Sequence[Path], sizes: Sequence[int], source_commit: str
) -> dict[str, object]:
    if len(paths) != len(sizes):
        raise ValueError("one metadata file is required per size")
    payloads = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
    common_fields = (
        "git_commit",
        "seed",
        "replica_counter_first",
        "replica_counter_last_exclusive",
        "samples_per_pair",
        "batches",
        "rng",
    )
    for field in common_fields:
        values = {json.dumps(payload.get(field), sort_keys=True) for payload in payloads}
        if len(values) != 1:
            raise ValueError(f"metadata field differs across sizes: {field}")
    if payloads[0]["git_commit"] != source_commit or len(source_commit) != 40:
        raise ValueError("metadata does not identify the declared full source commit")
    seen_sizes = []
    for payload in payloads:
        designs = payload.get("designs")
        if not isinstance(designs, list) or len(designs) != 1:
            raise ValueError("each metadata file must contain exactly one design")
        seen_sizes.append(int(designs[0]["N"]))
    if sorted(seen_sizes) != sorted(sizes):
        raise ValueError("metadata sizes do not match histogram sizes")
    return {
        field: payloads[0][field] for field in common_fields
    } | {
        "per_size_elapsed_seconds": {
            str(int(payload["designs"][0]["N"])): payload["elapsed_seconds"]
            for payload in payloads
        },
        "commands": [payload["command"] for payload in payloads],
    }


def aggregate_orientation(
    records: Dict[Key, Dict[str, object]], n: int, orientation: str, drop_batch: int | None
) -> dict[str, object]:
    selected = [
        records[key]
        for key in sorted(records)
        if key[0] == n and key[1] == orientation and key[2] != drop_batch
    ]
    if len(selected) < 2:
        raise ValueError("insufficient batches after deletion")
    return {
        "a": int(selected[0]["a"]),
        "b": int(selected[0]["b"]),
        "samples": sum(int(row["samples"]) for row in selected),
        "minus": add_histograms(selected, "minus"),
        "plus": add_histograms(selected, "plus"),
    }


def estimate(
    records: Dict[Key, Dict[str, object]], n: int, p: mp.mpf, drop_batch: int | None = None
) -> dict[str, float]:
    first = aggregate_orientation(records, n, "first", drop_batch)
    second = aggregate_orientation(records, n, "second", drop_batch)
    first_m = matching_value(n, first["samples"], first["minus"], first["plus"], p)
    second_m = matching_value(n, second["samples"], second["minus"], second["plus"], p)
    first_d = matching_derivative(n, first["samples"], first["minus"], first["plus"], p)
    second_d = matching_derivative(n, second["samples"], second["minus"], second["plus"], p)
    first_root = matching_root(n, first["samples"], first["minus"], first["plus"])
    second_root = matching_root(n, second["samples"], second["minus"], second["plus"])
    delta_m = first_m - second_m
    slope = (first_d + second_d) / 2
    root_gap = first_root - second_root
    delta_cos4 = cos4(first["a"], first["b"]) - cos4(second["a"], second["b"])
    a_m = mp.power(n, mp.mpf(13) / 8) * delta_m / delta_cos4
    b_value = mp.power(n, -mp.mpf(3) / 8) * slope
    return {
        "delta_M": float(delta_m),
        "mean_M_prime": float(slope),
        "root_gap": float(root_gap),
        "closure_C": float(-root_gap * slope / delta_m),
        "A_M": float(a_m),
        "B": float(b_value),
        "A_p": float(-n * n * root_gap / delta_cos4),
        "A_p_predicted_from_A_M_over_B": float(a_m / b_value),
        "delta_cos4": delta_cos4,
    }


def jackknife_covariance(xs: Sequence[float], ys: Sequence[float]) -> float:
    if len(xs) != len(ys) or len(xs) < 2:
        raise ValueError("jackknife vectors must have the same nontrivial length")
    mean_x = math.fsum(xs) / len(xs)
    mean_y = math.fsum(ys) / len(ys)
    return (len(xs) - 1) / len(xs) * math.fsum(
        (x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)
    )


def invert_2x2(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    a, b = matrix[0]
    c, d = matrix[1]
    determinant = a * d - b * c
    if determinant <= 0.0:
        raise ValueError("covariance matrix is not positive definite")
    return [[d / determinant, -b / determinant], [-c / determinant, a / determinant]]


def quadratic(values: Sequence[float], inverse: Sequence[Sequence[float]]) -> float:
    return math.fsum(
        values[i] * inverse[i][j] * values[j]
        for i in range(len(values)) for j in range(len(values))
    )


def score(
    records: Dict[Key, Dict[str, object]], p: mp.mpf, prediction: float,
    prediction_se: float,
) -> dict[str, object]:
    sizes = sorted({key[0] for key in records})
    if sizes != [65, 85]:
        raise ValueError(f"frozen primary sizes must be [65, 85], got {sizes}")
    batch_ids = {
        n: sorted({key[2] for key in records if key[0] == n}) for n in sizes
    }
    expected = list(range(len(batch_ids[65])))
    if any(batch_ids[n] != expected for n in sizes):
        raise ValueError("sizes must have the same complete zero-based batch ids")

    full = {n: estimate(records, n, p) for n in sizes}
    delete_one = {
        n: [estimate(records, n, p, drop_batch=batch) for batch in expected]
        for n in sizes
    }
    by_size: dict[str, object] = {}
    for n in sizes:
        covariance = [
            [
                jackknife_covariance(
                    [row[left] for row in delete_one[n]],
                    [row[right] for row in delete_one[n]],
                )
                for right in METRICS
            ]
            for left in METRICS
        ]
        by_size[str(n)] = {
            "estimate": full[n],
            "metric_order": list(METRICS),
            "jackknife_covariance": covariance,
            "standard_errors": {
                metric: math.sqrt(max(0.0, covariance[i][i]))
                for i, metric in enumerate(METRICS)
            },
        }

    a_p_covariance = [
        [
            jackknife_covariance(
                [row["A_p"] for row in delete_one[left_n]],
                [row["A_p"] for row in delete_one[right_n]],
            )
            for right_n in sizes
        ]
        for left_n in sizes
    ]
    observations = [full[n]["A_p"] for n in sizes]
    prediction_covariance = [
        [a_p_covariance[i][j] + prediction_se**2 for j in range(2)]
        for i in range(2)
    ]
    frozen_residuals = [value - prediction for value in observations]
    frozen_chi2 = quadratic(frozen_residuals, invert_2x2(prediction_covariance))
    zero_chi2 = quadratic(observations, invert_2x2(a_p_covariance))

    inverse = invert_2x2(a_p_covariance)
    inverse_times_one = [sum(row) for row in inverse]
    normalization = sum(inverse_times_one)
    weights = [value / normalization for value in inverse_times_one]
    free_amplitude = math.fsum(w * y for w, y in zip(weights, observations))
    free_residuals = [value - free_amplitude for value in observations]
    free_chi2 = quadratic(free_residuals, inverse)

    return {
        "schema": "frozen angular-normalized root amplitude score v1",
        "p_ref": mp.nstr(p, mp.mp.dps),
        "sizes": sizes,
        "batch_count": len(expected),
        "by_size": by_size,
        "A_p_cross_size_jackknife_covariance": a_p_covariance,
        "frozen_prediction": {
            "value": prediction,
            "source_standard_error": prediction_se,
            "residuals": frozen_residuals,
            "residual_covariance": prediction_covariance,
            "chi_square": frozen_chi2,
            "degrees_of_freedom": 2,
        },
        "free_common_amplitude": {
            "value": free_amplitude,
            "standard_error": math.sqrt(1.0 / normalization),
            "weights": weights,
            "internal_chi_square": free_chi2,
            "degrees_of_freedom": 1,
        },
        "zero_effect": {"chi_square": zero_chi2, "degrees_of_freedom": 2},
    }


def write_csv(path: Path, result: dict[str, object]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        fields = ["N", *METRICS, *(metric + "_se" for metric in METRICS)]
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        by_size = result["by_size"]
        assert isinstance(by_size, dict)
        for n in result["sizes"]:
            row = by_size[str(n)]
            writer.writerow({
                "N": n,
                **{metric: row["estimate"][metric] for metric in METRICS},
                **{metric + "_se": row["standard_errors"][metric] for metric in METRICS},
            })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histograms", nargs="+", required=True, type=Path)
    parser.add_argument("--moments", nargs="+", required=True, type=Path)
    parser.add_argument("--metadata", nargs="+", required=True, type=Path)
    parser.add_argument("--p", default="0.592746050790")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--prediction", type=float, required=True)
    parser.add_argument("--prediction-se", type=float, required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument("--binary-sha256", required=True)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    records = merge_inputs(args.histograms, args.moments)
    result = score(records, mp.mpf(args.p), args.prediction, args.prediction_se)
    provenance = validate_metadata(args.metadata, result["sizes"], args.source_commit)
    result["provenance"] = provenance | {
        "source_sha256": args.source_sha256,
        "binary_sha256": args.binary_sha256,
        "tracked_source_clean_before_and_after_build": True,
        "cross_size_rng_policy": "deliberate common stream with aligned batches and full covariance",
    }
    args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    write_csv(args.csv, result)


if __name__ == "__main__":
    main()
