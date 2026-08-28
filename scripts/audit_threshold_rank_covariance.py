#!/usr/bin/env python3
"""Audit cross-size covariance in aligned threshold-rank orientation batches.

The C++ threshold-rank engine may deliberately reuse one counter stream across
different sizes. When batch ids are aligned, this script reconstructs each
batch's orientation difference and slope, and uses delete-one jackknife
pseudo-values for the nonlinear root gap. It then emits covariance matrices of
the estimators across sizes and compares full-covariance and diagonal-only
held-out constant-amplitude tests for

    A_M = N^(13/8) Delta M / Delta cos(4 theta)
    A_p = -N^2 Delta p* / Delta cos(4 theta).

This audits existing aggregates. It never changes the frozen simulation, RNG,
or model-selection protocol.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import mpmath as mp

from analyze_threshold_rank_orientation import (
    add_histograms,
    cos4,
    evaluate_histogram,
    read_histograms,
)
from analyze_threshold_ranks import matching_root


Matrix = List[List[float]]
Vector = List[float]
Record = Dict[str, object]
RecordKey = Tuple[int, str, int]
Grouped = Dict[int, Dict[str, Dict[int, Record]]]


def _matvec(matrix: Matrix, vector: Vector) -> Vector:
    return [math.fsum(a * b for a, b in zip(row, vector)) for row in matrix]


def _solve(matrix: Matrix, vector: Vector) -> Vector:
    n = len(vector)
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError("linear system must be square")
    augmented = [list(map(float, matrix[i])) + [float(vector[i])] for i in range(n)]
    scale = max((abs(value) for row in matrix for value in row), default=0.0)
    tolerance = max(scale * 1e-13, 1e-300)
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= tolerance:
            raise ArithmeticError("singular or ill-resolved matrix")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        for entry in range(column, n + 1):
            augmented[column][entry] /= divisor
        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            for entry in range(column, n + 1):
                augmented[row][entry] -= factor * augmented[column][entry]
    return [augmented[row][-1] for row in range(n)]


def _inverse(matrix: Matrix) -> Matrix:
    n = len(matrix)
    columns = [
        _solve(matrix, [1.0 if row == column else 0.0 for row in range(n)])
        for column in range(n)
    ]
    return [[columns[column][row] for column in range(n)] for row in range(n)]


def _subset(
    matrix: Matrix,
    rows: Sequence[int],
    columns: Optional[Sequence[int]] = None,
) -> Matrix:
    if columns is None:
        columns = rows
    return [[matrix[i][j] for j in columns] for i in rows]


def _quadratic(vector: Vector, inverse_covariance: Matrix) -> float:
    return math.fsum(
        a * b for a, b in zip(vector, _matvec(inverse_covariance, vector))
    )


def covariance_of_mean(
    batch_vectors: Sequence[Sequence[float]],
) -> Tuple[Vector, Matrix]:
    """Return column means and covariance matrix of those means.

    Rows are aligned independent batches and columns are sizes/observables.
    """

    if len(batch_vectors) < 2:
        raise ValueError("at least two aligned batches are required")
    width = len(batch_vectors[0])
    if width == 0 or any(len(row) != width for row in batch_vectors):
        raise ValueError("batch matrix must be nonempty and rectangular")
    if any(not math.isfinite(value) for row in batch_vectors for value in row):
        raise ValueError("batch values must be finite")
    batches = len(batch_vectors)
    means = [
        math.fsum(row[column] for row in batch_vectors) / batches
        for column in range(width)
    ]
    covariance = [[0.0] * width for _ in range(width)]
    denominator = batches * (batches - 1)
    for i in range(width):
        for j in range(i, width):
            value = math.fsum(
                (row[i] - means[i]) * (row[j] - means[j])
                for row in batch_vectors
            ) / denominator
            covariance[i][j] = covariance[j][i] = value
    return means, covariance


def jackknife_pseudovalues(
    full_estimate: float,
    delete_one_estimates: Sequence[float],
) -> List[float]:
    """Convert delete-one estimates into standard jackknife pseudo-values."""

    batches = len(delete_one_estimates)
    if batches < 2:
        raise ValueError("at least two delete-one estimates are required")
    if not math.isfinite(full_estimate) or any(
        not math.isfinite(value) for value in delete_one_estimates
    ):
        raise ValueError("jackknife estimates must be finite")
    return [
        batches * full_estimate - (batches - 1) * value
        for value in delete_one_estimates
    ]


def correlations(covariance: Matrix) -> Matrix:
    output = [[0.0] * len(covariance) for _ in covariance]
    for i in range(len(covariance)):
        if covariance[i][i] <= 0:
            raise ValueError("covariance diagonal must be positive")
        for j in range(len(covariance)):
            output[i][j] = covariance[i][j] / math.sqrt(
                covariance[i][i] * covariance[j][j]
            )
    return output


def constant_heldout_audit(
    values: Sequence[float],
    covariance: Matrix,
    sizes: Sequence[int],
    training_sizes: Sequence[int],
    heldout_sizes: Sequence[int],
) -> Dict[str, object]:
    """Fit one constant on training sizes and score correlated held-out residuals."""

    if len(values) != len(sizes) or len(covariance) != len(sizes):
        raise ValueError("values/covariance/sizes dimensions disagree")
    indices = {size: index for index, size in enumerate(sizes)}
    if len(indices) != len(sizes):
        raise ValueError("sizes must be unique")
    try:
        training = [indices[size] for size in training_sizes]
        heldout = [indices[size] for size in heldout_sizes]
    except KeyError as exc:
        raise ValueError("training or held-out size is absent") from exc

    c_tt = _subset(covariance, training)
    inv_tt = _inverse(c_tt)
    ones_t = [1.0] * len(training)
    raw_weights = _matvec(inv_tt, ones_t)
    denominator = math.fsum(raw_weights)
    if denominator <= 0:
        raise ValueError("constant-fit information must be positive")
    weights = [value / denominator for value in raw_weights]
    amplitude = math.fsum(
        weight * values[index]
        for weight, index in zip(weights, training)
    )
    amplitude_variance = 1.0 / denominator

    c_hh = _subset(covariance, heldout)
    c_ht = _subset(covariance, heldout, training)
    covariance_with_amplitude = _matvec(c_ht, weights)
    residual_covariance = [
        [
            c_hh[i][j]
            + amplitude_variance
            - covariance_with_amplitude[i]
            - covariance_with_amplitude[j]
            for j in range(len(heldout))
        ]
        for i in range(len(heldout))
    ]
    residuals = [values[index] - amplitude for index in heldout]
    chi_square = _quadratic(residuals, _inverse(residual_covariance))
    return {
        "training_sizes": list(training_sizes),
        "heldout_sizes": list(heldout_sizes),
        "amplitude": amplitude,
        "amplitude_se": math.sqrt(amplitude_variance),
        "training_weights": weights,
        "heldout_observed": [values[index] for index in heldout],
        "heldout_residuals": residuals,
        "heldout_residual_covariance": residual_covariance,
        "heldout_chi_square": chi_square,
        "heldout_dof": len(heldout),
    }


def _orientation_batches(
    records: Mapping[RecordKey, Record],
) -> Tuple[List[int], List[int], Grouped]:
    sizes = sorted({key[0] for key in records})
    by_size: Grouped = {}
    common_batches: Optional[List[int]] = None
    for n in sizes:
        by_size[n] = {}
        for orientation in ("first", "second"):
            selected = {
                key[2]: records[key]
                for key in records
                if key[0] == n and key[1] == orientation
            }
            if not selected:
                raise ValueError("N={} has no {} batches".format(n, orientation))
            ids = sorted(selected)
            if ids != list(range(len(ids))):
                raise ValueError(
                    "N={} {} batch ids are not contiguous".format(n, orientation)
                )
            by_size[n][orientation] = selected
            if common_batches is None:
                common_batches = ids
            elif ids != common_batches:
                raise ValueError("all sizes/orientations must share aligned batch ids")
        for batch in common_batches or []:
            first_samples = int(by_size[n]["first"][batch]["samples"])
            second_samples = int(by_size[n]["second"][batch]["samples"])
            if first_samples != second_samples:
                raise ValueError("paired orientations must have equal batch samples")
    if common_batches is None:
        raise ValueError("no batches found")
    return sizes, common_batches, by_size


def _root_from_records(n: int, selected: Sequence[Record]) -> mp.mpf:
    if not selected:
        raise ValueError("cannot reconstruct a root from zero batches")
    minus = add_histograms(selected, "minus")
    plus = add_histograms(selected, "plus")
    samples = sum(int(row["samples"]) for row in selected)
    return matching_root(n, samples, minus, plus)


def _root_gap_pseudovalues(
    sizes: Sequence[int],
    batches: Sequence[int],
    by_size: Grouped,
) -> Tuple[Dict[Tuple[int, int], float], Dict[str, object]]:
    pseudo_by_size_batch: Dict[Tuple[int, int], float] = {}
    details: Dict[str, object] = {}
    for n in sizes:
        first_all = [by_size[n]["first"][batch] for batch in batches]
        second_all = [by_size[n]["second"][batch] for batch in batches]
        full_gap = float(
            _root_from_records(n, first_all) - _root_from_records(n, second_all)
        )
        delete_one: List[float] = []
        for omitted in batches:
            first_reduced = [
                by_size[n]["first"][batch]
                for batch in batches
                if batch != omitted
            ]
            second_reduced = [
                by_size[n]["second"][batch]
                for batch in batches
                if batch != omitted
            ]
            delete_one.append(
                float(
                    _root_from_records(n, first_reduced)
                    - _root_from_records(n, second_reduced)
                )
            )
        pseudovalues = jackknife_pseudovalues(full_gap, delete_one)
        for batch, value in zip(batches, pseudovalues):
            pseudo_by_size_batch[(n, batch)] = value
        bias_corrected = math.fsum(pseudovalues) / len(pseudovalues)
        details[str(n)] = {
            "full_estimate": full_gap,
            "bias_corrected_estimate": bias_corrected,
            "bias_correction": bias_corrected - full_gap,
            "delete_one_count": len(delete_one),
        }
    return pseudo_by_size_batch, details


def reconstruct_batch_metrics(
    records: Mapping[RecordKey, Record],
    p: mp.mpf,
) -> Tuple[List[int], List[Record], Dict[str, object]]:
    sizes, batches, by_size = _orientation_batches(records)
    root_pseudovalues, root_details = _root_gap_pseudovalues(
        sizes, batches, by_size
    )
    output: List[Record] = []
    for batch in batches:
        for n in sizes:
            first = by_size[n]["first"][batch]
            second = by_size[n]["second"][batch]
            first_m, first_d = evaluate_histogram(first, p)
            second_m, second_d = evaluate_histogram(second, p)
            delta_m = first_m - second_m
            mean_slope = (first_d + second_d) / 2
            root_gap_pseudovalue = root_pseudovalues[(n, batch)]
            delta_cos4 = cos4(
                int(first["a"]), int(first["b"])
            ) - cos4(int(second["a"]), int(second["b"]))
            if delta_cos4 == 0:
                raise ValueError("zero angular leverage in a batch")
            output.append(
                {
                    "N": n,
                    "batch": batch,
                    "samples": int(first["samples"]),
                    "a1": int(first["a"]),
                    "b1": int(first["b"]),
                    "a2": int(second["a"]),
                    "b2": int(second["b"]),
                    "delta_cos4": delta_cos4,
                    "delta_M": float(delta_m),
                    "mean_M_prime": float(mean_slope),
                    "root_gap_jackknife_pseudovalue": root_gap_pseudovalue,
                    "A_M": float(
                        n ** (13.0 / 8.0) * delta_m / delta_cos4
                    ),
                    "B": float(n ** (-3.0 / 8.0) * mean_slope),
                    "A_p_jackknife_pseudovalue": float(
                        -n * n * root_gap_pseudovalue / delta_cos4
                    ),
                }
            )
    return sizes, output, root_details


def _matrix_by_field(
    rows: Sequence[Record],
    sizes: Sequence[int],
    field: str,
) -> List[List[float]]:
    by_batch: Dict[int, Dict[int, float]] = {}
    for row in rows:
        by_batch.setdefault(int(row["batch"]), {})[int(row["N"])] = float(
            row[field]
        )
    output = []
    for batch in sorted(by_batch):
        if set(by_batch[batch]) != set(sizes):
            raise ValueError(
                "batch {} is incomplete for field {}".format(batch, field)
            )
        output.append([by_batch[batch][n] for n in sizes])
    return output


def _diagonal(covariance: Matrix) -> Matrix:
    return [
        [covariance[i][i] if i == j else 0.0 for j in range(len(covariance))]
        for i in range(len(covariance))
    ]


def audit(
    records: Mapping[RecordKey, Record],
    p: mp.mpf,
    training_sizes: Sequence[int],
    heldout_sizes: Sequence[int],
) -> Tuple[List[Record], Dict[str, object]]:
    sizes, batch_rows, root_details = reconstruct_batch_metrics(records, p)
    metric_fields = {
        "delta_M": "delta_M",
        "root_gap": "root_gap_jackknife_pseudovalue",
        "mean_M_prime": "mean_M_prime",
        "A_M": "A_M",
        "B": "B",
        "A_p": "A_p_jackknife_pseudovalue",
    }
    payload: Dict[str, object] = {
        "format_version": 2,
        "p": mp.nstr(p, mp.mp.dps),
        "sizes": sizes,
        "batch_count": len({int(row["batch"]) for row in batch_rows}),
        "nonlinear_estimator": {
            "root_gap_method": "delete_one_jackknife_pseudovalues",
            "by_N": root_details,
        },
        "metrics": {},
        "constant_amplitude_audits": {},
    }
    for metric, field in metric_fields.items():
        matrix = _matrix_by_field(batch_rows, sizes, field)
        means, covariance = covariance_of_mean(matrix)
        metric_payload = {
            "batch_field": field,
            "means": dict(zip(map(str, sizes), means)),
            "standard_errors": dict(
                zip(
                    map(str, sizes),
                    [math.sqrt(covariance[i][i]) for i in range(len(sizes))],
                )
            ),
            "covariance_of_means": covariance,
            "correlation_of_batch_values": correlations(covariance),
        }
        metrics_payload = payload["metrics"]
        assert isinstance(metrics_payload, dict)
        metrics_payload[metric] = metric_payload
        if metric in ("A_M", "A_p"):
            audits = payload["constant_amplitude_audits"]
            assert isinstance(audits, dict)
            audits[metric] = {
                "full_covariance": constant_heldout_audit(
                    means, covariance, sizes, training_sizes, heldout_sizes
                ),
                "diagonal_covariance": constant_heldout_audit(
                    means,
                    _diagonal(covariance),
                    sizes,
                    training_sizes,
                    heldout_sizes,
                ),
            }
    return batch_rows, payload


def write_batch_csv(path: Path, rows: Sequence[Record]) -> None:
    fields = [
        "N",
        "batch",
        "samples",
        "a1",
        "b1",
        "a2",
        "b2",
        "delta_cos4",
        "delta_M",
        "mean_M_prime",
        "root_gap_jackknife_pseudovalue",
        "A_M",
        "B",
        "A_p_jackknife_pseudovalue",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_covariance_csv(path: Path, payload: Mapping[str, object]) -> None:
    sizes = [int(value) for value in payload["sizes"]]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "metric",
                "N_i",
                "N_j",
                "covariance_of_means",
                "correlation",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        metrics = payload["metrics"]
        assert isinstance(metrics, dict)
        for metric, raw in metrics.items():
            assert isinstance(raw, dict)
            covariance = raw["covariance_of_means"]
            correlation = raw["correlation_of_batch_values"]
            for i, n_i in enumerate(sizes):
                for j in range(i, len(sizes)):
                    writer.writerow(
                        {
                            "metric": metric,
                            "N_i": n_i,
                            "N_j": sizes[j],
                            "covariance_of_means": covariance[i][j],
                            "correlation": correlation[i][j],
                        }
                    )


def write_challenge_inputs(
    output_dir: Path,
    batch_rows: Sequence[Record],
    payload: Mapping[str, object],
    seed_label: str,
) -> None:
    sizes = [int(value) for value in payload["sizes"]]
    representatives: Dict[int, Record] = {}
    for row in batch_rows:
        representatives.setdefault(int(row["N"]), row)
    metrics = payload["metrics"]
    assert isinstance(metrics, dict)
    delta = metrics["delta_M"]
    means = delta["means"]
    standard_errors = delta["standard_errors"]
    covariance = delta["covariance_of_means"]

    observation_path = output_dir / "delta_M_observations.csv"
    with observation_path.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "row_id",
            "N",
            "seed",
            "delta_M",
            "delta_M_se",
            "delta_cos4",
            "a1",
            "b1",
            "a2",
            "b2",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for n in sizes:
            row = representatives[n]
            writer.writerow(
                {
                    "row_id": "{}:{}".format(n, seed_label),
                    "N": n,
                    "seed": seed_label,
                    "delta_M": means[str(n)],
                    "delta_M_se": standard_errors[str(n)],
                    "delta_cos4": row["delta_cos4"],
                    "a1": row["a1"],
                    "b1": row["b1"],
                    "a2": row["a2"],
                    "b2": row["b2"],
                }
            )

    covariance_path = output_dir / "delta_M_covariance.csv"
    with covariance_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["row_id_i", "row_id_j", "covariance"],
            lineterminator="\n",
        )
        writer.writeheader()
        for i, n_i in enumerate(sizes):
            for j in range(i, len(sizes)):
                writer.writerow(
                    {
                        "row_id_i": "{}:{}".format(n_i, seed_label),
                        "row_id_j": "{}:{}".format(sizes[j], seed_label),
                        "covariance": covariance[i][j],
                    }
                )


def parse_sizes(text: str) -> Tuple[int, ...]:
    try:
        values = tuple(int(value) for value in text.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "sizes must be comma-separated integers"
        ) from exc
    if (
        not values
        or any(value <= 0 for value in values)
        or len(set(values)) != len(values)
    ):
        raise argparse.ArgumentTypeError("sizes must be unique positive integers")
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histograms", type=Path, required=True)
    parser.add_argument("--p", default="0.592746050790")
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument(
        "--training-sizes", type=parse_sizes, default=(65, 85, 130)
    )
    parser.add_argument(
        "--heldout-sizes", type=parse_sizes, default=(145, 170)
    )
    parser.add_argument("--seed-label", default="threshold-rank-coupled")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.dps < 30:
        raise SystemExit("--dps must be at least 30")
    if set(args.training_sizes) & set(args.heldout_sizes):
        raise SystemExit("training and held-out sizes must be disjoint")
    mp.mp.dps = args.dps
    p = mp.mpf(args.p)
    if not 0 < p < 1:
        raise SystemExit("--p must lie strictly between zero and one")

    records = read_histograms(args.histograms)
    batch_rows, payload = audit(
        records, p, args.training_sizes, args.heldout_sizes
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_batch_csv(args.output_dir / "batch_metrics.csv", batch_rows)
    write_covariance_csv(
        args.output_dir / "cross_size_covariance.csv", payload
    )
    write_challenge_inputs(
        args.output_dir, batch_rows, payload, args.seed_label
    )
    (args.output_dir / "summary.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )

    audits = payload["constant_amplitude_audits"]
    print(
        json.dumps(
            {
                "sizes": payload["sizes"],
                "batch_count": payload["batch_count"],
                "root_gap_method": payload["nonlinear_estimator"][
                    "root_gap_method"
                ],
                "A_M_full_heldout_chi_square": audits["A_M"][
                    "full_covariance"
                ]["heldout_chi_square"],
                "A_M_diagonal_heldout_chi_square": audits["A_M"][
                    "diagonal_covariance"
                ]["heldout_chi_square"],
                "A_p_full_heldout_chi_square": audits["A_p"][
                    "full_covariance"
                ]["heldout_chi_square"],
                "A_p_diagonal_heldout_chi_square": audits["A_p"][
                    "diagonal_covariance"
                ]["heldout_chi_square"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
