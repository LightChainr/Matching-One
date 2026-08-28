#!/usr/bin/env python3
"""Validate the committed threshold-rank covariance audit archive.

This checker is deliberately independent of the simulation and reconstruction
path.  It verifies the batch/geometry contract, covariance symmetry and
positive definiteness, reported standard errors/correlations, matrix
conditioning, jackknife centers, and held-out score metadata.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple


Matrix = List[List[float]]
Record = Dict[str, str]


def _as_float(value: object, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("{} must be numeric".format(label)) from exc
    if not math.isfinite(result):
        raise ValueError("{} must be finite".format(label))
    return result


def _matrix(value: object, label: str) -> Matrix:
    if not isinstance(value, list) or not value:
        raise ValueError("{} must be a nonempty matrix".format(label))
    output: Matrix = []
    width = None
    for row_index, raw_row in enumerate(value):
        if not isinstance(raw_row, list):
            raise ValueError("{} row {} is not a list".format(label, row_index))
        row = [
            _as_float(entry, "{}[{}]".format(label, row_index))
            for entry in raw_row
        ]
        if width is None:
            width = len(row)
        if not row or len(row) != width:
            raise ValueError("{} must be nonempty and rectangular".format(label))
        output.append(row)
    return output


def _inverse(matrix: Matrix) -> Matrix:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    augmented = [
        list(map(float, matrix[row]))
        + [1.0 if row == column else 0.0 for column in range(size)]
        for row in range(size)
    ]
    scale = max(abs(value) for row in matrix for value in row)
    tolerance = max(scale * 1e-14, 1e-300)
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= tolerance:
            raise ValueError("covariance matrix is singular or numerically unresolved")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        for entry in range(2 * size):
            augmented[column][entry] /= divisor
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            for entry in range(2 * size):
                augmented[row][entry] -= factor * augmented[column][entry]
    return [row[size:] for row in augmented]


def covariance_diagnostics(
    covariance: Matrix,
    label: str,
    max_condition: float,
) -> Dict[str, float]:
    """Require a finite symmetric positive-definite covariance matrix."""

    size = len(covariance)
    if size == 0 or any(len(row) != size for row in covariance):
        raise ValueError("{} must be nonempty and square".format(label))
    if any(not math.isfinite(value) for row in covariance for value in row):
        raise ValueError("{} contains non-finite entries".format(label))
    scale = max(abs(value) for row in covariance for value in row)
    if scale <= 0.0:
        raise ValueError("{} has no positive scale".format(label))
    symmetry_error = max(
        abs(covariance[i][j] - covariance[j][i])
        for i in range(size)
        for j in range(size)
    )
    if symmetry_error > max(scale * 1e-12, 1e-300):
        raise ValueError("{} is not symmetric".format(label))

    symmetric = [
        [0.5 * (covariance[i][j] + covariance[j][i]) for j in range(size)]
        for i in range(size)
    ]
    cholesky = [[0.0] * size for _ in range(size)]
    minimum_pivot = math.inf
    for i in range(size):
        for j in range(i + 1):
            value = symmetric[i][j] - math.fsum(
                cholesky[i][k] * cholesky[j][k] for k in range(j)
            )
            if i == j:
                if value <= 0.0 or not math.isfinite(value):
                    raise ValueError("{} is not positive definite".format(label))
                minimum_pivot = min(minimum_pivot, value)
                cholesky[i][j] = math.sqrt(value)
            else:
                cholesky[i][j] = value / cholesky[j][j]

    inverse = _inverse(symmetric)
    norm = max(math.fsum(abs(value) for value in row) for row in symmetric)
    inverse_norm = max(math.fsum(abs(value) for value in row) for row in inverse)
    condition = norm * inverse_norm
    if not math.isfinite(condition) or condition > max_condition:
        raise ValueError(
            "{} condition number {:.6g} exceeds {:.6g}".format(
                label, condition, max_condition
            )
        )
    return {
        "dimension": float(size),
        "infinity_norm_condition": condition,
        "minimum_cholesky_pivot": minimum_pivot,
        "maximum_symmetry_error": symmetry_error,
    }


def read_batch_rows(path: Path) -> List[Record]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise ValueError("batch metrics file is empty")
    return rows


def validate_batch_rows(rows: Sequence[Mapping[str, str]]) -> Dict[str, object]:
    required = {
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
    }
    sample_counts = set()
    geometries: Dict[int, set] = {}
    angular_leverage: Dict[int, set] = {}
    batches_by_size: Dict[int, set] = {}
    numeric_fields = required - {"N", "batch", "samples", "a1", "b1", "a2", "b2"}

    for index, row in enumerate(rows):
        missing = required - set(row)
        if missing:
            raise ValueError(
                "batch row {} is missing {}".format(index, sorted(missing))
            )
        try:
            n = int(row["N"])
            batch = int(row["batch"])
            samples = int(row["samples"])
            a1, b1 = int(row["a1"]), int(row["b1"])
            a2, b2 = int(row["a2"]), int(row["b2"])
        except (TypeError, ValueError) as exc:
            raise ValueError("batch row {} has invalid integer fields".format(index)) from exc
        if n <= 0 or batch < 0 or samples <= 0:
            raise ValueError("batch row {} has invalid N/batch/samples".format(index))
        if a1 * a1 + b1 * b1 != n or a2 * a2 + b2 * b2 != n:
            raise ValueError("batch row {} violates the Gaussian norm contract".format(index))
        if (a1, b1) == (a2, b2):
            raise ValueError("batch row {} has identical orientations".format(index))
        for field in numeric_fields:
            _as_float(row[field], "row {} {}".format(index, field))

        sample_counts.add(samples)
        geometries.setdefault(n, set()).add((a1, b1, a2, b2))
        angular_leverage.setdefault(n, set()).add(row["delta_cos4"])
        batches_by_size.setdefault(n, set()).add(batch)

    if len(sample_counts) != 1:
        raise ValueError("all aligned size/orientation batches must use one sample count")
    for n, values in geometries.items():
        if len(values) != 1:
            raise ValueError("N={} changes geometry across aligned batches".format(n))
    for n, values in angular_leverage.items():
        if len(values) != 1:
            raise ValueError("N={} changes angular leverage across batches".format(n))
    common_batches = None
    for n, values in batches_by_size.items():
        ordered = sorted(values)
        if ordered != list(range(len(ordered))):
            raise ValueError("N={} batch ids are not contiguous".format(n))
        if common_batches is None:
            common_batches = ordered
        elif ordered != common_batches:
            raise ValueError("all sizes must have the same aligned batch ids")
    assert common_batches is not None

    return {
        "sizes": sorted(batches_by_size),
        "batch_count": len(common_batches),
        "samples_per_size_batch": next(iter(sample_counts)),
        "geometry_by_N": {
            str(n): list(next(iter(geometries[n]))) for n in sorted(geometries)
        },
    }


def _mapping_by_size(value: object, sizes: Sequence[int], label: str) -> Dict[str, float]:
    if not isinstance(value, dict):
        raise ValueError("{} must be a mapping".format(label))
    expected = {str(size) for size in sizes}
    if set(value) != expected:
        raise ValueError("{} keys do not match sizes".format(label))
    return {key: _as_float(raw, "{} {}".format(label, key)) for key, raw in value.items()}


def validate_summary(
    summary: Mapping[str, object],
    max_condition: float,
) -> Dict[str, object]:
    if summary.get("format_version") != 2:
        raise ValueError("summary format_version must be 2")
    raw_sizes = summary.get("sizes")
    if not isinstance(raw_sizes, list) or not raw_sizes:
        raise ValueError("summary sizes must be a nonempty list")
    sizes = [int(value) for value in raw_sizes]
    if len(set(sizes)) != len(sizes) or any(value <= 0 for value in sizes):
        raise ValueError("summary sizes must be unique positive integers")
    batch_count = int(summary.get("batch_count", 0))
    if batch_count < 2:
        raise ValueError("summary batch_count must be at least two")

    raw_metrics = summary.get("metrics")
    if not isinstance(raw_metrics, dict) or not raw_metrics:
        raise ValueError("summary metrics must be a nonempty mapping")
    metric_reports: Dict[str, object] = {}
    for name, raw_metric in raw_metrics.items():
        if not isinstance(raw_metric, dict):
            raise ValueError("metric {} must be a mapping".format(name))
        means = _mapping_by_size(raw_metric.get("means"), sizes, "{} means".format(name))
        standard_errors = _mapping_by_size(
            raw_metric.get("standard_errors"), sizes, "{} standard_errors".format(name)
        )
        covariance = _matrix(
            raw_metric.get("covariance_of_means"),
            "{} covariance_of_means".format(name),
        )
        if len(covariance) != len(sizes) or any(
            len(row) != len(sizes) for row in covariance
        ):
            raise ValueError("metric {} covariance dimension disagrees with sizes".format(name))
        diagnostics = covariance_diagnostics(
            covariance, "metric {} covariance".format(name), max_condition
        )
        for index, size in enumerate(sizes):
            variance = covariance[index][index]
            reported = standard_errors[str(size)]
            if variance <= 0.0 or not math.isclose(
                reported * reported, variance, rel_tol=1e-9, abs_tol=1e-300
            ):
                raise ValueError("metric {} N={} SE disagrees with covariance".format(name, size))

        correlation = _matrix(
            raw_metric.get("correlation_of_batch_values"),
            "{} correlation".format(name),
        )
        if len(correlation) != len(sizes) or any(
            len(row) != len(sizes) for row in correlation
        ):
            raise ValueError("metric {} correlation dimension disagrees with sizes".format(name))
        maximum_correlation_error = 0.0
        for i in range(len(sizes)):
            for j in range(len(sizes)):
                expected = covariance[i][j] / math.sqrt(
                    covariance[i][i] * covariance[j][j]
                )
                maximum_correlation_error = max(
                    maximum_correlation_error, abs(correlation[i][j] - expected)
                )
        if maximum_correlation_error > 1e-10:
            raise ValueError("metric {} correlation disagrees with covariance".format(name))
        diagnostics["maximum_correlation_error"] = maximum_correlation_error
        diagnostics["maximum_absolute_mean"] = max(abs(value) for value in means.values())
        metric_reports[str(name)] = diagnostics

    nonlinear = summary.get("nonlinear_estimator")
    if not isinstance(nonlinear, dict) or nonlinear.get("root_gap_method") != "delete_one_jackknife_pseudovalues":
        raise ValueError("summary must declare delete-one jackknife root-gap pseudo-values")
    by_n = nonlinear.get("by_N")
    if not isinstance(by_n, dict) or set(by_n) != {str(size) for size in sizes}:
        raise ValueError("nonlinear by_N keys do not match sizes")
    root_metric = raw_metrics.get("root_gap")
    if not isinstance(root_metric, dict) or not isinstance(root_metric.get("means"), dict):
        raise ValueError("root_gap metric is missing")
    for size in sizes:
        detail = by_n[str(size)]
        if not isinstance(detail, dict):
            raise ValueError("root detail N={} is invalid".format(size))
        if int(detail.get("delete_one_count", 0)) != batch_count:
            raise ValueError("root detail N={} delete-one count disagrees".format(size))
        center = _as_float(detail.get("bias_corrected_estimate"), "root center")
        metric_center = _as_float(root_metric["means"][str(size)], "root mean")
        if not math.isclose(center, metric_center, rel_tol=1e-13, abs_tol=1e-300):
            raise ValueError("root pseudo-value center disagrees at N={}".format(size))

    audits = summary.get("constant_amplitude_audits")
    if not isinstance(audits, dict) or set(audits) != {"A_M", "A_p"}:
        raise ValueError("summary must contain A_M and A_p constant-amplitude audits")
    audit_reports: Dict[str, object] = {}
    for metric, modes in audits.items():
        if not isinstance(modes, dict) or set(modes) != {
            "full_covariance",
            "diagonal_covariance",
        }:
            raise ValueError("{} audit modes are incomplete".format(metric))
        audit_reports[metric] = {}
        for mode, raw in modes.items():
            if not isinstance(raw, dict):
                raise ValueError("{} {} audit is invalid".format(metric, mode))
            residual_covariance = _matrix(
                raw.get("heldout_residual_covariance"),
                "{} {} heldout residual covariance".format(metric, mode),
            )
            diagnostics = covariance_diagnostics(
                residual_covariance,
                "{} {} heldout residual covariance".format(metric, mode),
                max_condition,
            )
            chi_square = _as_float(raw.get("heldout_chi_square"), "heldout chi-square")
            dof = int(raw.get("heldout_dof", 0))
            if chi_square < 0.0 or dof != len(residual_covariance):
                raise ValueError("{} {} heldout score metadata is invalid".format(metric, mode))
            diagnostics["heldout_chi_square"] = chi_square
            diagnostics["heldout_dof"] = dof
            audit_reports[metric][mode] = diagnostics

    return {
        "format_version": 1,
        "sizes": sizes,
        "batch_count": batch_count,
        "max_condition_allowed": max_condition,
        "metric_covariance_diagnostics": metric_reports,
        "heldout_residual_diagnostics": audit_reports,
    }


def validate_archive(
    batch_metrics: Path,
    summary_path: Path,
    max_condition: float,
) -> Dict[str, object]:
    batch_report = validate_batch_rows(read_batch_rows(batch_metrics))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if not isinstance(summary, dict):
        raise ValueError("summary JSON must contain an object")
    summary_report = validate_summary(summary, max_condition)
    if batch_report["sizes"] != summary_report["sizes"]:
        raise ValueError("batch and summary sizes disagree")
    if batch_report["batch_count"] != summary_report["batch_count"]:
        raise ValueError("batch and summary batch counts disagree")
    return {
        "format_version": 1,
        "batch_contract": batch_report,
        "summary_contract": summary_report,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-metrics", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--max-condition", type=float, default=1e12)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.max_condition) or args.max_condition <= 1.0:
        raise SystemExit("--max-condition must be finite and greater than one")
    try:
        result = validate_archive(
            args.batch_metrics, args.summary, args.max_condition
        )
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        raise SystemExit(str(exc)) from exc
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
