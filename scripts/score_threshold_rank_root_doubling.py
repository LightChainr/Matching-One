#!/usr/bin/env python3
"""Score the no-fit Gaussian root-gap doubling relation from covariance audit output.

Input is ``summary.json`` from ``audit_threshold_rank_covariance.py``.  The
root-gap means and covariance there are based on delete-one jackknife
pseudo-values.  This script preserves multiplication-by-``1+i`` lineage order,
which reverses the repository's display order at N=130 and N=170, and tests

    Delta p*(2N) / Delta p*(N) = -1/4.

Both the measured full cross-size covariance and a diagonal-only diagnostic are
reported.  The diagonal result is a comparison, not a replacement.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple


Matrix = List[List[float]]
Vector = List[float]
LINEAGES = ((65, 130), (85, 170))
TARGET_RATIO = -0.25


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
            raise ArithmeticError("singular or ill-resolved residual covariance")
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


def _matvec(matrix: Matrix, vector: Vector) -> Vector:
    return [math.fsum(a * b for a, b in zip(row, vector)) for row in matrix]


def _quadratic(vector: Vector, covariance: Matrix) -> float:
    solution = _solve(covariance, vector)
    return math.fsum(a * b for a, b in zip(vector, solution))


def _diagonal(covariance: Matrix) -> Matrix:
    return [
        [covariance[i][i] if i == j else 0.0 for j in range(len(covariance))]
        for i in range(len(covariance))
    ]


def _transform_covariance(covariance: Matrix, transform: Matrix) -> Matrix:
    width = len(covariance)
    if any(len(row) != width for row in covariance):
        raise ValueError("root covariance must be square")
    if any(len(row) != width for row in transform):
        raise ValueError("lineage transform has the wrong width")
    return [
        [
            math.fsum(
                transform[i][a] * covariance[a][b] * transform[j][b]
                for a in range(width)
                for b in range(width)
            )
            for j in range(len(transform))
        ]
        for i in range(len(transform))
    ]


def _score_covariance(
    values: Vector,
    covariance: Matrix,
    sizes: Sequence[int],
) -> Dict[str, object]:
    index = {size: position for position, size in enumerate(sizes)}
    missing = sorted({size for lineage in LINEAGES for size in lineage} - set(index))
    if missing:
        raise ValueError("missing Gaussian doubling sizes: {}".format(missing))

    transform: Matrix = [[0.0] * len(sizes) for _ in LINEAGES]
    for row, (parent, child) in enumerate(LINEAGES):
        # child_lineage = -child_stored, so
        # residual = child_lineage - TARGET_RATIO*parent
        #          = -child_stored + 0.25*parent.
        transform[row][index[parent]] = -TARGET_RATIO
        transform[row][index[child]] = -1.0

    residuals = _matvec(transform, values)
    residual_covariance = _transform_covariance(covariance, transform)
    chi_square = _quadratic(residuals, residual_covariance)

    lineages = []
    for row, (parent, child) in enumerate(LINEAGES):
        p_index = index[parent]
        c_index = index[child]
        parent_gap = values[p_index]
        stored_child_gap = values[c_index]
        lineage_child_gap = -stored_child_gap
        if parent_gap == 0.0:
            raise ValueError("zero parent root gap prevents a ratio")
        ratio = lineage_child_gap / parent_gap
        d_parent = stored_child_gap / (parent_gap * parent_gap)
        d_child = -1.0 / parent_gap
        ratio_variance = (
            d_parent * d_parent * covariance[p_index][p_index]
            + d_child * d_child * covariance[c_index][c_index]
            + 2.0 * d_parent * d_child * covariance[p_index][c_index]
        )
        residual_se = math.sqrt(residual_covariance[row][row])
        lineages.append(
            {
                "parent_N": parent,
                "child_N": child,
                "parent_stored_root_gap": parent_gap,
                "child_stored_root_gap": stored_child_gap,
                "child_lineage_root_gap": lineage_child_gap,
                "observed_ratio": ratio,
                "ratio_se_delta_method": math.sqrt(max(0.0, ratio_variance)),
                "target_ratio": TARGET_RATIO,
                "fixed_prediction_residual": residuals[row],
                "residual_se": residual_se,
                "residual_z": residuals[row] / residual_se,
            }
        )

    return {
        "target_ratio": TARGET_RATIO,
        "lineage_order_note": (
            "N=130 and N=170 stored orientation order is reversed from 1+i genealogy"
        ),
        "lineages": lineages,
        "residual_covariance": residual_covariance,
        "joint_residual_chi_square": chi_square,
        "joint_degrees_of_freedom": 2,
    }


def score(summary: Mapping[str, object]) -> Dict[str, object]:
    sizes = [int(value) for value in summary["sizes"]]
    metrics = summary["metrics"]
    if not isinstance(metrics, dict) or "root_gap" not in metrics:
        raise ValueError("summary does not contain root-gap covariance")
    root = metrics["root_gap"]
    if not isinstance(root, dict):
        raise ValueError("invalid root-gap metric payload")
    means_map = root["means"]
    covariance = root["covariance_of_means"]
    if not isinstance(means_map, dict) or not isinstance(covariance, list):
        raise ValueError("invalid root-gap means or covariance")
    values = [float(means_map[str(size)]) for size in sizes]
    covariance = [list(map(float, row)) for row in covariance]
    if len(covariance) != len(sizes) or any(len(row) != len(sizes) for row in covariance):
        raise ValueError("root-gap covariance dimension does not match sizes")

    return {
        "format_version": 1,
        "source_summary_format_version": summary.get("format_version"),
        "root_gap_method": summary.get("nonlinear_estimator", {}).get(
            "root_gap_method"
        ) if isinstance(summary.get("nonlinear_estimator"), dict) else None,
        "full_cross_size_covariance": _score_covariance(values, covariance, sizes),
        "diagonal_cross_size_covariance": _score_covariance(
            values, _diagonal(covariance), sizes
        ),
    }


def write_outputs(
    result: Mapping[str, object],
    json_path: Path,
    csv_path: Path,
    report_path: Path,
) -> None:
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    rows = []
    for covariance_mode in (
        "full_cross_size_covariance",
        "diagonal_cross_size_covariance",
    ):
        audit = result[covariance_mode]
        assert isinstance(audit, dict)
        for row in audit["lineages"]:
            rows.append({"covariance_mode": covariance_mode, **row})
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    full = result["full_cross_size_covariance"]
    diagonal = result["diagonal_cross_size_covariance"]
    assert isinstance(full, dict) and isinstance(diagonal, dict)
    lines = [
        "# Threshold-rank Gaussian root-gap doubling audit",
        "",
        "Frozen lineage prediction: `Delta p*(2N) / Delta p*(N) = -1/4`.",
        "Child signs follow multiplication by `1+i`, not stored display order.",
        "",
        "| lineage | observed ratio | ratio SE | full-cov residual z | diagonal residual z |",
        "|---|---:|---:|---:|---:|",
    ]
    for full_row, diagonal_row in zip(full["lineages"], diagonal["lineages"]):
        lines.append(
            "| {}->{} | {:.8g} | {:.3g} | {:.3f} | {:.3f} |".format(
                full_row["parent_N"],
                full_row["child_N"],
                full_row["observed_ratio"],
                full_row["ratio_se_delta_method"],
                full_row["residual_z"],
                diagonal_row["residual_z"],
            )
        )
    lines.extend(
        [
            "",
            "Full-covariance joint residual chi-square: **{:.6g} / 2**.".format(
                full["joint_residual_chi_square"]
            ),
            "Diagonal diagnostic chi-square: **{:.6g} / 2**.".format(
                diagonal["joint_residual_chi_square"]
            ),
            "",
        ]
    )
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    result = score(summary)
    write_outputs(result, args.json, args.csv, args.report)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
