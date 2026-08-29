#!/usr/bin/env python3
"""Frozen scorer views for the Issue #225 norm-5 multi-radius pilot."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


DESIGNS = {
    "n325_first": (325, 17, 6),
    "n325_second": (325, 18, 1),
    "n425_first": (425, 16, 13),
    "n425_second": (425, 19, 8),
}
PAIRS = {
    325: ("n325_first", "n325_second"),
    425: ("n425_first", "n425_second"),
}
RADII = [2, 4, 7, 8]
INTEGER_FIELDS = (
    "n", "a", "b", "radius", "batch", "counter_first",
    "counter_last_exclusive", "samples", "primal_pivotal",
    "matching_pivotal", "primal_h4", "matching_h4", "h4_plus", "h4_minus",
)
SUM_FIELDS = (
    "samples", "primal_pivotal", "matching_pivotal", "primal_h4",
    "matching_h4", "h4_plus", "h4_minus",
)


def read_rows(path: Path):
    rows = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not {"label", *INTEGER_FIELDS}.issubset(reader.fieldnames):
            raise ValueError("batch schema mismatch")
        for raw in reader:
            label = raw["label"]
            if label not in DESIGNS:
                raise ValueError(f"unexpected design {label}")
            row = {field: int(raw[field]) for field in INTEGER_FIELDS}
            if (row["n"], row["a"], row["b"]) != DESIGNS[label]:
                raise ValueError(f"design arithmetic mismatch for {label}")
            if row["h4_plus"] != row["primal_h4"] + row["matching_h4"]:
                raise ValueError("plus identity failed")
            if row["h4_minus"] != row["primal_h4"] - row["matching_h4"]:
                raise ValueError("minus identity failed")
            key = label, row["radius"], row["batch"]
            if key in rows:
                raise ValueError(f"duplicate row {key}")
            rows[key] = row
    batches = sorted({key[2] for key in rows})
    expected = {(label, radius, batch) for label in DESIGNS for radius in RADII for batch in batches}
    if set(rows) != expected:
        raise ValueError("design/radius/batch grid mismatch")
    for label in DESIGNS:
        for batch in batches:
            first = rows[label, RADII[0], batch]
            for radius in RADII[1:]:
                row = rows[label, radius, batch]
                for field in (
                    "samples", "primal_pivotal", "matching_pivotal",
                    "counter_first", "counter_last_exclusive",
                ):
                    if row[field] != first[field]:
                        raise ValueError(f"cross-radius alignment failed for {label}, batch {batch}")
    return rows, batches


def aggregate(rows, label, radius, batches, omit=None):
    total = {field: 0 for field in SUM_FIELDS}
    for batch in batches:
        if batch == omit:
            continue
        row = rows[label, radius, batch]
        for field in SUM_FIELDS:
            total[field] += row[field]
    pivotal = total["primal_pivotal"] + total["matching_pivotal"]
    if pivotal == 0:
        raise ValueError(f"zero pivotal denominator for {label},R={radius}")
    n = DESIGNS[label][0]
    return {
        "mu0": n * pivotal / total["samples"],
        "A_plus": total["h4_plus"] / pivotal,
        "A_minus": total["h4_minus"] / pivotal,
        "pivotal_events": pivotal,
        "samples": total["samples"],
    }


def contrast_vectors(rows, batches, omit=None):
    points = {
        (label, radius): aggregate(rows, label, radius, batches, omit)
        for label in DESIGNS for radius in RADII
    }
    order = []
    vector = []
    contrasts = {}
    for n, (first, second) in PAIRS.items():
        for radius in RADII:
            contrasts[n, radius] = {}
            for channel in ("A_plus", "A_minus"):
                value = points[first, radius][channel] - points[second, radius][channel]
                contrasts[n, radius][channel] = value
                order.append(f"N{n}_R{radius}_Delta_{channel}")
                vector.append(value)
    shell_order = []
    shell_vector = []
    shells = []
    for n in PAIRS:
        for first_radius, second_radius in ((2, 4), (4, 8)):
            record = {"N": n, "R_first": first_radius, "R_second": second_radius}
            for channel in ("A_plus", "A_minus"):
                value = (contrasts[n, second_radius][channel] -
                         contrasts[n, first_radius][channel]) / math.log(2)
                record[f"per_log_{channel}"] = value
                shell_order.append(
                    f"N{n}_R{first_radius}_to_R{second_radius}_per_log_Delta_{channel}")
                shell_vector.append(value)
            shells.append(record)
    matched_order = []
    matched_vector = []
    for n, radius in ((325, 7), (425, 8)):
        for channel in ("A_plus", "A_minus"):
            matched_order.append(f"N{n}_R{radius}_Delta_{channel}")
            matched_vector.append(contrasts[n, radius][channel])
    return points, contrasts, order, vector, shells, shell_order, shell_vector, matched_order, matched_vector


def covariance(replicates):
    count = len(replicates)
    width = len(replicates[0])
    mean = [sum(row[i] for row in replicates) / count for i in range(width)]
    return [[
        (count - 1) / count * sum(
            (row[i] - mean[i]) * (row[j] - mean[j]) for row in replicates)
        for j in range(width)] for i in range(width)]


def solve(matrix, rhs):
    """Small dense partial-pivot solve for the frozen covariance scores."""
    size = len(rhs)
    augmented = [list(map(float, matrix[row])) + [float(rhs[row])] for row in range(size)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1e-18:
            raise ArithmeticError("singular covariance in GLS score")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                augmented[row][index] - factor * augmented[column][index]
                for index in range(size + 1)
            ]
    return [augmented[row][-1] for row in range(size)]


def chi_square_survival(value, degrees_of_freedom):
    if degrees_of_freedom == 2:
        return math.exp(-value / 2)
    if degrees_of_freedom == 3:
        return (math.erfc(math.sqrt(value / 2)) +
                math.sqrt(2 * value / math.pi) * math.exp(-value / 2))
    raise ValueError("only the frozen df=2/3 scores are supported")


def constant_gls(values, matrix):
    precision_one = solve(matrix, [1.0] * len(values))
    denominator = sum(precision_one)
    mean = sum(weight * value for weight, value in zip(precision_one, values)) / denominator
    residual = [value - mean for value in values]
    precision_residual = solve(matrix, residual)
    chi_square = sum(left * right for left, right in zip(residual, precision_residual))
    degrees_of_freedom = len(values) - 1
    return {
        "constant_GLS_amplitude": mean,
        "constant_GLS_SE": math.sqrt(1 / denominator),
        "chi_square": chi_square,
        "degrees_of_freedom": degrees_of_freedom,
        "chi_square_survival": chi_square_survival(chi_square, degrees_of_freedom),
    }


def subcovariance(full_order, full_cov, selected):
    indices = [full_order.index(label) for label in selected]
    return [[full_cov[i][j] for j in indices] for i in indices]


def analyze(batch_path: Path, metadata_path: Path, require_production: bool = True):
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("schema") != "matching-one/matching-multiradius-pivotal/v1":
        raise ValueError("metadata schema mismatch")
    if metadata.get("cutoff") != "euclidean" or metadata.get("radii") != RADII:
        raise ValueError("cutoff/radius contract mismatch")
    if [(row["label"], row["N"], row["a"], row["b"]) for row in metadata["designs"]] != [
        (label, *DESIGNS[label]) for label in DESIGNS
    ]:
        raise ValueError("metadata design order mismatch")
    if require_production and (metadata.get("samples_per_design") != 200000 or metadata.get("batches") != 200):
        raise ValueError("production scorer requires 200k samples and 200 batches")
    rows, batches = read_rows(batch_path)
    points, contrasts, order, vector, shells, shell_order, shell_vector, matched_order, matched_vector = contrast_vectors(rows, batches)
    delete_vectors = []
    delete_shells = []
    for omitted in batches:
        result = contrast_vectors(rows, batches, omitted)
        delete_vectors.append(result[3])
        delete_shells.append(result[6])
    full_cov = covariance(delete_vectors)
    shell_cov = covariance(delete_shells)
    same_r_rows = []
    for n in PAIRS:
        for radius in (2, 4, 8):
            row = {"N": n, "R": radius, "delta": radius / math.sqrt(n)}
            for channel in ("A_plus", "A_minus"):
                label = f"N{n}_R{radius}_Delta_{channel}"
                index = order.index(label)
                row[f"Delta_{channel}"] = vector[index]
                row[f"Delta_{channel}_SE"] = math.sqrt(max(full_cov[index][index], 0.0))
            same_r_rows.append(row)
    for row in shells:
        for channel in ("A_plus", "A_minus"):
            label = f"N{row['N']}_R{row['R_first']}_to_R{row['R_second']}_per_log_Delta_{channel}"
            index = shell_order.index(label)
            row[f"per_log_{channel}_SE"] = math.sqrt(max(shell_cov[index][index], 0.0))
    matched_cov = subcovariance(order, full_cov, matched_order)
    transform = [[1, 0, -1, 0], [0, 1, 0, -1]]
    matched_difference = [matched_vector[0] - matched_vector[2], matched_vector[1] - matched_vector[3]]
    matched_difference_cov = [[
        sum(transform[i][a] * matched_cov[a][b] * transform[j][b]
            for a in range(4) for b in range(4))
        for j in range(2)] for i in range(2)]
    shell_scores = {}
    for channel in ("A_plus", "A_minus"):
        indices = [
            index for index, label in enumerate(shell_order)
            if label.endswith(channel)
        ]
        values = [shell_vector[index] for index in indices]
        matrix = [[shell_cov[first][second] for second in indices] for first in indices]
        shell_scores[channel] = constant_gls(values, matrix)
    matched_precision_residual = solve(matched_difference_cov, matched_difference)
    matched_chi_square = sum(
        left * right for left, right in zip(matched_difference, matched_precision_residual))
    matched_score = {
        "null": "N325_R7 and N425_R8 orientation-contrast vectors are equal at matched delta",
        "chi_square": matched_chi_square,
        "degrees_of_freedom": 2,
        "chi_square_survival": chi_square_survival(matched_chi_square, 2),
    }
    return {
        "schema": "matching-one/norm5-multiradius-pivotal-score/v1",
        "source": {"batches": str(batch_path), "metadata": str(metadata_path)},
        "orientation_registry": {
            "order": "first_minus_second",
            "N325": {"first": [17, 6], "second": [18, 1], "Delta_cos4": "-16128/21125"},
            "N425": {"first": [16, 13], "second": [19, 8], "Delta_cos4": "-32256/36125"},
            "sign_reversal_between_sizes": False,
        },
        "same_R_UV": same_r_rows,
        "dyadic_shells": shells,
        "matched_delta": {
            "rows": [{"N": 325, "R": 7, "delta": 7 / math.sqrt(325)},
                     {"N": 425, "R": 8, "delta": 8 / math.sqrt(425)}],
            "relative_delta_mismatch": abs(7 / math.sqrt(325) - 8 / math.sqrt(425)) /
                                       ((7 / math.sqrt(325) + 8 / math.sqrt(425)) / 2),
            "order": matched_order,
            "point": matched_vector,
            "covariance": matched_cov,
            "N325_minus_N425": matched_difference,
            "N325_minus_N425_covariance": matched_difference_cov,
        },
        "contrast_vector": {"order": order, "point": vector, "covariance": full_cov},
        "shell_vector": {"order": shell_order, "point": shell_vector, "covariance": shell_cov},
        "hypothesis_scores": {
            "constant_shell_GLS": shell_scores,
            "matched_delta_equality": matched_score,
        },
        "interpretation": {
            "same_R": "UV comparison only; delta changes with N",
            "matched_delta": "N325,R7 versus N425,R8 holds delta to about 0.1 percent",
            "evidence_rule": "all radii/channels/sizes are one correlated block",
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batches", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-smoke", action="store_true")
    args = parser.parse_args()
    result = analyze(args.batches, args.metadata, not args.allow_smoke)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
