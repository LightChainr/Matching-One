#!/usr/bin/env python3
"""Frozen Issue #154 Phase-E local-singlet pilot score.

The protocol, candidates and promotion gate were frozen in commit 4f249dd.
This scorer consumes only the two locked batch tables.  It reports every
candidate; it does not select a row by its p value.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable


FIELDS = ("A", "E", "C", "W", "J_black", "J_white", "J_even", "J_odd")
CANDIDATES = (("A", "E", "C"), ("A", "E", "J_even"), ("A", "E", "J_odd"))
EXPECTED = {
    65: {
        "first": (8, 1), "second": (7, 4), "samples": 20000, "batches": 100,
        "seed": 202615465, "counter": (15465000000, 15465020000),
    },
    130: {
        "first": (11, 3), "second": (9, 7), "samples": 20000, "batches": 100,
        "seed": 2026154130, "counter": (15465200000, 15465220000),
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def cos4(a: int, b: int) -> float:
    norm = a * a + b * b
    return (a**4 - 6 * a * a * b * b + b**4) / (norm * norm)


def row_values(row: dict[str, str]) -> dict[str, float]:
    n = int(row["n"])
    samples = int(row["samples"])
    k1 = int(row["sum_k1"]) / samples
    k2 = int(row["sum_k2"]) / samples
    i0 = int(row["sum_i0"]) / samples
    i2 = int(row["sum_i2"]) / samples
    black = int(row["sum_black_axis_pairs"]) / (samples * 2 * n)
    white = int(row["sum_white_matching_axis_pairs"]) / (samples * 2 * n)
    return {
        "A": i2 - i0,
        "E": i2 + i0,
        "C": (k1 + k2) / (2 * (n + 1)) - 0.5,
        "W": (k2 - k1) / (n + 1),
        "J_black": black,
        "J_white": white,
        "J_even": 0.5 * (black + white),
        "J_odd": 0.5 * (black - white),
    }


def read_projected(path: Path, metadata_path: Path) -> tuple[int, list[list[float]], dict]:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    n = int(metadata["N"])
    expected = EXPECTED[n]
    if metadata["git_commit"] != "4f249dd74f8efb6a6175addd44d97b876c411c8f":
        raise ValueError(f"{path}: runner commit differs from freeze")
    for key, wanted in (("samples", expected["samples"]), ("batches", expected["batches"]),
                        ("seed", expected["seed"]),
                        ("replica_counter_first", expected["counter"][0]),
                        ("replica_counter_last_exclusive", expected["counter"][1])):
        if int(metadata[key]) != wanted:
            raise ValueError(f"{path}: {key} violates freeze")
    if not math.isclose(float(metadata["p_ref"]), 0.59274605079, rel_tol=0, abs_tol=1e-15):
        raise ValueError(f"{path}: p_ref violates freeze")

    rows: dict[tuple[int, str], dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["n"]) != n:
                raise ValueError("mixed N in batch file")
            rows[(int(row["batch"]), row["orientation"])] = row
    if len(rows) != 2 * expected["batches"]:
        raise ValueError("batch/orientation row count violates freeze")
    first = expected["first"]
    second = expected["second"]
    delta = cos4(*first) - cos4(*second)
    vectors = []
    for batch in range(expected["batches"]):
        one = rows[(batch, "first")]
        two = rows[(batch, "second")]
        if (int(one["a"]), int(one["b"])) != first or (int(two["a"]), int(two["b"])) != second:
            raise ValueError("orientation label violates freeze")
        if int(one["samples"]) != expected["samples"] // expected["batches"] or int(two["samples"]) != expected["samples"] // expected["batches"]:
            raise ValueError("per-batch samples violate freeze")
        x, y = row_values(one), row_values(two)
        vectors.append([(x[field] - y[field]) / delta for field in FIELDS])
    return n, vectors, metadata


def mean_cov(vectors: list[list[float]]) -> tuple[list[float], list[list[float]]]:
    count, dimension = len(vectors), len(vectors[0])
    mean = [sum(row[j] for row in vectors) / count for j in range(dimension)]
    covariance = [[0.0] * dimension for _ in range(dimension)]
    for j in range(dimension):
        for k in range(dimension):
            covariance[j][k] = sum(
                (row[j] - mean[j]) * (row[k] - mean[k]) for row in vectors
            ) / ((count - 1) * count)
    return mean, covariance


def solve(matrix: list[list[float]], vector: list[float]) -> list[float]:
    n = len(vector)
    augmented = [list(matrix[i]) + [vector[i]] for i in range(n)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1e-24:
            raise ValueError("singular frozen covariance block")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [a - factor * b for a, b in zip(augmented[row], augmented[column])]
    return [augmented[i][-1] for i in range(n)]


def quadratic(vector: list[float], covariance: list[list[float]]) -> float:
    inverse_times = solve(covariance, vector)
    return sum(a * b for a, b in zip(vector, inverse_times))


def subvector(mean: list[float], covariance: list[list[float]], fields: Iterable[str]):
    indices = [FIELDS.index(field) for field in fields]
    return ([mean[i] for i in indices], [[covariance[i][j] for j in indices] for i in indices])


def ray_objective(scale: float, first, second) -> float:
    first_mean, first_covariance = first
    second_mean, second_covariance = second
    residual = [y - scale * x for x, y in zip(first_mean, second_mean)]
    covariance = [[second_covariance[i][j] + scale * scale * first_covariance[i][j]
                   for j in range(len(residual))] for i in range(len(residual))]
    return quadratic(residual, covariance)


def fit_ray(first, second) -> tuple[float, float]:
    grid = [(-20.0 + 0.05 * index) for index in range(801)]
    values = [ray_objective(scale, first, second) for scale in grid]
    best = min(range(len(grid)), key=values.__getitem__)
    left = grid[max(0, best - 1)]
    right = grid[min(len(grid) - 1, best + 1)]
    golden = (math.sqrt(5.0) - 1.0) / 2.0
    x1 = right - golden * (right - left)
    x2 = left + golden * (right - left)
    f1, f2 = ray_objective(x1, first, second), ray_objective(x2, first, second)
    for _ in range(80):
        if f1 <= f2:
            right, x2, f2 = x2, x1, f1
            x1 = right - golden * (right - left)
            f1 = ray_objective(x1, first, second)
        else:
            left, x1, f1 = x1, x2, f2
            x2 = left + golden * (right - left)
            f2 = ray_objective(x2, first, second)
    scale = 0.5 * (left + right)
    return scale, ray_objective(scale, first, second)


def calculate(raw_dir: Path) -> dict:
    points = {}
    inputs = []
    for n in (65, 130):
        csv_path = raw_dir / f"n{n}_pilot.batches.csv"
        metadata_path = raw_dir / f"n{n}_pilot.metadata.json"
        actual_n, vectors, metadata = read_projected(csv_path, metadata_path)
        mean, covariance = mean_cov(vectors)
        points[actual_n] = (mean, covariance)
        inputs += [{"path": str(csv_path), "sha256": sha256(csv_path)},
                   {"path": str(metadata_path), "sha256": sha256(metadata_path)}]

    estimates = {}
    for n, (mean, covariance) in points.items():
        estimates[str(n)] = {
            field: {"estimate": mean[index], "se": math.sqrt(max(covariance[index][index], 0.0)),
                    "z": mean[index] / math.sqrt(covariance[index][index])}
            for index, field in enumerate(FIELDS)
        }
        estimates[str(n)]["covariance"] = covariance

    rays = {}
    for fields in CANDIDATES:
        first = subvector(*points[65], fields)
        second = subvector(*points[130], fields)
        scale, chi2 = fit_ray(first, second)
        rays["/".join(fields)] = {
            "fields": list(fields), "scale_N130_over_N65": scale,
            "chi_square": chi2, "df": 2, "p_value": math.exp(-0.5 * chi2),
        }

    baseline = rays["A/E/C"]["chi_square"]
    qualifying = []
    for row in ("J_even", "J_odd"):
        z_min = min(abs(estimates[str(n)][row]["z"]) for n in (65, 130))
        improvement = baseline - rays[f"A/E/{row}"]["chi_square"]
        if z_min >= 3.0 and improvement >= 4.0:
            qualifying.append(row)
    return {
        "issue": 154,
        "status": "pilot scored once against commit-4f249dd frozen candidates",
        "inputs": inputs,
        "projected_fields": list(FIELDS),
        "estimates": estimates,
        "common_ray_candidates": rays,
        "promotion_gate": {
            "rule": "same fixed local row has |z|>=3 at both N and common-ray chi2 improves by >=4 versus A/E/C",
            "qualifying_rows": qualifying,
            "decision": "extend_both_to_100k" if qualifying else "stop_at_20k",
        },
        "interpretation_boundary": "finite-lattice local connectivity singlet; not a continuum-energy identification",
    }


def report(payload: dict) -> str:
    lines = [
        "# Issue #154 Phase-E local-singlet pilot",
        "",
        "The frozen 20k same-stream pilot completed on N65 and N130. The new rows are local",
        "connectivity observables not determined by the occupation count K; they are not named as",
        "continuum energy fields.",
        "",
        "| N | P4[J_black] | P4[J_white] | P4[J_even] | P4[J_odd] |",
        "|---:|---:|---:|---:|---:|",
    ]
    for n in (65, 130):
        row = payload["estimates"][str(n)]
        def cell(field: str) -> str:
            return f"{row[field]['estimate']:.7g} +/- {row[field]['se']:.2g} (z={row[field]['z']:.2f})"
        lines.append(f"| {n} | {cell('J_black')} | {cell('J_white')} | {cell('J_even')} | {cell('J_odd')} |")
    lines += ["", "## Frozen common-ray comparison", "",
              "| candidate | scale N130/N65 | chi2 / 2 df | p |",
              "|---|---:|---:|---:|"]
    for name, row in payload["common_ray_candidates"].items():
        lines.append(f"| {name} | {row['scale_N130_over_N65']:.5g} | {row['chi_square']:.4g} | {row['p_value']:.3g} |")
    gate = payload["promotion_gate"]
    lines += ["", "## Frozen gate", "",
              f"Decision: **{gate['decision']}**. Qualifying fixed rows: `{gate['qualifying_rows']}`.",
              "Neither fixed local row resolves at both sizes, and neither improves the A/E/C",
              "common-ray score. This pilot therefore gives no production evidence that the",
              "radius-1 connectivity row replaces C or pins the E_top plane; the frozen extension",
              "is not run.",
              "All three candidates are reported; none was selected by its p value. N65 and N130",
              "are one parent-child lineage and constitute one scale comparison, not independent",
              "geometric votes. Full within-size batch covariance is retained in `score.json`."]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    payload = calculate(args.raw_dir)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.report.write_text(report(payload), encoding="utf-8")


if __name__ == "__main__":
    main()
