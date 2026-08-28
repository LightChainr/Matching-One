#!/usr/bin/env python3
"""Score the parameter-free Gaussian N -> 2N spin-4 doubling test."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


LINEAGES = ((65, 130), (85, 170))
EXPECTED = {
    65: ((8, 1), (7, 4), 1.0),
    85: ((9, 2), (7, 6), 1.0),
    130: ((11, 3), (9, 7), -1.0),
    170: ((13, 1), (11, 7), -1.0),
}
FACTOR = 2.0 ** (-13.0 / 8.0)


def mean_se(values: list[float]) -> tuple[float, float]:
    if len(values) < 2:
        raise ValueError("at least two batches are required")
    return statistics.fmean(values), statistics.stdev(values) / math.sqrt(len(values))


def covariance_of_means(xs: list[float], ys: list[float]) -> float:
    if len(xs) != len(ys) or len(xs) < 2:
        raise ValueError("paired covariance requires equal nontrivial samples")
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    return math.fsum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (
        len(xs) * (len(xs) - 1)
    )


def load_cross_batches(paths: list[Path]) -> dict[int, dict[int, float]]:
    grouped: dict[int, dict[int, float]] = {}
    signatures: dict[int, tuple[tuple[int, int], tuple[int, int]]] = {}
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row["channel"] != "cross":
                    continue
                n = int(row["n"])
                if n not in EXPECTED:
                    continue
                first = (int(row["a1"]), int(row["b1"]))
                second = (int(row["a2"]), int(row["b2"]))
                expected_first, expected_second, lineage_sign = EXPECTED[n]
                if (first, second) != (expected_first, expected_second):
                    raise ValueError(
                        f"N={n} orientation order {(first, second)} != "
                        f"{(expected_first, expected_second)}"
                    )
                signatures[n] = (first, second)
                batch = int(row["batch"])
                samples = int(row["samples"])
                if samples <= 0 or batch in grouped.setdefault(n, {}):
                    raise ValueError(f"invalid or duplicate N={n}, batch={batch}")
                stored_delta = (
                    int(row["first_primal_sum"])
                    - int(row["first_matching_sum"])
                    - int(row["second_primal_sum"])
                    + int(row["second_matching_sum"])
                ) / samples
                grouped[n][batch] = lineage_sign * stored_delta
    missing = sorted(set(EXPECTED) - set(grouped))
    if missing:
        raise ValueError(f"missing required sizes: {missing}")
    batch_sets = {n: set(values) for n, values in grouped.items()}
    reference = batch_sets[65]
    for n, batches in batch_sets.items():
        if batches != reference:
            raise ValueError(f"N={n} batch ids are not aligned")
    return grouped


def score(grouped: dict[int, dict[int, float]]) -> dict[str, object]:
    batches = sorted(grouped[65])
    lineages: list[dict[str, object]] = []
    residual_vectors: list[list[float]] = []
    for parent_n, child_n in LINEAGES:
        parent = [grouped[parent_n][batch] for batch in batches]
        child = [grouped[child_n][batch] for batch in batches]
        residual = [y + FACTOR * x for x, y in zip(parent, child)]
        parent_mean, parent_se = mean_se(parent)
        child_mean, child_se = mean_se(child)
        residual_mean, residual_se = mean_se(residual)
        cross_cov = covariance_of_means(parent, child)
        ratio = child_mean / parent_mean
        ratio_variance = (
            child_se**2 / parent_mean**2
            + child_mean**2 * parent_se**2 / parent_mean**4
            - 2.0 * child_mean * cross_cov / parent_mean**3
        )
        lineages.append(
            {
                "parent_N": parent_n,
                "child_N": child_n,
                "parent_delta_M": parent_mean,
                "parent_se": parent_se,
                "child_lineage_delta_M": child_mean,
                "child_se": child_se,
                "ratio": ratio,
                "ratio_se_delta_method": math.sqrt(max(0.0, ratio_variance)),
                "target_ratio": -FACTOR,
                "fixed_prediction_residual": residual_mean,
                "residual_se": residual_se,
                "residual_z": residual_mean / residual_se,
                "parent_child_batch_correlation": statistics.correlation(parent, child),
            }
        )
        residual_vectors.append(residual)

    r1, r2 = residual_vectors
    m1, _ = mean_se(r1)
    m2, _ = mean_se(r2)
    v1 = covariance_of_means(r1, r1)
    v2 = covariance_of_means(r2, r2)
    c12 = covariance_of_means(r1, r2)
    determinant = v1 * v2 - c12 * c12
    if determinant <= 0.0:
        raise ValueError("residual covariance is not positive definite")
    chi2 = (v2 * m1 * m1 - 2.0 * c12 * m1 * m2 + v1 * m2 * m2) / determinant
    return {
        "schema": "Gaussian doubling spin-4 test v1",
        "target_ratio": -FACTOR,
        "batches": len(batches),
        "lineages": lineages,
        "joint_residual_chi2": chi2,
        "joint_degrees_of_freedom": 2,
        "residual_mean_covariance": [[v1, c12], [c12, v2]],
    }


def write_outputs(result: dict[str, object], json_path: Path, csv_path: Path, report: Path) -> None:
    json_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    rows = result["lineages"]
    assert isinstance(rows, list)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    lines = [
        "# Parameter-free Gaussian doubling test",
        "",
        f"Frozen prediction: `Delta M(2N) / Delta M(N) = {result['target_ratio']:.15g}`.",
        "Lineage order is preserved after multiplication by `1+i`; this reverses",
        "the stored display order at N=130 and N=170.",
        "",
        "| lineage | ratio | ratio SE | fixed residual | residual SE | z |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['parent_N']}->{row['child_N']} | {row['ratio']:.8g} | "
            f"{row['ratio_se_delta_method']:.3g} | {row['fixed_prediction_residual']:.8g} | "
            f"{row['residual_se']:.3g} | {row['residual_z']:.3f} |"
        )
    lines.extend(
        [
            "",
            f"Joint fixed-prediction residual chi-square: **{result['joint_residual_chi2']:.4g} / 2 dof**.",
            "",
        ]
    )
    report.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batches", nargs="+", type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    result = score(load_cross_batches(args.batches))
    write_outputs(result, args.json, args.csv, args.report)


if __name__ == "__main__":
    main()
