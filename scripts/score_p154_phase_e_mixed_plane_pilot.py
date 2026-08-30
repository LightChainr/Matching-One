#!/usr/bin/env python3
"""Score the frozen P154 Phase-E local-singlet mixed-plane pilot."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from score_p154_local_singlet_pilot import cos4, fit_ray, mean_cov, sha256


FREEZE_COMMIT = "0578105d92d3822cb48f5c421bd23ff339295cc6"
FIELDS = ("A", "E", "C", "W", "B", "J_top", "J_bulk", "B_var")
CANDIDATES = (
    ("A", "E", "C"),
    ("A", "E", "J_bulk"),
    ("A", "E", "J_top"),
    ("A", "E", "B"),
)
EXPECTED = {
    65: {"first": (8, 1), "second": (7, 4), "samples": 20000, "batches": 100,
         "seed": 202615465, "counter": (15466000000, 15466020000)},
    130: {"first": (11, 3), "second": (9, 7), "samples": 20000, "batches": 100,
          "seed": 2026154130, "counter": (15466200000, 15466220000)},
}


def row_values(row: dict[str, str]) -> dict[str, float]:
    n = int(row["n"])
    samples = int(row["samples"])
    if samples <= 1:
        raise ValueError("each batch needs at least two replicas")
    i0 = int(row["sum_i0"]) / samples
    i2 = int(row["sum_i2"]) / samples
    k1 = int(row["sum_k1"]) / samples
    k2 = int(row["sum_k2"]) / samples
    scale = 4.0 * n
    even_sum = int(row["sum_black_axis_pairs"]) + int(row["sum_white_matching_axis_pairs"])
    b = even_sum / (samples * scale)
    b2 = int(row["sum_even_numerator_squared"]) / (samples * scale * scale)
    i0b = int(row["sum_i0_even_numerator"]) / (samples * scale)
    i2b = int(row["sum_i2_even_numerator"]) / (samples * scale)
    unbiased = samples / (samples - 1.0)
    j0 = unbiased * (i0b - i0 * b)
    j2 = unbiased * (i2b - i2 * b)
    b_var = unbiased * (b2 - b * b)
    return {
        "A": i2 - i0,
        "E": i2 + i0,
        "C": (k1 + k2) / (2.0 * (n + 1)) - 0.5,
        "W": (k2 - k1) / (n + 1.0),
        "B": b,
        "J_top": j2 - j0,
        "J_bulk": j2 + j0,
        "B_var": b_var,
    }


def read_projected(csv_path: Path, metadata_path: Path, expected_commit: str):
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    n = int(metadata["N"])
    expected = EXPECTED[n]
    if metadata["git_commit"] != expected_commit:
        raise ValueError(f"{csv_path}: runner commit differs from frozen commit")
    checks = {
        "samples": expected["samples"],
        "batches": expected["batches"],
        "seed": expected["seed"],
        "replica_counter_first": expected["counter"][0],
        "replica_counter_last_exclusive": expected["counter"][1],
    }
    for key, wanted in checks.items():
        if int(metadata[key]) != wanted:
            raise ValueError(f"{csv_path}: {key} violates freeze")
    rows = {}
    with csv_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows[(int(row["batch"]), row["orientation"])] = row
    if len(rows) != 2 * expected["batches"]:
        raise ValueError(f"{csv_path}: incomplete batch/orientation table")
    delta = cos4(*expected["first"]) - cos4(*expected["second"])
    vectors = []
    for batch in range(expected["batches"]):
        first = rows[(batch, "first")]
        second = rows[(batch, "second")]
        if (int(first["a"]), int(first["b"])) != expected["first"]:
            raise ValueError("first orientation violates freeze")
        if (int(second["a"]), int(second["b"])) != expected["second"]:
            raise ValueError("second orientation violates freeze")
        x, y = row_values(first), row_values(second)
        vectors.append([(x[field] - y[field]) / delta for field in FIELDS])
    return n, vectors, metadata


def subvector(mean, covariance, fields):
    indices = [FIELDS.index(field) for field in fields]
    return ([mean[i] for i in indices],
            [[covariance[i][j] for j in indices] for i in indices])


def calculate(raw_dir: Path, expected_commit: str) -> dict:
    points = {}
    inputs = []
    for n in (65, 130):
        csv_path = raw_dir / f"n{n}_mixed.batches.csv"
        metadata_path = raw_dir / f"n{n}_mixed.metadata.json"
        actual_n, vectors, _ = read_projected(csv_path, metadata_path, expected_commit)
        points[actual_n] = mean_cov(vectors)
        inputs += [
            {"path": str(csv_path), "sha256": sha256(csv_path)},
            {"path": str(metadata_path), "sha256": sha256(metadata_path)},
        ]
    estimates = {}
    for n, (mean, covariance) in points.items():
        estimates[str(n)] = {
            field: {
                "estimate": mean[index],
                "se": math.sqrt(max(covariance[index][index], 0.0)),
                "z": mean[index] / math.sqrt(covariance[index][index]),
            }
            for index, field in enumerate(FIELDS)
        }
        estimates[str(n)]["covariance"] = covariance
    rays = {}
    for fields in CANDIDATES:
        first = subvector(*points[65], fields)
        second = subvector(*points[130], fields)
        scale, chi2 = fit_ray(first, second)
        rays["/".join(fields)] = {
            "fields": list(fields),
            "scale_N130_over_N65": scale,
            "chi_square": chi2,
            "df": 2,
            "p_value": math.exp(-0.5 * chi2),
        }
    improvement = rays["A/E/C"]["chi_square"] - rays["A/E/J_bulk"]["chi_square"]
    resolved = min(abs(estimates[str(n)]["J_bulk"]["z"]) for n in (65, 130)) >= 3.0
    return {
        "issue": 154,
        "phase": "E local-singlet mixed plane",
        "freeze_commit": expected_commit,
        "inputs": inputs,
        "projected_fields": list(FIELDS),
        "estimates": estimates,
        "common_ray_candidates": rays,
        "primary_decision": {
            "baseline": "A/E/C",
            "candidate": "A/E/J_bulk",
            "chi_square_improvement": improvement,
            "J_bulk_resolved_both_sizes": resolved,
            "decision": "mixed_local_plane_advances" if resolved and improvement >= 4.0
                        else "mixed_local_plane_not_selected_at_20k",
        },
        "interpretation_boundary": "finite-lattice matching-even mixed local-singlet response; not a continuum energy identity",
    }


def make_report(payload: dict) -> str:
    lines = [
        "# Issue #154 Phase-E mixed local-singlet pilot",
        "",
        "The pilot stores `B`, `B^2`, `I0*B`, and `I2*B` on the same stream,",
        "so `J_top=Cov(I2-I0,B)` and `J_bulk=Cov(I2+I0,B)` are directly scoreable.",
        "",
        "| N | P4[B] | P4[J_top] | P4[J_bulk] | P4[Var(B)] |",
        "|---:|---:|---:|---:|---:|",
    ]
    for n in (65, 130):
        row = payload["estimates"][str(n)]
        def cell(field: str) -> str:
            value = row[field]
            return f"{value['estimate']:.7g} +/- {value['se']:.2g} (z={value['z']:.2f})"
        lines.append(f"| {n} | {cell('B')} | {cell('J_top')} | {cell('J_bulk')} | {cell('B_var')} |")
    lines += ["", "## Common-plane score", "",
              "| candidate | scale N130/N65 | chi2 / 2 df | p |",
              "|---|---:|---:|---:|"]
    for name, row in payload["common_ray_candidates"].items():
        lines.append(f"| {name} | {row['scale_N130_over_N65']:.5g} | "
                     f"{row['chi_square']:.5g} | {row['p_value']:.4g} |")
    decision = payload["primary_decision"]
    lines += ["", "## Decision", "",
              f"`{decision['decision']}`. The A/E/J_bulk improvement over A/E/C is "
              f"{decision['chi_square_improvement']:.5g}; resolved at both sizes: "
              f"{decision['J_bulk_resolved_both_sizes']}.", "",
              "This is a finite-volume mixed-observation decision. It does not name `B` as the",
              "continuum thermal energy or turn a surviving plane into a field identity.", ""]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--expected-commit", default=FREEZE_COMMIT)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    payload = calculate(args.raw_dir, args.expected_commit)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.report.write_text(make_report(payload), encoding="utf-8")
    print(json.dumps(payload["primary_decision"], sort_keys=True))


if __name__ == "__main__":
    main()
