#!/usr/bin/env python3
"""Analyze aligned multi-radius plus/minus pivotal-H4 batch statistics."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


INTEGER_FIELDS = (
    "n", "a", "b", "radius", "batch", "counter_first",
    "counter_last_exclusive", "samples", "sum_score_t",
    "sum_score_lambda", "sum_global_twice", "global_twice_score_t",
    "global_twice_score_lambda", "black_pivotal", "white_pivotal",
    "black_h4", "white_h4", "h4_plus", "h4_minus",
    "h4_plus_score_t", "h4_plus_score_lambda", "h4_minus_score_t",
    "h4_minus_score_lambda",
)
COMMON_FIELDS = (
    "n", "a", "b", "batch", "counter_first", "counter_last_exclusive",
    "samples", "sum_score_t", "sum_score_lambda", "sum_global_twice",
    "global_twice_score_t", "global_twice_score_lambda", "black_pivotal",
    "white_pivotal",
)
SUM_FIELDS = tuple(field for field in INTEGER_FIELDS if field not in {"n", "a", "b", "radius", "batch", "counter_first", "counter_last_exclusive"})


def read_rows(path: Path) -> tuple[dict[tuple[int, int, int], dict[str, int]], list[int], list[int], list[int]]:
    rows: dict[tuple[int, int, int], dict[str, int]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not set(INTEGER_FIELDS).issubset(reader.fieldnames):
            raise ValueError("multi-radius batch schema mismatch")
        for raw in reader:
            row = {field: int(raw[field]) for field in INTEGER_FIELDS}
            key = row["n"], row["radius"], row["batch"]
            if key in rows:
                raise ValueError(f"duplicate row {key}")
            if row["h4_plus"] != row["black_h4"] + row["white_h4"]:
                raise ValueError(f"plus identity failed at {key}")
            if row["h4_minus"] != row["black_h4"] - row["white_h4"]:
                raise ValueError(f"minus identity failed at {key}")
            rows[key] = row
    sizes = sorted({key[0] for key in rows})
    radii = sorted({key[1] for key in rows})
    batches = sorted({key[2] for key in rows})
    if len(sizes) < 1 or len(radii) < 2 or len(batches) < 2:
        raise ValueError("need at least one size, two radii and two batches")
    expected = {(n, radius, batch) for n in sizes for radius in radii for batch in batches}
    if set(rows) != expected:
        raise ValueError("size/radius/batch grid is incomplete")
    for n in sizes:
        for batch in batches:
            reference = rows[n, radii[0], batch]
            for radius in radii[1:]:
                row = rows[n, radius, batch]
                if any(row[field] != reference[field] for field in COMMON_FIELDS):
                    raise ValueError(f"common-field alignment failed for N={n}, batch={batch}")
    return rows, sizes, radii, batches


def aggregate(rows: dict[tuple[int, int, int], dict[str, int]], n: int, radius: int,
              batches: list[int], omit: int | None = None) -> dict[str, int]:
    total = {field: 0 for field in SUM_FIELDS}
    for batch in batches:
        if batch == omit:
            continue
        row = rows[n, radius, batch]
        for field in SUM_FIELDS:
            total[field] += row[field]
    return total


def coordinates(total: dict[str, int], n: int, radius: int) -> dict[str, float]:
    pivotal = total["black_pivotal"] + total["white_pivotal"]
    if pivotal == 0 or total["samples"] == 0:
        raise ValueError(f"zero pivotal/sample denominator for N={n}, R={radius}")
    return {
        "N": n,
        "R": radius,
        "delta": radius / math.sqrt(n),
        "mu0": n * pivotal / total["samples"],
        "mu4_plus": n * total["h4_plus"] / total["samples"],
        "mu4_minus": n * total["h4_minus"] / total["samples"],
        "A_plus": total["h4_plus"] / pivotal,
        "A_minus": total["h4_minus"] / pivotal,
        "black_h4": total["black_h4"],
        "white_h4": total["white_h4"],
        "pivotal_events": pivotal,
        "samples": total["samples"],
    }


def vectors(rows, sizes, radii, batches, omit=None):
    points = {
        (n, radius): coordinates(aggregate(rows, n, radius, batches, omit), n, radius)
        for n in sizes for radius in radii
    }
    base_order = [f"N{n}_R{radius}_{channel}" for n in sizes for radius in radii for channel in ("A_plus", "A_minus")]
    base = [points[n, radius][channel] for n in sizes for radius in radii for channel in ("A_plus", "A_minus")]
    shell_order: list[str] = []
    shell: list[float] = []
    shell_rows: list[dict[str, float | int]] = []
    for n in sizes:
        for first, second in zip(radii, radii[1:]):
            log_ratio = math.log(second / first)
            record: dict[str, float | int] = {"N": n, "R_first": first, "R_second": second, "log_ratio": log_ratio}
            for channel in ("A_plus", "A_minus"):
                increment = points[n, second][channel] - points[n, first][channel]
                normalized = increment / log_ratio
                record[f"delta_{channel}"] = increment
                record[f"per_log_{channel}"] = normalized
                shell_order.append(f"N{n}_R{first}_to_R{second}_per_log_{channel}")
                shell.append(normalized)
            shell_rows.append(record)
    return points, base_order, base, shell_order, shell, shell_rows


def jackknife_covariance(replicates: list[list[float]]) -> list[list[float]]:
    count = len(replicates)
    width = len(replicates[0])
    means = [sum(row[index] for row in replicates) / count for index in range(width)]
    scale = (count - 1) / count
    return [
        [scale * sum((row[i] - means[i]) * (row[j] - means[j]) for row in replicates)
         for j in range(width)]
        for i in range(width)
    ]


def analyze(batch_path: Path, metadata_path: Path) -> dict[str, object]:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("schema") != "matching-one/c4-multiradius-pivotal/v1":
        raise ValueError("metadata schema mismatch")
    rows, sizes, radii, batches = read_rows(batch_path)
    if metadata.get("radii") != radii:
        raise ValueError("metadata/CSV radius mismatch")
    points, base_order, base, shell_order, shell, shell_rows = vectors(
        rows, sizes, radii, batches)
    delete_base: list[list[float]] = []
    delete_shell: list[list[float]] = []
    for omitted in batches:
        _, _, base_rep, _, shell_rep, _ = vectors(rows, sizes, radii, batches, omitted)
        delete_base.append(base_rep)
        delete_shell.append(shell_rep)
    base_cov = jackknife_covariance(delete_base)
    shell_cov = jackknife_covariance(delete_shell)
    base_se = [math.sqrt(max(base_cov[i][i], 0.0)) for i in range(len(base))]
    shell_se = [math.sqrt(max(shell_cov[i][i], 0.0)) for i in range(len(shell))]
    point_rows = []
    for n in sizes:
        for radius in radii:
            point = dict(points[n, radius])
            for channel in ("A_plus", "A_minus"):
                index = base_order.index(f"N{n}_R{radius}_{channel}")
                point[f"{channel}_SE"] = base_se[index]
            point_rows.append(point)
    for record in shell_rows:
        n = record["N"]
        first = record["R_first"]
        second = record["R_second"]
        for channel in ("A_plus", "A_minus"):
            label = f"N{n}_R{first}_to_R{second}_per_log_{channel}"
            record[f"per_log_{channel}_SE"] = shell_se[shell_order.index(label)]
    return {
        "schema": "matching-one/c4-multiradius-pivotal-analysis/v1",
        "source": {"batches": str(batch_path), "metadata": str(metadata_path)},
        "design": {
            "sizes": sizes,
            "radii": radii,
            "batches": len(batches),
            "radius_semantics": "fixed R uses lattice units; fixed delta requires R proportional to sqrt(N)",
            "same_block_warning": "all radii/channels/sizes are correlated views of one counter stream",
        },
        "points": point_rows,
        "base_vector": {"order": base_order, "point": base, "covariance": base_cov},
        "shell_increments": shell_rows,
        "shell_vector": {"order": shell_order, "point": shell, "covariance": shell_cov},
        "interpretation": {
            "fixed_R": "UV limit: delta=R/sqrt(N) tends to zero with lattice spacing",
            "fixed_delta": "two-cutoff limit: choose R(N)=round(delta*sqrt(N)) before N grows",
            "shell_test": "compare adjacent per-log increments jointly; do not count shells independently",
            "pilot_scope": "engineering smoke only; no Jordan or exponent claim",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batches", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.batches, args.metadata)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
