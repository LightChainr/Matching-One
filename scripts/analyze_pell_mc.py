#!/usr/bin/env python3
"""Analyze batch aggregates from ``src/pell_matching_mc.cpp``.

The engine evaluates the wrapping-side matching function at ``p_ref-h``,
``p_ref``, and ``p_ref+h`` with common random numbers.  This script reports
batch-jackknife uncertainty for the point values, central derivative, linear
root, and paired axis/diamond combinations.  It consumes aggregate counts;
individual configurations are neither needed nor written.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Optional, Union


def read_metadata(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("metadata must contain a JSON object")
    return value


def read_batches(path: Path) -> tuple[list[int], list[float], dict[tuple[int, float], dict[str, int]]]:
    rows: dict[tuple[int, float], dict[str, int]] = {}
    batches: set[int] = set()
    probabilities: set[float] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"batch", "samples", "p", "axis_sum", "diamond_sum"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"missing CSV fields: {sorted(required - set(reader.fieldnames or []))}")
        for raw in reader:
            batch = int(raw["batch"])
            p = float(raw["p"])
            key = (batch, p)
            if key in rows:
                raise ValueError(f"duplicate batch/probability row: {key}")
            rows[key] = {
                "samples": int(raw["samples"]),
                "axis": int(raw["axis_sum"]),
                "diamond": int(raw["diamond_sum"]),
            }
            batches.add(batch)
            probabilities.add(p)
    ordered_batches = sorted(batches)
    ordered_p = sorted(probabilities)
    if len(ordered_batches) < 2:
        raise ValueError("at least two batches are required")
    if len(ordered_p) != 3:
        raise ValueError(f"expected exactly three probabilities, found {len(ordered_p)}")
    expected = {(batch, p) for batch in ordered_batches for p in ordered_p}
    if set(rows) != expected:
        raise ValueError("CSV does not contain a complete batch x probability grid")
    sizes = {rows[key]["samples"] for key in rows}
    if len(sizes) != 1:
        raise ValueError("jackknife requires equal batch sizes")
    return ordered_batches, ordered_p, rows


def estimate_grid(
    batches: list[int],
    probabilities: list[float],
    rows: dict[tuple[int, float], dict[str, int]],
    omitted: Optional[int],
) -> dict[str, list[float]]:
    result = {"axis": [], "diamond": [], "difference": []}
    for p in probabilities:
        selected = [batch for batch in batches if batch != omitted]
        samples = sum(rows[(batch, p)]["samples"] for batch in selected)
        axis = sum(rows[(batch, p)]["axis"] for batch in selected) / samples
        diamond = sum(rows[(batch, p)]["diamond"] for batch in selected) / samples
        result["axis"].append(axis)
        result["diamond"].append(diamond)
        result["difference"].append(diamond - axis)
    return result


def derived(grid: dict[str, list[float]], probabilities: list[float], metadata: dict[str, object]) -> dict[str, float]:
    p_minus, p_ref, p_plus = probabilities
    h_left = p_ref - p_minus
    h_right = p_plus - p_ref
    tolerance = 64 * math.ulp(max(abs(p_ref), 1.0))
    if abs(h_left - h_right) > tolerance:
        raise ValueError("probability grid is not symmetric around p_ref")
    h = (h_left + h_right) / 2
    output: dict[str, float] = {}
    roots: dict[str, float] = {}
    for channel in ("axis", "diamond"):
        y_minus, y_zero, y_plus = grid[channel]
        derivative = (y_plus - y_minus) / (2 * h)
        curvature = (y_plus - 2 * y_zero + y_minus) / (h * h)
        if derivative == 0:
            raise ValueError(f"zero central derivative for {channel}; increase samples or h")
        root = p_ref - y_zero / derivative
        roots[channel] = root
        output[f"{channel}_M_p_ref"] = y_zero
        output[f"{channel}_derivative"] = derivative
        output[f"{channel}_curvature"] = curvature
        output[f"{channel}_linear_root"] = root

    output["orientation_root_gap_diamond_minus_axis"] = roots["diamond"] - roots["axis"]
    output["simple_mean_root"] = (roots["axis"] + roots["diamond"]) / 2
    axis_length = float(metadata["axis_physical_period"])
    diamond_length = float(metadata["diamond_physical_period"])
    axis_weight = axis_length**4 / (axis_length**4 + diamond_length**4)
    diamond_weight = 1.0 - axis_weight
    output["L4_axis_weight"] = axis_weight
    output["L4_diamond_weight"] = diamond_weight
    output["L4_weighted_root"] = axis_weight * roots["axis"] + diamond_weight * roots["diamond"]
    output["paired_difference_M_p_ref"] = grid["difference"][1]
    output["paired_difference_derivative"] = (
        grid["difference"][2] - grid["difference"][0]
    ) / (2 * h)
    return output


def jackknife(
    full: dict[str, float],
    leave_one_out: list[dict[str, float]],
) -> dict[str, dict[str, float]]:
    count = len(leave_one_out)
    summary: dict[str, dict[str, float]] = {}
    for name, estimate in full.items():
        values = [row[name] for row in leave_one_out]
        center = math.fsum(values) / count
        se = math.sqrt((count - 1) / count * math.fsum((value - center) ** 2 for value in values))
        summary[name] = {"estimate": estimate, "jackknife_se": se}
    return summary


def point_summaries(
    batches: list[int],
    probabilities: list[float],
    rows: dict[tuple[int, float], dict[str, int]],
) -> list[dict[str, Union[float, str]]]:
    output: list[dict[str, Union[float, str]]] = []
    for point, p in enumerate(probabilities):
        for channel in ("axis", "diamond", "difference"):
            values = []
            for batch in batches:
                row = rows[(batch, p)]
                numerator = (
                    row["axis"] if channel == "axis" else
                    row["diamond"] if channel == "diamond" else
                    row["diamond"] - row["axis"]
                )
                values.append(numerator / row["samples"])
            mean = math.fsum(values) / len(values)
            se = math.sqrt(
                math.fsum((value - mean) ** 2 for value in values)
                / (len(values) * (len(values) - 1))
            )
            output.append({"point": point, "p": p, "channel": channel, "estimate": mean, "batch_se": se})
    return output


def write_summary_csv(path: Path, summary: dict[str, dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["metric", "estimate", "jackknife_se"],
            lineterminator="\n",
        )
        writer.writeheader()
        for metric, values in summary.items():
            writer.writerow({"metric": metric, **values})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", type=Path, required=True, help="*.batches.csv from the C++ engine")
    parser.add_argument("--metadata", type=Path, required=True, help="*.metadata.json from the C++ engine")
    parser.add_argument("--json", type=Path, required=True, help="analysis JSON output")
    parser.add_argument("--csv", type=Path, required=True, help="flat derived-metric CSV output")
    args = parser.parse_args()

    metadata = read_metadata(args.metadata)
    batches, probabilities, rows = read_batches(args.batches)
    full_grid = estimate_grid(batches, probabilities, rows, omitted=None)
    full = derived(full_grid, probabilities, metadata)
    leave_one_out = [
        derived(estimate_grid(batches, probabilities, rows, omitted=batch), probabilities, metadata)
        for batch in batches
    ]
    summary = jackknife(full, leave_one_out)
    points = point_summaries(batches, probabilities, rows)
    payload = {
        "analysis": "three-point central-difference with delete-one-batch jackknife",
        "source_batches": str(args.batches),
        "source_metadata": str(args.metadata),
        "batch_count": len(batches),
        "probabilities": probabilities,
        "point_estimates": points,
        "derived": summary,
        "caveats": [
            "linear_root uses only M(p_ref) and the central finite-difference derivative",
            "finite-difference truncation error is not included in the jackknife standard error",
            "fixed-p discovery scans are not a replacement for a full Newman-Ziff production run",
        ],
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_summary_csv(args.csv, summary)

    for name in (
        "axis_linear_root", "diamond_linear_root",
        "orientation_root_gap_diamond_minus_axis", "simple_mean_root", "L4_weighted_root",
    ):
        row = summary[name]
        print(f"{name}: {row['estimate']:.12g} +/- {row['jackknife_se']:.3g}")
    print(f"wrote {args.json}")
    print(f"wrote {args.csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
