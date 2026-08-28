#!/usr/bin/env python3
"""Summarize batch output from ``src/gaussian_orientation_mc.cpp``.

The four sectors are formed after pairing, so every reported orientation
difference retains the common-random-number covariance.  ``even`` is
``(R_G+R_hat)/2`` and ``odd`` is ``(R_G-R_hat)/2``; ``matching_function`` is
twice ``odd`` and is included to avoid hidden factor-of-two conventions.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Tuple


BatchKey = Tuple[int, str, int]


def mean_se(values: List[float]) -> Tuple[float, float]:
    if len(values) < 2:
        raise ValueError("at least two batches are required")
    mean = math.fsum(values) / len(values)
    variance = math.fsum((value - mean) ** 2 for value in values)
    return mean, math.sqrt(variance / (len(values) * (len(values) - 1)))


def cos4(a: int, b: int) -> float:
    a2 = a * a
    b2 = b * b
    return (a2 * a2 - 6 * a2 * b2 + b2 * b2) / float((a2 + b2) ** 2)


def read_rows(path: Path) -> Tuple[Dict[BatchKey, Dict[str, int]], Dict[int, Tuple[int, int, int, int]]]:
    rows: Dict[BatchKey, Dict[str, int]] = {}
    designs: Dict[int, Tuple[int, int, int, int]] = {}
    required = {
        "n", "batch", "samples", "p_ref", "channel", "a1", "b1", "a2", "b2",
        "first_primal_sum", "first_matching_sum", "second_primal_sum", "second_matching_sum",
    }
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        if not required.issubset(fields):
            raise ValueError("missing CSV fields: " + ", ".join(sorted(required - fields)))
        for raw in reader:
            n = int(raw["n"])
            batch = int(raw["batch"])
            channel = raw["channel"]
            if channel not in ("either", "cross"):
                raise ValueError("unknown channel: " + channel)
            key = (n, channel, batch)
            if key in rows:
                raise ValueError("duplicate N/channel/batch row: " + repr(key))
            design = (int(raw["a1"]), int(raw["b1"]), int(raw["a2"]), int(raw["b2"]))
            if n in designs and designs[n] != design:
                raise ValueError("inconsistent representations for N=" + str(n))
            if design[0] ** 2 + design[1] ** 2 != n or design[2] ** 2 + design[3] ** 2 != n:
                raise ValueError("representation does not have advertised N")
            designs[n] = design
            samples = int(raw["samples"])
            if samples <= 0:
                raise ValueError("batch sample count must be positive")
            values = {
                "samples": samples,
                "first_primal": int(raw["first_primal_sum"]),
                "first_matching": int(raw["first_matching_sum"]),
                "second_primal": int(raw["second_primal_sum"]),
                "second_matching": int(raw["second_matching_sum"]),
            }
            if any(values[name] < 0 or values[name] > samples for name in values if name != "samples"):
                raise ValueError("indicator sum lies outside [0,samples]")
            rows[key] = values
    if not rows:
        raise ValueError("batch CSV is empty")
    for n in designs:
        for channel in ("either", "cross"):
            batches = sorted(key[2] for key in rows if key[:2] == (n, channel))
            if len(batches) < 2 or batches != list(range(len(batches))):
                raise ValueError("batches must be a complete zero-based range")
    return rows, designs


def sector_values(row: Dict[str, int], prefix: str) -> Dict[str, float]:
    samples = row["samples"]
    primal = row[prefix + "_primal"] / samples
    matching = row[prefix + "_matching"] / samples
    return {
        "primal": primal,
        "matching": matching,
        "even": 0.5 * (primal + matching),
        "odd": 0.5 * (primal - matching),
        "matching_function": primal - matching,
    }


def analyze(rows: Dict[BatchKey, Dict[str, int]], designs: Dict[int, Tuple[int, int, int, int]]) -> List[Dict[str, object]]:
    output: List[Dict[str, object]] = []
    for n in sorted(designs):
        a1, b1, a2, b2 = designs[n]
        delta_cos4 = cos4(a1, b1) - cos4(a2, b2)
        if delta_cos4 == 0:
            raise ValueError("orientation pair has zero delta cos(4 theta)")
        for channel in ("either", "cross"):
            keys = sorted((key for key in rows if key[:2] == (n, channel)), key=lambda key: key[2])
            for sector in ("primal", "matching", "even", "odd", "matching_function"):
                first_batches: List[float] = []
                second_batches: List[float] = []
                difference_batches: List[float] = []
                for key in keys:
                    first = sector_values(rows[key], "first")[sector]
                    second = sector_values(rows[key], "second")[sector]
                    first_batches.append(first)
                    second_batches.append(second)
                    difference_batches.append(first - second)
                first_mean, first_se = mean_se(first_batches)
                second_mean, second_se = mean_se(second_batches)
                difference, difference_se = mean_se(difference_batches)
                normalized = difference / delta_cos4
                normalized_se = difference_se / abs(delta_cos4)
                exponent = 13.0 / 8.0 if sector in ("odd", "matching_function") else 1.0
                scale = n ** exponent
                output.append({
                    "N": n,
                    "channel": channel,
                    "sector": sector,
                    "first_rep": [a1, b1],
                    "second_rep": [a2, b2],
                    "cos4_first": cos4(a1, b1),
                    "cos4_second": cos4(a2, b2),
                    "delta_cos4_first_minus_second": delta_cos4,
                    "first_estimate": first_mean,
                    "first_batch_se": first_se,
                    "second_estimate": second_mean,
                    "second_batch_se": second_se,
                    "difference_first_minus_second": difference,
                    "difference_batch_se": difference_se,
                    "difference_z": difference / difference_se if difference_se else None,
                    "normalized_by_delta_cos4": normalized,
                    "normalized_batch_se": normalized_se,
                    "hypothesis_N_exponent": exponent,
                    "hypothesis_scaled_amplitude": scale * normalized,
                    "hypothesis_scaled_batch_se": scale * normalized_se,
                })
    return output


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    scalar_fields = [name for name in rows[0] if name not in ("first_rep", "second_rep")]
    fields = ["a1", "b1", "a2", "b2"] + scalar_fields
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            flat = dict(row)
            first = flat.pop("first_rep")
            second = flat.pop("second_rep")
            flat.update({"a1": first[0], "b1": first[1], "a2": second[0], "b2": second[1]})
            writer.writerow(flat)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()

    with args.metadata.open(encoding="utf-8") as handle:
        metadata = json.load(handle)
    if not isinstance(metadata, dict):
        raise ValueError("metadata must contain an object")
    rows, designs = read_rows(args.batches)
    summaries = analyze(rows, designs)
    payload = {
        "analysis": "paired batch means for same-N Gaussian orientations",
        "source_batches": str(args.batches),
        "source_metadata": str(args.metadata),
        "metadata": metadata,
        "definitions": {
            "even": "(R_primal + R_white_matching)/2",
            "odd": "(R_primal - R_white_matching)/2",
            "matching_function": "R_primal - R_white_matching = 2*odd",
            "orientation_difference": "first representation minus second representation",
        },
        "scaling_diagnostics": {
            "primal_matching_even": "N * difference / delta_cos4 (L^-2 hypothesis)",
            "odd_matching_function": "N^(13/8) * difference / delta_cos4 (L^-13/4 hypothesis)",
        },
        "caveats": [
            "batch standard errors quantify Monte Carlo noise, not finite-size systematic error",
            "a resolved matching-even harmonic is required before interpreting a null matching-odd difference",
            "either and cross are separate prespecified wrapping conventions, not independent replicas",
        ],
        "summaries": summaries,
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_csv(args.csv, summaries)

    for row in summaries:
        if row["sector"] not in ("even", "matching_function"):
            continue
        z = row["difference_z"]
        z_text = "undefined" if z is None else "{:.2f}".format(z)
        print(
            "N={N} {channel} {sector}: diff={difference_first_minus_second:.6g} "
            "+/- {difference_batch_se:.3g}, z={z}, scaled={hypothesis_scaled_amplitude:.6g}"
            .format(z=z_text, **row)
        )
    print("wrote " + str(args.json))
    print("wrote " + str(args.csv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
