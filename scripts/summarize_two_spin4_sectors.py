#!/usr/bin/env python3
"""Summarize fixed-exponent spin-4 amplitudes from a same-N analysis CSV.

Expected input is the long-form output of analyze_gaussian_orientation_mc.py,
containing channel/sector rows and `hypothesis_scaled_amplitude` columns.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable


def weighted_constant(rows: Iterable[tuple[float, float]]) -> dict[str, float]:
    values = list(rows)
    if not values:
        raise ValueError("no rows selected")
    if any(se <= 0 or not math.isfinite(se) for _value, se in values):
        raise ValueError("all standard errors must be finite and positive")
    weights = [1.0 / se**2 for _value, se in values]
    mean = math.fsum(w * value for w, (value, _se) in zip(weights, values)) / math.fsum(weights)
    standard_error = math.sqrt(1.0 / math.fsum(weights))
    chi_square = math.fsum(((value - mean) / se) ** 2 for value, se in values)
    return {
        "mean": mean,
        "standard_error": standard_error,
        "chi_square": chi_square,
        "degrees_of_freedom": len(values) - 1,
        "z_from_zero": mean / standard_error,
    }


def load(path: Path, channel: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    even = []
    matching = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "N", "channel", "sector", "hypothesis_N_exponent",
            "hypothesis_scaled_amplitude", "hypothesis_scaled_batch_se",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("missing columns: " + ", ".join(sorted(missing)))
        for raw in reader:
            if raw["channel"] != channel:
                continue
            sector = raw["sector"]
            if sector not in ("even", "matching_function"):
                continue
            row = {
                "N": int(raw["N"]),
                "amplitude": float(raw["hypothesis_scaled_amplitude"]),
                "standard_error": float(raw["hypothesis_scaled_batch_se"]),
                "N_exponent": float(raw["hypothesis_N_exponent"]),
            }
            (even if sector == "even" else matching).append(row)
    even.sort(key=lambda row: row["N"])
    matching.sort(key=lambda row: row["N"])
    if [row["N"] for row in even] != [row["N"] for row in matching]:
        raise ValueError("even and matching-function sizes do not match")
    return even, matching


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--channel", default="either")
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    even, matching = load(args.csv, args.channel)
    even_fit = weighted_constant((row["amplitude"], row["standard_error"]) for row in even)
    matching_fit = weighted_constant(
        (row["amplitude"], row["standard_error"]) for row in matching
    )

    # D=M/2. The raw ratio D/S scales as N^[-(13/8-1)] = N^-5/8.
    a_d = matching_fit["mean"] / 2.0
    a_s = even_fit["mean"]
    ratio_prefactor = a_d / a_s
    crossover_N = ratio_prefactor ** (1.0 / (13.0 / 8.0 - 1.0))

    payload = {
        "channel": args.channel,
        "sizes": [row["N"] for row in even],
        "matching_even_model": "P4[S] = A_I N^-1",
        "matching_even": even_fit,
        "matching_function_model": "P4[M] = A_M N^-13/8",
        "matching_function": matching_fit,
        "D_over_S_prediction": {
            "formula": "P4[D]/P4[S] = (A_M/2/A_I) N^-5/8",
            "prefactor": ratio_prefactor,
            "crossover_N": crossover_N,
            "crossover_linear_length": math.sqrt(crossover_N),
        },
        "rows": {"even": even, "matching_function": matching},
    }
    text = json.dumps(payload, indent=2)
    print(text)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
