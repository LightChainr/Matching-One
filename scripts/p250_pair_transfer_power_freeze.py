#!/usr/bin/env python3
"""Freeze a pair-transfer sample size from the existing P250 10k stream."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from score_z5_charged_multiseparation import denominator_rows, read_batches


GRID = (10_000, 20_000, 40_000, 80_000)
SEPARATIONS = (1, 2, 3)
MINIMUM_Z = 5.0


def freeze(path: Path) -> dict:
    rows = read_batches(path)
    baseline_samples = sum(row["samples"] for row in rows)
    if baseline_samples != 10_000:
        raise ValueError("power source must be the frozen 10k stream")
    source = {}
    for separation in SEPARATIONS:
        channels = denominator_rows(rows, separation)
        source[str(separation)] = {
            "channel_abs_z": {name: value["abs_z"] for name, value in channels.items()},
            "minimum_real_abs_z": min(value["abs_z"] for value in channels.values()),
        }
    table = []
    for samples in GRID:
        projected = {
            separation: source[str(separation)]["minimum_real_abs_z"]
            * math.sqrt(samples / baseline_samples)
            for separation in SEPARATIONS
        }
        table.append({
            "samples": samples,
            "projected_minimum_real_abs_z": projected,
            "qualifies": min(projected.values()) >= MINIMUM_Z,
        })
    selected = next(row["samples"] for row in table if row["qualifies"])
    return {
        "schema": "matching-one/p250-pair-transfer-power-freeze/v1",
        "source": {
            "path": str(path),
            "samples": baseline_samples,
            "uses": "Hermitian two-point rows only; these equal Re(T_r) exactly",
            "excluded": ["cubic support", "cubic phase", "transfer-shape residuals"],
        },
        "primary_separations": list(SEPARATIONS),
        "minimum_real_abs_z": MINIMUM_Z,
        "source_resolution": source,
        "grid": table,
        "selection_rule": "first grid point with projected weakest Re(T) z>=5 at every d=1,2,3",
        "selected_samples": selected,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batches", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(freeze(args.batches), indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
