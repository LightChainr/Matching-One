#!/usr/bin/env python3
"""Stream a marked-birth archive into frozen two-time rank sufficient statistics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


SCHEMA = "matching-one/p334-two-time-kernel-summary/v1"
Z_GRID = (-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5)
ORIENTATIONS = ("first", "second")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def layers(n: int, k0: int) -> list[int]:
    width = n ** (5.0 / 8.0)
    return [math.floor(k0 + z * width + 0.5) for z in Z_GRID]


def empty_row() -> dict[str, object]:
    size = len(Z_GRID)
    return {
        "samples": 0,
        "sum_r": [0] * size,
        "sum_rr": [[0] * size for _ in range(size)],
        "sum_f1": [0] * size,
        "sum_f2": [0] * size,
        "sum_joint": [[0] * size for _ in range(size)],
    }


def summarize(path: Path, n: int, k0: int, expected_sha256: str) -> dict[str, object]:
    actual_sha256 = sha256(path)
    if actual_sha256 != expected_sha256:
        raise ValueError(f"raw hash changed: {actual_sha256}")
    grid = layers(n, k0)
    rows: dict[tuple[str, int], dict[str, object]] = {}
    per_batch: dict[tuple[str, int], int] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"n", "orientation", "batch", "samples", "k1", "k2", "count"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("marked-birth header changed")
        for source in reader:
            orientation = source["orientation"]
            if int(source["n"]) != n or orientation not in ORIENTATIONS:
                raise ValueError("size/orientation contract changed")
            batch = int(source["batch"])
            key = orientation, batch
            samples = int(source["samples"])
            if key in per_batch and per_batch[key] != samples:
                raise ValueError("samples-per-batch changed within a batch")
            per_batch[key] = samples
            k1, k2, count = int(source["k1"]), int(source["k2"]), int(source["count"])
            if not 1 <= k1 <= k2 <= n or count <= 0:
                raise ValueError("birth support changed")
            row = rows.setdefault(key, empty_row())
            row["samples"] += count
            rank = [int(k1 <= layer) + int(k2 <= layer) for layer in grid]
            for i, value in enumerate(rank):
                row["sum_r"][i] += count * value
                row["sum_f1"][i] += count * int(k1 <= grid[i])
                row["sum_f2"][i] += count * int(k2 <= grid[i])
                for j, other in enumerate(rank):
                    row["sum_rr"][i][j] += count * value * other
                    row["sum_joint"][i][j] += count * int(k1 <= grid[i] and k2 <= grid[j])
    if len(rows) != 200:
        raise ValueError(f"expected 2x100 batches, found {len(rows)}")
    output_rows = []
    for (orientation, batch), row in sorted(rows.items()):
        if row["samples"] != per_batch[(orientation, batch)]:
            raise ValueError(f"{orientation} batch {batch} count mismatch")
        output_rows.append({"orientation": orientation, "batch": batch, **row})
    return {
        "schema": SCHEMA,
        "source": str(path),
        "source_sha256": actual_sha256,
        "N": n,
        "k0": k0,
        "z_grid": list(Z_GRID),
        "layers": grid,
        "batches": output_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--k0", type=int, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = summarize(args.input, args.n, args.k0, args.sha256)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
