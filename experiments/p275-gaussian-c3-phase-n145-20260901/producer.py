#!/usr/bin/env python3
"""Independent paired N145 producer for the three-model P275 phase gate."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import json
from math import sqrt
from pathlib import Path
import sys
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from integer_period_torus import Matrix, integer_torus_geometry  # noqa: E402
from pinson_arguin_primitive import (  # noqa: E402
    engine_to_paper,
    primitive_probability_direct,
)
from square_bond_primitive_pilot import (  # noqa: E402
    TARGET_LINES,
    _splitmix64,
    classify_bond_mask,
)


G1: Matrix = ((12, -1), (1, 12))
G2: Matrix = ((9, -8), (8, 9))
DESIGNS = (("g1_12_plus_i", G1), ("g2_9_plus_8i", G2))
WEIGHTS = {
    "l0": (1.0, 0.0),
    "l1": (-0.5, -sqrt(3.0) / 2.0),
    "l2": (-0.5, +sqrt(3.0) / 2.0),
}
PRODUCTION_SAMPLES = 5_000_000
PRODUCTION_BATCHES = 100
PRODUCTION_SEED = 20_260_901_277


@dataclass(frozen=True)
class PairedBatch:
    batch: int
    replica_first: int
    samples: int
    z1_re: float
    z1_im: float
    z2_re: float
    z2_im: float


def validate_geometries() -> None:
    for identifier, matrix in DESIGNS:
        geometry = integer_torus_geometry(matrix)
        if geometry.n != 145 or len(geometry.primal_edges) != 290:
            raise ValueError(f"{identifier} is not the declared N145/290-edge torus")
    # g2*conj(g1)=116+87i, hence exp(i delta)=(4+3i)/5.
    if (9 + 8j) * (12 - 1j) != 116 + 87j:
        raise AssertionError("N145 Gaussian rotation arithmetic changed")


def _counter_mask(seed: int, replica: int, bits: int) -> int:
    words = (bits + 63) // 64
    value = 0
    for word in range(words):
        counter = replica * words + word
        draw = _splitmix64(seed ^ _splitmix64(counter + 1))
        value |= draw << (64 * word)
    return value & ((1 << bits) - 1)


def _run_batch(task: tuple[int, int, int, int]) -> PairedBatch:
    batch, replica_first, samples, seed = task
    geometries = [integer_torus_geometry(matrix) for _, matrix in DESIGNS]
    edge_count = len(geometries[0].primal_edges)
    sums = [[0.0, 0.0], [0.0, 0.0]]
    for replica in range(replica_first, replica_first + samples):
        mask = _counter_mask(seed, replica, edge_count)
        for index, geometry in enumerate(geometries):
            category, _ = classify_bond_mask(geometry, mask)
            if category == "invariant_failure":
                raise RuntimeError(
                    f"homology invariant failure in batch {batch}, replica {replica}"
                )
            weight = WEIGHTS.get(category)
            if weight is not None:
                sums[index][0] += weight[0]
                sums[index][1] += weight[1]
    return PairedBatch(
        batch,
        replica_first,
        samples,
        sums[0][0] / samples,
        sums[0][1] / samples,
        sums[1][0] / samples,
        sums[1][1] / samples,
    )


def primitive_c3_baseline(*, dps: int) -> tuple[float, float]:
    import mpmath as mp

    probabilities = [
        primitive_probability_direct(*engine_to_paper(line), mp.j, dps=dps)
        for _, line in TARGET_LINES
    ]
    re = probabilities[0] - (probabilities[1] + probabilities[2]) / 2
    im = mp.sqrt(3) * (probabilities[2] - probabilities[1]) / 2
    return float(re), float(im)


def run_batches(
    *, samples: int, batches: int, seed: int, workers: int
) -> list[PairedBatch]:
    validate_geometries()
    if samples <= 0 or batches <= 4 or samples % batches:
        raise ValueError("samples must be positive and divisible by batches>4")
    if workers <= 0:
        raise ValueError("workers must be positive")
    per_batch = samples // batches
    tasks = [(batch, batch * per_batch, per_batch, seed) for batch in range(batches)]
    if workers == 1:
        return [_run_batch(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(_run_batch, tasks))


def subtract_baseline(
    rows: Iterable[PairedBatch], baseline: tuple[float, float]
) -> list[PairedBatch]:
    re, im = baseline
    return [
        PairedBatch(
            row.batch,
            row.replica_first,
            row.samples,
            row.z1_re - re,
            row.z1_im - im,
            row.z2_re - re,
            row.z2_im - im,
        )
        for row in rows
    ]


def write_batches(path: Path, rows: Sequence[PairedBatch]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(path)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ("batch", "replica_first", "samples", "z1_re", "z1_im", "z2_re", "z2_im")
        )
        for row in sorted(rows, key=lambda item: item.batch):
            writer.writerow(
                (
                    row.batch,
                    row.replica_first,
                    row.samples,
                    format(row.z1_re, ".17g"),
                    format(row.z1_im, ".17g"),
                    format(row.z2_re, ".17g"),
                    format(row.z2_im, ".17g"),
                )
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("smoke", "production"), default="smoke")
    parser.add_argument("--samples", type=int, default=2_000)
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--seed", type=int, default=PRODUCTION_SEED)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.mode == "production" and (
        args.samples,
        args.batches,
        args.seed,
    ) != (PRODUCTION_SAMPLES, PRODUCTION_BATCHES, PRODUCTION_SEED):
        raise ValueError("production mode must use the frozen samples/batches/seed")
    baseline = primitive_c3_baseline(dps=args.dps)
    rows = subtract_baseline(
        run_batches(
            samples=args.samples,
            batches=args.batches,
            seed=args.seed,
            workers=args.workers,
        ),
        baseline,
    )
    write_batches(args.output, rows)
    metadata_path = args.output.with_suffix(args.output.suffix + ".metadata.json")
    if metadata_path.exists():
        raise FileExistsError(metadata_path)
    metadata_path.write_text(
        json.dumps(
            {
                "schema": "matching-one/p275-gaussian-c3-n145-paired-batches/v1",
                "mode": args.mode,
                "samples": args.samples,
                "batches": args.batches,
                "seed": args.seed,
                "workers": args.workers,
                "baseline_re_im": list(baseline),
                "period_matrices": {
                    identifier: [list(row) for row in matrix]
                    for identifier, matrix in DESIGNS
                },
                "retained_scientific_coordinates": [
                    "z1_re", "z1_im", "z2_re", "z2_im"
                ],
                "status": "smoke_only" if args.mode == "smoke" else "heldout_production",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"mode": args.mode, "samples": args.samples, "batches": args.batches, "output": str(args.output)}))


if __name__ == "__main__":
    main()
