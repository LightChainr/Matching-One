#!/usr/bin/env python3
"""Single authorized 5k paired anchor pilot for Issue 418."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

from norm5_chiral_fixedp_mc import P_FIXED, splitmix64
from z5_projective_leg_bivariate_mc import (
    GRID,
    gauge_charged_rows,
    label,
    pair_value,
    rotation_gauges,
    translated_points,
)
from z5_projective_leg_cross_scale_mc import (
    CHILD_ORDER,
    CHARGES,
    HANDS,
    PARENT_GEOMETRY,
    charged_rows,
    contexts,
    counter_uniform,
)


SCHEMA = "matching-one/p418-anchor-paired-pilot/v1"
ROOT = Path(__file__).resolve().parents[1]
FREEZE = ROOT / "analysis/p418_anchor_paired_pilot_freeze.json"
DEFAULT_OUTPUT = ROOT / "results/local-20260830/P418-anchor-paired-5k/response_5k.json"
ESTIMATORS = ("current", "independent", "full")
MASK64 = (1 << 64) - 1


def observable_order() -> tuple[str, ...]:
    return tuple(
        f"{estimator}__{label(a, b)}_r{charge}_{hand}_{part}"
        for estimator in ESTIMATORS
        for a, b in GRID
        for hand in HANDS
        for charge in CHARGES
        for part in ("re", "im")
    )


FIELD_ORDER = observable_order()


def independent_anchor(seed: int, replica: int) -> int:
    key = seed ^ splitmix64(replica + 0x418505101A11)
    return splitmix64(key & MASK64) % PARENT_GEOMETRY.n


def targets(origin: int) -> dict[tuple[int, int], int]:
    x, y = PARENT_GEOMETRY.coordinates[origin]
    return {
        displacement: PARENT_GEOMETRY.vertex((x + displacement[0], y + displacement[1]))
        for displacement in GRID
    }


TARGETS_BY_ORIGIN = tuple(targets(origin) for origin in range(PARENT_GEOMETRY.n))
_GAUGES = None


def cached_gauges():
    global _GAUGES
    if _GAUGES is None:
        _GAUGES = rotation_gauges()
    return _GAUGES


def replica_observables(field_seed: int, independent_seed: int, replica: int) -> tuple[dict[str, float], int, int, bytes]:
    field = [counter_uniform(field_seed, replica, site) < P_FIXED for site in range(CHILD_ORDER)]
    current, _ = translated_points(field_seed, replica)
    independent = independent_anchor(independent_seed, replica)
    anchor_targets = {
        "current": TARGETS_BY_ORIGIN[current],
        "independent": TARGETS_BY_ORIGIN[independent],
    }
    gauges = cached_gauges()
    output = {name: 0.0 for name in FIELD_ORDER}
    parents = range(PARENT_GEOMETRY.n)
    for hand, context in zip(HANDS, contexts()):
        rows = gauge_charged_rows(charged_rows(context, field, parents), parents, gauges[hand])
        for displacement in GRID:
            for charge in CHARGES:
                for estimator, origin in (("current", current), ("independent", independent)):
                    value = pair_value(rows, origin, anchor_targets[estimator][displacement], charge)
                    prefix = f"{estimator}__{label(*displacement)}_r{charge}_{hand}_"
                    output[prefix + "re"] = value.real
                    output[prefix + "im"] = value.imag
                full = sum(
                    pair_value(rows, origin, TARGETS_BY_ORIGIN[origin][displacement], charge)
                    for origin in parents
                ) / PARENT_GEOMETRY.n
                prefix = f"full__{label(*displacement)}_r{charge}_{hand}_"
                output[prefix + "re"] = full.real
                output[prefix + "im"] = full.imag
    return output, current, independent, bytes(field)


def _run_batch(task: tuple[int, int, int, int, int]) -> dict[str, object]:
    batch, start, samples, field_seed, independent_seed = task
    sums = {name: 0.0 for name in FIELD_ORDER}
    field_digest = hashlib.sha256()
    current_digest = hashlib.sha256()
    independent_digest = hashlib.sha256()
    anchor_matches = 0
    for replica in range(start, start + samples):
        values, current, independent, field = replica_observables(field_seed, independent_seed, replica)
        field_digest.update(field)
        current_digest.update(current.to_bytes(2, "little"))
        independent_digest.update(independent.to_bytes(2, "little"))
        anchor_matches += int(current == independent)
        for name, value in values.items():
            sums[name] += value
    return {
        "batch": batch,
        "replica_first": start,
        "samples": samples,
        "field_sha256": field_digest.hexdigest(),
        "current_anchor_sha256": current_digest.hexdigest(),
        "independent_anchor_sha256": independent_digest.hexdigest(),
        "anchor_matches": anchor_matches,
        **sums,
    }


def validate_freeze(args: argparse.Namespace) -> dict[str, object]:
    freeze = json.loads(FREEZE.read_text())
    if freeze["status"] != "authorized_single_5k_paired_pilot":
        raise ValueError("paired pilot is not authorized")
    observed = {
        "samples": args.samples,
        "batches": args.batches,
        "workers": args.workers,
        "p": args.p,
        "field_seed": args.field_seed,
        "replica_offset": args.replica_offset,
        "independent_anchor_seed": args.independent_anchor_seed,
    }
    for key, value in observed.items():
        if freeze[key] != value:
            raise ValueError(f"run differs from freeze for {key}")
    if args.samples != 5_000 or args.samples % args.batches:
        raise ValueError("the only authorized run is 5k with equal batches")
    return freeze


def run(args: argparse.Namespace) -> list[dict[str, object]]:
    per_batch = args.samples // args.batches
    tasks = [
        (
            batch,
            args.replica_offset + batch * per_batch,
            per_batch,
            args.field_seed,
            args.independent_anchor_seed,
        )
        for batch in range(args.batches)
    ]
    if args.workers == 1:
        return [_run_batch(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        return list(pool.map(_run_batch, tasks))


def write_batches(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=5_000)
    parser.add_argument("--batches", type=int, default=100)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--p", type=float, default=P_FIXED)
    parser.add_argument("--field-seed", type=int, default=25050510120261130)
    parser.add_argument("--independent-anchor-seed", type=int, default=41850510120260830)
    parser.add_argument("--replica-offset", type=int, default=0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--batches-output", type=Path)
    args = parser.parse_args()
    freeze = validate_freeze(args)
    rows = run(args)
    batches = args.batches_output or args.output.with_suffix(".batches.csv")
    write_batches(batches, rows)
    payload = {
        "schema": SCHEMA,
        "status": "completed_single_authorized_5k_paired_pilot",
        "issues": [418, 250],
        "freeze": str(FREEZE.relative_to(ROOT)),
        "run": {
            key: getattr(args, key)
            for key in ("samples", "batches", "workers", "p", "field_seed", "independent_anchor_seed", "replica_offset")
        },
        "estimators": freeze["estimators"],
        "coordinates": [list(point) for point in GRID],
        "batch_rows_retain_full_cross_estimator_channel_covariance": True,
        "anchor_matches": sum(int(row["anchor_matches"]) for row in rows),
        "batches_output": str(batches.resolve().relative_to(ROOT)),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
