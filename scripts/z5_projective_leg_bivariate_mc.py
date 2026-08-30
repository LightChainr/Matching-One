#!/usr/bin/env python3
"""Fresh C4-closed mixed-displacement projective-leg pair stream for P250."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

from norm5_chiral_fixedp_mc import PHASES, P_FIXED, splitmix64
from z5_projective_leg_cross_scale_mc import (
    CHILD_ORDER,
    HANDS,
    CHARGES,
    PARENT_GEOMETRY,
    contexts,
    counter_uniform,
    charged_rows,
    exact_gate as cross_scale_exact_gate,
)
from z5_charged_threepoint_mc import covariance_of_mean


RADIUS = 4
GRID = tuple(
    (a, b)
    for radius in range(RADIUS + 1)
    for a in range(-radius, radius + 1)
    for b in range(-radius, radius + 1)
    if abs(a) + abs(b) == radius
)
FIRST_QUADRANT = tuple((a, b) for a, b in GRID if a >= 0 and b >= 0)
SCHEMA = "matching-one/z5-projective-leg-bivariate-response/v1"
DEFAULT_CAP = 2_000
MASK64 = (1 << 64) - 1


def label(a: int, b: int) -> str:
    return f"a{'p' if a >= 0 else 'm'}{abs(a)}_b{'p' if b >= 0 else 'm'}{abs(b)}"


def field_order() -> tuple[str, ...]:
    return tuple(
        f"{label(a, b)}_r{charge}_{hand}_{part}"
        for a, b in GRID
        for hand in HANDS
        for charge in CHARGES
        for part in ("re", "im")
    )


FIELD_ORDER = field_order()


def translated_points(seed: int, replica: int) -> tuple[int, dict[tuple[int, int], int]]:
    translation_index = splitmix64(seed ^ splitmix64(replica + 0x250505101B1A)) % PARENT_GEOMETRY.n
    x0, y0 = PARENT_GEOMETRY.coordinates[translation_index]
    points = {
        displacement: PARENT_GEOMETRY.vertex((x0 + displacement[0], y0 + displacement[1]))
        for displacement in GRID
    }
    if len(set(points.values())) != len(GRID):
        raise AssertionError("mixed-displacement diamond collapsed")
    return translation_index, points


def pair_value(rows: Mapping[int, Mapping[int, complex]], origin: int, target: int, charge: int) -> complex:
    return rows[origin][charge] * rows[target][(-charge) % 5]


def rotate(point: tuple[int, int]) -> tuple[int, int]:
    return -point[1], point[0]


def _rotation_gauge_for_context(context, multiplier: int) -> tuple[list[int], list[int], int]:
    """Trivialize the affine C4 fiber cocycle.

    Raw section labels transform as ``j' = multiplier*j + shift(parent)``.
    A site phase ``zeta^(r*gauge(parent))`` removes the shift when
    ``gauge(Rp)=shift(p)+multiplier*gauge(p)`` modulo five.
    """
    rotation = []
    shifts = []
    coordinate_failures = 0
    for parent_index, representative in enumerate(PARENT_GEOMETRY.coordinates):
        rotated_parent = PARENT_GEOMETRY.vertex(rotate(representative))
        rotation.append(rotated_parent)
        rotated_zero = context.geometry.vertex(rotate(context.field_coordinates[5 * parent_index]))
        candidates = [
            fiber for fiber in range(5)
            if context.field_to_vertex[5 * rotated_parent + fiber] == rotated_zero
        ]
        if len(candidates) != 1:
            raise AssertionError("C4 affine fiber shift is not unique")
        shift = candidates[0]
        shifts.append(shift)
        for fiber in range(5):
            field = 5 * parent_index + fiber
            rotated_point = rotate(context.field_coordinates[field])
            target = 5 * rotated_parent + (shift + multiplier * fiber) % 5
            coordinate_failures += int(context.geometry.vertex(rotated_point) != context.field_to_vertex[target])

    gauge: list[int | None] = [None] * PARENT_GEOMETRY.n
    for start in range(PARENT_GEOMETRY.n):
        if gauge[start] is not None:
            continue
        orbit = []
        point = start
        while point not in orbit:
            if gauge[point] is not None:
                raise AssertionError("C4 orbits unexpectedly overlap")
            orbit.append(point)
            point = rotation[point]
        if point != start:
            raise AssertionError("C4 orbit did not close at its start")
        solution = None
        for initial in range(5):
            trial = {start: initial}
            point = start
            valid = True
            for _ in orbit:
                target = rotation[point]
                value = (shifts[point] + multiplier * trial[point]) % 5
                if target in trial and trial[target] != value:
                    valid = False
                    break
                trial[target] = value
                point = target
            if valid:
                solution = trial
                break
        if solution is None:
            raise AssertionError("C4 affine cocycle has no covariant gauge")
        for parent_index in orbit:
            gauge[parent_index] = solution[parent_index]
    if any(value is None for value in gauge):
        raise AssertionError("C4 gauge is incomplete")
    gauge_values = [int(value) for value in gauge]
    equation_failures = sum(
        (gauge_values[rotation[parent]] - shifts[parent] - multiplier * gauge_values[parent]) % 5 != 0
        for parent in range(PARENT_GEOMETRY.n)
    )
    return gauge_values, shifts, coordinate_failures + equation_failures


def rotation_gauges() -> dict[str, list[int]]:
    plus, minus = contexts()
    plus_gauge, _, plus_failures = _rotation_gauge_for_context(plus, 3)
    minus_gauge, _, minus_failures = _rotation_gauge_for_context(minus, 2)
    if plus_failures or minus_failures:
        raise AssertionError("C4 gauge construction failed")
    return {"plus": plus_gauge, "minus": minus_gauge}


def gauge_charged_rows(rows, parent_indices: Sequence[int], gauge: Sequence[int]):
    return {
        parent: {
            charge: complex(*PHASES[(charge * gauge[parent]) % 5]) * rows[parent][charge]
            for charge in range(1, 5)
        }
        for parent in parent_indices
    }


def _rotation_fiber_gate() -> dict:
    multipliers = {"plus": 3, "minus": 2}
    failures = 0
    shifts_payload = {}
    gauge_payload = {}
    for hand, context in zip(HANDS, contexts()):
        gauge, shifts, row_failures = _rotation_gauge_for_context(context, multipliers[hand])
        failures += row_failures
        shifts_payload[hand] = shifts
        gauge_payload[hand] = gauge
    channel_maps = {
        "plus": {"r1": "conjugate(r2)", "r2": "r1"},
        "minus": {"r1": "r2", "r2": "conjugate(r1)"},
    }
    # Both maps have order four on the realified (r1,r2) pair.
    return {
        "fiber_multipliers_mod5": multipliers,
        "raw_affine_shifts_mod5": shifts_payload,
        "covariant_gauge_mod5": gauge_payload,
        "gauge_equation": "t(Rx)=s(x)+k*t(x) mod 5",
        "channel_maps": channel_maps,
        "channel_map_fourth_power": "identity",
        "failures": failures,
        "passed": failures == 0,
    }


def exact_gate() -> dict:
    base = cross_scale_exact_gate()
    points = {point: PARENT_GEOMETRY.vertex(point) for point in GRID}
    rotation = _rotation_fiber_gate()
    required = {(1, 1), (2, 1), (1, 2)}
    passed = (
        base["passed"]
        and len(GRID) == 41
        and len(set(points.values())) == len(GRID)
        and required.issubset(GRID)
        and all(rotate(point) in GRID for point in GRID)
        and rotation["passed"]
    )
    if not passed:
        raise AssertionError("bivariate projective-leg exact gate failed")
    return {
        "cross_scale_gate": base,
        "grid": [list(point) for point in GRID],
        "first_quadrant": [list(point) for point in FIRST_QUADRANT],
        "grid_points": len(GRID),
        "distinct_parent_vertices": len(set(points.values())),
        "required_commuting_gate_points": [[1, 1], [2, 1], [1, 2]],
        "C4_closed": True,
        "rotation_fiber_gate": rotation,
        "coordinates": len(FIELD_ORDER),
        "passed": True,
    }


def _run_batch(task: tuple[int, int, int, float, int]) -> dict:
    batch, start, samples, p, seed = task
    plus, minus = contexts()
    sums = {name: 0.0 for name in FIELD_ORDER}
    field_digest = hashlib.sha256()
    translation_digest = hashlib.sha256()
    gauges = rotation_gauges()
    for replica in range(start, start + samples):
        field = [counter_uniform(seed, replica, site) < p for site in range(CHILD_ORDER)]
        field_digest.update(bytes(field))
        translation_index, points = translated_points(seed, replica)
        translation_digest.update(translation_index.to_bytes(2, "little"))
        unique = sorted(set(points.values()))
        origin = points[(0, 0)]
        for hand, context in (("plus", plus), ("minus", minus)):
            rows = gauge_charged_rows(charged_rows(context, field, unique), unique, gauges[hand])
            for (a, b), target in points.items():
                for charge in CHARGES:
                    value = pair_value(rows, origin, target, charge)
                    prefix = f"{label(a, b)}_r{charge}_{hand}_"
                    sums[prefix + "re"] += value.real
                    sums[prefix + "im"] += value.imag
    return {
        "batch": batch,
        "replica_first": start,
        "samples": samples,
        "field_sha256": field_digest.hexdigest(),
        "translation_sha256": translation_digest.hexdigest(),
        **sums,
    }


def summarize(rows: Sequence[Mapping[str, object]]) -> dict:
    batch_means = [[float(row[name]) / int(row["samples"]) for name in FIELD_ORDER] for row in rows]
    return {
        "order": list(FIELD_ORDER),
        "point": [sum(row[index] for row in batch_means) / len(batch_means) for index in range(len(FIELD_ORDER))],
        "covariance_of_mean": covariance_of_mean(batch_means),
        "batch_replicates": len(rows),
    }


def run(samples: int, batches: int, workers: int, p: float, seed: int, replica_offset: int, *, cap: int):
    if samples <= 0 or batches < 2 or samples % batches:
        raise ValueError("samples must be divisible by batches>=2")
    if samples > cap:
        raise ValueError(f"bivariate run exceeds authorized cap {cap}")
    per_batch = samples // batches
    tasks = [(batch, replica_offset + batch * per_batch, per_batch, p, seed) for batch in range(batches)]
    if workers == 1:
        rows = [_run_batch(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(_run_batch, tasks))
    return rows, summarize(rows)


def validate_manifest(args, manifest: Mapping[str, object]) -> int:
    if manifest.get("status") != "authorized_fresh_bivariate_state":
        raise ValueError("bivariate manifest is not authorized")
    observed = {
        "samples": args.samples,
        "batches": args.batches,
        "workers": args.workers,
        "p": args.p,
        "seed": args.seed,
        "replica_offset": args.replica_offset,
        "replica_last_exclusive": args.replica_offset + args.samples,
    }
    for key, value in observed.items():
        if manifest["run"].get(key) != value:
            raise ValueError(f"run differs from manifest for {key}")
    if manifest.get("grid") != [list(point) for point in GRID]:
        raise ValueError("mixed-displacement grid changed")
    return args.samples


def write_batches(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=2_000)
    parser.add_argument("--batches", type=int, default=40)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--p", type=float, default=P_FIXED)
    parser.add_argument("--seed", type=int, default=25050510120261130)
    parser.add_argument("--replica-offset", type=int, default=0)
    parser.add_argument("--production-manifest", type=Path)
    parser.add_argument("--exact-gate", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batches-output", type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.exact_gate:
        args.output.write_text(json.dumps(exact_gate(), indent=2) + "\n")
        return 0
    cap = DEFAULT_CAP
    manifest = None
    if args.production_manifest:
        manifest = json.loads(args.production_manifest.read_text())
        cap = validate_manifest(args, manifest)
    rows, analysis = run(args.samples, args.batches, args.workers, args.p, args.seed, args.replica_offset, cap=cap)
    batches_path = args.batches_output or args.output.with_suffix(".batches.csv")
    write_batches(batches_path, rows)
    payload = {
        "schema": SCHEMA,
        "status": "fresh_bivariate_state" if manifest else "bivariate_state_smoke",
        "issues": [249, 250, 255],
        "exact_gate": exact_gate(),
        "run": {
            "samples": args.samples,
            "batches": args.batches,
            "workers": args.workers,
            "p": args.p,
            "seed": args.seed,
            "replica_offset": args.replica_offset,
            "replica_last_exclusive": args.replica_offset + args.samples,
            "batches_output": str(batches_path),
        },
        "manifest_runner_commit": manifest.get("runner_commit") if manifest else None,
        "observable": {
            "operator": "C4-covariant-gauge neutral projective-leg charged pair G_r(a,b)",
            "parent_order": PARENT_GEOMETRY.n,
            "child_order": CHILD_ORDER,
            "grid": [list(point) for point in GRID],
            "charges": list(CHARGES),
            "hands": list(HANDS),
            "cubic_fields": [],
            "gauge": "exact affine-fiber C4 cocycle trivialization from exact_gate",
        },
        "analysis": analysis,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
