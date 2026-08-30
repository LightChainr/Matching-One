#!/usr/bin/env python3
"""Frozen radius-five shell acquisition for the P250 sector-morphism test."""

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
    gauge_charged_rows,
    label,
    pair_value,
    rotate,
    rotation_gauges,
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


RADIUS = 5
SHELL = tuple(
    (a, b) for a in range(-RADIUS, RADIUS + 1) for b in range(-RADIUS, RADIUS + 1)
    if abs(a) + abs(b) == RADIUS
)
SOURCE_NEW = ((5, 0), (4, 1), (3, 2))
SCHEMA = "matching-one/z5-projective-leg-radius5-morphism-response/v1"
DEFAULT_CAP = 2_000
MASK64 = (1 << 64) - 1


def reflect(point: tuple[int, int]) -> tuple[int, int]:
    return point[0], -point[1]


def c4_orbit(point: tuple[int, int]) -> set[tuple[int, int]]:
    output = set()
    for _ in range(4):
        output.add(point)
        point = rotate(point)
    return output


def d4_closure(points: Sequence[tuple[int, int]]) -> set[tuple[int, int]]:
    return set().union(*(c4_orbit(point) | c4_orbit(reflect(point)) for point in points))


def field_order() -> tuple[str, ...]:
    return tuple(
        f"{label(a, b)}_r{charge}_{hand}_{part}"
        for a, b in SHELL for hand in HANDS for charge in CHARGES for part in ("re", "im")
    )


FIELD_ORDER = field_order()


def exact_gate() -> dict:
    vertices = {point: PARENT_GEOMETRY.vertex(point) for point in SHELL}
    closure = d4_closure(SOURCE_NEW)
    passed = (
        len(SHELL) == 20
        and len(set(vertices.values())) == 20
        and closure == set(SHELL)
        and all(rotate(point) in SHELL and reflect(point) in SHELL for point in SHELL)
    )
    if not passed:
        raise AssertionError("radius-five morphism geometry gate failed")
    return {
        "radius": RADIUS,
        "shell": [list(point) for point in SHELL],
        "shell_points": len(SHELL),
        "distinct_parent_vertices": len(set(vertices.values())),
        "source_degree5_endpoints": [list(point) for point in SOURCE_NEW],
        "source_C4_plus_Alexander_reflection_closure": [list(point) for point in sorted(closure)],
        "closure_is_full_shell": closure == set(SHELL),
        "minimality": "The three new endpoints in the u=(3,0) shifted annihilator generate all 20 shell points under C4 and Alexander reflection.",
        "hands": list(HANDS),
        "charges": list(CHARGES),
        "complex_rows_per_batch": len(SHELL) * len(HANDS) * len(CHARGES),
        "real_coordinates_per_batch": len(FIELD_ORDER),
        "passed": True,
    }


def translated_points(seed: int, replica: int) -> tuple[int, int, dict[tuple[int, int], int]]:
    translation_index = splitmix64(seed ^ splitmix64(replica + 0x250505101250)) % PARENT_GEOMETRY.n
    x0, y0 = PARENT_GEOMETRY.coordinates[translation_index]
    origin = PARENT_GEOMETRY.vertex((x0, y0))
    points = {
        displacement: PARENT_GEOMETRY.vertex((x0 + displacement[0], y0 + displacement[1]))
        for displacement in SHELL
    }
    if len(set(points.values())) != len(SHELL) or origin in points.values():
        raise AssertionError("radius-five translated shell collapsed")
    return translation_index, origin, points


def _run_batch(task: tuple[int, int, int, float, int]) -> dict:
    batch, start, samples, p, seed = task
    plus, minus = contexts()
    gauges = rotation_gauges()
    sums = {name: 0.0 for name in FIELD_ORDER}
    field_digest = hashlib.sha256()
    translation_digest = hashlib.sha256()
    for replica in range(start, start + samples):
        field = [counter_uniform(seed, replica, site) < p for site in range(CHILD_ORDER)]
        field_digest.update(bytes(field))
        translation_index, origin, points = translated_points(seed, replica)
        translation_digest.update(translation_index.to_bytes(2, "little"))
        unique = sorted({origin, *points.values()})
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


def validate_manifest(args, manifest: Mapping[str, object]) -> int:
    if not manifest.get("execution_authorized") or manifest.get("status") != "authorized_fresh_radius5_morphism":
        raise ValueError("radius-five morphism execution is not authorized")
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
    if manifest.get("geometry", {}).get("shell") != [list(point) for point in SHELL]:
        raise ValueError("radius-five shell changed")
    return args.samples


def run(samples: int, batches: int, workers: int, p: float, seed: int, replica_offset: int, cap: int):
    if samples <= 0 or batches < 2 or samples % batches:
        raise ValueError("samples must be divisible by batches>=2")
    if samples > cap:
        raise ValueError(f"radius-five run exceeds authorized cap {cap}")
    per_batch = samples // batches
    tasks = [(batch, replica_offset + batch * per_batch, per_batch, p, seed) for batch in range(batches)]
    if workers == 1:
        return [_run_batch(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(_run_batch, tasks))


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
    parser.add_argument("--seed", type=int, default=25050510120261250)
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
    if args.production_manifest is None:
        raise ValueError("a production manifest is required; no unfrozen smoke is permitted")
    manifest = json.loads(args.production_manifest.read_text())
    cap = validate_manifest(args, manifest)
    rows = run(args.samples, args.batches, args.workers, args.p, args.seed, args.replica_offset, cap)
    batches_path = args.batches_output or args.output.with_suffix(".batches.csv")
    write_batches(batches_path, rows)
    payload = {
        "schema": SCHEMA,
        "status": "fresh_radius5_morphism",
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
        "manifest_runner_commit": manifest.get("runner_commit"),
        "observable": {
            "operator": "C4-gauged projective-leg radius-five morphism shell",
            "shell": [list(point) for point in SHELL],
            "hands": list(HANDS),
            "charges": list(CHARGES),
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
