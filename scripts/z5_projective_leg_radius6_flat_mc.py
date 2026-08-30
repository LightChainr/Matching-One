#!/usr/bin/env python3
"""Frozen minimal degree-six acquisition for the P250 R2 flat-extension test."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

from norm5_chiral_fixedp_mc import P_FIXED, splitmix64
from z5_projective_leg_bivariate_mc import gauge_charged_rows, label, pair_value, rotation_gauges
from z5_projective_leg_cross_scale_mc import (
    CHILD_ORDER,
    CHARGES,
    PARENT_GEOMETRY,
    charged_rows,
    contexts,
    counter_uniform,
)


DEGREE3 = ((3, 0), (2, 1), (1, 2), (0, 3))
SOURCE_ENDPOINTS = tuple((6 - index, index) for index in range(7))


def selected_r2_alexander(point: tuple[int, int]) -> tuple[int, int]:
    """Alexander y-reflection followed by the selected R2 half-turn."""
    return -point[0], point[1]


TARGET_ENDPOINTS = tuple(selected_r2_alexander(point) for point in SOURCE_ENDPOINTS)
POINTS_BY_HAND = {"plus": SOURCE_ENDPOINTS, "minus": TARGET_ENDPOINTS}
ACQUIRED_ENDPOINTS = tuple(sorted(set(SOURCE_ENDPOINTS) | set(TARGET_ENDPOINTS)))
SCHEMA = "matching-one/z5-projective-leg-radius6-flat-response/v1"
DEFAULT_CAP = 2_000


def add(first: tuple[int, int], second: tuple[int, int]) -> tuple[int, int]:
    return first[0] + second[0], first[1] + second[1]


def field_order() -> tuple[str, ...]:
    return tuple(
        f"{label(a, b)}_r{charge}_{hand}_{part}"
        for hand, points in POINTS_BY_HAND.items()
        for a, b in points
        for charge in CHARGES
        for part in ("re", "im")
    )


FIELD_ORDER = field_order()


def exact_gate() -> dict:
    degree_six_sums = {add(first, second) for first in DEGREE3 for second in DEGREE3}
    vertices = {point: PARENT_GEOMETRY.vertex(point) for point in ACQUIRED_ENDPOINTS}
    origin = PARENT_GEOMETRY.vertex((0, 0))
    full_shell = {
        (a, b) for a in range(-6, 7) for b in range(-6, 7) if abs(a) + abs(b) == 6
    }
    passed = (
        degree_six_sums == set(SOURCE_ENDPOINTS)
        and {selected_r2_alexander(point) for point in SOURCE_ENDPOINTS} == set(TARGET_ENDPOINTS)
        and all(selected_r2_alexander(selected_r2_alexander(point)) == point for point in ACQUIRED_ENDPOINTS)
        and len(ACQUIRED_ENDPOINTS) == 13
        and len(vertices) == len(set(vertices.values()))
        and origin not in vertices.values()
        and set(ACQUIRED_ENDPOINTS) < full_shell
    )
    if not passed:
        raise AssertionError("radius-six flat-extension geometry gate failed")
    return {
        "selected_gauge": "Alexander reflection followed by R2 and coefficient conjugation",
        "spatial_map": "phi(a,b)=(-a,b)",
        "degree3_monomials": [list(point) for point in DEGREE3],
        "source_degree6_endpoints": [list(point) for point in SOURCE_ENDPOINTS],
        "target_degree6_endpoints": [list(point) for point in TARGET_ENDPOINTS],
        "acquired_endpoint_union": [list(point) for point in ACQUIRED_ENDPOINTS],
        "acquired_spatial_points": len(ACQUIRED_ENDPOINTS),
        "full_radius6_shell_points": len(full_shell),
        "omitted_full_shell_points": len(full_shell) - len(ACQUIRED_ENDPOINTS),
        "distinct_parent_vertices": len(set(vertices.values())),
        "source_is_exact_degree3_plus_degree3_sumset": degree_six_sums == set(SOURCE_ENDPOINTS),
        "selected_map_is_involution": True,
        "hands": {hand: [list(point) for point in points] for hand, points in POINTS_BY_HAND.items()},
        "charges": list(CHARGES),
        "complex_coordinates_per_batch": sum(len(points) for points in POINTS_BY_HAND.values()) * len(CHARGES),
        "real_coordinates_per_batch": len(FIELD_ORDER),
        "minimality": "Every acquired endpoint is a distinct entry of the missing degree3-by-degree3 Hankel block in the fixed plus/R2-minus gauges; the other 11 radius-six shell points do not enter either block.",
        "passed": True,
    }


def translated_points(seed: int, replica: int) -> tuple[int, int, dict[tuple[int, int], int]]:
    translation_index = splitmix64(seed ^ splitmix64(replica + 0x250606101260)) % PARENT_GEOMETRY.n
    x0, y0 = PARENT_GEOMETRY.coordinates[translation_index]
    origin = PARENT_GEOMETRY.vertex((x0, y0))
    points = {
        displacement: PARENT_GEOMETRY.vertex((x0 + displacement[0], y0 + displacement[1]))
        for displacement in ACQUIRED_ENDPOINTS
    }
    if len(set(points.values())) != len(points) or origin in points.values():
        raise AssertionError("radius-six translated endpoint set collapsed")
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
            for displacement in POINTS_BY_HAND[hand]:
                target = points[displacement]
                for charge in CHARGES:
                    value = pair_value(rows, origin, target, charge)
                    prefix = f"{label(*displacement)}_r{charge}_{hand}_"
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
    if not manifest.get("execution_authorized") or manifest.get("status") != "authorized_fresh_radius6_flat_extension":
        raise ValueError("radius-six flat-extension execution is not authorized")
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
    if manifest.get("geometry", {}).get("points_by_hand") != {
        hand: [list(point) for point in points] for hand, points in POINTS_BY_HAND.items()
    }:
        raise ValueError("radius-six endpoint contract changed")
    return args.samples


def run(samples: int, batches: int, workers: int, p: float, seed: int, replica_offset: int, cap: int):
    if samples <= 0 or batches < 2 or samples % batches:
        raise ValueError("samples must be divisible by batches>=2")
    if samples > cap:
        raise ValueError(f"radius-six run exceeds authorized cap {cap}")
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
    parser.add_argument("--seed", type=int, default=25060610120261250)
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
        "status": "fresh_radius6_flat_extension",
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
            "operator": "minimal R2-gauged projective-leg degree-six flat-extension endpoints",
            "points_by_hand": {hand: [list(point) for point in points] for hand, points in POINTS_BY_HAND.items()},
            "charges": list(CHARGES),
        },
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
