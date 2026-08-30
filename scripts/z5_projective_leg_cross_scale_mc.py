#!/usr/bin/env python3
"""Fresh N505 projective-leg pair row for the P250 cross-scale test."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from math import gcd, isqrt
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

from integer_period_torus import gaussian_integer_torus, integer_torus_geometry, matrix_product, matrix_vector
from norm5_chiral_fixedp_mc import M_MINUS, M_PLUS, P_FIXED, splitmix64
from z5_charged_threepoint_mc import covariance_of_mean, dft_charges
from z5_projective_leg_multiseparation_mc import ProjectiveLegIndex
from z5_projective_leg_pair_transfer_mc import pair_rows


PARENT = (10, 1)
PARENT_ORDER = 101
PARENT_MATRIX = ((10, -1), (1, 10))
CHILD_PLUS = matrix_product(PARENT_MATRIX, M_PLUS)
CHILD_MINUS = matrix_product(PARENT_MATRIX, M_MINUS)
CHILD_ORDER = 505
SEPARATIONS = (1, 2, 3, 4, 5)
HANDS = ("plus", "minus")
CHARGES = (1, 2)
SCHEMA = "matching-one/z5-projective-leg-cross-scale-response/v1"
DEFAULT_CAP = 2_000
MASK64 = (1 << 64) - 1
PARENT_GEOMETRY = gaussian_integer_torus(*PARENT)


def field_order() -> tuple[str, ...]:
    return tuple(
        f"d{separation}_{row}{charge}_{hand}_{part}"
        for separation in SEPARATIONS
        for hand in HANDS
        for charge in CHARGES
        for row in ("T", "A")
        for part in ("re", "im")
    )


FIELD_ORDER = field_order()


class ScaleCoverContext:
    def __init__(self, periods, *, name: str) -> None:
        self.name = name
        self.geometry = integer_torus_geometry(periods, name=name)
        representatives = PARENT_GEOMETRY.coordinates
        if self.geometry.n != CHILD_ORDER or len(representatives) != PARENT_ORDER:
            raise ValueError("wrong N101 to N505 cover")
        self.field_coordinates = []
        self.field_to_vertex = []
        for representative in representatives:
            for fiber in range(5):
                offset = matrix_vector(PARENT_MATRIX, (fiber, 0))
                point = representative[0] + offset[0], representative[1] + offset[1]
                self.field_coordinates.append(point)
                self.field_to_vertex.append(self.geometry.vertex(point))
        if len(set(self.field_to_vertex)) != CHILD_ORDER:
            raise AssertionError("parent/fiber labels do not cover the child")

    def active_from_field(self, field: Sequence[bool]) -> list[bool]:
        active = [False] * CHILD_ORDER
        for field_index, value in enumerate(field):
            active[self.field_to_vertex[field_index]] = bool(value)
        return active


_CONTEXTS: tuple[ScaleCoverContext, ScaleCoverContext] | None = None


def contexts() -> tuple[ScaleCoverContext, ScaleCoverContext]:
    global _CONTEXTS
    if _CONTEXTS is None:
        _CONTEXTS = (
            ScaleCoverContext(CHILD_PLUS, name="N505-parent10+i-times2+i"),
            ScaleCoverContext(CHILD_MINUS, name="N505-parent10+i-times2-i"),
        )
    return _CONTEXTS


def counter_uniform(seed: int, replica: int, site: int) -> float:
    key = (
        seed
        ^ splitmix64(CHILD_ORDER)
        ^ splitmix64(replica + 0xD1B54A32D192ED03)
        ^ splitmix64(site + 0x94D049BB133111EB)
    )
    return (splitmix64(key & MASK64) >> 11) * (2.0**-53)


def translated_points(seed: int, replica: int) -> tuple[int, dict[int, tuple[int, int, int]]]:
    translation_index = splitmix64(seed ^ splitmix64(replica + 0x250505101)) % PARENT_ORDER
    x0, y0 = PARENT_GEOMETRY.coordinates[translation_index]
    origin = PARENT_GEOMETRY.vertex((x0, y0))
    points = {
        separation: (
            origin,
            PARENT_GEOMETRY.vertex((x0 + separation, y0)),
            PARENT_GEOMETRY.vertex((x0, y0 + separation)),
        )
        for separation in SEPARATIONS
    }
    if any(len(set(row)) != 3 for row in points.values()):
        raise AssertionError("cross-scale anchors collapsed")
    return translation_index, points


def charged_rows(context: ScaleCoverContext, field: Sequence[bool], parent_indices: Sequence[int]):
    index = ProjectiveLegIndex(context.geometry, context.active_from_field(field))
    return {
        parent_index: dft_charges([
            index.scalar(context.field_to_vertex[5 * parent_index + fiber])
            for fiber in range(5)
        ])
        for parent_index in parent_indices
    }


def primitive_oblique_candidate_audit() -> dict:
    candidates = []
    for norm in range(66, PARENT_ORDER + 1):
        rows = []
        for b in range(1, isqrt(norm) + 1):
            a2 = norm - b * b
            a = isqrt(a2)
            if a >= b and a * a == a2 and gcd(a, b) == 1:
                rows.append([a, b])
        if rows:
            candidates.append({
                "norm": norm,
                "representatives": rows,
                "strict_axis_separations_before_half_period": int((norm ** 0.5) // 2),
            })
    return {
        "scope": "primitive oblique Gaussian parents above norm 65",
        "candidates": candidates,
        "selected": [PARENT_ORDER, *PARENT],
        "reason": "first candidate with five strict unit-axis separations; norm 100 is real/nonprimitive and collapses the two hands to a reflection pair",
    }


def exact_gate() -> dict:
    plus, minus = contexts()
    projection_failures = 0
    deck_failures = 0
    deck_step = matrix_vector(PARENT_MATRIX, (1, 0))
    for context in (plus, minus):
        for field_index, point in enumerate(context.field_coordinates):
            projection_failures += int(PARENT_GEOMETRY.vertex(point) != field_index // 5)
            translated = point[0] + deck_step[0], point[1] + deck_step[1]
            target = 5 * (field_index // 5) + (field_index % 5 + 1) % 5
            deck_failures += int(context.geometry.vertex(translated) != context.field_to_vertex[target])
    anchor_failures = 0
    displacement_vertices = set()
    for separation in SEPARATIONS:
        row = (
            PARENT_GEOMETRY.vertex((separation, 0)),
            PARENT_GEOMETRY.vertex((-separation, 0)),
            PARENT_GEOMETRY.vertex((0, separation)),
            PARENT_GEOMETRY.vertex((0, -separation)),
        )
        anchor_failures += int(len(set(row)) != 4)
        displacement_vertices.update(row)
    passed = (
        projection_failures == 0 and deck_failures == 0 and anchor_failures == 0
        and len(displacement_vertices) == 4 * len(SEPARATIONS)
    )
    if not passed:
        raise AssertionError("N101 to N505 exact gate failed")
    return {
        "candidate_audit": primitive_oblique_candidate_audit(),
        "parent": list(PARENT),
        "parent_periods": PARENT_MATRIX,
        "parent_order": PARENT_ORDER,
        "plus_child_periods": CHILD_PLUS,
        "minus_child_periods": CHILD_MINUS,
        "child_order": CHILD_ORDER,
        "separations": list(SEPARATIONS),
        "distinct_signed_axis_displacements": len(displacement_vertices),
        "projection_failures": projection_failures,
        "deck_failures": deck_failures,
        "anchor_failures": anchor_failures,
        "coordinates": len(FIELD_ORDER),
        "passed": True,
    }


def _run_batch(task: tuple[int, int, int, float, int]) -> dict:
    batch, start, samples, p, seed = task
    plus, minus = contexts()
    sums = {name: 0.0 for name in FIELD_ORDER}
    field_digest = hashlib.sha256()
    translation_digest = hashlib.sha256()
    for replica in range(start, start + samples):
        field = [counter_uniform(seed, replica, site) < p for site in range(CHILD_ORDER)]
        field_digest.update(bytes(field))
        translation_index, anchors = translated_points(seed, replica)
        translation_digest.update(translation_index.to_bytes(2, "little"))
        unique = sorted({index for row in anchors.values() for index in row})
        for hand, context in (("plus", plus), ("minus", minus)):
            charged = charged_rows(context, field, unique)
            for separation, points in anchors.items():
                for (row, charge), value in pair_rows(charged, points).items():
                    sums[f"d{separation}_{row}{charge}_{hand}_re"] += value.real
                    sums[f"d{separation}_{row}{charge}_{hand}_im"] += value.imag
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
        "point": [sum(row[j] for row in batch_means) / len(batch_means) for j in range(len(FIELD_ORDER))],
        "covariance_of_mean": covariance_of_mean(batch_means),
        "batch_replicates": len(rows),
    }


def run(samples: int, batches: int, workers: int, p: float, seed: int, replica_offset: int, *, cap: int):
    if samples <= 0 or batches < 2 or samples % batches:
        raise ValueError("samples must be divisible by batches>=2")
    if samples > cap:
        raise ValueError(f"cross-scale run exceeds authorized cap {cap}")
    per_batch = samples // batches
    tasks = [(batch, replica_offset + batch * per_batch, per_batch, p, seed) for batch in range(batches)]
    if workers == 1:
        rows = [_run_batch(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(_run_batch, tasks))
    return rows, summarize(rows)


def validate_manifest(args, manifest: Mapping[str, object]) -> int:
    if manifest.get("status") != "authorized_fresh_cross_scale_transfer":
        raise ValueError("cross-scale manifest is not authorized")
    observed = {
        "samples": args.samples, "batches": args.batches, "workers": args.workers,
        "p": args.p, "seed": args.seed, "replica_offset": args.replica_offset,
        "replica_last_exclusive": args.replica_offset + args.samples,
    }
    for key, value in observed.items():
        if manifest["run"].get(key) != value:
            raise ValueError(f"run differs from manifest for {key}")
    if manifest.get("target", {}).get("parent") != list(PARENT):
        raise ValueError("manifest target parent changed")
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
    parser.add_argument("--seed", type=int, default=25050510120260830)
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
    if args.production_manifest:
        cap = validate_manifest(args, json.loads(args.production_manifest.read_text()))
    rows, analysis = run(args.samples, args.batches, args.workers, args.p, args.seed, args.replica_offset, cap=cap)
    batches_path = args.batches_output or args.output.with_suffix(".batches.csv")
    write_batches(batches_path, rows)
    payload = {
        "schema": SCHEMA,
        "status": "fresh_cross_scale_transfer" if args.production_manifest else "cross_scale_smoke",
        "issues": [250],
        "exact_gate": exact_gate(),
        "run": {
            "samples": args.samples, "batches": args.batches, "workers": args.workers,
            "p": args.p, "seed": args.seed, "replica_offset": args.replica_offset,
            "replica_last_exclusive": args.replica_offset + args.samples,
            "batches_output": str(batches_path),
        },
        "observable": {
            "operator": "projective-leg rank-one root membership",
            "parent": list(PARENT), "parent_order": PARENT_ORDER,
            "child_order": CHILD_ORDER, "separations": list(SEPARATIONS),
            "rows": ["complex axis-average T", "complex axis-difference A"],
            "cubic_fields": [],
        },
        "analysis": analysis,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
