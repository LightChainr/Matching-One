#!/usr/bin/env python3
"""Pair-only Z5 projective-leg transfer stream for Issue 250."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

from norm5_chiral_fixedp_mc import P_FIXED, contexts, counter_uniform
from z5_charged_threepoint_mc import PARENT_GEOMETRY, covariance_of_mean, parent_anchor_indices
from z5_projective_leg_multiseparation_mc import ProjectiveLegIndex, exact_gate as leg_exact_gate
from z5_charged_threepoint_mc import dft_charges


SEPARATIONS = (1, 2, 3, 4, 5, 6)
HANDS = ("plus", "minus")
CHARGES = (1, 2)
SCHEMA = "matching-one/z5-projective-leg-pair-transfer-response/v1"
DEFAULT_CAP = 2_000


def field_order() -> tuple[str, ...]:
    output = []
    for separation in SEPARATIONS:
        for hand in HANDS:
            for charge in CHARGES:
                for row in ("T", "A"):
                    output.extend(
                        (f"d{separation}_{row}{charge}_{hand}_re",
                         f"d{separation}_{row}{charge}_{hand}_im")
                    )
    return tuple(output)


FIELD_ORDER = field_order()


def translated_points(seed: int, replica: int) -> tuple[int, dict[int, tuple[int, int, int]]]:
    translation_index, base = parent_anchor_indices(seed, replica)
    origin = base[0]
    x0, y0 = PARENT_GEOMETRY.coordinates[origin]
    points = {
        separation: (
            origin,
            PARENT_GEOMETRY.vertex((x0 + separation, y0)),
            PARENT_GEOMETRY.vertex((x0, y0 + separation)),
        )
        for separation in SEPARATIONS
    }
    if any(len(set(row)) != 3 for row in points.values()):
        raise AssertionError("pair-transfer anchors collapsed")
    return translation_index, points


def charged_rows(context, field: Sequence[bool], parent_indices: Sequence[int]):
    active = context.active_from_field(field)
    index = ProjectiveLegIndex(context.geometry, active)
    rows = {}
    for parent_index in parent_indices:
        fiber_values = [
            index.scalar(context.field_to_vertex[5 * parent_index + fiber])
            for fiber in range(5)
        ]
        rows[parent_index] = dft_charges(fiber_values)
    return rows


def pair_rows(rows: Mapping[int, Mapping[int, complex]], anchors: Sequence[int]):
    origin, xpoint, ypoint = (rows[index] for index in anchors)
    output = {}
    for charge in CHARGES:
        conjugate = (-charge) % 5
        xvalue = origin[charge] * xpoint[conjugate]
        yvalue = origin[charge] * ypoint[conjugate]
        output[("T", charge)] = (xvalue + yvalue) / 2.0
        output[("A", charge)] = (xvalue - yvalue) / 2.0
    return output


def exact_gate() -> dict:
    base = leg_exact_gate()
    collapse_failures = 0
    for translation in PARENT_GEOMETRY.coordinates:
        x0, y0 = translation
        for separation in SEPARATIONS:
            anchors = (
                PARENT_GEOMETRY.vertex(translation),
                PARENT_GEOMETRY.vertex((x0 + separation, y0)),
                PARENT_GEOMETRY.vertex((x0, y0 + separation)),
            )
            collapse_failures += int(len(set(anchors)) != 3)
    probe = {
        index: dft_charges((index, 1 - index, 2 * index - 3, 4 - index, index - 2))
        for index in range(3)
    }
    rows = pair_rows(probe, (0, 1, 2))
    finite_failures = sum(
        int(not (value.real == value.real and value.imag == value.imag))
        for value in rows.values()
    )
    passed = base["passed"] and collapse_failures == 0 and finite_failures == 0
    if not passed:
        raise AssertionError("pair-transfer exact gate failed")
    return {
        "projective_leg_gate": base,
        "separations": list(SEPARATIONS),
        "parent_translations_checked": PARENT_GEOMETRY.n,
        "anchor_rows_checked": PARENT_GEOMETRY.n * len(SEPARATIONS),
        "anchor_collapse_failures": collapse_failures,
        "pair_definition": "T_r=(O_r(0)O_-r(dx)+O_r(0)O_-r(dy))/2; A_r uses the axis difference",
        "coordinates": len(FIELD_ORDER),
        "finite_probe_failures": finite_failures,
        "passed": True,
    }


def _run_batch(task: tuple[int, int, int, float, int]) -> dict:
    batch, start, samples, p, seed = task
    plus, minus, _ = contexts()
    sums = {name: 0.0 for name in FIELD_ORDER}
    field_digest = hashlib.sha256()
    translation_digest = hashlib.sha256()
    for replica in range(start, start + samples):
        field = [counter_uniform(seed, replica, site) < p for site in range(325)]
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
    batch_means = [
        [float(row[name]) / int(row["samples"]) for name in FIELD_ORDER]
        for row in rows
    ]
    point = [
        sum(row[index] for row in batch_means) / len(batch_means)
        for index in range(len(FIELD_ORDER))
    ]
    return {
        "order": list(FIELD_ORDER),
        "point": point,
        "covariance_of_mean": covariance_of_mean(batch_means),
        "batch_replicates": len(rows),
    }


def run(
    samples: int,
    batches: int,
    workers: int,
    p: float,
    seed: int,
    replica_offset: int,
    *,
    sample_cap: int = DEFAULT_CAP,
):
    if samples <= 0 or batches < 2 or samples % batches:
        raise ValueError("samples must be divisible by batches>=2")
    if samples > sample_cap:
        raise ValueError(f"pair-transfer run exceeds authorized cap {sample_cap}")
    per_batch = samples // batches
    tasks = [
        (batch, replica_offset + batch * per_batch, per_batch, p, seed)
        for batch in range(batches)
    ]
    if workers == 1:
        rows = [_run_batch(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(_run_batch, tasks))
    return rows, summarize(rows)


def write_batches(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def validate_manifest(args, manifest: Mapping[str, object]) -> int:
    if manifest.get("status") != "authorized_fresh_pair_transfer":
        raise ValueError("pair-transfer manifest is not authorized")
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
            raise ValueError(f"run differs from pair-transfer manifest for {key}")
    if manifest.get("separations") != list(SEPARATIONS):
        raise ValueError("manifest separation list changed")
    return int(manifest["run"]["samples"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=2_000)
    parser.add_argument("--batches", type=int, default=40)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--p", type=float, default=P_FIXED)
    parser.add_argument("--seed", type=int, default=25025033720260831)
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
    sample_cap = DEFAULT_CAP
    if args.production_manifest is not None:
        sample_cap = validate_manifest(args, json.loads(args.production_manifest.read_text()))
    rows, analysis = run(
        args.samples, args.batches, args.workers, args.p, args.seed,
        args.replica_offset, sample_cap=sample_cap,
    )
    batches_path = args.batches_output or args.output.with_suffix(".batches.csv")
    write_batches(batches_path, rows)
    payload = {
        "schema": SCHEMA,
        "status": "fresh_pair_transfer" if args.production_manifest else "pair_transfer_smoke",
        "issues": [250],
        "exact_gate": exact_gate(),
        "run": {
            "samples": args.samples, "batches": args.batches,
            "workers": args.workers, "p": args.p, "seed": args.seed,
            "replica_offset": args.replica_offset,
            "replica_last_exclusive": args.replica_offset + args.samples,
            "batches_output": str(batches_path),
        },
        "observable": {
            "operator": "projective-leg rank-one root membership",
            "separations": list(SEPARATIONS),
            "charges": list(CHARGES),
            "hands": list(HANDS),
            "rows": ["complex axis-average T", "complex axis-difference A"],
            "cubic_fields": [],
        },
        "analysis": analysis,
        "production_manifest": str(args.production_manifest) if args.production_manifest else None,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
