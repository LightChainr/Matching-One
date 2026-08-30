#!/usr/bin/env python3
"""Minimal N325 Z5-charged three-anchor local-H4 stream for Issue #250."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
from typing import Sequence

from marked_pivotal_h4_reference import marked_pair
from norm5_chiral_fixedp_mc import (
    PARENT,
    PHASES,
    P_FIXED,
    contexts,
    counter_uniform,
    mapping_gate,
    splitmix64,
)
from integer_period_torus import gaussian_integer_torus


SMOKE_CAP = 20_000
PRODUCTION_ID = "P250-N325-Z5-charged-threepoint-v1"
ANCHOR_OFFSETS = ((0, 0), (1, 0), (0, 1))
PRIMARY_COMPLEX_ORDER = ("C113_plus", "C113_minus", "C122_plus", "C122_minus")
CONJUGATE_COMPLEX_ORDER = ("C244_plus", "C244_minus", "C334_plus", "C334_minus")
NONNEUTRAL_COMPLEX_ORDER = ("C111_plus", "C111_minus", "C112_plus", "C112_minus")
JOINT_COMPLEX_ORDER = PRIMARY_COMPLEX_ORDER + CONJUGATE_COMPLEX_ORDER + NONNEUTRAL_COMPLEX_ORDER
PRIMARY_REAL_ORDER = tuple(
    coordinate for name in PRIMARY_COMPLEX_ORDER for coordinate in (f"{name}_re", f"{name}_im")
)
JOINT_REAL_ORDER = tuple(
    coordinate for name in JOINT_COMPLEX_ORDER for coordinate in (f"{name}_re", f"{name}_im")
)
PARENT_GEOMETRY = gaussian_integer_torus(*PARENT)


def covariance_of_mean(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    count = len(rows)
    means = [sum(row[j] for row in rows) / count for j in range(len(rows[0]))]
    return [
        [
            sum((row[i] - means[i]) * (row[j] - means[j]) for row in rows)
            / (count * (count - 1))
            for j in range(len(means))
        ]
        for i in range(len(means))
    ]


def parent_anchor_indices(seed: int, replica: int) -> tuple[int, tuple[int, int, int]]:
    translation_index = splitmix64(seed ^ splitmix64(replica + 0x250113122)) % PARENT_GEOMETRY.n
    translation = PARENT_GEOMETRY.coordinates[translation_index]
    anchors = tuple(
        PARENT_GEOMETRY.vertex((translation[0] + offset[0], translation[1] + offset[1]))
        for offset in ANCHOR_OFFSETS
    )
    if len(set(anchors)) != 3:
        raise AssertionError("three-anchor pattern collapsed")
    return translation_index, anchors


def local_h4_at_root(context, active: Sequence[bool], root_field: int, radius: int) -> int:
    shifted = context.shifted_to_root(active, root_field)
    pair = marked_pair(context.geometry, shifted, radius)
    return pair["primal"]["h4"] - pair["matching"]["h4"]


def dft_charges(real_fiber_values: Sequence[int]) -> dict[int, complex]:
    if len(real_fiber_values) != 5:
        raise ValueError("Z5 local row requires exactly five fiber values")
    return {
        charge: sum(
            complex(*PHASES[(-charge * fiber) % 5]) * real_fiber_values[fiber]
            for fiber in range(5)
        ) / 5.0
        for charge in range(1, 5)
    }


def hand_channels(context, field: Sequence[bool], anchors: Sequence[int], radius: int) -> dict[str, complex]:
    active = context.active_from_field(field)
    rows = []
    for parent_index in anchors:
        local = [
            local_h4_at_root(context, active, 5 * parent_index + fiber, radius)
            for fiber in range(5)
        ]
        rows.append(dft_charges(local))
    return {
        "C113": rows[0][1] * rows[1][1] * rows[2][3],
        "C122": rows[0][1] * rows[1][2] * rows[2][2],
        "C244": rows[0][4] * rows[1][4] * rows[2][2],
        "C334": rows[0][4] * rows[1][3] * rows[2][3],
        "C111": rows[0][1] * rows[1][1] * rows[2][1],
        "C112": rows[0][1] * rows[1][1] * rows[2][2],
    }


def exact_mapping_gate() -> dict:
    gate = mapping_gate()
    plus, minus, _ = contexts()
    projection_failures = 0
    anchor_failures = 0
    for translation_index, translation in enumerate(PARENT_GEOMETRY.coordinates):
        anchors = tuple(
            PARENT_GEOMETRY.vertex((translation[0] + offset[0], translation[1] + offset[1]))
            for offset in ANCHOR_OFFSETS
        )
        anchor_failures += int(len(set(anchors)) != 3)
        for context in (plus, minus):
            for parent_index in anchors:
                for fiber in range(5):
                    field = 5 * parent_index + fiber
                    point = context.field_coordinates[field]
                    projection_failures += int(PARENT_GEOMETRY.vertex(point) != parent_index)
    probe = dft_charges((0, 1, -2, 3, -4))
    conjugacy_residual = max(abs(probe[4] - probe[1].conjugate()), abs(probe[3] - probe[2].conjugate()))
    if projection_failures or anchor_failures or conjugacy_residual > 1e-14:
        raise AssertionError("charged three-anchor exact gate failed")
    return {
        "base_cover_gate": gate,
        "anchor_offsets": [list(row) for row in ANCHOR_OFFSETS],
        "parent_translations_checked": PARENT_GEOMETRY.n,
        "root_labels_checked": PARENT_GEOMETRY.n * 3 * 5 * 2,
        "projection_failures": projection_failures,
        "anchor_failures": anchor_failures,
        "DFT_convention": "O_r=(1/5) sum_k zeta5^(-r k) local_H4(k)",
        "DFT_probe_conjugacy_residual": conjugacy_residual,
        "charges": {
            "neutral": {"A": [1, 1, 3], "B": [1, 2, 2]},
            "conjugates": {"Abar": [4, 4, 2], "Bbar": [4, 3, 3]},
            "nonneutral_controls": {"N111": [1, 1, 1], "N112": [1, 1, 2]},
        },
        "passed": True,
    }


def _run_batch(task: tuple[int, int, int, float, int, int]) -> dict:
    batch, start, samples, p, seed, radius = task
    plus, minus, _ = contexts()
    sums = {name: 0.0 for name in JOINT_REAL_ORDER}
    field_digest = hashlib.sha256()
    translation_digest = hashlib.sha256()
    conjugacy_max = 0.0
    for replica in range(start, start + samples):
        field = [counter_uniform(seed, replica, site) < p for site in range(325)]
        field_digest.update(bytes(field))
        translation_index, anchors = parent_anchor_indices(seed, replica)
        translation_digest.update(translation_index.to_bytes(2, "little"))
        for hand, context in (("plus", plus), ("minus", minus)):
            channels = hand_channels(context, field, anchors, radius)
            conjugacy_max = max(
                conjugacy_max,
                abs(channels["C244"] - channels["C113"].conjugate()),
                abs(channels["C334"] - channels["C122"].conjugate()),
            )
            for channel in ("C113", "C122", "C244", "C334", "C111", "C112"):
                value = channels[channel]
                sums[f"{channel}_{hand}_re"] += value.real
                sums[f"{channel}_{hand}_im"] += value.imag
    return {
        "batch": batch,
        "replica_first": start,
        "samples": samples,
        "field_sha256": field_digest.hexdigest(),
        "translation_sha256": translation_digest.hexdigest(),
        "conjugacy_max_abs": conjugacy_max,
        **sums,
    }


def _complex_from_row(row: dict, name: str, samples: int) -> complex:
    return complex(row[f"{name}_re"], row[f"{name}_im"]) / samples


def summarize(batches: Sequence[dict]) -> dict:
    joint_rows = [
        [batch[name] / batch["samples"] for name in JOINT_REAL_ORDER]
        for batch in batches
    ]
    joint_point = [sum(row[j] for row in joint_rows) / len(joint_rows) for j in range(len(JOINT_REAL_ORDER))]
    joint_covariance = covariance_of_mean(joint_rows)
    primary_indices = [JOINT_REAL_ORDER.index(name) for name in PRIMARY_REAL_ORDER]
    primary_point = [joint_point[index] for index in primary_indices]
    primary_covariance = [
        [joint_covariance[i][j] for j in primary_indices] for i in primary_indices
    ]

    total_samples = sum(batch["samples"] for batch in batches)
    totals = {name: sum(batch[name] for batch in batches) for name in PRIMARY_REAL_ORDER}
    closure_delete_one = []
    for batch in batches:
        kept = total_samples - batch["samples"]
        values = {
            name: complex(
                (totals[f"{name}_re"] - batch[f"{name}_re"]) / kept,
                (totals[f"{name}_im"] - batch[f"{name}_im"]) / kept,
            )
            for name in PRIMARY_COMPLEX_ORDER
        }
        closure_delete_one.append(
            values["C113_plus"] * values["C122_minus"]
            - values["C113_minus"] * values["C122_plus"]
        )
    means = [
        complex(primary_point[2*j], primary_point[2*j+1])
        for j in range(len(PRIMARY_COMPLEX_ORDER))
    ]
    closure = means[0] * means[3] - means[1] * means[2]
    jackknife_mean = sum(closure_delete_one) / len(closure_delete_one)
    factor = (len(batches) - 1) / len(batches)
    closure_covariance = [
        [
            factor * sum(
                ((value.real, value.imag)[i] - (jackknife_mean.real, jackknife_mean.imag)[i])
                * ((value.real, value.imag)[j] - (jackknife_mean.real, jackknife_mean.imag)[j])
                for value in closure_delete_one
            )
            for j in range(2)
        ]
        for i in range(2)
    ]
    nonneutral_indices = [
        JOINT_REAL_ORDER.index(f"{name}_{part}")
        for name in NONNEUTRAL_COMPLEX_ORDER
        for part in ("re", "im")
    ]
    return {
        "primary_order": list(PRIMARY_REAL_ORDER),
        "primary_point": primary_point,
        "primary_covariance_of_mean": primary_covariance,
        "joint_order": list(JOINT_REAL_ORDER),
        "joint_point": joint_point,
        "joint_covariance_of_mean": joint_covariance,
        "closure": {
            "relation": "C113_plus*C122_minus-C113_minus*C122_plus",
            "point_re_im": [closure.real, closure.imag],
            "delete_one_covariance_re_im": closure_covariance,
            "delete_one_replicates": len(closure_delete_one),
        },
        "nonneutral_controls": {
            "order": [JOINT_REAL_ORDER[index] for index in nonneutral_indices],
            "point": [joint_point[index] for index in nonneutral_indices],
            "covariance_of_mean": [
                [joint_covariance[i][j] for j in nonneutral_indices]
                for i in nonneutral_indices
            ],
        },
        "conjugacy_max_abs": max(batch["conjugacy_max_abs"] for batch in batches),
        "conjugacy_relations": [
            "C244_hand=conj(C113_hand)", "C334_hand=conj(C122_hand)"
        ],
    }


def run(samples: int, batches: int, workers: int, p: float, seed: int, radius: int, replica_offset: int):
    if samples <= 0 or batches < 2 or samples % batches:
        raise ValueError("samples must be divisible by batches>=2")
    per_batch = samples // batches
    tasks = [
        (batch, replica_offset + batch * per_batch, per_batch, p, seed, radius)
        for batch in range(batches)
    ]
    if workers == 1:
        rows = [_run_batch(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(_run_batch, tasks))
    return rows, summarize(rows)


def write_batches(path: Path, rows: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def validate_gate(args, batches_path: Path) -> tuple[str, dict | None]:
    if args.samples <= SMOKE_CAP:
        if args.production_manifest is not None:
            raise ValueError("a production manifest cannot label a smoke run")
        return "engineering_variance_smoke", None
    if args.production_manifest is None:
        raise ValueError("runs above 20k require an authorized frozen manifest")
    manifest = json.loads(args.production_manifest.read_text())
    actual = {
        "samples": args.samples,
        "batches": args.batches,
        "workers": args.workers,
        "p": args.p,
        "seed": args.seed,
        "radius": args.radius,
        "replica_offset": args.replica_offset,
        "output": str(args.output),
        "batches_output": str(batches_path),
    }
    if (
        manifest.get("schema") != "matching-one/z5-charged-threepoint-production/v1"
        or manifest.get("production_id") != PRODUCTION_ID
        or not manifest.get("production_authorized")
        or manifest.get("run") != actual
    ):
        raise ValueError("CLI differs from the authorized frozen acquisition")
    return "production_under_frozen_manifest", manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=2_000)
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--p", type=float, default=P_FIXED)
    parser.add_argument("--seed", type=int, default=25011312220260830)
    parser.add_argument("--radius", type=int, default=1)
    parser.add_argument("--replica-offset", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batches-output", type=Path)
    parser.add_argument("--production-manifest", type=Path)
    parser.add_argument("--exact-gate", action="store_true")
    args = parser.parse_args()
    if args.exact_gate:
        args.output.write_text(json.dumps(exact_mapping_gate(), indent=2) + "\n")
        return 0
    batches_path = args.batches_output or args.output.with_suffix(".batches.csv")
    status, manifest = validate_gate(args, batches_path)
    gate = exact_mapping_gate()
    rows, analysis = run(
        args.samples, args.batches, args.workers, args.p, args.seed, args.radius,
        args.replica_offset,
    )
    write_batches(batches_path, rows)
    payload = {
        "schema": "matching-one/z5-charged-threepoint-response/v1",
        "status": status,
        "production_id": manifest.get("production_id") if manifest else PRODUCTION_ID,
        "manifest_runner_commit": manifest.get("runner_commit") if manifest else None,
        "issues": [250],
        "mapping_gate": gate,
        "run": {
            "samples": args.samples,
            "batches": args.batches,
            "workers": args.workers,
            "p": args.p,
            "seed": args.seed,
            "radius": args.radius,
            "replica_offset": args.replica_offset,
            "replica_last_exclusive": args.replica_offset + args.samples,
            "batches_output": str(batches_path),
        },
        "observable": {
            "anchors": [list(row) for row in ANCHOR_OFFSETS],
            "anchor_translation": "counter-derived uniform parent translation per replica",
            "local_row": "primal-minus-matching landing-pivotal H4 at all five fibers",
            "fiber_DFT": "O_r=(1/5) sum_k zeta5^(-r k) local_H4(k)",
            "primary_channels": list(PRIMARY_COMPLEX_ORDER),
            "local_H4_evaluations_per_replica": 30,
            "same_field_hands": ["(8+i)(2+i)", "(8+i)(2-i)"],
        },
        "analysis": analysis,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
