#!/usr/bin/env python3
"""Small N325 multi-separation charged-cubic stream for Issue 250."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

from norm5_chiral_fixedp_mc import P_FIXED, contexts, counter_uniform, mapping_gate
from z5_charged_threepoint_mc import (
    PARENT_GEOMETRY,
    covariance_of_mean,
    dft_charges,
    local_h4_at_root,
    parent_anchor_indices,
)


SEPARATIONS = (1, 2, 3)
HANDS = ("plus", "minus")
CHANNELS = ("C113", "C122", "C244", "C334")
PAIR_ROWS = ("G1", "G2", "V1", "V2")
SMOKE_CAP = 5_000
SCHEMA = "matching-one/z5-charged-multiseparation-response/v1"


def field_order() -> tuple[str, ...]:
    rows = []
    for separation in SEPARATIONS:
        for hand in HANDS:
            for channel in CHANNELS:
                rows.extend(
                    (f"d{separation}_{channel}_{hand}_re", f"d{separation}_{channel}_{hand}_im")
                )
            rows.extend(f"d{separation}_{name}_{hand}" for name in PAIR_ROWS)
    return tuple(rows)


FIELD_ORDER = field_order()


def translated_anchors(seed: int, replica: int) -> tuple[int, dict[int, tuple[int, int, int]]]:
    translation_index, base = parent_anchor_indices(seed, replica)
    origin = base[0]
    x0, y0 = PARENT_GEOMETRY.coordinates[origin]
    anchors = {
        d: (
            origin,
            PARENT_GEOMETRY.vertex((x0 + d, y0)),
            PARENT_GEOMETRY.vertex((x0, y0 + d)),
        )
        for d in SEPARATIONS
    }
    if any(len(set(row)) != 3 for row in anchors.values()):
        raise AssertionError("multi-separation anchor triangle collapsed")
    return translation_index, anchors


def charged_rows(context, field: Sequence[bool], parent_indices: Sequence[int], radius: int):
    active = context.active_from_field(field)
    output = {}
    for parent_index in parent_indices:
        local = [
            local_h4_at_root(context, active, 5 * parent_index + fiber, radius)
            for fiber in range(5)
        ]
        output[parent_index] = dft_charges(local)
    return output


def hermitian_pair(origin: Mapping[int, complex], xrow: Mapping[int, complex],
                   yrow: Mapping[int, complex], charge: int) -> complex:
    conjugate = (-charge) % 5
    return (
        origin[charge] * xrow[conjugate]
        + origin[conjugate] * xrow[charge]
        + origin[charge] * yrow[conjugate]
        + origin[conjugate] * yrow[charge]
    ) / 4.0


def observables_for_hand(rows: Mapping[int, Mapping[int, complex]], anchors: Sequence[int]):
    origin, xpoint, ypoint = (rows[index] for index in anchors)
    channels = {
        "C113": origin[1] * xpoint[1] * ypoint[3],
        "C122": origin[1] * xpoint[2] * ypoint[2],
        "C244": origin[4] * xpoint[4] * ypoint[2],
        "C334": origin[4] * xpoint[3] * ypoint[3],
    }
    pair1 = hermitian_pair(origin, xpoint, ypoint, 1)
    pair2 = hermitian_pair(origin, xpoint, ypoint, 2)
    variance1 = sum(abs(row[1]) ** 2 for row in (origin, xpoint, ypoint)) / 3.0
    variance2 = sum(abs(row[2]) ** 2 for row in (origin, xpoint, ypoint)) / 3.0
    return channels, {"G1": pair1, "G2": pair2, "V1": variance1, "V2": variance2}


def exact_mapping_gate() -> dict:
    plus, minus, _ = contexts()
    base_gate = mapping_gate()
    collapse_failures = 0
    projection_failures = 0
    for translation in PARENT_GEOMETRY.coordinates:
        x0, y0 = translation
        for separation in SEPARATIONS:
            anchors = (
                PARENT_GEOMETRY.vertex(translation),
                PARENT_GEOMETRY.vertex((x0 + separation, y0)),
                PARENT_GEOMETRY.vertex((x0, y0 + separation)),
            )
            collapse_failures += int(len(set(anchors)) != 3)
            for context in (plus, minus):
                for parent_index in anchors:
                    for fiber in range(5):
                        child_point = context.field_coordinates[5 * parent_index + fiber]
                        projection_failures += int(
                            PARENT_GEOMETRY.vertex(child_point) != parent_index
                        )
    if collapse_failures or projection_failures:
        raise AssertionError("multi-separation mapping gate failed")
    return {
        "base_cover_gate": base_gate,
        "separations": list(SEPARATIONS),
        "anchor_rule": "(0,0),(d,0),(0,d) after uniform parent translation",
        "translations_checked": PARENT_GEOMETRY.n,
        "root_labels_checked": PARENT_GEOMETRY.n * len(SEPARATIONS) * 3 * 5 * 2,
        "collapse_failures": collapse_failures,
        "projection_failures": projection_failures,
        "charged_connectedness": (
            "113 and 122 have no neutral proper charge subset; their raw third moments "
            "are connected cumulants by exact Z5 charge conservation"
        ),
        "passed": True,
    }


def _run_batch(task: tuple[int, int, int, float, int, int]) -> dict:
    batch, start, samples, p, seed, radius = task
    plus, minus, _ = contexts()
    sums = {name: 0.0 for name in FIELD_ORDER}
    field_digest = hashlib.sha256()
    translation_digest = hashlib.sha256()
    pair_imaginary_max = 0.0
    conjugacy_max = 0.0
    for replica in range(start, start + samples):
        field = [counter_uniform(seed, replica, site) < p for site in range(325)]
        field_digest.update(bytes(field))
        translation_index, triangles = translated_anchors(seed, replica)
        translation_digest.update(translation_index.to_bytes(2, "little"))
        unique = sorted({index for triangle in triangles.values() for index in triangle})
        for hand, context in (("plus", plus), ("minus", minus)):
            rows = charged_rows(context, field, unique, radius)
            for separation, anchors in triangles.items():
                channels, pairs = observables_for_hand(rows, anchors)
                conjugacy_max = max(
                    conjugacy_max,
                    abs(channels["C244"] - channels["C113"].conjugate()),
                    abs(channels["C334"] - channels["C122"].conjugate()),
                )
                for channel, value in channels.items():
                    sums[f"d{separation}_{channel}_{hand}_re"] += value.real
                    sums[f"d{separation}_{channel}_{hand}_im"] += value.imag
                for name, value in pairs.items():
                    if isinstance(value, complex):
                        pair_imaginary_max = max(pair_imaginary_max, abs(value.imag))
                        value = value.real
                    sums[f"d{separation}_{name}_{hand}"] += value
    return {
        "batch": batch,
        "replica_first": start,
        "samples": samples,
        "field_sha256": field_digest.hexdigest(),
        "translation_sha256": translation_digest.hexdigest(),
        "pair_imaginary_max": pair_imaginary_max,
        "conjugacy_max_abs": conjugacy_max,
        **sums,
    }


def summarize(rows: Sequence[dict]) -> dict:
    batch_means = [[row[name] / row["samples"] for name in FIELD_ORDER] for row in rows]
    point = [sum(row[index] for row in batch_means) / len(batch_means) for index in range(len(FIELD_ORDER))]
    return {
        "order": list(FIELD_ORDER),
        "point": point,
        "covariance_of_mean": covariance_of_mean(batch_means),
        "batch_replicates": len(rows),
        "pair_imaginary_max": max(row["pair_imaginary_max"] for row in rows),
        "DFT_conjugacy_max_abs": max(row["conjugacy_max_abs"] for row in rows),
    }


def run(samples: int, batches: int, workers: int, p: float, seed: int, radius: int,
        replica_offset: int):
    if samples <= 0 or batches < 2 or samples % batches:
        raise ValueError("samples must be divisible by batches>=2")
    if samples > SMOKE_CAP:
        raise ValueError("this branch is capped at a 5k model-development smoke")
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=2_000)
    parser.add_argument("--batches", type=int, default=40)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--p", type=float, default=P_FIXED)
    parser.add_argument("--seed", type=int, default=25011312220260901)
    parser.add_argument("--radius", type=int, default=1)
    parser.add_argument("--replica-offset", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batches-output", type=Path)
    parser.add_argument("--exact-gate", action="store_true")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.exact_gate:
        args.output.write_text(json.dumps(exact_mapping_gate(), indent=2) + "\n")
        return 0
    rows, analysis = run(
        args.samples, args.batches, args.workers, args.p, args.seed, args.radius,
        args.replica_offset,
    )
    batches_path = args.batches_output or args.output.with_suffix(".batches.csv")
    write_batches(batches_path, rows)
    payload = {
        "schema": SCHEMA,
        "status": "local_model_development_smoke",
        "issues": [250],
        "mapping_gate": exact_mapping_gate(),
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
            "separations": list(SEPARATIONS),
            "triangles": "(0,0),(d,0),(0,d)",
            "two_point": (
                "G_r(d)=1/4 sum_(axis x,y and reversed charges) "
                "O_r(0) O_-r(d axis); exactly real"
            ),
            "local_variance": "V_r(d)=mean over the three anchors of |O_r|^2",
            "separation_normalization": (
                "Omega113=C113/sqrt(|G1|^2 |G2|); "
                "Omega122=C122/sqrt(|G1| |G2|^2)"
            ),
        },
        "analysis": analysis,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
