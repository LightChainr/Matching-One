#!/usr/bin/env python3
"""Common-field N112 rho-child acquisition for delta E_top and primitive H4.

The primary observer is the Alexander-even rank coordinate

    E_top = 1_rank0 + 1_rank2.

Its continuum-subtracted three-child vector is transformed to the nontrivial
complex C3 character.  Primitive H4 is retained in the same stream only for a
secondary two-observer ray determinant; it is not rescored as H4/H8/H12.
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

from integer_period_torus import (
    determinant, integer_torus_geometry, matrix_product, matrix_vector,
)
from square_bond_primitive_pilot import classify_bond_mask


MASK64 = (1 << 64) - 1
PARENT_MATRIX = ((8, 4), (0, 7))
EMBEDDINGS = (
    ("2omega", ((1, 0), (0, 2))),
    ("omega_over_2", ((2, 0), (0, 1))),
    ("omega_plus_1_over_2", ((2, 1), (0, 1))),
)
CHILDREN = tuple(
    (name, matrix_product(PARENT_MATRIX, embedding))
    for name, embedding in EMBEDDINGS
)
CHILD_ORDER = tuple(name for name, _ in CHILDREN)
EDGE_COUNT = 224
SMOKE_CAP = 5_000


def splitmix64(value):
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return (value ^ (value >> 31)) & MASK64


def counter_mask(seed, replica, bits):
    output = 0
    replica_key = splitmix64(replica + 0xD1B54A32D192ED03)
    for word in range((bits + 63) // 64):
        value = splitmix64(seed ^ replica_key ^ splitmix64(word + 0x94D049BB133111EB))
        output |= value << (64 * word)
    return output & ((1 << bits) - 1)


def physical_phase(matrix, line):
    lifted = matrix_vector(matrix, line)
    value = complex(lifted[0], lifted[1])
    if not value:
        raise ValueError("zero primitive line")
    return (value / abs(value)) ** 4


def exact_gate():
    geometries = [integer_torus_geometry(matrix) for _, matrix in CHILDREN]
    return {
        "parent_N": abs(determinant(PARENT_MATRIX)),
        "child_order": list(CHILD_ORDER),
        "children": [
            {"id": name, "period_matrix_rows": [list(row) for row in matrix],
             "N": geometry.n, "bonds": len(geometry.primal_edges)}
            for (name, matrix), geometry in zip(CHILDREN, geometries)
        ],
        "common_field": "same counter-derived 224-bit bond mask in child edge order",
        "passed": (
            abs(determinant(PARENT_MATRIX)) == 56
            and all(geometry.n == 112 for geometry in geometries)
            and all(len(geometry.primal_edges) == EDGE_COUNT for geometry in geometries)
        ),
    }


def empty_row(batch, first, samples):
    row = {"batch": batch, "replica_first": first, "samples": samples}
    for name in CHILD_ORDER:
        for category in ("rank0", "rank1", "rank2", "invalid"):
            row[f"{name}_{category}"] = 0
        row[f"{name}_H4_re"] = 0.0
        row[f"{name}_H4_im"] = 0.0
    return row


def run_batch(task):
    batch, start, samples, seed = task
    geometries = [integer_torus_geometry(matrix) for _, matrix in CHILDREN]
    row = empty_row(batch, start, samples)
    digest = hashlib.sha256()
    for replica in range(start, start + samples):
        mask = counter_mask(seed, replica, EDGE_COUNT)
        digest.update(mask.to_bytes(EDGE_COUNT // 8, "little"))
        for (name, matrix), geometry in zip(CHILDREN, geometries):
            category, line = classify_bond_mask(geometry, mask)
            if category == "invariant_failure":
                row[f"{name}_invalid"] += 1
            elif line is None:
                row[f"{name}_{category}"] += 1
            else:
                row[f"{name}_rank1"] += 1
                value = physical_phase(matrix, line)
                row[f"{name}_H4_re"] += value.real
                row[f"{name}_H4_im"] += value.imag
    row["common_field_sha256"] = digest.hexdigest()
    return row


def covariance_of_mean(rows):
    count = len(rows)
    means = [sum(row[j] for row in rows) / count for j in range(len(rows[0]))]
    return [
        [sum((row[i] - means[i]) * (row[j] - means[j]) for row in rows)
         / (count * (count - 1)) for j in range(len(means))]
        for i in range(len(means))
    ]


def summarize(rows):
    batch_vectors = []
    for row in rows:
        samples = int(row["samples"])
        vector = []
        for name in CHILD_ORDER:
            vector.extend([
                (int(row[f"{name}_rank0"]) + int(row[f"{name}_rank2"])) / samples,
                float(row[f"{name}_H4_re"]) / samples,
                float(row[f"{name}_H4_im"]) / samples,
            ])
        batch_vectors.append(vector)
    point = [sum(row[j] for row in batch_vectors) / len(batch_vectors)
             for j in range(len(batch_vectors[0]))]
    invalid = {name: sum(int(row[f"{name}_invalid"]) for row in rows)
               for name in CHILD_ORDER}
    return {
        "coordinate_order": [
            coordinate for name in CHILD_ORDER
            for coordinate in (f"{name}_Etop", f"{name}_H4_re", f"{name}_H4_im")
        ],
        "point": point,
        "full_covariance_of_mean_9x9": covariance_of_mean(batch_vectors),
        "invalid_totals": invalid,
        "all_invariant_failures_zero": all(value == 0 for value in invalid.values()),
    }


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=2_000_000)
    parser.add_argument("--batches", type=int, default=100)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--seed", type=int, default=2672751123001)
    parser.add_argument("--replica-offset", type=int, default=27_500_000_000)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batches-output", type=Path, required=True)
    parser.add_argument("--git-commit", default="unknown")
    parser.add_argument("--environment", default="unknown")
    args = parser.parse_args()
    if args.samples <= 0 or args.batches < 20 or args.samples % args.batches:
        raise ValueError("samples must be divisible by batches>=20")
    gate = exact_gate()
    if not gate["passed"]:
        raise AssertionError("rho-child exact geometry gate failed")
    manifest = None
    if args.samples > SMOKE_CAP:
        if args.manifest is None:
            raise ValueError("production above smoke cap requires frozen manifest")
        manifest = json.loads(args.manifest.read_text())
        expected = manifest["run"]
        for key in ("samples", "batches", "workers", "seed", "replica_offset"):
            if getattr(args, key.replace("replica_offset", "replica_offset")) != expected[key]:
                raise ValueError(f"manifest mismatch: {key}")
    per_batch = args.samples // args.batches
    tasks = [(batch, args.replica_offset + batch * per_batch, per_batch, args.seed)
             for batch in range(args.batches)]
    started = time.perf_counter()
    if args.workers == 1:
        rows = [run_batch(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            rows = list(pool.map(run_batch, tasks))
    rows.sort(key=lambda row: row["batch"])
    batch_hash = write_csv(args.batches_output, rows)
    payload = {
        "schema": "matching-one/p267-rho-child-etop-c3-run/v1",
        "run": {"samples": args.samples, "batches": args.batches,
                "workers": args.workers, "seed": args.seed,
                "replica_offset": args.replica_offset,
                "elapsed_seconds": time.perf_counter() - started,
                "git_commit": args.git_commit, "environment": args.environment,
                "python": sys.version.split()[0], "platform": platform.platform()},
        "geometry_gate": gate,
        "batch_sha256": batch_hash,
        "summary": summarize(rows),
        "observer_boundary": (
            "E_top is the Alexander-even rank coordinate. Primitive H4 is retained only for "
            "the same-stream observer-ray determinant, not harmonic voting."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
