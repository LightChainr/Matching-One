#!/usr/bin/env python3
"""Common-field square-bond pilot for the three degree-2 rho children.

The measured row is the named primitive-homology automorphic harmonic

    H4 = E[1_{rank=1} (P ell / |P ell|)^4].

It resolves polarization within rank one and is therefore not a function of
``A_top=P2-P0``.  One counter-derived 2N-bit bond field is shared by all three
index-two child quotients in the frozen order 2*tau, tau/2, (tau+1)/2.
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
from typing import Sequence

from integer_period_torus import Matrix, determinant, integer_torus_geometry, matrix_product, matrix_vector
from square_bond_primitive_pilot import classify_bond_mask


MASK64 = (1 << 64) - 1
PARENT_ID = "pell_Dplus1_N56"
PARENT_MATRIX: Matrix = ((8, 4), (0, 7))
CHILD_EMBEDDINGS: tuple[tuple[str, Matrix], ...] = (
    ("2omega", ((1, 0), (0, 2))),
    ("omega_over_2", ((2, 0), (0, 1))),
    ("omega_plus_1_over_2", ((2, 1), (0, 1))),
)
CHILD_DESIGNS: tuple[tuple[str, Matrix], ...] = tuple(
    (name, matrix_product(PARENT_MATRIX, embedding))
    for name, embedding in CHILD_EMBEDDINGS
)
PRIMARY_ORDER = tuple(
    coordinate
    for name, _ in CHILD_DESIGNS
    for coordinate in (f"{name}_re", f"{name}_im")
)
PRODUCTION_ID = "P267-rho-C3-primitive-H4-N112-v1"
SMOKE_CAP = 5_000


def splitmix64(value: int) -> int:
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return (value ^ (value >> 31)) & MASK64


def counter_mask(seed: int, replica: int, bits: int) -> int:
    """Partition-independent counter field, shared byte-for-byte by children."""

    output = 0
    replica_key = splitmix64(replica + 0xD1B54A32D192ED03)
    for word in range((bits + 63) // 64):
        value = splitmix64(seed ^ replica_key ^ splitmix64(word + 0x94D049BB133111EB))
        output |= value << (64 * word)
    return output & ((1 << bits) - 1)


def physical_phase(matrix: Matrix, line: tuple[int, int], spin: int = 4) -> complex:
    lifted = matrix_vector(matrix, line)
    value = complex(lifted[0], lifted[1])
    if not value:
        raise ValueError("zero winding line")
    return (value / abs(value)) ** spin


def child_gate() -> dict[str, object]:
    parent_n = abs(determinant(PARENT_MATRIX))
    geometries = [integer_torus_geometry(matrix) for _, matrix in CHILD_DESIGNS]
    edge_counts = [len(geometry.primal_edges) for geometry in geometries]
    if parent_n != 56 or any(geometry.n != 112 for geometry in geometries):
        raise AssertionError("frozen parent/child determinants changed")
    if len(set(edge_counts)) != 1 or edge_counts[0] != 224:
        raise AssertionError("three-child common bond field is not 224 bits")
    return {
        "parent_id": PARENT_ID,
        "parent_matrix_rows": [list(row) for row in PARENT_MATRIX],
        "parent_N": parent_n,
        "child_order": [name for name, _ in CHILD_DESIGNS],
        "children": [
            {
                "id": name,
                "embedding_in_parent_basis": [list(row) for row in embedding],
                "period_matrix_rows": [list(row) for row in matrix],
                "N": abs(determinant(matrix)),
                "bonds": edge_count,
            }
            for (name, embedding), (_, matrix), edge_count in zip(
                CHILD_EMBEDDINGS, CHILD_DESIGNS, edge_counts
            )
        ],
        "common_field": "same counter-derived 224-bit vector in deterministic primal-edge order",
        "passed": True,
    }


def tiny_oracle() -> dict[str, object]:
    """Exhaust all masks on the three index-two children of the N=4 control."""

    parent: Matrix = ((2, 1), (0, 2))
    designs = [matrix_product(parent, embedding) for _, embedding in CHILD_EMBEDDINGS]
    geometries = [integer_torus_geometry(matrix) for matrix in designs]
    if any(len(geometry.primal_edges) != 16 for geometry in geometries):
        raise AssertionError("tiny oracle must have 16 bonds per child")
    invalid = [0, 0, 0]
    rank1 = [0, 0, 0]
    harmonic = [0j, 0j, 0j]
    for mask in range(1 << 16):
        for index, (geometry, matrix) in enumerate(zip(geometries, designs)):
            category, line = classify_bond_mask(geometry, mask)
            invalid[index] += int(category == "invariant_failure")
            if line is not None:
                rank1[index] += 1
                harmonic[index] += physical_phase(matrix, line)
    digest_payload = ";".join(
        f"{rank1[index]}:{harmonic[index].real:.12f}:{harmonic[index].imag:.12f}"
        for index in range(3)
    )
    return {
        "parent_matrix_rows": [list(row) for row in parent],
        "configurations": 1 << 16,
        "invalid_counts": invalid,
        "rank1_counts": rank1,
        "harmonic_sums_re_im": [[value.real, value.imag] for value in harmonic],
        "digest_sha256": hashlib.sha256(digest_payload.encode()).hexdigest(),
        "passed": invalid == [0, 0, 0],
    }


def _run_batch(task: tuple[int, int, int, int]) -> dict[str, object]:
    batch, start, samples, seed = task
    geometries = [integer_torus_geometry(matrix) for _, matrix in CHILD_DESIGNS]
    edge_count = len(geometries[0].primal_edges)
    sums = {name: 0.0 for name in PRIMARY_ORDER}
    counts = {
        f"{name}_{category}": 0
        for name, _ in CHILD_DESIGNS
        for category in ("rank0", "rank1", "rank2", "invalid")
    }
    field_digest = hashlib.sha256()
    for replica in range(start, start + samples):
        mask = counter_mask(seed, replica, edge_count)
        field_digest.update(mask.to_bytes((edge_count + 7) // 8, "little"))
        for (name, matrix), geometry in zip(CHILD_DESIGNS, geometries):
            category, line = classify_bond_mask(geometry, mask)
            if category == "invariant_failure":
                counts[f"{name}_invalid"] += 1
            elif line is None:
                counts[f"{name}_{category}"] += 1
            else:
                counts[f"{name}_rank1"] += 1
                value = physical_phase(matrix, line)
                sums[f"{name}_re"] += value.real
                sums[f"{name}_im"] += value.imag
    return {
        "batch": batch,
        "replica_first": start,
        "samples": samples,
        "common_field_sha256": field_digest.hexdigest(),
        **sums,
        **counts,
    }


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


def summarize(batches: Sequence[dict[str, object]]) -> dict[str, object]:
    rows = [
        [float(batch[name]) / int(batch["samples"]) for name in PRIMARY_ORDER]
        for batch in batches
    ]
    point = [sum(row[j] for row in rows) / len(rows) for j in range(len(PRIMARY_ORDER))]
    covariance = covariance_of_mean(rows)
    totals = {
        key: sum(int(batch[key]) for batch in batches)
        for key in batches[0]
        if key.endswith(("_rank0", "_rank1", "_rank2", "_invalid"))
    }
    return {
        "primary_order": list(PRIMARY_ORDER),
        "lattice_H4_point_re_im": point,
        "full_common_field_covariance_6x6": covariance,
        "category_totals": totals,
        "all_invariant_failures_zero": all(
            value == 0 for key, value in totals.items() if key.endswith("_invalid")
        ),
    }


def run(samples: int, batches: int, workers: int, seed: int) -> tuple[list[dict], dict]:
    if samples <= 0 or batches < 2 or samples % batches:
        raise ValueError("samples must be divisible by batches>=2")
    per_batch = samples // batches
    tasks = [(batch, batch * per_batch, per_batch, seed) for batch in range(batches)]
    if workers == 1:
        output = [_run_batch(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            output = list(pool.map(_run_batch, tasks))
    return output, summarize(output)


def write_batches(path: Path, rows: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=2_000)
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--seed", type=int, default=267156112)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--batches-output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps({"child_gate": child_gate(), "tiny_oracle": tiny_oracle()}, indent=2))
        return 0
    if args.output is None or args.batches_output is None:
        raise SystemExit("--output and --batches-output are required")
    status = "engineering_smoke"
    manifest = None
    if args.samples > SMOKE_CAP:
        if args.manifest is None:
            raise ValueError("runs above the smoke cap require a manifest")
        manifest = json.loads(args.manifest.read_text())
        expected = manifest["acquisition"]
        actual = {
            "samples": args.samples,
            "batches": args.batches,
            "workers": args.workers,
            "seed": args.seed,
        }
        if (
            manifest.get("production_id") != PRODUCTION_ID
            or not manifest.get("production_authorized")
            or expected != actual
        ):
            raise ValueError("CLI differs from the authorized frozen acquisition")
        status = "production_under_frozen_manifest"
    batches, summary = run(args.samples, args.batches, args.workers, args.seed)
    payload = {
        "schema": "matching-one/rho-child-primitive-h4-mc/v1",
        "status": status,
        "production_id": PRODUCTION_ID,
        "observable": "delta primitive-homology automorphic H4; baseline scored separately",
        "not_A_top": True,
        "p": "1/2 square-bond",
        "samples": args.samples,
        "batches": args.batches,
        "workers": args.workers,
        "seed": args.seed,
        "child_gate": child_gate(),
        "summary": summary,
        "manifest_runner_commit": manifest.get("runner_commit") if manifest else None,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    write_batches(args.batches_output, batches)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
