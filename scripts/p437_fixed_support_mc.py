#!/usr/bin/env python3
"""Fixed five-bond positive topology-energy endpoint, fresh local 20k."""
from __future__ import annotations
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time

from integer_period_torus import integer_torus_geometry
from p437_high_pass_mc import MATRICES, random_mask
from square_bond_primitive_pilot import classify_bond_mask

SUPPORT = (0, 28, 56, 84, 112)
SUPPORT_MASK = sum(1 << site for site in SUPPORT)
VERTICES = tuple((sum(1 << site for j, site in enumerate(SUPPORT) if (subset >> j) & 1),
                  (-1) ** (5 - bin(subset).count("1"))) for subset in range(32))
DENOMINATOR = 192 ** 2


def child_differences(outside_mask, geometries):
    base = outside_mask & ~SUPPORT_MASK
    sums = [0, 0, 0]
    for vertex, sign in VERTICES:
        mask = base | vertex
        for child, geometry in enumerate(geometries):
            category, _ = classify_bond_mask(geometry, mask)
            if category == "invariant_failure":
                raise AssertionError("topology invariant failure")
            sums[child] += sign * int(category in ("rank0", "rank2"))
    return tuple(sums)


def energy_numerator(differences):
    a, b, c = differences
    real, imaginary = 2 * a - b - c, c - b
    return real * real + 3 * imaginary * imaginary


def run_batch(task):
    batch, first, count, seed = task
    geometries = [integer_torus_geometry(matrix) for matrix in MATRICES]
    classes = Counter()
    digest = hashlib.sha256()
    began = time.process_time()
    for replica in range(first, first + count):
        mask = random_mask(seed, replica, 0) & ~SUPPORT_MASK
        digest.update(mask.to_bytes(28, "little"))
        classes[child_differences(mask, geometries)] += 1
    return {"batch": batch, "replica_first": first, "samples": count,
            "cpu_seconds": time.process_time() - began, "outside_sha256": digest.hexdigest(),
            "nonzero": sum(count for value, count in classes.items() if energy_numerator(value)),
            "energy_numerator_sum": sum(energy_numerator(value) * count for value, count in classes.items()),
            "classes": [{"child_difference_numerators": value, "count": count,
                         "energy_numerator": energy_numerator(value)} for value, count in sorted(classes.items())]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    run = manifest["run"]
    if manifest["support"] != list(SUPPORT) or run["samples"] != 20000 or run["batches"] != 100:
        raise ValueError("frozen support/size mismatch")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    per_batch = run["samples"] // run["batches"]
    tasks = [(b, run["replica_offset"] + b * per_batch, per_batch, run["seed"]) for b in range(run["batches"])]
    began = time.perf_counter()
    with ProcessPoolExecutor(max_workers=run["workers"]) as pool:
        rows = list(pool.map(run_batch, tasks))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    batch_path = args.output_dir / "batches.json"
    batch_path.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")
    metadata = {"schema": "matching-one/p437-fixed-support-run/v1", "run": run,
                "git_commit": revision, "platform": platform.platform(),
                "wall_seconds": time.perf_counter() - began, "cpu_seconds": sum(row["cpu_seconds"] for row in rows),
                "batch_sha256": hashlib.sha256(batch_path.read_bytes()).hexdigest(),
                "manifest_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
                "support": SUPPORT, "energy_denominator": DENOMINATOR}
    (args.output_dir / "run.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metadata))


if __name__ == "__main__":
    main()
