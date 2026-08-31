#!/usr/bin/env python3
"""Fixed 20k six-noise-level N112 C3 spectral-energy pilot (#437)."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import json
from math import sqrt
from pathlib import Path
import platform
import sys
import time

from integer_period_torus import integer_torus_geometry
from observer_bandwidth_k_centered_euler import euler_observer_values, k_center
from square_bond_primitive_pilot import classify_bond_mask

MATRICES = (((8, 8), (0, 14)), ((16, 4), (0, 7)), ((16, 12), (0, 7)))
COEFFICIENTS = (1, -31, 310, -1240, 1984, -1024)
BITS = 224
MASK64 = (1 << 64) - 1
MASK = (1 << BITS) - 1
NAMES = tuple(f"l{level}_{name}" for level in range(6) for name in
              ("F_re", "F_im", "energy_re", "energy_im", "euler_energy", "degree5_energy")) + ("unfiltered_variance", "euler_variance")


def splitmix64(value):
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return (value ^ (value >> 31)) & MASK64


def random_mask(seed, replica, domain):
    key = splitmix64(replica) ^ splitmix64(seed) ^ splitmix64(domain + 0xD1B54A32D192ED03)
    return sum(splitmix64(key ^ splitmix64(word)) << (64 * word) for word in range(4)) & MASK


def coupled_masks(seed, replica):
    base, replacement = (random_mask(seed, replica, d) for d in (0, 1))
    retention = MASK
    levels = [base]
    for level in range(1, 6):
        retention &= random_mask(seed, replica, level + 1)
        levels.append((base & retention) | (replacement & (MASK ^ retention)))
    return levels, replacement


def c3_etop(mask, geometries):
    values = []
    for geometry in geometries:
        category, _ = classify_bond_mask(geometry, mask)
        if category == "invariant_failure":
            raise AssertionError("topology invariant failure")
        values.append(int(category in ("rank0", "rank2")))
    return complex((2 * values[0] - values[1] - values[2]) / 6,
                   sqrt(3) * (values[2] - values[1]) / 6)


def degree5(mask):
    return -1 if bin(mask & 31).count("1") % 2 == 0 else 1


def non_low_degree_certificate():
    """A 32-point fifth finite difference, exact in Q(i sqrt(3))."""
    geometries = [integer_torus_geometry(matrix) for matrix in MATRICES]
    free_edges, fixed_edges = (0, 28, 56, 84, 112), (140, 168, 196)
    fixed = sum(1 << edge for edge in fixed_edges)
    real_numerator = imaginary_numerator = 0
    for subset in range(32):
        mask = fixed | sum(1 << edge for j, edge in enumerate(free_edges) if (subset >> j) & 1)
        values = [int(classify_bond_mask(g, mask)[0] in ("rank0", "rank2")) for g in geometries]
        sign = (-1) ** (5 - bin(subset).count("1"))
        real_numerator += sign * (2 * values[0] - values[1] - values[2])
        imaginary_numerator += sign * (values[2] - values[1])
    if (real_numerator, imaginary_numerator) != (-2, 0):
        raise AssertionError("fifth-degree topology witness changed")
    return {"free_edges": free_edges, "fixed_open_edges": fixed_edges, "all_other_edges": "closed",
            "point_checks": 32, "mixed_difference_re": "-1/3", "mixed_difference_im": "0",
            "conclusion": "F is not degree<=4; high-pass self-energy is strictly positive under full-support Bernoulli(1/2)"}


def run_batch(task):
    batch, start, count, seed = task
    geometries = [integer_torus_geometry(matrix) for matrix in MATRICES]
    values, _ = euler_observer_values(3)
    euler = {mask: float(value) for mask, value in k_center(values, 9).items()}
    sums = [0.0] * len(NAMES)
    digest = hashlib.sha256()
    base_cpu = 0.0
    began = time.process_time()
    for replica in range(start, start + count):
        levels, replacement = coupled_masks(seed, replica)
        digest.update(levels[0].to_bytes(28, "little"))
        base_started = time.process_time()
        source = c3_etop(levels[0], geometries)
        independent = c3_etop(replacement, geometries)
        base_cpu += time.process_time() - base_started
        source_euler, source_five = euler[levels[0] & 511], degree5(levels[0])
        for level, mask in enumerate(levels):
            target = source if level == 0 else c3_etop(mask, geometries)
            response = source.conjugate() * target
            row = (target.real, target.imag, response.real, response.imag,
                   source_euler * euler[mask & 511], source_five * degree5(mask))
            for index, value in enumerate(row):
                sums[6 * level + index] += value
        sums[-2] += abs(source - independent) ** 2 / 2
        sums[-1] += (source_euler - euler[replacement & 511]) ** 2 / 2
    return {"batch": batch, "replica_first": start, "samples": count,
            "cpu_seconds": time.process_time() - began, "comparator_classification_cpu_seconds": base_cpu,
            "base_sha256": digest.hexdigest(), **dict(zip(NAMES, sums))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--git-commit", required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    run = manifest["run"]
    if manifest["coefficients"] != list(COEFFICIENTS) or run["samples"] > 20000:
        raise ValueError("frozen six-level pilot contract mismatch")
    if run["samples"] % run["batches"]:
        raise ValueError("equal-size batches required")
    per_batch = run["samples"] // run["batches"]
    tasks = [(b, run["replica_offset"] + b * per_batch, per_batch, run["seed"])
             for b in range(run["batches"])]
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=run["workers"]) as pool:
        rows = list(pool.map(run_batch, tasks))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    batch_path = args.output_dir / "batches.csv"
    with batch_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    report = {"schema": "matching-one/p437-high-pass-run/v1", "run": run,
              "non_low_degree_certificate": non_low_degree_certificate(),
              "git_commit": args.git_commit, "manifest_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
              "batch_sha256": hashlib.sha256(batch_path.read_bytes()).hexdigest(),
              "wall_seconds": time.perf_counter() - started,
              "cpu_seconds": sum(row["cpu_seconds"] for row in rows),
              "comparator_classification_cpu_seconds": sum(row["comparator_classification_cpu_seconds"] for row in rows),
              "coordinate_order": NAMES, "python": sys.version, "platform": platform.platform()}
    (args.output_dir / "run.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: report[key] for key in ("wall_seconds", "cpu_seconds", "batch_sha256")}))


if __name__ == "__main__":
    main()
