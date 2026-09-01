#!/usr/bin/env python3
"""Minimal Phase-A Monte Carlo for the Issue #234 triangular log pair.

The engine generates only Bernoulli site configurations.  Camia--Feng cluster
signs are integrated out analytically by ``triangular_energy_logpair_stats``.
Every independently seeded batch archives the frozen integer sufficient
statistics and yields an unbiased same-stream (LL,LD,DD) U-statistic.
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import json
import math
import os
from pathlib import Path
import platform
import random
import sys
from typing import Iterable, Sequence

from triangular_energy_logpair_stats import (
    MOMENT_ORDER,
    PAIR_ORDER,
    PRODUCT_ORDER,
    SufficientSums,
    black_cluster_roots,
    covariance_of_mean,
    translation_averaged_configuration_sums,
    triangular_edges,
    unbiased_moment_vector,
)


MASK64 = (1 << 64) - 1
DELTA_DENOMINATORS = (8, 12, 16)


@dataclass(frozen=True)
class BatchArchive:
    batch: int
    sufficient: SufficientSums


def splitmix64(value: int) -> int:
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return value ^ (value >> 31)


def delta_value(denominator: int) -> float:
    if denominator not in DELTA_DENOMINATORS:
        raise ValueError(f"delta denominator must be one of {DELTA_DENOMINATORS}")
    return 1.0 / (denominator * math.sqrt(2.0))


def delta_radius(length: int, denominator: int) -> int:
    if length <= 0 or length % 2:
        raise ValueError("L must be a positive even integer")
    scaled = delta_value(denominator) * length
    radius = math.floor(scaled + 0.5)
    if radius <= 1 or 2 * radius >= length:
        raise ValueError("the declared L/delta pair does not resolve a bilocal radius")
    return radius


def _run_batch(task: tuple[int, int, int, int, int]) -> BatchArchive:
    length, denominator, batch, samples, seed = task
    radius = delta_radius(length, denominator)
    vertex_count = length * length
    edges = triangular_edges(length, length)
    displacement = (length // 2, length // 2)
    rng = random.Random(splitmix64(seed ^ splitmix64(batch + 1)))
    sums = SufficientSums.empty()
    for _ in range(samples):
        mask = rng.getrandbits(vertex_count)
        roots = black_cluster_roots(length, length, mask, edges)
        pair_sums, four_sums = translation_averaged_configuration_sums(
            roots,
            length,
            length,
            delta_radius=radius,
            center_displacement=displacement,
        )
        sums.add(pair_sums, four_sums)
    return BatchArchive(batch=batch, sufficient=sums)


def run_batches(
    *,
    length: int,
    denominator: int,
    samples: int,
    batches: int,
    seed: int,
    workers: int,
) -> list[BatchArchive]:
    delta_radius(length, denominator)
    if batches < 2 or samples <= 0 or samples % batches:
        raise ValueError("samples must be positive and divisible by batches>=2")
    if workers <= 0:
        raise ValueError("workers must be positive")
    per_batch = samples // batches
    if per_batch < 2:
        raise ValueError("the unbiased U-statistic needs at least two samples per batch")
    tasks = [
        (length, denominator, batch, per_batch, seed) for batch in range(batches)
    ]
    if workers == 1:
        archives = [_run_batch(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            archives = list(executor.map(_run_batch, tasks))
    return sorted(archives, key=lambda row: row.batch)


def _diagnostic(
    estimate: Sequence[float], covariance: Sequence[Sequence[float]]
) -> dict[str, object]:
    ll, ld, dd = estimate
    determinant = ll * dd - ld * ld
    determinant_gradient = (dd, -2 * ld, ll)
    determinant_variance = sum(
        determinant_gradient[i] * covariance[i][j] * determinant_gradient[j]
        for i in range(3)
        for j in range(3)
    )
    payload: dict[str, object] = {
        "gram_determinant_LL_DD_minus_LD_squared": determinant,
        "gram_determinant_standard_error_delta_method": math.sqrt(
            max(0.0, determinant_variance)
        ),
    }
    determinant_se = float(payload["gram_determinant_standard_error_delta_method"])
    payload["gram_determinant_signed_z"] = (
        determinant / determinant_se if determinant_se else None
    )
    if ld:
        invariant = ll * dd / (ld * ld) - 1.0
        gradient = (
            dd / (ld * ld),
            -2 * ll * dd / (ld**3),
            ll / (ld * ld),
        )
        variance = sum(
            gradient[i] * covariance[i][j] * gradient[j]
            for i in range(3)
            for j in range(3)
        )
        payload.update(
            {
                "rescaling_invariant_J": invariant,
                "rescaling_invariant_J_definition": "LL*DD/LD^2-1",
                "rescaling_invariant_J_standard_error_delta_method": math.sqrt(
                    max(0.0, variance)
                ),
            }
        )
    else:
        payload.update(
            {
                "rescaling_invariant_J": None,
                "rescaling_invariant_J_definition": "LL*DD/LD^2-1",
                "rescaling_invariant_J_standard_error_delta_method": None,
            }
        )
    return payload


def analyze_archives(
    archives: Sequence[BatchArchive], *, length: int, denominator: int
) -> dict[str, object]:
    if len(archives) < 2:
        raise ValueError("at least two batch archives are required")
    samples_per_batch = archives[0].sufficient.samples
    if samples_per_batch < 2 or any(
        row.sufficient.samples != samples_per_batch for row in archives
    ):
        raise ValueError("all batches must have the same sample count >=2")
    placements = length * length
    block_vectors = [
        [float(value) for value in unbiased_moment_vector(row.sufficient, placements=placements)]
        for row in archives
    ]
    estimate = [
        math.fsum(row[index] for row in block_vectors) / len(block_vectors)
        for index in range(3)
    ]
    covariance = covariance_of_mean(block_vectors)
    standard_error = [math.sqrt(max(0.0, covariance[i][i])) for i in range(3)]
    radius = delta_radius(length, denominator)
    delta = delta_value(denominator)
    return {
        "schema": "matching-one/p234-triangular-energy-logpair-phaseA/v1",
        "issue": 234,
        "status": "raw_parent_pair_phaseA_no_universal_coefficient",
        "geometry": {
            "lattice": "triangular_site",
            "p": 0.5,
            "L": length,
            "sites": length * length,
            "a": 1.0 / length,
            "period_basis_lattice_coordinates": [[length, 0], [0, length]],
            "center_displacement": [length // 2, length // 2],
            "translations_per_configuration": placements,
            "delta_formula": f"1/({denominator}*sqrt(2))",
            "delta": delta,
            "bilocal_radius_lattice_units": radius,
            "realized_delta": radius / length,
            "realized_delta_error": radius / length - delta,
        },
        "monte_carlo": {
            "batches": len(archives),
            "samples_per_batch": samples_per_batch,
            "total_samples": len(archives) * samples_per_batch,
            "random_cluster_sign_draws": 0,
            "sign_treatment": "exact conditional integration over cluster signs",
            "centering": "unbiased cross-configuration order-2 U-statistic",
        },
        "moment_order": list(MOMENT_ORDER),
        "estimate": estimate,
        "standard_error": standard_error,
        "covariance_of_mean": covariance,
        "block_estimates": block_vectors,
        "two_field_matrix": [[estimate[0], estimate[1]], [estimate[1], estimate[2]]],
        "jordan_diagnostics": _diagnostic(estimate, covariance),
        "scientific_boundary": [
            "These are raw E_a/E_a_delta correlations, before pi_a and field normalization.",
            "J is invariant under separate nonzero rescalings of the two fields, but one small run is only a pipeline smoke test.",
            "Take L to infinity at fixed declared delta before comparing different delta rows.",
            "No value of kappa=C1*CL/C2 is assumed or estimated.",
        ],
    }


def write_archives(path: Path, archives: Sequence[BatchArchive]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = (
        ["batch", "samples"]
        + [f"sum_{name}" for name in PAIR_ORDER]
        + [f"sum_product_{name}" for name in PRODUCT_ORDER]
        + [f"sum_four_{name}" for name in PRODUCT_ORDER]
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for archive in archives:
            sums = archive.sufficient
            row = {"batch": archive.batch, "samples": sums.samples}
            row.update({f"sum_{name}": sums.pair_sums[name] for name in PAIR_ORDER})
            row.update(
                {
                    f"sum_product_{name}": sums.within_products[name]
                    for name in PRODUCT_ORDER
                }
            )
            row.update(
                {
                    f"sum_four_{name}": sums.four_spin_sums[name]
                    for name in PRODUCT_ORDER
                }
            )
            writer.writerow(row)


def read_archives(path: Path) -> list[BatchArchive]:
    archives = []
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            sums = SufficientSums(
                samples=int(raw["samples"]),
                pair_sums={name: int(raw[f"sum_{name}"]) for name in PAIR_ORDER},
                within_products={
                    name: int(raw[f"sum_product_{name}"]) for name in PRODUCT_ORDER
                },
                four_spin_sums={
                    name: int(raw[f"sum_four_{name}"]) for name in PRODUCT_ORDER
                },
            )
            archives.append(BatchArchive(batch=int(raw["batch"]), sufficient=sums))
    archives.sort(key=lambda row: row.batch)
    if [row.batch for row in archives] != list(range(len(archives))):
        raise ValueError("batch ids must be contiguous from zero")
    return archives


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run = subparsers.add_parser("run", help="sample, archive, and analyze")
    run.add_argument("--L", type=int, required=True)
    run.add_argument("--delta-denominator", type=int, choices=DELTA_DENOMINATORS, required=True)
    run.add_argument("--samples", type=int, required=True)
    run.add_argument("--batches", type=int, required=True)
    run.add_argument("--seed", type=int, required=True)
    run.add_argument("--workers", type=int, default=max(1, os.cpu_count() or 1))
    run.add_argument("--output-prefix", type=Path, required=True)
    analyze = subparsers.add_parser("analyze", help="rebuild JSON from archived sums")
    analyze.add_argument("--L", type=int, required=True)
    analyze.add_argument("--delta-denominator", type=int, choices=DELTA_DENOMINATORS, required=True)
    analyze.add_argument("--batches-csv", type=Path, required=True)
    analyze.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.command == "run":
        archives = run_batches(
            length=args.L,
            denominator=args.delta_denominator,
            samples=args.samples,
            batches=args.batches,
            seed=args.seed,
            workers=args.workers,
        )
        csv_path = Path(str(args.output_prefix) + ".batches.csv")
        json_path = Path(str(args.output_prefix) + ".json")
        write_archives(csv_path, archives)
        payload = analyze_archives(
            archives, length=args.L, denominator=args.delta_denominator
        )
        payload["run"] = {
            "seed": args.seed,
            "workers": args.workers,
            "python": sys.version,
            "platform": platform.platform(),
            "argv": sys.argv,
            "batch_archive": str(csv_path),
        }
        _write_json(json_path, payload)
        print(json_path)
    else:
        archives = read_archives(args.batches_csv)
        payload = analyze_archives(
            archives, length=args.L, denominator=args.delta_denominator
        )
        payload["reanalyzed_from"] = str(args.batches_csv)
        _write_json(args.output_json, payload)
        print(args.output_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
