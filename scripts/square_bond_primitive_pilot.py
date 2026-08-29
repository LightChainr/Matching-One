#!/usr/bin/env python3
"""Exact and fixed-p square-bond primitive-sector controls for issue #156.

The experiment is intentionally narrow: one N=4 exhaustive oracle and two
200k-replica Pell pilots at p=1/2.  It reports the three shortest primitive
rank-1 sectors, their Pinson--Arguin continuum residuals, and the real C3
contrasts C, Q, and S.  It does not fit an exponent.
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import json
from math import gcd, sqrt
import os
from pathlib import Path
import random
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import mpmath as mp

from integer_period_torus import (
    IntegerHomologyUnionFind,
    IntegerTorusGeometry,
    Matrix,
    integer_torus_geometry,
)
from pinson_arguin_primitive import (
    engine_to_paper,
    primitive_probability_direct,
)


MASK64 = (1 << 64) - 1
TARGET_LINES: Tuple[Tuple[str, Tuple[int, int]], ...] = (
    ("l0", (1, 0)),
    ("l1", (0, 1)),
    ("l2", (1, -1)),
)
CATEGORIES: Tuple[str, ...] = (
    "rank0",
    "l0",
    "l1",
    "l2",
    "rank1_other",
    "rank2",
    "invariant_failure",
)
PILOT_DESIGNS: Tuple[Tuple[str, Matrix], ...] = (
    ("pell_Dminus2_N30", ((6, 3), (0, 5))),
    ("pell_Dplus1_N56", ((8, 4), (0, 7))),
)
EXACT_DESIGN: Tuple[str, Matrix] = (
    "pell_fundamental_N4",
    ((2, 1), (0, 2)),
)


@dataclass(frozen=True)
class BatchResult:
    design: str
    batch: int
    samples: int
    counts: Dict[str, int]


def _canonical_line(vector: Tuple[int, int]) -> Tuple[int, int]:
    first, second = vector
    divisor = gcd(abs(first), abs(second))
    if divisor == 0:
        raise ValueError("zero has no primitive winding line")
    first //= divisor
    second //= divisor
    if first < 0 or (first == 0 and second < 0):
        first, second = -first, -second
    return first, second


def classify_bond_mask(
    geometry: IntegerTorusGeometry,
    mask: int,
) -> Tuple[str, Optional[Tuple[int, int]]]:
    """Classify the full open-bond homology image of one configuration."""

    union_find = IntegerHomologyUnionFind(geometry.n, geometry.periods)
    remaining = mask
    while remaining:
        lowest = remaining & -remaining
        edge_index = lowest.bit_length() - 1
        edge = geometry.primal_edges[edge_index]
        union_find.add_edge(edge.i, edge.j, edge.dx, edge.dy)
        remaining ^= lowest

    roots = {union_find.find(vertex)[0] for vertex in range(geometry.n)}
    directions = set()
    for root in roots:
        component = union_find.component(root)
        if component.rank == 2:
            return "rank2", None
        if component.rank == 1:
            directions.add(_canonical_line(component.basis[0]))

    if not directions:
        return "rank0", None
    if len(directions) != 1:
        return "invariant_failure", None
    direction = next(iter(directions))
    for name, line in TARGET_LINES:
        if direction == line:
            return name, direction
    return "rank1_other", direction


def exact_oracle() -> dict[str, object]:
    identifier, matrix = EXACT_DESIGN
    geometry = integer_torus_geometry(matrix)
    edge_count = len(geometry.primal_edges)
    if edge_count > 24:
        raise ValueError("exact oracle exceeds the 24-bond safety limit")
    counts = {category: 0 for category in CATEGORIES}
    other_lines: Dict[Tuple[int, int], int] = {}
    for mask in range(1 << edge_count):
        category, direction = classify_bond_mask(geometry, mask)
        counts[category] += 1
        if category == "rank1_other":
            assert direction is not None
            other_lines[direction] = other_lines.get(direction, 0) + 1
    return {
        "design": identifier,
        "period_matrix_rows": [list(row) for row in matrix],
        "N_vertices": geometry.n,
        "N_bonds": edge_count,
        "p": "1/2",
        "configurations": 1 << edge_count,
        "counts": counts,
        "rank1_other_lines": [
            {"engine_winding": list(line), "count": count}
            for line, count in sorted(other_lines.items())
        ],
        "passed": counts
        == {
            "rank0": 75,
            "l0": 57,
            "l1": 24,
            "l2": 24,
            "rank1_other": 1,
            "rank2": 75,
            "invariant_failure": 0,
        },
    }


def _splitmix64(value: int) -> int:
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return value ^ (value >> 31)


def _run_batch(task: Tuple[str, Matrix, int, int, int]) -> BatchResult:
    identifier, matrix, batch, samples, seed = task
    geometry = integer_torus_geometry(matrix)
    rng = random.Random(_splitmix64(seed ^ _splitmix64(batch + 1)))
    edge_count = len(geometry.primal_edges)
    counts = {category: 0 for category in CATEGORIES}
    for _ in range(samples):
        category, _ = classify_bond_mask(geometry, rng.getrandbits(edge_count))
        counts[category] += 1
    return BatchResult(identifier, batch, samples, counts)


def run_pilot_batches(
    *,
    samples: int,
    batches: int,
    seed: int,
    workers: int,
) -> List[BatchResult]:
    if samples <= 0 or batches <= 1 or samples % batches:
        raise ValueError("samples must be positive and divisible by batches>1")
    per_batch = samples // batches
    tasks = []
    for design_index, (identifier, matrix) in enumerate(PILOT_DESIGNS):
        design_seed = _splitmix64(seed ^ _splitmix64(design_index + 101))
        tasks.extend(
            (identifier, matrix, batch, per_batch, design_seed)
            for batch in range(batches)
        )
    if workers == 1:
        return [_run_batch(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(_run_batch, tasks))


def _tau_from_matrix(matrix: Matrix, *, dps: int = 80) -> mp.mpc:
    (a, b), (c, d) = matrix
    with mp.workdps(dps + 20):
        first = mp.mpc(a, c)
        second = mp.mpc(b, d)
        tau = second / first
        if mp.im(tau) <= 0:
            raise ValueError(
                "period matrix must define a positive-orientation modulus"
            )
        return +tau


def _covariance_of_mean(rows: Sequence[Sequence[float]]) -> List[List[float]]:
    count = len(rows)
    if count <= 1:
        raise ValueError("at least two batches are required for covariance")
    width = len(rows[0])
    means = [sum(row[index] for row in rows) / count for index in range(width)]
    return [
        [
            sum(
                (row[first] - means[first]) * (row[second] - means[second])
                for row in rows
            )
            / (count * (count - 1))
            for second in range(width)
        ]
        for first in range(width)
    ]


def _matrix_covariance(
    transform: Sequence[Sequence[float]],
    covariance: Sequence[Sequence[float]],
) -> List[List[float]]:
    return [
        [
            sum(
                transform[first][i] * covariance[i][j] * transform[second][j]
                for i in range(3)
                for j in range(3)
            )
            for second in range(3)
        ]
        for first in range(3)
    ]


def analyze_design(
    identifier: str,
    matrix: Matrix,
    rows: Sequence[BatchResult],
    *,
    dps: int,
) -> dict[str, object]:
    geometry = integer_torus_geometry(matrix)
    selected = sorted(
        (row for row in rows if row.design == identifier), key=lambda row: row.batch
    )
    if not selected:
        raise ValueError(f"no batches for {identifier}")
    total_samples = sum(row.samples for row in selected)
    totals = {
        category: sum(row.counts[category] for row in selected)
        for category in CATEGORIES
    }
    batch_probabilities = [
        [row.counts[name] / row.samples for name in ("l0", "l1", "l2")]
        for row in selected
    ]
    probabilities = [totals[name] / total_samples for name in ("l0", "l1", "l2")]
    covariance = _covariance_of_mean(batch_probabilities)

    tau = _tau_from_matrix(matrix, dps=dps)
    high_precision_baselines = [
        primitive_probability_direct(
            *engine_to_paper(line), tau, dps=dps
        )
        for _, line in TARGET_LINES
    ]
    baselines = [float(value) for value in high_precision_baselines]
    residuals = [value - baseline for value, baseline in zip(probabilities, baselines)]
    root_three_over_two = sqrt(3) / 2
    transform = (
        (1.0, -0.5, -0.5),
        (0.0, -root_three_over_two, root_three_over_two),
        (1.0, 1.0, 1.0),
    )
    names = ("C_nontrivial_real", "Q_reflection_null", "S_scalar")
    contrasts = [
        sum(transform[row][column] * residuals[column] for column in range(3))
        for row in range(3)
    ]
    contrast_covariance = _matrix_covariance(transform, covariance)
    contrast_payload = {}
    for index, name in enumerate(names):
        standard_error = sqrt(max(0.0, contrast_covariance[index][index]))
        contrast_payload[name] = {
            "value": contrasts[index],
            "standard_error": standard_error,
            "z": contrasts[index] / standard_error if standard_error else None,
        }
    return {
        "design": identifier,
        "period_matrix_rows": [list(row) for row in matrix],
        "N_vertices": geometry.n,
        "N_bonds": len(geometry.primal_edges),
        "tau_real": float(mp.re(tau)),
        "tau_imag": float(mp.im(tau)),
        "p": 0.5,
        "samples": total_samples,
        "batches": len(selected),
        "category_counts": totals,
        "target_sector_order": [name for name, _ in TARGET_LINES],
        "engine_windings": [list(line) for _, line in TARGET_LINES],
        "paper_types": [list(engine_to_paper(line)) for _, line in TARGET_LINES],
        "probabilities": probabilities,
        "continuum_baselines": baselines,
        "continuum_baselines_50dps": [
            mp.nstr(value, 50) for value in high_precision_baselines
        ],
        "residuals": residuals,
        "probability_covariance_of_mean": covariance,
        "contrasts": contrast_payload,
        "contrast_order": list(names),
        "contrast_covariance_of_mean": contrast_covariance,
    }


def build_result(
    batch_rows: Sequence[BatchResult],
    *,
    samples: int,
    batches: int,
    seed: int,
    dps: int,
) -> dict[str, object]:
    return {
        "schema": "p156-square-bond-primitive-pilot-v1",
        "issue": 156,
        "model": "square_bond_percolation",
        "p": "1/2",
        "seed": seed,
        "samples_per_design": samples,
        "batches_per_design": batches,
        "exact_oracle": exact_oracle(),
        "convention": {
            "engine_winding": "(u,v)=u*omega1+v*omega2",
            "paper_type": "{a,b}=a*omega1-b*omega2",
            "map": "(u,v)->{u,-v}, primitive saturated and unoriented",
            "target_engine_lines": [list(line) for _, line in TARGET_LINES],
        },
        "pilot": [
            analyze_design(identifier, matrix, batch_rows, dps=dps)
            for identifier, matrix in PILOT_DESIGNS
        ],
        "interpretation_boundary": [
            "C, Q, and S are continuum-baseline-subtracted real contrasts.",
            "Q is a reflection-null convention/control statistic.",
            "The two nontrivial complex C3 DFT modes are conjugate and counted once through C and Q.",
            "Two pilot sizes do not support an exponent fit or an H4-versus-H8 claim.",
        ],
    }


def write_batches(path: Path, rows: Sequence[BatchResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("design", "batch", "samples") + CATEGORIES)
        for row in sorted(rows, key=lambda item: (item.design, item.batch)):
            writer.writerow(
                (row.design, row.batch, row.samples)
                + tuple(row.counts[category] for category in CATEGORIES)
            )


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=200_000)
    parser.add_argument("--batches", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260829)
    parser.add_argument("--workers", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--output-prefix", type=Path)
    args = parser.parse_args(list(argv) if argv is not None else None)
    exact = exact_oracle()
    if not exact["passed"]:
        raise SystemExit("N=4 exact oracle failed")
    rows = run_pilot_batches(
        samples=args.samples,
        batches=args.batches,
        seed=args.seed,
        workers=args.workers,
    )
    payload = build_result(
        rows,
        samples=args.samples,
        batches=args.batches,
        seed=args.seed,
        dps=args.dps,
    )
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output_prefix is None:
        print(text, end="")
    else:
        args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
        Path(str(args.output_prefix) + ".json").write_text(text, encoding="utf-8")
        write_batches(Path(str(args.output_prefix) + ".batches.csv"), rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
