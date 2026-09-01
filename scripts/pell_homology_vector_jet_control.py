#!/usr/bin/env python3
"""Small square-bond full-homology Taylor-jet control for issues #156/#159.

The two period matrices are the first useful Pell approximants on opposite
sides of the Eisenstein modulus.  At each matrix, common random numbers are
used at ``p0-h``, ``p0`` and ``p0+h``.  Every observed primitive winding line
is retained separately, together with rank zero and rank two, so the output is
a probability-vector value, first derivative and second derivative rather
than a scalar H4 vote.
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from math import sqrt
import os
from pathlib import Path
import platform
import random
import time
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import mpmath as mp

from integer_period_torus import Matrix, integer_torus_geometry
from pinson_arguin_primitive import engine_to_paper, primitive_probability_direct
from square_bond_primitive_pilot import (
    PILOT_DESIGNS,
    TARGET_LINES,
    classify_bond_mask,
)


MASK64 = (1 << 64) - 1
POINTS = ("minus", "center", "plus")
JET_ORDERS = ("value", "d_dp", "d2_dp2")


@dataclass(frozen=True)
class BatchResult:
    design: str
    batch: int
    samples: int
    counts: Dict[str, Dict[str, int]]


def _splitmix64(value: int) -> int:
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return value ^ (value >> 31)


def category_key(category: str, direction: Optional[Tuple[int, int]]) -> str:
    if category in {"l0", "l1", "l2", "rank1_other"}:
        if direction is None:
            raise ValueError("rank-one category must retain its winding line")
        return f"rank1:{direction[0]},{direction[1]}"
    return category


def category_sort_key(category: str) -> Tuple[int, int, int]:
    if category == "rank0":
        return (0, 0, 0)
    if category.startswith("rank1:"):
        first, second = (int(value) for value in category.split(":", 1)[1].split(","))
        return (1, first, second)
    if category == "rank2":
        return (2, 0, 0)
    if category == "invariant_failure":
        return (3, 0, 0)
    raise ValueError(f"unknown category {category!r}")


def _masks_from_uniforms(
    rng: random.Random,
    edge_count: int,
    thresholds: Sequence[int],
) -> Tuple[int, ...]:
    masks = [0] * len(thresholds)
    for edge in range(edge_count):
        value = rng.getrandbits(64)
        bit = 1 << edge
        for index, threshold in enumerate(thresholds):
            if value < threshold:
                masks[index] |= bit
    return tuple(masks)


def _run_batch(
    task: Tuple[str, Matrix, int, int, int, float, float],
) -> BatchResult:
    identifier, matrix, batch, samples, seed, p0, h = task
    geometry = integer_torus_geometry(matrix)
    rng = random.Random(_splitmix64(seed ^ _splitmix64(batch + 1)))
    thresholds = tuple(int((p0 + offset * h) * (1 << 64)) for offset in (-1, 0, 1))
    counts: Dict[str, Dict[str, int]] = {point: {} for point in POINTS}
    for _ in range(samples):
        masks = _masks_from_uniforms(rng, len(geometry.primal_edges), thresholds)
        for point, mask in zip(POINTS, masks):
            category, direction = classify_bond_mask(geometry, mask)
            key = category_key(category, direction)
            counts[point][key] = counts[point].get(key, 0) + 1
    return BatchResult(identifier, batch, samples, counts)


def run_batches(
    *,
    samples: int,
    batches: int,
    seed: int,
    workers: int,
    p0: float,
    h: float,
) -> List[BatchResult]:
    if samples <= 0 or batches <= 1 or samples % batches:
        raise ValueError("samples must be positive and divisible by batches>1")
    if not (0.0 < p0 - h < p0 < p0 + h < 1.0):
        raise ValueError("p0-h, p0 and p0+h must lie strictly in (0,1)")
    per_batch = samples // batches
    tasks = []
    for design_index, (identifier, matrix) in enumerate(PILOT_DESIGNS):
        design_seed = _splitmix64(seed ^ _splitmix64(design_index + 401))
        tasks.extend(
            (identifier, matrix, batch, per_batch, design_seed, p0, h)
            for batch in range(batches)
        )
    if workers == 1:
        return [_run_batch(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return list(executor.map(_run_batch, tasks))


def _tau_from_matrix(matrix: Matrix, *, dps: int) -> mp.mpc:
    (a, b), (c, d) = matrix
    with mp.workdps(dps + 20):
        tau = mp.mpc(b, d) / mp.mpc(a, c)
        if mp.im(tau) <= 0:
            raise ValueError("period matrix must define positive orientation")
        return +tau


def _covariance_of_mean(rows: Sequence[Sequence[float]]) -> List[List[float]]:
    count = len(rows)
    if count <= 1:
        raise ValueError("at least two batches are required")
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


def _jet_from_probabilities(
    minus: Sequence[float], center: Sequence[float], plus: Sequence[float], h: float
) -> Tuple[List[float], List[float], List[float]]:
    return (
        list(center),
        [(right - left) / (2.0 * h) for left, right in zip(minus, plus)],
        [
            (right - 2.0 * middle + left) / (h * h)
            for left, middle, right in zip(minus, center, plus)
        ],
    )


def _batch_jet(
    row: BatchResult, support: Sequence[str], h: float
) -> Tuple[List[float], List[float], List[float]]:
    points = [
        [row.counts[point].get(category, 0) / row.samples for category in support]
        for point in POINTS
    ]
    return _jet_from_probabilities(points[0], points[1], points[2], h)


def _line_from_key(category: str) -> Optional[Tuple[int, int]]:
    if not category.startswith("rank1:"):
        return None
    return tuple(int(value) for value in category.split(":", 1)[1].split(","))  # type: ignore[return-value]


def _three_line_contrasts(
    batch_jets: Sequence[Tuple[List[float], List[float], List[float]]],
    support: Sequence[str],
    baselines: Mapping[str, float],
    samples: int,
) -> Dict[str, object]:
    keys = [category_key(name, line) for name, line in TARGET_LINES]
    positions = [support.index(key) for key in keys]
    root_three_over_two = sqrt(3.0) / 2.0
    transform = (
        (1.0, -0.5, -0.5),
        (0.0, -root_three_over_two, root_three_over_two),
        (1.0, 1.0, 1.0),
    )
    names = ("C_nontrivial_real", "Q_reflection_null", "S_scalar")
    result: Dict[str, object] = {}
    for order_index, order in enumerate(JET_ORDERS):
        rows = []
        for jet in batch_jets:
            values = [jet[order_index][position] for position in positions]
            if order == "value":
                values = [value - baselines[key] for value, key in zip(values, keys)]
            rows.append(
                [
                    sum(transform[row][column] * values[column] for column in range(3))
                    for row in range(3)
                ]
            )
        means = [sum(row[index] for row in rows) / len(rows) for index in range(3)]
        covariance = _covariance_of_mean(rows)
        payload = {}
        for index, name in enumerate(names):
            standard_error = sqrt(max(0.0, covariance[index][index]))
            z = means[index] / standard_error if standard_error else None
            required = None
            if z not in (None, 0.0):
                required = int(round(samples * (5.0 / abs(z)) ** 2))
            payload[name] = {
                "estimate": means[index],
                "standard_error": standard_error,
                "z": z,
                "samples_for_abs_z5_at_same_effect": required,
            }
        result[order] = {"coordinates": payload, "covariance_of_mean": covariance}
    return result


def analyze_design(
    identifier: str,
    matrix: Matrix,
    rows: Sequence[BatchResult],
    *,
    h: float,
    dps: int,
) -> Dict[str, object]:
    selected = sorted(
        (row for row in rows if row.design == identifier), key=lambda row: row.batch
    )
    if not selected:
        raise ValueError(f"no batches for {identifier}")
    samples = sum(row.samples for row in selected)
    support = sorted(
        {
            category
            for row in selected
            for point in POINTS
            for category in row.counts[point]
        }
        | {"rank0", "rank2", "invariant_failure"},
        key=category_sort_key,
    )
    totals = {
        point: {
            category: sum(row.counts[point].get(category, 0) for row in selected)
            for category in support
        }
        for point in POINTS
    }
    probabilities = {
        point: [totals[point][category] / samples for category in support]
        for point in POINTS
    }
    jet = _jet_from_probabilities(
        probabilities["minus"], probabilities["center"], probabilities["plus"], h
    )
    batch_jets = [_batch_jet(row, support, h) for row in selected]
    covariances = [
        _covariance_of_mean([row[order] for row in batch_jets]) for order in range(3)
    ]
    tau = _tau_from_matrix(matrix, dps=dps)
    baselines: Dict[str, float] = {}
    baseline_text: Dict[str, str] = {}
    for category in support:
        line = _line_from_key(category)
        if line is None:
            continue
        value = primitive_probability_direct(*engine_to_paper(line), tau, dps=dps)
        baselines[category] = float(value)
        baseline_text[category] = mp.nstr(value, 40)

    coordinates = []
    for index, category in enumerate(support):
        estimates = {}
        for order_index, order in enumerate(JET_ORDERS):
            standard_error = sqrt(max(0.0, covariances[order_index][index][index]))
            estimate = jet[order_index][index]
            estimates[order] = {
                "estimate": estimate,
                "standard_error": standard_error,
                "z": estimate / standard_error if standard_error else None,
            }
        coordinates.append(
            {
                "category": category,
                "continuum_baseline_at_p0": baseline_text.get(category),
                "center_residual": (
                    jet[0][index] - baselines[category]
                    if category in baselines
                    else None
                ),
                "jet": estimates,
            }
        )
    geometry = integer_torus_geometry(matrix)
    return {
        "design": identifier,
        "period_matrix_rows": [list(row) for row in matrix],
        "N_vertices": geometry.n,
        "N_bonds": len(geometry.primal_edges),
        "tau_real": float(mp.re(tau)),
        "tau_imag": float(mp.im(tau)),
        "samples": samples,
        "batches": len(selected),
        "support": support,
        "probabilities": probabilities,
        "coordinates": coordinates,
        "jet_covariance_of_mean": {
            order: covariance for order, covariance in zip(JET_ORDERS, covariances)
        },
        "three_shortest_line_contrasts": _three_line_contrasts(
            batch_jets, support, baselines, samples
        ),
        "count_conservation": {
            point: sum(totals[point].values()) == samples for point in POINTS
        },
        "invariant_failures": {
            point: totals[point].get("invariant_failure", 0) for point in POINTS
        },
    }


def exact_n4_jet(p0: float = 0.5) -> Dict[str, object]:
    identifier = "pell_fundamental_N4"
    matrix: Matrix = ((2, 1), (0, 2))
    geometry = integer_torus_geometry(matrix)
    edge_count = len(geometry.primal_edges)
    polynomials: Dict[str, List[int]] = {}
    for mask in range(1 << edge_count):
        category, direction = classify_bond_mask(geometry, mask)
        key = category_key(category, direction)
        coefficients = polynomials.setdefault(key, [0] * (edge_count + 1))
        coefficients[mask.bit_count()] += 1
    support = sorted(polynomials, key=category_sort_key)
    q0 = 1.0 - p0
    coordinates = []
    for category in support:
        value = first = second = 0.0
        for occupied, multiplicity in enumerate(polynomials[category]):
            if not multiplicity:
                continue
            weight = p0**occupied * q0 ** (edge_count - occupied)
            log_first = occupied / p0 - (edge_count - occupied) / q0
            log_second = (
                log_first * log_first
                - occupied / (p0 * p0)
                - (edge_count - occupied) / (q0 * q0)
            )
            value += multiplicity * weight
            first += multiplicity * weight * log_first
            second += multiplicity * weight * log_second
        coordinates.append(
            {
                "category": category,
                "edge_count_polynomial": polynomials[category],
                "jet": {"value": value, "d_dp": first, "d2_dp2": second},
            }
        )
    return {
        "design": identifier,
        "period_matrix_rows": [list(row) for row in matrix],
        "N_vertices": geometry.n,
        "N_bonds": edge_count,
        "p0": p0,
        "configurations": 1 << edge_count,
        "coordinates": coordinates,
        "sum_rules": {
            order: sum(row["jet"][order] for row in coordinates)
            for order in JET_ORDERS
        },
    }


def write_batches(path: Path, rows: Sequence[BatchResult]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("design", "batch", "samples", "point", "category", "count"))
        for row in sorted(rows, key=lambda item: (item.design, item.batch)):
            for point in POINTS:
                for category in sorted(row.counts[point], key=category_sort_key):
                    writer.writerow(
                        (row.design, row.batch, row.samples, point, category, row.counts[point][category])
                    )


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=240_000)
    parser.add_argument("--batches", type=int, default=48)
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--workers", type=int, default=min(14, max(1, os.cpu_count() or 1)))
    parser.add_argument("--p0", type=float, default=0.5)
    parser.add_argument("--h", type=float, default=0.01)
    parser.add_argument("--dps", type=int, default=70)
    parser.add_argument("--source-commit", default="unknown")
    parser.add_argument("--output-prefix", type=Path)
    args = parser.parse_args(list(argv) if argv is not None else None)
    started = time.perf_counter()
    rows = run_batches(
        samples=args.samples,
        batches=args.batches,
        seed=args.seed,
        workers=args.workers,
        p0=args.p0,
        h=args.h,
    )
    payload = {
        "schema": "matching-one/pell-homology-vector-jet-control/v1",
        "issues": [156, 159],
        "status": "small_square_bond_control_complete",
        "model": "critical_square_bond_percolation",
        "p_grid": [args.p0 - args.h, args.p0, args.p0 + args.h],
        "finite_difference": {
            "d_dp": "(P(p0+h)-P(p0-h))/(2h)",
            "d2_dp2": "(P(p0+h)-2P(p0)+P(p0-h))/h^2",
            "common_random_numbers": True,
        },
        "seed": args.seed,
        "samples_per_design": args.samples,
        "batches_per_design": args.batches,
        "exact_n4_oracle": exact_n4_jet(args.p0),
        "pilot": [
            analyze_design(identifier, matrix, rows, h=args.h, dps=args.dps)
            for identifier, matrix in PILOT_DESIGNS
        ],
        "runtime": {
            "source_commit": args.source_commit,
            "utc_finished": datetime.now(timezone.utc).isoformat(),
            "wall_seconds": time.perf_counter() - started,
            "workers": args.workers,
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "interpretation_boundary": [
            "This is the requested small square-bond control, not N418/N780 production.",
            "The vector retains every observed primitive winding line plus rank zero and rank two.",
            "The first and second p-jets use a fixed three-point finite-difference contract with common random numbers.",
            "The control tests vector organization at both first useful Pell moduli; it does not identify H4 versus H8 by itself.",
        ],
    }
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
