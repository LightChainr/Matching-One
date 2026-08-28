#!/usr/bin/env python3
"""Bidirectional Newman--Ziff threshold-rank reference implementation.

For a permutation ``v_1,...,v_N``, define ``B_k={v_1,...,v_k}`` and its
white complement ``W_k``.  This module freezes the off-by-one convention

* ``K_plus  = min{k in 1..N: B_k has primal cross wrapping}``;
* ``K_minus = min{k in 1..N: W_k has lost matching cross wrapping}``.

The reverse matching sweep first crosses after ``r`` white additions, hence
``K_minus = N-r+1``.  For the matching pair, ``K_minus <= K_plus`` is checked
for every sample.  Integer marginal and sparse joint counts are sufficient to
reconstruct the canonical cross-channel matching curve without rerunning the
connectivity calculation.

This Python code prioritizes auditable conventions and tiny exact tests.  It
is not the million-sample production kernel.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from integer_period_torus import (
    IntegerHomologyUnionFind,
    IntegerTorusGeometry,
    axis_integer_torus,
    diamond_integer_torus,
    gaussian_integer_torus,
    integer_torus_geometry,
)
from matched_torus_reference import Edge


MASK64 = (1 << 64) - 1


class SplitMix64:
    """Small reproducible 64-bit generator used for counter-derived streams."""

    def __init__(self, state: int) -> None:
        self.state = state & MASK64

    def next_u64(self) -> int:
        self.state = (self.state + 0x9E3779B97F4A7C15) & MASK64
        value = self.state
        value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
        value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
        return (value ^ (value >> 31)) & MASK64

    def randbelow(self, bound: int) -> int:
        if bound <= 0:
            raise ValueError("bound must be positive")
        limit = (1 << 64) - ((1 << 64) % bound)
        while True:
            value = self.next_u64()
            if value < limit:
                return value % bound


def _mix64(value: int) -> int:
    generator = SplitMix64(value)
    return generator.next_u64()


def counter_permutation(n: int, seed: int, counter: int) -> Tuple[int, ...]:
    """Return a deterministic Fisher--Yates permutation for one sample counter."""

    if n < 0 or counter < 0:
        raise ValueError("n and counter must be nonnegative")
    stream_key = _mix64((seed & MASK64) ^ _mix64(counter + 0xD1B54A32D192ED03))
    generator = SplitMix64(stream_key)
    values = list(range(n))
    for stop in range(n - 1, 0, -1):
        other = generator.randbelow(stop + 1)
        values[stop], values[other] = values[other], values[stop]
    return tuple(values)


def _incident_edges(n: int, edges: Iterable[Edge]) -> Tuple[Tuple[Edge, ...], ...]:
    incident: List[List[Edge]] = [[] for _ in range(n)]
    for edge in edges:
        incident[edge.i].append(edge)
        if edge.j != edge.i:
            incident[edge.j].append(edge)
    return tuple(tuple(values) for values in incident)


def _first_cross_rank_validated(
    geometry: IntegerTorusGeometry,
    permutation: Sequence[int],
    incident: Sequence[Sequence[Edge]],
) -> int:
    active = [False] * geometry.n
    union_find = IntegerHomologyUnionFind(geometry.n, geometry.periods)

    for rank, vertex in enumerate(permutation, start=1):
        active[vertex] = True
        for edge in incident[vertex]:
            if active[edge.i] and active[edge.j]:
                union_find.add_edge(edge.i, edge.j, edge.dx, edge.dy)
        if union_find.component(vertex).cross:
            return rank
    raise AssertionError("fully occupied periodic graph did not cross wrap")


class ThresholdRankEngine:
    """Cache geometry adjacency while evaluating many independent permutations."""

    def __init__(self, geometry: IntegerTorusGeometry) -> None:
        self.geometry = geometry
        self.primal_incident = _incident_edges(geometry.n, geometry.primal_edges)
        self.matching_incident = _incident_edges(geometry.n, geometry.matching_edges)

    def first_cross_rank(
        self, permutation: Sequence[int], *, matching: bool = False
    ) -> int:
        if (
            len(permutation) != self.geometry.n
            or set(permutation) != set(range(self.geometry.n))
        ):
            raise ValueError("permutation must contain every vertex exactly once")
        incident = self.matching_incident if matching else self.primal_incident
        return _first_cross_rank_validated(self.geometry, permutation, incident)

    def threshold_ranks(self, permutation: Sequence[int]) -> Tuple[int, int]:
        if (
            len(permutation) != self.geometry.n
            or set(permutation) != set(range(self.geometry.n))
        ):
            raise ValueError("permutation must contain every vertex exactly once")
        k_plus = _first_cross_rank_validated(
            self.geometry, permutation, self.primal_incident
        )
        reverse_white_rank = _first_cross_rank_validated(
            self.geometry, tuple(reversed(permutation)), self.matching_incident
        )
        k_minus = self.geometry.n - reverse_white_rank + 1
        if k_minus > k_plus:
            raise AssertionError(
                "matching cross thresholds are out of order: "
                f"K_minus={k_minus}, K_plus={k_plus}"
            )
        return k_minus, k_plus


def first_cross_rank(
    geometry: IntegerTorusGeometry,
    permutation: Sequence[int],
    *,
    matching: bool = False,
) -> int:
    """Return the 1-based occupation rank at which cross wrapping first occurs."""

    return ThresholdRankEngine(geometry).first_cross_rank(
        permutation, matching=matching
    )


def threshold_ranks(
    geometry: IntegerTorusGeometry, permutation: Sequence[int]
) -> Tuple[int, int]:
    """Return ``(K_minus,K_plus)`` under the module's frozen convention."""

    return ThresholdRankEngine(geometry).threshold_ranks(permutation)


@dataclass
class ThresholdRankCounts:
    n: int
    kminus: List[int] = field(init=False)
    kplus: List[int] = field(init=False)
    joint: Dict[Tuple[int, int], int] = field(default_factory=dict)
    sample_count: int = 0
    sum_kminus: int = 0
    sum_kplus: int = 0
    sum_kminus2: int = 0
    sum_kplus2: int = 0
    sum_product: int = 0
    sum_gap: int = 0
    sum_gap2: int = 0

    def __post_init__(self) -> None:
        if self.n <= 0:
            raise ValueError("n must be positive")
        self.kminus = [0] * (self.n + 1)
        self.kplus = [0] * (self.n + 1)

    def add(self, k_minus: int, k_plus: int) -> None:
        if not 1 <= k_minus <= k_plus <= self.n:
            raise ValueError("require 1 <= K_minus <= K_plus <= N")
        self.kminus[k_minus] += 1
        self.kplus[k_plus] += 1
        pair = (k_minus, k_plus)
        self.joint[pair] = self.joint.get(pair, 0) + 1
        self.sample_count += 1
        self.sum_kminus += k_minus
        self.sum_kplus += k_plus
        self.sum_kminus2 += k_minus * k_minus
        self.sum_kplus2 += k_plus * k_plus
        self.sum_product += k_minus * k_plus
        gap = k_plus - k_minus
        self.sum_gap += gap
        self.sum_gap2 += gap * gap

    def moments(self) -> Dict[str, int]:
        return {
            "sum_K_minus": self.sum_kminus,
            "sum_K_plus": self.sum_kplus,
            "sum_K_minus_squared": self.sum_kminus2,
            "sum_K_plus_squared": self.sum_kplus2,
            "sum_K_minus_times_K_plus": self.sum_product,
            "sum_gap": self.sum_gap,
            "sum_gap_squared": self.sum_gap2,
        }


def simulate(
    geometry: IntegerTorusGeometry,
    samples: int,
    seed: int,
    counter_start: int = 0,
) -> ThresholdRankCounts:
    if samples <= 0 or counter_start < 0:
        raise ValueError("samples must be positive and counter_start nonnegative")
    counts = ThresholdRankCounts(geometry.n)
    engine = ThresholdRankEngine(geometry)
    for counter in range(counter_start, counter_start + samples):
        permutation = counter_permutation(geometry.n, seed, counter)
        counts.add(*engine.threshold_ranks(permutation))
    return counts


def enumerate_exact(geometry: IntegerTorusGeometry) -> ThresholdRankCounts:
    """Enumerate all N! permutations for tiny regression systems."""

    if geometry.n > 9:
        raise ValueError("exact permutation enumeration is limited to N<=9")
    counts = ThresholdRankCounts(geometry.n)
    engine = ThresholdRankEngine(geometry)
    for permutation in itertools.permutations(range(geometry.n)):
        counts.add(*engine.threshold_ranks(permutation))
    return counts


def _write_histogram(path: Path, name: str, values: Sequence[int]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow([name, "count"])
        for rank, count in enumerate(values):
            if count:
                writer.writerow([rank, count])


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_counts(
    output_dir: Path,
    geometry: IntegerTorusGeometry,
    counts: ThresholdRankCounts,
    *,
    seed: int,
    counter_start: int,
    elapsed_seconds: float,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    minus_path = output_dir / "kminus_hist.csv"
    plus_path = output_dir / "kplus_hist.csv"
    joint_path = output_dir / "joint_hist.csv"
    _write_histogram(minus_path, "K_minus", counts.kminus)
    _write_histogram(plus_path, "K_plus", counts.kplus)
    with joint_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["K_minus", "K_plus", "count"])
        for (k_minus, k_plus), count in sorted(counts.joint.items()):
            writer.writerow([k_minus, k_plus, count])

    metadata = {
        "format_version": 1,
        "geometry": geometry.name,
        "N": geometry.n,
        "period_matrix_row_major": [list(row) for row in geometry.periods.matrix],
        "period_matrix_convention": "columns are the two lifted period vectors",
        "primal_forward_vectors": sorted(
            {tuple((edge.dx, edge.dy)) for edge in geometry.primal_edges}
        ),
        "matching_forward_vectors": sorted(
            {tuple((edge.dx, edge.dy)) for edge in geometry.matching_edges}
        ),
        "channel": "cross_rank_2_single_component",
        "off_by_one_convention": {
            "K_plus": "first k in 1..N where B_k primal-cross-wraps",
            "K_minus": "first k in 1..N where W_k matching-cross no longer wraps",
            "reverse_conversion": "K_minus = N - reverse_white_first_cross_rank + 1",
        },
        "sample_count": counts.sample_count,
        "rng": {
            "algorithm": "counter-derived SplitMix64 streams + unbiased Fisher-Yates",
            "seed_u64": seed & MASK64,
            "sample_counter_start_inclusive": counter_start,
            "sample_counter_stop_exclusive": counter_start + counts.sample_count,
        },
        "first_second_joint_integer_moments": counts.moments(),
        "elapsed_seconds": elapsed_seconds,
    }
    metadata_path = output_dir / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    checksum_paths = (minus_path, plus_path, joint_path, metadata_path)
    (output_dir / "checksums.sha256").write_text(
        "".join(f"{_sha256(path)}  {path.name}\n" for path in checksum_paths),
        encoding="utf-8",
    )


def _geometry_from_args(args: argparse.Namespace) -> IntegerTorusGeometry:
    if args.geometry == "axis":
        if args.L is None:
            raise SystemExit("--geometry axis requires --L")
        return axis_integer_torus(args.L)
    if args.geometry == "diamond":
        if args.L is None:
            raise SystemExit("--geometry diamond requires --L")
        return diamond_integer_torus(args.L)
    if args.geometry == "gaussian":
        if args.a is None or args.b is None:
            raise SystemExit("--geometry gaussian requires --a and --b")
        return gaussian_integer_torus(args.a, args.b)
    if args.matrix is None:
        raise SystemExit("--geometry matrix requires --matrix A B C D")
    a, b, c, d = args.matrix
    return integer_torus_geometry(((a, b), (c, d)))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--geometry", choices=("axis", "diamond", "gaussian", "matrix"), required=True
    )
    parser.add_argument("--L", type=int)
    parser.add_argument("--a", type=int)
    parser.add_argument("--b", type=int)
    parser.add_argument("--matrix", type=int, nargs=4, metavar=("A", "B", "C", "D"))
    parser.add_argument("--samples", type=int, required=True)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--counter-start", type=int, default=0)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    geometry = _geometry_from_args(args)
    started = time.perf_counter()
    counts = simulate(geometry, args.samples, args.seed, args.counter_start)
    elapsed = time.perf_counter() - started
    write_counts(
        args.output_dir,
        geometry,
        counts,
        seed=args.seed,
        counter_start=args.counter_start,
        elapsed_seconds=elapsed,
    )
    print(
        f"geometry={geometry.name} N={geometry.n} samples={counts.sample_count} "
        f"joint_bins={len(counts.joint)} elapsed_seconds={elapsed:.6f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
