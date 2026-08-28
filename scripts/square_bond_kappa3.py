#!/usr/bin/env python3
"""Exact and CPU Monte Carlo kappa_3 control for square bond percolation.

The configuration observable is the difference between the event that open
primal bonds wrap a square torus and the event that closed dual bonds wrap the
dual torus.  Its expectation is the finite-volume matching function.  At the
exact threshold p=1/2, ``kappa3_half_score`` converts one Bernoulli ensemble
into estimates of the first and third derivatives without finite differences.

Monte Carlo work is divided into independently seeded blocks.  The block
seeds, rather than worker scheduling, determine the random stream, so changing
``--workers`` does not change the output.  The same blocks provide a
delete-one-block jackknife estimate of the ratio uncertainty.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import platform
import random
import sys
from typing import Iterable, Optional, Sequence, Union

from kappa3_half_score import first_score, third_weight
from torus_homology import HomologyUnionFind


MASK64 = (1 << 64) - 1


@dataclass(frozen=True)
class BondPair:
    """One primal bond and the dual bond crossing it."""

    primal: tuple[int, int, int, int]
    dual: tuple[int, int, int, int]


@dataclass(frozen=True)
class BlockSums:
    samples: int
    observable_sum: int
    observable_square_sum: int
    first_sum: int
    first_square_sum: int
    third_sum: int
    third_square_sum: int
    first_third_sum: int


def _vertex_id(length: int, x: int, y: int) -> int:
    return (y % length) * length + (x % length)


def square_bond_pairs(length: int) -> tuple[BondPair, ...]:
    """Return the 2 L^2 primal bonds paired with their crossing dual bonds."""

    if length < 2:
        raise ValueError("L must be at least 2")
    pairs: list[BondPair] = []
    for y in range(length):
        for x in range(length):
            # Horizontal primal bond; the crossing dual bond is vertical.
            pairs.append(
                BondPair(
                    primal=(
                        _vertex_id(length, x, y),
                        _vertex_id(length, x + 1, y),
                        1,
                        0,
                    ),
                    dual=(
                        _vertex_id(length, x, y - 1),
                        _vertex_id(length, x, y),
                        0,
                        1,
                    ),
                )
            )
            # Vertical primal bond; the crossing dual bond is horizontal.
            pairs.append(
                BondPair(
                    primal=(
                        _vertex_id(length, x, y),
                        _vertex_id(length, x, y + 1),
                        0,
                        1,
                    ),
                    dual=(
                        _vertex_id(length, x - 1, y),
                        _vertex_id(length, x, y),
                        1,
                        0,
                    ),
                )
            )
    return tuple(pairs)


def _has_wrapping(union_find: HomologyUnionFind) -> bool:
    return any(
        union_find.parent[root] == root and bool(union_find.basis[root])
        for root in range(len(union_find.parent))
    )


def wrapping_difference(length: int, mask: int, pairs: Sequence[BondPair]) -> int:
    """Return open-primal wrapping minus closed-dual wrapping."""

    vertex_count = length * length
    primal = HomologyUnionFind(vertex_count, (length, length))
    dual = HomologyUnionFind(vertex_count, (length, length))
    for bond_index, pair in enumerate(pairs):
        if (mask >> bond_index) & 1:
            primal.add_edge(*pair.primal)
        else:
            dual.add_edge(*pair.dual)
    return int(_has_wrapping(primal)) - int(_has_wrapping(dual))


def _fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def exact_estimate(length: int) -> dict[str, object]:
    """Exhaustively enumerate a tiny square torus at p=1/2."""

    pairs = square_bond_pairs(length)
    bond_count = len(pairs)
    if bond_count > 24:
        raise ValueError(
            f"exact enumeration for L={length} requires 2^{bond_count} states; "
            "the safety limit is 24 bonds"
        )

    configuration_count = 1 << bond_count
    d_sum = 0
    d2_sum = 0
    first_sum = 0
    first2_sum = 0
    third_sum = 0
    third2_sum = 0
    first_third_sum = 0
    for mask in range(configuration_count):
        occupied = bin(mask).count("1")
        observable = wrapping_difference(length, mask, pairs)
        first = first_score(bond_count, occupied) * observable
        third = third_weight(bond_count, occupied) * observable
        d_sum += observable
        d2_sum += observable * observable
        first_sum += first
        first2_sum += first * first
        third_sum += third
        third2_sum += third * third
        first_third_sum += first * third

    first_derivative = Fraction(first_sum, configuration_count)
    third_derivative = Fraction(third_sum, configuration_count)
    if first_derivative == 0:
        raise ZeroDivisionError("exact first derivative is zero")
    kappa3 = third_derivative / first_derivative**3

    def population_variance(total: int, square_total: int) -> Fraction:
        mean = Fraction(total, configuration_count)
        return Fraction(square_total, configuration_count) - mean**2

    covariance = (
        Fraction(first_third_sum, configuration_count)
        - first_derivative * third_derivative
    )
    return {
        "model": "square_bond_square_torus",
        "mode": "exact",
        "L": length,
        "vertices": length * length,
        "bernoulli_variables": bond_count,
        "configurations": configuration_count,
        "mean_observable": float(Fraction(d_sum, configuration_count)),
        "first_derivative": float(first_derivative),
        "third_derivative": float(third_derivative),
        "kappa3": float(kappa3),
        "exact": {
            "mean_observable": _fraction_text(Fraction(d_sum, configuration_count)),
            "first_derivative": _fraction_text(first_derivative),
            "third_derivative": _fraction_text(third_derivative),
            "kappa3": _fraction_text(kappa3),
        },
        "per_configuration_score_variance": {
            "first": float(population_variance(first_sum, first2_sum)),
            "third": float(population_variance(third_sum, third2_sum)),
            "first_third_covariance": float(covariance),
        },
        "jackknife": None,
    }


def _splitmix64(value: int) -> int:
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return (value ^ (value >> 31)) & MASK64


def _simulate_block(task: tuple[int, int, int]) -> BlockSums:
    length, samples, seed = task
    pairs = square_bond_pairs(length)
    bond_count = len(pairs)
    rng = random.Random(seed)
    d_sum = d2_sum = 0
    first_sum = first2_sum = 0
    third_sum = third2_sum = first_third_sum = 0
    for _ in range(samples):
        mask = rng.getrandbits(bond_count)
        occupied = bin(mask).count("1")
        observable = wrapping_difference(length, mask, pairs)
        first = first_score(bond_count, occupied) * observable
        third = third_weight(bond_count, occupied) * observable
        d_sum += observable
        d2_sum += observable * observable
        first_sum += first
        first2_sum += first * first
        third_sum += third
        third2_sum += third * third
        first_third_sum += first * third
    return BlockSums(
        samples=samples,
        observable_sum=d_sum,
        observable_square_sum=d2_sum,
        first_sum=first_sum,
        first_square_sum=first2_sum,
        third_sum=third_sum,
        third_square_sum=third2_sum,
        first_third_sum=first_third_sum,
    )


def _combine(blocks: Iterable[BlockSums]) -> BlockSums:
    values = tuple(blocks)
    return BlockSums(
        **{
            field: sum(getattr(value, field) for value in values)
            for field in BlockSums.__dataclass_fields__
        }
    )


def _sample_variance(total: int, square_total: int, samples: int) -> float:
    if samples < 2:
        return math.nan
    return (square_total - total * total / samples) / (samples - 1)


def _ratio(first_mean: float, third_mean: float) -> float:
    if first_mean == 0.0:
        raise ZeroDivisionError("estimated first derivative is zero")
    return third_mean / first_mean**3


def _jackknife(
    total: BlockSums, blocks: Sequence[BlockSums]
) -> dict[str, Union[float, int]]:
    if len(blocks) < 2:
        raise ValueError("jackknife requires at least two blocks")
    if len({block.samples for block in blocks}) != 1:
        raise ValueError("delete-one-block jackknife requires equal block sizes")

    full_first = total.first_sum / total.samples
    full_third = total.third_sum / total.samples
    full_ratio = _ratio(full_first, full_third)
    leave_one_out: list[float] = []
    for block in blocks:
        remaining = total.samples - block.samples
        first = (total.first_sum - block.first_sum) / remaining
        third = (total.third_sum - block.third_sum) / remaining
        leave_one_out.append(_ratio(first, third))

    count = len(leave_one_out)
    leave_mean = math.fsum(leave_one_out) / count
    variance = (count - 1) / count * math.fsum(
        (value - leave_mean) ** 2 for value in leave_one_out
    )
    standard_error = math.sqrt(variance)
    bias_corrected = count * full_ratio - (count - 1) * leave_mean
    return {
        "blocks": count,
        "block_samples": blocks[0].samples,
        "raw_estimate": full_ratio,
        "bias_corrected_estimate": bias_corrected,
        "estimated_bias": full_ratio - bias_corrected,
        "variance": variance,
        "standard_error": standard_error,
        "normal_95_low": bias_corrected - 1.96 * standard_error,
        "normal_95_high": bias_corrected + 1.96 * standard_error,
    }


def monte_carlo_estimate(
    length: int,
    samples: int,
    blocks: int,
    seed: int,
    workers: int,
) -> dict[str, object]:
    """Estimate square-bond kappa_3 with independent reproducible blocks."""

    if samples <= 0:
        raise ValueError("samples must be positive")
    if blocks < 2 or samples % blocks:
        raise ValueError("samples must be divisible by at least two blocks")
    if workers <= 0:
        raise ValueError("workers must be positive")
    square_bond_pairs(length)  # Validate before creating worker processes.

    per_block = samples // blocks
    tasks = [
        (length, per_block, _splitmix64(seed + block_index))
        for block_index in range(blocks)
    ]
    if workers == 1:
        block_results = [_simulate_block(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            block_results = list(executor.map(_simulate_block, tasks))

    total = _combine(block_results)
    first_mean = total.first_sum / samples
    third_mean = total.third_sum / samples
    kappa3 = _ratio(first_mean, third_mean)
    first_variance = _sample_variance(
        total.first_sum, total.first_square_sum, samples
    )
    third_variance = _sample_variance(
        total.third_sum, total.third_square_sum, samples
    )
    covariance = (
        total.first_third_sum - total.first_sum * total.third_sum / samples
    ) / (samples - 1)
    observable_variance = _sample_variance(
        total.observable_sum, total.observable_square_sum, samples
    )
    return {
        "model": "square_bond_square_torus",
        "mode": "monte_carlo",
        "L": length,
        "vertices": length * length,
        "bernoulli_variables": 2 * length * length,
        "samples": samples,
        "blocks": blocks,
        "seed": seed,
        "block_seeds": [task[2] for task in tasks],
        "mean_observable": total.observable_sum / samples,
        "mean_observable_standard_error": math.sqrt(observable_variance / samples),
        "first_derivative": first_mean,
        "first_derivative_variance": first_variance / samples,
        "first_derivative_standard_error": math.sqrt(first_variance / samples),
        "third_derivative": third_mean,
        "third_derivative_variance": third_variance / samples,
        "third_derivative_standard_error": math.sqrt(third_variance / samples),
        "first_third_score_covariance": covariance,
        "kappa3": kappa3,
        "jackknife": _jackknife(total, block_results),
        "block_sums": [asdict(block) for block in block_results],
    }


def _parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("exact", "monte-carlo"), required=True)
    parser.add_argument("--sizes", nargs="+", type=int, required=True, metavar="L")
    parser.add_argument("--samples", type=int, default=10000)
    parser.add_argument("--blocks", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260828)
    parser.add_argument("--workers", type=int, default=max(1, os.cpu_count() or 1))
    parser.add_argument("--output", type=Path)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parse_args(argv)
    if len(set(args.sizes)) != len(args.sizes):
        raise ValueError("sizes must not contain duplicates")

    results: list[dict[str, object]] = []
    for length in args.sizes:
        if args.mode == "exact":
            results.append(exact_estimate(length))
        else:
            geometry_seed = _splitmix64(args.seed + length * 1_000_003)
            results.append(
                monte_carlo_estimate(
                    length,
                    args.samples,
                    args.blocks,
                    geometry_seed,
                    args.workers,
                )
            )

    payload = {
        "schema_version": 1,
        "estimator": "p=1/2 Bernoulli score with wrapping-difference observable",
        "random_generator": "Python random.Random (MT19937), independent SplitMix64 block seeds",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "implemented_models": ["square_bond_square_torus"],
        "unimplemented_models": [
            "union_jack_site_square_torus",
            "triangular_site_rhombic_torus",
        ],
        "results": results,
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    sys.stdout.write(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
