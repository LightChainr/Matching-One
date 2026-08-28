#!/usr/bin/env python3
"""Exact-marginal order-statistic coupling pilot for Gaussian graph covers.

For every degree-Q cover fiber, start from Q iid child priorities.  If U_(k)
is the k-th order statistic, then

    V_k = F_Beta(k,Q+1-k)(U_(k))

is exactly Uniform(0,1).  Disjoint fibers are independent, so every fixed k
gives an iid parent priority field and hence an exact uniform parent
Newman--Ziff permutation.  The Q parent fields are correlated with each other;
they are Q exact *marginal* couplings, not Q mutually independent fields.

The executable pilot uses the frozen norm-5 N=65 -> 325 lineage and compares
the frozen H4 semigroup residual with an independent-parent baseline.  It is a
small method pilot, not a physics score and not production evidence.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
import time
from pathlib import Path
from typing import Iterable, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from gaussian_covering_map import (  # noqa: E402
    CoveringMap,
    GaussianPair,
    canonical_cover,
)
from integer_period_torus import gaussian_integer_torus  # noqa: E402
from threshold_rank_nz import ThresholdRankEngine  # noqa: E402


P_REF = 0.592746050790
H4_NORM5_RAW_RATIO = (-14.0 / 25.0) * 5.0 ** (-13.0 / 8.0)


def beta_order_cdf(order: int, degree: int, value: float) -> float:
    """CDF of Beta(order, degree+1-order), for integer parameters.

    This binomial-tail polynomial avoids a scipy dependency and is exact up to
    ordinary floating-point evaluation for the small degrees used here.
    """

    if not 1 <= order <= degree:
        raise ValueError("require 1 <= order <= degree")
    if not 0.0 <= value <= 1.0:
        raise ValueError("value must lie in [0,1]")
    return sum(
        math.comb(degree, count)
        * value**count
        * (1.0 - value) ** (degree - count)
        for count in range(order, degree + 1)
    )


def orderstat_parent_fields(
    fiber_uniforms: Sequence[Sequence[float]],
) -> tuple[tuple[float, ...], ...]:
    """Return Q exact-marginal iid parent fields from child fiber uniforms.

    Input is indexed ``[parent_label][kernel_index]``.  Output is indexed
    ``[order-1][parent_label]``.  Fields at different order indices are
    generally dependent and must not be described as mutually iid.
    """

    if not fiber_uniforms:
        raise ValueError("at least one fiber is required")
    degree = len(fiber_uniforms[0])
    if degree == 0 or any(len(row) != degree for row in fiber_uniforms):
        raise ValueError("all nonempty fibers must have the same degree")
    output = [[] for _ in range(degree)]
    for fiber in fiber_uniforms:
        ordered = sorted(fiber)
        for index, value in enumerate(ordered, start=1):
            output[index - 1].append(beta_order_cdf(index, degree, value))
    return tuple(tuple(field) for field in output)


def priorities_to_permutation(priorities: Sequence[float]) -> tuple[int, ...]:
    if len(set(priorities)) != len(priorities):
        raise ValueError("priorities must be distinct")
    return tuple(sorted(range(len(priorities)), key=priorities.__getitem__))


def child_priorities(
    cover: CoveringMap, fiber_uniforms: Sequence[Sequence[float]]
) -> tuple[float, ...]:
    if len(fiber_uniforms) != cover.parent.n:
        raise ValueError("one uniform fiber is required per parent label")
    if any(len(row) != cover.degree for row in fiber_uniforms):
        raise ValueError("fiber width must equal cover degree")
    priorities = [0.0] * cover.child.n
    for parent_label, values in enumerate(fiber_uniforms):
        for kernel_index, child_label in enumerate(cover.fiber(parent_label)):
            priorities[child_label] = values[kernel_index]
    return tuple(priorities)


def binomial_probabilities(n: int, p: float) -> tuple[float, ...]:
    probabilities = [0.0] * (n + 1)
    probabilities[0] = (1.0 - p) ** n
    ratio = p / (1.0 - p)
    for count in range(n):
        probabilities[count + 1] = (
            probabilities[count] * (n - count) / (count + 1) * ratio
        )
    total = sum(probabilities)
    return tuple(value / total for value in probabilities)


def matching_value_from_ranks(
    k_minus: int, k_plus: int, binomial: Sequence[float]
) -> float:
    """Rao--Blackwellized fixed-p matching value for one permutation."""

    return sum(binomial[k_plus:]) - sum(binomial[:k_minus])


def orientation_difference(
    permutation: Sequence[int],
    first: ThresholdRankEngine,
    second: ThresholdRankEngine,
    binomial: Sequence[float],
) -> float:
    first_ranks = first.threshold_ranks(permutation)
    second_ranks = second.threshold_ranks(permutation)
    return matching_value_from_ranks(*first_ranks, binomial) - matching_value_from_ranks(
        *second_ranks, binomial
    )


def _variance(values: Sequence[float]) -> float:
    return statistics.variance(values) if len(values) >= 2 else math.nan


def _covariance(first: Sequence[float], second: Sequence[float]) -> float:
    return statistics.covariance(first, second)


def _solve(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> list[float]:
    """Small dense solve with pivoting; sufficient for Q<=5 pilot weights."""

    size = len(rhs)
    augmented = [list(matrix[row]) + [rhs[row]] for row in range(size)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1e-15:
            raise ValueError("singular pilot covariance matrix")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                left - factor * right
                for left, right in zip(augmented[row], augmented[column])
            ]
    return [augmented[row][-1] for row in range(size)]


def constrained_residual_weights(
    child: Sequence[float],
    parents: Sequence[Sequence[float]],
    multiplier: float,
    ridge_fraction: float = 1e-8,
) -> tuple[float, ...]:
    """Train sum-to-one weights minimizing Var(child-r*w'parent)."""

    if multiplier == 0.0:
        raise ValueError("multiplier must be nonzero")
    degree = len(parents)
    covariance = [
        [_covariance(parents[row], parents[column]) for column in range(degree)]
        for row in range(degree)
    ]
    scale = sum(covariance[index][index] for index in range(degree)) / degree
    for index in range(degree):
        covariance[index][index] += ridge_fraction * max(scale, 1e-30)
    target = [_covariance(parents[index], child) / multiplier for index in range(degree)]
    inverse_target = _solve(covariance, target)
    inverse_one = _solve(covariance, [1.0] * degree)
    correction = (1.0 - sum(inverse_target)) / sum(inverse_one)
    return tuple(
        inverse_target[index] + correction * inverse_one[index]
        for index in range(degree)
    )


def _weighted_parent(
    parent_fields: Sequence[Sequence[float]], weights: Sequence[float], index: int
) -> float:
    return sum(weight * parent[index] for weight, parent in zip(weights, parent_fields))


def _paired_variance_comparison(
    candidate: Sequence[float], baseline: Sequence[float]
) -> dict[str, float]:
    """Paired first-order uncertainty for a variance difference."""

    candidate_mean = statistics.mean(candidate)
    baseline_mean = statistics.mean(baseline)
    influence = [
        (left - candidate_mean) ** 2 - (right - baseline_mean) ** 2
        for left, right in zip(candidate, baseline)
    ]
    difference = statistics.mean(influence)
    standard_error = statistics.stdev(influence) / math.sqrt(len(influence))
    return {
        "variance_difference_vs_independent": difference,
        "variance_difference_se": standard_error,
        "variance_difference_z": difference / standard_error,
    }


def run_pilot(samples: int, train_samples: int, seed: int) -> tuple[dict, list[dict]]:
    if not 10 <= train_samples < samples:
        raise ValueError("require 10 <= train_samples < samples")
    parent_pairs = (GaussianPair(8, 1), GaussianPair(7, 4))
    child_pairs = (GaussianPair(17, 6), GaussianPair(18, 1))
    covers = tuple(
        canonical_cover(parent, child)
        for parent, child in zip(parent_pairs, child_pairs)
    )
    degree = covers[0].degree
    if any(cover.degree != degree for cover in covers):
        raise AssertionError("lineage cover degrees disagree")

    parent_engines = tuple(
        ThresholdRankEngine(gaussian_integer_torus(pair.a, pair.b))
        for pair in parent_pairs
    )
    child_engines = tuple(
        ThresholdRankEngine(gaussian_integer_torus(pair.a, pair.b))
        for pair in child_pairs
    )
    parent_binomial = binomial_probabilities(parent_pairs[0].n, P_REF)
    child_binomial = binomial_probabilities(child_pairs[0].n, P_REF)
    generator = random.Random(seed)

    child_values: list[float] = []
    parent_values: list[list[float]] = [[] for _ in range(degree)]
    independent_values: list[float] = []
    rows: list[dict] = []
    started = time.perf_counter()
    for replica in range(samples):
        fibers = tuple(
            tuple(generator.random() for _ in range(degree))
            for _ in range(parent_pairs[0].n)
        )
        child_permutations = tuple(
            priorities_to_permutation(child_priorities(cover, fibers))
            for cover in covers
        )
        child_value = (
            matching_value_from_ranks(
                *child_engines[0].threshold_ranks(child_permutations[0]),
                child_binomial,
            )
            - matching_value_from_ranks(
                *child_engines[1].threshold_ranks(child_permutations[1]),
                child_binomial,
            )
        )
        fields = orderstat_parent_fields(fibers)
        per_order = [
            orientation_difference(
                priorities_to_permutation(field),
                parent_engines[0],
                parent_engines[1],
                parent_binomial,
            )
            for field in fields
        ]
        independent_priorities = tuple(
            generator.random() for _ in range(parent_pairs[0].n)
        )
        independent = orientation_difference(
            priorities_to_permutation(independent_priorities),
            parent_engines[0],
            parent_engines[1],
            parent_binomial,
        )
        child_values.append(child_value)
        independent_values.append(independent)
        for order, value in enumerate(per_order):
            parent_values[order].append(value)
        rows.append(
            {
                "replica": replica,
                "split": "train" if replica < train_samples else "evaluation",
                "child_delta_M": child_value,
                "independent_parent_delta_M": independent,
                **{
                    f"order_{order + 1}_parent_delta_M": value
                    for order, value in enumerate(per_order)
                },
            }
        )

    weights = constrained_residual_weights(
        child_values[:train_samples],
        [values[:train_samples] for values in parent_values],
        H4_NORM5_RAW_RATIO,
    )
    evaluation = range(train_samples, samples)
    baseline = [
        child_values[index] - H4_NORM5_RAW_RATIO * independent_values[index]
        for index in evaluation
    ]
    equal = [
        child_values[index]
        - H4_NORM5_RAW_RATIO
        * sum(values[index] for values in parent_values)
        / degree
        for index in evaluation
    ]
    trained = [
        child_values[index]
        - H4_NORM5_RAW_RATIO * _weighted_parent(parent_values, weights, index)
        for index in evaluation
    ]
    single = [
        [
            child_values[index] - H4_NORM5_RAW_RATIO * parent_values[order][index]
            for index in evaluation
        ]
        for order in range(degree)
    ]
    baseline_variance = _variance(baseline)
    elapsed = time.perf_counter() - started
    evaluation_parent_means = [
        statistics.mean(values[train_samples:]) for values in parent_values
    ]
    evaluation_child = child_values[train_samples:]
    evaluation_independent = independent_values[train_samples:]
    equal_parent = [
        sum(values[index] for values in parent_values) / degree
        for index in evaluation
    ]
    trained_parent = [
        _weighted_parent(parent_values, weights, index) for index in evaluation
    ]
    summary = {
        "schema": "norm5 cover order-statistic coupling pilot v1",
        "status": "method pilot; not a physics score or production evidence",
        "samples": samples,
        "train_samples": train_samples,
        "evaluation_samples": samples - train_samples,
        "seed": seed,
        "p_ref": P_REF,
        "frozen_raw_h4_ratio": H4_NORM5_RAW_RATIO,
        "lineage": {
            "parent": [[pair.a, pair.b] for pair in parent_pairs],
            "child": [[pair.a, pair.b] for pair in child_pairs],
            "cover_degree": degree,
            "cover_units": [cover.t for cover in covers],
            "direction_maps": [cover.direction_map() for cover in covers],
            "homology_matrices": [cover.homology_matrix() for cover in covers],
        },
        "trained_weights": list(weights),
        "weight_sum": sum(weights),
        "evaluation_parent_means": evaluation_parent_means,
        "evaluation_parent_mean_span": max(evaluation_parent_means)
        - min(evaluation_parent_means),
        "evaluation": {
            "independent_baseline": {
                "mean": statistics.mean(baseline),
                "variance": baseline_variance,
                "variance_ratio_vs_independent": 1.0,
            },
            "equal_orderstat_average": {
                "mean": statistics.mean(equal),
                "variance": _variance(equal),
                "variance_ratio_vs_independent": _variance(equal) / baseline_variance,
                "variance_gain_vs_independent": baseline_variance / _variance(equal),
                "child_parent_covariance": _covariance(evaluation_child, equal_parent),
                "child_parent_correlation": statistics.correlation(
                    evaluation_child, equal_parent
                ),
                **_paired_variance_comparison(equal, baseline),
            },
            "trained_orderstat_average": {
                "mean": statistics.mean(trained),
                "variance": _variance(trained),
                "variance_ratio_vs_independent": _variance(trained) / baseline_variance,
                "variance_gain_vs_independent": baseline_variance / _variance(trained),
                "child_parent_covariance": _covariance(evaluation_child, trained_parent),
                "child_parent_correlation": statistics.correlation(
                    evaluation_child, trained_parent
                ),
                **_paired_variance_comparison(trained, baseline),
            },
            "single_order_fields": [
                {
                    "order": order + 1,
                    "mean": statistics.mean(values),
                    "variance": _variance(values),
                    "variance_ratio_vs_independent": _variance(values)
                    / baseline_variance,
                    "child_parent_covariance": _covariance(
                        evaluation_child, parent_values[order][train_samples:]
                    ),
                    "child_parent_correlation": statistics.correlation(
                        evaluation_child, parent_values[order][train_samples:]
                    ),
                    **_paired_variance_comparison(values, baseline),
                }
                for order, values in enumerate(single)
            ],
        },
        "elapsed_seconds": elapsed,
    }
    return summary, rows


def write_outputs(output_dir: Path, summary: dict, rows: Iterable[dict]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "pilot.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    rows = list(rows)
    with (output_dir / "replicas.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=400)
    parser.add_argument("--train-samples", type=int, default=200)
    parser.add_argument("--seed", type=int, default=2026086701)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summary, rows = run_pilot(args.samples, args.train_samples, args.seed)
    write_outputs(args.output_dir, summary, rows)
    print(json.dumps(summary["evaluation"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
