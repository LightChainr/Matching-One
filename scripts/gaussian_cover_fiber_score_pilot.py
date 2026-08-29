#!/usr/bin/env python3
"""Fixed-p norm-two Gaussian-cover fiber-score pilot for Issue #226.

The child sites are cyclically labelled by ``j``.  For the cover obtained by
multiplication by ``1+i``, the two sites above a parent site are ``j`` and
``j+N``.  This makes the deck-trivial and deck-detail score directions exact:

    S_+ = sum_j (X_j-p),    S_- = sum_{j<N} (X_j-X_{j+N}).

No child observations are collapsed to a fabricated parent configuration.
The script only measures likelihood-score responses on the child torus.
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from fractions import Fraction
import json
from math import sqrt
from pathlib import Path
import sys
from typing import Callable, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))

from integer_period_torus import (  # noqa: E402
    IntegerTorusGeometry,
    classify_configuration,
    gaussian_integer_torus,
)


MASK64 = (1 << 64) - 1
DEFAULT_P = 0.592746050790
PRIMITIVE_NAMES = (
    "lineage_8_1_trivial",
    "lineage_8_1_detail",
    "lineage_7_4_trivial",
    "lineage_7_4_detail",
)
DERIVED_NAMES = (
    "global_matching_odd_trivial",
    "global_matching_odd_detail",
    "orientation_H4_trivial",
    "orientation_H4_detail",
)
DERIVED_TRANSFORM = (
    (0.5, 0.0, 0.5, 0.0),
    (0.0, 0.5, 0.0, 0.5),
    (1.0, 0.0, -1.0, 0.0),
    (0.0, 1.0, 0.0, -1.0),
)


@dataclass(frozen=True)
class Lineage:
    parent: tuple[int, int]
    child: tuple[int, int]
    parent_order: int
    parent_multiplier: int


LINEAGES = (
    Lineage((8, 1), (7, 9), 65, 29),
    Lineage((7, 4), (3, 11), 65, 24),
)


def splitmix64(value: int) -> int:
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return (value ^ (value >> 31)) & MASK64


def counter_uniform(seed: int, n: int, replica: int, site: int) -> float:
    key = (
        seed
        ^ splitmix64(n)
        ^ splitmix64(replica + 0xD1B54A32D192ED03)
        ^ splitmix64(site + 0x94D049BB133111EB)
    )
    return (splitmix64(key) >> 11) * (2.0 ** -53)


def cyclic_label(a: int, b: int, coordinate: tuple[int, int]) -> int:
    """Cyclic quotient label used by ``gaussian_orientation_mc.cpp``."""

    x, y = coordinate
    return (a * x + b * y) % (a * a + b * b)


def mask_from_labels(
    geometry: IntegerTorusGeometry, a: int, b: int, labels: Sequence[bool]
) -> list[bool]:
    if len(labels) != geometry.n:
        raise ValueError("label field has wrong order")
    return [labels[cyclic_label(a, b, point)] for point in geometry.coordinates]


def fiber_parent_label(lineage: Lineage, child_label: int) -> int:
    return lineage.parent_multiplier * (child_label % lineage.parent_order) % lineage.parent_order


def validate_lineage(lineage: Lineage) -> None:
    a, b = lineage.parent
    c, d = lineage.child
    n = a * a + b * b
    if lineage.parent_order != n or c != a - b or d != a + b:
        raise ValueError("lineage is not multiplication by 1+i")
    if (lineage.parent_multiplier * c - a) % n:
        raise ValueError("child-to-parent multiplier fails x generator")
    if (lineage.parent_multiplier * d - b) % n:
        raise ValueError("child-to-parent multiplier fails y generator")
    for child_label in range(2 * n):
        if fiber_parent_label(lineage, child_label) != fiber_parent_label(
            lineage, (child_label + n) % (2 * n)
        ):
            raise ValueError("deck translate does not preserve parent label")


def matching_odd_channels(
    geometry: IntegerTorusGeometry,
    a: int,
    b: int,
    black_labels: Sequence[bool],
) -> tuple[int, int, int]:
    black = mask_from_labels(geometry, a, b, black_labels)
    white = [not value for value in black]
    primal, _ = classify_configuration(geometry, black)
    matching, _ = classify_configuration(geometry, white, matching=True)
    return (
        int(primal.cross) - int(matching.cross),
        int(primal.direction_0) - int(matching.direction_0),
        int(primal.direction_1) - int(matching.direction_1),
    )


def child_observables(
    geometry: IntegerTorusGeometry,
    a: int,
    b: int,
    black_labels: Sequence[bool],
) -> tuple[int, int]:
    cross, direction_0, direction_1 = matching_odd_channels(
        geometry, a, b, black_labels
    )
    return cross, direction_0 - direction_1


def score_values(black_labels: Sequence[bool], p: float) -> tuple[float, int]:
    n = len(black_labels)
    if n % 2:
        raise ValueError("norm-two cover must have even order")
    half = n // 2
    trivial = sum(int(value) - p for value in black_labels)
    detail = sum(
        int(black_labels[j]) - int(black_labels[j + half]) for j in range(half)
    )
    return trivial, detail


def _run_batch(task: tuple[int, int, int, float, int]) -> dict[str, float | int]:
    batch_index, start, samples, p, seed = task
    geometries = [gaussian_integer_torus(*lineage.child) for lineage in LINEAGES]
    sums: dict[str, float | int] = {
        "batch": batch_index,
        "start_replica": start,
        "samples": samples,
        "sum_m1": 0,
        "sum_m2": 0,
        "sum_score_trivial": 0.0,
        "sum_score_detail": 0,
        "sum_m1_score_trivial": 0.0,
        "sum_m1_score_detail": 0,
        "sum_m2_score_trivial": 0.0,
        "sum_m2_score_detail": 0,
    }
    n = geometries[0].n
    for replica in range(start, start + samples):
        labels = [counter_uniform(seed, n, replica, site) < p for site in range(n)]
        trivial, detail = score_values(labels, p)
        m1 = child_observables(geometries[0], *LINEAGES[0].child, labels)[0]
        m2 = child_observables(geometries[1], *LINEAGES[1].child, labels)[0]
        sums["sum_m1"] += m1
        sums["sum_m2"] += m2
        sums["sum_score_trivial"] += trivial
        sums["sum_score_detail"] += detail
        sums["sum_m1_score_trivial"] += m1 * trivial
        sums["sum_m1_score_detail"] += m1 * detail
        sums["sum_m2_score_trivial"] += m2 * trivial
        sums["sum_m2_score_detail"] += m2 * detail
    return sums


def covariance_of_mean(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    count = len(rows)
    width = len(rows[0])
    means = [sum(row[j] for row in rows) / count for j in range(width)]
    if count < 2:
        return [[0.0] * width for _ in range(width)]
    return [
        [
            sum((row[i] - means[i]) * (row[j] - means[j]) for row in rows)
            / (count * (count - 1))
            for j in range(width)
        ]
        for i in range(width)
    ]


def transform_vector(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> list[float]:
    return [sum(row[j] * vector[j] for j in range(len(vector))) for row in matrix]


def transform_covariance(
    matrix: Sequence[Sequence[float]], covariance: Sequence[Sequence[float]]
) -> list[list[float]]:
    width = len(matrix)
    return [
        [
            sum(
                matrix[i][a] * covariance[a][b] * matrix[j][b]
                for a in range(width)
                for b in range(width)
            )
            for j in range(width)
        ]
        for i in range(width)
    ]


def summarize_batches(batches: Sequence[dict[str, float | int]], p: float) -> dict:
    pq = p * (1.0 - p)
    response_rows = []
    for batch in batches:
        scale = float(batch["samples"]) * pq
        response_rows.append(
            [
                float(batch["sum_m1_score_trivial"]) / scale,
                float(batch["sum_m1_score_detail"]) / scale,
                float(batch["sum_m2_score_trivial"]) / scale,
                float(batch["sum_m2_score_detail"]) / scale,
            ]
        )
    primitive = [sum(row[j] for row in response_rows) / len(response_rows) for j in range(4)]
    primitive_cov = covariance_of_mean(response_rows)
    derived = transform_vector(DERIVED_TRANSFORM, primitive)
    derived_cov = transform_covariance(DERIVED_TRANSFORM, primitive_cov)
    return {
        "primitive_order": list(PRIMITIVE_NAMES),
        "primitive_point": primitive,
        "primitive_standard_error": [sqrt(max(0.0, primitive_cov[j][j])) for j in range(4)],
        "primitive_covariance_of_mean": primitive_cov,
        "derived_order": list(DERIVED_NAMES),
        "derived_point": derived,
        "derived_standard_error": [sqrt(max(0.0, derived_cov[j][j])) for j in range(4)],
        "derived_covariance_of_mean": derived_cov,
        "batch_response_rows": response_rows,
    }


def run_pilot(samples: int, batches: int, p: float, seed: int, workers: int) -> tuple[list[dict], dict]:
    if samples <= 0 or batches <= 1 or samples % batches:
        raise ValueError("samples must be positive and divisible by at least two batches")
    for lineage in LINEAGES:
        validate_lineage(lineage)
    per_batch = samples // batches
    tasks = [
        (batch, batch * per_batch, per_batch, p, seed) for batch in range(batches)
    ]
    if workers == 1:
        rows = [_run_batch(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(_run_batch, tasks))
    return rows, summarize_batches(rows, p)


def _active_from_mask(mask: int, n: int) -> list[bool]:
    return [bool(mask & (1 << site)) for site in range(n)]


def _fraction_expectation(
    observable: Callable[[Sequence[bool]], int],
    base_p: Fraction,
    epsilon: Fraction,
    signs: Sequence[int],
) -> Fraction:
    total = Fraction(0)
    n = len(signs)
    probabilities = [base_p + epsilon * sign for sign in signs]
    for mask in range(1 << n):
        active = _active_from_mask(mask, n)
        weight = Fraction(1)
        for state, probability in zip(active, probabilities):
            weight *= probability if state else 1 - probability
        total += observable(active) * weight
    return total


def _lagrange_derivative_weights(nodes: Sequence[Fraction]) -> list[Fraction]:
    weights = []
    for j, node in enumerate(nodes):
        polynomial = [Fraction(1)]
        denominator = Fraction(1)
        for k, other in enumerate(nodes):
            if k == j:
                continue
            next_polynomial = [Fraction(0)] * (len(polynomial) + 1)
            for degree, coefficient in enumerate(polynomial):
                next_polynomial[degree] -= other * coefficient
                next_polynomial[degree + 1] += coefficient
            polynomial = next_polynomial
            denominator *= node - other
        weights.append(polynomial[1] / denominator)
    return weights


def exact_tiny_oracle() -> dict:
    """N=5 -> 10 exact score identity and independent polynomial derivative."""

    lineage = Lineage((2, 1), (1, 3), 5, 2)
    validate_lineage(lineage)
    geometry = gaussian_integer_torus(*lineage.child)
    p = Fraction(2, 5)
    signs_by_mode = {
        "trivial": [1] * 10,
        "detail": [1] * 5 + [-1] * 5,
    }

    def cross(active: Sequence[bool]) -> int:
        return child_observables(geometry, *lineage.child, active)[0]

    def orientation(active: Sequence[bool]) -> int:
        return child_observables(geometry, *lineage.child, active)[1]

    observables = {"matching_odd_cross": cross, "directional_matching_odd_H4": orientation}
    nodes = [Fraction(index, 100) for index in range(-5, 6)]
    weights = _lagrange_derivative_weights(nodes)
    checks = {}
    for observable_name, observable in observables.items():
        for mode, signs in signs_by_mode.items():
            score_numerator = Fraction(0)
            for mask in range(1 << 10):
                active = _active_from_mask(mask, 10)
                occupied = sum(active)
                weight = p**occupied * (1 - p) ** (10 - occupied)
                raw_score = sum(sign * (int(state) - p) for sign, state in zip(signs, active))
                score_numerator += observable(active) * raw_score * weight
            score_derivative = score_numerator / (p * (1 - p))
            finite_difference = sum(
                coefficient * _fraction_expectation(observable, p, epsilon, signs)
                for coefficient, epsilon in zip(weights, nodes)
            )
            key = f"{observable_name}:{mode}"
            checks[key] = {
                "score_derivative": str(score_derivative),
                "polynomial_finite_difference": str(finite_difference),
                "equal": score_derivative == finite_difference,
            }
            if score_derivative != finite_difference:
                raise AssertionError(f"tiny score identity failed for {key}")
    return {
        "parent": [2, 1],
        "child": [1, 3],
        "parent_order": 5,
        "child_order": 10,
        "child_to_parent_multiplier_mod_5": 2,
        "p": str(p),
        "finite_difference_nodes": [str(value) for value in nodes],
        "checks": checks,
    }


def write_batches(path: Path, batches: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(batches[0]))
        writer.writeheader()
        writer.writerows(batches)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=20_000)
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--p", type=float, default=DEFAULT_P)
    parser.add_argument("--seed", type=int, default=226_20260829)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batches-output", type=Path)
    args = parser.parse_args()
    if not 0.0 < args.p < 1.0:
        raise ValueError("p must lie strictly between zero and one")
    if args.samples > 20_000:
        raise ValueError("this pilot is frozen to at most 20,000 samples")
    batches, summary = run_pilot(args.samples, args.batches, args.p, args.seed, args.workers)
    batch_path = args.batches_output or args.output.with_suffix(".batches.csv")
    write_batches(batch_path, batches)
    payload = {
        "schema": "matching-one.gaussian-cover-fiber-score-pilot.v1",
        "issue": 226,
        "status": "smoke_only_not_production_evidence",
        "design": {
            "parent_order": 65,
            "child_order": 130,
            "cover_multiplier": "1+i",
            "deck_group": "Z/2",
            "lineages": [
                {
                    "parent": list(lineage.parent),
                    "child": list(lineage.child),
                    "child_to_parent_multiplier_mod_65": lineage.parent_multiplier,
                    "fiber": "{j,j+65}",
                }
                for lineage in LINEAGES
            ],
            "trivial_score": "sum_j (X_j-p)",
            "detail_score": "sum_{j=0}^{64} (X_j-X_{j+65})",
            "primitive_observable": "primal.cross(X)-matching.cross(1-X)",
            "global_observable": "(M_(7,9)+M_(3,11))/2",
            "orientation_H4_observable": "M_(7,9)-M_(3,11)",
            "parent_configuration_constructed": False,
        },
        "run": {
            "samples": args.samples,
            "batches": args.batches,
            "p": args.p,
            "seed": args.seed,
            "workers": args.workers,
            "batch_sufficient_statistics": str(batch_path),
        },
        "exact_score_oracle": exact_tiny_oracle(),
        "response": summary,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["response"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
