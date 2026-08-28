#!/usr/bin/env python3
"""Euler and local-motif control variates for square-site matching.

For a complementary black-primal/white-matching configuration on a periodic
square lattice, the uncentered configuration identity used here is

    N_B - N_W = q + V - E + F0,

where ``q`` is the cross-wrapping difference, ``V`` is occupied sites, ``E``
is occupied nearest-neighbour edges, and ``F0`` is fully occupied elementary
faces.  Subtracting their known Bernoulli means gives the matching-function
identity in issue #34.

The Monte Carlo mode learns coefficients for exact zero-mean Euler and local
motif controls on an independent pilot, freezes them, and reports performance
on a fresh evaluation stream relative to the better of ``q`` and the centered
cluster-number estimator.  It is a correctness/power reference, not a
production Newman-Ziff kernel.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
import math
from pathlib import Path
import platform
import random
from typing import Dict, Iterable, Optional, Sequence, Tuple

from control_variate_estimator import (
    FrozenZeroMeanControls,
    sample_covariance,
)
from integer_period_torus import (
    IntegerTorusGeometry,
    axis_integer_torus,
    classify_configuration,
    diamond_integer_torus,
    gaussian_integer_torus,
)


MASK64 = (1 << 64) - 1
CONTROL_NAMES = ("Z_V", "Z_E", "Z_F0", "Z_diag_11", "Z_corner_3")
CONTROL_BASES = (
    ("euler", (0, 1, 2)),
    ("euler_plus_diagonal", (0, 1, 2, 3)),
    ("euler_plus_local_motifs", (0, 1, 2, 3, 4)),
)


def falling(value: int, order: int) -> int:
    """Return the falling factorial ``(value)_order``."""

    if order < 0:
        raise ValueError("falling-factorial order must be nonnegative")
    result = 1
    for offset in range(order):
        result *= value - offset
    return result


@dataclass(frozen=True)
class EulerConvention:
    geometry: IntegerTorusGeometry
    faces: Tuple[Tuple[int, ...], ...]
    diagonal_11: Tuple[Tuple[int, ...], ...]
    corner_3: Tuple[Tuple[int, ...], ...]


@dataclass(frozen=True)
class EulerObservables:
    occupied_sites: int
    occupied_nn_edges: int
    occupied_faces: int
    occupied_diagonal_11: int
    occupied_corner_3: int
    black_clusters: int
    white_matching_clusters: int
    q_cross: int

    @property
    def cluster_difference(self) -> int:
        return self.black_clusters - self.white_matching_clusters


def _motif_embeddings(
    geometry: IntegerTorusGeometry,
    offsets: Sequence[Tuple[int, int]],
    name: str,
) -> Tuple[Tuple[int, ...], ...]:
    embeddings = []
    for x, y in geometry.coordinates:
        vertices = tuple(
            geometry.vertex((x + dx, y + dy)) for dx, dy in offsets
        )
        if len(set(vertices)) != len(vertices):
            raise ValueError(
                f"{geometry.name}: motif {name} identifies distinct declared "
                "vertices on this quotient"
            )
        embeddings.append(vertices)
    return tuple(embeddings)


def prepare_convention(geometry: IntegerTorusGeometry) -> EulerConvention:
    """Freeze the local square-cell and auxiliary motif conventions."""

    return EulerConvention(
        geometry=geometry,
        faces=_motif_embeddings(
            geometry, ((0, 0), (1, 0), (0, 1), (1, 1)), "F0"
        ),
        # One declared diagonal orientation.  Its N embeddings have exact
        # Bernoulli mean N p^2 and fixed-K mean N (K)_2/(N)_2.
        diagonal_11=_motif_embeddings(
            geometry, ((0, 0), (1, 1)), "diagonal_11"
        ),
        # One declared right-angle three-site corner per square anchor.
        corner_3=_motif_embeddings(
            geometry, ((0, 0), (1, 0), (0, 1)), "corner_3"
        ),
    )


def _occupied_embeddings(
    active: Sequence[bool], embeddings: Iterable[Sequence[int]]
) -> int:
    return sum(all(active[vertex] for vertex in motif) for motif in embeddings)


def configuration_observables(
    convention: EulerConvention, active: Sequence[bool]
) -> EulerObservables:
    geometry = convention.geometry
    if len(active) != geometry.n:
        raise ValueError("active mask length does not match geometry")
    black_channels, black_components = classify_configuration(geometry, active)
    white = [not value for value in active]
    white_channels, white_components = classify_configuration(
        geometry, white, matching=True
    )
    return EulerObservables(
        occupied_sites=sum(active),
        occupied_nn_edges=sum(
            active[edge.i] and active[edge.j] for edge in geometry.primal_edges
        ),
        occupied_faces=_occupied_embeddings(active, convention.faces),
        occupied_diagonal_11=_occupied_embeddings(
            active, convention.diagonal_11
        ),
        occupied_corner_3=_occupied_embeddings(active, convention.corner_3),
        black_clusters=len(black_components),
        white_matching_clusters=len(white_components),
        q_cross=int(black_channels.cross) - int(white_channels.cross),
    )


def euler_identity_residual(observables: EulerObservables) -> int:
    """Return zero iff the uncentered Euler identity holds."""

    right = (
        observables.q_cross
        + observables.occupied_sites
        - observables.occupied_nn_edges
        + observables.occupied_faces
    )
    return observables.cluster_difference - right


def canonical_channels(
    observables: EulerObservables, site_count: int, p: float
) -> Tuple[float, list[float]]:
    """Return ``q`` and canonically centered zero-mean controls."""

    return float(observables.q_cross), [
        observables.occupied_sites - site_count * p,
        observables.occupied_nn_edges - 2 * site_count * p**2,
        observables.occupied_faces - site_count * p**4,
        observables.occupied_diagonal_11 - site_count * p**2,
        observables.occupied_corner_3 - site_count * p**3,
    ]


def centered_cluster_estimator(
    observables: EulerObservables, site_count: int, p: float
) -> float:
    chi = p - 2 * p**2 + p**4
    return observables.cluster_difference - site_count * chi


def microcanonical_control_means(
    site_count: int, occupied: int
) -> Dict[str, Fraction]:
    """Exact conditional motif means under uniform fixed-K sampling."""

    if not 0 <= occupied <= site_count:
        raise ValueError("occupied count lies outside [0,N]")

    def mean(embeddings: int, motif_sites: int) -> Fraction:
        denominator = falling(site_count, motif_sites)
        if denominator == 0:
            raise ValueError("motif has more distinct sites than the quotient")
        return Fraction(
            embeddings * falling(occupied, motif_sites), denominator
        )

    return {
        "V": Fraction(occupied, 1),
        "E": mean(2 * site_count, 2),
        "F0": mean(site_count, 4),
        "diag_11": mean(site_count, 2),
        "corner_3": mean(site_count, 3),
    }


def _fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def exhaustive_validation(geometry: IntegerTorusGeometry) -> dict[str, object]:
    """Verify identity and all conditional means on every tiny configuration."""

    if geometry.n > 16:
        raise ValueError("exhaustive validation is limited to N<=16")
    convention = prepare_convention(geometry)
    fields = {
        "V": "occupied_sites",
        "E": "occupied_nn_edges",
        "F0": "occupied_faces",
        "diag_11": "occupied_diagonal_11",
        "corner_3": "occupied_corner_3",
    }
    counts = [0] * (geometry.n + 1)
    sums = {
        name: [0] * (geometry.n + 1) for name in fields
    }
    maximum_residual = 0
    for mask in range(1 << geometry.n):
        active = [bool((mask >> vertex) & 1) for vertex in range(geometry.n)]
        observed = configuration_observables(convention, active)
        residual = euler_identity_residual(observed)
        maximum_residual = max(maximum_residual, abs(residual))
        occupied = observed.occupied_sites
        counts[occupied] += 1
        for name, attribute in fields.items():
            sums[name][occupied] += int(getattr(observed, attribute))

    conditional = []
    for occupied, configuration_count in enumerate(counts):
        expected = microcanonical_control_means(geometry.n, occupied)
        row_observed = {
            name: Fraction(sums[name][occupied], configuration_count)
            for name in fields
        }
        residuals = {
            name: row_observed[name] - expected[name] for name in fields
        }
        conditional.append(
            {
                "K": occupied,
                "configurations": configuration_count,
                "observed_means": {
                    name: _fraction_text(value)
                    for name, value in row_observed.items()
                },
                "expected_means": {
                    name: _fraction_text(value)
                    for name, value in expected.items()
                },
                "centered_sum_residuals": {
                    name: _fraction_text(value)
                    for name, value in residuals.items()
                },
            }
        )
        if any(residuals.values()):
            raise AssertionError(
                f"{geometry.name}: fixed-K motif mean failed at K={occupied}"
            )
    if maximum_residual:
        raise AssertionError(
            f"{geometry.name}: Euler identity maximum residual {maximum_residual}"
        )
    return {
        "geometry": geometry.name,
        "period_matrix_columns": [list(row) for row in geometry.periods.matrix],
        "N": geometry.n,
        "configurations": 1 << geometry.n,
        "q_convention": "black_cross_wrap_minus_white_matching_cross_wrap",
        "edge_count_convention": "2*N lifted positive NN edges",
        "face_convention": "N anchors with corners (0,0),(1,0),(0,1),(1,1)",
        "euler_identity_max_abs_residual": maximum_residual,
        "all_fixed_K_centered_sums_zero": True,
        "conditional_means": conditional,
    }


def exact_validation_suite() -> dict[str, object]:
    geometries = (
        axis_integer_torus(3),
        diamond_integer_torus(2),
        gaussian_integer_torus(2, 1),
    )
    return {
        "schema_version": 1,
        "identity": "N_B-N_W=q+V-E+F0",
        "geometries": [exhaustive_validation(geometry) for geometry in geometries],
    }


def _splitmix64(value: int) -> int:
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return (value ^ (value >> 31)) & MASK64


def _sample_rows(
    convention: EulerConvention, p: float, samples: int, seed: int
) -> Tuple[list[float], list[float], list[list[float]]]:
    rng = random.Random(seed)
    topology = []
    cluster = []
    controls = []
    for _ in range(samples):
        active = [rng.random() < p for _ in range(convention.geometry.n)]
        observed = configuration_observables(convention, active)
        if euler_identity_residual(observed):
            raise AssertionError("Euler identity failed during Monte Carlo")
        q, centered = canonical_channels(observed, convention.geometry.n, p)
        topology.append(q)
        cluster.append(
            centered_cluster_estimator(observed, convention.geometry.n, p)
        )
        controls.append(centered)
    return topology, cluster, controls


def _mean_variance(values: Sequence[float]) -> Tuple[float, float]:
    mean = math.fsum(values) / len(values)
    variance = math.fsum((value - mean) ** 2 for value in values) / (
        len(values) - 1
    )
    return mean, variance


def pilot_frozen_monte_carlo(
    geometry: IntegerTorusGeometry,
    p: float,
    pilot_samples: int,
    evaluation_samples: int,
    seed: int,
    ridge: float = 0.0,
) -> dict[str, object]:
    if not 0.0 <= p <= 1.0:
        raise ValueError("p must lie in [0,1]")
    if pilot_samples < 2 or evaluation_samples < 2:
        raise ValueError("pilot and evaluation each require at least two samples")
    convention = prepare_convention(geometry)
    pilot_seed = _splitmix64(seed)
    evaluation_seed = _splitmix64(seed + 1)
    pilot_q, _pilot_cluster, pilot_controls = _sample_rows(
        convention, p, pilot_samples, pilot_seed
    )
    evaluation_q, evaluation_cluster, evaluation_controls = _sample_rows(
        convention, p, evaluation_samples, evaluation_seed
    )
    q_mean, q_variance = _mean_variance(evaluation_q)
    cluster_mean, cluster_variance = _mean_variance(evaluation_cluster)
    best_name, best_variance = min(
        (("q_cross", q_variance), ("D_cluster", cluster_variance)),
        key=lambda item: item[1],
    )
    joint_rows = [
        [q, cluster] + controls
        for q, cluster, controls in zip(
            evaluation_q, evaluation_cluster, evaluation_controls
        )
    ]
    control_means = [
        math.fsum(row[index] for row in evaluation_controls) / evaluation_samples
        for index in range(len(CONTROL_NAMES))
    ]
    hierarchy: Dict[str, dict[str, object]] = {}
    frozen_full: Optional[FrozenZeroMeanControls] = None
    adjusted_full: Optional[dict[str, object]] = None
    for basis_name, indices in CONTROL_BASES:
        names = tuple(CONTROL_NAMES[index] for index in indices)
        pilot_subset = [
            [row[index] for index in indices] for row in pilot_controls
        ]
        evaluation_subset = [
            [row[index] for index in indices] for row in evaluation_controls
        ]
        frozen = FrozenZeroMeanControls.fit(
            names, pilot_q, pilot_subset, ridge
        )
        adjusted = frozen.evaluate(evaluation_q, evaluation_subset)
        adjusted_variance = float(adjusted["sample_variance"])
        hierarchy[basis_name] = {
            "control_names": list(names),
            "frozen_coefficients": list(frozen.coefficients),
            "coefficient_l1_norm": math.fsum(
                abs(value) for value in frozen.coefficients
            ),
            "coefficient_l2_norm": math.sqrt(
                math.fsum(value * value for value in frozen.coefficients)
            ),
            "applied_diagonal_ridge": frozen.applied_ridge,
            "evaluation": adjusted,
            "variance_reduction_vs_best_single": (
                best_variance / adjusted_variance if adjusted_variance else None
            ),
        }
        if basis_name == "euler_plus_local_motifs":
            frozen_full = frozen
            adjusted_full = adjusted
    if frozen_full is None or adjusted_full is None:
        raise AssertionError("full control basis is missing")
    adjusted_variance = float(adjusted_full["sample_variance"])
    return {
        "schema_version": 1,
        "method": "independent_pilot_frozen_zero_mean_controls",
        "wrapping_channel_policy": (
            "retain q_cross only; configuration-identical wrapping-difference "
            "channels are excluded rather than ridge-combined"
        ),
        "geometry": geometry.name,
        "period_matrix_columns": [list(row) for row in geometry.periods.matrix],
        "N": geometry.n,
        "p": p,
        "pilot_samples": pilot_samples,
        "evaluation_samples": evaluation_samples,
        "base_seed": seed,
        "pilot_seed": pilot_seed,
        "evaluation_seed": evaluation_seed,
        "rng": "Python random.Random (MT19937)",
        "control_names": list(CONTROL_NAMES),
        "frozen_coefficients": list(frozen_full.coefficients),
        "coefficient_l1_norm": math.fsum(
            abs(value) for value in frozen_full.coefficients
        ),
        "coefficient_l2_norm": math.sqrt(
            math.fsum(value * value for value in frozen_full.coefficients)
        ),
        "applied_diagonal_ridge": frozen_full.applied_ridge,
        "pilot_joint_covariance_q_then_controls": [
            list(row) for row in frozen_full.pilot_covariance
        ],
        "evaluation_joint_covariance_q_cluster_then_controls": sample_covariance(
            joint_rows
        ),
        "evaluation_control_means": dict(zip(CONTROL_NAMES, control_means)),
        "predeclared_control_hierarchy": hierarchy,
        "estimators": {
            "q_cross": {
                "mean": q_mean,
                "sample_variance": q_variance,
                "standard_error": math.sqrt(q_variance / evaluation_samples),
            },
            "D_cluster": {
                "mean": cluster_mean,
                "sample_variance": cluster_variance,
                "standard_error": math.sqrt(
                    cluster_variance / evaluation_samples
                ),
            },
            "pilot_frozen_control_variate": adjusted_full,
        },
        "best_single_estimator": best_name,
        "best_single_sample_variance": best_variance,
        "variance_reduction_vs_best_single": (
            best_variance / adjusted_variance if adjusted_variance else None
        ),
        "python": platform.python_version(),
        "platform": platform.platform(),
    }


def _geometry_from_args(args: argparse.Namespace) -> IntegerTorusGeometry:
    if args.geometry == "axis":
        if args.L is None:
            raise ValueError("axis geometry requires --L")
        return axis_integer_torus(args.L)
    if args.geometry == "diamond":
        if args.L is None:
            raise ValueError("diamond geometry requires --L")
        return diamond_integer_torus(args.L)
    if args.a is None or args.b is None:
        raise ValueError("gaussian geometry requires --a and --b")
    return gaussian_integer_torus(args.a, args.b)


def _parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-suite", action="store_true")
    parser.add_argument("--geometry", choices=("axis", "diamond", "gaussian"))
    parser.add_argument("--L", type=int)
    parser.add_argument("--a", type=int)
    parser.add_argument("--b", type=int)
    parser.add_argument("--p", type=float, default=0.59274605079210)
    parser.add_argument("--pilot-samples", type=int, default=20000)
    parser.add_argument("--evaluation-samples", type=int, default=80000)
    parser.add_argument("--seed", type=int, default=20260828)
    parser.add_argument("--ridge", type=float, default=0.0)
    parser.add_argument("--output", type=Path)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parse_args(argv)
    if args.exact_suite:
        payload = exact_validation_suite()
    else:
        if args.geometry is None:
            raise SystemExit("Monte Carlo mode requires --geometry")
        try:
            geometry = _geometry_from_args(args)
            payload = pilot_frozen_monte_carlo(
                geometry,
                args.p,
                args.pilot_samples,
                args.evaluation_samples,
                args.seed,
                args.ridge,
            )
        except (ValueError, ArithmeticError) as exc:
            raise SystemExit(str(exc)) from exc
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
