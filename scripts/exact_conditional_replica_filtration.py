#!/usr/bin/env python3
"""Exact N=5/N=10 conditional-replica martingale oracle for Issue #256."""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
import json
from math import atan2, comb, degrees
from pathlib import Path
from typing import Iterable, Sequence

from c4_self_matching_exact import c4_self_matching_torus
from exact_boolean_noise_semigroup import (
    P,
    biased_fourier_moments,
    mask_probability,
    popcount,
    truth_tables,
)
from integer_period_torus import IntegerTorusGeometry, gaussian_integer_torus


Q = Fraction
VECTORS = {
    "topology_pair": ("primal_cross", "matching_complement_cross"),
    "symmetry_pair": ("orientation_difference", "matching_odd_cross"),
}


def text(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def matrix_text(matrix: Sequence[Sequence[Q]]) -> list[list[str]]:
    return [[text(value) for value in row] for row in matrix]


def zero_matrix(width: int) -> list[list[Q]]:
    return [[Q(0) for _ in range(width)] for _ in range(width)]


def add_matrix(left, right):
    return [
        [left[i][j] + right[i][j] for j in range(len(left))]
        for i in range(len(left))
    ]


def subtract_matrix(left, right):
    return [
        [left[i][j] - right[i][j] for j in range(len(left))]
        for i in range(len(left))
    ]


def outer(left: Sequence[Q], right: Sequence[Q]) -> list[list[Q]]:
    return [[first * second for second in right] for first in left]


def scale_matrix(matrix, scalar: Q):
    return [[value * scalar for value in row] for row in matrix]


def state_key(mask: int, blocks: Sequence[Sequence[int]]) -> tuple[int, ...]:
    return tuple(sum(bool(mask & (1 << vertex)) for vertex in block) for block in blocks)


def group_by_partition(n: int, blocks: Sequence[Sequence[int]]) -> dict[tuple[int, ...], list[int]]:
    groups = defaultdict(list)
    for mask in range(1 << n):
        groups[state_key(mask, blocks)].append(mask)
    return dict(groups)


def validate_hierarchy(n: int, hierarchy: Sequence[dict]) -> None:
    universe = set(range(n))
    previous = [list(range(n))]
    for level in hierarchy:
        blocks = level["blocks"]
        flattened = [vertex for block in blocks for vertex in block]
        if set(flattened) != universe or len(flattened) != n:
            raise ValueError(f"{level['id']}: blocks do not partition the sites")
        for block in blocks:
            if not any(set(block).issubset(parent) for parent in map(set, previous)):
                raise ValueError(f"{level['id']}: partition is not nested")
        previous = blocks


def hierarchy_for(geometry: IntegerTorusGeometry) -> list[dict]:
    n = geometry.n
    if n == 5:
        return [
            {"id": "F0_total_occupation", "blocks": [list(range(5))]},
            {"id": "F1_marked_origin_vs_rest", "blocks": [[0], [1, 2, 3, 4]]},
            {"id": "F2_full_configuration", "blocks": [[site] for site in range(5)]},
        ]
    if n == 10:
        return [
            {"id": "F0_total_occupation", "blocks": [list(range(10))]},
            {
                "id": "F1_norm2_parent_fiber_counts",
                "blocks": [[site, site + 5] for site in range(5)],
            },
            {"id": "F2_full_configuration", "blocks": [[site] for site in range(10)]},
        ]
    raise ValueError("oracle is frozen to N=5 and N=10")


def conditional_means(groups, values: Sequence[Sequence[Q]]) -> dict[tuple[int, ...], list[Q]]:
    width = len(values[0])
    return {
        key: [sum(values[mask][column] for mask in masks) / len(masks) for column in range(width)]
        for key, masks in groups.items()
    }


def projection_by_mask(n: int, blocks, means) -> list[list[Q]]:
    return [means[state_key(mask, blocks)] for mask in range(1 << n)]


def expectation_outer(
    first: Sequence[Sequence[Q]],
    second: Sequence[Sequence[Q]],
    probabilities: Sequence[Q],
) -> list[list[Q]]:
    width = len(first[0])
    result = zero_matrix(width)
    for mask, probability in enumerate(probabilities):
        result = add_matrix(result, scale_matrix(outer(first[mask], second[mask]), probability))
    return result


def direct_conditional_replica_C(groups, values, probabilities) -> list[list[Q]]:
    """Enumerate two replicas conditional on one declared partition state."""

    width = len(values[0])
    result = zero_matrix(width)
    for masks in groups.values():
        group_probability = sum(probabilities[mask] for mask in masks)
        coefficient = group_probability / (len(masks) * len(masks))
        for first in masks:
            for second in masks:
                result = add_matrix(
                    result,
                    scale_matrix(outer(values[first], values[second]), coefficient),
                )
    return result


def four_replica_increment(
    n: int,
    coarse_blocks,
    fine_blocks,
    fine_means,
    probabilities,
) -> list[list[Q]]:
    """Collapsed exact four-replica sum over two independent refinements."""

    width = len(next(iter(fine_means.values())))
    result = zero_matrix(width)
    coarse_groups = group_by_partition(n, coarse_blocks)
    fine_groups = group_by_partition(n, fine_blocks)
    fine_probability = {
        key: sum(probabilities[mask] for mask in masks) for key, masks in fine_groups.items()
    }
    for coarse_key, masks in coarse_groups.items():
        coarse_probability = sum(probabilities[mask] for mask in masks)
        refinements = sorted({state_key(mask, fine_blocks) for mask in masks})
        for first_key in refinements:
            first_weight = fine_probability[first_key] / coarse_probability
            for second_key in refinements:
                second_weight = fine_probability[second_key] / coarse_probability
                difference = [
                    fine_means[first_key][column] - fine_means[second_key][column]
                    for column in range(width)
                ]
                result = add_matrix(
                    result,
                    scale_matrix(
                        outer(difference, difference),
                        coarse_probability * first_weight * second_weight / 2,
                    ),
                )
    return result


def rank_psd_payload(matrix: Sequence[Sequence[Q]]) -> dict:
    if len(matrix) != 2:
        raise ValueError("rank payload is frozen to two-coordinate vectors")
    a, b = matrix[0]
    c, d = matrix[1]
    if b != c:
        raise ArithmeticError("increment covariance is not symmetric")
    determinant = a * d - b * b
    if a < 0 or d < 0 or determinant < 0:
        raise ArithmeticError("increment covariance is not positive semidefinite")
    rank = 2 if determinant > 0 else (1 if a != 0 or b != 0 or d != 0 else 0)
    angle = None
    if rank:
        angle = degrees(0.5 * atan2(float(2 * b), float(a - d)))
    return {
        "rank": rank,
        "trace": text(a + d),
        "determinant": text(determinant),
        "psd_exact": True,
        "leading_axis_angle_degrees": angle,
        "offdiagonal_sign": "positive" if b > 0 else ("negative" if b < 0 else "zero"),
    }


def krawtchouk_value(n: int, degree: int, occupied: int, p: Q) -> Q:
    q = 1 - p
    return sum(
        Q(comb(occupied, chosen) * comb(n - occupied, degree - chosen))
        * q**chosen
        * (-p) ** (degree - chosen)
        for chosen in range(degree + 1)
        if chosen <= occupied and degree - chosen <= n - occupied
    )


def krawtchouk_control(
    names: Sequence[str],
    tables: dict[str, list[int]],
    centered: Sequence[Sequence[Q]],
    probabilities: Sequence[Q],
    moments: dict[str, list[Q]],
    C0,
    p: Q,
) -> dict:
    n = len(probabilities).bit_length() - 1
    pq = p * (1 - p)
    coefficients = []
    reconstruction = [[Q(0) for _ in names] for _ in range(n + 1)]
    covariance = zero_matrix(len(names))
    for degree in range(n + 1):
        norm = Q(comb(n, degree)) * pq**degree
        radial = []
        fourier = []
        for column, name in enumerate(names):
            numerator = sum(
                probabilities[mask]
                * centered[mask][column]
                * krawtchouk_value(n, degree, popcount(mask), p)
                for mask in range(1 << n)
            )
            coefficient = numerator / norm
            radial.append(coefficient)
            subset_sum = sum(
                moments[name][subset]
                for subset in range(1, 1 << n)
                if popcount(subset) == degree
            )
            if degree == 0:
                subset_sum = Q(0)
            fourier_coefficient = subset_sum / norm
            if coefficient != fourier_coefficient:
                raise ArithmeticError("Krawtchouk/Fourier radial coefficient mismatch")
            fourier.append(fourier_coefficient)
            for occupied in range(n + 1):
                reconstruction[occupied][column] += coefficient * krawtchouk_value(
                    n, degree, occupied, p
                )
        covariance = add_matrix(covariance, scale_matrix(outer(radial, radial), norm))
        coefficients.append(
            {
                "degree": degree,
                "norm": text(norm),
                "conditional_projection_coefficients": [text(value) for value in radial],
                "boolean_subset_average_coefficients": [text(value) for value in fourier],
            }
        )
    if covariance != C0:
        raise ArithmeticError("Krawtchouk covariance does not reproduce C0")
    total_groups = group_by_partition(n, [list(range(n))])
    radial_means = conditional_means(total_groups, centered)
    for occupied in range(n + 1):
        if reconstruction[occupied] != radial_means[(occupied,)]:
            raise ArithmeticError("Krawtchouk reconstruction of m0 failed")
    return {
        "basis": "K_r(K)=sum_j binom(K,j)binom(N-K,r-j)(1-p)^j(-p)^(r-j)",
        "coefficients": coefficients,
        "m0_reconstruction_by_occupation": {
            str(occupied): [text(value) for value in reconstruction[occupied]]
            for occupied in range(n + 1)
        },
        "C0_from_krawtchouk": matrix_text(covariance),
        "matches_PR245_p_biased_subset_moments": True,
    }


def analyze_vector(
    geometry: IntegerTorusGeometry,
    names: Sequence[str],
    tables: dict[str, list[int]],
    moments,
    p: Q,
) -> dict:
    n = geometry.n
    probabilities = [mask_probability(mask, n, p) for mask in range(1 << n)]
    means = [
        sum(probabilities[mask] * tables[name][mask] for mask in range(1 << n))
        for name in names
    ]
    centered = [
        [Q(tables[name][mask]) - means[column] for column, name in enumerate(names)]
        for mask in range(1 << n)
    ]
    hierarchy = hierarchy_for(geometry)
    validate_hierarchy(n, hierarchy)
    projections = []
    level_means = []
    groups_by_level = []
    for level in hierarchy:
        groups = group_by_partition(n, level["blocks"])
        conditional = conditional_means(groups, centered)
        groups_by_level.append(groups)
        level_means.append(conditional)
        projections.append(projection_by_mask(n, level["blocks"], conditional))

    zero_projection = [[Q(0) for _ in names] for _ in range(1 << n)]
    previous_projection = zero_projection
    previous_C = zero_matrix(len(names))
    increments = []
    gammas = []
    for level_index, level in enumerate(hierarchy):
        projection = projections[level_index]
        D = [
            [projection[mask][column] - previous_projection[mask][column] for column in range(len(names))]
            for mask in range(1 << n)
        ]
        C = expectation_outer(projection, projection, probabilities)
        gamma = expectation_outer(D, D, probabilities)
        if gamma != subtract_matrix(C, previous_C):
            raise ArithmeticError("Gamma != C_j-C_(j-1)")
        replica_C = direct_conditional_replica_C(
            groups_by_level[level_index], centered, probabilities
        )
        if replica_C != C:
            raise ArithmeticError("conditional-replica C identity failed")
        coarse_blocks = [] if level_index == 0 else hierarchy[level_index - 1]["blocks"]
        four_replica = four_replica_increment(
            n,
            coarse_blocks,
            level["blocks"],
            level_means[level_index],
            probabilities,
        )
        if four_replica != gamma:
            raise ArithmeticError("four-replica increment identity failed")

        histogram = defaultdict(lambda: {"configurations": 0, "probability": Q(0)})
        for mask, value in enumerate(D):
            key = tuple(value)
            histogram[key]["configurations"] += 1
            histogram[key]["probability"] += probabilities[mask]
        state_rows = {}
        for key, masks in groups_by_level[level_index].items():
            representative = masks[0]
            state_rows[",".join(map(str, key))] = {
                "multiplicity": len(masks),
                "probability": text(sum(probabilities[mask] for mask in masks)),
                "m_j": [text(value) for value in level_means[level_index][key]],
                "D_j": [text(value) for value in D[representative]],
            }
        increments.append(D)
        gammas.append(gamma)
        level["state_count"] = len(groups_by_level[level_index])
        level["C_j"] = matrix_text(C)
        level["C_j_conditional_replica"] = matrix_text(replica_C)
        level["Gamma_j"] = matrix_text(gamma)
        level["Gamma_j_four_replica"] = matrix_text(four_replica)
        level["rank_psd"] = rank_psd_payload(gamma)
        level["conditional_state_table"] = state_rows
        level["increment_histogram"] = [
            {
                "D_j": [text(value) for value in key],
                "configurations": row["configurations"],
                "probability": text(row["probability"]),
            }
            for key, row in sorted(histogram.items())
        ]
        previous_projection = projection
        previous_C = C

    orthogonality = {}
    for first, second in combinations(range(len(increments)), 2):
        cross = expectation_outer(increments[first], increments[second], probabilities)
        if cross != zero_matrix(len(names)):
            raise ArithmeticError("martingale increments are not orthogonal")
        orthogonality[f"D{first}_x_D{second}"] = matrix_text(cross)
    total_covariance = expectation_outer(centered, centered, probabilities)
    telescope = zero_matrix(len(names))
    for gamma in gammas:
        telescope = add_matrix(telescope, gamma)
    if telescope != total_covariance:
        raise ArithmeticError("Gamma telescope failed")

    total_trace = sum(total_covariance[index][index] for index in range(len(names)))
    for level, gamma in zip(hierarchy, gammas):
        gamma_trace = sum(gamma[index][index] for index in range(len(names)))
        C = [[Q(value) for value in row] for row in level["C_j"]]
        cumulative_trace = sum(C[index][index] for index in range(len(names)))
        level["trace_fraction_of_CovY"] = text(gamma_trace / total_trace)
        level["cumulative_predictable_trace_fraction"] = text(
            cumulative_trace / total_trace
        )

    krawtchouk = krawtchouk_control(
        names, tables, centered, probabilities, moments, gammas[0], p
    )
    angles = [level["rank_psd"]["leading_axis_angle_degrees"] for level in hierarchy]
    rotations = []
    previous_angle = None
    previous_id = None
    for level, angle in zip(hierarchy, angles):
        if angle is None:
            continue
        if previous_angle is not None:
            raw = abs(angle - previous_angle) % 180.0
            rotations.append(
                {
                    "from": previous_id,
                    "to": level["id"],
                    "principal_axis_rotation_degrees": min(raw, 180.0 - raw),
                }
            )
        previous_angle = angle
        previous_id = level["id"]
    return {
        "observable_order": list(names),
        "means": [text(value) for value in means],
        "levels": hierarchy,
        "martingale_orthogonality": orthogonality,
        "covariance_Y": matrix_text(total_covariance),
        "sum_Gamma": matrix_text(telescope),
        "telescoping_exact": True,
        "krawtchouk_F0_control": krawtchouk,
        "leading_axis_angles_degrees": angles,
        "consecutive_principal_axis_rotations": rotations,
    }


def analyze_geometry(geometry: IntegerTorusGeometry, p: Q) -> dict:
    tables = truth_tables(geometry)
    moments, _means = biased_fourier_moments(tables, geometry.n, p)
    vectors = {
        label: analyze_vector(geometry, names, tables, moments, p)
        for label, names in VECTORS.items()
    }
    return {
        "geometry": geometry.name,
        "N": geometry.n,
        "p": text(p),
        "configurations": 1 << geometry.n,
        "vectors": vectors,
    }


def render() -> dict:
    geometries = [gaussian_integer_torus(2, 1), c4_self_matching_torus(3, 1)]
    payload = {
        "schema": "matching-one.exact-conditional-replica-filtration.v1",
        "issue": 256,
        "status": "exact_finite_volume_martingale_oracle",
        "p": text(P),
        "geometries": [analyze_geometry(geometry, P) for geometry in geometries],
        "exact_identities": [
            "m_j=E[Y|F_j] and D_j=m_j-m_(j-1)",
            "C_j=E[m_j m_j^T]=E[Y(X)Y(X')^T|conditional replicas]",
            "Gamma_j=E[D_jD_j^T]=C_j-C_(j-1) is PSD",
            "Cov(Y)=sum_j Gamma_j",
            "Gamma_j equals the symmetric four-replica refinement identity",
            "F0=sigma(K) equals the Krawtchouk/radial projection from PR245 moments",
        ],
        "mechanism_conjecture": {
            "status": "new_risky_conjecture_from_exact_N10_structure",
            "statement": (
                "The odd-shell anomaly is a rotating two-dimensional primal/matching "
                "topology plane: radial, norm-2-fiber, and singleton increments retain "
                "rank two while their covariance axis changes sign/angle. At large N, "
                "a stable rotating rank-2 plane should replace any constant scalar shell law."
            ),
            "falsifier": (
                "At larger N the topology-pair Gamma increments become rank one in one "
                "fixed basis, or their leading subspace/axis has no cross-size stability."
            ),
            "representation_warning": (
                "The orientation_difference/matching_odd pair is exactly cross-orthogonal "
                "on these controls; simultaneous rank two there is a direct sum of irreps, "
                "not coordinate rotation. Rotation claims must be made within a symmetry block."
            ),
        },
    }
    return payload


def main(argv: Iterable[str] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(list(argv) if argv is not None else None)
    payload = render()
    output = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(output, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
