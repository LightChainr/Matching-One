#!/usr/bin/env python3
"""Score the projective Jordan pencil across fixed-delta P234 sizes.

For a canonical logarithmic two-point form, normalizing each Gram matrix by
its mixed entry turns the cross-size pencil into a lower unipotent matrix.
The lower entry is free; the other three entries form a normalization-free
three-degree-of-freedom gate.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Callable, Sequence

import mpmath


Vector = Sequence[float]
Matrix = list[list[float]]


def gram(vector: Vector, relative_field_scale: float = 1.0) -> Matrix:
    ll, ld, dd = (float(value) for value in vector)
    if ld == 0.0:
        raise ValueError("LD must be nonzero for projective normalization")
    if relative_field_scale == 0.0:
        raise ValueError("relative field scale must be nonzero")
    return [
        [relative_field_scale * ll / ld, 1.0],
        [1.0, dd / (relative_field_scale * ld)],
    ]


def determinant(value: Matrix) -> float:
    return value[0][0] * value[1][1] - value[0][1] * value[1][0]


def inverse2(value: Matrix) -> Matrix:
    det = determinant(value)
    if det == 0.0:
        raise ValueError("projective source Gram matrix is singular")
    return [
        [value[1][1] / det, -value[0][1] / det],
        [-value[1][0] / det, value[0][0] / det],
    ]


def multiply2(left: Matrix, right: Matrix) -> Matrix:
    return [
        [sum(left[i][k] * right[k][j] for k in range(2)) for j in range(2)]
        for i in range(2)
    ]


def frobenius(value: Matrix) -> float:
    return math.sqrt(sum(item * item for row in value for item in row))


def pencil(
    first: Vector,
    second: Vector,
    first_relative_scale: float = 1.0,
    second_relative_scale: float = 1.0,
) -> Matrix:
    return multiply2(
        gram(second, second_relative_scale),
        inverse2(gram(first, first_relative_scale)),
    )


def gate_vector(
    first: Vector,
    second: Vector,
    first_relative_scale: float = 1.0,
    second_relative_scale: float = 1.0,
) -> list[float]:
    transfer = pencil(first, second, first_relative_scale, second_relative_scale)
    return [transfer[0][0] - 1.0, transfer[0][1], transfer[1][1] - 1.0]


def flow_rates(
    vectors: Sequence[Vector],
    sizes: Sequence[float],
    relative_scales: Sequence[float],
) -> list[float]:
    if not (len(vectors) == len(sizes) == len(relative_scales)):
        raise ValueError("vectors, sizes, and relative scales must align")
    answer = []
    for index in range(len(vectors) - 1):
        transfer = pencil(
            vectors[index],
            vectors[index + 1],
            relative_scales[index],
            relative_scales[index + 1],
        )
        answer.append(transfer[1][0] / math.log(sizes[index + 1] / sizes[index]))
    return answer


def _block_covariance(first: Matrix, second: Matrix) -> Matrix:
    answer = [[0.0] * 6 for _ in range(6)]
    for i in range(3):
        for j in range(3):
            answer[i][j] = float(first[i][j])
            answer[i + 3][j + 3] = float(second[i][j])
    return answer


def _many_block_covariance(blocks: Sequence[Matrix]) -> Matrix:
    size = 3 * len(blocks)
    answer = [[0.0] * size for _ in range(size)]
    for block, covariance in enumerate(blocks):
        for i in range(3):
            for j in range(3):
                answer[3 * block + i][3 * block + j] = float(covariance[i][j])
    return answer


def numerical_jacobian(function: Callable[[Vector], list[float]], point: Vector) -> Matrix:
    base = list(float(value) for value in point)
    rows = len(function(base))
    answer = [[0.0] * len(base) for _ in range(rows)]
    for column, value in enumerate(base):
        step = 1e-6 * max(abs(value), 1e-8)
        plus = base[:]
        minus = base[:]
        plus[column] += step
        minus[column] -= step
        f_plus = function(plus)
        f_minus = function(minus)
        for row in range(rows):
            answer[row][column] = (f_plus[row] - f_minus[row]) / (2.0 * step)
    return answer


def propagate(jacobian: Matrix, covariance: Matrix) -> Matrix:
    return [
        [
            sum(
                jacobian[i][a] * covariance[a][b] * jacobian[j][b]
                for a in range(len(covariance))
                for b in range(len(covariance))
            )
            for j in range(len(jacobian))
        ]
        for i in range(len(jacobian))
    ]


def solve(matrix: Matrix, vector: Vector) -> list[float]:
    n = len(vector)
    augmented = [list(matrix[i]) + [float(vector[i])] for i in range(n)]
    scale = max(abs(value) for row in matrix for value in row)
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= 1e-12 * max(1.0, scale):
            raise ValueError("gate covariance is numerically singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                augmented[row][item] - factor * augmented[column][item]
                for item in range(n + 1)
            ]
    return [augmented[i][-1] for i in range(n)]


def chi_square_sf(value: float, degrees: int) -> float:
    return float(mpmath.gammainc(degrees / 2, value / 2, mpmath.inf, regularized=True))


def pair_score(
    first: dict,
    second: dict,
    first_relative_scale: float,
    second_relative_scale: float,
) -> dict:
    first_point = [float(value) for value in first["estimate"]]
    second_point = [float(value) for value in second["estimate"]]
    joined = first_point + second_point
    covariance = _block_covariance(
        first["covariance_of_mean"], second["covariance_of_mean"]
    )
    function = lambda values: gate_vector(
        values[:3], values[3:], first_relative_scale, second_relative_scale
    )
    gate = function(joined)
    jacobian = numerical_jacobian(function, joined)
    gate_covariance = propagate(jacobian, covariance)
    try:
        solved = solve(gate_covariance, gate)
        chi_square = sum(gate[i] * solved[i] for i in range(3))
        survival = chi_square_sf(chi_square, 3)
    except ValueError:
        chi_square = None
        survival = None

    transfer = pencil(
        first_point, second_point, first_relative_scale, second_relative_scale
    )
    trace = transfer[0][0] + transfer[1][1]
    det = determinant(transfer)
    first_gram = gram(first_point, first_relative_scale)
    condition = frobenius(first_gram) * frobenius(inverse2(first_gram))
    first_geometry = first["geometry"]
    second_geometry = second["geometry"]
    return {
        "L_pair": [first_geometry["L"], second_geometry["L"]],
        "declared_delta": [first_geometry["delta"], second_geometry["delta"]],
        "realized_delta": [
            first_geometry["realized_delta"],
            second_geometry["realized_delta"],
        ],
        "relative_field_scales_L_over_D": [
            first_relative_scale,
            second_relative_scale,
        ],
        "projective_transfer_M2_times_M1_inverse": transfer,
        "source_projective_condition_number_frobenius": condition,
        "free_log_mixing_entry": transfer[1][0],
        "similarity_invariants": {
            "trace": trace,
            "determinant": det,
            "discriminant": trace * trace - 4.0 * det,
        },
        "lower_unipotent_gate": {
            "order": ["T00_minus_1", "T01", "T11_minus_1"],
            "point": gate,
            "covariance_delta_method": gate_covariance,
            "chi_square": chi_square,
            "degrees_of_freedom": 3,
            "chi_square_survival": survival,
            "information_warning": "With exact LL=0 at both sizes this gate is an algebraic identity for arbitrary DD drift; it is not an independent constant-flow test.",
        },
    }


def flow_score(payloads: Sequence[dict], scales: Sequence[float]) -> dict:
    sizes = [float(row["geometry"]["L"]) for row in payloads]
    vectors = [[float(value) for value in row["estimate"]] for row in payloads]
    joined = [value for vector in vectors for value in vector]
    covariance = _many_block_covariance(
        [row["covariance_of_mean"] for row in payloads]
    )

    def function(values: Vector) -> list[float]:
        split = [values[3 * i : 3 * i + 3] for i in range(len(payloads))]
        return flow_rates(split, sizes, scales)

    rates = function(joined)
    jacobian = numerical_jacobian(function, joined)
    rate_covariance = propagate(jacobian, covariance)
    ones = [1.0] * len(rates)
    try:
        inverse_y = solve(rate_covariance, rates)
        inverse_one = solve(rate_covariance, ones)
        fitted = sum(inverse_y) / sum(inverse_one)
        residual = [value - fitted for value in rates]
        inverse_residual = solve(rate_covariance, residual)
        chi_square = sum(a * b for a, b in zip(residual, inverse_residual))
        degrees = len(rates) - 1
        survival = chi_square_sf(chi_square, degrees) if degrees else None
        fitted_se = math.sqrt(1.0 / sum(inverse_one))
    except ValueError:
        fitted = fitted_se = chi_square = survival = None
        degrees = len(rates) - 1
    return {
        "adjacent_L_pairs": [[int(sizes[i]), int(sizes[i + 1])] for i in range(len(sizes) - 1)],
        "mixing_rate_T10_per_log_L": rates,
        "covariance_delta_method": rate_covariance,
        "constant_rate_GLS": {
            "estimate": fitted,
            "standard_error": fitted_se,
            "chi_square": chi_square,
            "degrees_of_freedom": degrees,
            "chi_square_survival": survival,
        },
    }


def render(inputs: Sequence[Path], relative_scales: Sequence[float]) -> dict:
    if len(inputs) != len(relative_scales):
        raise ValueError("one relative field scale is required for every input")
    rows = [
        (json.loads(path.read_text(encoding="utf-8")), float(scale))
        for path, scale in zip(inputs, relative_scales)
    ]
    rows.sort(key=lambda row: row[0]["geometry"]["L"])
    payloads = [row[0] for row in rows]
    scales = [row[1] for row in rows]
    if len(payloads) < 2:
        raise ValueError("at least two P234 size blocks are required")
    denominators = {row["geometry"]["delta_formula"] for row in payloads}
    if len(denominators) != 1:
        raise ValueError("all size blocks must use the same declared fixed delta")
    return {
        "schema": "matching-one.p234-projective-jordan-pencil.v1",
        "issue": 234,
        "input_sizes": [row["geometry"]["L"] for row in payloads],
        "relative_field_scales_L_over_D": scales,
        "pairs": [
            pair_score(left, right, left_scale, right_scale)
            for left, right, left_scale, right_scale in zip(
                payloads, payloads[1:], scales, scales[1:]
            )
        ],
        "projective_flow": flow_score(payloads, scales),
        "exact_target": {
            "canonical_gram": "M(t)/LD = [[0,1],[1,d0+k*t]]",
            "pencil": "M(t2) M(t1)^-1 = [[1,0],[k*(t2-t1),1]]",
            "gate": "T00=1, T01=0, T11=1; T10 is an unfitted log-mixing coordinate",
            "constant_generator": "T10/log(L2/L1) is the same on every adjacent interval",
        },
        "scope": [
            "Each supplied scale is A_L/A_D for physical fields A_L*L_raw and A_D*D_raw; LD normalization alone cannot remove a size-dependent relative field gauge.",
            "The raw Phase-A blocks do not yet determine A_D because pi_a is absent, so no production pencil score is valid until those scale ratios are supplied.",
            "The two-size lower-unipotent gate collapses to the per-size LL-null condition; three or more calibrated sizes are required for new flow information.",
            "Delta-method covariance assumes the independently seeded size blocks are independent.",
            "A near-singular source Gram matrix makes the pencil weakly identified and is reported by its condition number.",
            "Passing is a projective Jordan-shape result, not a universal kappa measurement.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument(
        "--relative-scales",
        nargs="+",
        required=True,
        type=float,
        help="A_L/A_D for each input, in the same order as the input paths",
    )
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = render(args.inputs, args.relative_scales)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
