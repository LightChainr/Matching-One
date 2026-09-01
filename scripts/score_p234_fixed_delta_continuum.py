#!/usr/bin/env python3
"""Score the fixed-delta continuum collapse of the P234 parent log pair.

The Phase-A batch archive already contains the connection probability between
the two endpoints of each bilocal insertion.  At one fixed physical delta,
this probability is proportional to pi_a**2, with an L-independent factor.
It therefore supplies the missing *relative* field normalization without a
new Monte Carlo run.  It does not fix a universal amplitude or compare two
different delta values.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Sequence

import mpmath


Vector = Sequence[float]
Matrix = list[list[float]]


def covariance_of_mean(rows: Sequence[Vector]) -> Matrix:
    if len(rows) < 2:
        raise ValueError("at least two batches are required")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("batch vectors must have equal width")
    means = [sum(float(row[j]) for row in rows) / len(rows) for j in range(width)]
    return [
        [
            sum(
                (float(row[i]) - means[i]) * (float(row[j]) - means[j])
                for row in rows
            )
            / (len(rows) * (len(rows) - 1))
            for j in range(width)
        ]
        for i in range(width)
    ]


def transform(point: Vector, alpha: float, beta_numerator: float) -> list[float]:
    ll, ld, dd, connection = (float(value) for value in point)
    if connection <= 0.0:
        raise ValueError("bilocal endpoint connection probability must be positive")
    return [
        alpha * alpha * ll,
        alpha * beta_numerator * ld / connection,
        beta_numerator * beta_numerator * dd / (connection * connection),
    ]


def transform_jacobian(point: Vector, alpha: float, beta_numerator: float) -> Matrix:
    _, ld, dd, connection = (float(value) for value in point)
    return [
        [alpha * alpha, 0.0, 0.0, 0.0],
        [
            0.0,
            alpha * beta_numerator / connection,
            0.0,
            -alpha * beta_numerator * ld / connection**2,
        ],
        [
            0.0,
            0.0,
            beta_numerator**2 / connection**2,
            -2.0 * beta_numerator**2 * dd / connection**3,
        ],
    ]


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
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError("linear system must be square")
    augmented = [list(float(x) for x in matrix[i]) + [float(vector[i])] for i in range(n)]
    scale = max(abs(value) for row in matrix for value in row)
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= 1e-12 * max(scale, 1e-300):
            raise ValueError("covariance is numerically singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                augmented[row][j] - factor * augmented[column][j]
                for j in range(n + 1)
            ]
    return [augmented[i][-1] for i in range(n)]


def chi_square_sf(value: float, degrees: int) -> float:
    return float(mpmath.gammainc(degrees / 2, value / 2, mpmath.inf, regularized=True))


def block_diagonal(blocks: Sequence[Matrix]) -> Matrix:
    width = sum(len(block) for block in blocks)
    answer = [[0.0] * width for _ in range(width)]
    offset = 0
    for block in blocks:
        for i, row in enumerate(block):
            for j, value in enumerate(row):
                answer[offset + i][offset + j] = float(value)
        offset += len(block)
    return answer


def quadratic_form(vector: Vector, covariance: Matrix) -> float:
    solved = solve(covariance, vector)
    return sum(float(a) * b for a, b in zip(vector, solved))


def _gls(points: Sequence[Vector], covariances: Sequence[Matrix], design_columns: Sequence[Vector], names: Sequence[str], model: str) -> dict:
    y = [float(value) for point in points for value in point]
    covariance = block_diagonal(covariances)
    inverse_columns = []
    for column in design_columns:
        inverse_columns.append(solve(covariance, column))
    inverse_y = solve(covariance, y)
    normal = [
        [sum(design_columns[i][k] * inverse_columns[j][k] for k in range(len(y))) for j in range(len(design_columns))]
        for i in range(len(design_columns))
    ]
    rhs = [sum(column[k] * inverse_y[k] for k in range(len(y))) for column in design_columns]
    coefficients = solve(normal, rhs)
    fitted = [
        sum(coefficients[j] * design_columns[j][i] for j in range(len(coefficients)))
        for i in range(len(y))
    ]
    residual = [value - target for value, target in zip(y, fitted)]
    chi_square = quadratic_form(residual, covariance)
    degrees = len(y) - len(coefficients)
    coefficient_covariance = [
        solve(normal, [1.0 if i == j else 0.0 for i in range(len(coefficients))])
        for j in range(len(coefficients))
    ]
    coefficient_covariance = [
        [coefficient_covariance[j][i] for j in range(len(coefficients))]
        for i in range(len(coefficients))
    ]
    return {
        "model": model,
        "coefficient_order": list(names),
        "coefficients": coefficients,
        "coefficient_covariance": coefficient_covariance,
        "chi_square": chi_square,
        "degrees_of_freedom": degrees,
        "chi_square_survival": chi_square_sf(chi_square, degrees),
    }


def gls_parent_pair(points: Sequence[Vector], covariances: Sequence[Matrix]) -> dict:
    """Fit [zero bottom-bottom, constant mixed, constant top-top]."""
    width = 3 * len(points)
    design = [
        [1.0 if i % 3 == coordinate else 0.0 for i in range(width)]
        for coordinate in (1, 2)
    ]
    return _gls(
        points,
        covariances,
        design,
        ["constant_LD", "constant_DD"],
        "normalized_LL=0; normalized_LD and normalized_DD constant at every finite L",
    )


def gls_linear_continuum(sizes: Sequence[float], points: Sequence[Vector], covariances: Sequence[Matrix]) -> dict:
    """Fit the continuum target with one leading analytic 1/L correction."""
    if len(sizes) != len(points):
        raise ValueError("sizes and points must align")
    columns = [[0.0] * (3 * len(points)) for _ in range(5)]
    for block, size in enumerate(sizes):
        inverse_size = 1.0 / float(size)
        columns[0][3 * block] = inverse_size
        columns[1][3 * block + 1] = 1.0
        columns[2][3 * block + 1] = inverse_size
        columns[3][3 * block + 2] = 1.0
        columns[4][3 * block + 2] = inverse_size
    return _gls(
        points,
        covariances,
        columns,
        ["LL_1_over_L", "LD_continuum", "LD_1_over_L", "DD_continuum", "DD_1_over_L"],
        "normalized_LL=c_LL/L; normalized_LD=C_LD+c_LD/L; normalized_DD=C_DD+c_DD/L",
    )


def mesh_aware_log_partner_diagnostic(rows: Sequence[dict]) -> dict:
    """Resolve signed radius rounding after the frozen L192 reveal.

    This is deliberately labelled post-reveal.  It uses the natural realized
    cutoff in the bilocal prefactor, then gives each matrix coordinate a
    leading signed delta error in addition to the analytic 1/L correction.
    """
    points = []
    covariances = []
    for row in rows:
        ratio = (row["realized_delta"] / row["declared_delta"]) ** (-25.0 / 24.0)
        scaling = [1.0, ratio, ratio * ratio]
        points.append([value * scale for value, scale in zip(row["normalized_point"], scaling)])
        covariance = row["normalized_covariance_delta_method"]
        covariances.append(
            [
                [covariance[i][j] * scaling[i] * scaling[j] for j in range(3)]
                for i in range(3)
            ]
        )
    width = 3 * len(rows)
    columns = [[0.0] * width for _ in range(8)]
    for block, row in enumerate(rows):
        inverse_size = 1.0 / float(row["L"])
        delta_error = row["realized_delta"] - row["declared_delta"]
        columns[0][3 * block] = inverse_size
        columns[1][3 * block] = delta_error
        columns[2][3 * block + 1] = 1.0
        columns[3][3 * block + 1] = inverse_size
        columns[4][3 * block + 1] = delta_error
        columns[5][3 * block + 2] = 1.0
        columns[6][3 * block + 2] = inverse_size
        columns[7][3 * block + 2] = delta_error
    score = _gls(
        points,
        covariances,
        columns,
        [
            "LL_1_over_L",
            "LL_delta_error",
            "LD_continuum",
            "LD_1_over_L",
            "LD_delta_error",
            "DD_continuum",
            "DD_1_over_L",
            "DD_delta_error",
        ],
        "realized-cutoff fields with 1/L and signed (delta_realized-delta_declared) corrections",
    )
    mixed = score["coefficients"][2]
    dd_derivative = score["coefficients"][7]
    declared_delta = float(rows[0]["declared_delta"])
    kappa = -declared_delta * dd_derivative / (2.0 * mixed)
    gradient_mixed = declared_delta * dd_derivative / (2.0 * mixed**2)
    gradient_derivative = -declared_delta / (2.0 * mixed)
    covariance = score["coefficient_covariance"]
    kappa_variance = (
        gradient_mixed**2 * covariance[2][2]
        + gradient_derivative**2 * covariance[7][7]
        + 2.0 * gradient_mixed * gradient_derivative * covariance[2][7]
    )
    return {
        "status": "post_reveal_mechanism_diagnostic_added_after_the_frozen_L192_score",
        "natural_realized_cutoff_points": points,
        "score": score,
        "cutoff_shear_proxy": {
            "definition": "kappa_proxy=-delta*(d DD/d delta)/(2*LD) from B_delta=hat_phi-kappa*log(2delta)*phi",
            "estimate": kappa,
            "standard_error_delta_method": math.sqrt(max(0.0, kappa_variance)),
        },
    }


def read_size(json_path: Path) -> dict:
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    batch_path = json_path.with_suffix(".batches.csv")
    with batch_path.open(newline="", encoding="utf-8") as handle:
        batch_rows = list(csv.DictReader(handle))
    estimates = [[float(value) for value in row] for row in payload["block_estimates"]]
    if len(estimates) != len(batch_rows):
        raise ValueError(f"batch count mismatch for {json_path}")
    L = int(payload["geometry"]["L"])
    volume = L * L
    joint_rows = []
    for estimate, raw in zip(estimates, batch_rows):
        samples = int(raw["samples"])
        connection = (int(raw["sum_D1"]) + int(raw["sum_D2"])) / (2.0 * samples * volume)
        joint_rows.append(estimate + [connection])
    point = [sum(row[j] for row in joint_rows) / len(joint_rows) for j in range(4)]
    raw_covariance = covariance_of_mean(joint_rows)
    alpha = L ** 1.25 / math.log(L)
    declared_delta = float(payload["geometry"]["delta"])
    realized_delta = float(payload["geometry"]["realized_delta"])
    beta_numerator = (2.0 * declared_delta) ** (-25.0 / 24.0) * (2.0 * realized_delta) ** (-5.0 / 24.0)
    normalized = transform(point, alpha, beta_numerator)
    normalized_covariance = propagate(
        transform_jacobian(point, alpha, beta_numerator), raw_covariance
    )
    return {
        "L": L,
        "declared_delta": payload["geometry"]["delta"],
        "realized_delta": payload["geometry"]["realized_delta"],
        "input_json": str(json_path),
        "input_batches": str(batch_path),
        "batches": len(batch_rows),
        "connection_probability_proxy_for_pi_a_squared": point[3],
        "connection_probability_standard_error": math.sqrt(raw_covariance[3][3]),
        "local_normalization_alpha_L_5_over_4_over_log_L": alpha,
        "bilocal_normalization_numerator": beta_numerator,
        "normalized_order": [
            "alpha^2_LL",
            "alpha*beta_numerator*LD/p_connection",
            "beta_numerator^2*DD/p_connection^2",
        ],
        "normalized_point": normalized,
        "normalized_standard_error": [math.sqrt(max(0.0, normalized_covariance[i][i])) for i in range(3)],
        "normalized_covariance_delta_method": normalized_covariance,
    }


def render(paths: Sequence[Path]) -> dict:
    sizes = sorted((read_size(path) for path in paths), key=lambda row: row["L"])
    if len(sizes) < 2:
        raise ValueError("at least two sizes are required")
    deltas = {row["declared_delta"] for row in sizes}
    if len(deltas) != 1:
        raise ValueError("the connection proxy is only comparable at one fixed declared delta")
    finite_score = gls_parent_pair(
        [row["normalized_point"] for row in sizes],
        [row["normalized_covariance_delta_method"] for row in sizes],
    )
    continuum_score = gls_linear_continuum(
        [row["L"] for row in sizes],
        [row["normalized_point"] for row in sizes],
        [row["normalized_covariance_delta_method"] for row in sizes],
    )
    return {
        "schema": "matching-one.p234-fixed-delta-continuum.v1",
        "issue": 234,
        "status": "production score",
        "sizes": sizes,
        "joint_parent_pair_continuum_linear_in_1_over_L": continuum_score,
        "finite_L_constant_diagnostic": finite_score,
        "post_reveal_mesh_aware_log_partner_diagnostic": mesh_aware_log_partner_diagnostic(sizes),
        "derivation": {
            "paper_normalizations": "phi_a=a^(-5/4)|log a|^(-1)E_a; eta_a^delta=pi_a^(-2)E_a^delta",
            "connection_proxy": "pi_a^(-2)=constant*(2*delta_realized)^(-5/24)/p_connection*(1+o(1)).",
            "bilocal_factor": "beta_numerator=(2*delta_declared)^(-25/24)*(2*delta_realized)^(-5/24)",
            "scored_vector": "[alpha^2 LL, alpha*beta_numerator*LD/p_connection, beta_numerator^2*DD/p_connection^2]",
        },
        "scope": [
            "The endpoint-connection proxy and LL/LD/DD come from the same batch stream; their covariance is retained.",
            "The unknown connection-to-pi amplitude is L-independent only at one fixed declared delta, so it is absorbed into the two fitted constants.",
            "This score tests the a->0 continuum collapse at fixed physical geometry. It does not estimate kappa, a universal amplitude, or a Jordan logarithmic coupling.",
            "A coordinate-dilation Jordan shear is not a log(L) drift in a correctly renormalized fixed-geometry ultraviolet limit.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = render(args.inputs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
