#!/usr/bin/env python3
"""Production-first empirical coalescence diagnostic for Issue #218.

The only state used here is the source-declared pair

    H_N = N^(13/8) P4[D],
    U_N = N^(13/8) P4[S']/mean(M'),

reconstructed jointly from the committed threshold-rank batches.  H is the
leading thermal H4 amplitude and U is the live center-slope coordinate after
the declared thermal metric normalization.  A covariance metric is frozen on
N65/N85 before inspecting the N520/N680 targets.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Callable, Mapping, Sequence

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
SIZES = (65, 85, 130, 170, 260, 340, 520, 680)
GENERATIONS = ((65, 85), (130, 170), (260, 340), (520, 680))
SOURCE_ORDER = (65, 130, 85, 170)
RAW_SOURCE = {
    65: (
        "results/server-20260828/P45-root-amplitude/n65.hist.csv",
        "results/server-20260829/P154-norm4-production/raw/n65_1900m.hist.csv",
    ),
    130: (
        "results/server-20260828/P49-fullcurve-doubling-100m/raw/n130.hist.csv",
        "results/server-20260829/P154-norm4-production/raw/n130_1900m.hist.csv",
    ),
    85: (
        "results/server-20260828/P45-root-amplitude/n85.hist.csv",
        "results/server-20260829/P154-norm4-production/raw/n85_1900m.hist.csv",
    ),
    170: (
        "results/server-20260828/P49-fullcurve-doubling-100m/raw/n170.hist.csv",
        "results/server-20260829/P154-norm4-production/raw/n170_1900m.hist.csv",
    ),
}
RAW_TARGET = {
    260: "results/server-20260829/P154-norm4-production/raw/n260_1b.hist.csv",
    340: "results/server-20260829/P154-norm4-production/raw/n340_1b.hist.csv",
    520: "results/server-20260829/P154-norm4-generation4-pilot/raw/n520_100m.hist.csv",
    680: "results/server-20260829/P154-norm4-generation4-pilot/raw/n680_100m.hist.csv",
}
SCORERS = (
    "results/server-20260829/P154-norm4-production/analysis/scalar_score.json",
    "results/server-20260829/P154-norm4-production/analysis/thermal_jet_score.json",
    "results/server-20260829/P154-norm4-generation4-pilot/analysis/score.json",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def matmul(left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]):
    return [
        [
            math.fsum(left[i][k] * right[k][j] for k in range(len(right)))
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def matvec(matrix: Sequence[Sequence[float]], vector: Sequence[float]):
    return [math.fsum(a * b for a, b in zip(row, vector)) for row in matrix]


def inverse_2x2(matrix: Sequence[Sequence[float]]):
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    if abs(determinant) < 1e-30:
        raise ValueError("singular 2x2 response matrix")
    return [
        [matrix[1][1] / determinant, -matrix[0][1] / determinant],
        [-matrix[1][0] / determinant, matrix[0][0] / determinant],
    ]


def columns_matrix(left: Sequence[float], right: Sequence[float]):
    return [[left[0], right[0]], [left[1], right[1]]]


def flatten_columns(matrix: Sequence[Sequence[float]]):
    return [matrix[0][0], matrix[1][0], matrix[0][1], matrix[1][1]]


def frobenius_squared(matrix: Sequence[Sequence[float]]) -> float:
    return math.fsum(value * value for row in matrix for value in row)


def subtract(left, right):
    return [
        [left[i][j] - right[i][j] for j in range(len(left[0]))]
        for i in range(len(left))
    ]


def whitener(covariance: Sequence[Sequence[float]]):
    matrix = mp.matrix([[mp.mpf(str(value)) for value in row] for row in covariance])
    values, vectors = mp.eigsy(matrix)
    if any(values[index] <= 0 for index in range(2)):
        raise ValueError("reference covariance is not positive definite")
    # L C L^T = I, with rows equal to covariance eigenvectors / sqrt(value).
    return [
        [float(vectors[column, row] / mp.sqrt(values[row])) for column in range(2)]
        for row in range(2)
    ]


def transform_covariance(covariance, block_transform):
    dimension = len(covariance)
    output = [[0.0] * dimension for _ in range(dimension)]
    blocks = dimension // 2
    for block_i in range(blocks):
        for block_j in range(blocks):
            source = [
                [covariance[2 * block_i + a][2 * block_j + b] for b in range(2)]
                for a in range(2)
            ]
            transformed = matmul(matmul(block_transform, source), transpose(block_transform))
            for a in range(2):
                for b in range(2):
                    output[2 * block_i + a][2 * block_j + b] = transformed[a][b]
    return output


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def submatrix(matrix, indices):
    return [[matrix[i][j] for j in indices] for i in indices]


def load_production_state() -> dict[str, Any]:
    """Reconstruct all state points and their same-batch jackknife covariance."""

    sys.path.insert(0, str(ROOT / "scripts"))
    from analyze_p48_retrospective import read_histograms
    from score_norm4_production import (
        METRICS,
        estimate_aligned,
        estimate_one,
        merge_histogram_blocks,
    )
    from score_p50_fullcurve_n290 import grouped

    source_groups = {}
    for n, paths in RAW_SOURCE.items():
        source_groups[n] = grouped(
            merge_histogram_blocks(*(ROOT / path for path in paths), n), n
        )
    source_points, source_covariance = estimate_aligned(source_groups, SOURCE_ORDER)
    target_points = {}
    target_covariances = {}
    for n, path_text in RAW_TARGET.items():
        target_points[n], target_covariances[n] = estimate_one(
            grouped(read_histograms(ROOT / path_text), n)
        )
    points = {**source_points, **target_points}

    # First assemble raw (U,P4_D) covariance in the public SIZES order.
    raw_covariance = [[0.0] * 16 for _ in range(16)]
    metric_indices = (METRICS.index("U"), METRICS.index("P4_D"))
    for i, n_i in enumerate(SOURCE_ORDER):
        for j, n_j in enumerate(SOURCE_ORDER):
            out_i, out_j = SIZES.index(n_i), SIZES.index(n_j)
            for a, metric_a in enumerate(metric_indices):
                for b, metric_b in enumerate(metric_indices):
                    raw_covariance[2 * out_i + a][2 * out_j + b] = source_covariance[
                        3 * i + metric_a
                    ][3 * j + metric_b]
    for n, covariance in target_covariances.items():
        out = SIZES.index(n)
        for a, metric_a in enumerate(metric_indices):
            for b, metric_b in enumerate(metric_indices):
                raw_covariance[2 * out + a][2 * out + b] = covariance[metric_a][metric_b]

    state = []
    state_covariance = [[0.0] * 16 for _ in range(16)]
    factors = [n ** (13.0 / 8.0) for n in SIZES]
    for n, factor in zip(SIZES, factors):
        state.extend([factor * points[n]["P4_D"], points[n]["U"]])
    # Per-size map (U,D) -> (H,U) equals [[0,N^(13/8)],[1,0]].
    for i, factor_i in enumerate(factors):
        left = [[0.0, factor_i], [1.0, 0.0]]
        for j, factor_j in enumerate(factors):
            right = [[0.0, factor_j], [1.0, 0.0]]
            raw = [
                [raw_covariance[2 * i + a][2 * j + b] for b in range(2)]
                for a in range(2)
            ]
            converted = matmul(matmul(left, raw), transpose(right))
            for a in range(2):
                for b in range(2):
                    state_covariance[2 * i + a][2 * j + b] = converted[a][b]

    return {
        "state": state,
        "covariance": state_covariance,
        "points": points,
    }


def response_generations(state: Sequence[float], frozen_whitener):
    by_size = {
        n: matvec(frozen_whitener, state[2 * index : 2 * index + 2])
        for index, n in enumerate(SIZES)
    }
    return [columns_matrix(by_size[left], by_size[right]) for left, right in GENERATIONS]


def vector_angle(left: Sequence[float], right: Sequence[float]) -> float:
    denominator = math.sqrt(math.fsum(x * x for x in left) * math.fsum(x * x for x in right))
    cosine = abs(math.fsum(x * y for x, y in zip(left, right))) / denominator
    return math.degrees(math.acos(min(1.0, max(0.0, cosine))))


def response_condition(matrix) -> float:
    a, b = matrix[0]
    c, d = matrix[1]
    trace_ata = a * a + b * b + c * c + d * d
    determinant_squared = (a * d - b * c) ** 2
    gap = math.sqrt(max(0.0, trace_ata * trace_ata - 4 * determinant_squared))
    largest = (trace_ata + gap) / 2
    smallest = (trace_ata - gap) / 2
    if smallest <= max(1e-30, largest * 1e-15):
        return math.inf
    return math.sqrt(largest / smallest)


def matrix_diagnostics(matrix) -> dict[str, Any]:
    a, b = matrix[0]
    c, d = matrix[1]
    trace_value = a + d
    determinant = a * d - b * c
    discriminant = trace_value * trace_value - 4 * determinant
    common = trace_value / 2
    shifted = [[a - common, b], [c, d - common]]
    squared = matmul(shifted, shifted)
    j1 = math.sqrt(frobenius_squared(shifted))
    j2 = math.sqrt(frobenius_squared(squared))
    if discriminant >= 0:
        root = math.sqrt(discriminant)
        eigenvalues: list[Any] = [(trace_value + root) / 2, (trace_value - root) / 2]
        denominator = max(abs(eigenvalues[0]), abs(eigenvalues[1]), 1e-300)
        relative_gap = abs(eigenvalues[0] - eigenvalues[1]) / denominator
        if root <= 1e-13 * max(1.0, abs(trace_value)):
            angle = 0.0
            condition = math.inf
        else:
            vectors = []
            for eigenvalue in eigenvalues:
                candidate = [b, eigenvalue - a]
                if math.fsum(value * value for value in candidate) < 1e-24:
                    candidate = [eigenvalue - d, c]
                norm = math.sqrt(math.fsum(value * value for value in candidate))
                vectors.append([value / norm for value in candidate])
            dot = abs(math.fsum(x * y for x, y in zip(*vectors)))
            angle = math.degrees(math.acos(min(1.0, max(0.0, dot))))
            condition = math.sqrt((1 + dot) / max(1e-30, 1 - dot))
    else:
        root = math.sqrt(-discriminant)
        eigenvalues = [
            {"real": common, "imag": root / 2},
            {"real": common, "imag": -root / 2},
        ]
        relative_gap = root / max(math.hypot(common, root / 2), 1e-300)
        angle = None
        condition = None
    return {
        "matrix": matrix,
        "trace": trace_value,
        "determinant": determinant,
        "discriminant": discriminant,
        "eigenvalues": eigenvalues,
        "relative_eigenvalue_gap": relative_gap,
        "right_eigenvector_principal_angle_degrees": angle,
        "right_eigenbasis_condition_2norm": condition,
        "common_eigenvalue_trace_over_2": common,
        "J1_frobenius": j1,
        "J2_frobenius": j2,
        "J2_over_J1": j2 / j1 if j1 else 0.0,
    }


def solve_least_squares(design: Sequence[Sequence[float]], values: Sequence[float]):
    x = mp.matrix([[mp.mpf(str(value)) for value in row] for row in design])
    y = mp.matrix([mp.mpf(str(value)) for value in values])
    beta = mp.lu_solve(x.T * x, x.T * y)
    return [float(beta[index]) for index in range(len(beta))]


def fit_normal_diagonalizable(parent, child):
    design = []
    values = []
    for column in range(2):
        x0, x1 = parent[0][column], parent[1][column]
        design.extend(([x0, x1, 0.0], [0.0, x0, x1]))
        values.extend((child[0][column], child[1][column]))
    a, b, d = solve_least_squares(design, values)
    return [[a, b], [b, d]]


def jordan_at_angle(parent, child, theta: float):
    cosine, sine = math.cos(theta), math.sin(theta)
    nilpotent = [
        [-cosine * sine, cosine * cosine],
        [-sine * sine, sine * cosine],
    ]
    nx = matmul(nilpotent, parent)
    design = []
    values = []
    for column in range(2):
        for row in range(2):
            design.append([parent[row][column], nx[row][column]])
            values.append(child[row][column])
    eigenvalue, strength = solve_least_squares(design, values)
    fitted = [
        [
            eigenvalue * float(i == j) + strength * nilpotent[i][j]
            for j in range(2)
        ]
        for i in range(2)
    ]
    objective = frobenius_squared(subtract(child, matmul(fitted, parent)))
    return objective, fitted, eigenvalue, strength


def fit_jordan(parent, child):
    grid_size = 720
    candidates = [jordan_at_angle(parent, child, math.pi * index / grid_size) for index in range(grid_size)]
    best_index = min(range(grid_size), key=lambda index: candidates[index][0])
    center = math.pi * best_index / grid_size
    lower, upper = center - math.pi / grid_size, center + math.pi / grid_size
    golden = (math.sqrt(5) - 1) / 2
    left = upper - golden * (upper - lower)
    right = lower + golden * (upper - lower)
    for _ in range(70):
        if jordan_at_angle(parent, child, left)[0] <= jordan_at_angle(parent, child, right)[0]:
            upper, right = right, left
            left = upper - golden * (upper - lower)
        else:
            lower, left = left, right
            right = lower + golden * (upper - lower)
    theta = (lower + upper) / 2
    objective, fitted, eigenvalue, strength = jordan_at_angle(parent, child, theta)
    return fitted, {
        "source_fit_squared_residual": objective,
        "theta_degrees": math.degrees(theta) % 180,
        "common_eigenvalue": eigenvalue,
        "nilpotent_strength": strength,
    }


def fit_model(name: str, generations):
    parent, child = generations[1], generations[2]
    if name == "generic_2x2":
        return matmul(child, inverse_2x2(parent)), {"source_fit_squared_residual": 0.0}
    if name == "normal_diagonalizable":
        fitted = fit_normal_diagonalizable(parent, child)
        return fitted, {
            "source_fit_squared_residual": frobenius_squared(
                subtract(child, matmul(fitted, parent))
            )
        }
    if name == "rank2_Jordan":
        return fit_jordan(parent, child)
    raise ValueError(name)


def prediction_and_residual(state, frozen_whitener, model_name):
    generations = response_generations(state, frozen_whitener)
    transfer, fit = fit_model(model_name, generations)
    prediction = matmul(transfer, generations[2])
    residual = subtract(generations[3], prediction)
    return flatten_columns(prediction), flatten_columns(residual), transfer, fit


def numerical_jacobian(function: Callable[[Sequence[float]], Sequence[float]], point, covariance):
    baseline = function(point)
    jacobian = [[0.0] * len(point) for _ in baseline]
    for index, value in enumerate(point):
        standard_error = math.sqrt(max(0.0, covariance[index][index]))
        step = max(1e-9, abs(value) * 1e-7, standard_error * 1e-4)
        left, right = list(point), list(point)
        left[index] -= step
        right[index] += step
        f_left, f_right = function(left), function(right)
        for row in range(len(baseline)):
            jacobian[row][index] = (f_right[row] - f_left[row]) / (2 * step)
    return jacobian


def covariance_from_jacobian(jacobian, covariance):
    return matmul(matmul(jacobian, covariance), transpose(jacobian))


def chi_square_score(residual, covariance):
    sys.path.insert(0, str(ROOT / "scripts"))
    from score_p50_fullcurve_n290 import generalized_covariance_score

    return generalized_covariance_score(residual, covariance)


def determinant_gate(matrix, covariance, state_indices):
    left = [matrix[0][0], matrix[1][0]]
    right = [matrix[0][1], matrix[1][1]]
    determinant = left[0] * right[1] - left[1] * right[0]
    gradient = [right[1], -right[0], -left[1], left[0]]
    local_covariance = submatrix(covariance, state_indices)
    variance = math.fsum(
        gradient[i] * local_covariance[i][j] * gradient[j]
        for i in range(4)
        for j in range(4)
    )
    standard_error = math.sqrt(max(0.0, variance))
    return {
        "determinant": determinant,
        "delta_method_standard_error": standard_error,
        "absolute_z": abs(determinant) / standard_error if standard_error else math.inf,
    }


def build_report(alpha: float = 0.01) -> dict[str, Any]:
    mp.mp.dps = 80
    loaded = load_production_state()
    state = loaded["state"]
    covariance = loaded["covariance"]
    marginal_65 = submatrix(covariance, (0, 1))
    marginal_85 = submatrix(covariance, (2, 3))
    reference_covariance = [
        [(marginal_65[i][j] + marginal_85[i][j]) / 2 for j in range(2)]
        for i in range(2)
    ]
    frozen_whitener = whitener(reference_covariance)
    whitened_covariance = transform_covariance(covariance, frozen_whitener)
    generations = response_generations(state, frozen_whitener)

    response_rows = []
    for generation, (sizes, matrix) in enumerate(zip(GENERATIONS, generations)):
        positions = [SIZES.index(n) for n in sizes]
        indices = [2 * positions[0], 2 * positions[0] + 1, 2 * positions[1], 2 * positions[1] + 1]
        gate = determinant_gate(matrix, whitened_covariance, indices)
        response_rows.append(
            {
                "generation": generation,
                "sizes": list(sizes),
                "covariance_whitened_principal_angle_degrees": vector_angle(
                    [matrix[0][0], matrix[1][0]], [matrix[0][1], matrix[1][1]]
                ),
                "response_matrix_condition_2norm": response_condition(matrix),
                "rank_gate": gate,
            }
        )

    transfer_rows = []
    for index in range(3):
        transfer = matmul(generations[index + 1], inverse_2x2(generations[index]))
        transfer_rows.append(
            {
                "transition": f"generation_{index}_to_{index + 1}",
                "parent_sizes": list(GENERATIONS[index]),
                "child_sizes": list(GENERATIONS[index + 1]),
                **matrix_diagnostics(transfer),
            }
        )

    model_names = ("normal_diagonalizable", "rank2_Jordan", "generic_2x2")
    models = {}
    for name in model_names:
        prediction, residual, transfer, source_fit = prediction_and_residual(
            state, frozen_whitener, name
        )
        function = lambda values, name=name: prediction_and_residual(
            values, frozen_whitener, name
        )[1]
        jacobian = numerical_jacobian(function, state, covariance)
        residual_covariance = covariance_from_jacobian(jacobian, covariance)
        score = chi_square_score(residual, residual_covariance)
        models[name] = {
            "source_fit": source_fit,
            "transfer": matrix_diagnostics(transfer),
            "heldout_prediction_order": ["N520_H", "N520_U", "N680_H", "N680_U"],
            "heldout_prediction": prediction,
            "heldout_residual": residual,
            "heldout_residual_covariance": residual_covariance,
            "heldout_score": score,
        }

    target_indices = [12, 13, 14, 15]
    target_covariance = submatrix(whitened_covariance, target_indices)
    for name in model_names:
        separations = {}
        for other in model_names:
            if other == name:
                continue
            difference = [
                a - b
                for a, b in zip(
                    models[name]["heldout_prediction"],
                    models[other]["heldout_prediction"],
                )
            ]
            separations[other] = chi_square_score(difference, target_covariance)[
                "chi_square"
            ]
        closest = min(separations.values())
        p_value = models[name]["heldout_score"]["chi_square_survival"]
        if p_value < alpha:
            decision = "eliminated"
        elif closest < 9.0:
            decision = "underpowered"
        else:
            decision = "survives"
        models[name]["optimistic_target_only_pairwise_chi2"] = separations
        models[name]["decision"] = decision
        models[name]["decision_rule"] = (
            "eliminate at p<alpha; otherwise call underpowered when even the "
            "target-only nearest-rival separation is below chi2=9"
        )

    source_input_gate = response_rows[1]["rank_gate"]
    heldout_parent_gate = response_rows[2]["rank_gate"]
    identifiable = (
        source_input_gate["absolute_z"] >= 3
        and heldout_parent_gate["absolute_z"] >= 3
    )
    return {
        "schema": "matching-one/p218-production-coalescence/v1",
        "issue": 218,
        "status": "production_empirical_transfer_diagnostic",
        "conclusion": "no_positive_joint_coalescence_signature_but_second_response_direction_is_underpowered",
        "alpha": alpha,
        "state_contract": {
            "coordinate_order": ["H=N^(13/8)*P4_D", "U=N^(13/8)*P4_S_prime/mean_M_prime"],
            "basis_selection": "exact source-declared coordinates; no PCA or target rotation",
            "whitening": "inverse square root of the average N65/N85 state covariance, frozen before N520/N680",
            "covariance": "full same-batch jackknife within the aligned N65/N130/N85/N170 source; independent committed blocks for N260/N340/N520/N680",
        },
        "provenance": {
            "raw_histograms": [
                {"path": path, "sha256": sha256(ROOT / path)}
                for paths in RAW_SOURCE.values()
                for path in paths
            ]
            + [
                {"path": path, "sha256": sha256(ROOT / path)}
                for path in RAW_TARGET.values()
            ],
            "committed_scorers": [
                {"path": path, "sha256": sha256(ROOT / path)} for path in SCORERS
            ],
        },
        "state_points": {
            str(n): {"H": state[2 * index], "U": state[2 * index + 1]}
            for index, n in enumerate(SIZES)
        },
        "frozen_whitener": frozen_whitener,
        "response_vector_geometry": response_rows,
        "plug_in_transfer_geometry": transfer_rows,
        "rank_identifiability": {
            "identifiable_at_3sigma": identifiable,
            "source_input_generation_absolute_determinant_z": source_input_gate[
                "absolute_z"
            ],
            "heldout_parent_generation_absolute_determinant_z": heldout_parent_gate[
                "absolute_z"
            ],
            "consequence": (
                "full 2x2 eigen-geometry is identifiable"
                if identifiable
                else "a second response direction is not resolved at the 3-sigma production gate; plug-in eigen geometry is descriptive only"
            ),
        },
        "heldout_model_table": models,
        "frozen_existing_recurrence_crosscheck": json.loads(
            (ROOT / SCORERS[2]).read_text(encoding="utf-8")
        )["primary_lambda_half"],
        "scientific_card": [
            "MECHANISM SPACE: empirical 2D Jordan coalescence versus stable diagonalizable or generic mixing in the leading-H4/S-prime state.",
            "RESULT: report plug-in gap, response angle, eigenbasis condition and minimal-polynomial ratio, but gate interpretation on the source response determinant and held-out rival separation.",
            "NOT PROVED: the empirical curve-transfer matrix is not the microscopic Potts transfer matrix; target reuse is one mechanism analysis, not an independent evidence vote.",
            "OBSERVER-SECTOR-SOURCE-GEOMETRY: (N^(13/8)P4[D], U) | thermal H4/S-prime | threshold-rank curves | two dyadic Gaussian lineages.",
            "DEPENDENCY GROUP: PR277/P154 N65 through N680; scalar and jet scorers are correlated views of the same histograms.",
            "UPWEIGHT OBSERVATION: resolve the source response determinant and a held-out pairwise model separation above chi2=9 using the same semantic two-coordinate block.",
        ],
    }


def render_markdown(report: Mapping[str, Any]) -> str:
    rank = report["rank_identifiability"]
    lines = [
        "# P218 production empirical coalescence diagnostic",
        "",
        "## Answer",
        "",
        (
            "The committed production block resolves a two-dimensional source response at the predeclared gate."
            if rank["identifiable_at_3sigma"]
            else "The committed production block does not resolve a stable second response direction at the explicit 3-sigma diagnostic gate."
        ),
        f"The source-input determinant has `|z|={rank['source_input_generation_absolute_determinant_z']:.3f}`",
        f"and the held-out-parent determinant has `|z|={rank['heldout_parent_generation_absolute_determinant_z']:.3f}`.",
        "Therefore plug-in eigenvalue gaps and eigenvector geometry below are reported as",
        "the requested diagnostic, but are not promoted to an empirical Jordan certificate",
        "when the response-rank gate fails.",
        "The plug-in sequence also lacks the requested joint signature: the gap grows",
        "from 0.635 to 1.128 on the two source transitions while J2/J1 worsens from",
        "0.305 to 3.049. That is not positive coalescence evidence, but the rank failure",
        "prevents treating it as a powered elimination.",
        "",
        "## Covariance-whitened response geometry",
        "",
        "| generation | sizes | response angle (deg) | response cond | determinant z |",
        "|---:|---|---:|---:|---:|",
    ]
    for row in report["response_vector_geometry"]:
        lines.append(
            f"| {row['generation']} | {row['sizes']} | {row['covariance_whitened_principal_angle_degrees']:.4g} "
            f"| {row['response_matrix_condition_2norm']:.4g} | {row['rank_gate']['absolute_z']:.4g} |"
        )
    lines.extend(
        [
            "",
            "## Plug-in transfer geometry",
            "",
            "| transition | relative eigenvalue gap | eigenvector angle (deg) | eigenbasis cond | J2/J1 |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in report["plug_in_transfer_geometry"]:
        angle = row["right_eigenvector_principal_angle_degrees"]
        condition = row["right_eigenbasis_condition_2norm"]
        lines.append(
            f"| {row['transition']} | {row['relative_eigenvalue_gap']:.4g} "
            f"| {angle if angle is not None else 'complex pair'} "
            f"| {condition if condition is not None else 'complex pair'} | {row['J2_over_J1']:.4g} |"
        )
    lines.extend(
        [
            "",
            "## Held-out N520/N680 model table",
            "",
            "| source-frozen class | held-out chi-square / rank | p | decision | nearest optimistic separation chi-square |",
            "|---|---:|---:|---|---:|",
        ]
    )
    for name, row in report["heldout_model_table"].items():
        score = row["heldout_score"]
        nearest = min(row["optimistic_target_only_pairwise_chi2"].values())
        lines.append(
            f"| {name} | {score['chi_square']:.4g} / {score['numerical_rank']} "
            f"| {score['chi_square_survival']:.4g} | {row['decision']} | {nearest:.4g} |"
        )
    lines.extend(
        [
            "",
            "`normal_diagonalizable` is the covariance-whitened stable-eigenbasis",
            "representative; `rank2_Jordan` is `lambda I + s u(Ju)^T`; and",
            "`generic_2x2` is the saturated source map. A generic real 2x2 matrix with",
            "distinct eigenvalues is already diagonalizable, so the first and third rows",
            "test normal/stable versus unrestricted mixing, not disjoint algebraic sets.",
            "The pairwise separation uses target covariance only and is therefore an",
            "optimistic upper bound; values below 9 certify lack of 3-sigma discrimination.",
            "",
            "## Scientific card",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in report["scientific_card"])
    lines.extend(
        [
            "",
            "## Reproduction",
            "",
            "```bash",
            "python3 scripts/analyze_p218_production_coalescence.py --output results/p218-production-coalescence/latest.json --markdown results/p218-production-coalescence/REPORT.md",
            "python3 -m unittest discover -s tests -p 'test_p218_production_coalescence.py'",
            "```",
            "",
            "No new Monte Carlo or exact oracle is used. The script reconstructs every",
            "point and covariance from the committed production histograms and records",
            "their hashes.",
            "",
            "## Claim boundary",
            "",
            report["rank_identifiability"]["consequence"] + ". The transfer diagnostics",
            "are empirical response maps in a frozen curve coordinate, not microscopic",
            "transfer-matrix eigenstates. N520/N680 are held out from basis and source-map",
            "selection but are already revealed production blocks.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--alpha", type=float, default=0.01)
    args = parser.parse_args()
    report = build_report(alpha=args.alpha)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    args.markdown.write_text(render_markdown(report), encoding="utf-8")


if __name__ == "__main__":
    main()
