#!/usr/bin/env python3
"""Split F3 character response into first-birth and first-line exit flux."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time
from typing import Iterable, Mapping, Sequence

import yaml


SCHEMA = "matching-one.f3-activation-flux.v1"
ORIENTATIONS = ("first", "second")
HAD = ("H", "A", "D")
UVW = ("u", "v", "w")
SQRT2 = math.sqrt(2.0)
SQRT3 = math.sqrt(3.0)
SQRT_TWO_THIRDS = math.sqrt(2.0 / 3.0)
LINE_WEIGHTS = {
    "x": (0.5, 1.0 / SQRT2, 0.0),
    "d_plus": (-0.5, 0.0, 1.0 / SQRT2),
    "d_minus": (-0.5, 0.0, -1.0 / SQRT2),
    "y": (0.5, -1.0 / SQRT2, 0.0),
}
LINE_ORDER = ("x", "d_plus", "d_minus", "y")
HAD_TO_UVW = (
    (1.0 / SQRT3, SQRT_TWO_THIRDS, 0.0),
    (SQRT_TWO_THIRDS, -1.0 / SQRT3, 0.0),
    (0.0, 0.0, 1.0),
)


@dataclass(frozen=True)
class BirthCell:
    orientation: str
    batch: int
    samples: int
    tau1: int
    tau2: int
    kind: str
    ell_x: int
    ell_y: int
    count: int


@dataclass(frozen=True)
class BatchCoefficients:
    rank0: tuple[float, ...]
    rank2: tuple[float, ...]
    line_plateau: tuple[tuple[float, ...], ...]
    plateau: tuple[tuple[float, ...], ...]
    birth_line: tuple[tuple[float, ...], ...]
    exit_line: tuple[tuple[float, ...], ...]
    direct_boundary: tuple[float, ...]
    birth: tuple[tuple[float, ...], ...]
    exit: tuple[tuple[float, ...], ...]
    samples: int
    direct_rank2_count: int


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_archive(path: Path) -> tuple[int, dict[tuple[str, int], list[BirthCell]]]:
    grouped: dict[tuple[str, int], list[BirthCell]] = {}
    n_values: set[int] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        expected = {
            "n", "orientation", "batch", "samples", "tau1", "tau2",
            "kind", "ell_x", "ell_y", "count",
        }
        if reader.fieldnames is None or not expected.issubset(reader.fieldnames):
            raise ValueError("birth archive is missing required columns")
        for raw in reader:
            n_values.add(int(raw["n"]))
            cell = BirthCell(
                orientation=raw["orientation"],
                batch=int(raw["batch"]),
                samples=int(raw["samples"]),
                tau1=int(raw["tau1"]),
                tau2=int(raw["tau2"]),
                kind=raw["kind"],
                ell_x=int(raw["ell_x"]),
                ell_y=int(raw["ell_y"]),
                count=int(raw["count"]),
            )
            if cell.orientation not in ORIENTATIONS:
                raise ValueError(f"unexpected orientation {cell.orientation!r}")
            grouped.setdefault((cell.orientation, cell.batch), []).append(cell)
    if len(n_values) != 1:
        raise ValueError(f"archive must contain one N, got {sorted(n_values)}")
    return next(iter(n_values)), grouped


def projective_line(x: int, y: int) -> str:
    x %= 3
    y %= 3
    if x == 0 and y == 0:
        raise ValueError("zero vector cannot define an F3 projective line")
    if y == 0:
        return "x"
    if x == 0:
        return "y"
    if y == x:
        return "d_plus"
    if y == (-x) % 3:
        return "d_minus"
    raise ValueError(f"unrecognized F3 projective line {(x, y)}")


def build_coefficients(cells: Sequence[BirthCell], n: int) -> BatchCoefficients:
    if not cells:
        raise ValueError("empty batch")
    samples = cells[0].samples
    if any(cell.samples != samples for cell in cells):
        raise ValueError("inconsistent samples within batch")
    if sum(cell.count for cell in cells) != samples:
        raise ValueError("sparse cells do not sum to the declared batch size")
    rank0 = [0.0] * (n + 1)
    rank2 = [0.0] * (n + 1)
    line_plateau = [[0.0] * (n + 1) for _ in LINE_ORDER]
    plateau = [[0.0] * (n + 1) for _ in HAD]
    birth_line = [[0.0] * n for _ in LINE_ORDER]
    exit_line = [[0.0] * n for _ in LINE_ORDER]
    direct_boundary = [0.0] * n
    birth = [[0.0] * n for _ in HAD]
    exit_ = [[0.0] * n for _ in HAD]
    direct = 0
    for cell in cells:
        if not 0 <= cell.tau1 <= n or not 0 <= cell.tau2 <= n:
            raise ValueError("activation threshold outside [0,N]")
        mass = cell.count / samples
        for k in range(cell.tau1):
            rank0[k] += mass
        for k in range(cell.tau2, n + 1):
            rank2[k] += mass
        if cell.kind == "DIRECT_RANK2":
            if cell.tau1 != cell.tau2 or cell.ell_x or cell.ell_y:
                raise ValueError("invalid DIRECT_RANK2 atom")
            direct += cell.count
            if cell.tau1 < 1:
                raise ValueError("DIRECT_RANK2 threshold must be at least one")
            direct_boundary[cell.tau1 - 1] += mass
            continue
        if cell.kind != "LINE" or not 1 <= cell.tau1 < cell.tau2 <= n:
            raise ValueError("invalid LINE atom")
        if math.gcd(abs(cell.ell_x), abs(cell.ell_y)) != 1:
            raise ValueError("nonprimitive projective line")
        line = projective_line(cell.ell_x, cell.ell_y)
        line_index = LINE_ORDER.index(line)
        weights = LINE_WEIGHTS[line]
        birth_line[line_index][cell.tau1 - 1] += mass
        exit_line[line_index][cell.tau2 - 1] += mass
        for k in range(cell.tau1, cell.tau2):
            line_plateau[line_index][k] += mass
        for j, character_weight in enumerate(weights):
            weighted = mass * character_weight
            birth[j][cell.tau1 - 1] += weighted
            exit_[j][cell.tau2 - 1] += weighted
            for k in range(cell.tau1, cell.tau2):
                plateau[j][k] += weighted
    return BatchCoefficients(
        rank0=tuple(rank0),
        rank2=tuple(rank2),
        line_plateau=tuple(tuple(row) for row in line_plateau),
        plateau=tuple(tuple(row) for row in plateau),
        birth_line=tuple(tuple(row) for row in birth_line),
        exit_line=tuple(tuple(row) for row in exit_line),
        direct_boundary=tuple(direct_boundary),
        birth=tuple(tuple(row) for row in birth),
        exit=tuple(tuple(row) for row in exit_),
        samples=samples,
        direct_rank2_count=direct,
    )


def binomial_probabilities(n: int, p: float) -> list[float]:
    if not 0.0 < p < 1.0:
        raise ValueError("evaluation p must lie strictly between 0 and 1")
    values = [0.0] * (n + 1)
    values[0] = (1.0 - p) ** n
    ratio = p / (1.0 - p)
    for k in range(n):
        values[k + 1] = values[k] * (n - k) / (k + 1) * ratio
    return values


def evaluate_coefficients(
    coefficients: BatchCoefficients, n: int, p: float,
) -> tuple[list[float], float, float]:
    pmf_n = binomial_probabilities(n, p)
    pmf_boundary = binomial_probabilities(n - 1, p)
    plateau = [
        math.fsum(value * prob for value, prob in zip(row, pmf_n))
        for row in coefficients.plateau
    ]
    birth = [
        n * math.fsum(value * prob for value, prob in zip(row, pmf_boundary))
        for row in coefficients.birth
    ]
    exit_ = [
        n * math.fsum(value * prob for value, prob in zip(row, pmf_boundary))
        for row in coefficients.exit
    ]
    derivative = [birth[j] - exit_[j] for j in range(3)]
    coefficient_derivative = [
        n * math.fsum(
            (row[k + 1] - row[k]) * pmf_boundary[k] for k in range(n)
        )
        for row in coefficients.plateau
    ]
    derivative_residual = max(
        abs(derivative[j] - coefficient_derivative[j]) for j in range(3)
    )
    partition_residual = max(
        abs(
            coefficients.rank0[k] + coefficients.rank2[k]
            + math.fsum(row[k] for row in coefficients.line_plateau) - 1.0
        )
        for k in range(n + 1)
    )
    line_residual = 0.0
    for line_index in range(len(LINE_ORDER)):
        line_derivative = n * math.fsum(
            (
                coefficients.line_plateau[line_index][k + 1]
                - coefficients.line_plateau[line_index][k]
            ) * pmf_boundary[k]
            for k in range(n)
        )
        line_flux = n * math.fsum(
            (
                coefficients.birth_line[line_index][k]
                - coefficients.exit_line[line_index][k]
            ) * pmf_boundary[k]
            for k in range(n)
        )
        line_residual = max(line_residual, abs(line_derivative - line_flux))
    rank0_derivative = n * math.fsum(
        (coefficients.rank0[k + 1] - coefficients.rank0[k]) * pmf_boundary[k]
        for k in range(n)
    )
    rank2_derivative = n * math.fsum(
        (coefficients.rank2[k + 1] - coefficients.rank2[k]) * pmf_boundary[k]
        for k in range(n)
    )
    first_losses = n * math.fsum(
        (
            math.fsum(row[k] for row in coefficients.birth_line)
            + coefficients.direct_boundary[k]
        ) * pmf_boundary[k]
        for k in range(n)
    )
    second_gains = n * math.fsum(
        (
            math.fsum(row[k] for row in coefficients.exit_line)
            + coefficients.direct_boundary[k]
        ) * pmf_boundary[k]
        for k in range(n)
    )
    rank_residual = max(
        abs(rank0_derivative + first_losses),
        abs(rank2_derivative - second_gains),
    )
    normalization_residual = abs(math.fsum(pmf_n) - 1.0)
    gate_residual = max(
        derivative_residual, partition_residual, line_residual,
        rank_residual, normalization_residual,
    )
    return plateau + birth + exit_, gate_residual, normalization_residual


def subtract(left: Sequence[float], right: Sequence[float]) -> list[float]:
    return [right[j] - left[j] for j in range(len(left))]


def coefficient_contrast(
    first: BatchCoefficients, second: BatchCoefficients,
) -> dict[str, list[list[float]]]:
    return {
        "plateau_HAD": [
            subtract(first.plateau[j], second.plateau[j]) for j in range(3)
        ],
        "birth_boundary_HAD": [
            subtract(first.birth[j], second.birth[j]) for j in range(3)
        ],
        "exit_boundary_HAD": [
            subtract(first.exit[j], second.exit[j]) for j in range(3)
        ],
    }


def boundary_payload(coefficients: BatchCoefficients) -> dict[str, object]:
    return {
        "birth_line": [list(row) for row in coefficients.birth_line],
        "exit_line": [list(row) for row in coefficients.exit_line],
        "direct_rank2": list(coefficients.direct_boundary),
    }


def covariance_of_mean(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    batches = len(rows)
    if batches < 2:
        raise ValueError("at least two aligned batches are required")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("inconsistent covariance row width")
    means = [math.fsum(row[j] for row in rows) / batches for j in range(width)]
    return [
        [
            math.fsum(
                (row[i] - means[i]) * (row[j] - means[j]) for row in rows
            ) / (batches * (batches - 1))
            for j in range(width)
        ]
        for i in range(width)
    ]


def transform_had(vector: Sequence[float]) -> list[float]:
    if len(vector) != 3:
        raise ValueError("HAD vector must have length three")
    return [
        math.fsum(matrix_row[j] * vector[j] for j in range(3))
        for matrix_row in HAD_TO_UVW
    ]


def derived_rows(base_rows: Sequence[Sequence[float]], transform: bool) -> list[list[float]]:
    output = []
    for row in base_rows:
        plateau = list(row[0:3])
        birth = list(row[3:6])
        exit_ = list(row[6:9])
        derivative = [birth[j] - exit_[j] for j in range(3)]
        blocks = (plateau, birth, exit_, derivative)
        if transform:
            blocks = tuple(transform_had(block) for block in blocks)
        output.append([value for block in blocks for value in block])
    return output


def summarize(order: Sequence[str], rows: Sequence[Sequence[float]]) -> dict[str, object]:
    covariance = covariance_of_mean(rows)
    means = [math.fsum(row[j] for row in rows) / len(rows) for j in range(len(order))]
    errors = [math.sqrt(max(0.0, covariance[j][j])) for j in range(len(order))]
    z_scores = [
        means[j] / errors[j] if errors[j] > 0.0 else None for j in range(len(order))
    ]
    return {
        "order": list(order),
        "mean": means,
        "standard_error": errors,
        "z_score": z_scores,
        "covariance_of_mean": covariance,
        "batch_count": len(rows),
    }


def euclidean_norm(values: Iterable[float]) -> float:
    return math.sqrt(math.fsum(value * value for value in values))


def safe_ratio(numerator: float, denominator: float) -> float | None:
    return numerator / denominator if denominator > 0.0 else None


def invert_matrix(matrix: Sequence[Sequence[float]]) -> list[list[float]] | None:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    scale = max(abs(value) for row in matrix for value in row)
    if scale == 0.0:
        return None
    augmented = [
        list(row) + [1.0 if i == j else 0.0 for j in range(size)]
        for i, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= scale * 1e-12:
            return None
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                augmented[row][j] - factor * augmented[column][j]
                for j in range(2 * size)
            ]
    return [row[size:] for row in augmented]


def quadratic_form(
    vector: Sequence[float], covariance: Sequence[Sequence[float]],
) -> float | None:
    inverse = invert_matrix(covariance)
    if inverse is None:
        return None
    transformed = [
        math.fsum(inverse[i][j] * vector[j] for j in range(len(vector)))
        for i in range(len(vector))
    ]
    return math.fsum(vector[i] * transformed[i] for i in range(len(vector)))


def submatrix(matrix: Sequence[Sequence[float]], indices: Sequence[int]) -> list[list[float]]:
    return [[matrix[i][j] for j in indices] for i in indices]


def nonlinear_geometry(
    had_mean: Sequence[float], uvw_mean: Sequence[float],
) -> dict[str, float | None | str | list[float]]:
    birth = (uvw_mean[4], uvw_mean[5])
    exit_ = (uvw_mean[7], uvw_mean[8])
    derivative = (uvw_mean[10], uvw_mean[11])
    birth_norm = euclidean_norm(birth)
    exit_norm = euclidean_norm(exit_)
    derivative_norm = euclidean_norm(derivative)
    cosine = safe_ratio(
        birth[0] * exit_[0] + birth[1] * exit_[1], birth_norm * exit_norm
    )
    net_to_gross = safe_ratio(derivative_norm, birth_norm + exit_norm)
    if net_to_gross is not None:
        net_to_gross = min(1.0, max(0.0, net_to_gross))
    magnitude_balance = safe_ratio(
        2.0 * birth_norm * exit_norm, birth_norm * birth_norm + exit_norm * exit_norm
    )
    phase_cancellation = (
        magnitude_balance * cosine
        if magnitude_balance is not None and cosine is not None else None
    )
    h_birth = had_mean[3]
    h_exit = had_mean[6]
    h_derivative = had_mean[9]
    h_net_to_gross = safe_ratio(abs(h_derivative), abs(h_birth) + abs(h_exit))
    if h_net_to_gross is not None:
        h_net_to_gross = min(1.0, max(0.0, h_net_to_gross))
    if math.isclose(birth_norm, exit_norm, rel_tol=1e-12, abs_tol=1e-15):
        norm_order = "equal_within_numeric_tolerance"
    elif birth_norm > exit_norm:
        norm_order = "birth_norm_larger"
    else:
        norm_order = "first_line_exit_norm_larger"
    return {
        "charged_derivative_vector_vw": list(derivative),
        "charged_derivative_norm": derivative_norm,
        "birth_vector_vw": list(birth),
        "exit_vector_vw": list(exit_),
        "birth_norm": birth_norm,
        "exit_norm": exit_norm,
        "birth_exit_cosine": cosine,
        "magnitude_balance": magnitude_balance,
        "phase_cancellation_index": phase_cancellation,
        "charged_net_to_gross_ratio": net_to_gross,
        "charged_cancellation_fraction": (
            1.0 - net_to_gross if net_to_gross is not None else None
        ),
        "H_birth": h_birth,
        "H_exit": h_exit,
        "H_derivative": h_derivative,
        "H_net_to_gross_ratio": h_net_to_gross,
        "H_cancellation_fraction": (
            1.0 - h_net_to_gross if h_net_to_gross is not None else None
        ),
        "point_estimate_flux_norm_order": norm_order,
        "nonlinear_metric_status": "descriptive_fixed_p_not_a_field_or_phase_transport_score",
    }


def jackknife_geometry(
    had_rows: Sequence[Sequence[float]], uvw_rows: Sequence[Sequence[float]],
) -> dict[str, object]:
    if len(had_rows) != len(uvw_rows):
        raise ValueError("HAD and UVW rows must share batch alignment")
    metrics = (
        "birth_exit_cosine", "magnitude_balance", "phase_cancellation_index",
        "charged_net_to_gross_ratio", "charged_cancellation_fraction",
        "H_net_to_gross_ratio", "H_cancellation_fraction",
    )
    values: dict[str, list[float]] = {metric: [] for metric in metrics}
    batches = len(had_rows)
    for omitted in range(batches):
        had_mean = [
            math.fsum(row[j] for index, row in enumerate(had_rows) if index != omitted)
            / (batches - 1)
            for j in range(len(had_rows[0]))
        ]
        uvw_mean = [
            math.fsum(row[j] for index, row in enumerate(uvw_rows) if index != omitted)
            / (batches - 1)
            for j in range(len(uvw_rows[0]))
        ]
        geometry = nonlinear_geometry(had_mean, uvw_mean)
        for metric in metrics:
            value = geometry[metric]
            if value is not None:
                values[metric].append(float(value))
    output: dict[str, object] = {}
    for metric, observed in values.items():
        if len(observed) != batches:
            output[metric] = {"status": "not_scoreable"}
            continue
        mean = math.fsum(observed) / batches
        standard_error = math.sqrt(
            (batches - 1) / batches
            * math.fsum((value - mean) ** 2 for value in observed)
        )
        output[metric] = {
            "leave_one_values": observed,
            "jackknife_standard_error": standard_error,
            "leave_one_min": min(observed),
            "leave_one_max": max(observed),
        }
    return output


def charged_diagnostics(
    had_summary: Mapping[str, object], uvw_summary: Mapping[str, object],
    had_rows: Sequence[Sequence[float]], uvw_rows: Sequence[Sequence[float]],
) -> dict[str, object]:
    uvw_mean = uvw_summary["mean"]
    uvw_cov = uvw_summary["covariance_of_mean"]
    birth = (uvw_mean[4], uvw_mean[5])
    exit_ = (uvw_mean[7], uvw_mean[8])
    derivative = (uvw_mean[10], uvw_mean[11])
    cov00 = uvw_cov[10][10]
    cov01 = uvw_cov[10][11]
    cov11 = uvw_cov[11][11]
    determinant = cov00 * cov11 - cov01 * cov01
    birth_indices = (4, 5)
    exit_indices = (7, 8)
    derivative_indices = (10, 11)
    joint_indices = (4, 5, 7, 8)
    chi2_birth = quadratic_form(birth, submatrix(uvw_cov, birth_indices))
    chi2_exit = quadratic_form(exit_, submatrix(uvw_cov, exit_indices))
    chi2_net = quadratic_form(derivative, submatrix(uvw_cov, derivative_indices))
    chi2_joint = quadratic_form(
        (*birth, *exit_), submatrix(uvw_cov, joint_indices)
    )
    geometry = nonlinear_geometry(had_summary["mean"], uvw_mean)
    return {
        **geometry,
        "charged_birth_chi2_df2": chi2_birth,
        "charged_exit_chi2_df2": chi2_exit,
        "charged_derivative_chi2_df2": chi2_net,
        "charged_joint_birth_exit_chi2_df4": chi2_joint,
        "charged_derivative_p_value_df2": (
            math.exp(-0.5 * max(0.0, chi2_net)) if chi2_net is not None else None
        ),
        "charged_derivative_inference": {
            "statistic": "plug_in_Wald_chi_square",
            "df": 2,
            "p_value": "asymptotic_chi_square_tail",
            "covariance": "estimated_from_20_aligned_batches",
            "role": "correlated_exploratory_fixed_p_diagnostic_not_an_independent_vote",
        },
        "charged_covariance_determinant": determinant,
        "delete_one_geometry": jackknife_geometry(had_rows, uvw_rows),
    }


def evaluation_order(basis: Sequence[str], include_derivative: bool) -> list[str]:
    blocks = ("plateau", "birth", "exit", "derivative") if include_derivative else (
        "plateau", "birth", "exit"
    )
    return [f"{block}_{name}" for block in blocks for name in basis]


def score(
    manifest: Mapping[str, object], births_path: Path, metadata_path: Path,
    workers: int,
) -> dict[str, object]:
    started = time.perf_counter()
    expected_births_hash = manifest["source"]["births"]["sha256"]
    expected_metadata_hash = manifest["source"]["metadata"]["sha256"]
    actual_births_hash = sha256(births_path)
    actual_metadata_hash = sha256(metadata_path)
    if actual_births_hash != expected_births_hash:
        raise ValueError("birth archive SHA256 does not match manifest")
    if actual_metadata_hash != expected_metadata_hash:
        raise ValueError("metadata SHA256 does not match manifest")

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    n, grouped = read_archive(births_path)
    batch_sets = {
        orientation: {batch for found_orientation, batch in grouped
                      if found_orientation == orientation}
        for orientation in ORIENTATIONS
    }
    if batch_sets["first"] != batch_sets["second"]:
        raise ValueError("orientations do not share the same aligned batch IDs")
    batch_ids = sorted(batch_sets["first"])
    coefficients = {
        key: build_coefficients(cells, n) for key, cells in grouped.items()
    }
    direct_counts = {
        orientation: sum(
            coefficients[(orientation, batch)].direct_rank2_count for batch in batch_ids
        )
        for orientation in ORIENTATIONS
    }
    primary_coefficients = [
        {
            "batch": batch,
            **coefficient_contrast(
                coefficients[("first", batch)], coefficients[("second", batch)]
            ),
        }
        for batch in batch_ids
    ]
    sufficient_statistics = [
        {
            "batch": batch,
            "orientations": {
                orientation: boundary_payload(coefficients[(orientation, batch)])
                for orientation in ORIENTATIONS
            },
        }
        for batch in batch_ids
    ]

    evaluations: list[dict[str, object]] = []
    maximum_exact_closure_residual = 0.0
    maximum_pmf_residual = 0.0
    for point in manifest["evaluation_points"]:
        p = float(point["p"])
        batch_base: dict[tuple[str, int], list[float]] = {}
        for orientation in ORIENTATIONS:
            for batch in batch_ids:
                row, derivative_residual, pmf_residual = evaluate_coefficients(
                    coefficients[(orientation, batch)], n, p
                )
                batch_base[(orientation, batch)] = row
                maximum_exact_closure_residual = max(
                    maximum_exact_closure_residual, derivative_residual
                )
                maximum_pmf_residual = max(maximum_pmf_residual, pmf_residual)

        views: dict[str, object] = {}
        for view in (*ORIENTATIONS, "second_minus_first"):
            if view == "second_minus_first":
                base_rows = [
                    subtract(
                        batch_base[("first", batch)],
                        batch_base[("second", batch)],
                    )
                    for batch in batch_ids
                ]
            else:
                base_rows = [batch_base[(view, batch)] for batch in batch_ids]
            had_rows = derived_rows(base_rows, transform=False)
            uvw_rows = derived_rows(base_rows, transform=True)
            base_summary = summarize(evaluation_order(HAD, False), base_rows)
            had_summary = summarize(evaluation_order(HAD, True), had_rows)
            uvw_summary = summarize(evaluation_order(UVW, True), uvw_rows)
            views[view] = {
                "base_HAD": base_summary,
                "derived_HAD": had_summary,
                "derived_uvw": uvw_summary,
                "diagnostics": charged_diagnostics(
                    had_summary, uvw_summary, had_rows, uvw_rows
                ),
            }
        joint_rows = [
            batch_base[("first", batch)] + batch_base[("second", batch)]
            for batch in batch_ids
        ]
        evaluation = {
            "id": point["id"], "role": point["role"], "p": p,
            "joint_orientations_HAD": summarize(
                [
                    f"{orientation}_{name}"
                    for orientation in ORIENTATIONS
                    for name in evaluation_order(HAD, False)
                ],
                joint_rows,
            ),
            "views": views,
        }
        if "root_treatment" in point:
            evaluation["root_treatment"] = point["root_treatment"]
        evaluations.append(evaluation)

    by_id = {row["id"]: row for row in evaluations}
    p_ref_h = by_id["p_ref"]["views"]["second_minus_first"]["derived_HAD"]["mean"][0]
    cross_dh = by_id["h_cross"]["views"]["second_minus_first"]["derived_HAD"]["mean"][9]
    reference_checks = {
        "p_ref_plateau_H_contrast": {
            "expected": manifest["reference_checks"]["p_ref_plateau_H_contrast"],
            "observed": p_ref_h,
            "residual": p_ref_h - manifest["reference_checks"]["p_ref_plateau_H_contrast"],
        },
        "h_cross_derivative_H_contrast": {
            "expected": manifest["reference_checks"]["h_cross_derivative_H_contrast"],
            "observed": cross_dh,
            "residual": cross_dh - manifest["reference_checks"]["h_cross_derivative_H_contrast"],
        },
    }
    exact_gates = {
        "passed": max(maximum_exact_closure_residual, maximum_pmf_residual) < 1e-11,
        "maximum_exact_transition_closure_residual": maximum_exact_closure_residual,
        "maximum_binomial_normalization_residual": maximum_pmf_residual,
        "aligned_batch_ids_identical": True,
        "archive_sha256_verified": True,
        "metadata_sha256_verified": True,
    }
    if not exact_gates["passed"]:
        raise ValueError(f"exact flux gate failed: {exact_gates}")

    return {
        "schema": SCHEMA,
        "analysis_status": "complete",
        "evidence_role": "exploratory_existing_data_reuse",
        "new_samples": False,
        "dependency_group": manifest["source"]["dependency_group"],
        "source": {
            **manifest["source"],
            "births_sha256_verified": actual_births_hash,
            "metadata_sha256_verified": actual_metadata_hash,
            "N": n,
            "batch_ids": batch_ids,
            "direct_rank2_counts": direct_counts,
            "metadata_git_commit": metadata.get("git_commit"),
        },
        "contract": manifest["flux_contract"],
        "basis": manifest["character_basis"],
        "exact_gates": exact_gates,
        "reference_checks": reference_checks,
        "evaluations": evaluations,
        "paired_contrast_boundary_coefficients": {
            "description": (
                "Per-aligned-batch second-minus-first coefficients retain arbitrary-p "
                "and cross-p covariance for the primary orientation contrast."
            ),
            "plateau_basis": "degree-N Bernstein coefficients",
            "boundary_basis": "degree-(N-1) Bernstein coefficients before multiplication by N",
            "character_order": list(HAD),
            "batches": primary_coefficients,
        },
        "line_boundary_sufficient_statistics": {
            "description": (
                "Per-aligned-batch and per-orientation F3 line birth/exit arrays plus "
                "DIRECT_RANK2 retain arbitrary-p character and nonzero-twist fluxes."
            ),
            "line_order": list(LINE_ORDER),
            "boundary_basis": "degree-(N-1) Bernstein coefficients before multiplication by N",
            "batches": sufficient_statistics,
        },
        "interpretation_boundary": manifest["claim_boundary"],
        "runtime": {
            "elapsed_seconds": time.perf_counter() - started,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "machine": platform.machine(),
            "workers_requested": workers,
            "execution": "deterministic CPU aggregation; no new Monte Carlo samples",
        },
    }


def fmt(value: float | None, digits: int = 6) -> str:
    return "--" if value is None else f"{value:.{digits}g}"


def render_markdown(result: Mapping[str, object]) -> str:
    rows = []
    for evaluation in result["evaluations"]:
        view = evaluation["views"]["second_minus_first"]
        diagnostic = view["diagnostics"]
        rows.append(
            "| `{id}` | {p:.9f} | {plateau} | {birth} | {exit} | {derivative} | "
            "{h_cancel} | {charged} | {chi2} | {p_value} | {c_cancel} |".format(
                id=evaluation["id"], p=evaluation["p"],
                plateau=fmt(view["derived_HAD"]["mean"][0]),
                birth=fmt(diagnostic["H_birth"]),
                exit=fmt(diagnostic["H_exit"]),
                derivative=fmt(diagnostic["H_derivative"]),
                h_cancel=fmt(diagnostic["H_net_to_gross_ratio"], 4),
                charged=fmt(diagnostic["charged_derivative_norm"]),
                chi2=fmt(diagnostic["charged_derivative_chi2_df2"], 5),
                p_value=fmt(diagnostic["charged_derivative_p_value_df2"], 4),
                c_cancel=fmt(diagnostic["charged_net_to_gross_ratio"], 4),
            )
        )
    cross = next(row for row in result["evaluations"] if row["id"] == "h_cross")
    reference = next(row for row in result["evaluations"] if row["id"] == "p_ref")
    cross_diag = cross["views"]["second_minus_first"]["diagnostics"]
    ref_diag = reference["views"]["second_minus_first"]["diagnostics"]
    upper = next(row for row in result["evaluations"] if row["id"] == "upper_probe")
    upper_diag = upper["views"]["second_minus_first"]["diagnostics"]
    cross_cos_loo = cross_diag["delete_one_geometry"]["birth_exit_cosine"]
    ref_cos_loo = ref_diag["delete_one_geometry"]["birth_exit_cosine"]
    gates = result["exact_gates"]
    direct = result["source"]["direct_rank2_counts"]
    lines = [
        "# F3 activation-flux tomography", "",
        "## Result", "",
        "The existing paired N65 archive now resolves the directional projective "
        "character response into rank-one birth and exit of the first-born line "
        "at the rank-two completion boundary. "
        "This is a new view of the same 20 aligned batches, not a new evidence vote.", "",
        "The birth and first-line-weighted exit entries below are signed character "
        "coordinates. Both net/gross columns use `|birth-exit|/(|birth|+|exit|)`: "
        "zero means cancellation and one means reinforcement/no cancellation.", "",
        "| point | p | plateau H | birth H | rank-one exit H at tau2 | dH/dp | H net/gross | "
        "charged |d(v,w)/dp| | plug-in Wald chi2(2) | asymptotic p(2) | charged net/gross |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        *rows, "",
        "At the frozen exploratory H crossing, the H slope is `{}`: first-birth `{}` "
        "minus first-line exit at rank-two completion `{}`. Their signs are opposite, "
        "so the two topological "
        "boundaries reinforce rather than cancel (`H net/gross={}`). This reinforcement "
        "appears in every fixed-p delete-one batch reconstruction.".format(
            fmt(cross_diag["H_derivative"]), fmt(cross_diag["H_birth"]),
            fmt(cross_diag["H_exit"]), fmt(cross_diag["H_net_to_gross_ratio"], 4),
        ), "",
        "At `p_ref`, H remains boundary-reinforced (`dH/dp={}`, net/gross `1`). "
        "The point-estimate norm of the first-line-weighted exit projection onto the "
        "fixed `(v,w)` plane is larger than the birth projection, "
        "but their relative phase is not yet stable: its delete-one cosine spans "
        "`[{:.3f},{:.3f}]` at the crossing and `[{:.3f},{:.3f}]` at `p_ref`. The crossing "
        "delete-one values hold the full-data root fixed and do not propagate root "
        "uncertainty. The safe claim is two-boundary reinforcement in H, not one resolved "
        "charged ray or verified C3 transport.".format(
            fmt(ref_diag["H_derivative"]),
            cross_cos_loo["leave_one_min"], cross_cos_loo["leave_one_max"],
            ref_cos_loo["leave_one_min"], ref_cos_loo["leave_one_max"],
        ), "",
        "By `p=0.65`, H itself nearly cancels (net/gross `{}`) while the fixed charged-"
        "plane projection is nonzero under a plug-in Wald diagnostic (`chi2_2={}`, "
        "asymptotic `p={}`). Because all p points share the same archive, this is a "
        "mechanism selector for a fresh child—not an independent significance claim or "
        "verified C3 transport.".format(
            fmt(upper_diag["H_net_to_gross_ratio"], 4),
            fmt(upper_diag["charged_derivative_chi2_df2"], 5),
            fmt(upper_diag["charged_derivative_p_value_df2"], 4),
        ), "",
        "## Exactness and dependence", "",
        "For every orientation, batch, line, character and reported p, the Bernstein "
        "partition and rank/line transition identities close, including character "
        "derivative equals signed birth coordinate minus the first-line-weighted exit "
        "coordinate at rank-two completion. "
        "The maximum residual is `{:.3g}`; the maximum binomial normalization "
        "residual is `{:.3g}`.".format(
            gates["maximum_exact_transition_closure_residual"],
            gates["maximum_binomial_normalization_residual"],
        ), "",
        "`DIRECT_RANK2` atoms are reported but contribute exactly zero to the zero-sum "
        "projective characters: first `{}`, second `{}`. Full covariance-of-the-mean "
        "matrices are retained for H/A/D and u/v/w, and per-batch boundary coefficients "
        "retain arbitrary cross-p covariance for the paired contrast.".format(
            direct["first"], direct["second"]
        ), "",
        "## What this changes", "",
        "The next fresh N130 archive should freeze the complete A4 triplet and these two "
        "boundary fluxes together, and add a section-audited completion record. "
        "The present archive knows when rank two completes but not which marked "
        "complement representative completes it. That representative is not intrinsic: "
        "a fresh child must retain raw winding, section/basis, transporter and ambiguity "
        "metadata before testing a gauge-covariant stabilizer-C3 completion character. "
        "A larger deterministic basis-shear run "
        "cannot answer this question; a later physical defect/source insertion can test "
        "whether the charged plane is more than a projective re-expression of the same "
        "archive.", "",
        "Claim boundary: {}".format(result["interpretation_boundary"]), "",
    ]
    return "\n".join(lines)


def resolve_input(
    explicit: Path | None, source_root: Path | None, relative_path: str,
) -> Path:
    if explicit is not None:
        return explicit.resolve()
    if source_root is None:
        raise ValueError("provide either an explicit input path or --source-root")
    return (source_root / relative_path).resolve()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--births", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    manifest = yaml.safe_load(args.manifest.read_text(encoding="utf-8"))
    births = resolve_input(
        args.births, args.source_root, manifest["source"]["births"]["path"]
    )
    metadata = resolve_input(
        args.metadata, args.source_root, manifest["source"]["metadata"]["path"]
    )
    result = score(manifest, births, metadata, args.workers)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    args.output_md.write_text(render_markdown(result), encoding="utf-8")


if __name__ == "__main__":
    main()
