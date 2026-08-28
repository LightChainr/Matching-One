#!/usr/bin/env python3
"""Audit cross-size covariance in aligned threshold-rank orientation batches.

The C++ threshold-rank engine may deliberately reuse one counter stream across
different sizes. When batch ids are aligned, this script reconstructs each
batch's orientation difference and slope, and uses delete-one jackknife
pseudo-values for the nonlinear root gap. It then emits covariance matrices of
the estimators across sizes and compares full-covariance and diagonal-only
held-out constant-amplitude tests for

    A_M = N^(13/8) Delta M / Delta cos(4 theta)
    A_p = -N^2 Delta p* / Delta cos(4 theta).

Coupled scores are fail-closed: every aligned batch must share one sample
count, and same batch ids are treated as coupled only after metadata proves
that the sizes share RNG schema, seed, and replica-counter ranges. Scientific
quadratic forms use Cholesky or an SVD/eigen pseudoinverse, never Gauss-Jordan
alone. The plug-in chi-square is labelled asymptotic and is accompanied by
Hotelling/F and batch-bootstrap calibration.

This audits existing aggregates. It never changes the frozen simulation, RNG,
or model-selection protocol.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import mpmath as mp

from analyze_threshold_rank_orientation import (
    add_histograms,
    cos4,
    evaluate_histogram,
    read_histograms,
)
from analyze_threshold_ranks import matching_root


Matrix = List[List[float]]
Vector = List[float]
Record = Dict[str, object]
RecordKey = Tuple[int, str, int]
Grouped = Dict[int, Dict[str, Dict[int, Record]]]

FORMAT_VERSION = 3
DEFAULT_EIGENVALUE_CUTOFF = 1e-12
EIGENVALUE_CUTOFF_GRID = (0.0, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6)
DEFAULT_BOOTSTRAP_REPLICATES = 2000
DEFAULT_BOOTSTRAP_SEED = 20260828
CONDITION_SWITCH_TO_SVD = 1.0 / DEFAULT_EIGENVALUE_CUTOFF


def _matvec(matrix: Matrix, vector: Vector) -> Vector:
    return [math.fsum(a * b for a, b in zip(row, vector)) for row in matrix]


def _symmetrize(matrix: Matrix) -> Matrix:
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    return [
        [0.5 * (matrix[i][j] + matrix[j][i]) for j in range(n)]
        for i in range(n)
    ]


def _solve(matrix: Matrix, vector: Vector) -> Vector:
    """Gauss-Jordan solve retained only as a numerical cross-check."""

    n = len(vector)
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError("linear system must be square")
    augmented = [list(map(float, matrix[i])) + [float(vector[i])] for i in range(n)]
    scale = max((abs(value) for row in matrix for value in row), default=0.0)
    tolerance = max(scale * 1e-13, 1e-300)
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= tolerance:
            raise ArithmeticError("singular or ill-resolved matrix")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        for entry in range(column, n + 1):
            augmented[column][entry] /= divisor
        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            for entry in range(column, n + 1):
                augmented[row][entry] -= factor * augmented[column][entry]
    return [augmented[row][-1] for row in range(n)]


def _inverse(matrix: Matrix) -> Matrix:
    n = len(matrix)
    columns = [
        _solve(matrix, [1.0 if row == column else 0.0 for row in range(n)])
        for column in range(n)
    ]
    return [[columns[column][row] for column in range(n)] for row in range(n)]


def _cholesky_factor(matrix: Matrix) -> Matrix:
    symmetric = _symmetrize(matrix)
    n = len(symmetric)
    factor = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            value = symmetric[i][j] - math.fsum(
                factor[i][k] * factor[j][k] for k in range(j)
            )
            if i == j:
                if value <= 0.0 or not math.isfinite(value):
                    raise ArithmeticError("matrix is not positive definite")
                factor[i][j] = math.sqrt(value)
            else:
                factor[i][j] = value / factor[j][j]
    return factor


def _cholesky_solve(factor: Matrix, vector: Vector) -> Vector:
    n = len(vector)
    down = [0.0] * n
    for i in range(n):
        down[i] = (
            vector[i] - math.fsum(factor[i][j] * down[j] for j in range(i))
        ) / factor[i][i]
    up = [0.0] * n
    for i in range(n - 1, -1, -1):
        up[i] = (
            down[i]
            - math.fsum(factor[j][i] * up[j] for j in range(i + 1, n))
        ) / factor[i][i]
    return up


def cholesky_inverse(matrix: Matrix) -> Matrix:
    factor = _cholesky_factor(matrix)
    n = len(matrix)
    columns = [
        _cholesky_solve(factor, [1.0 if row == column else 0.0 for row in range(n)])
        for column in range(n)
    ]
    return [[columns[column][row] for column in range(n)] for row in range(n)]


def _jacobi_eigh(matrix: Matrix) -> Tuple[Vector, Matrix]:
    """Return ascending eigenvalues and columns-as-eigenvectors for a symmetric matrix."""

    symmetric = _symmetrize(matrix)
    n = len(symmetric)
    vectors = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    if n == 1:
        return [symmetric[0][0]], vectors
    scale = max((abs(value) for row in symmetric for value in row), default=0.0)
    if scale == 0.0:
        return [0.0] * n, vectors
    work = [list(row) for row in symmetric]
    tolerance = max(scale * 1e-15, 1e-300)
    for _ in range(max(32, 12 * n * n)):
        pivot_i = 0
        pivot_j = 1
        largest = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                value = abs(work[i][j])
                if value > largest:
                    largest = value
                    pivot_i, pivot_j = i, j
        if largest <= tolerance:
            break
        app = work[pivot_i][pivot_i]
        aqq = work[pivot_j][pivot_j]
        apq = work[pivot_i][pivot_j]
        theta = 0.5 * (aqq - app) / apq
        sign = -1.0 if theta < 0.0 else 1.0
        tangent = sign / (abs(theta) + math.sqrt(1.0 + theta * theta))
        cosine = 1.0 / math.sqrt(1.0 + tangent * tangent)
        sine = tangent * cosine
        work[pivot_i][pivot_i] = app - tangent * apq
        work[pivot_j][pivot_j] = aqq + tangent * apq
        work[pivot_i][pivot_j] = work[pivot_j][pivot_i] = 0.0
        for k in range(n):
            if k == pivot_i or k == pivot_j:
                continue
            aik = work[k][pivot_i]
            ajk = work[k][pivot_j]
            rotated_i = cosine * aik - sine * ajk
            rotated_j = sine * aik + cosine * ajk
            work[k][pivot_i] = work[pivot_i][k] = rotated_i
            work[k][pivot_j] = work[pivot_j][k] = rotated_j
        for k in range(n):
            vip = vectors[k][pivot_i]
            viq = vectors[k][pivot_j]
            vectors[k][pivot_i] = cosine * vip - sine * viq
            vectors[k][pivot_j] = sine * vip + cosine * viq
    eigenvalues = [work[i][i] for i in range(n)]
    order = sorted(range(n), key=lambda index: eigenvalues[index])
    eigenvalues = [eigenvalues[index] for index in order]
    vectors = [[row[index] for index in order] for row in vectors]
    return eigenvalues, vectors


def _absolute_cutoff(eigenvalues: Sequence[float], relative_cutoff: float) -> float:
    maximum = max((abs(value) for value in eigenvalues), default=0.0)
    if maximum == 0.0:
        return 0.0
    return max(relative_cutoff * maximum, 0.0)


def covariance_eigenstructure(
    matrix: Matrix,
    relative_cutoff: float = DEFAULT_EIGENVALUE_CUTOFF,
) -> Dict[str, object]:
    eigenvalues, eigenvectors = _jacobi_eigh(matrix)
    cutoff = _absolute_cutoff(eigenvalues, relative_cutoff)
    retained = [value for value in eigenvalues if abs(value) > cutoff]
    maximum = max((abs(value) for value in eigenvalues), default=0.0)
    minimum_abs = min((abs(value) for value in retained), default=0.0)
    if not retained:
        condition = math.inf
    else:
        condition = maximum / minimum_abs if minimum_abs > 0.0 else math.inf
    return {
        "eigenvalues": eigenvalues,
        "eigenvectors": eigenvectors,
        "relative_cutoff": relative_cutoff,
        "absolute_cutoff": cutoff,
        "effective_rank": len(retained),
        "dimension": len(eigenvalues),
        "minimum_eigenvalue": min(eigenvalues) if eigenvalues else 0.0,
        "maximum_eigenvalue": max(eigenvalues) if eigenvalues else 0.0,
        "condition_number": condition,
        "truncated_eigenvalue_count": len(eigenvalues) - len(retained),
        "negative_eigenvalue_count": sum(1 for value in eigenvalues if value < -cutoff),
    }


def _pseudoinverse_from_eigh(
    eigenvalues: Sequence[float],
    eigenvectors: Matrix,
    absolute_cutoff: float,
) -> Matrix:
    n = len(eigenvalues)
    scales = [
        (1.0 / value) if abs(value) > absolute_cutoff else 0.0
        for value in eigenvalues
    ]
    inverse = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            value = math.fsum(
                eigenvectors[i][k] * scales[k] * eigenvectors[j][k]
                for k in range(n)
            )
            inverse[i][j] = inverse[j][i] = value
    return inverse


def stable_inverse(
    matrix: Matrix,
    relative_cutoff: float = DEFAULT_EIGENVALUE_CUTOFF,
) -> Tuple[Matrix, Dict[str, object]]:
    """Invert a covariance with Cholesky first and SVD/eigen pseudoinverse fallback.

    Scientific scores use this path. The retained Gauss-Jordan helper is only a
    cross-check on well-conditioned problems.
    """

    diagnostics = covariance_eigenstructure(matrix, relative_cutoff)
    rank = int(diagnostics["effective_rank"])
    if rank == 0:
        raise ArithmeticError("covariance has no retained eigenvalues")
    condition = float(diagnostics["condition_number"])
    use_svd = (
        rank < int(diagnostics["dimension"])
        or int(diagnostics["negative_eigenvalue_count"]) > 0
        or not math.isfinite(condition)
        or condition > CONDITION_SWITCH_TO_SVD
    )
    if not use_svd:
        try:
            inverse = cholesky_inverse(matrix)
            diagnostics["solver"] = "cholesky"
            return inverse, diagnostics
        except (ArithmeticError, ValueError):
            use_svd = True
    inverse = _pseudoinverse_from_eigh(
        diagnostics["eigenvalues"],  # type: ignore[arg-type]
        diagnostics["eigenvectors"],  # type: ignore[arg-type]
        float(diagnostics["absolute_cutoff"]),
    )
    diagnostics["solver"] = "svd_pseudoinverse"
    return inverse, diagnostics


def _public_eigenstructure(diagnostics: Mapping[str, object]) -> Dict[str, object]:
    condition = diagnostics["condition_number"]
    return {
        "eigenvalues": list(diagnostics["eigenvalues"]),  # type: ignore[arg-type]
        "relative_cutoff": diagnostics["relative_cutoff"],
        "absolute_cutoff": diagnostics["absolute_cutoff"],
        "effective_rank": diagnostics["effective_rank"],
        "dimension": diagnostics["dimension"],
        "minimum_eigenvalue": diagnostics["minimum_eigenvalue"],
        "maximum_eigenvalue": diagnostics["maximum_eigenvalue"],
        "condition_number": condition if math.isfinite(float(condition)) else None,
        "condition_number_infinite": not math.isfinite(float(condition)),
        "truncated_eigenvalue_count": diagnostics["truncated_eigenvalue_count"],
        "negative_eigenvalue_count": diagnostics["negative_eigenvalue_count"],
        "solver": diagnostics.get("solver"),
    }


def _subset(
    matrix: Matrix,
    rows: Sequence[int],
    columns: Optional[Sequence[int]] = None,
) -> Matrix:
    if columns is None:
        columns = rows
    return [[matrix[i][j] for j in columns] for i in rows]


def _quadratic(vector: Vector, inverse_covariance: Matrix) -> float:
    return math.fsum(
        a * b for a, b in zip(vector, _matvec(inverse_covariance, vector))
    )


def truncation_sensitivity(
    vector: Vector,
    matrix: Matrix,
    cutoffs: Sequence[float] = EIGENVALUE_CUTOFF_GRID,
) -> List[Dict[str, object]]:
    eigenvalues, eigenvectors = _jacobi_eigh(matrix)
    rows = []
    for relative in cutoffs:
        cutoff = _absolute_cutoff(eigenvalues, relative)
        inverse = _pseudoinverse_from_eigh(eigenvalues, eigenvectors, cutoff)
        rank = sum(1 for value in eigenvalues if abs(value) > cutoff)
        rows.append(
            {
                "relative_cutoff": relative,
                "absolute_cutoff": cutoff,
                "effective_rank": rank,
                "plugin_chi_square": _quadratic(vector, inverse),
            }
        )
    return rows


def _regularized_betainc(a: float, b: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    return float(mp.betainc(a, b, 0, x, regularized=True))


def chi_square_survival(statistic: float, dof: int) -> float:
    if statistic <= 0.0:
        return 1.0
    if dof <= 0:
        raise ValueError("chi-square degrees of freedom must be positive")
    return float(mp.gammainc(dof / 2.0, statistic / 2.0, mp.inf, regularized=True))


def f_survival(statistic: float, df1: float, df2: float) -> float:
    if statistic <= 0.0:
        return 1.0
    if df1 <= 0.0 or df2 <= 0.0:
        raise ValueError("F degrees of freedom must be positive")
    x = (df1 * statistic) / (df1 * statistic + df2)
    return 1.0 - _regularized_betainc(df1 / 2.0, df2 / 2.0, x)


def hotelling_calibration(
    plugin_chi_square: float,
    dimension: int,
    batch_count: int,
) -> Dict[str, object]:
    """Finite-batch Hotelling/F calibration of a plug-in mean quadratic form.

    covariance_of_mean already returns Cov(mean), so the plug-in chi-square is
    Hotelling's T^2. This is a diagnostic, not a paper-facing p-value.
    """

    note = (
        "finite-batch calibration of the plug-in quadratic form; "
        "not a paper-facing p-value"
    )
    if dimension < 1 or batch_count <= dimension:
        return {
            "applicable": False,
            "reason": "need batch_count > held-out dimension",
            "note": note,
        }
    t_square = plugin_chi_square
    f_statistic = ((batch_count - dimension) / (dimension * (batch_count - 1))) * t_square
    df1 = float(dimension)
    df2 = float(batch_count - dimension)
    return {
        "applicable": True,
        "t_square": t_square,
        "f_statistic": f_statistic,
        "df": [dimension, batch_count - dimension],
        "f_survival": f_survival(f_statistic, df1, df2),
        "plugin_chi_square_survival": chi_square_survival(plugin_chi_square, dimension),
        "note": note,
    }


def covariance_of_mean(
    batch_vectors: Sequence[Sequence[float]],
) -> Tuple[Vector, Matrix]:
    """Return column means and covariance matrix of those means.

    Rows are aligned independent batches and columns are sizes/observables.
    This is an unweighted equal-batch estimator; callers must assert a common
    sample count rather than silently mixing unequal batch weights.
    """

    if len(batch_vectors) < 2:
        raise ValueError("at least two aligned batches are required")
    width = len(batch_vectors[0])
    if width == 0 or any(len(row) != width for row in batch_vectors):
        raise ValueError("batch matrix must be nonempty and rectangular")
    if any(not math.isfinite(value) for row in batch_vectors for value in row):
        raise ValueError("batch values must be finite")
    batches = len(batch_vectors)
    means = [
        math.fsum(row[column] for row in batch_vectors) / batches
        for column in range(width)
    ]
    covariance = [[0.0] * width for _ in range(width)]
    denominator = batches * (batches - 1)
    for i in range(width):
        for j in range(i, width):
            value = math.fsum(
                (row[i] - means[i]) * (row[j] - means[j])
                for row in batch_vectors
            ) / denominator
            covariance[i][j] = covariance[j][i] = value
    return means, covariance


def jackknife_pseudovalues(
    full_estimate: float,
    delete_one_estimates: Sequence[float],
) -> List[float]:
    """Convert delete-one estimates into standard jackknife pseudo-values."""

    batches = len(delete_one_estimates)
    if batches < 2:
        raise ValueError("at least two delete-one estimates are required")
    if not math.isfinite(full_estimate) or any(
        not math.isfinite(value) for value in delete_one_estimates
    ):
        raise ValueError("jackknife estimates must be finite")
    return [
        batches * full_estimate - (batches - 1) * value
        for value in delete_one_estimates
    ]


def correlations(covariance: Matrix) -> Matrix:
    output = [[0.0] * len(covariance) for _ in covariance]
    for i in range(len(covariance)):
        if covariance[i][i] <= 0:
            raise ValueError("covariance diagonal must be positive")
        for j in range(len(covariance)):
            output[i][j] = covariance[i][j] / math.sqrt(
                covariance[i][i] * covariance[j][j]
            )
    return output


def constant_heldout_audit(
    values: Sequence[float],
    covariance: Matrix,
    sizes: Sequence[int],
    training_sizes: Sequence[int],
    heldout_sizes: Sequence[int],
    batch_count: Optional[int] = None,
    eigenvalue_cutoff: float = DEFAULT_EIGENVALUE_CUTOFF,
) -> Dict[str, object]:
    """Fit one constant on training sizes and score correlated held-out residuals.

    The quadratic form is a plug-in/asymptotic score: the covariance is treated
    as known. When batch_count is supplied, Hotelling/F calibration is attached.
    Inverses use Cholesky or an SVD/eigen pseudoinverse.
    """

    if len(values) != len(sizes) or len(covariance) != len(sizes):
        raise ValueError("values/covariance/sizes dimensions disagree")
    indices = {size: index for index, size in enumerate(sizes)}
    if len(indices) != len(sizes):
        raise ValueError("sizes must be unique")
    try:
        training = [indices[size] for size in training_sizes]
        heldout = [indices[size] for size in heldout_sizes]
    except KeyError as exc:
        raise ValueError("training or held-out size is absent") from exc

    c_tt = _subset(covariance, training)
    inv_tt, training_eigen = stable_inverse(c_tt, eigenvalue_cutoff)
    ones_t = [1.0] * len(training)
    raw_weights = _matvec(inv_tt, ones_t)
    denominator = math.fsum(raw_weights)
    if denominator <= 0:
        raise ValueError("constant-fit information must be positive")
    weights = [value / denominator for value in raw_weights]
    amplitude = math.fsum(
        weight * values[index]
        for weight, index in zip(weights, training)
    )
    amplitude_variance = 1.0 / denominator

    c_hh = _subset(covariance, heldout)
    c_ht = _subset(covariance, heldout, training)
    covariance_with_amplitude = _matvec(c_ht, weights)
    residual_covariance = [
        [
            c_hh[i][j]
            + amplitude_variance
            - covariance_with_amplitude[i]
            - covariance_with_amplitude[j]
            for j in range(len(heldout))
        ]
        for i in range(len(heldout))
    ]
    residuals = [values[index] - amplitude for index in heldout]
    inv_resid, residual_eigen = stable_inverse(
        residual_covariance, eigenvalue_cutoff
    )
    chi_square = _quadratic(residuals, inv_resid)
    gauss_jordan_crosscheck = None
    try:
        gauss_jordan_crosscheck = _quadratic(
            residuals, _inverse(residual_covariance)
        )
    except (ArithmeticError, ValueError):
        gauss_jordan_crosscheck = None

    payload: Dict[str, object] = {
        "training_sizes": list(training_sizes),
        "heldout_sizes": list(heldout_sizes),
        "amplitude": amplitude,
        "amplitude_se": math.sqrt(amplitude_variance),
        "training_weights": weights,
        "heldout_observed": [values[index] for index in heldout],
        "heldout_residuals": residuals,
        "heldout_residual_covariance": residual_covariance,
        "heldout_chi_square": chi_square,
        "heldout_chi_square_kind": "plugin_asymptotic",
        "heldout_dof": len(heldout),
        "score_solver": residual_eigen.get("solver"),
        "training_covariance_eigenstructure": _public_eigenstructure(training_eigen),
        "heldout_residual_eigenstructure": _public_eigenstructure(residual_eigen),
        "eigenvalue_truncation_sensitivity": truncation_sensitivity(
            residuals, residual_covariance
        ),
        "gauss_jordan_crosscheck_chi_square": gauss_jordan_crosscheck,
    }
    if batch_count is not None:
        payload["finite_batch_calibration"] = {
            "batch_count": batch_count,
            "plugin_chi_square": chi_square,
            "plugin_chi_square_kind": "asymptotic_known_covariance",
            "plugin_chi_square_survival": chi_square_survival(chi_square, len(heldout)),
            "hotelling": hotelling_calibration(chi_square, len(heldout), batch_count),
        }
    return payload


def bootstrap_heldout_calibration(
    batch_vectors: Sequence[Sequence[float]],
    sizes: Sequence[int],
    training_sizes: Sequence[int],
    heldout_sizes: Sequence[int],
    observed_chi_square: float,
    replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
    seed: int = DEFAULT_BOOTSTRAP_SEED,
    eigenvalue_cutoff: float = DEFAULT_EIGENVALUE_CUTOFF,
    diagonal_only: bool = False,
) -> Dict[str, object]:
    if replicates < 10:
        raise ValueError("bootstrap replicates must be at least 10")
    rng = random.Random(seed)
    batches = len(batch_vectors)
    scores: List[float] = []
    failures = 0
    for _ in range(replicates):
        sample = [batch_vectors[rng.randrange(batches)] for _ in range(batches)]
        means, covariance = covariance_of_mean(sample)
        if diagonal_only:
            covariance = _diagonal(covariance)
        try:
            result = constant_heldout_audit(
                means,
                covariance,
                sizes,
                training_sizes,
                heldout_sizes,
                batch_count=batches,
                eigenvalue_cutoff=eigenvalue_cutoff,
            )
        except (ArithmeticError, ValueError):
            failures += 1
            continue
        scores.append(float(result["heldout_chi_square"]))
    if len(scores) < 10:
        raise ArithmeticError("bootstrap produced too few finite held-out scores")
    scores.sort()

    def quantile(level: float) -> float:
        if not 0.0 <= level <= 1.0:
            raise ValueError("quantile level must lie in [0, 1]")
        position = level * (len(scores) - 1)
        left = int(math.floor(position))
        right = int(math.ceil(position))
        if left == right:
            return scores[left]
        weight = position - left
        return (1.0 - weight) * scores[left] + weight * scores[right]

    return {
        "replicates_requested": replicates,
        "replicates_used": len(scores),
        "failures": failures,
        "seed": seed,
        "mean": math.fsum(scores) / len(scores),
        "quantiles": {
            "0.5": quantile(0.5),
            "0.9": quantile(0.9),
            "0.95": quantile(0.95),
            "0.99": quantile(0.99),
        },
        "observed_plugin_chi_square": observed_chi_square,
        "observed_survival": sum(
            1 for value in scores if value >= observed_chi_square
        )
        / len(scores),
        "diagonal_only": diagonal_only,
        "note": (
            "uncentered batch bootstrap of the held-out plug-in chi-square, "
            "measuring sampling variability around the observed sample; "
            "not a null p-value and not paper-facing"
        ),
    }


def sidecar_metadata_path(histograms: Path) -> Path:
    name = histograms.name
    if name.endswith(".hist.csv"):
        return histograms.with_name(name[: -len(".hist.csv")] + ".metadata.json")
    if name.endswith(".csv"):
        return histograms.with_suffix(".metadata.json")
    return histograms.with_name(histograms.name + ".metadata.json")


def validate_coupling_metadata(
    metadata: Mapping[str, object],
    records: Mapping[RecordKey, Record],
    seed_label: str,
    expected_seed: Optional[int] = None,
) -> Dict[str, object]:
    """Prove that aligned batch ids share one RNG schema, seed, and counters.

    Equal batch ids are not themselves evidence of cross-size coupling.
    """

    if not str(seed_label).strip():
        raise ValueError("seed-label is required to name a validated coupling contract")
    required = (
        "engine",
        "rng",
        "seed",
        "replica_counter_first",
        "replica_counter_last_exclusive",
        "batches",
        "samples_per_pair",
    )
    missing = [key for key in required if key not in metadata]
    if missing:
        raise ValueError("coupling metadata missing {}".format(missing))
    rng = str(metadata["rng"]).strip()
    engine = str(metadata["engine"]).strip()
    if not rng or not engine:
        raise ValueError("metadata rng schema and engine must be nonempty")
    seed = int(metadata["seed"])
    counter_first = int(metadata["replica_counter_first"])
    counter_last = int(metadata["replica_counter_last_exclusive"])
    batches = int(metadata["batches"])
    samples_per_pair = int(metadata["samples_per_pair"])
    if expected_seed is not None and seed != expected_seed:
        raise ValueError(
            "metadata seed {} disagrees with declared seed {}".format(
                seed, expected_seed
            )
        )
    if batches < 2:
        raise ValueError("metadata batches must be at least two")
    if samples_per_pair <= 0 or samples_per_pair % batches != 0:
        raise ValueError("samples_per_pair must be positive and divisible by batches")
    samples_per_batch = samples_per_pair // batches
    if counter_last != counter_first + samples_per_pair:
        raise ValueError(
            "replica counter interval must equal samples_per_pair; "
            "sizes cannot share a declared coupling stream"
        )
    if counter_last <= counter_first:
        raise ValueError("replica counter range must be nonempty")

    designs = metadata.get("designs")
    if not isinstance(designs, list) or not designs:
        raise ValueError("metadata designs are required to validate cross-size coupling")
    design_sizes = []
    for index, design in enumerate(designs):
        if not isinstance(design, dict) or "N" not in design:
            raise ValueError("metadata design {} is invalid".format(index))
        if any(key in design for key in ("seed", "rng", "replica_counter_first")):
            raise ValueError(
                "per-size RNG overrides are present; refuse to treat batch ids as coupled"
            )
        design_sizes.append(int(design["N"]))
    if len(set(design_sizes)) != len(design_sizes):
        raise ValueError("metadata designs contain duplicate sizes")

    sizes, hist_batches, _grouped = _orientation_batches(records)
    if set(sizes) - set(design_sizes):
        raise ValueError("histogram sizes are not covered by metadata designs")
    if len(hist_batches) != batches:
        raise ValueError("histogram batch count disagrees with metadata")
    histogram_samples = {
        int(record["samples"]) for record in records.values()
    }
    if histogram_samples != {samples_per_batch}:
        raise ValueError("histogram samples disagree with metadata samples_per_batch")

    per_batch_ranges = [
        {
            "batch": batch,
            "replica_counter_first": counter_first + batch * samples_per_batch,
            "replica_counter_last_exclusive": counter_first
            + (batch + 1) * samples_per_batch,
        }
        for batch in hist_batches
    ]
    return {
        "seed_label": seed_label,
        "validated": True,
        "engine": engine,
        "rng": rng,
        "seed": seed,
        "replica_counter_first": counter_first,
        "replica_counter_last_exclusive": counter_last,
        "batches": batches,
        "samples_per_pair": samples_per_pair,
        "samples_per_batch": samples_per_batch,
        "sizes": sizes,
        "design_sizes": design_sizes,
        "coupling": str(metadata.get("coupling", "")),
        "shared_rng_schema": True,
        "shared_seed": True,
        "shared_counter_ranges": True,
        "per_batch_counter_range_formula": (
            "replica_counter_first + batch * samples_per_batch"
        ),
        "per_batch_counter_ranges": per_batch_ranges,
    }


def _orientation_batches(
    records: Mapping[RecordKey, Record],
) -> Tuple[List[int], List[int], Grouped]:
    sizes = sorted({key[0] for key in records})
    by_size: Grouped = {}
    common_batches: Optional[List[int]] = None
    sample_counts = set()
    for n in sizes:
        by_size[n] = {}
        for orientation in ("first", "second"):
            selected = {
                key[2]: records[key]
                for key in records
                if key[0] == n and key[1] == orientation
            }
            if not selected:
                raise ValueError("N={} has no {} batches".format(n, orientation))
            ids = sorted(selected)
            if ids != list(range(len(ids))):
                raise ValueError(
                    "N={} {} batch ids are not contiguous".format(n, orientation)
                )
            by_size[n][orientation] = selected
            if common_batches is None:
                common_batches = ids
            elif ids != common_batches:
                raise ValueError("all sizes/orientations must share aligned batch ids")
        for batch in common_batches or []:
            first_samples = int(by_size[n]["first"][batch]["samples"])
            second_samples = int(by_size[n]["second"][batch]["samples"])
            if first_samples != second_samples:
                raise ValueError("paired orientations must have equal batch samples")
            sample_counts.add(first_samples)
    if common_batches is None:
        raise ValueError("no batches found")
    if len(sample_counts) != 1:
        raise ValueError(
            "all aligned sizes and batches must share one sample count; "
            "covariance_of_mean is an unweighted equal-batch estimator"
        )
    return sizes, common_batches, by_size


def _root_from_records(n: int, selected: Sequence[Record]) -> mp.mpf:
    if not selected:
        raise ValueError("cannot reconstruct a root from zero batches")
    minus = add_histograms(selected, "minus")
    plus = add_histograms(selected, "plus")
    samples = sum(int(row["samples"]) for row in selected)
    return matching_root(n, samples, minus, plus)


def _root_gap_pseudovalues(
    sizes: Sequence[int],
    batches: Sequence[int],
    by_size: Grouped,
) -> Tuple[Dict[Tuple[int, int], float], Dict[str, object]]:
    pseudo_by_size_batch: Dict[Tuple[int, int], float] = {}
    details: Dict[str, object] = {}
    for n in sizes:
        first_all = [by_size[n]["first"][batch] for batch in batches]
        second_all = [by_size[n]["second"][batch] for batch in batches]
        full_gap = float(
            _root_from_records(n, first_all) - _root_from_records(n, second_all)
        )
        delete_one: List[float] = []
        for omitted in batches:
            first_reduced = [
                by_size[n]["first"][batch]
                for batch in batches
                if batch != omitted
            ]
            second_reduced = [
                by_size[n]["second"][batch]
                for batch in batches
                if batch != omitted
            ]
            delete_one.append(
                float(
                    _root_from_records(n, first_reduced)
                    - _root_from_records(n, second_reduced)
                )
            )
        pseudovalues = jackknife_pseudovalues(full_gap, delete_one)
        for batch, value in zip(batches, pseudovalues):
            pseudo_by_size_batch[(n, batch)] = value
        bias_corrected = math.fsum(pseudovalues) / len(pseudovalues)
        details[str(n)] = {
            "full_estimate": full_gap,
            "bias_corrected_estimate": bias_corrected,
            "bias_correction": bias_corrected - full_gap,
            "delete_one_count": len(delete_one),
        }
    return pseudo_by_size_batch, details


def reconstruct_batch_metrics(
    records: Mapping[RecordKey, Record],
    p: mp.mpf,
) -> Tuple[List[int], List[Record], Dict[str, object]]:
    sizes, batches, by_size = _orientation_batches(records)
    root_pseudovalues, root_details = _root_gap_pseudovalues(
        sizes, batches, by_size
    )
    output: List[Record] = []
    for batch in batches:
        for n in sizes:
            first = by_size[n]["first"][batch]
            second = by_size[n]["second"][batch]
            first_m, first_d = evaluate_histogram(first, p)
            second_m, second_d = evaluate_histogram(second, p)
            delta_m = first_m - second_m
            mean_slope = (first_d + second_d) / 2
            root_gap_pseudovalue = root_pseudovalues[(n, batch)]
            delta_cos4 = cos4(
                int(first["a"]), int(first["b"])
            ) - cos4(int(second["a"]), int(second["b"]))
            if delta_cos4 == 0:
                raise ValueError("zero angular leverage in a batch")
            output.append(
                {
                    "N": n,
                    "batch": batch,
                    "samples": int(first["samples"]),
                    "a1": int(first["a"]),
                    "b1": int(first["b"]),
                    "a2": int(second["a"]),
                    "b2": int(second["b"]),
                    "delta_cos4": delta_cos4,
                    "delta_M": float(delta_m),
                    "mean_M_prime": float(mean_slope),
                    "root_gap_jackknife_pseudovalue": root_gap_pseudovalue,
                    "A_M": float(
                        n ** (13.0 / 8.0) * delta_m / delta_cos4
                    ),
                    "B": float(n ** (-3.0 / 8.0) * mean_slope),
                    "A_p_jackknife_pseudovalue": float(
                        -n * n * root_gap_pseudovalue / delta_cos4
                    ),
                }
            )
    return sizes, output, root_details


def _matrix_by_field(
    rows: Sequence[Record],
    sizes: Sequence[int],
    field: str,
) -> List[List[float]]:
    by_batch: Dict[int, Dict[int, float]] = {}
    for row in rows:
        by_batch.setdefault(int(row["batch"]), {})[int(row["N"])] = float(
            row[field]
        )
    output = []
    for batch in sorted(by_batch):
        if set(by_batch[batch]) != set(sizes):
            raise ValueError(
                "batch {} is incomplete for field {}".format(batch, field)
            )
        output.append([by_batch[batch][n] for n in sizes])
    return output


def _diagonal(covariance: Matrix) -> Matrix:
    return [
        [covariance[i][i] if i == j else 0.0 for j in range(len(covariance))]
        for i in range(len(covariance))
    ]


def audit(
    records: Mapping[RecordKey, Record],
    p: mp.mpf,
    training_sizes: Sequence[int],
    heldout_sizes: Sequence[int],
    coupling_contract: Optional[Mapping[str, object]] = None,
    bootstrap_replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
    bootstrap_seed: int = DEFAULT_BOOTSTRAP_SEED,
    eigenvalue_cutoff: float = DEFAULT_EIGENVALUE_CUTOFF,
) -> Tuple[List[Record], Dict[str, object]]:
    sizes, batch_rows, root_details = reconstruct_batch_metrics(records, p)
    sample_counts = {int(row["samples"]) for row in batch_rows}
    if len(sample_counts) != 1:
        raise ValueError(
            "all aligned sizes and batches must share one sample count; "
            "covariance_of_mean is an unweighted equal-batch estimator"
        )
    samples_per_batch = next(iter(sample_counts))
    batch_count = len({int(row["batch"]) for row in batch_rows})
    metric_fields = {
        "delta_M": "delta_M",
        "root_gap": "root_gap_jackknife_pseudovalue",
        "mean_M_prime": "mean_M_prime",
        "A_M": "A_M",
        "B": "B",
        "A_p": "A_p_jackknife_pseudovalue",
    }
    payload: Dict[str, object] = {
        "format_version": FORMAT_VERSION,
        "p": mp.nstr(p, mp.mp.dps),
        "sizes": sizes,
        "batch_count": batch_count,
        "equal_batch_weight_contract": {
            "estimator": "unweighted_equal_batch_mean",
            "asserted": True,
            "samples_per_batch": samples_per_batch,
        },
        "coupling_contract": dict(coupling_contract) if coupling_contract else None,
        "nonlinear_estimator": {
            "root_gap_method": "delete_one_jackknife_pseudovalues",
            "by_N": root_details,
        },
        "metrics": {},
        "constant_amplitude_audits": {},
        "score_numerics": {
            "inverse": "cholesky_then_svd_pseudoinverse",
            "eigenvalue_relative_cutoff": eigenvalue_cutoff,
            "heldout_chi_square_kind": "plugin_asymptotic",
        },
    }
    for metric, field in metric_fields.items():
        matrix = _matrix_by_field(batch_rows, sizes, field)
        means, covariance = covariance_of_mean(matrix)
        eigen = covariance_eigenstructure(covariance, eigenvalue_cutoff)
        metric_payload = {
            "batch_field": field,
            "means": dict(zip(map(str, sizes), means)),
            "standard_errors": dict(
                zip(
                    map(str, sizes),
                    [math.sqrt(covariance[i][i]) for i in range(len(sizes))],
                )
            ),
            "covariance_of_means": covariance,
            "correlation_of_batch_values": correlations(covariance),
            "eigenstructure": _public_eigenstructure(eigen),
        }
        metrics_payload = payload["metrics"]
        assert isinstance(metrics_payload, dict)
        metrics_payload[metric] = metric_payload
        if metric in ("A_M", "A_p"):
            full = constant_heldout_audit(
                means,
                covariance,
                sizes,
                training_sizes,
                heldout_sizes,
                batch_count=batch_count,
                eigenvalue_cutoff=eigenvalue_cutoff,
            )
            diagonal = constant_heldout_audit(
                means,
                _diagonal(covariance),
                sizes,
                training_sizes,
                heldout_sizes,
                batch_count=batch_count,
                eigenvalue_cutoff=eigenvalue_cutoff,
            )
            full_calibration = full.get("finite_batch_calibration")
            if isinstance(full_calibration, dict):
                full_calibration["bootstrap"] = bootstrap_heldout_calibration(
                    matrix,
                    sizes,
                    training_sizes,
                    heldout_sizes,
                    float(full["heldout_chi_square"]),
                    replicates=bootstrap_replicates,
                    seed=bootstrap_seed,
                    eigenvalue_cutoff=eigenvalue_cutoff,
                )
            diagonal_calibration = diagonal.get("finite_batch_calibration")
            if isinstance(diagonal_calibration, dict):
                diagonal_calibration["bootstrap"] = bootstrap_heldout_calibration(
                    matrix,
                    sizes,
                    training_sizes,
                    heldout_sizes,
                    float(diagonal["heldout_chi_square"]),
                    replicates=bootstrap_replicates,
                    seed=bootstrap_seed,
                    eigenvalue_cutoff=eigenvalue_cutoff,
                    diagonal_only=True,
                )
            audits = payload["constant_amplitude_audits"]
            assert isinstance(audits, dict)
            audits[metric] = {
                "full_covariance": full,
                "diagonal_covariance": diagonal,
            }
    return batch_rows, payload


def write_batch_csv(path: Path, rows: Sequence[Record]) -> None:
    fields = [
        "N",
        "batch",
        "samples",
        "a1",
        "b1",
        "a2",
        "b2",
        "delta_cos4",
        "delta_M",
        "mean_M_prime",
        "root_gap_jackknife_pseudovalue",
        "A_M",
        "B",
        "A_p_jackknife_pseudovalue",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_covariance_csv(path: Path, payload: Mapping[str, object]) -> None:
    sizes = [int(value) for value in payload["sizes"]]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "metric",
                "N_i",
                "N_j",
                "covariance_of_means",
                "correlation",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        metrics = payload["metrics"]
        assert isinstance(metrics, dict)
        for metric, raw in metrics.items():
            assert isinstance(raw, dict)
            covariance = raw["covariance_of_means"]
            correlation = raw["correlation_of_batch_values"]
            for i, n_i in enumerate(sizes):
                for j in range(i, len(sizes)):
                    writer.writerow(
                        {
                            "metric": metric,
                            "N_i": n_i,
                            "N_j": sizes[j],
                            "covariance_of_means": covariance[i][j],
                            "correlation": correlation[i][j],
                        }
                    )


def write_challenge_inputs(
    output_dir: Path,
    batch_rows: Sequence[Record],
    payload: Mapping[str, object],
    seed_label: str,
) -> None:
    sizes = [int(value) for value in payload["sizes"]]
    representatives: Dict[int, Record] = {}
    for row in batch_rows:
        representatives.setdefault(int(row["N"]), row)
    metrics = payload["metrics"]
    assert isinstance(metrics, dict)
    delta = metrics["delta_M"]
    means = delta["means"]
    standard_errors = delta["standard_errors"]
    covariance = delta["covariance_of_means"]

    observation_path = output_dir / "delta_M_observations.csv"
    with observation_path.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "row_id",
            "N",
            "seed",
            "delta_M",
            "delta_M_se",
            "delta_cos4",
            "a1",
            "b1",
            "a2",
            "b2",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for n in sizes:
            row = representatives[n]
            writer.writerow(
                {
                    "row_id": "{}:{}".format(n, seed_label),
                    "N": n,
                    "seed": seed_label,
                    "delta_M": means[str(n)],
                    "delta_M_se": standard_errors[str(n)],
                    "delta_cos4": row["delta_cos4"],
                    "a1": row["a1"],
                    "b1": row["b1"],
                    "a2": row["a2"],
                    "b2": row["b2"],
                }
            )

    covariance_path = output_dir / "delta_M_covariance.csv"
    with covariance_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["row_id_i", "row_id_j", "covariance"],
            lineterminator="\n",
        )
        writer.writeheader()
        for i, n_i in enumerate(sizes):
            for j in range(i, len(sizes)):
                writer.writerow(
                    {
                        "row_id_i": "{}:{}".format(n_i, seed_label),
                        "row_id_j": "{}:{}".format(sizes[j], seed_label),
                        "covariance": covariance[i][j],
                    }
                )


def parse_sizes(text: str) -> Tuple[int, ...]:
    try:
        values = tuple(int(value) for value in text.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "sizes must be comma-separated integers"
        ) from exc
    if (
        not values
        or any(value <= 0 for value in values)
        or len(set(values)) != len(values)
    ):
        raise argparse.ArgumentTypeError("sizes must be unique positive integers")
    return values


def _json_ready(value: object) -> object:
    if isinstance(value, float):
        if math.isfinite(value):
            return value
        if math.isinf(value):
            return "Infinity" if value > 0 else "-Infinity"
        return "NaN"
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    return value


def load_metadata(path: Path) -> Dict[str, object]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("coupling metadata must be a JSON object")
    return raw


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histograms", type=Path, required=True)
    parser.add_argument(
        "--metadata",
        type=Path,
        help=(
            "run metadata proving shared RNG schema, seed, and replica-counter "
            "ranges; default is the histogram sidecar .metadata.json"
        ),
    )
    parser.add_argument("--p", default="0.592746050790")
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument(
        "--training-sizes", type=parse_sizes, default=(65, 85, 130)
    )
    parser.add_argument(
        "--heldout-sizes", type=parse_sizes, default=(145, 170)
    )
    parser.add_argument("--seed-label", default="threshold-rank-coupled")
    parser.add_argument(
        "--expected-seed",
        type=int,
        help="if set, metadata seed must equal this value",
    )
    parser.add_argument(
        "--bootstrap-replicates",
        type=int,
        default=DEFAULT_BOOTSTRAP_REPLICATES,
    )
    parser.add_argument(
        "--bootstrap-seed",
        type=int,
        default=DEFAULT_BOOTSTRAP_SEED,
    )
    parser.add_argument(
        "--eigenvalue-cutoff",
        type=float,
        default=DEFAULT_EIGENVALUE_CUTOFF,
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.dps < 30:
        raise SystemExit("--dps must be at least 30")
    if set(args.training_sizes) & set(args.heldout_sizes):
        raise SystemExit("training and held-out sizes must be disjoint")
    if args.bootstrap_replicates < 10:
        raise SystemExit("--bootstrap-replicates must be at least 10")
    if not math.isfinite(args.eigenvalue_cutoff) or args.eigenvalue_cutoff < 0:
        raise SystemExit("--eigenvalue-cutoff must be a finite nonnegative number")
    mp.mp.dps = args.dps
    p = mp.mpf(args.p)
    if not 0 < p < 1:
        raise SystemExit("--p must lie strictly between zero and one")

    records = read_histograms(args.histograms)
    metadata_path = args.metadata or sidecar_metadata_path(args.histograms)
    if not metadata_path.is_file():
        raise SystemExit(
            "refusing to treat aligned batch ids as coupled without metadata "
            "that proves shared RNG schema, seed, and counter ranges: {}".format(
                metadata_path
            )
        )
    try:
        coupling_contract = validate_coupling_metadata(
            load_metadata(metadata_path),
            records,
            args.seed_label,
            expected_seed=args.expected_seed,
        )
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    coupling_contract["metadata_path"] = str(metadata_path)

    batch_rows, payload = audit(
        records,
        p,
        args.training_sizes,
        args.heldout_sizes,
        coupling_contract=coupling_contract,
        bootstrap_replicates=args.bootstrap_replicates,
        bootstrap_seed=args.bootstrap_seed,
        eigenvalue_cutoff=args.eigenvalue_cutoff,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_batch_csv(args.output_dir / "batch_metrics.csv", batch_rows)
    write_covariance_csv(
        args.output_dir / "cross_size_covariance.csv", payload
    )
    write_challenge_inputs(
        args.output_dir, batch_rows, payload, args.seed_label
    )
    (args.output_dir / "summary.json").write_text(
        json.dumps(_json_ready(payload), indent=2) + "\n", encoding="utf-8"
    )

    audits = payload["constant_amplitude_audits"]
    print(
        json.dumps(
            _json_ready(
                {
                    "sizes": payload["sizes"],
                    "batch_count": payload["batch_count"],
                    "root_gap_method": payload["nonlinear_estimator"][
                        "root_gap_method"
                    ],
                    "equal_batch_weight_contract": payload[
                        "equal_batch_weight_contract"
                    ],
                    "coupling_validated": True,
                    "A_M_full_heldout_chi_square": audits["A_M"][
                        "full_covariance"
                    ]["heldout_chi_square"],
                    "A_M_diagonal_heldout_chi_square": audits["A_M"][
                        "diagonal_covariance"
                    ]["heldout_chi_square"],
                    "A_p_full_heldout_chi_square": audits["A_p"][
                        "full_covariance"
                    ]["heldout_chi_square"],
                    "A_p_diagonal_heldout_chi_square": audits["A_p"][
                        "diagonal_covariance"
                    ]["heldout_chi_square"],
                    "A_M_full_score_solver": audits["A_M"]["full_covariance"][
                        "score_solver"
                    ],
                    "A_p_full_score_solver": audits["A_p"]["full_covariance"][
                        "score_solver"
                    ],
                }
            ),
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
