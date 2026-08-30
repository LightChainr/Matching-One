#!/usr/bin/env python3
"""Model-free multivariate Hankel rank lower bound for P250/P249/P255."""

from __future__ import annotations

import argparse
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from scipy.stats import chi2, f

from score_z5_projective_leg_bivariate_state import means, pair, read_batches
from score_z5_projective_leg_cross_scale import jackknife_covariance
from score_z5_projective_leg_pair_transfer import CHANNELS


MONOMIALS_1 = ((0, 0), (1, 0), (0, 1))
MONOMIALS_2 = ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2))
MATRIX_GROUPS = {
    "plus_charge1": (("plus", 1),),
    "plus_charge2": (("plus", 2),),
    "minus_charge1": (("minus", 1),),
    "minus_charge2": (("minus", 2),),
    "plus_block": (("plus", 1), ("plus", 2)),
    "minus_block": (("minus", 1), ("minus", 2)),
    "shared_block": CHANNELS,
}
RANKS = (1, 2, 3, 4, 5)
EIGEN_RELATIVE_CUTOFF = 1e-10


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add(first: tuple[int, int], second: tuple[int, int]) -> tuple[int, int]:
    return first[0] + second[0], first[1] + second[1]


def hankel_matrix(
    values: Mapping[str, float], channels: Sequence[tuple[str, int]],
    monomials: Sequence[tuple[int, int]] = MONOMIALS_2,
) -> np.ndarray:
    return np.asarray([
        [pair(values, add(left, right), channel) for right in monomials]
        for channel in channels for left in monomials
    ], dtype=complex)


def row_labels(channels: Sequence[tuple[str, int]], monomials=MONOMIALS_2) -> list[dict]:
    return [
        {"hand": channel[0], "charge": channel[1], "monomial": list(monomial)}
        for channel in channels for monomial in monomials
    ]


def maximum_volume_pivot(matrix: np.ndarray, rank: int) -> dict:
    """Choose a stable deterministic coordinate chart for rank<=r."""
    best: tuple[float, tuple[int, ...], tuple[int, ...]] | None = None
    for columns in combinations(range(matrix.shape[1]), rank):
        for rows in combinations(range(matrix.shape[0]), rank):
            volume = float(abs(np.linalg.det(matrix[np.ix_(rows, columns)])))
            candidate = (volume, tuple(rows), tuple(columns))
            if best is None or candidate[0] > best[0]:
                best = candidate
    if best is None or not best[0] > 0.0:
        raise ValueError("no invertible pivot chart")
    pivot = matrix[np.ix_(best[1], best[2])]
    return {
        "rows": best[1],
        "columns": best[2],
        "abs_determinant": best[0],
        "condition_number": float(np.linalg.cond(pivot)),
    }


def schur_complement(matrix: np.ndarray, pivot: Mapping[str, object]) -> np.ndarray:
    rows = tuple(pivot["rows"])
    columns = tuple(pivot["columns"])
    other_rows = tuple(index for index in range(matrix.shape[0]) if index not in rows)
    other_columns = tuple(index for index in range(matrix.shape[1]) if index not in columns)
    p = matrix[np.ix_(rows, columns)]
    q = matrix[np.ix_(rows, other_columns)]
    r = matrix[np.ix_(other_rows, columns)]
    s = matrix[np.ix_(other_rows, other_columns)]
    return s - r @ np.linalg.solve(p, q)


def realify(matrix: np.ndarray) -> list[float]:
    output = []
    for value in matrix.ravel():
        output.extend((float(value.real), float(value.imag)))
    return output


def covariance_score(point: Sequence[float], deleted: Sequence[Sequence[float]]) -> dict:
    """Correlation-normalized pseudoinverse plus finite-batch Hotelling score."""
    covariance = np.asarray(jackknife_covariance(deleted), dtype=float)
    point_array = np.asarray(point, dtype=float)
    scales = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    if np.any(scales <= 0.0):
        raise ValueError("rank residual has a zero covariance diagonal")
    correlation = covariance / scales[:, None] / scales[None, :]
    correlation = 0.5 * (correlation + correlation.T)
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    largest = float(eigenvalues[-1])
    kept = eigenvalues > EIGEN_RELATIVE_CUTOFF * largest
    if not np.any(kept):
        raise ValueError("rank residual covariance has no resolved mode")
    normalized = point_array / scales
    projections = eigenvectors[:, kept].T @ normalized
    statistic = float(np.sum(projections**2 / eigenvalues[kept]))
    degrees = int(np.count_nonzero(kept))
    batches = len(deleted)
    if degrees >= batches:
        hotelling_f = math.inf
        hotelling_p = 0.0
        hotelling_denominator_df = batches - degrees
    else:
        hotelling_f = (batches - degrees) * statistic / (degrees * (batches - 1))
        hotelling_denominator_df = batches - degrees
        hotelling_p = float(f.sf(hotelling_f, degrees, hotelling_denominator_df))
    return {
        "residual": list(point),
        "covariance": covariance.tolist(),
        "eigen_relative_cutoff": EIGEN_RELATIVE_CUTOFF,
        "resolved_covariance_modes": degrees,
        "discarded_covariance_modes": len(point) - degrees,
        "correlation_eigenvalue_min_kept": float(eigenvalues[kept][0]),
        "correlation_eigenvalue_max": largest,
        "asymptotic_chi_square": statistic,
        "asymptotic_degrees_of_freedom": degrees,
        "asymptotic_survival_p": float(chi2.sf(statistic, degrees)),
        "finite_batch_hotelling_F": hotelling_f,
        "finite_batch_numerator_df": degrees,
        "finite_batch_denominator_df": hotelling_denominator_df,
        "finite_batch_survival_p": hotelling_p,
    }


def rank_score(
    matrix: np.ndarray, deleted_matrices: Sequence[np.ndarray], rank: int,
) -> dict:
    pivot = maximum_volume_pivot(matrix, rank)
    point = realify(schur_complement(matrix, pivot))
    deleted = [realify(schur_complement(row, pivot)) for row in deleted_matrices]
    score = covariance_score(point, deleted)
    other_rows = [index for index in range(matrix.shape[0]) if index not in pivot["rows"]]
    other_columns = [index for index in range(matrix.shape[1]) if index not in pivot["columns"]]
    return {
        "null": f"rank(H)<={rank}",
        "pivot": {
            "rows": list(pivot["rows"]),
            "columns": list(pivot["columns"]),
            "abs_determinant": pivot["abs_determinant"],
            "condition_number": pivot["condition_number"],
        },
        "Schur_shape": [len(other_rows), len(other_columns)],
        "vanishing_minor_chart": "S-R P^-1 Q; locally complete rank-r determinantal constraints",
        "score": score,
    }


def decomposition_gate(values: Mapping[str, float]) -> dict:
    by_sum: dict[tuple[int, int], list[tuple[tuple[int, int], tuple[int, int]]]] = {}
    for left in MONOMIALS_2:
        for right in MONOMIALS_2:
            by_sum.setdefault(add(left, right), []).append((left, right))
    maximum = 0.0
    for channel in CHANNELS:
        for target, decompositions in by_sum.items():
            rows = [pair(values, add(left, right), channel) for left, right in decompositions]
            maximum = max(maximum, max(abs(value - rows[0]) for value in rows))
    return {
        "monomial_sums": len(by_sum),
        "maximum_decomposition_difference": maximum,
        "passed": maximum == 0.0,
        "scope": "This is an exact endpoint-label/Hankel construction gate, not an ordered-path TxTy versus TyTx observation.",
    }


def singular_payload(matrix: np.ndarray) -> list[float]:
    values = np.linalg.svd(matrix, compute_uv=False)
    return [float(value) for value in values]


def score(batches: Sequence[dict], manifest: Mapping[str, object]) -> dict:
    path = Path(manifest["input_batches"])
    if sha256(path) != manifest["input_batches_sha256"]:
        raise ValueError("bivariate batch hash changed")
    values = means(batches)
    deleted_values = [means(batches, index) for index in range(len(batches))]
    alpha = float(manifest["decision_alpha"])
    groups = {}
    for name, channels in MATRIX_GROUPS.items():
        matrix = hankel_matrix(values, channels)
        deleted_matrices = [hankel_matrix(row, channels) for row in deleted_values]
        ranks = {str(rank): rank_score(matrix, deleted_matrices, rank) for rank in RANKS}
        rejected = [
            rank for rank in RANKS
            if ranks[str(rank)]["score"]["finite_batch_survival_p"] < alpha
        ]
        lower_bound = max(rejected, default=0) + 1
        order1_flat = ranks["3"]["score"]["finite_batch_survival_p"] >= alpha
        groups[name] = {
            "channels": [list(channel) for channel in channels],
            "shape": list(matrix.shape),
            "row_labels": row_labels(channels),
            "column_monomials": [list(point) for point in MONOMIALS_2],
            "singular_values_descriptive_only": singular_payload(matrix),
            "rank_nulls": ranks,
            "rank_lower_bound_at_alpha": lower_bound,
            "order1_flat_extension_rank_le_3_survives": order1_flat,
        }
    shared_lower_bound = groups["shared_block"]["rank_lower_bound_at_alpha"]
    if shared_lower_bound >= 6:
        decision = "shared_multivariate_Hankel_rank_at_least_6"
    else:
        decision = f"shared_multivariate_Hankel_rank_at_least_{shared_lower_bound}"
    return {
        "schema": "matching-one/z5-projective-leg-hankel-rank-score/v1",
        "status": "model_free_existing_data_reanalysis",
        "monomial_basis_degree_le_2": [list(point) for point in MONOMIALS_2],
        "Hankel_decomposition_gate": decomposition_gate(values),
        "groups": groups,
        "decision_alpha": alpha,
        "primary_p_value": "finite_batch_hotelling_survival_p",
        "decision": decision,
        "flat_extension_scope": {
            "order1": "rank H_degree2 = rank H_degree1 can only close at rank<=3 and is tested by the rank-3 null",
            "order2": "not identifiable: testing rank H_degree3 = rank H_degree2 requires moments through total degree six; this stream stops at degree four",
            "ordered_path": "not identifiable: the stream records G(a,b) by endpoint displacement and does not separately record xy and yx paths",
        },
        "claim_boundary": [
            "The determinantal rank constraints assume only a shared path-independent bivariate moment sequence; they do not assume diagonalizability and therefore include commuting Jordan realizations.",
            "A shared-block lower bound applies to a common transfer state with channel-specific left functionals/amplitudes.",
            "Failure of order-one flatness means more than three states are needed; it is not by itself evidence of path or context dependence.",
            "Current endpoint-only data cannot test TxTy against TyTx. Ordered-path rows or moments through degree six are the precise missing observables for those stronger claims.",
        ],
    }


def render(result: Mapping[str, object]) -> str:
    lines = [
        "# P250 model-free multivariate Hankel rank score",
        "",
        "| group | rank<=1 p | rank<=2 p | rank<=3 p | rank<=4 p | rank<=5 p | lower bound |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for name, row in result["groups"].items():
        probabilities = [row["rank_nulls"][str(rank)]["score"]["finite_batch_survival_p"] for rank in RANKS]
        lines.append(
            f"| {name} | " + " | ".join(f"{value:.6g}" for value in probabilities)
            + f" | {row['rank_lower_bound_at_alpha']} |"
        )
    lines += [
        "",
        f"Decision: `{result['decision']}`.",
        "",
        "The primary probabilities use the finite-400-batch Hotelling correction; full residual covariances and asymptotic chi-square scores are in JSON.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    result = score(read_batches(Path(manifest["input_batches"])), manifest)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
