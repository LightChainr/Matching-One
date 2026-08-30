#!/usr/bin/env python3
"""Frozen rank ladder and R2 ideal-bridge score for the P250 degree-six shell."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from scipy.stats import chi2, f

from score_z5_projective_leg_bivariate_state import means as old_means, pair as old_pair, read_batches as read_old
from score_z5_projective_leg_cross_scale import jackknife_covariance
from score_z5_projective_leg_radius5_morphism import new_means as radius5_means, new_pair as radius5_pair, read_new as read_radius5
from z5_projective_leg_radius6_flat_mc import FIELD_ORDER, SCHEMA, selected_r2_alexander


MONOMIALS_2 = ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2))
DEGREE3 = ((3, 0), (2, 1), (1, 2), (0, 3))
MONOMIALS_3 = MONOMIALS_2 + DEGREE3
RANKS = (5, 6, 7, 8, 9)
EIGEN_RELATIVE_CUTOFF = 1e-10


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_new(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return [
            {
                key: int(value) if key in {"batch", "replica_first", "samples"}
                else value if key in {"field_sha256", "translation_sha256"}
                else float(value)
                for key, value in row.items()
            }
            for row in csv.DictReader(handle)
        ]


def new_means(rows: Sequence[dict], excluded: int | None = None) -> dict[str, float]:
    kept = [row for index, row in enumerate(rows) if index != excluded]
    samples = sum(row["samples"] for row in kept)
    return {key: sum(row[key] for row in kept) / samples for key in FIELD_ORDER}


def new_pair(values: Mapping[str, float], point: tuple[int, int], channel: tuple[str, int]) -> complex:
    from z5_projective_leg_bivariate_mc import label

    hand, charge = channel
    prefix = f"{label(*point)}_r{charge}_{hand}_"
    return complex(values[prefix + "re"], values[prefix + "im"])


def physical(point: tuple[int, int], hand: str) -> tuple[int, int]:
    return point if hand == "plus" else selected_r2_alexander(point)


def moment(
    old4: Mapping[str, float], old5: Mapping[str, float], fresh6: Mapping[str, float],
    abstract_point: tuple[int, int], hand: str, charge: int,
) -> complex:
    point = physical(abstract_point, hand)
    degree = abs(point[0]) + abs(point[1])
    if degree <= 4:
        return old_pair(old4, point, (hand, charge))
    if degree == 5:
        return radius5_pair(old5, point, (hand, charge))
    if degree == 6:
        return new_pair(fresh6, point, (hand, charge))
    raise ValueError(f"point {point} lies outside the degree-six frozen domain")


def hankel(
    old4: Mapping[str, float], old5: Mapping[str, float], fresh6: Mapping[str, float], hand: str,
) -> np.ndarray:
    return np.asarray([
        [
            moment(old4, old5, fresh6, (left[0] + right[0], left[1] + right[1]), hand, charge)
            for right in MONOMIALS_3
        ]
        for charge in (1, 2) for left in MONOMIALS_3
    ], dtype=complex)


def schur(matrix: np.ndarray, rows: Sequence[int], columns: Sequence[int]) -> np.ndarray:
    other_rows = tuple(index for index in range(matrix.shape[0]) if index not in rows)
    other_columns = tuple(index for index in range(matrix.shape[1]) if index not in columns)
    p = matrix[np.ix_(rows, columns)]
    q = matrix[np.ix_(rows, other_columns)]
    r = matrix[np.ix_(other_rows, columns)]
    s = matrix[np.ix_(other_rows, other_columns)]
    return s - r @ np.linalg.solve(p, q)


def realify(array: np.ndarray) -> list[float]:
    output = []
    for value in array.ravel():
        output.extend((float(value.real), float(value.imag)))
    return output


def covariance_sum(groups: Sequence[Sequence[Sequence[float]]]) -> np.ndarray:
    return sum((np.asarray(jackknife_covariance(group), dtype=float) for group in groups), start=0.0)


def score_vector(point: Sequence[float], deleted_groups: Sequence[Sequence[Sequence[float]]]) -> dict:
    covariance = covariance_sum(deleted_groups)
    point_array = np.asarray(point, dtype=float)
    scales = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    if np.any(scales <= 0.0):
        raise ValueError("flat-extension residual has zero covariance scale")
    correlation = covariance / scales[:, None] / scales[None, :]
    correlation = 0.5 * (correlation + correlation.T)
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    keep = eigenvalues > EIGEN_RELATIVE_CUTOFF * eigenvalues[-1]
    projection = eigenvectors[:, keep].T @ (point_array / scales)
    statistic = float(np.sum(projection**2 / eigenvalues[keep]))
    degrees = int(np.count_nonzero(keep))
    batches = min(len(group) for group in deleted_groups)
    if degrees >= batches:
        hotelling = float("inf")
        finite_p = 0.0
    else:
        hotelling = (batches - degrees) * statistic / (degrees * (batches - 1))
        finite_p = float(f.sf(hotelling, degrees, batches - degrees))
    return {
        "residual": list(point),
        "covariance": covariance.tolist(),
        "resolved_covariance_modes": degrees,
        "asymptotic_chi_square": statistic,
        "asymptotic_survival_p": float(chi2.sf(statistic, degrees)),
        "finite_batch_Hotelling_F": hotelling,
        "finite_batch_survival_p": finite_p,
    }


def fixed_rank_residual(matrix: np.ndarray, pivot: Mapping[str, Sequence[int]]) -> list[float]:
    return realify(schur(matrix, tuple(pivot["rows"]), tuple(pivot["columns"])))


def kernel_projector(matrix: np.ndarray, rank: int = 5) -> np.ndarray:
    _, _, right_adjoint = np.linalg.svd(matrix, full_matrices=False)
    kernel = right_adjoint[rank:].conjugate().T
    return kernel @ kernel.conjugate().T


def ideal_bridge_residual(plus: np.ndarray, minus: np.ndarray) -> list[float]:
    return realify(kernel_projector(plus).conjugate() - kernel_projector(minus))


def score(old4_rows: Sequence[dict], old5_rows: Sequence[dict], new6_rows: Sequence[dict], manifest: Mapping[str, object]) -> dict:
    for key, hash_key in (("old4_batches", "old4_sha256"), ("old5_batches", "old5_sha256"), ("new6_batches", "new6_sha256")):
        if sha256(Path(manifest[key])) != manifest[hash_key]:
            raise ValueError(f"{key} hash changed")

    old4 = old_means(old4_rows)
    old5 = radius5_means(old5_rows)
    new6 = new_means(new6_rows)
    old4_deleted = [old_means(old4_rows, index) for index in range(len(old4_rows))]
    old5_deleted = [radius5_means(old5_rows, index) for index in range(len(old5_rows))]
    new6_deleted = [new_means(new6_rows, index) for index in range(len(new6_rows))]

    full = {hand: hankel(old4, old5, new6, hand) for hand in ("plus", "minus")}
    deleted_by_source = {
        "old4": [
            {hand: hankel(row, old5, new6, hand) for hand in ("plus", "minus")}
            for row in old4_deleted
        ],
        "old5": [
            {hand: hankel(old4, row, new6, hand) for hand in ("plus", "minus")}
            for row in old5_deleted
        ],
        "new6": [
            {hand: hankel(old4, old5, row, hand) for hand in ("plus", "minus")}
            for row in new6_deleted
        ],
    }

    alpha = float(manifest["decision_alpha"])
    hands = {}
    for hand in ("plus", "minus"):
        ranks = {}
        for rank in RANKS:
            pivot = manifest["frozen_pivots"][hand][str(rank)]
            point = fixed_rank_residual(full[hand], pivot)
            deleted = [
                [fixed_rank_residual(row[hand], pivot) for row in deleted_by_source[source]]
                for source in ("old4", "old5", "new6")
            ]
            ranks[str(rank)] = {
                "null": f"rank(H3)<={rank}",
                "pivot": pivot,
                "Schur_shape": [20 - rank, 10 - rank],
                "score": score_vector(point, deleted),
            }
        compatible = [rank for rank in RANKS if ranks[str(rank)]["score"]["finite_batch_survival_p"] >= alpha]
        hands[hand] = {
            "matrix_shape": [20, 10],
            "rank_nulls": ranks,
            "rank_lower_bound_at_alpha": min(compatible, default=10),
            "rank5_flat_extension_survives": ranks["5"]["score"]["finite_batch_survival_p"] >= alpha,
        }

    flat = all(row["rank5_flat_extension_survives"] for row in hands.values())
    if flat:
        bridge_point = ideal_bridge_residual(full["plus"], full["minus"])
        bridge_deleted = [
            [ideal_bridge_residual(row["plus"], row["minus"]) for row in deleted_by_source[source]]
            for source in ("old4", "old5", "new6")
        ]
        bridge = score_vector(bridge_point, bridge_deleted)
        if bridge["finite_batch_survival_p"] >= alpha:
            decision = "rank5_flat_and_R2_conjugate_ideal_bridge_survives"
        else:
            decision = "rank5_flat_but_R2_conjugate_ideal_bridge_rejected"
    else:
        bridge = "LOCKED_RANK5_FLAT_EXTENSION_FAILED"
        decision = "rank5_flat_extension_rejected"

    return {
        "schema": "matching-one/z5-projective-leg-radius6-flat-score/v1",
        "status": "fresh_degree6_parameter_free_reveal",
        "abstract_monomials_degree_le_3": [list(point) for point in MONOMIALS_3],
        "fixed_gauge": "plus identity; minus Alexander R2; coefficient-space comparison uses conjugation",
        "hand_rank_ladders": hands,
        "R2_conjugate_kernel_projector_bridge": bridge,
        "decision_alpha": alpha,
        "decision": decision,
        "claim_boundary": [
            "The five-map family is not rescored; R2 plus conjugation is fixed before the degree-six stream.",
            "Rank-five flatness is a path-independent commuting truncated-moment claim through degree six, not an ordered-path TxTy versus TyTx measurement.",
            "If rank five fails, the frozen rank ladder localizes the minimum compatible H3 rank up to ten without naming a transfer algebra.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    result = score(
        read_old(Path(manifest["old4_batches"])),
        read_radius5(Path(manifest["old5_batches"])),
        read_new(Path(manifest["new6_batches"])),
        manifest,
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    lines = ["# P250 radius-six flat-extension score", "", f"Decision: `{result['decision']}`.", ""]
    args.markdown.write_text("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
