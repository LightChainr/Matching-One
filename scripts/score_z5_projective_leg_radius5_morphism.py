#!/usr/bin/env python3
"""Frozen support-first scorer for the P250 radius-five morphism acquisition."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from scipy.stats import chi2, f

from score_z5_projective_leg_annihilator_bridge import (
    annihilator_line,
    bridge_score,
    projective_residual,
    transformed_basis,
)
from score_z5_projective_leg_bivariate_state import means as old_means, pair as old_pair, read_batches as read_old
from score_z5_projective_leg_cross_scale import jackknife_covariance
from score_z5_projective_leg_hankel_rank import MONOMIALS_2
from z5_projective_leg_bivariate_mc import label
from z5_projective_leg_radius5_morphism_mc import FIELD_ORDER, SCHEMA


ROWS_3 = ((3, 0), (2, 1), (1, 2), (0, 3))
CHANNELS = (("plus", 1), ("plus", 2), ("minus", 1), ("minus", 2))
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
    hand, charge = channel
    prefix = f"{label(*point)}_r{charge}_{hand}_"
    return complex(values[prefix + "re"], values[prefix + "im"])


def moment(
    old_values: Mapping[str, float], new_values: Mapping[str, float],
    point: tuple[int, int], channel: tuple[str, int],
) -> complex:
    degree = abs(point[0]) + abs(point[1])
    if degree == 5:
        return new_pair(new_values, point, channel)
    if degree <= 4:
        return old_pair(old_values, point, channel)
    raise ValueError(f"point {point} lies outside the frozen old+new domain")


def extension_matrix(
    old_values: Mapping[str, float], new_values: Mapping[str, float], hand: str,
    basis: Sequence[tuple[int, int]], rows3: Sequence[tuple[int, int]],
) -> np.ndarray:
    old_rows = np.asarray([
        [old_pair(old_values, (u[0] + v[0], u[1] + v[1]), (hand, charge)) for v in basis]
        for charge in (1, 2) for u in basis
    ], dtype=complex)
    new_rows = np.asarray([
        [moment(old_values, new_values, (u[0] + v[0], u[1] + v[1]), (hand, charge)) for v in basis]
        for charge in (1, 2) for u in rows3
    ], dtype=complex)
    return np.vstack((old_rows, new_rows))


def extension_residual(
    line: np.ndarray, old_values: Mapping[str, float], new_values: Mapping[str, float], hand: str,
    basis: Sequence[tuple[int, int]], rows3: Sequence[tuple[int, int]],
) -> list[float]:
    output = []
    for charge in (1, 2):
        for u in rows3:
            value = sum(
                line[index] * moment(old_values, new_values, (u[0] + v[0], u[1] + v[1]), (hand, charge))
                for index, v in enumerate(basis)
            )
            output.extend((float(value.real), float(value.imag)))
    return output


def covariance_sum(first: Sequence[Sequence[float]], second: Sequence[Sequence[float]]) -> np.ndarray:
    return np.asarray(first, dtype=float) + np.asarray(second, dtype=float)


def score_from_covariance(point: Sequence[float], covariance: np.ndarray, batches: int) -> dict:
    point_array = np.asarray(point, dtype=float)
    scales = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    if np.any(scales <= 0.0):
        raise ValueError("zero residual covariance scale")
    correlation = covariance / scales[:, None] / scales[None, :]
    correlation = 0.5 * (correlation + correlation.T)
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    keep = eigenvalues > EIGEN_RELATIVE_CUTOFF * eigenvalues[-1]
    normalized = point_array / scales
    projection = eigenvectors[:, keep].T @ normalized
    statistic = float(np.sum(projection**2 / eigenvalues[keep]))
    degrees = int(np.count_nonzero(keep))
    hotelling = (batches - degrees) * statistic / (degrees * (batches - 1))
    return {
        "residual": list(point),
        "covariance": covariance.tolist(),
        "resolved_covariance_modes": degrees,
        "asymptotic_chi_square": statistic,
        "asymptotic_survival_p": float(chi2.sf(statistic, degrees)),
        "finite_batch_Hotelling_F": hotelling,
        "finite_batch_survival_p": float(f.sf(hotelling, degrees, batches - degrees)),
    }


def independent_score(point, deleted_old, deleted_new) -> dict:
    covariance = covariance_sum(jackknife_covariance(deleted_old), jackknife_covariance(deleted_new))
    return score_from_covariance(point, covariance, min(len(deleted_old), len(deleted_new)))


def transformed_rows(alexander: bool, rotation_power: int) -> tuple[tuple[int, int], ...]:
    base = transformed_basis(alexander_reflection=alexander, rotation_power=rotation_power)
    transform = dict(zip(MONOMIALS_2, base))
    # The same linear D4 action is determined by images of x and y.
    ximage = transform[(1, 0)]
    yimage = transform[(0, 1)]
    return tuple((u[0] * ximage[0] + u[1] * yimage[0], u[0] * ximage[1] + u[1] * yimage[1]) for u in ROWS_3)


def score(old_rows: Sequence[dict], new_rows: Sequence[dict], manifest: Mapping[str, object]) -> dict:
    old_path = Path(manifest["old_batches"])
    new_path = Path(manifest["new_batches"])
    if sha256(old_path) != manifest["old_batches_sha256"] or sha256(new_path) != manifest["new_batches_sha256"]:
        raise ValueError("morphism batch hash changed")
    old_full = old_means(old_rows)
    new_full = new_means(new_rows)
    old_deleted = [old_means(old_rows, index) for index in range(len(old_rows))]
    new_deleted = [new_means(new_rows, index) for index in range(len(new_rows))]
    basis0 = transformed_basis(alexander_reflection=False, rotation_power=0)

    old_lines = {
        hand: annihilator_line(extension_matrix(old_full, new_full, hand, basis0, ()))[0]
        for hand in ("plus", "minus")
    }
    extension = {}
    for hand in ("plus", "minus"):
        point = extension_residual(old_lines[hand], old_full, new_full, hand, basis0, ROWS_3)
        deleted_from_old = [
            extension_residual(
                annihilator_line(extension_matrix(row, new_full, hand, basis0, ()))[0],
                row, new_full, hand, basis0, ROWS_3,
            ) for row in old_deleted
        ]
        deleted_from_new = [
            extension_residual(old_lines[hand], old_full, row, hand, basis0, ROWS_3)
            for row in new_deleted
        ]
        extension[hand] = independent_score(point, deleted_from_old, deleted_from_new)

    augmented = {}
    augmented_old_deleted = {}
    augmented_new_deleted = {}
    for hand in ("plus", "minus"):
        augmented[hand] = annihilator_line(extension_matrix(old_full, new_full, hand, basis0, ROWS_3))[0]
        augmented_old_deleted[hand] = [
            annihilator_line(extension_matrix(row, new_full, hand, basis0, ROWS_3))[0] for row in old_deleted
        ]
        augmented_new_deleted[hand] = [
            annihilator_line(extension_matrix(old_full, row, hand, basis0, ROWS_3))[0] for row in new_deleted
        ]

    maps = {}
    candidates = [("identity_conjugation", False, 0)] + [
        (f"Alexander_R{power}_conjugation", True, power) for power in range(4)
    ]
    for name, alexander, power in candidates:
        basis = transformed_basis(alexander_reflection=alexander, rotation_power=power)
        rows3 = transformed_rows(alexander, power)
        minus = annihilator_line(extension_matrix(old_full, new_full, "minus", basis, rows3))[0]
        minus_old_deleted = [
            annihilator_line(extension_matrix(row, new_full, "minus", basis, rows3))[0] for row in old_deleted
        ]
        minus_new_deleted = [
            annihilator_line(extension_matrix(old_full, row, "minus", basis, rows3))[0] for row in new_deleted
        ]
        plus = augmented["plus"].conjugate()
        # Use the same projective pivot for the point and both independent jackknife sources.
        shared = np.minimum(np.abs(plus), np.abs(minus))
        pivot = int(np.argmax(shared))
        point = projective_residual(plus, minus, pivot)
        deleted_from_old = [
            projective_residual(first.conjugate(), second, pivot)
            for first, second in zip(augmented_old_deleted["plus"], minus_old_deleted)
        ]
        deleted_from_new = [
            projective_residual(first.conjugate(), second, pivot)
            for first, second in zip(augmented_new_deleted["plus"], minus_new_deleted)
        ]
        maps[name] = {
            "alexander_reflection": alexander,
            "rotation_power": power,
            "comparison_pivot": pivot,
            "score": independent_score(point, deleted_from_old, deleted_from_new),
        }

    alpha = float(manifest["decision_alpha"])
    support = all(row["finite_batch_survival_p"] >= alpha for row in extension.values())
    survivors = [name for name, row in maps.items() if row["score"]["finite_batch_survival_p"] >= alpha]
    if not support:
        decision = "five_plus_five_direct_sum_extension_fails"
    elif survivors == ["identity_conjugation"]:
        decision = "identity_conjugation_identified"
    elif "identity_conjugation" not in survivors and any(name.startswith("Alexander") for name in survivors):
        decision = "Alexander_family_selected_over_identity"
    elif not survivors:
        decision = "general_five_plus_five_bridge_only"
    else:
        decision = "parameter_free_bridge_remains_nonidentified"
    return {
        "schema": "matching-one/z5-projective-leg-radius5-morphism-score/v1",
        "status": "fresh_radius5_support_first_reveal",
        "general_five_plus_five_extension_gate": extension,
        "parameter_free_maps": maps if support else "LOCKED_SUPPORT_FAILED",
        "surviving_maps": survivors if support else [],
        "decision_alpha": alpha,
        "decision": decision,
        "claim_boundary": [
            "The fresh shell first tests whether the two old annihilators extend; map scores are interpreted only after both support gates pass.",
            "A general 5+5 direct sum imposes no cross-hand line equality and is represented by the two independent extension gates.",
            "The old and fresh streams use independent seeds; their delete-one covariance contributions are added before scoring.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    result = score(read_old(Path(manifest["old_batches"])), read_new(Path(manifest["new_batches"])), manifest)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    lines = ["# P250 radius-five morphism score", "", f"Decision: `{result['decision']}`.", ""]
    args.markdown.write_text("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
