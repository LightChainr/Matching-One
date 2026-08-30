#!/usr/bin/env python3
"""Score the frozen P250 rank-eight R2-conjugate kernel-plane bridge."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

import score_z5_projective_leg_radius6_flat as base


SCHEMA = "matching-one/p250-rank8-r2-kernel-bridge-score/v1"
RANK = 8


def checked_inputs(manifest: Mapping[str, object]) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for name, entry in manifest["inputs"].items():
        path = Path(entry["path"])
        if base.sha256(path) != entry["sha256"]:
            raise ValueError(f"immutable input hash changed: {name}")
        paths[name] = path
    return paths


def check_upstream_gate(path: Path, alpha: float) -> dict[str, object]:
    upstream = json.loads(path.read_text())
    if "Alexander R2" not in upstream["fixed_gauge"]:
        raise ValueError("upstream score does not use the selected Alexander R2 gauge")
    summary: dict[str, object] = {}
    for hand in ("plus", "minus"):
        row = upstream["hand_rank_ladders"][hand]
        rank7_p = float(row["rank_nulls"]["7"]["score"]["finite_batch_survival_p"])
        rank8_p = float(row["rank_nulls"]["8"]["score"]["finite_batch_survival_p"])
        if int(row["rank_lower_bound_at_alpha"]) != RANK:
            raise ValueError(f"{hand} does not have the frozen rank-eight lower bound")
        if not rank7_p < alpha:
            raise ValueError(f"{hand} rank<=7 was not eliminated")
        if not rank8_p >= alpha:
            raise ValueError(f"{hand} rank<=8 was eliminated")
        summary[hand] = {
            "rank7_finite_batch_p": rank7_p,
            "rank8_finite_batch_p": rank8_p,
            "rank_lower_bound_at_alpha": RANK,
        }
    return summary


def kernel_basis(matrix: np.ndarray, rank: int = RANK) -> np.ndarray:
    _, _, right_adjoint = np.linalg.svd(matrix, full_matrices=False)
    return right_adjoint[rank:].conjugate().T


def kernel_projector(matrix: np.ndarray, rank: int = RANK) -> np.ndarray:
    kernel = kernel_basis(matrix, rank)
    return kernel @ kernel.conjugate().T


def bridge_matrix(plus: np.ndarray, minus: np.ndarray) -> np.ndarray:
    return kernel_projector(plus).conjugate() - kernel_projector(minus)


def bridge_residual(plus: np.ndarray, minus: np.ndarray) -> list[float]:
    return base.realify(bridge_matrix(plus, minus))


def descriptive_geometry(plus: np.ndarray, minus: np.ndarray) -> dict[str, object]:
    plus_values = np.linalg.svd(plus, compute_uv=False)
    minus_values = np.linalg.svd(minus, compute_uv=False)
    plus_kernel = kernel_basis(plus).conjugate()
    minus_kernel = kernel_basis(minus)
    cosines = np.linalg.svd(plus_kernel.conjugate().T @ minus_kernel, compute_uv=False)
    cosines = np.clip(cosines, 0.0, 1.0)
    angles = np.arccos(cosines)
    return {
        "plus_singular_values": plus_values.tolist(),
        "minus_R2_singular_values": minus_values.tolist(),
        "plus_rank8_to_rank9_singular_ratio": float(plus_values[7] / plus_values[8]),
        "minus_rank8_to_rank9_singular_ratio": float(minus_values[7] / minus_values[8]),
        "kernel_plane_principal_cosines": cosines.tolist(),
        "kernel_plane_principal_angles_radians": angles.tolist(),
        "kernel_projector_frobenius_distance": float(np.linalg.norm(bridge_matrix(plus, minus))),
    }


def matrices_from_means(
    old4: Mapping[str, float], old5: Mapping[str, float], fresh6: Mapping[str, float]
) -> dict[str, np.ndarray]:
    return {
        hand: base.hankel(old4, old5, fresh6, hand)
        for hand in ("plus", "minus")
    }


def score(manifest: Mapping[str, object]) -> dict[str, object]:
    paths = checked_inputs(manifest)
    alpha = float(manifest["statistic"]["alpha"])
    upstream = check_upstream_gate(paths["upstream_score"], alpha)

    old4_rows = base.read_old(paths["old4_batches"])
    old5_rows = base.read_radius5(paths["old5_batches"])
    fresh6_rows = base.read_new(paths["fresh6_batches"])
    if {len(old4_rows), len(old5_rows), len(fresh6_rows)} != {400}:
        raise ValueError("each dependency block must contain exactly 400 batches")

    old4 = base.old_means(old4_rows)
    old5 = base.radius5_means(old5_rows)
    fresh6 = base.new_means(fresh6_rows)
    full = matrices_from_means(old4, old5, fresh6)

    old4_deleted = [base.old_means(old4_rows, index) for index in range(len(old4_rows))]
    old5_deleted = [base.radius5_means(old5_rows, index) for index in range(len(old5_rows))]
    fresh6_deleted = [base.new_means(fresh6_rows, index) for index in range(len(fresh6_rows))]
    deleted_by_source = {
        "old4": [matrices_from_means(row, old5, fresh6) for row in old4_deleted],
        "old5": [matrices_from_means(old4, row, fresh6) for row in old5_deleted],
        "fresh6": [matrices_from_means(old4, old5, row) for row in fresh6_deleted],
    }

    point = bridge_residual(full["plus"], full["minus"])
    deleted = [
        [bridge_residual(row["plus"], row["minus"]) for row in deleted_by_source[source]]
        for source in ("old4", "old5", "fresh6")
    ]
    bridge_score = base.score_vector(point, deleted)
    finite_p = float(bridge_score["finite_batch_survival_p"])
    decision = (
        manifest["decision"]["p_below_alpha"]
        if finite_p < alpha
        else manifest["decision"]["p_at_or_above_alpha"]
    )
    return {
        "schema": SCHEMA,
        "status": "existing_data_rank8_bridge_reveal",
        "freeze_manifest_sha256": base.sha256(Path("analysis/p250_rank8_r2_kernel_bridge_freeze.json")),
        "dependency_group": manifest["dependency_group"],
        "upstream_gate": upstream,
        "matrix_shape_per_hand": list(full["plus"].shape),
        "conditional_rank": RANK,
        "kernel_dimension": full["plus"].shape[1] - RANK,
        "fixed_gauge": "plus identity; minus Alexander R2; compare coefficient spaces by conjugation",
        "primary_null": manifest["primary_null"],
        "decision_alpha": alpha,
        "bridge_score": bridge_score,
        "descriptive_geometry": descriptive_geometry(full["plus"], full["minus"]),
        "decision": decision,
        "claim_boundary": manifest["claim_boundary"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        default="analysis/p250_rank8_r2_kernel_bridge_freeze.json",
        type=Path,
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    result = score(manifest)
    payload = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
