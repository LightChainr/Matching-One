#!/usr/bin/env python3
"""Separate single-radius cone failure from cross-radius spectrum sharing for P418."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Callable, Optional, Sequence

import numpy as np

import score_p406_spatial_fourier_cone as p406
import score_p418_crt_degauging as p418


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/huawei-20260830/P418-radius-resolved-elimination/score.json"
P418_SCORE = ROOT / "results/huawei-20260830/P418-crt-degauging/score.json"
SCHEMA = "matching-one/p418-radius-resolved-elimination/v1"
RADIUS_NAMES = ("radius4", "radius5", "radius6")


def quantiles(values: Sequence[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {
        "q01": float(np.quantile(array, 0.01)),
        "q50": float(np.quantile(array, 0.50)),
        "q99": float(np.quantile(array, 0.99)),
    }


def pvalue(observed: float, reference: Sequence[float]) -> float:
    return (1 + sum(value >= observed for value in reference)) / (len(reference) + 1)


def fit_family(
    blocks: list[dict[str, object]],
    design_builder: Callable[[int, list[tuple[int, int]]], np.ndarray],
    bootstrap: int,
    seed: int,
) -> dict[str, Any]:
    whitened = [
        p406.whiten(block["values"], design_builder(index, block["coordinates"]))
        for index, block in enumerate(blocks)
    ]
    matrices = [row["X"] for row in whitened]
    vectors = [row["y"] for row in whitened]
    lengths = [len(vector) for vector in vectors]
    starts = np.cumsum([0] + lengths)
    shared_matrix = np.vstack(matrices)
    shared_vector = np.concatenate(vectors)

    shared_weights, shared_distance = p406.fit_nonnegative(shared_matrix, shared_vector)
    shared_fitted = shared_matrix @ shared_weights
    separate = []
    separate_fitted = []
    separate_distance = 0.0
    for name, matrix, vector, white in zip(RADIUS_NAMES, matrices, vectors, whitened):
        weights, distance = p406.fit_nonnegative(matrix, vector)
        separate_distance += distance
        separate_fitted.append(matrix @ weights)
        separate.append({
            "radius": name,
            "distance_squared": distance,
            "resolved_modes": len(vector),
            "design_rank": int(np.linalg.matrix_rank(matrix, tol=1e-10)),
            "positive_weights": int(np.count_nonzero(weights > max(weights.max() * 1e-10, 1e-14))),
        })

    observed_penalty = shared_distance - separate_distance
    rng = np.random.default_rng(seed)
    noises = rng.standard_normal((bootstrap, len(shared_vector)))
    shared_reference = []
    penalty_reference = []
    radius_reference = [[] for _ in blocks]
    for noise in noises:
        _, common_value = p406.fit_nonnegative(shared_matrix, shared_fitted + noise)
        shared_reference.append(common_value)
        split_value = 0.0
        for index, matrix in enumerate(matrices):
            piece = noise[starts[index] : starts[index + 1]]
            _, value = p406.fit_nonnegative(matrix, separate_fitted[index] + piece)
            radius_reference[index].append(value)
            # Under the common-spectrum null, use the common fitted center for
            # both nested fits so the sharing penalty has the correct null.
            common_piece = shared_fitted[starts[index] : starts[index + 1]]
            _, nested_value = p406.fit_nonnegative(matrix, common_piece + piece)
            split_value += nested_value
        penalty_reference.append(common_value - split_value)

    for row, reference in zip(separate, radius_reference):
        row["bootstrap_p"] = pvalue(row["distance_squared"], reference)
        row["bootstrap_quantiles"] = quantiles(reference)
        row["decision"] = "cone_not_rejected" if row["bootstrap_p"] >= 0.01 else "cone_rejected"

    penalty_p = pvalue(observed_penalty, penalty_reference)
    return {
        "per_radius_separate_spectra": separate,
        "sum_separate_distance_squared": separate_distance,
        "shared_spectrum": {
            "distance_squared": shared_distance,
            "bootstrap_p": pvalue(shared_distance, shared_reference),
            "bootstrap_quantiles": quantiles(shared_reference),
            "decision": "cone_not_rejected" if pvalue(shared_distance, shared_reference) >= 0.01 else "cone_rejected",
        },
        "cross_radius_sharing_penalty": {
            "distance_squared": observed_penalty,
            "bootstrap_p_under_shared_spectrum": penalty_p,
            "bootstrap_quantiles_under_shared_spectrum": quantiles(penalty_reference),
            "decision": "shared_spectrum_rejected" if penalty_p < 0.01 else "shared_spectrum_not_rejected",
        },
    }


def build_result(paths: list[tuple[Path, str]], bootstrap: int, seed: int) -> dict[str, Any]:
    archives = [p406.read_block(path, digest) for path, digest in paths]
    exact = p418.exact_section_and_masks()
    published = json.loads(P418_SCORE.read_text(encoding="utf-8"))
    channels = {}
    maximum_shared_replay_error = 0.0
    maximum_shared_p_error = 0.0

    for hand_index, hand in enumerate(p406.HANDS):
        for charge in p406.CHARGES:
            key = f"{hand}_r{charge}"
            blocks = [archive[(hand, charge)] for archive in archives]
            channel_seed = seed + 1000 * hand_index + charge
            raw = fit_family(
                blocks,
                lambda _index, coordinates: p406.design(coordinates),
                bootstrap,
                channel_seed,
            )
            masked = fit_family(
                blocks,
                lambda _index, coordinates: p418.masked_design(exact, hand, charge, coordinates),
                bootstrap,
                channel_seed,
            )
            expected = published["channels"][key]
            maximum_shared_replay_error = max(
                maximum_shared_replay_error,
                abs(raw["shared_spectrum"]["distance_squared"] - expected["raw_cone"]["distance_squared"]),
                abs(masked["shared_spectrum"]["distance_squared"] - expected["masked_cone"]["distance_squared"]),
            )
            maximum_shared_p_error = max(
                maximum_shared_p_error,
                abs(raw["shared_spectrum"]["bootstrap_p"] - expected["raw_cone"]["bootstrap_p"]),
                abs(masked["shared_spectrum"]["bootstrap_p"] - expected["masked_cone"]["bootstrap_p"]),
            )
            per_radius = []
            for raw_row, masked_row in zip(
                raw["per_radius_separate_spectra"], masked["per_radius_separate_spectra"]
            ):
                per_radius.append({
                    "radius": raw_row["radius"],
                    "raw_cone": raw_row,
                    "masked_cone": masked_row,
                    "mask_penalty_over_raw": masked_row["distance_squared"] - raw_row["distance_squared"],
                })
            channels[key] = {
                "per_radius": per_radius,
                "raw_family": raw,
                "masked_family": masked,
            }

    if maximum_shared_replay_error > 1e-9 or maximum_shared_p_error > 1e-15:
        raise AssertionError("published shared-spectrum P418 score did not replay")

    single_radius_mask_rejections = [
        f"{key}:{row['radius']}"
        for key, channel in channels.items()
        for row in channel["per_radius"]
        if row["masked_cone"]["decision"] == "cone_rejected"
    ]
    shared_mask_rejections = [
        key for key, channel in channels.items()
        if channel["masked_family"]["shared_spectrum"]["decision"] == "cone_rejected"
    ]
    sharing_penalty_rejections = [
        key for key, channel in channels.items()
        if channel["masked_family"]["cross_radius_sharing_penalty"]["decision"] == "shared_spectrum_rejected"
    ]
    if single_radius_mask_rejections:
        mechanism = "masked_positive_stationarity_or_transport_already_fails_within_at_least_one_radius"
    elif shared_mask_rejections and sharing_penalty_rejections:
        mechanism = "each_radius_survives_but_one_shared_cross_radius_masked_spectrum_fails"
    else:
        mechanism = "radius_decomposition_is_mixed_or_inconclusive"
    return {
        "schema": SCHEMA,
        "status": "archived_production_radius_resolved_once",
        "issues": [418, 406, 250],
        "new_monte_carlo": False,
        "inputs": [{"path": str(path.resolve().relative_to(ROOT)), "sha256": digest} for path, digest in paths],
        "bootstrap": {"replicates": bootstrap, "seed": seed, "decision_alpha": 0.01},
        "shared_score_replay": {
            "source": str(P418_SCORE.relative_to(ROOT)),
            "maximum_distance_error": maximum_shared_replay_error,
            "maximum_p_error": maximum_shared_p_error,
            "passed": True,
        },
        "channels": channels,
        "decision": {
            "mechanism": mechanism,
            "single_radius_mask_rejections": single_radius_mask_rejections,
            "shared_mask_rejections": shared_mask_rejections,
            "masked_sharing_penalty_rejections": sharing_penalty_rejections,
        },
        "claim_boundary": (
            "This separates per-radius cone compatibility from one shared spectrum using the pinned archive "
            "means and complete within-radius covariance. It adds no samples and identifies no state count, "
            "field, Jordan block or ordered-memory mechanism."
        ),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius4", type=Path, required=True)
    parser.add_argument("--radius4-sha256", required=True)
    parser.add_argument("--radius5", type=Path, required=True)
    parser.add_argument("--radius5-sha256", required=True)
    parser.add_argument("--radius6", type=Path, required=True)
    parser.add_argument("--radius6-sha256", required=True)
    parser.add_argument("--bootstrap", type=int, default=250)
    parser.add_argument("--seed", type=int, default=40610120260830)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    result = build_result([
        (args.radius4, args.radius4_sha256),
        (args.radius5, args.radius5_sha256),
        (args.radius6, args.radius6_sha256),
    ], args.bootstrap, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["decision"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
