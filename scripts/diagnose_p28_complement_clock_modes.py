#!/usr/bin/env python3
"""Project P28 birth-clock curvature into exact complement even/odd modes."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import mpmath as mp

from diagnose_p28_birth_clock_mixture import component_profiles
from diagnose_p28_threshold_profile_tail_rejection import canonical_modes
from score_p28_threshold_profile_tails import (
    DEFAULT_MANIFEST,
    ROOT,
    correlated_chi_square,
    jackknife_covariance,
    load_histogram,
    require,
    sha256,
)


EXACT_ORACLE = ROOT / "results/digital-alexander-filtration/latest.json"
REQUIRED_IDENTITIES = {
    "K_minus^G(pi)+K_plus^Ghat(reverse(pi))=N+1",
    "K_plus^G(pi)+K_minus^Ghat(reverse(pi))=N+1",
}


def residual(values: Sequence[float], modes: Sequence[Sequence[float]]) -> list[float]:
    return [math.fsum(weight * value for weight, value in zip(mode, values)) for mode in modes]


def complement_coordinates(
    profiles: Mapping[str, Any], z_grid: Sequence[float]
) -> tuple[dict[str, list[float]], list[dict[str, Any]], list[dict[str, Any]]]:
    modes = canonical_modes(z_grid, "stretched_4_over_3")
    output = {"even": [], "odd": []}
    labels = []
    pairs = []
    scale = math.sqrt(2.0)
    # Exact reflection turns K1^G-right into K2^Ghat-left and K1^G-left into K2^Ghat-right.
    pair_map = (
        ("second_birth_left", "right", "left"),
        ("second_birth_right", "left", "right"),
    )
    for orientation in ("first", "second"):
        for mapped_tail, k1_side, k2_side in pair_map:
            matching = residual(profiles["K1"][orientation]["log_density"][k1_side], modes)
            primal = residual(profiles["K2"][orientation]["log_density"][k2_side], modes)
            pair_even = [(left + right) / scale for left, right in zip(matching, primal)]
            pair_odd = [(left - right) / scale for left, right in zip(matching, primal)]
            output["even"].extend(pair_even)
            output["odd"].extend(pair_odd)
            for mode in range(2, 5):
                labels.append({"orientation": orientation, "mapped_tail": mapped_tail, "curvature_mode": mode})
            primal_norm = math.sqrt(math.fsum(value * value for value in primal))
            matching_norm = math.sqrt(math.fsum(value * value for value in matching))
            dot = math.fsum(left * right for left, right in zip(matching, primal))
            ratio = dot / math.fsum(value * value for value in primal)
            pairs.append({
                "orientation": orientation,
                "mapped_tail": mapped_tail,
                "matching_residual": matching,
                "primal_residual": primal,
                "euclidean_cosine": dot / (matching_norm * primal_norm),
                "matching_over_primal_amplitude": ratio,
                "relative_shape_mismatch_after_amplitude": math.sqrt(
                    math.fsum((left - ratio * right) ** 2 for left, right in zip(matching, primal))
                ) / matching_norm,
            })
    return output, labels, pairs


def closure_bases(pairs: Sequence[Mapping[str, Any]]) -> list[list[list[float]]]:
    bases = []
    for pair in pairs:
        primal = pair["primal_residual"]
        norm = math.sqrt(math.fsum(value * value for value in primal))
        direction = [value / norm for value in primal]
        rows: list[list[float]] = []
        for index in range(3):
            value = [1.0 if row == index else 0.0 for row in range(3)]
            for basis in [direction] + rows:
                projection = math.fsum(left * right for left, right in zip(value, basis))
                value = [left - projection * right for left, right in zip(value, basis)]
            length = math.sqrt(math.fsum(item * item for item in value))
            if length > 1e-10:
                rows.append([item / length for item in value])
            if len(rows) == 2:
                break
        require(len(rows) == 2, "failed to construct amplitude-free pair basis")
        bases.append(rows)
    return bases


def shape_closure_vector(
    pairs: Sequence[Mapping[str, Any]], bases: Sequence[Sequence[Sequence[float]]]
) -> list[float]:
    output = []
    for pair, rows in zip(pairs, bases):
        matching = pair["matching_residual"]
        primal = pair["primal_residual"]
        ratio = math.fsum(left * right for left, right in zip(matching, primal)) / math.fsum(
            value * value for value in primal
        )
        difference = [left - ratio * right for left, right in zip(matching, primal)]
        output.extend(
            math.fsum(weight * value for weight, value in zip(row, difference)) for row in rows
        )
    return output


def subset_score(
    vector: Sequence[float], covariance: Sequence[Sequence[float]], indices: Sequence[int], cutoff: float
) -> dict[str, Any]:
    return correlated_chi_square(
        [vector[index] for index in indices],
        [[covariance[left][right] for right in indices] for left in indices],
        cutoff,
    )


def diagnose_archive(archive: Mapping[str, Any], z_grid: Sequence[float], cutoff: float) -> dict[str, Any]:
    full_profiles = component_profiles(archive, z_grid)
    full, labels, pairs = complement_coordinates(full_profiles, z_grid)
    bases = closure_bases(pairs)
    full_joint = full["even"] + full["odd"]
    full_closure = shape_closure_vector(pairs, bases)
    deleted_joint = []
    deleted_closure = []
    for omitted in archive["batches"]:
        profiles = component_profiles(archive, z_grid, omitted)
        coordinates, _, deleted_pairs = complement_coordinates(profiles, z_grid)
        deleted_joint.append(coordinates["even"] + coordinates["odd"])
        deleted_closure.append(shape_closure_vector(deleted_pairs, bases))
    covariance = jackknife_covariance(deleted_joint)
    closure_covariance = jackknife_covariance(deleted_closure)
    width = len(full["even"])
    even_indices = list(range(width))
    odd_indices = list(range(width, 2 * width))
    scores = {
        "even": subset_score(full_joint, covariance, even_indices, cutoff),
        "odd": subset_score(full_joint, covariance, odd_indices, cutoff),
    }
    by_coordinate = {"even": {}, "odd": {}}
    for coordinate, offset in (("even", 0), ("odd", width)):
        for field in ("mapped_tail", "curvature_mode", "orientation"):
            by_coordinate[coordinate][field] = {}
            for value in sorted({str(label[field]) for label in labels}):
                indices = [
                    offset + index for index, label in enumerate(labels) if str(label[field]) == value
                ]
                by_coordinate[coordinate][field][value] = subset_score(
                    full_joint, covariance, indices, cutoff
                )
    closure_score = correlated_chi_square(full_closure, closure_covariance, cutoff)
    maximum_reconstruction_error = 0.0
    scale = math.sqrt(2.0)
    for pair_index, pair in enumerate(pairs):
        start = pair_index * 3
        for mode in range(3):
            even = full["even"][start + mode]
            odd = full["odd"][start + mode]
            maximum_reconstruction_error = max(
                maximum_reconstruction_error,
                abs((even + odd) / scale - pair["matching_residual"][mode]),
                abs((even - odd) / scale - pair["primal_residual"][mode]),
            )
    return {
        "N": archive["N"],
        "coordinate_scores": scores,
        "marginal_decomposition": by_coordinate,
        "pair_geometry": pairs,
        "amplitude_free_shape_closure": closure_score,
        "exact_coordinate_reconstruction_gate": {
            "maximum_absolute_error": maximum_reconstruction_error,
            "passed": maximum_reconstruction_error < 1e-12,
        },
    }


def combine(blocks: Sequence[Mapping[str, Any]], key: str) -> dict[str, Any]:
    chi_square = math.fsum(block["coordinate_scores"][key]["chi_square"] for block in blocks)
    df = sum(block["coordinate_scores"][key]["effective_df"] for block in blocks)
    probability = mp.gammainc(mp.mpf(df) / 2, chi_square / 2, mp.inf, regularized=True)
    return {"chi_square": chi_square, "effective_df": df, "survival_probability": float(probability)}


def combine_closure(blocks: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    chi_square = math.fsum(block["amplitude_free_shape_closure"]["chi_square"] for block in blocks)
    df = sum(block["amplitude_free_shape_closure"]["effective_df"] for block in blocks)
    probability = mp.gammainc(mp.mpf(df) / 2, chi_square / 2, mp.inf, regularized=True)
    return {"chi_square": chi_square, "effective_df": df, "survival_probability": float(probability)}


def run(manifest_path: Path) -> dict[str, Any]:
    exact = json.loads(EXACT_ORACLE.read_text(encoding="utf-8"))
    require(REQUIRED_IDENTITIES <= set(exact["exact_identities"]), "exact complement identities missing")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    z_grid = [float(value) for value in manifest["tail_design"]["absolute_z_grid"]]
    cutoff = float(manifest["tail_design"]["covariance_relative_eigenvalue_cutoff"])
    blocks = []
    for row in manifest["heldout_archives"]:
        path = ROOT / row["path"]
        require(sha256(path) == row["sha256"], f"checksum mismatch: {row['path']}")
        blocks.append(diagnose_archive(load_histogram(path), z_grid, cutoff))
    return {
        "schema": "matching-one/p28-complement-clock-curvature/v1",
        "issue": 28,
        "related_issue": 337,
        "status": "post_reveal_mechanism_diagnostic_no_tail_refit",
        "exact_oracle": {
            "path": str(EXACT_ORACLE.relative_to(ROOT)),
            "sha256": sha256(EXACT_ORACLE),
            "identities": sorted(REQUIRED_IDENTITIES),
            "derived_density_map": {
                "K1_G_right_z": "K2_Ghat_left_z",
                "K1_G_left_z": "K2_Ghat_right_z",
            },
            "boundary": (
                "The identity exchanges G with its matching graph Ghat. It does not assert that "
                "primal and matching curvature profiles are equal."
            ),
        },
        "coordinate_definition": {
            "even": "(mapped K1^G tail + K2^G tail)/sqrt(2)",
            "odd": "(mapped K1^G tail - K2^G tail)/sqrt(2)",
            "mapped_left_pair": "K1-right encodes K2^Ghat-left; compare with K2^G-left",
            "mapped_right_pair": "K1-left encodes K2^Ghat-right; compare with K2^G-right",
        },
        "global_coordinate_scores": {"even": combine(blocks, "even"), "odd": combine(blocks, "odd")},
        "global_amplitude_free_shape_closure": combine_closure(blocks),
        "blocks": blocks,
        "interpretation_boundary": (
            "Odd=0 would identify a common primal/matching curvature residual. Nonzero odd or failed "
            "amplitude-free collinearity shows distinct large-deviation shapes on the measured window; "
            "neither outcome identifies a continuum operator."
        ),
        "new_monte_carlo": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = run(args.manifest.resolve())
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
