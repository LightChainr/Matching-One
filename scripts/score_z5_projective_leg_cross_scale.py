#!/usr/bin/env python3
"""Frozen cross-scale sine-law versus exponential score for Issue 250."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from score_norm5_chiral_phase import matvec, quadratic
from score_z5_charged_threepoint import chi_square_survival, covariance_precision, zero_score
from score_z5_projective_leg_pair_transfer import CHANNELS, read_batches, means, transfer
from z5_projective_leg_cross_scale_mc import SCHEMA


SOURCE_L = math.sqrt(65.0)
TARGET_L = math.sqrt(101.0)
SOURCE_DISTANCES = (2, 3)
TARGET_DISTANCES = (2, 3, 4, 5)
FIXED_ALPHA = 1.25


def jackknife_covariance(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    count = len(rows)
    center = [sum(row[j] for row in rows) / count for j in range(len(rows[0]))]
    factor = (count - 1) / count
    return [[
        factor * sum((row[i] - center[i]) * (row[j] - center[j]) for row in rows)
        for j in range(len(center))
    ] for i in range(len(center))]


def add_covariance(first, second):
    return [[first[i][j] + second[i][j] for j in range(len(first))] for i in range(len(first))]


def decay_vector(values: Mapping[str, float], distances: Sequence[int]) -> list[float]:
    output = []
    for hand, charge in CHANNELS:
        origin = transfer(values, 1, hand, charge)
        if abs(origin) <= 0.0:
            raise ValueError("zero d1 transfer")
        for distance in distances:
            value = transfer(values, distance, hand, charge)
            if abs(value) <= 0.0:
                raise ValueError("zero target transfer")
            output.append(-math.log(abs(value / origin)))
    return output


def phase_vector(values: Mapping[str, float], distances: Sequence[int]) -> list[float]:
    output = []
    for hand, charge in CHANNELS:
        origin = transfer(values, 1, hand, charge)
        for distance in distances:
            ratio = transfer(values, distance, hand, charge) / origin
            output.append(math.atan2(ratio.imag, ratio.real))
    return output


def sine_design(length: float, distances: Sequence[int]) -> list[float]:
    base = math.sin(math.pi / length)
    row = [math.log(math.sin(math.pi * distance / length) / base) for distance in distances]
    return row * len(CHANNELS)


def exponential_design(distances: Sequence[int]) -> list[float]:
    return [float(distance - 1) for distance in distances] * len(CHANNELS)


def scalar_gls(point, covariance, design) -> dict:
    precision = covariance_precision(covariance)
    precision_point = matvec(precision, point)
    precision_design = matvec(precision, design)
    denominator = sum(design[i] * precision_design[i] for i in range(len(design)))
    beta = sum(design[i] * precision_point[i] for i in range(len(design))) / denominator
    residual = [point[i] - beta * design[i] for i in range(len(point))]
    chi_square = quadratic(residual, precision)
    return {
        "coefficient": beta,
        "standard_error": math.sqrt(1.0 / denominator),
        "residual": residual,
        "chi_square": chi_square,
        "degrees_of_freedom": len(point) - 1,
        "survival_p": chi_square_survival(chi_square, len(point) - 1),
    }


def source_freeze(source_batches: Sequence[dict]) -> dict:
    point = decay_vector(means(source_batches), SOURCE_DISTANCES)
    deleted = [
        decay_vector(means(source_batches, omitted), SOURCE_DISTANCES)
        for omitted in range(len(source_batches))
    ]
    covariance = jackknife_covariance(deleted)
    cylinder = scalar_gls(point, covariance, sine_design(SOURCE_L, SOURCE_DISTANCES))
    exponential = scalar_gls(point, covariance, exponential_design(SOURCE_DISTANCES))
    return {
        "source_parent_order": 65,
        "source_length": SOURCE_L,
        "fit_distances": list(SOURCE_DISTANCES),
        "common_cylinder_alpha": cylinder,
        "common_exponential_mass": exponential,
    }


def batch_resolution(batches: Sequence[dict]) -> dict:
    output = {}
    for distance in (1, *TARGET_DISTANCES):
        for hand, charge in CHANNELS:
            key = f"d{distance}_T{charge}_{hand}"
            rows = [[row[key + "_re"] / row["samples"], row[key + "_im"] / row["samples"]] for row in batches]
            point = [sum(row[j] for row in rows) / len(rows) for j in range(2)]
            covariance = [[
                sum((row[i] - point[i]) * (row[j] - point[j]) for row in rows)
                / (len(rows) * (len(rows) - 1))
                for j in range(2)
            ] for i in range(2)]
            output[key] = {
                "point_re_im": point,
                "covariance": covariance,
                "real_abs_z": abs(point[0]) / math.sqrt(covariance[0][0]),
            }
    return output


def propagated_model_score(
    target_point,
    target_deleted,
    target_design,
    coefficient,
    source_deleted_coefficients=None,
):
    residual = [target_point[i] - coefficient * target_design[i] for i in range(len(target_point))]
    target_residual_deleted = [
        [row[i] - coefficient * target_design[i] for i in range(len(row))]
        for row in target_deleted
    ]
    covariance = jackknife_covariance(target_residual_deleted)
    if source_deleted_coefficients is not None:
        source_residual_deleted = [
            [target_point[i] - value * target_design[i] for i in range(len(target_point))]
            for value in source_deleted_coefficients
        ]
        covariance = add_covariance(covariance, jackknife_covariance(source_residual_deleted))
    full = zero_score(residual, covariance)
    heldout_indices = [len(TARGET_DISTANCES) * channel + 3 for channel in range(len(CHANNELS))]
    heldout = [residual[index] for index in heldout_indices]
    heldout_covariance = [[covariance[i][j] for j in heldout_indices] for i in heldout_indices]
    return {
        "residual_order": [
            f"{hand}_r{charge}_d{distance}"
            for hand, charge in CHANNELS for distance in TARGET_DISTANCES
        ],
        "residual": residual,
        "covariance": covariance,
        "all_target_distances_score": full,
        "heldout_d5": {
            "order": [f"{hand}_r{charge}_d5" for hand, charge in CHANNELS],
            "residual": heldout,
            "covariance": heldout_covariance,
            "zero_score": zero_score(heldout, heldout_covariance),
        },
    }


def score(payload: dict, target_batches: Sequence[dict], source_batches: Sequence[dict], manifest: Mapping[str, object]) -> dict:
    if payload.get("schema") != SCHEMA or not payload["exact_gate"]["passed"]:
        raise ValueError("wrong or failed cross-scale response")
    for key, expected in manifest["run"].items():
        if payload["run"].get(key) != expected:
            raise ValueError(f"run differs from manifest for {key}")
    frozen = source_freeze(source_batches)
    expected_freeze = manifest["source_fit"]
    if abs(frozen["common_cylinder_alpha"]["coefficient"] - expected_freeze["common_cylinder_alpha"]) > 1e-12:
        raise ValueError("source-fitted common alpha changed")
    if abs(frozen["common_exponential_mass"]["coefficient"] - expected_freeze["common_exponential_mass"]) > 1e-12:
        raise ValueError("source-fitted common mass changed")
    raw = batch_resolution(target_batches)
    minimum_z = min(row["real_abs_z"] for row in raw.values())
    gate = minimum_z >= manifest["transfer_gate"]["minimum_real_abs_z"]
    if not gate:
        models = {"status": "locked_target_transfer_gate_failed", "computed": False}
    else:
        target_point = decay_vector(means(target_batches), TARGET_DISTANCES)
        target_deleted = [
            decay_vector(means(target_batches, omitted), TARGET_DISTANCES)
            for omitted in range(len(target_batches))
        ]
        source_sine = sine_design(SOURCE_L, SOURCE_DISTANCES)
        source_exp = exponential_design(SOURCE_DISTANCES)
        source_deleted_alpha = []
        source_deleted_mass = []
        for omitted in range(len(source_batches)):
            kept = [row for index, row in enumerate(source_batches) if index != omitted]
            point = decay_vector(means(kept), SOURCE_DISTANCES)
            nested = [decay_vector(means(kept, inner), SOURCE_DISTANCES) for inner in range(len(kept))]
            covariance = jackknife_covariance(nested)
            source_deleted_alpha.append(scalar_gls(point, covariance, source_sine)["coefficient"])
            source_deleted_mass.append(scalar_gls(point, covariance, source_exp)["coefficient"])
        target_sine = sine_design(TARGET_L, TARGET_DISTANCES)
        target_exp = exponential_design(TARGET_DISTANCES)
        fixed = propagated_model_score(target_point, target_deleted, target_sine, FIXED_ALPHA)
        fixed.update({"alpha": FIXED_ALPHA, "meaning": "fixed thermal 2x=5/4 cylinder law"})
        free = propagated_model_score(
            target_point, target_deleted, target_sine,
            frozen["common_cylinder_alpha"]["coefficient"], source_deleted_alpha,
        )
        free.update({
            "alpha": frozen["common_cylinder_alpha"]["coefficient"],
            "alpha_source_standard_error": frozen["common_cylinder_alpha"]["standard_error"],
            "meaning": "common alpha fitted only on source N65 d1-d3",
        })
        exponential = propagated_model_score(
            target_point, target_deleted, target_exp,
            frozen["common_exponential_mass"]["coefficient"], source_deleted_mass,
        )
        exponential.update({
            "mass": frozen["common_exponential_mass"]["coefficient"],
            "mass_source_standard_error": frozen["common_exponential_mass"]["standard_error"],
            "meaning": "lattice-unit finite correlation mass fitted only on source N65 d1-d3",
        })
        phase_point = phase_vector(means(target_batches), TARGET_DISTANCES)
        phase_deleted = [phase_vector(means(target_batches, omitted), TARGET_DISTANCES) for omitted in range(len(target_batches))]
        phase_covariance = jackknife_covariance(phase_deleted)
        models = {
            "status": "cross_scale_models_revealed", "computed": True,
            "fixed_thermal_sine": fixed,
            "source_fitted_sine": free,
            "single_exponential": exponential,
            "target_deck_phase_zero": {
                "point": phase_point, "covariance": phase_covariance,
                "zero_score": zero_score(phase_point, phase_covariance),
            },
        }
    return {
        "schema": "matching-one/z5-projective-leg-cross-scale-score/v1",
        "status": "fresh_cross_scale_reveal",
        "source_fit": frozen,
        "target_parent_order": 101,
        "target_length": TARGET_L,
        "minimum_d1_d5_real_abs_z": minimum_z,
        "target_transfer_gate_passed": gate,
        "raw_target_transfer": raw,
        "models": models,
        "heldout_contract": "all N101 distances are geometry-held-out; d5 is additionally reported as the decisive antipodal separation",
        "cubic_fields_used": [],
        "claim_boundary": [
            "The fixed 5/4 target was frozen before N101 collection.",
            "The free alpha and exponential mass use only N65 d1-d3 and propagate source uncertainty.",
            "Cylinder compatibility is a propagator-shape result, not a local-primary or OPE identification.",
        ],
    }


def render(result: Mapping[str, object]) -> str:
    lines = [
        "# P250 cross-scale projective-leg propagator", "",
        f"Target d1-d5 minimum real z: `{result['minimum_d1_d5_real_abs_z']}`; gate `{result['target_transfer_gate_passed']}`.", "",
    ]
    if result["models"]["computed"]:
        lines += ["| model | all target chi2/df | p | held-out d5 chi2/df | p |", "|---|---:|---:|---:|---:|"]
        for key in ("fixed_thermal_sine", "source_fitted_sine", "single_exponential"):
            row = result["models"][key]
            full = row["all_target_distances_score"]
            held = row["heldout_d5"]["zero_score"]
            lines.append(f"| {key} | {full['chi_square']:.6g}/{full['degrees_of_freedom']} | {full['survival_p']:.6g} | {held['chi_square']:.6g}/{held['degrees_of_freedom']} | {held['survival_p']:.6g} |")
        phase = result["models"]["target_deck_phase_zero"]["zero_score"]
        lines += ["", f"Deck-phase zero score: `{phase['chi_square']}/{phase['degrees_of_freedom']}`, p `{phase['survival_p']}`."]
    else:
        lines += [f"Models locked: `{result['models']['status']}`."]
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("response", type=Path)
    parser.add_argument("target_batches", type=Path)
    parser.add_argument("source_batches", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = score(
        json.loads(args.response.read_text()), read_batches(args.target_batches),
        read_batches(args.source_batches), json.loads(args.manifest.read_text()),
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
