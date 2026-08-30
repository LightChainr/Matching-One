#!/usr/bin/env python3
"""Freeze the N680 heldout discriminator from the three-generation recurrence."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Mapping, Optional, Sequence


SAMPLES = 120_000_000
BATCHES = 80
SEED = 202608337680


def preregister(recurrence_path: Path, n340_path: Path) -> dict[str, object]:
    recurrence = json.loads(recurrence_path.read_text(encoding="utf-8"))
    n340 = json.loads(n340_path.read_text(encoding="utf-8"))
    if recurrence.get("schema") != "matching-one/P337-three-generation-H4-recurrence/v1":
        raise ValueError("unexpected recurrence schema")
    if n340.get("schema") != "matching-one/P337-N340-second-child-score/v1":
        raise ValueError("unexpected N340 schema")

    covectors = [Fraction(-4633, 7225), Fraction(6887, 7225)]
    dc = float(covectors[1] - covectors[0])
    transform = [
        [-1.0 / dc, 1.0 / dc],
        [2.0 * float(covectors[1]) / dc, -2.0 * float(covectors[0]) / dc],
    ]
    rmodel = recurrence["two_mode_recurrence"]
    comparators = recurrence["comparators"]
    frozen = [
        ("two_mode_recurrence", rmodel["N680_H4_amplitude_prediction"], rmodel["N680_prediction_standard_error_delta"] ** 2, "primary recurrence target"),
        ("single_frozen_lambda0", comparators["single_frozen_lambda0"]["N680_H4_amplitude_prediction"], comparators["single_frozen_lambda0"]["N680_prediction_standard_error"] ** 2, "frozen nominal H4 single mode"),
        ("single_free_lambda", comparators["single_free_lambda"]["N680_H4_amplitude_prediction"], comparators["single_free_lambda"]["N680_prediction_standard_error"] ** 2, "same-lineage descriptive free transfer"),
        ("scale_neutral", comparators["scale_neutral"]["N680_H4_amplitude_prediction"], comparators["scale_neutral"]["N680_prediction_standard_error"] ** 2, "no-decay control"),
    ]
    models = []
    for name, amplitude, variance, role in frozen:
        models.append({
            "name": name, "role": role,
            "H4_amplitude": amplitude,
            "H4_amplitude_target_variance": variance,
            "absolute_K_A": [amplitude * float(c) for c in covectors],
            "absolute_target_covariance": [
                [variance * float(left) * float(right) for right in covectors]
                for left in covectors
            ],
            "pair_second_minus_first": amplitude * dc,
            "pair_sign": "negative",
        })

    n340_variance = n340["decomposition"]["measurement_covariance"][0][0]
    projected_variance = n340_variance * n340["source"]["samples_per_shape"] / SAMPLES
    projected_se = math.sqrt(projected_variance)
    values = {model["name"]: model["H4_amplitude"] for model in models}
    gaps = {
        "two_mode_vs_free_single": abs(values["two_mode_recurrence"] - values["single_free_lambda"]),
        "two_mode_vs_fixed_lambda0": abs(values["two_mode_recurrence"] - values["single_frozen_lambda0"]),
        "two_mode_vs_scale_neutral": abs(values["two_mode_recurrence"] - values["scale_neutral"]),
    }
    required_3sigma = n340["source"]["samples_per_shape"] * (
        3.0 * math.sqrt(n340_variance) / gaps["two_mode_vs_free_single"]
    ) ** 2
    return {
        "schema": "matching-one/P337-N680-heldout-preregistration/v1",
        "status": "frozen before any N680 projective-birth generation or inspection",
        "source": {
            "recurrence_commit": "4024a7c",
            "recurrence": str(recurrence_path),
            "N340_variance_source": str(n340_path),
            "p_ref": n340["source"]["p_ref"],
        },
        "geometry": {
            "N": 680,
            "first": [22, 14], "second": [26, 2],
            "first_period_matrix": [[22, -14], [14, 22]],
            "second_period_matrix": [[26, -2], [2, 26]],
            "parent_N340": [[18, 4], [14, 12]],
            "raw_multiplier_images": [[14, 22], [2, 26]],
            "canonicalization": "D4-equivalent canonical representatives only",
            "common_multiplier": [1, 1],
            "H4_covectors_exact": [str(value) for value in covectors],
            "parent_H4_covectors_exact": ["4633/7225", "-6887/7225"],
            "exact_angle_flip": True,
            "Smith_classes": [[2, 340], [2, 340]],
            "F3_reduction_nondegenerate": True,
            "charged_projective_scalar": {"value": 0.5, "frozen_amplitude_target": 0.0},
        },
        "observable": {
            "definition": "K_A=p(1-p)Jminus_A/W_A=d_eta log W_A",
            "vector_order": ["K_A_first_22_plus_14i", "K_A_second_26_plus_2i"],
            "decomposition": {
                "H4_amplitude": "(K_second-K_first)/(c_second-c_first)",
                "A_projective_scalar": "2(c_second*K_first-c_first*K_second)/(c_second-c_first)",
                "transform": transform,
            },
        },
        "frozen_models_in_scoring_order": models,
        "production": {
            "samples_per_shape": SAMPLES, "batches": BATCHES, "seed": SEED,
            "replica_offset": 0, "threads": 16,
            "engine": "1714141 projective-birth integer-period runner",
            "machine": {
                "name": "DevEnvC_HZsCM6", "id": "033945d8bf8b47a7acf475c595169e07",
                "class": "Kunpeng AArch64, 16 vCPU, 32 GiB",
                "selection": "idle at freeze audit",
            },
            "full_batch_covariance_required": True,
            "archive_policy": "retain exact raw remotely; commit lossless compressed raw only if below GitHub hard limit",
        },
        "power": {
            "variance_source": "N340 12M/shape observed H4-amplitude variance",
            "projected_measurement_standard_error": projected_se,
            "fixed_amplitude_gaps": gaps,
            "projected_measurement_z_for_gaps": {key: value / projected_se for key, value in gaps.items()},
            "required_samples_for_3sigma_two_mode_vs_free_single": required_3sigma,
            "why_120m": "smallest round design above the N340-variance 3-sigma requirement; 80M would project only 2.47 sigma",
            "target_uncertainty_note": "power is for fixed point forecasts; source-fit uncertainty is retained separately in predictive scoring",
        },
        "scoring_contract": {
            "primary_coordinate": "H4_amplitude",
            "fixed_model_order": [model["name"] for model in models],
            "report_measurement_only_and_source_uncertainty_aware_residuals": True,
            "projective_scalar_zero_control": True,
            "pair_expected_sign": "negative",
            "no_model_refit": True,
            "no_exponent_fit": True,
            "no_post_reveal_basis_change": True,
        },
        "claim_boundary": "one preregistered N680 same-lineage heldout discriminator; model point forecasts and projective zero control only",
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    power = payload["power"]
    lines = [
        "# Preregistration: N680 same-lineage heldout", "",
        "Frozen before generating or inspecting N680 data. N680 `(22+14i,26+2i)` is the exact next `1+i` child of N340; H4 flips negative again and the projective scalar remains a zero-amplitude control.", "",
        "Fixed H4-amplitude forecasts:", "",
    ]
    for model in payload["frozen_models_in_scoring_order"]:
        lines.append(f"- `{model['name']}`: `{model['H4_amplitude']:+.12g}` (target SE `{math.sqrt(model['H4_amplitude_target_variance']):.3g}`)")
    lines.extend([
        "",
        f"N340 variance implies 3-sigma two-mode/free-single separation at {power['required_samples_for_3sigma_two_mode_vs_free_single']:.0f} samples/shape. The frozen 120M design projects SE `{power['projected_measurement_standard_error']:.4g}` and separations `{power['projected_measurement_z_for_gaps']['two_mode_vs_free_single']:.3f}` (two/free), `{power['projected_measurement_z_for_gaps']['two_mode_vs_fixed_lambda0']:.3f}` (two/fixed), and `{power['projected_measurement_z_for_gaps']['two_mode_vs_scale_neutral']:.3f}` (two/neutral).", "",
        "Power refers to fixed point forecasts. Source-fit uncertainty remains in a separate predictive score; no model, exponent, or basis may be refit after reveal.", "",
    ])
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recurrence", type=Path, required=True)
    parser.add_argument("--n340", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    payload = preregister(args.recurrence, args.n340)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
