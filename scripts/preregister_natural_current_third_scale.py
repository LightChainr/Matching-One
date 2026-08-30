#!/usr/bin/env python3
"""Freeze the N145 held-out natural-current transfer test."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping, Optional, Sequence


SOURCE_N = 85
TARGET_N = 145
H4_EXPONENT = -13.0 / 8.0


def preregister(n65_prereg_path: Path, n85_score_path: Path) -> dict[str, object]:
    n65 = json.loads(n65_prereg_path.read_text(encoding="utf-8"))
    n85 = json.loads(n85_score_path.read_text(encoding="utf-8"))
    if n65.get("schema") != "matching-one/P337-natural-current-scale-preregistration/v1":
        raise ValueError("unexpected N65 target source")
    if n85.get("schema") != "matching-one/P337-natural-current-scale-N85-score/v1":
        raise ValueError("unexpected N85 result source")

    delta65 = n65["natural_coordinate"]["value"][2]
    var65 = n65["natural_coordinate"]["covariance"][2][2]
    delta85 = n85["natural_coordinate"]["value"][2]
    var85 = n85["natural_coordinate"]["covariance"][2][2]
    h4_ratio = (TARGET_N / SOURCE_N) ** H4_EXPONENT
    log_scale_power = math.log(TARGET_N / SOURCE_N) / math.log(SOURCE_N / 65.0)
    descriptive_ratio = delta85 / delta65
    effective_target = delta85 * descriptive_ratio ** log_scale_power
    derivative_85 = (1.0 + log_scale_power) * effective_target / delta85
    derivative_65 = -log_scale_power * effective_target / delta65
    effective_variance = derivative_85**2 * var85 + derivative_65**2 * var65
    targets = {
        "zero": {"value": 0.0, "fit_variance": 0.0, "tier": "primary"},
        "source_fitted_scale_neutral": {
            "value": delta85, "fit_variance": var85, "tier": "primary",
        },
        "source_fitted_project_H4": {
            "value": h4_ratio * delta85,
            "fit_variance": h4_ratio**2 * var85,
            "tier": "primary", "exponent_in_area_N": H4_EXPONENT,
            "scale_ratio": h4_ratio,
        },
        "secondary_post_reveal_effective_transfer": {
            "value": effective_target,
            "fit_variance": effective_variance,
            "tier": "secondary",
            "construction": (
                "Delta85*(Delta85/Delta65)^[log(145/85)/log(85/65)]; "
                "descriptive ratio transport only, not a promoted exponent model"
            ),
            "descriptive_ratio_65_to_85": descriptive_ratio,
            "log_scale_power": log_scale_power,
        },
    }

    target_values = {name: row["value"] for name, row in targets.items()}
    pairs = [
        ("zero_vs_scale_neutral", "zero", "source_fitted_scale_neutral"),
        ("zero_vs_H4", "zero", "source_fitted_project_H4"),
        ("zero_vs_effective", "zero", "secondary_post_reveal_effective_transfer"),
        ("scale_neutral_vs_H4", "source_fitted_scale_neutral", "source_fitted_project_H4"),
        ("scale_neutral_vs_effective", "source_fitted_scale_neutral", "secondary_post_reveal_effective_transfer"),
        ("H4_vs_effective", "source_fitted_project_H4", "secondary_post_reveal_effective_transfer"),
    ]
    source_samples = n85["source"]["samples_per_shape"]
    source_se = n85["natural_coordinate"]["standard_error"][2]

    def required(name_a: str, name_b: str, z: float) -> float:
        gap = abs(target_values[name_a] - target_values[name_b])
        return source_samples * (z * source_se / gap) ** 2

    geometry = {
        "N": 145,
        "first": [12, 1], "second": [9, 8],
        "first_period_matrix": [[12, -1], [1, 12]],
        "second_period_matrix": [[9, -8], [8, 9]],
        "audit": {
            "equal_area": 12 * 12 + 1 == 9 * 9 + 8 * 8 == 145,
            "primitive_representations": math.gcd(12, 1) == math.gcd(9, 8) == 1,
            "smith_invariants": [[1, 145], [1, 145]],
            "nonassociate_under_Gaussian_units_and_conjugation": {12, 1} != {9, 8},
            "F3_reduction_nondegenerate": 145 % 3 != 0,
            "same_square_modulus": True,
            "legal": True,
        },
    }
    return {
        "schema": "matching-one/P337-natural-current-third-scale-preregistration/v1",
        "status": "frozen before any N145 projective-birth generation or inspection",
        "source": {
            "N65_archive_commit": "1714141",
            "N85_result_commit": "8783977",
            "N65_preregistration": str(n65_prereg_path),
            "N85_score": str(n85_score_path),
            "p_ref": n85["source"].get("p_ref", n65["source"]["p_ref"]),
        },
        "observable": {
            "definition": "K_A=p(1-p)Jminus_A/W_A=d_eta log W_A",
            "orientation_contrast": "Delta_K_A=K_A(second)-K_A(first)",
            "normalization_changed": False,
        },
        "geometry": geometry,
        "frozen_targets_at_N145": targets,
        "design": {
            "samples_per_shape": 2400000,
            "batches": 40,
            "seed": 202608337145,
            "replica_offset": 0,
            "engine": "1714141 projective-birth runner",
            "variance_source": "N85 independent 200k block only",
            "projected_standard_error": source_se * math.sqrt(source_samples / 2400000),
            "required_samples_3sigma": {
                label: required(a, b, 3.0) for label, a, b in pairs
            },
            "required_samples_5sigma": {
                label: required(a, b, 5.0) for label, a, b in pairs
            },
            "why_2_4m": (
                "exceeds the N85-variance 3-sigma H4-versus-secondary-effective "
                "requirement; no larger block is justified for this held-out discriminator"
            ),
        },
        "scoring_contract": {
            "primary_order": ["zero", "source_fitted_scale_neutral", "source_fitted_project_H4"],
            "secondary_order": ["secondary_post_reveal_effective_transfer"],
            "measurement_only_residual": "(observed-target)^2/Var_N145",
            "predictive_residual": "(observed-target)^2/(Var_N145+fit_variance)",
            "effective_transfer_boundary": (
                "secondary diagnostic derived only from the already revealed N65/N85 ratio; "
                "it cannot replace or refit the project H4 target"
            ),
            "scientific_question": (
                "does the faster N85 attenuation bend back toward project H4 curvature, "
                "continue along the descriptive effective transfer, or disappear"
            ),
        },
        "claim_boundary": "held-out N145 targets only; no N145 values and no fitted new exponent",
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    targets = payload["frozen_targets_at_N145"]
    design = payload["design"]
    lines = [
        "# Preregistration: third-scale natural charged current", "",
        "Frozen before generating or inspecting any N145 projective-birth data.", "",
        "The primitive equal-area pair `(12+i,9+8i)` passes the Gaussian legality audit: both norms are 145, both representations and Smith quotients are primitive, they are not associates/conjugates, and F3 reduction is nondegenerate.", "",
        "`K_A=p(1-p)Jminus_A/W_A` and the orientation order are unchanged.", "",
        "Frozen N145 targets:", "",
        f"- zero: `0`",
        f"- N85-fitted scale-neutral: `{targets['source_fitted_scale_neutral']['value']:.12g}`",
        f"- N85-fitted project H4: `{targets['source_fitted_project_H4']['value']:.12g}`",
        f"- secondary N65-to-N85 descriptive-ratio transfer: `{targets['secondary_post_reveal_effective_transfer']['value']:.12g}`", "",
        "The fourth target is explicitly secondary and is not a new exponent model.", "",
        f"Design: `{design['samples_per_shape']}` samples/shape, `{design['batches']}` aligned batches, seed `{design['seed']}`. N85 variance projects SE `{design['projected_standard_error']:.4g}`; this exceeds the 3-sigma H4-versus-effective gap requirement.", "",
        "Primary reporting keeps zero, scale-neutral and project H4 in their original order. The post-reveal effective transfer is reported separately to diagnose continuing fast attenuation versus curvature.", "",
    ]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n65-preregistration", type=Path, required=True)
    parser.add_argument("--n85-score", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    payload = preregister(args.n65_preregistration, args.n85_score)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
