#!/usr/bin/env python3
"""Freeze the N170 exact angle-flip charged-current production."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Mapping, Optional, Sequence


AREA_POWER = -13.0 / 8.0
TARGET_N = 170
SAMPLES = 8_000_000
BATCHES = 80
SEED = 202608337170


def chi4_real(a: int, b: int) -> Fraction:
    n = a * a + b * b
    return Fraction(a**4 - 6 * a * a * b * b + b**4, n * n)


def preregister(crosswalk_path: Path, n145_score_path: Path) -> dict[str, object]:
    crosswalk = json.loads(crosswalk_path.read_text(encoding="utf-8"))
    n145 = json.loads(n145_score_path.read_text(encoding="utf-8"))
    if crosswalk.get("schema") != "matching-one/P337-natural-current-geometry-crosswalk/v1":
        raise ValueError("unexpected geometry crosswalk schema")
    if n145.get("schema") != "matching-one/P337-natural-current-third-scale-N145-score/v1":
        raise ValueError("unexpected N145 variance source")

    h4 = next(model for model in crosswalk["models"]
              if model["name"] == "one_H4_geometry_covector")
    beta = h4["parameters"][0]
    beta_variance = h4["parameter_covariance"][0][0]
    child = ((11, 7), (13, 1))
    parent = ((9, 2), (7, 6))
    covectors = [chi4_real(*row) for row in child]
    parent_covectors = [chi4_real(*row) for row in parent]
    if any(child_value != -parent_value
           for child_value, parent_value in zip(covectors, parent_covectors)):
        raise AssertionError("N170 is not the exact H4 angle-flip child of N85")
    radial = TARGET_N ** AREA_POWER
    absolute_targets = [beta * radial * float(value) for value in covectors]
    gradients = [radial * float(value) for value in covectors]
    fit_covariance = [[beta_variance * left * right for right in gradients]
                      for left in gradients]
    delta_covector = float(covectors[1] - covectors[0])
    pair_target = absolute_targets[1] - absolute_targets[0]
    h4_amplitude_target = beta * radial
    h4_amplitude_fit_variance = beta_variance * radial * radial
    transform = [
        [-1.0 / delta_covector, 1.0 / delta_covector],
        [2.0 * float(covectors[1]) / delta_covector,
         -2.0 * float(covectors[0]) / delta_covector],
    ]

    source_covariance = [row[:2] for row in n145["natural_coordinate"]["covariance"][:2]]
    scale = n145["source"]["samples_per_shape"] / SAMPLES

    def projected_variance(weights: Sequence[float]) -> float:
        return scale * math.fsum(
            weights[i] * source_covariance[i][j] * weights[j]
            for i in range(2) for j in range(2)
        )

    pair_weights = [-1.0, 1.0]
    pair_se = math.sqrt(projected_variance(pair_weights))
    h4_se = math.sqrt(projected_variance(transform[0]))
    scalar_se = math.sqrt(projected_variance(transform[1]))
    required_5sigma = n145["source"]["samples_per_shape"] * (
        5.0 * n145["natural_coordinate"]["standard_error"][2] / abs(pair_target)
    ) ** 2

    return {
        "schema": "matching-one/P337-N170-angle-flip-preregistration/v1",
        "status": "frozen before any N170 projective-birth generation or inspection",
        "source": {
            "geometry_crosswalk_commit": "186d72a",
            "geometry_crosswalk": str(crosswalk_path),
            "H4_model_fit_scales": [65, 85],
            "N145_role": "variance source only; not used to fit the frozen H4 direction",
            "p_ref": n145["source"]["p_ref"],
        },
        "geometry": {
            "N": TARGET_N,
            "first": [11, 7], "second": [13, 1],
            "first_period_matrix": [[11, -7], [7, 11]],
            "second_period_matrix": [[13, -1], [1, 13]],
            "parent_N85": [[9, 2], [7, 6]],
            "common_multiplier": [1, 1],
            "H4_covectors_exact": [str(value) for value in covectors],
            "parent_H4_covectors_exact": [str(value) for value in parent_covectors],
            "exact_angle_flip": True,
            "charged_projective_scalar": {
                "value": 0.5,
                "transport": "unchanged from parent because q_A^2=(u+H_F3)/2",
            },
            "tau": "i",
            "F3_reduction_nondegenerate": TARGET_N % 3 != 0,
        },
        "observable": {
            "definition": "K_A=p(1-p)Jminus_A/W_A=d_eta log W_A",
            "vector_order": ["K_A_first_11_plus_7i", "K_A_second_13_plus_i"],
            "decomposition": {
                "H4_amplitude": "(K_second-K_first)/(c_second-c_first)",
                "A_projective_scalar": "2(c_second*K_first-c_first*K_second)/(c_second-c_first)",
                "transform": transform,
            },
        },
        "frozen_H4_only_prediction": {
            "model_parameter_beta": beta,
            "model_parameter_variance": beta_variance,
            "area_power": AREA_POWER,
            "absolute_K_A": absolute_targets,
            "absolute_fit_covariance": fit_covariance,
            "pair_second_minus_first": pair_target,
            "pair_sign": "negative",
            "H4_amplitude": h4_amplitude_target,
            "H4_amplitude_fit_variance": h4_amplitude_fit_variance,
            "A_projective_scalar": 0.0,
        },
        "primary_contrasts": {
            "geometry_sign_flip_vs_scalar": {
                "statistic": "Delta_K_A=K_second-K_first",
                "H4_target": pair_target,
                "scalar_or_geometry_blind_target": 0.0,
                "decision_reading": "negative resolves the exact angle flip; zero is scalar/geometry-blind",
            },
            "curvature_vs_projective_split": {
                "H4_amplitude_target": h4_amplitude_target,
                "A_projective_scalar_target": 0.0,
                "decision_reading": (
                    "departure along H4 amplitude is scale curvature; departure in the "
                    "orthogonal A scalar is a charged/projective common mode"
                ),
            },
        },
        "production": {
            "samples_per_shape": SAMPLES,
            "batches": BATCHES,
            "seed": SEED,
            "replica_offset": 0,
            "threads": 16,
            "engine": "1714141 projective-birth integer-period runner",
            "machine": {
                "name": "DevEnvC_HZsCM6",
                "id": "033945d8bf8b47a7acf475c595169e07",
                "class": "Kunpeng AArch64, 16 vCPU, 32 GiB",
                "selection": "idle at preregistration audit; no compute process",
            },
            "full_batch_covariance_required": True,
        },
        "power": {
            "variance_source": "N145 2.4M/shape full paired covariance",
            "projected_pair_standard_error": pair_se,
            "projected_H4_amplitude_standard_error": h4_se,
            "projected_A_scalar_standard_error": scalar_se,
            "projected_pair_z_at_H4_target": abs(pair_target) / pair_se,
            "required_samples_for_5sigma_pair_vs_zero": required_5sigma,
            "A_scalar_5sigma_sensitivity": 5.0 * scalar_se,
            "why_8m": (
                "exceeds the covariance-projected 5-sigma requirement for the parameter-free "
                "pair sign-flip versus scalar zero contrast"
            ),
        },
        "scoring_contract": {
            "primary": [
                "full two-component H4 predictive quadratic with model-fit covariance",
                "Delta_K_A against H4 target and zero",
                "H4_amplitude and A_projective_scalar with full transformed covariance",
            ],
            "no_H4_H8_vote": True,
            "no_exponent_fit": True,
            "no_post_reveal_basis_change": True,
        },
        "claim_boundary": (
            "N170 held-out angle-flip test only; distinguishes geometry direction from "
            "scalar/common-mode and localizes residual curvature without reopening harmonic voting"
        ),
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    prediction = payload["frozen_H4_only_prediction"]
    power = payload["power"]
    lines = [
        "# Preregistration: N170 exact angle-flip natural current", "",
        "Frozen before generating or inspecting any N170 projective-birth data.", "",
        "N170 `(11+7i,13+i)` is the exact `1+i` child of N85. The reflection-even H4 covector flips sign for both orientations, while the charged scalar `q_A^2=(u+H_F3)/2` remains fixed.", "",
        "The H4-only model from commit `186d72a` freezes:", "",
        f"- `K_A(first)={prediction['absolute_K_A'][0]:+.12g}`",
        f"- `K_A(second)={prediction['absolute_K_A'][1]:+.12g}`",
        f"- `Delta_K_A(second-first)={prediction['pair_second_minus_first']:+.12g}`",
        f"- H4 amplitude `{prediction['H4_amplitude']:+.12g}`",
        "- charged/projective scalar `0`", "",
        "The primary parameter-free discriminator is the sign-flipped negative pair contrast versus scalar/common-mode zero. Residual in the H4 amplitude coordinate is called curvature; residual in the orthogonal A-scalar coordinate is called projective/common mode.", "",
        f"Design: 8M samples/shape, 80 aligned batches, seed `202608337170`, on idle Huawei `DevEnvC_HZsCM6`. N145 covariance projects pair SE `{power['projected_pair_standard_error']:.4g}` and expected H4-versus-zero `z={power['projected_pair_z_at_H4_target']:.3f}`.", "",
        "No H4/H8 vote, exponent fit, or post-reveal basis change is allowed.", "",
    ]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crosswalk", type=Path, required=True)
    parser.add_argument("--n145-score", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    payload = preregister(args.crosswalk, args.n145_score)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
