#!/usr/bin/env python3
"""Freeze the N340 second-child charged-current heldout production."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Mapping, Optional, Sequence


AREA_POWER = -13.0 / 8.0
TARGET_N = 340
SAMPLES = 12_000_000
BATCHES = 80
SEED = 202608337340


def chi4_real(a: int, b: int) -> Fraction:
    n = a * a + b * b
    return Fraction(a**4 - 6 * a * a * b * b + b**4, n * n)


def projected_variance(weights: Sequence[float], covariance: Sequence[Sequence[float]]) -> float:
    return math.fsum(
        weights[i] * covariance[i][j] * weights[j]
        for i in range(len(weights)) for j in range(len(weights))
    )


def preregister(n85_path: Path, n170_path: Path) -> dict[str, object]:
    n85 = json.loads(n85_path.read_text(encoding="utf-8"))
    n170 = json.loads(n170_path.read_text(encoding="utf-8"))
    if n85.get("schema") != "matching-one/P337-natural-current-scale-N85-score/v1":
        raise ValueError("unexpected N85 score schema")
    if n170.get("schema") != "matching-one/P337-N170-angle-flip-score/v1":
        raise ValueError("unexpected N170 score schema")

    parent = ((11, 7), (13, 1))
    child = ((18, 4), (14, 12))
    parent_covectors = [chi4_real(*row) for row in parent]
    covectors = [chi4_real(*row) for row in child]
    if any(child_value != -parent_value
           for child_value, parent_value in zip(covectors, parent_covectors)):
        raise AssertionError("N340 is not the exact second angle-flip child of N170")
    dc = float(covectors[1] - covectors[0])
    transform = [
        [-1.0 / dc, 1.0 / dc],
        [2.0 * float(covectors[1]) / dc, -2.0 * float(covectors[0]) / dc],
    ]

    n85_values = n85["natural_coordinate"]["value"][:2]
    n85_covariance = [row[:2] for row in n85["natural_coordinate"]["covariance"][:2]]
    n85_covectors = [chi4_real(9, 2), chi4_real(7, 6)]
    n85_dc = float(n85_covectors[1] - n85_covectors[0])
    n85_weights = [-1.0 / n85_dc, 1.0 / n85_dc]
    amplitude_85 = math.fsum(a * b for a, b in zip(n85_weights, n85_values))
    amplitude_85_variance = projected_variance(n85_weights, n85_covariance)

    split170 = n170["curvature_projective_decomposition"]
    amplitude_170 = split170["observed"][0]
    amplitude_170_variance = split170["measurement_covariance"][0][0]
    effective_ratio = amplitude_170 / amplitude_85
    nominal_ratio = 2.0 ** AREA_POWER
    target_amplitudes = {
        "nominal_area_H4": nominal_ratio * amplitude_170,
        "observed_N85_to_N170_effective": amplitude_170 * effective_ratio,
        "scale_neutral": amplitude_170,
    }
    target_variances = {
        "nominal_area_H4": nominal_ratio * nominal_ratio * amplitude_170_variance,
        "observed_N85_to_N170_effective": (
            (2.0 * effective_ratio) ** 2 * amplitude_170_variance
            + effective_ratio**4 * amplitude_85_variance
        ),
        "scale_neutral": amplitude_170_variance,
    }

    models = []
    for name in ("nominal_area_H4", "observed_N85_to_N170_effective", "scale_neutral"):
        amplitude = target_amplitudes[name]
        variance = target_variances[name]
        vector = [amplitude * float(c) for c in covectors]
        vector_covariance = [
            [variance * float(left) * float(right) for right in covectors]
            for left in covectors
        ]
        models.append({
            "name": name,
            "role": (
                "primary asymptotic H4 target" if name == "nominal_area_H4" else
                "secondary descriptive two-state/slow-state continuation" if name == "observed_N85_to_N170_effective" else
                "fixed scale-neutral comparator"
            ),
            "H4_amplitude": amplitude,
            "H4_amplitude_target_variance": variance,
            "absolute_K_A": vector,
            "absolute_target_covariance": vector_covariance,
            "pair_second_minus_first": amplitude * dc,
            "pair_sign": "positive",
        })

    measurement_var_projection = amplitude_170_variance * 8_000_000 / SAMPLES
    gaps = {
        "nominal_vs_effective": abs(target_amplitudes["nominal_area_H4"] - target_amplitudes["observed_N85_to_N170_effective"]),
        "effective_vs_neutral": abs(target_amplitudes["observed_N85_to_N170_effective"] - target_amplitudes["scale_neutral"]),
        "nominal_vs_neutral": abs(target_amplitudes["nominal_area_H4"] - target_amplitudes["scale_neutral"]),
    }
    return {
        "schema": "matching-one/P337-N340-second-child-preregistration/v1",
        "status": "frozen before any N340 projective-birth generation or inspection",
        "source": {
            "N85_score": str(n85_path),
            "N170_score": str(n170_path),
            "N170_reveal_commit": "0db21b7",
            "p_ref": n170["source"]["p_ref"],
            "independence": "N85 and N170 are independent seed/counter blocks",
        },
        "geometry": {
            "N": TARGET_N,
            "first": [18, 4], "second": [14, 12],
            "first_period_matrix": [[18, -4], [4, 18]],
            "second_period_matrix": [[14, -12], [12, 14]],
            "raw_multiplier_images": [[4, 18], [12, 14]],
            "canonicalization": "D4-equivalent canonical representatives used by the runner",
            "parent_N170": [list(row) for row in parent],
            "common_multiplier": [1, 1],
            "H4_covectors_exact": [str(value) for value in covectors],
            "parent_H4_covectors_exact": [str(value) for value in parent_covectors],
            "exact_angle_flip": True,
            "Smith_classes": [[2, 170], [2, 170]],
            "charged_projective_scalar": {
                "value": 0.5,
                "frozen_amplitude_target": 0.0,
                "transport": "unchanged along the Gaussian 1+i lineage",
            },
            "F3_reduction_nondegenerate": TARGET_N % 3 != 0,
        },
        "observable": {
            "definition": "K_A=p(1-p)Jminus_A/W_A=d_eta log W_A",
            "vector_order": ["K_A_first_18_plus_4i", "K_A_second_14_plus_12i"],
            "decomposition": {
                "H4_amplitude": "(K_second-K_first)/(c_second-c_first)",
                "A_projective_scalar": "2(c_second*K_first-c_first*K_second)/(c_second-c_first)",
                "transform": transform,
            },
        },
        "frozen_models_in_scoring_order": models,
        "frozen_source_amplitudes": {
            "N85_H4_amplitude": amplitude_85,
            "N85_measurement_variance": amplitude_85_variance,
            "N170_H4_amplitude": amplitude_170,
            "N170_measurement_variance": amplitude_170_variance,
            "nominal_area_ratio": nominal_ratio,
            "observed_N85_to_N170_ratio": effective_ratio,
            "ratio_is_secondary_descriptive": True,
        },
        "production": {
            "samples_per_shape": SAMPLES, "batches": BATCHES, "seed": SEED,
            "replica_offset": 0, "threads": 16,
            "engine": "1714141 projective-birth integer-period runner",
            "machine": {
                "name": "DevEnvC_HZsCM6",
                "id": "033945d8bf8b47a7acf475c595169e07",
                "class": "Kunpeng AArch64, 16 vCPU, 32 GiB",
            },
            "full_batch_covariance_required": True,
        },
        "power": {
            "variance_source": "N170 8M/shape H4-amplitude covariance",
            "projected_measurement_standard_error": math.sqrt(measurement_var_projection),
            "fixed_amplitude_gaps": gaps,
            "projected_measurement_z_for_gaps": {
                key: value / math.sqrt(measurement_var_projection) for key, value in gaps.items()
            },
            "why_12m": "resolves nominal versus effective at >5 SE and effective versus neutral at >4 SE under N170 variance",
        },
        "scoring_contract": {
            "primary_coordinate": "H4_amplitude",
            "fixed_model_order": [model["name"] for model in models],
            "projective_scalar_zero_control": True,
            "report_measurement_only_and_source_uncertainty_aware_residuals": True,
            "no_exponent_fit": True,
            "no_post_reveal_model_or_basis_change": True,
            "no_H4_H8_vote": True,
        },
        "claim_boundary": "one heldout second child in the same Gaussian lineage; discriminates fixed scale transfers and scalar leakage, not a fitted continuum exponent",
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    models = payload["frozen_models_in_scoring_order"]
    source = payload["frozen_source_amplitudes"]
    power = payload["power"]
    lines = [
        "# Preregistration: N340 same-lineage second child", "",
        "Frozen before generating or inspecting N340 data. Multiplication of N170 `(11+7i,13+i)` by `1+i`, followed only by D4 canonicalization, gives N340 `(18+4i,14+12i)`. The exact H4 covector flips again; the projective scalar remains a frozen zero-amplitude control.", "",
        "Primary coordinate: `A_H=(K_second-K_first)/(c_second-c_first)`. Fixed targets:", "",
    ]
    for model in models:
        lines.append(f"- `{model['name']}`: `A_H={model['H4_amplitude']:+.12g}`, pair `{model['pair_second_minus_first']:+.12g}`")
    lines.extend([
        "",
        f"The secondary effective transfer is the already revealed N85-to-N170 ratio `{source['observed_N85_to_N170_ratio']:.9g}`; it is not an exponent fit. N85 target uncertainty is retained.", "",
        f"Design: 12M/shape, 80 aligned batches, seed `202608337340`, HZsCM6. N170 variance projects H4-amplitude SE `{power['projected_measurement_standard_error']:.4g}`; fixed target gaps are {power['projected_measurement_z_for_gaps']['nominal_vs_effective']:.2f}, {power['projected_measurement_z_for_gaps']['effective_vs_neutral']:.2f}, and {power['projected_measurement_z_for_gaps']['nominal_vs_neutral']:.2f} measurement SE.", "",
        "Scoring order is frozen: nominal area H4, observed effective continuation, scale-neutral. Both measurement-only and source-uncertainty-aware residuals are reported. No exponent fit, basis change, or H4/H8 revote is allowed.", "",
    ])
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n85-score", type=Path, required=True)
    parser.add_argument("--n170-score", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    payload = preregister(args.n85_score, args.n170_score)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
