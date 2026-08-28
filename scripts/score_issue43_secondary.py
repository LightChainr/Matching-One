#!/usr/bin/env python3
"""Score the frozen Issue #43 secondary hypotheses without reading run raw data.

The input from Issue #43 is the JSON already produced by
``score_issue43_full_curve.py``.  This script deliberately has no histogram,
moments, or metadata arguments: the primary scorer remains the only consumer
of N=185/265 production sufficient statistics.

The fixed reporting order is:

1. reference the already-scored original x=21/4 two-sector primary;
2. score the frozen x=17/4 DeltaM competitor, without refitting;
3. reference the already-scored zero benchmark;
4. report the predeclared shared H4+H12 comparison as NOT SCORABLE unless a
   pre-target amplitude vector and its source covariance exist (none is
   frozen in the repository artifact set at this commit);
5. only then, optionally score the #72 P48 S-prime frozen models.

The invalidated wrong-Kac-branch V_<1,3> artifact is never scored.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import score_p48_sprime_prospective as p48


TARGET_SIZES = (185, 265)
PRIMARY_PREDICTION_SHA256 = (
    "a370e79a10854341fac3ee75e8c518dbf3533e8c077cba2c2ec1018178144f44"
)
X17_ARTIFACT_SHA256 = (
    "941af010cc146c76e26985ecf3edf58f0df28d987fc79c03725ebc21f64964f5"
)
P48_ARTIFACT_SHA256 = (
    "0d44228ae117f94cb1f99d1e2727eb47390aae950c3ae70c21dd8bc5a09454ae"
)

# Frozen in predictions/x17_spin4_competitor_20260828.yaml before target reveal.
X17_MEAN = (0.0002868415917648024, 0.0002799185461225267)
X17_SOURCE_SE = (0.000012974176398815644, 0.000012661039050681047)

# Exact pre-target harmonic audit columns.  They establish leverage only; they
# do not supply the missing frozen A12 coefficient or its covariance.
H12_OVER_H4 = (-1.0255132676434948, 0.1505079050060446)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("input must be a JSON object")
    return payload


def validate_primary(payload: Mapping[str, object]) -> tuple[list[float], list[float]]:
    if payload.get("protocol") != "Issue #43 prospective N=185/265 two-spin4 full-curve score":
        raise ValueError("input is not the frozen Issue #43 primary-score output")
    if payload.get("prediction_artifact_sha256") != PRIMARY_PREDICTION_SHA256:
        raise ValueError("primary frozen prediction hash mismatch")
    scores = payload.get("scores")
    if not isinstance(scores, dict) or "DeltaM" not in scores or "DeltaS" not in scores:
        raise ValueError("primary output lacks both DeltaM and DeltaS scores")
    delta_m = scores["DeltaM"]
    if not isinstance(delta_m, dict):
        raise ValueError("primary DeltaM score is malformed")
    observed = [float(value) for value in delta_m.get("observed", ())]
    sampling_se = [float(value) for value in delta_m.get("sampling_se", ())]
    if len(observed) != 2 or len(sampling_se) != 2:
        raise ValueError("primary DeltaM must contain N=185/265 values")
    if any(not math.isfinite(value) for value in observed):
        raise ValueError("primary DeltaM observations must be finite")
    if any(not math.isfinite(value) or value <= 0.0 for value in sampling_se):
        raise ValueError("primary DeltaM sampling errors must be finite and positive")
    return observed, sampling_se


def quadratic_2(vector: Sequence[float], covariance: Sequence[Sequence[float]]) -> float:
    a = float(covariance[0][0])
    b = float(covariance[0][1])
    d = float(covariance[1][1])
    determinant = a * d - b * b
    if determinant <= 0.0:
        raise ValueError("score covariance is not positive definite")
    x, y = map(float, vector)
    return (d * x * x - 2.0 * b * x * y + a * y * y) / determinant


def copied_primary_stage(primary: Mapping[str, object]) -> dict:
    scores = primary["scores"]
    return {
        "order": 1,
        "name": "original_x21_H4_two_sector",
        "status": "ALREADY_SCORED_BY_PRIMARY",
        "target_refit_parameters": 0,
        "DeltaM": {
            "chi_square": scores["DeltaM"]["target_chi_square"],
            "df": scores["DeltaM"]["target_df"],
        },
        "DeltaS": {
            "chi_square": scores["DeltaS"]["target_chi_square"],
            "df": scores["DeltaS"]["target_df"],
        },
    }


def score_x17(
    observed: Sequence[float], sampling_se: Sequence[float], artifact_path: Path
) -> dict:
    if sha256(artifact_path) != X17_ARTIFACT_SHA256:
        raise ValueError("frozen x=17/4 prediction artifact hash mismatch")
    residual = [float(observed[i]) - X17_MEAN[i] for i in range(2)]
    covariance = [
        [
            (float(sampling_se[i]) ** 2 if i == j else 0.0)
            + X17_SOURCE_SE[i] * X17_SOURCE_SE[j]
            for j in range(2)
        ]
        for i in range(2)
    ]
    chi_square = quadratic_2(residual, covariance)
    return {
        "order": 2,
        "name": "x17_over_4_H4_adversarial_radial",
        "status": "SCORED_FROZEN_NO_REFIT",
        "law": "DeltaM=A17*DeltaCos4*N^(-9/8)",
        "sizes": list(TARGET_SIZES),
        "observed": list(map(float, observed)),
        "sampling_se": list(map(float, sampling_se)),
        "frozen_mean": list(X17_MEAN),
        "source_amplitude_only_se": list(X17_SOURCE_SE),
        "source_error_correlation": [[1.0, 1.0], [1.0, 1.0]],
        "residual": residual,
        "residual_covariance": covariance,
        "marginal_signed_z": [
            residual[i] / math.sqrt(covariance[i][i]) for i in range(2)
        ],
        "chi_square": chi_square,
        "df": 2,
        "chi_square_survival_df2": math.exp(-0.5 * chi_square),
        "target_refit_parameters": 0,
        "prediction_artifact": str(artifact_path),
        "prediction_artifact_sha256": X17_ARTIFACT_SHA256,
    }


def copied_zero_stage(primary: Mapping[str, object]) -> dict:
    scores = primary["scores"]
    return {
        "order": 3,
        "name": "zero_effect",
        "status": "ALREADY_SCORED_BY_PRIMARY",
        "target_refit_parameters": 0,
        "DeltaM": {
            "chi_square": scores["DeltaM"]["zero_chi_square"],
            "df": scores["DeltaM"]["zero_df"],
        },
        "DeltaS": {
            "chi_square": scores["DeltaS"]["zero_chi_square"],
            "df": scores["DeltaS"]["zero_df"],
        },
    }


def h4_h12_stage() -> dict:
    return {
        "order": 4,
        "name": "predeclared_shared_H4_plus_H12",
        "status": "NOT_SCORABLE",
        "sizes": list(TARGET_SIZES),
        "exact_DeltaCos12_over_DeltaCos4": list(H12_OVER_H4),
        "reason": (
            "The pre-target record freezes exact H12 harmonic columns but no complete "
            "shared H4/H12 amplitude vector and no source-prediction covariance. "
            "Fitting A12 to N185/N265 here would be a forbidden target refit."
        ),
        "required_to_unlock": (
            "A pre-target artifact containing both frozen amplitudes, their full source "
            "covariance, radial powers, and an immutable hash."
        ),
        "target_refit_parameters": 0,
    }


def p48_stage(target_path: Path | None, artifact_path: Path) -> dict:
    base = {
        "order": 5,
        "name": "issue72_P48_S_prime",
        "model_order": [
            "pure_power_baseline",
            "rank2_log_primary_correction",
            "analytic_inverse_N_competitor",
            "zero_effect",
        ],
        "global_guard": "score only after Issue #43 orders 1-4 are reported",
        "target_refit_parameters_per_model": 0,
    }
    if sha256(artifact_path) != P48_ARTIFACT_SHA256:
        raise ValueError("frozen P48 prediction artifact hash mismatch")
    base["prediction_artifact"] = str(artifact_path)
    base["prediction_artifact_sha256"] = P48_ARTIFACT_SHA256
    if target_path is None:
        base.update({
            "status": "READY_AWAITING_DERIVATIVE_TARGET",
            "reason": "Supply the aggregated P4[S_prime] target JSON; do not supply production raw.",
        })
        return base
    result = p48.score(p48.read_json(target_path), p48.read_yaml(artifact_path))
    observed_order = [row["name"] for row in result["results"]]
    if observed_order != base["model_order"]:
        raise ValueError("P48 scorer model order differs from the frozen secondary order")
    base.update({"status": "SCORED_FROZEN_NO_REFIT", "score": result})
    return base


def score(
    primary: Mapping[str, object],
    x17_artifact: Path,
    p48_artifact: Path,
    p48_target: Path | None = None,
) -> dict:
    observed, sampling_se = validate_primary(primary)
    stages = [
        copied_primary_stage(primary),
        score_x17(observed, sampling_se, x17_artifact),
        copied_zero_stage(primary),
        h4_h12_stage(),
        p48_stage(p48_target, p48_artifact),
    ]
    return {
        "protocol": "Issue #43 frozen secondary scoring ledger",
        "status": "fixed secondary scores only; no target refit",
        "sizes": list(TARGET_SIZES),
        "stage_order": [stage["name"] for stage in stages],
        "stages": stages,
        "excluded_models": [
            {
                "name": "V_<1,3>_N^-4/3",
                "status": "EXCLUDED_INVALIDATED_WRONG_KAC_BRANCH",
                "scored": False,
            }
        ],
        "raw_data_boundary": (
            "This scorer consumes the primary score JSON and an optional aggregated "
            "P48 target JSON. It has no interface for N185/N265 raw histograms or moments."
        ),
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary-score", type=Path, required=True)
    parser.add_argument("--p48-target", type=Path)
    parser.add_argument(
        "--x17-artifact",
        type=Path,
        default=root / "predictions/x17_spin4_competitor_20260828.yaml",
    )
    parser.add_argument(
        "--p48-artifact",
        type=Path,
        default=root / "predictions/p48_sprime_correction_20260828.yaml",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = score(
        read_json(args.primary_score),
        args.x17_artifact,
        args.p48_artifact,
        args.p48_target,
    )
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
