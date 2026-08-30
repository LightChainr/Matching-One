#!/usr/bin/env python3
"""Power gate for the Issue 158 norm-10 Gaussian root character.

This is a source-only calculation.  It reads the revealed P45 parent score and
P57 child variance, but never reads N=650/N=850 target data.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from statistics import NormalDist
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_P45 = ROOT / "results/server-20260828/P45-root-amplitude/score.json"
DEFAULT_P57 = ROOT / "results/server-20260829/P57-norm5-500m/primary_score.json"
DEFAULT_METADATA = {
    325: ROOT / "results/server-20260829/P57-norm5-500m/raw/n325_500m.metadata.json",
    425: ROOT / "results/server-20260829/P57-norm5-500m/raw/n425_500m.metadata.json",
}

Q = 10
ROOT_FACTOR = Fraction(7, 1250)
DELTA_COS4_FACTOR = ROOT_FACTOR * Q * Q
FIXED_P_FACTOR = float(DELTA_COS4_FACTOR) * Q ** (-13 / 8)
SLOPE_FACTOR = Q ** (3 / 8)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def required_child_samples(
    signal: float,
    child_variance_coefficient: float,
    source_variance: float,
    z_design: float,
) -> float:
    """Target samples for a Gaussian known-direction contrast.

    ``child_variance_coefficient`` is ``n * Var(child estimate)`` and
    ``source_variance`` is the already-acquired parent contribution after
    multiplication by the character factor.
    """

    budget = (abs(signal) / z_design) ** 2 - source_variance
    if budget <= 0:
        return math.inf
    return child_variance_coefficient / budget


def _quadratic_2x2(values: Sequence[float], covariance: Sequence[Sequence[float]]) -> float:
    a, b = covariance[0]
    c, d = covariance[1]
    determinant = a * d - b * c
    _require(determinant > 0, "covariance must be positive definite")
    return (
        d * values[0] ** 2
        - (b + c) * values[0] * values[1]
        + a * values[1] ** 2
    ) / determinant


def joint_equal_depth_samples(
    signals: Sequence[float],
    child_coefficients: Sequence[float],
    source_covariance: Sequence[Sequence[float]],
    z_design: float,
) -> float:
    """Equal target depth needed for the two-lineage matched-direction score."""

    target = z_design * z_design

    def information(samples: float) -> float:
        covariance = [list(row) for row in source_covariance]
        covariance[0][0] += child_coefficients[0] / samples
        covariance[1][1] += child_coefficients[1] / samples
        return _quadratic_2x2(signals, covariance)

    low, high = 1.0, 1.0
    while information(high) < target and high < 1e18:
        high *= 2
    if high >= 1e18 and information(high) < target:
        return math.inf
    for _ in range(120):
        middle = (low + high) / 2
        if information(middle) < target:
            low = middle
        else:
            high = middle
    return high


def build_plan(
    p45: Mapping[str, Any],
    p57: Mapping[str, Any],
    metadata: Mapping[int, Mapping[str, Any]],
    *,
    alpha: float = 0.01,
    power: float = 0.80,
) -> dict[str, Any]:
    _require(p45.get("schema") == "frozen angular-normalized root amplitude score v1", "wrong P45 score")
    _require(p57.get("schema") == "norm5 frozen harmonic primary score v1", "wrong P57 score")
    _require(p57.get("size_order") == [65, 85, 325, 425], "P57 size order drift")
    _require(0 < alpha < 1 and 0.5 < power < 1, "invalid design probability")
    _require(math.isclose(FIXED_P_FACTOR / SLOPE_FACTOR, float(ROOT_FACTOR), rel_tol=1e-14), "slope conversion drift")

    z_design = NormalDist().inv_cdf(1 - alpha / 2) + NormalDist().inv_cdf(power)
    lineages = []
    child_coefficients = []
    signals = []
    source_covariance_delta = [[0.0, 0.0], [0.0, 0.0]]
    parent_sizes = (65, 85)
    variance_proxy_sizes = (325, 425)
    observation_covariance = [
        [float(value) for value in row] for row in p57["observation_covariance"]
    ]

    for i, (parent_n, proxy_n) in enumerate(zip(parent_sizes, variance_proxy_sizes)):
        parent = p45["by_size"][str(parent_n)]
        estimate = parent["estimate"]
        errors = parent["standard_errors"]
        samples = int(metadata[proxy_n]["samples_per_pair"])
        elapsed = float(metadata[proxy_n]["elapsed_seconds"])
        child_delta_variance_coefficient = observation_covariance[i + 2][i + 2] * samples
        predicted_delta_m = FIXED_P_FACTOR * float(estimate["delta_M"])
        predicted_slope = SLOPE_FACTOR * float(estimate["mean_M_prime"])
        predicted_root_gap = float(ROOT_FACTOR) * float(estimate["root_gap"])
        converted_root_gap = -predicted_delta_m / predicted_slope
        root_variance_coefficient = child_delta_variance_coefficient / predicted_slope**2
        source_delta_variance = (FIXED_P_FACTOR * float(errors["delta_M"])) ** 2
        source_root_variance = (float(ROOT_FACTOR) * float(errors["root_gap"])) ** 2
        delta_samples = required_child_samples(
            predicted_delta_m, child_delta_variance_coefficient,
            source_delta_variance, z_design,
        )
        root_samples = required_child_samples(
            predicted_root_gap, root_variance_coefficient,
            source_root_variance, z_design,
        )
        relative_precision_samples = {}
        for relative in (0.20, 0.10, 0.05):
            required = required_child_samples(
                predicted_delta_m, child_delta_variance_coefficient,
                source_delta_variance, 1 / relative,
            )
            relative_precision_samples[str(relative)] = None if math.isinf(required) else required

        parent_delta_coefficient = float(errors["delta_M"]) ** 2 * int(p45["provenance"]["samples_per_pair"])
        perfect_crn_lower_coefficient = (
            math.sqrt(child_delta_variance_coefficient)
            - abs(FIXED_P_FACTOR) * math.sqrt(parent_delta_coefficient)
        ) ** 2
        perfect_crn_reduction = 1 - perfect_crn_lower_coefficient / child_delta_variance_coefficient
        perfect_crn_samples = perfect_crn_lower_coefficient / (abs(predicted_delta_m) / z_design) ** 2

        proxy_rate = samples / elapsed
        target_n = parent_n * Q
        target_rate_8_thread = proxy_rate * proxy_n / target_n
        lineages.append({
            "parent_N": parent_n,
            "target_N": target_n,
            "variance_proxy_N": proxy_n,
            "leading_targets": {
                "delta_M_at_fixed_p": predicted_delta_m,
                "mean_M_prime": predicted_slope,
                "root_gap": predicted_root_gap,
                "root_gap_from_minus_delta_M_over_slope": converted_root_gap,
                "root_linearization_relative_mismatch": abs(converted_root_gap / predicted_root_gap - 1),
            },
            "per_replica_variance_coefficients": {
                "delta_M": child_delta_variance_coefficient,
                "root_gap_via_slope": root_variance_coefficient,
            },
            "archived_parent_standard_errors_after_transport": {
                "delta_M": math.sqrt(source_delta_variance),
                "root_gap": math.sqrt(source_root_variance),
            },
            "no_parent_child_crn": {
                "target_samples_from_delta_M": delta_samples,
                "target_samples_from_root_gap": root_samples,
                "expected_z_at_500M_from_delta_M": abs(predicted_delta_m) / math.sqrt(
                    child_delta_variance_coefficient / 500_000_000 + source_delta_variance
                ),
                "target_samples_for_relative_standard_error": relative_precision_samples,
                "estimated_8_thread_hours_from_P57_inverse_N_rate": delta_samples / target_rate_8_thread / 3600,
            },
            "unattainable_perfect_parent_child_crn_bound": {
                "maximum_variance_reduction_fraction": perfect_crn_reduction,
                "equal_depth_samples_lower_bound": perfect_crn_samples,
                "assumption": "fresh parent rerun, equal replica depth, absolute correlation one",
            },
        })
        child_coefficients.append(child_delta_variance_coefficient)
        signals.append(predicted_delta_m)

    for i in range(2):
        for j in range(2):
            source_covariance_delta[i][j] = FIXED_P_FACTOR**2 * observation_covariance[i][j]
    joint_samples = joint_equal_depth_samples(
        signals, child_coefficients, source_covariance_delta, z_design,
    )
    return {
        "schema": "matching-one/gaussian-commuting-square-power-gate/v1",
        "status": "source_only_power_no_target_data",
        "design": {
            "two_sided_alpha": alpha,
            "power": power,
            "normal_z_sum": z_design,
            "interpretation": "power to distinguish the leading character target from a zero child contrast, not precision on the exact rational ratio",
        },
        "character_conversion": {
            "norm": Q,
            "delta_cos4_child_over_parent": str(DELTA_COS4_FACTOR),
            "fixed_p_delta_M_child_over_parent": FIXED_P_FACTOR,
            "slope_child_over_parent": SLOPE_FACTOR,
            "root_child_over_parent": str(ROOT_FACTOR),
            "identity": "root_factor = fixed_p_factor / slope_factor",
        },
        "lineages": lineages,
        "joint_equal_target_depth": {
            "samples_per_lineage": joint_samples,
            "total_target_samples": 2 * joint_samples,
            "score": "one known-direction contrast against zero using both independent target streams and archived parent covariance",
            "warning": "this does not establish per-lineage commuting-square closure",
        },
        "runner_reuse": {
            "engine": "src/threshold_rank_integer_period_mc.cpp",
            "N650_matrices": [[[23, -11], [11, 23]], [[19, -17], [17, 19]]],
            "N850_matrices": [[[29, -3], [3, 29]], [[27, -11], [11, 27]]],
            "same_N_orientation_crn_already_present": True,
            "parent_child_cover_crn_present": False,
            "direct_and_two_step_paths_are_the_same_target_quotient": True,
            "scorer_gap": "score_angular_root_amplitude.py hard-codes [65,85]; parameterize sizes/provenance before production",
        },
        "crn_conclusion": {
            "classification": "power_no_go_for_new_norm10_covering_crn",
            "reason": "the transported parent coefficient is only 0.01328 in fixed-p M (0.0056 in roots), so Cauchy-Schwarz caps even perfect fresh parent-child CRN at less than 2.4 percent variance reduction",
            "minimum_new_coupling_if_still_requested": "ten-fiber child uniforms plus an audited symmetric measure-preserving parent rank map, with joint batch moments; pilot only, because the absolute-correlation bound already rules out material gain",
        },
        "source_paths": {
            "P45": str(DEFAULT_P45.relative_to(ROOT)),
            "P57": str(DEFAULT_P57.relative_to(ROOT)),
            "P57_metadata": {str(key): str(path.relative_to(ROOT)) for key, path in DEFAULT_METADATA.items()},
        },
        "contains_target_data": False,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p45", type=Path, default=DEFAULT_P45)
    parser.add_argument("--p57", type=Path, default=DEFAULT_P57)
    parser.add_argument("--alpha", type=float, default=0.01)
    parser.add_argument("--power", type=float, default=0.80)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args(argv)
    p45 = json.loads(args.p45.read_text(encoding="utf-8"))
    p57 = json.loads(args.p57.read_text(encoding="utf-8"))
    metadata = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in DEFAULT_METADATA.items()
    }
    result = build_plan(p45, p57, metadata, alpha=args.alpha, power=args.power)
    payload = json.dumps(result, allow_nan=False, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
