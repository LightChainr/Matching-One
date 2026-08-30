#!/usr/bin/env python3
"""Score only propagation and cubic-interface gates for the P250 leg row."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

from score_z5_charged_multiseparation import (
    denominator_rows,
    jackknife,
    read_batches,
)
from score_z5_charged_threepoint import covariance_precision, zero_score
from z5_charged_multiseparation_mc import SEPARATIONS
from z5_projective_leg_multiseparation_mc import SCHEMA


def covariance_nondegenerate(covariance: Sequence[Sequence[float]]) -> tuple[bool, str]:
    diagonal = [float(covariance[index][index]) for index in range(len(covariance))]
    if len(covariance) != 8 or any(len(row) != 8 for row in covariance):
        return False, "not_8x8"
    if any(not math.isfinite(value) or value <= 0.0 for value in diagonal):
        return False, "nonpositive_or_nonfinite_diagonal"
    try:
        covariance_precision(covariance)
    except ValueError as error:
        return False, f"singular_correlation:{error}"
    return True, "positive_diagonal_and_invertible_correlation"


def score(payload: dict, batches: Sequence[dict]) -> dict:
    if payload.get("schema") != SCHEMA or not payload["exact_gate"]["passed"]:
        raise ValueError("wrong or failed projective-leg response")
    separations = {}
    ready_count = 0
    for separation in SEPARATIONS:
        denominators = denominator_rows(batches, separation)
        pair_ready = all(value["abs_z"] >= 2.0 for value in denominators.values())
        local_point, local_covariance, _ = jackknife(
            batches, separation, "local_variance"
        )
        nondegenerate, reason = covariance_nondegenerate(local_covariance)
        support = zero_score(local_point, local_covariance) if nondegenerate else None
        usable = pair_ready and nondegenerate
        ready_count += int(usable)
        separations[str(separation)] = {
            "two_point_denominators": denominators,
            "minimum_two_point_abs_z": min(value["abs_z"] for value in denominators.values()),
            "all_four_two_point_denominators_ready": pair_ready,
            "local_variance_normalized_cubic_point": local_point,
            "local_variance_normalized_cubic_covariance": local_covariance,
            "cubic_interface_nondegenerate": nondegenerate,
            "cubic_interface_reason": reason,
            "descriptive_cubic_support_zero_score": support,
            "usable_for_promotion": usable,
            "phase_score_computed": False,
        }
    return {
        "schema": "matching-one/z5-projective-leg-propagation-score/v1",
        "status": "frozen_low_sample_operator_propagation_smoke",
        "primary_insertion": payload["observable"]["primary_insertion"],
        "usable_separation_count": ready_count,
        "promotion_gate_at_least_two": ready_count >= 2,
        "separations": separations,
        "decision": (
            "projective_leg_operator_passes_low_sample_propagation_gate"
            if ready_count >= 2 else
            "projective_leg_operator_does_not_pass_low_sample_propagation_gate"
        ),
        "next_selector": (
            "freeze a fresh production-sized support-first block without changing the operator"
            if ready_count >= 2 else
            "do not increase samples until a different root-marked operator is defined"
        ),
        "claim_boundary": [
            "The 2k stream tests operator propagation and interface rank only.",
            "No phase statistic was computed or used.",
            "Passing does not identify a continuum field or OPE coefficient.",
        ],
    }


def render(result: dict) -> str:
    lines = [
        "# P250 projective-leg propagation smoke", "",
        "No phase score is computed.  The only promotion gate is resolved charged two-point propagation plus a nondegenerate cubic interface.", "",
        "| d | minimum pair z | all four pairs ready | cubic interface | usable |",
        "|---:|---:|---|---|---|",
    ]
    for separation, row in result["separations"].items():
        lines.append(
            f"| {separation} | {row['minimum_two_point_abs_z']:.3g} | "
            f"{row['all_four_two_point_denominators_ready']} | "
            f"{row['cubic_interface_nondegenerate']} | {row['usable_for_promotion']} |"
        )
    lines.extend([
        "", f"Decision: `{result['decision']}`.", "",
        f"Next selector: {result['next_selector']}", "",
    ])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("response", type=Path)
    parser.add_argument("batches", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = score(json.loads(args.response.read_text()), read_batches(args.batches))
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
