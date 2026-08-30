#!/usr/bin/env python3
"""Apply the frozen support-first P250 radius selector."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from norm5_chiral_fixedp_mc import contexts
from score_z5_charged_threepoint import covariance_precision
from z5_charged_multiseparation_mc import SCHEMA, SEPARATIONS


CANDIDATE_ORDER = ("R1", "R2", "R3", "R4")
FIXED_RUN = {
    "samples": 4000,
    "batches": 40,
    "p": 0.59274605079,
    "seed": 25011312220260901,
    "replica_offset": 0,
    "replica_last_exclusive": 4000,
}


def read_batches(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def annulus_injectivity(radius: int) -> dict:
    plus, minus, _ = contexts()
    points = [
        (x, y)
        for y in range(-radius, radius + 1)
        for x in range(-radius, radius + 1)
        if (x, y) != (0, 0)
    ]
    hands = {}
    for name, context in (("plus", plus), ("minus", minus)):
        images = [context.geometry.vertex(point) for point in points]
        hands[name] = {
            "points": len(points),
            "unique_images": len(set(images)),
            "injective": len(set(images)) == len(points),
        }
    return {"radius": radius, "hands": hands, "passed": all(row["injective"] for row in hands.values())}


def covariance_nondegenerate(covariance: Sequence[Sequence[float]]) -> tuple[bool, str]:
    if len(covariance) != 8 or any(len(row) != 8 for row in covariance):
        return False, "not_8x8"
    diagonal = [float(covariance[index][index]) for index in range(8)]
    if any(not math.isfinite(value) or value <= 0.0 for value in diagonal):
        return False, "nonpositive_or_nonfinite_diagonal"
    try:
        covariance_precision(covariance)
    except ValueError as error:
        return False, f"singular_correlation:{error}"
    return True, "positive_diagonal_and_invertible_correlation"


def common_counter_certificate(candidates: Mapping[str, dict]) -> dict:
    reference = candidates[CANDIDATE_ORDER[0]]
    reference_rows = reference["batches"]
    failures = []
    for name in CANDIDATE_ORDER:
        row = candidates[name]
        run = row["response"]["run"]
        for key, expected in FIXED_RUN.items():
            if run.get(key) != expected:
                failures.append(f"{name}:run:{key}")
        if run.get("radius") != row["radius"]:
            failures.append(f"{name}:run:radius")
        if len(row["batches"]) != len(reference_rows):
            failures.append(f"{name}:batch_count")
            continue
        for index, (observed, frozen) in enumerate(zip(row["batches"], reference_rows)):
            for key in ("batch", "replica_first", "samples", "field_sha256", "translation_sha256"):
                if observed[key] != frozen[key]:
                    failures.append(f"{name}:batch{index}:{key}")
    return {
        "reference": "R1",
        "checked_fields": ["batch", "replica_first", "samples", "field_sha256", "translation_sha256"],
        "failures": failures,
        "passed": not failures,
    }


def select(candidates: Mapping[str, dict]) -> dict:
    if tuple(candidates) != CANDIDATE_ORDER:
        raise ValueError("candidate order differs from frozen R1,R2,R3,R4")
    certificate = common_counter_certificate(candidates)
    if not certificate["passed"]:
        raise ValueError("common-counter certificate failed")
    rows = {}
    for name in CANDIDATE_ORDER:
        candidate = candidates[name]
        radius = candidate["radius"]
        response = candidate["response"]
        score = candidate["score"]
        if response.get("schema") != SCHEMA or not response["mapping_gate"]["passed"]:
            raise ValueError(f"{name}: response gate failed")
        injection = annulus_injectivity(radius)
        if not injection["passed"]:
            raise ValueError(f"{name}: local annulus is not injective")
        separations = {}
        usable_count = 0
        for separation in SEPARATIONS:
            observed = score["separations"][str(separation)]
            nondegenerate, reason = covariance_nondegenerate(
                observed["local_variance_normalized_covariance"]
            )
            two_point_ready = bool(observed["two_point_ready_abs_z_ge_2"])
            usable = two_point_ready and nondegenerate
            usable_count += int(usable)
            support = observed["cubic_support_zero_score"]
            separations[str(separation)] = {
                "minimum_two_point_abs_z": min(
                    value["abs_z"] for value in observed["two_point_denominators"].values()
                ),
                "all_four_two_point_denominators_ready": two_point_ready,
                "cubic_covariance_nondegenerate": nondegenerate,
                "cubic_covariance_reason": reason,
                "usable_for_production_selector": usable,
                "descriptive_cubic_support_chi_square": support["chi_square"],
                "descriptive_cubic_support_df": support["degrees_of_freedom"],
                "descriptive_cubic_support_p": support["survival_p"],
                "phase_score_used_for_selection": False,
            }
        rows[name] = {
            "radius": radius,
            "annulus_injectivity": injection,
            "usable_separation_count": usable_count,
            "candidate_pass_at_least_two": usable_count >= 2,
            "separations": separations,
        }
    advancing = [name for name in CANDIDATE_ORDER if rows[name]["candidate_pass_at_least_two"]]
    selected = advancing[0] if advancing else None
    return {
        "schema": "matching-one/p250-radius-selector-result/v1",
        "status": "common_counter_geometry_selector_complete",
        "common_counter_certificate": certificate,
        "selection_rule": {
            "candidate_pass": "at least two separations have all four abs(G z)>=2 and nondegenerate cubic covariance",
            "tie_break": "first in frozen R1,R2,R3,R4 order",
            "phase_policy": "phase scores ignored",
        },
        "candidates": rows,
        "advancing_candidates": advancing,
        "selected_candidate": selected,
        "decision": (
            f"advance_{selected}" if selected else
            "no_production_candidate_in_radius_1_to_4_local_landing_h4_family"
        ),
        "next_selector": (
            "freeze a different charged insertion, such as a leg-defect or mesoscopic row, before acquiring more cubic replicas"
            if selected is None else "freeze a fresh production stream for the selected radius"
        ),
        "claim_boundary": [
            "The scan is a 4k common-counter geometry selector, not a production result.",
            "A small cubic support p-value at an unresolved two-point separation cannot select a candidate.",
            "No phase statistic was used in the selection.",
            "The result concerns only radii 1 through 4 of the current local landing-H4 insertion.",
        ],
    }


def render(result: dict) -> str:
    lines = [
        "# P250 common-counter radius selector", "",
        "Primary gate: all four two-point denominators plus a nondegenerate cubic covariance at at least two separations.", "",
        "| radius | d=1 min z / usable | d=2 min z / usable | d=3 min z / usable | pass |",
        "|---:|---:|---:|---:|---|",
    ]
    for name, row in result["candidates"].items():
        cells = []
        for separation in ("1", "2", "3"):
            item = row["separations"][separation]
            cells.append(f"{item['minimum_two_point_abs_z']:.3g} / {item['usable_for_production_selector']}")
        lines.append(f"| {name} | {cells[0]} | {cells[1]} | {cells[2]} | {row['candidate_pass_at_least_two']} |")
    lines.extend([
        "", f"Decision: `{result['decision']}`.", "",
        "The R4/d3 cubic zero score is descriptive only: its two-point denominator fails the frozen gate, so it cannot select R4.", "",
        f"Next selector: {result['next_selector']}", "",
    ])
    return "\n".join(lines)


def parse_candidate(value: str) -> tuple[str, int, Path, Path, Path]:
    fields = value.split(",")
    if len(fields) != 5:
        raise argparse.ArgumentTypeError("candidate must be NAME,RADIUS,RESPONSE,BATCHES,SCORE")
    return fields[0], int(fields[1]), Path(fields[2]), Path(fields[3]), Path(fields[4])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", action="append", type=parse_candidate, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    candidates = {}
    for name, radius, response_path, batches_path, score_path in args.candidate:
        candidates[name] = {
            "radius": radius,
            "response": json.loads(response_path.read_text()),
            "batches": read_batches(batches_path),
            "score": json.loads(score_path.read_text()),
        }
    result = select(candidates)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
