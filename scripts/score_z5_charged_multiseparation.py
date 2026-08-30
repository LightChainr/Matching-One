#!/usr/bin/env python3
"""Support-first score for the P250 multi-separation charged cubic."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from score_z5_charged_threepoint import zero_score
from z5_charged_multiseparation_mc import HANDS, SCHEMA, SEPARATIONS


PRIMARY = ("C113_plus", "C113_minus", "C122_plus", "C122_minus")


def read_batches(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return [
            {
                key: int(value) if key in {"batch", "replica_first", "samples"} else value
                if key in {"field_sha256", "translation_sha256"}
                else float(value)
                for key, value in row.items()
            }
            for row in csv.DictReader(handle)
        ]


def _complex(means: Mapping[str, float], separation: int, channel: str, hand: str) -> complex:
    prefix = f"d{separation}_{channel}_{hand}"
    return complex(means[prefix + "_re"], means[prefix + "_im"])


def means(rows: Sequence[dict], excluded: int | None = None) -> dict[str, float]:
    kept = [row for index, row in enumerate(rows) if index != excluded]
    samples = sum(row["samples"] for row in kept)
    fields = [key for key in kept[0] if key.startswith("d")]
    return {field: sum(row[field] for row in kept) / samples for field in fields}


def normalized_vector(values: Mapping[str, float], separation: int, mode: str) -> list[float]:
    output = []
    for channel, hand in (("C113", "plus"), ("C113", "minus"),
                          ("C122", "plus"), ("C122", "minus")):
        g1 = abs(values[f"d{separation}_G1_{hand}"])
        g2 = abs(values[f"d{separation}_G2_{hand}"])
        v1 = values[f"d{separation}_V1_{hand}"]
        v2 = values[f"d{separation}_V2_{hand}"]
        if mode == "separation":
            denominator = g1 * math.sqrt(g2) if channel == "C113" else math.sqrt(g1) * g2
        elif mode == "local_variance":
            denominator = v1 * math.sqrt(v2) if channel == "C113" else math.sqrt(v1) * v2
        else:
            raise ValueError("unknown normalization")
        if denominator <= 0.0:
            raise ValueError("normalization denominator is not positive")
        value = _complex(values, separation, channel, hand) / denominator
        output.extend((value.real, value.imag))
    return output


def jackknife(rows: Sequence[dict], separation: int, mode: str):
    point = normalized_vector(means(rows), separation, mode)
    delete_one = [normalized_vector(means(rows, index), separation, mode) for index in range(len(rows))]
    center = [sum(row[j] for row in delete_one) / len(delete_one) for j in range(len(point))]
    factor = (len(rows) - 1) / len(rows)
    covariance = [
        [
            factor * sum((row[i] - center[i]) * (row[j] - center[j]) for row in delete_one)
            for j in range(len(point))
        ]
        for i in range(len(point))
    ]
    return point, covariance, delete_one


def closure(vector: Sequence[float]) -> complex:
    a_plus = complex(vector[0], vector[1])
    a_minus = complex(vector[2], vector[3])
    b_plus = complex(vector[4], vector[5])
    b_minus = complex(vector[6], vector[7])
    return a_plus * b_minus - a_minus * b_plus


def closure_score(point, delete_one):
    value = closure(point)
    replicates = [closure(row) for row in delete_one]
    center = sum(replicates) / len(replicates)
    factor = (len(replicates) - 1) / len(replicates)
    covariance = [
        [
            factor * sum(
                ((row.real, row.imag)[i] - (center.real, center.imag)[i])
                * ((row.real, row.imag)[j] - (center.real, center.imag)[j])
                for row in replicates
            )
            for j in range(2)
        ]
        for i in range(2)
    ]
    try:
        result = zero_score([value.real, value.imag], covariance)
    except ValueError as error:
        result = {"status": "singular_smoke_covariance", "error": str(error)}
    return {"point_re_im": [value.real, value.imag], "covariance": covariance, **result}


def denominator_rows(rows: Sequence[dict], separation: int):
    batch_values = [
        {
            f"{name}_{hand}": row[f"d{separation}_{name}_{hand}"] / row["samples"]
            for hand in HANDS for name in ("G1", "G2")
        }
        for row in rows
    ]
    output = {}
    for key in batch_values[0]:
        point = sum(row[key] for row in batch_values) / len(batch_values)
        variance = sum((row[key] - point) ** 2 for row in batch_values) / (
            len(batch_values) * (len(batch_values) - 1)
        )
        se = math.sqrt(variance)
        output[key] = {"point": point, "standard_error": se, "abs_z": abs(point) / se if se else math.inf}
    return output


def score(payload: dict, rows: Sequence[dict]) -> dict:
    if payload.get("schema") != SCHEMA or not payload["mapping_gate"]["passed"]:
        raise ValueError("wrong or failed multi-separation response")
    results = {}
    for separation in SEPARATIONS:
        denominators = denominator_rows(rows, separation)
        denominator_ready = all(row["abs_z"] >= 2.0 for row in denominators.values())
        sep_point, sep_covariance, sep_delete_one = jackknife(rows, separation, "separation")
        local_point, local_covariance, _ = jackknife(rows, separation, "local_variance")
        try:
            support = zero_score(local_point, local_covariance)
            support_status = "detected" if support["survival_p"] < 0.05 else "not_detected"
        except ValueError as error:
            support = {"status": "singular_smoke_covariance", "error": str(error)}
            support_status = "not_detected"
        diagnostic_closure = closure_score(sep_point, sep_delete_one)
        phase_status = (
            "interpretable" if denominator_ready and support_status == "detected"
            else "not_interpretable_until_nonzero_support"
        )
        results[str(separation)] = {
            "two_point_denominators": denominators,
            "two_point_ready_abs_z_ge_2": denominator_ready,
            "separation_normalized_order": [part for name in PRIMARY for part in (name + "_re", name + "_im")],
            "separation_normalized_point": sep_point,
            "separation_normalized_covariance": sep_covariance,
            "local_variance_normalized_point": local_point,
            "local_variance_normalized_covariance": local_covariance,
            "cubic_support_zero_score": {**support, "decision_at_0.05": support_status},
            "phase_closure": {**diagnostic_closure, "decision_status": phase_status},
        }
    return {
        "schema": "matching-one/z5-charged-multiseparation-score/v1",
        "status": "local_model_development_smoke",
        "score_order": [
            "two_point_denominator_nonzero",
            "local_variance_normalized_cubic_support_vs_zero",
            "conditional_separation_normalized_phase_closure",
        ],
        "separations": results,
        "connected_cumulant_identity": payload["mapping_gate"]["charged_connectedness"],
        "archive_reanalysis_boundary": (
            "be80f25 retains only batch-aggregated cubics; P226 retains a different global "
            "one-point marked row. Neither contains same-replica local pair products, "
            "separation labels, or cubic-pair cross-covariance."
        ),
        "claim_boundary": [
            "This is a <=5k plumbing and variance smoke, not production evidence.",
            "A phase closure is not scientifically interpreted unless nonzero cubic support and stable two-point denominators are both established.",
            "The transported deck basis fixes the remaining charged phase; the positive denominator cancels field magnitudes only.",
        ],
    }


def render(result: dict) -> str:
    lines = [
        "# P250 separation-normalized charged cubic smoke", "",
        "Score order: two-point support, cubic support, then conditional phase closure.", "",
        "| d | min two-point | cubic support chi2/8 | p | phase status |", "|---:|---:|---:|---:|---|",
    ]
    for separation, row in result["separations"].items():
        minimum = min(value["abs_z"] for value in row["two_point_denominators"].values())
        support = row["cubic_support_zero_score"]
        chi = support.get("chi_square", "singular")
        pvalue = support.get("survival_p", "--")
        lines.append(
            f"| {separation} | {minimum:.3g} sigma | {chi} | {pvalue} | {row['phase_closure']['decision_status']} |"
        )
    lines.extend(["", "## Boundary", "", *[f"- {line}" for line in result["claim_boundary"]], ""])
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
