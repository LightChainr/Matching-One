#!/usr/bin/env python3
"""Amplitude-free exponential/power transfer score for the P250 pair row."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from score_z5_charged_threepoint import zero_score
from z5_projective_leg_pair_transfer_mc import CHARGES, HANDS, SCHEMA, SEPARATIONS


CHANNELS = tuple((hand, charge) for hand in HANDS for charge in CHARGES)
PRIMARY_SEPARATIONS = (1, 2, 3)


def read_batches(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return [
            {
                key: int(value) if key in {"batch", "replica_first", "samples"}
                else value if key in {"field_sha256", "translation_sha256"}
                else float(value)
                for key, value in row.items()
            }
            for row in csv.DictReader(handle)
        ]


def means(rows: Sequence[dict], excluded: int | None = None) -> dict[str, float]:
    kept = [row for index, row in enumerate(rows) if index != excluded]
    samples = sum(row["samples"] for row in kept)
    fields = [key for key in kept[0] if key.startswith("d")]
    return {key: sum(row[key] for row in kept) / samples for key in fields}


def transfer(values: Mapping[str, float], separation: int, hand: str, charge: int) -> complex:
    prefix = f"d{separation}_T{charge}_{hand}_"
    return complex(values[prefix + "re"], values[prefix + "im"])


def wrap_phase(value: float) -> float:
    return math.atan2(math.sin(value), math.cos(value))


def channel_metrics(values: Mapping[str, float], hand: str, charge: int) -> dict[str, float]:
    first = transfer(values, 1, hand, charge)
    second = transfer(values, 2, hand, charge)
    third = transfer(values, 3, hand, charge)
    if min(abs(first), abs(second), abs(third)) <= 0.0:
        raise ValueError("zero transfer denominator")
    ratio12 = second / first
    ratio23 = third / second
    mass12 = -math.log(abs(ratio12))
    mass23 = -math.log(abs(ratio23))
    eta12 = mass12 / math.log(2.0)
    eta23 = mass23 / math.log(1.5)
    phase12 = math.atan2(ratio12.imag, ratio12.real)
    phase23 = math.atan2(ratio23.imag, ratio23.real)
    return {
        "ratio12_re": ratio12.real, "ratio12_im": ratio12.imag,
        "ratio23_re": ratio23.real, "ratio23_im": ratio23.imag,
        "mass12": mass12, "mass23": mass23,
        "eta12": eta12, "eta23": eta23,
        "phase12": phase12, "phase23": phase23,
        "exponential_mass_delta": mass23 - mass12,
        "power_eta_delta": eta23 - eta12,
        "phase_step_delta": wrap_phase(phase23 - phase12),
    }


def metric_vectors(values: Mapping[str, float]) -> tuple[list[float], list[float], list[float], dict]:
    exponential = []
    power = []
    phase = []
    details = {}
    for hand, charge in CHANNELS:
        row = channel_metrics(values, hand, charge)
        name = f"{hand}_r{charge}"
        details[name] = row
        exponential.extend((row["exponential_mass_delta"], row["phase_step_delta"]))
        power.extend((row["power_eta_delta"], row["phase_step_delta"]))
        phase.extend((row["phase12"], row["phase23"]))
    return exponential, power, phase, details


def jackknife_covariance(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    count = len(rows)
    center = [sum(row[index] for row in rows) / count for index in range(len(rows[0]))]
    factor = (count - 1) / count
    return [
        [
            factor * sum(
                (row[first] - center[first]) * (row[second] - center[second])
                for row in rows
            )
            for second in range(len(center))
        ]
        for first in range(len(center))
    ]


def raw_transfer_rows(batches: Sequence[dict]) -> dict:
    output = {}
    batch_count = len(batches)
    for separation in SEPARATIONS:
        for hand, charge in CHANNELS:
            key = f"d{separation}_T{charge}_{hand}"
            rows = [
                (batch[key + "_re"] / batch["samples"],
                 batch[key + "_im"] / batch["samples"])
                for batch in batches
            ]
            point = [sum(row[index] for row in rows) / batch_count for index in range(2)]
            covariance = [
                [
                    sum((row[i] - point[i]) * (row[j] - point[j]) for row in rows)
                    / (batch_count * (batch_count - 1))
                    for j in range(2)
                ]
                for i in range(2)
            ]
            real_se = math.sqrt(covariance[0][0])
            output[key] = {
                "point_re_im": point,
                "covariance": covariance,
                "real_abs_z": abs(point[0]) / real_se if real_se else math.inf,
                "complex_zero_score": zero_score(point, covariance),
            }
    return output


def score(payload: dict, batches: Sequence[dict], manifest: Mapping[str, object]) -> dict:
    if payload.get("schema") != SCHEMA or not payload["exact_gate"]["passed"]:
        raise ValueError("wrong or failed pair-transfer response")
    for key, expected in manifest["run"].items():
        if payload["run"].get(key) != expected:
            raise ValueError(f"run differs from manifest for {key}")
    raw = raw_transfer_rows(batches)
    minimum_z = min(
        raw[f"d{separation}_T{charge}_{hand}"]["real_abs_z"]
        for separation in PRIMARY_SEPARATIONS for hand, charge in CHANNELS
    )
    transfer_gate = minimum_z >= float(manifest["transfer_gate"]["minimum_real_abs_z"])
    if not transfer_gate:
        models = {
            "status": "locked_primary_transfer_gate_failed",
            "computed": False,
        }
    else:
        full_values = means(batches)
        exp_point, power_point, phase_point, details = metric_vectors(full_values)
        exp_deleted, power_deleted, phase_deleted = [], [], []
        for omitted in range(len(batches)):
            exp_row, power_row, phase_row, _ = metric_vectors(means(batches, omitted))
            exp_deleted.append(exp_row)
            power_deleted.append(power_row)
            phase_deleted.append(phase_row)
        exp_covariance = jackknife_covariance(exp_deleted)
        power_covariance = jackknife_covariance(power_deleted)
        phase_covariance = jackknife_covariance(phase_deleted)
        phase_delta_indices = [2 * index + 1 for index in range(len(CHANNELS))]
        phase_delta = [exp_point[index] for index in phase_delta_indices]
        phase_delta_covariance = [
            [exp_covariance[i][j] for j in phase_delta_indices]
            for i in phase_delta_indices
        ]
        models = {
            "status": "amplitude_free_transfer_models_revealed",
            "computed": True,
            "channel_metrics": details,
            "exponential": {
                "null": "constant complex adjacent eigenvalue: mass12=mass23 and phase12=phase23",
                "residual_order": [
                    coordinate for hand, charge in CHANNELS
                    for coordinate in (f"{hand}_r{charge}_mass_delta", f"{hand}_r{charge}_phase_delta")
                ],
                "residual": exp_point,
                "covariance": exp_covariance,
                "zero_score": zero_score(exp_point, exp_covariance),
            },
            "power": {
                "null": "constant effective eta with a constant deck phase step",
                "residual_order": [
                    coordinate for hand, charge in CHANNELS
                    for coordinate in (f"{hand}_r{charge}_eta_delta", f"{hand}_r{charge}_phase_delta")
                ],
                "residual": power_point,
                "covariance": power_covariance,
                "zero_score": zero_score(power_point, power_covariance),
            },
            "deck_phase": {
                "order": [
                    coordinate for hand, charge in CHANNELS
                    for coordinate in (f"{hand}_r{charge}_phase12", f"{hand}_r{charge}_phase23")
                ],
                "point": phase_point,
                "covariance": phase_covariance,
                "nonzero_score": zero_score(phase_point, phase_covariance),
                "constant_step_score": zero_score(phase_delta, phase_delta_covariance),
            },
        }
    return {
        "schema": "matching-one/z5-projective-leg-pair-transfer-score/v1",
        "status": "fresh_pair_only_transfer_reveal",
        "score_order": [
            "d1_d3_real_transfer_resolution",
            "amplitude_free_exponential_residual",
            "amplitude_free_power_residual",
            "deck_character_phase",
        ],
        "minimum_primary_real_abs_z": minimum_z,
        "primary_transfer_gate_passed": transfer_gate,
        "raw_transfer": raw,
        "models": models,
        "cubic_fields_used": [],
        "claim_boundary": [
            "Only pair-transfer rows enter this score.",
            "The exponential and power tests are amplitude-free local transfer-shape diagnostics over d1,d2,d3.",
            "Failure of both simple shapes permits mixtures or finite-torus transfer spectra; it is not absence of propagation.",
        ],
    }


def render(result: Mapping[str, object]) -> str:
    lines = [
        "# P250 projective-leg pair-only transfer spectrum", "",
        f"Primary d1-d3 minimum real z: `{result['minimum_primary_real_abs_z']}`; gate `{result['primary_transfer_gate_passed']}`.", "",
        "| d | plus r1 | plus r2 | minus r1 | minus r2 |", "|---:|---:|---:|---:|---:|",
    ]
    for separation in SEPARATIONS:
        cells = [
            result["raw_transfer"][f"d{separation}_T{charge}_{hand}"]["real_abs_z"]
            for hand, charge in CHANNELS
        ]
        lines.append(f"| {separation} | " + " | ".join(f"{value:.3g}" for value in cells) + " |")
    if result["models"]["computed"]:
        exponential = result["models"]["exponential"]["zero_score"]
        power = result["models"]["power"]["zero_score"]
        phase = result["models"]["deck_phase"]
        lines.extend([
            "", "## Amplitude-free shape scores", "",
            f"- exponential: `{exponential['chi_square']}/{exponential['degrees_of_freedom']}`, p `{exponential['survival_p']}`",
            f"- power: `{power['chi_square']}/{power['degrees_of_freedom']}`, p `{power['survival_p']}`",
            f"- deck phase nonzero: `{phase['nonzero_score']['chi_square']}/{phase['nonzero_score']['degrees_of_freedom']}`, p `{phase['nonzero_score']['survival_p']}`",
            f"- deck phase constant-step: `{phase['constant_step_score']['chi_square']}/{phase['constant_step_score']['degrees_of_freedom']}`, p `{phase['constant_step_score']['survival_p']}`",
        ])
    else:
        lines.extend(["", f"Models: `{result['models']['status']}`."])
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("response", type=Path)
    parser.add_argument("batches", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = score(
        json.loads(args.response.read_text()), read_batches(args.batches),
        json.loads(args.manifest.read_text()),
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
