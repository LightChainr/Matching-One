#!/usr/bin/env python3
"""Archive-first CRT de-gauging score for Issue 418."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Sequence

import numpy as np

import score_p406_spatial_fourier_cone as p406
from norm5_chiral_fixedp_mc import PHASES
from z5_projective_leg_bivariate_mc import contexts, rotate, rotation_gauges
from z5_projective_leg_cross_scale_mc import PARENT_GEOMETRY, PARENT_MATRIX


SCHEMA = "matching-one/p418-crt-degauging-score/v1"
GROUP_ORDER = 101
CRT_MULTIPLIER = 405
DECK_ORDER = 5
ROOT = Path(__file__).resolve().parents[1]
FREEZE = ROOT / "analysis/p418_crt_degauging_freeze.json"
P406_SCORE = ROOT / "results/huawei-20260830/P406-spatial-fourier-cone/score.json"
DEFAULT_OUTPUT = ROOT / "results/huawei-20260830/P418-crt-degauging/score.json"


def add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return left[0] + right[0], left[1] + right[1]


def scale(multiplier: int, point: tuple[int, int]) -> tuple[int, int]:
    return multiplier * point[0], multiplier * point[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def phase(exponent: int) -> complex:
    return complex(*PHASES[exponent % DECK_ORDER])


def residue(point: tuple[int, int]) -> int:
    return (point[0] - 10 * point[1]) % GROUP_ORDER


def exact_section_and_masks() -> dict[str, Any]:
    parent_coordinates = PARENT_GEOMETRY.coordinates
    residue_to_parent = {residue(point): index for index, point in enumerate(parent_coordinates)}
    if len(residue_to_parent) != GROUP_ORDER:
        raise AssertionError("parent residue coordinate is not bijective")
    gauges = rotation_gauges()
    hands: dict[str, Any] = {}
    deck_step = (PARENT_MATRIX[0][0], PARENT_MATRIX[1][0])
    for hand, context in zip(p406.HANDS, contexts()):
        section_points = [scale(CRT_MULTIPLIER, point) for point in parent_coordinates]
        section_vertices = [context.geometry.vertex(point) for point in section_points]
        projection_failures = sum(
            PARENT_GEOMETRY.vertex(point) != index
            for index, point in enumerate(section_points)
        )
        homomorphism_failures = 0
        for left_index, left_point in enumerate(parent_coordinates):
            for right_index, right_point in enumerate(parent_coordinates):
                parent_sum = PARENT_GEOMETRY.vertex(add(left_point, right_point))
                child_sum = context.geometry.vertex(add(section_points[left_index], section_points[right_index]))
                homomorphism_failures += int(child_sum != section_vertices[parent_sum])
        rotation_failures = 0
        for index, point in enumerate(parent_coordinates):
            parent_rotated = PARENT_GEOMETRY.vertex(rotate(point))
            child_rotated = context.geometry.vertex(rotate(section_points[index]))
            rotation_failures += int(child_rotated != section_vertices[parent_rotated])
        deck_annihilation_failures = int(
            context.geometry.vertex(scale(CRT_MULTIPLIER, deck_step))
            != context.geometry.vertex((0, 0))
        )
        offsets = []
        offset_failures = 0
        for parent_index, child_vertex in enumerate(section_vertices):
            matches = [
                fiber
                for fiber in range(DECK_ORDER)
                if context.field_to_vertex[DECK_ORDER * parent_index + fiber] == child_vertex
            ]
            if len(matches) != 1:
                offset_failures += 1
                offsets.append(-1)
            else:
                offsets.append(matches[0])
        if offset_failures:
            raise AssertionError("CRT fiber offset is not unique")
        residual_phase = [
            (gauges[hand][index] - offsets[index]) % DECK_ORDER
            for index in range(GROUP_ORDER)
        ]
        masks = {}
        minimum_abs = math.inf
        maximum_abs = 0.0
        zero_count = 0
        for charge in p406.CHARGES:
            rows = []
            complex_values = []
            for displacement_residue in range(GROUP_ORDER):
                displacement = parent_coordinates[residue_to_parent[displacement_residue]]
                counts = [0] * DECK_ORDER
                for parent_index, origin in enumerate(parent_coordinates):
                    target = PARENT_GEOMETRY.vertex(add(origin, displacement))
                    exponent = charge * (residual_phase[parent_index] - residual_phase[target])
                    counts[exponent % DECK_ORDER] += 1
                value = sum(count * phase(index) for index, count in enumerate(counts)) / GROUP_ORDER
                zero = len(set(counts)) == 1
                zero_count += int(zero)
                minimum_abs = min(minimum_abs, abs(value))
                maximum_abs = max(maximum_abs, abs(value))
                complex_values.append(value)
                rows.append(
                    {
                        "residue": displacement_residue,
                        "phase_counts": counts,
                        "real": value.real,
                        "imag": value.imag,
                        "abs": abs(value),
                        "exact_zero": zero,
                    }
                )
            spectrum = np.fft.fft(np.asarray(complex_values)) / GROUP_ORDER
            masks[f"r{charge}"] = {
                "values": rows,
                "A0_residual": abs(complex_values[0] - 1.0),
                "minimum_abs": min(abs(value) for value in complex_values),
                "minimum_abs_residue": int(np.argmin(np.abs(complex_values))),
                "maximum_abs": max(abs(value) for value in complex_values),
                "zero_count": sum(row["exact_zero"] for row in rows),
                "mask_fourier_min_real": float(np.min(spectrum.real)),
                "mask_fourier_max_abs_imag": float(np.max(np.abs(spectrum.imag))),
            }
        hands[hand] = {
            "section_vertices": section_vertices,
            "fiber_offset_b": offsets,
            "runner_rotation_gauge_t": gauges[hand],
            "residual_phase_u": residual_phase,
            "fiber_offset_counts": [offsets.count(value) for value in range(DECK_ORDER)],
            "residual_phase_counts": [residual_phase.count(value) for value in range(DECK_ORDER)],
            "gates": {
                "projection_failures": projection_failures,
                "homomorphism_pairs_checked": GROUP_ORDER * GROUP_ORDER,
                "homomorphism_failures": homomorphism_failures,
                "rotation_sites_checked": GROUP_ORDER,
                "rotation_failures": rotation_failures,
                "deck_annihilation_failures": deck_annihilation_failures,
                "fiber_offset_failures": offset_failures,
                "passed": not any(
                    (
                        projection_failures,
                        homomorphism_failures,
                        rotation_failures,
                        deck_annihilation_failures,
                        offset_failures,
                    )
                ),
            },
            "masks": masks,
        }
    witnesses = {
        "plus_r1": {"residue": 27, "phase_counts": [9, 30, 14, 20, 28]},
        "plus_r2": {"residue": 33, "phase_counts": [9, 30, 14, 20, 28]},
        "minus_r1": {"residue": 14, "phase_counts": [28, 8, 24, 15, 26]},
        "minus_r2": {"residue": 39, "phase_counts": [28, 8, 24, 15, 26]},
    }
    witness_failures = 0
    for channel, expected in witnesses.items():
        hand, charge_label = channel.split("_")
        observed = hands[hand]["masks"][charge_label]["values"][expected["residue"]]["phase_counts"]
        witness_failures += int(observed != expected["phase_counts"])
    if witness_failures:
        raise AssertionError("Issue 418 attenuation witness changed")
    return {
        "group_order": GROUP_ORDER,
        "deck_order": DECK_ORDER,
        "CRT_multiplier": CRT_MULTIPLIER,
        "arithmetic": {
            "405_mod_101": CRT_MULTIPLIER % GROUP_ORDER,
            "405_mod_5": CRT_MULTIPLIER % DECK_ORDER,
            "gcd_101_5": math.gcd(GROUP_ORDER, DECK_ORDER),
            "unique_homomorphic_section": True,
            "uniqueness_reason": "Hom(Z/101,Z/5)=0 because gcd(101,5)=1",
        },
        "hands": hands,
        "attenuation_witnesses": witnesses,
        "attenuation_witness_failures": witness_failures,
        "passed": all(row["gates"]["passed"] for row in hands.values()) and not witness_failures,
    }


def mask_complex(exact: dict[str, Any], hand: str, charge: int, coordinate: tuple[int, int]) -> complex:
    row = exact["hands"][hand]["masks"][f"r{charge}"]["values"][residue(coordinate)]
    return complex(row["real"], row["imag"])


def masked_design(exact: dict[str, Any], hand: str, charge: int, coordinates: list[tuple[int, int]]) -> np.ndarray:
    raw = p406.design(coordinates)
    output = raw.copy()
    for index, coordinate in enumerate(coordinates):
        value = mask_complex(exact, hand, charge, coordinate)
        real = raw[2 * index].copy()
        imag = raw[2 * index + 1].copy()
        output[2 * index] = value.real * real - value.imag * imag
        output[2 * index + 1] = value.imag * real + value.real * imag
    return output


def quantiles(values: np.ndarray) -> dict[str, float]:
    return {
        "q01": float(np.quantile(values, 0.01)),
        "q50": float(np.quantile(values, 0.50)),
        "q99": float(np.quantile(values, 0.99)),
    }


def prediction_envelope(
    exact: dict[str, Any], hand: str, charge: int, bootstrap_weights: np.ndarray, observed_residues: set[int]
) -> dict[str, Any]:
    coordinates_by_residue = [
        next(point for point in PARENT_GEOMETRY.coordinates if residue(point) == value)
        for value in range(GROUP_ORDER)
    ]
    predictions = masked_design(exact, hand, charge, coordinates_by_residue) @ bootstrap_weights.T
    rows = []
    observed_widths = []
    heldout_widths = []
    for value in range(GROUP_ORDER):
        real = quantiles(predictions[2 * value])
        imag = quantiles(predictions[2 * value + 1])
        width = math.hypot(real["q99"] - real["q01"], imag["q99"] - imag["q01"])
        is_observed = value in observed_residues
        (observed_widths if is_observed else heldout_widths).append(width)
        rows.append(
            {
                "residue": value,
                "observed_in_archive": is_observed,
                "real": real,
                "imag": imag,
                "complex_envelope_width": width,
            }
        )
    return {
        "status": "model_conditional_on_rejected_cone",
        "accepted_model_envelope": False,
        "residues": rows,
        "summary": {
            "observed_residues": len(observed_widths),
            "heldout_residues": len(heldout_widths),
            "median_observed_width": float(np.median(observed_widths)),
            "maximum_observed_width": max(observed_widths),
            "median_heldout_width": float(np.median(heldout_widths)),
            "maximum_heldout_width": max(heldout_widths),
        },
    }


def score_channel(
    exact: dict[str, Any],
    blocks: list[dict[str, object]],
    hand: str,
    charge: int,
    bootstrap: int,
    seed: int,
) -> dict[str, Any]:
    raw_whitened = []
    masked_whitened = []
    observed_residues: set[int] = set()
    for block in blocks:
        coordinates = block["coordinates"]
        observed_residues.update(residue(coordinate) for coordinate in coordinates)
        raw_whitened.append(p406.whiten(block["values"], p406.design(coordinates)))
        masked_whitened.append(
            p406.whiten(block["values"], masked_design(exact, hand, charge, coordinates))
        )
    raw_matrix = np.vstack([row["X"] for row in raw_whitened])
    masked_matrix = np.vstack([row["X"] for row in masked_whitened])
    vector = np.concatenate([row["y"] for row in raw_whitened])
    raw_weights, raw_statistic = p406.fit_nonnegative(raw_matrix, vector)
    masked_weights, masked_statistic = p406.fit_nonnegative(masked_matrix, vector)
    masked_unconstrained, _, masked_rank, _ = np.linalg.lstsq(masked_matrix, vector, rcond=1e-10)
    masked_floor = float(np.sum((masked_matrix @ masked_unconstrained - vector) ** 2))
    raw_fitted = raw_matrix @ raw_weights
    masked_fitted = masked_matrix @ masked_weights
    rng = np.random.default_rng(seed)
    noises = rng.standard_normal((bootstrap, len(vector)))
    raw_reference = []
    masked_reference = []
    masked_bootstrap_weights = []
    for noise in noises:
        raw_star, raw_value = p406.fit_nonnegative(raw_matrix, raw_fitted + noise)
        masked_star, masked_value = p406.fit_nonnegative(masked_matrix, masked_fitted + noise)
        del raw_star
        raw_reference.append(raw_value)
        masked_reference.append(masked_value)
        masked_bootstrap_weights.append(masked_star)
    raw_p = (1 + sum(value >= raw_statistic for value in raw_reference)) / (bootstrap + 1)
    masked_p = (1 + sum(value >= masked_statistic for value in masked_reference)) / (bootstrap + 1)
    bootstrap_weights = np.asarray(masked_bootstrap_weights)
    weight_q01 = np.quantile(bootstrap_weights, 0.01, axis=0)
    weight_q99 = np.quantile(bootstrap_weights, 0.99, axis=0)
    weight_width = weight_q99 - weight_q01
    block_contributions = []
    raw_residual = raw_fitted - vector
    masked_residual = masked_fitted - vector
    start = 0
    for name, raw_row, masked_row in zip(("radius4", "radius5", "radius6"), raw_whitened, masked_whitened):
        length = len(raw_row["y"])
        block_contributions.append(
            {
                "block": name,
                "resolved_modes": length,
                "raw_distance_squared": float(np.sum(raw_residual[start : start + length] ** 2)),
                "masked_distance_squared": float(np.sum(masked_residual[start : start + length] ** 2)),
            }
        )
        start += length
    prediction = prediction_envelope(exact, hand, charge, bootstrap_weights, observed_residues)
    masked_reference_array = np.asarray(masked_reference)
    decision = "masked_cone_not_rejected" if masked_p >= 0.01 else "masked_cone_rejected"
    return {
        "raw_cone": {
            "distance_squared": raw_statistic,
            "bootstrap_p": raw_p,
            "bootstrap_quantiles": quantiles(np.asarray(raw_reference)),
        },
        "masked_cone": {
            "distance_squared": masked_statistic,
            "bootstrap_p": masked_p,
            "bootstrap_quantiles": quantiles(masked_reference_array),
            "unconstrained_distance_squared": masked_floor,
            "cone_increment_over_unconstrained": masked_statistic - masked_floor,
            "decision": decision,
        },
        "increment_over_raw_cone": masked_statistic - raw_statistic,
        "distance_ratio_to_raw": masked_statistic / raw_statistic,
        "masked_unconstrained_floor_minus_raw_cone": masked_floor - raw_statistic,
        "design": {
            "resolved_whitened_coordinates": len(vector),
            "rank": int(masked_rank),
            "spectral_coordinates": GROUP_ORDER,
            "linear_nullity": GROUP_ORDER - int(masked_rank),
            "structurally_unique_spectrum": False,
            "interpretation": "rank/nullity describe the finite spatial design, not physical state count",
        },
        "spectral_mass_bootstrap": {
            "status": "model_conditional_on_rejected_cone" if decision == "masked_cone_rejected" else "model_conditional",
            "total_mass": quantiles(bootstrap_weights.sum(axis=1)),
            "per_frequency_interval_width": {
                "minimum": float(np.min(weight_width)),
                "median": float(np.median(weight_width)),
                "maximum": float(np.max(weight_width)),
                "frequencies_with_q01_positive": int(np.sum(weight_q01 > 0)),
            },
        },
        "prediction_envelope": prediction,
        "block_contributions": block_contributions,
        "observed_mask_abs": {
            "minimum": min(
                abs(mask_complex(exact, hand, charge, coordinate))
                for block in blocks
                for coordinate in block["coordinates"]
            ),
            "maximum": max(
                abs(mask_complex(exact, hand, charge, coordinate))
                for block in blocks
                for coordinate in block["coordinates"]
            ),
        },
        "boundary": (
            "Bootstrap spectral masses and predictions are conditional diagnostics. A rejected mask-times-cone "
            "has no accepted CRT-degauged prediction envelope."
        ),
    }


def build_result(paths: list[tuple[Path, str]], bootstrap: int, seed: int) -> dict[str, Any]:
    freeze = json.loads(FREEZE.read_text(encoding="utf-8"))
    if freeze["status"] != "frozen_before_archive_reveal":
        raise AssertionError("P418 protocol is not frozen")
    exact = exact_section_and_masks()
    archives = [p406.read_block(path, digest) for path, digest in paths]
    published = json.loads(P406_SCORE.read_text(encoding="utf-8"))
    channels = {}
    raw_replay_max_distance_error = 0.0
    raw_replay_max_p_error = 0.0
    for hand_index, hand in enumerate(p406.HANDS):
        for charge in p406.CHARGES:
            key = f"{hand}_r{charge}"
            blocks = [archive[(hand, charge)] for archive in archives]
            row = score_channel(exact, blocks, hand, charge, bootstrap, seed + 1000 * hand_index + charge)
            raw_replay_max_distance_error = max(
                raw_replay_max_distance_error,
                abs(row["raw_cone"]["distance_squared"] - published["channels"][key]["minimum_cone_distance_squared"]),
            )
            raw_replay_max_p_error = max(
                raw_replay_max_p_error,
                abs(row["raw_cone"]["bootstrap_p"] - published["channels"][key]["bootstrap_p"]),
            )
            channels[key] = row
    if raw_replay_max_distance_error > 1e-9 or raw_replay_max_p_error > 1e-15:
        raise AssertionError("P406 raw cone replay changed")
    all_rejected = all(
        row["masked_cone"]["decision"] == "masked_cone_rejected"
        for row in channels.values()
    )
    return {
        "schema": SCHEMA,
        "status": "archived_production_data_scored_once",
        "issues": [418, 406, 250],
        "freeze": {
            "path": str(FREEZE.relative_to(ROOT)),
            "sha256": sha256(FREEZE),
            "commit": "5977ce3",
        },
        "new_monte_carlo": False,
        "inputs": [
            {"path": str(path.resolve().relative_to(ROOT)), "sha256": digest}
            for path, digest in paths
        ],
        "exact_CRT_section_and_masks": exact,
        "raw_P406_replay": {
            "score_sha256": sha256(P406_SCORE),
            "maximum_distance_error": raw_replay_max_distance_error,
            "maximum_p_error": raw_replay_max_p_error,
            "passed": True,
        },
        "bootstrap": {
            "replicates": bootstrap,
            "seed": seed,
            "same_channel_noise_draws_for_raw_and_masked": True,
        },
        "channels": channels,
        "decision": "mask_times_positive_fourier_cone_rejected_in_all_channels" if all_rejected else "mixed_channel_decision",
        "mechanism_update": (
            "The raw positive Fourier cone survives, but the stricter exact CRT-section mask-times-positive cone "
            "is rejected in every archived hand/charge channel. The deterministic phase mask therefore does not "
            "explain the raw endpoint complexity under the stated section/stationarity/observable-transport contract."
        ),
        "first_followup": (
            "Audit the ensemble factorization against the runner's one-anchor sampling and the translation transport "
            "of the projective-leg root observable before assigning the residual to a physical state mechanism."
        ),
        "claim_boundary": freeze["claim_boundary"],
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--radius4", type=Path, required=True)
    parser.add_argument("--radius4-sha256", required=True)
    parser.add_argument("--radius5", type=Path, required=True)
    parser.add_argument("--radius5-sha256", required=True)
    parser.add_argument("--radius6", type=Path, required=True)
    parser.add_argument("--radius6-sha256", required=True)
    parser.add_argument("--bootstrap", type=int, default=250)
    parser.add_argument("--seed", type=int, default=40610120260830)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    result = build_result(
        [
            (args.radius4, args.radius4_sha256),
            (args.radius5, args.radius5_sha256),
            (args.radius6, args.radius6_sha256),
        ],
        args.bootstrap,
        args.seed,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "decision": result["decision"],
                "channels": {
                    key: {
                        "raw": row["raw_cone"]["distance_squared"],
                        "masked": row["masked_cone"]["distance_squared"],
                        "p": row["masked_cone"]["bootstrap_p"],
                    }
                    for key, row in result["channels"].items()
                },
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
