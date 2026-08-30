#!/usr/bin/env python3
"""Frozen bivariate common-state rank and algebra score for P250/P249/P255."""

from __future__ import annotations

import argparse
import cmath
import csv
import hashlib
from itertools import permutations
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from score_z5_charged_threepoint import zero_score
from score_z5_projective_leg_cross_scale import jackknife_covariance
from score_z5_projective_leg_pair_transfer import CHANNELS
from score_z5_projective_leg_state_dimension import solve_complex
from z5_projective_leg_bivariate_mc import GRID, SCHEMA, label, rotate


AXIS_TRAIN = ((0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (0, 3))
MIXED_COMMUTING_GATE = ((1, 1), (2, 1), (1, 2))
DEGREE4_HOLDOUT = ((4, 0), (3, 1), (2, 2), (1, 3), (0, 4))
ROTATION_REPRESENTATIVES = ((1, 0), (2, 0), (3, 0), (4, 0), (1, 1), (2, 1), (3, 1), (2, 2))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    fields = [key for key in kept[0] if key.startswith(("ap", "am"))]
    return {key: sum(row[key] for row in kept) / samples for key in fields}


def pair(values: Mapping[str, float], point: tuple[int, int], channel: tuple[str, int]) -> complex:
    hand, charge = channel
    prefix = f"{label(*point)}_r{charge}_{hand}_"
    return complex(values[prefix + "re"], values[prefix + "im"])


def fit_axis_recurrence(values: Mapping[str, float], axis: str, rank: int) -> list[complex]:
    design = []
    observed = []
    for channel in CHANNELS:
        for end in range(rank, 4):
            point = (end, 0) if axis == "x" else (0, end)
            observed.append(pair(values, point, channel))
            design.append([
                pair(values, (end - offset, 0) if axis == "x" else (0, end - offset), channel)
                for offset in range(1, rank + 1)
            ])
    normal = [[
        sum(row[i].conjugate() * row[j] for row in design)
        for j in range(rank)
    ] for i in range(rank)]
    rhs = [sum(row[i].conjugate() * value for row, value in zip(design, observed)) for i in range(rank)]
    return solve_complex(normal, rhs)


def polynomial(coefficients: Sequence[complex], value: complex) -> complex:
    rank = len(coefficients)
    return value**rank - sum(coefficients[index] * value ** (rank - 1 - index) for index in range(rank))


def recurrence_roots(coefficients: Sequence[complex]) -> list[complex]:
    rank = len(coefficients)
    if rank == 1:
        return [coefficients[0]]
    radius = 1.0 + max(abs(value) for value in coefficients)
    roots = [radius * cmath.exp(2j * math.pi * (index + 0.25) / rank) for index in range(rank)]
    for _ in range(500):
        updated = []
        largest = 0.0
        for index, root in enumerate(roots):
            denominator = math.prod(root - other for other_index, other in enumerate(roots) if other_index != index)
            if abs(denominator) < 1e-20:
                denominator += 1e-20
            value = root - polynomial(coefficients, root) / denominator
            updated.append(value)
            largest = max(largest, abs(value - root))
        roots = updated
        if largest < 1e-14:
            break
    if max(abs(polynomial(coefficients, root)) for root in roots) > 1e-8:
        raise ValueError("recurrence roots did not converge")
    return sorted(roots, key=lambda value: (-abs(value), math.atan2(value.imag, value.real)))


def fit_amplitudes(
    values: Mapping[str, float], xroots: Sequence[complex], yroots: Sequence[complex], points: Sequence[tuple[int, int]]
) -> dict[tuple[str, int], list[complex]]:
    design = [[xroots[index] ** a * yroots[index] ** b for index in range(len(xroots))] for a, b in points]
    normal = [[
        sum(row[i].conjugate() * row[j] for row in design)
        for j in range(len(xroots))
    ] for i in range(len(xroots))]
    output = {}
    for channel in CHANNELS:
        observed = [pair(values, point, channel) for point in points]
        rhs = [sum(row[i].conjugate() * value for row, value in zip(design, observed)) for i in range(len(xroots))]
        output[channel] = solve_complex(normal, rhs)
    return output


def predict(
    amplitudes: Mapping[tuple[str, int], Sequence[complex]], xroots: Sequence[complex], yroots: Sequence[complex],
    point: tuple[int, int], channel: tuple[str, int],
) -> complex:
    a, b = point
    return sum(
        amplitudes[channel][index] * xroots[index] ** a * yroots[index] ** b
        for index in range(len(xroots))
    )


def fit_model(values: Mapping[str, float], rank: int) -> dict:
    xcoefficients = fit_axis_recurrence(values, "x", rank)
    ycoefficients = fit_axis_recurrence(values, "y", rank)
    xroots = recurrence_roots(xcoefficients)
    yroots_unpaired = recurrence_roots(ycoefficients)
    candidates = []
    for permutation in permutations(range(rank)):
        yroots = [yroots_unpaired[index] for index in permutation]
        amplitudes = fit_amplitudes(values, xroots, yroots, AXIS_TRAIN)
        objective = sum(
            abs(pair(values, point, channel) - predict(amplitudes, xroots, yroots, point, channel)) ** 2
            for channel in CHANNELS for point in AXIS_TRAIN
        )
        candidates.append((objective, permutation, yroots, amplitudes))
    objective, permutation, yroots, amplitudes = min(candidates, key=lambda row: (row[0], row[1]))
    return {
        "rank": rank,
        "xcoefficients": xcoefficients,
        "ycoefficients": ycoefficients,
        "xroots": xroots,
        "yroots": yroots,
        "yroot_permutation": list(permutation),
        "amplitudes": amplitudes,
        "axis_objective": objective,
    }


def model_residual(values: Mapping[str, float], model: Mapping[str, object], points: Sequence[tuple[int, int]]) -> list[float]:
    output = []
    for channel in CHANNELS:
        for point in points:
            residual = pair(values, point, channel) - predict(
                model["amplitudes"], model["xroots"], model["yroots"], point, channel
            )
            output.extend((residual.real, residual.imag))
    return output


def characteristic_residual(model: Mapping[str, object]) -> list[float]:
    output = []
    for xvalue, yvalue in zip(model["xcoefficients"], model["ycoefficients"]):
        residual = xvalue - yvalue
        output.extend((residual.real, residual.imag))
    return output


def vector_score(point: Sequence[float], deleted: Sequence[Sequence[float]]) -> dict:
    covariance = jackknife_covariance(deleted)
    return {"residual": list(point), "covariance": covariance, "zero_score": zero_score(point, covariance)}


def complex_payload(values: Sequence[complex]) -> list[dict]:
    return [{
        "re": value.real,
        "im": value.imag,
        "abs": abs(value),
        "phase": math.atan2(value.imag, value.real),
    } for value in values]


def rank_score(values: Mapping[str, float], deleted_values: Sequence[Mapping[str, float]], rank: int) -> dict:
    model = fit_model(values, rank)
    deleted_models = [fit_model(row, rank) for row in deleted_values]
    mixed = vector_score(
        model_residual(values, model, MIXED_COMMUTING_GATE),
        [model_residual(row, fitted, MIXED_COMMUTING_GATE) for row, fitted in zip(deleted_values, deleted_models)],
    )
    holdout = vector_score(
        model_residual(values, model, DEGREE4_HOLDOUT),
        [model_residual(row, fitted, DEGREE4_HOLDOUT) for row, fitted in zip(deleted_values, deleted_models)],
    )
    characteristic = vector_score(
        characteristic_residual(model),
        [characteristic_residual(fitted) for fitted in deleted_models],
    )
    return {
        "rank": rank,
        "fit_points": [list(point) for point in AXIS_TRAIN],
        "x_characteristic_coefficients": complex_payload(model["xcoefficients"]),
        "y_characteristic_coefficients": complex_payload(model["ycoefficients"]),
        "joint_eigenpairs": [
            {"x": complex_payload([model["xroots"][index]])[0], "y": complex_payload([model["yroots"][index]])[0]}
            for index in range(rank)
        ],
        "yroot_permutation": model["yroot_permutation"],
        "axis_fit_objective": model["axis_objective"],
        "commuting_mixed_gate": {"points": [list(point) for point in MIXED_COMMUTING_GATE], **mixed},
        "degree4_heldout": {"points": [list(point) for point in DEGREE4_HOLDOUT], **holdout},
        "C4_similarity_necessary_characteristic_score": characteristic,
        "algebra_scope": "The joint diagonal realization enforces [Tx,Ty]=0; the mixed gate tests that commuting product out of axis fit.",
    }


def rotated_expected(values: Mapping[str, float], point: tuple[int, int], channel: tuple[str, int]) -> complex:
    hand, charge = channel
    if hand == "plus" and charge == 1:
        return pair(values, point, ("plus", 2)).conjugate()
    if hand == "plus" and charge == 2:
        return pair(values, point, ("plus", 1))
    if hand == "minus" and charge == 1:
        return pair(values, point, ("minus", 2))
    if hand == "minus" and charge == 2:
        return pair(values, point, ("minus", 1)).conjugate()
    raise ValueError("unknown hand-charge channel")


def rotation_residual(values: Mapping[str, float]) -> list[float]:
    output = []
    for point in ROTATION_REPRESENTATIVES:
        rotated = rotate(point)
        for channel in CHANNELS:
            residual = pair(values, rotated, channel) - rotated_expected(values, point, channel)
            output.extend((residual.real, residual.imag))
    return output


def validate(payload: Mapping[str, object], manifest: Mapping[str, object]) -> None:
    if payload.get("schema") != SCHEMA or not payload.get("exact_gate", {}).get("passed"):
        raise ValueError("wrong or failed bivariate response")
    if payload.get("manifest_runner_commit") != manifest.get("runner_commit"):
        raise ValueError("runner commit differs from manifest")
    for key, value in manifest["run"].items():
        if payload["run"].get(key) != value:
            raise ValueError(f"run differs from manifest for {key}")


def score(payload: Mapping[str, object], batches: Sequence[dict], manifest: Mapping[str, object]) -> dict:
    path = Path(manifest["input_batches"])
    if sha256(path) != manifest["input_batches_sha256"]:
        raise ValueError("bivariate batch hash changed")
    validate(payload, manifest)
    values = means(batches)
    deleted_values = [means(batches, index) for index in range(len(batches))]
    ranks = {str(rank): rank_score(values, deleted_values, rank) for rank in (1, 2, 3)}
    rotation = vector_score(rotation_residual(values), [rotation_residual(row) for row in deleted_values])
    alpha = float(manifest["decision_alpha"])
    minimal = None
    for rank in (1, 2, 3):
        row = ranks[str(rank)]
        if (
            row["commuting_mixed_gate"]["zero_score"]["survival_p"] >= alpha
            and row["degree4_heldout"]["zero_score"]["survival_p"] >= alpha
        ):
            minimal = rank
            break
    decision = f"minimal_commuting_common_rank_{minimal}" if minimal else "no_commuting_common_rank_le_3"
    return {
        "schema": "matching-one/z5-projective-leg-bivariate-state-score/v1",
        "status": "fresh_frozen_reveal",
        "rank_order": [1, 2, 3],
        "rank_models": ranks,
        "minimal_rank": minimal,
        "decision": decision,
        "C4_covariant_rotation_score": {
            "representatives": [list(point) for point in ROTATION_REPRESENTATIVES],
            "channel_map": payload["exact_gate"]["rotation_fiber_gate"]["channel_maps"],
            "channel_map_fourth_power": "identity",
            **rotation,
        },
        "algebra_decision": {
            "commutator": "tested as held-out mixed-moment closure of an axis-fitted commuting realization",
            "Ty_equals_RTxRinv": "characteristic-polynomial equality is a necessary score; neutral moments do not uniquely identify R up to the state centralizer",
            "R4": "the exact realified channel map has fourth power identity and is scored on signed C4 displacement orbits",
            "D5": "not identifiable: the neutral charged pair cancels simultaneous one-step deck phase",
        },
        "decision_alpha": alpha,
        "claim_boundary": [
            "The first passing rank is the minimal common commuting realization on the frozen stencil, not a count of continuum primary fields.",
            "Failure through rank three leaves higher rank, non-diagonal/Jordan state, and finite-torus image structure open.",
            "A passing characteristic score is necessary but not sufficient to construct a unique latent R.",
            "No scalar phase label is selected by this score.",
        ],
    }


def render(result: Mapping[str, object]) -> str:
    lines = [
        "# P250 bivariate common-state score",
        "",
        "| rank | mixed commuting gate | p | degree-4 heldout | p | C4 characteristic p |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for rank in (1, 2, 3):
        row = result["rank_models"][str(rank)]
        mixed = row["commuting_mixed_gate"]["zero_score"]
        holdout = row["degree4_heldout"]["zero_score"]
        c4 = row["C4_similarity_necessary_characteristic_score"]["zero_score"]
        lines.append(
            f"| {rank} | {mixed['chi_square']:.6g}/{mixed['degrees_of_freedom']} | {mixed['survival_p']:.6g} "
            f"| {holdout['chi_square']:.6g}/{holdout['degrees_of_freedom']} | {holdout['survival_p']:.6g} | {c4['survival_p']:.6g} |"
        )
    rotation = result["C4_covariant_rotation_score"]["zero_score"]
    lines += [
        "",
        f"Signed C4 covariant orbit score: `{rotation['chi_square']:.6g}/{rotation['degrees_of_freedom']}`, p `{rotation['survival_p']:.6g}`.",
        f"Decision: `{result['decision']}`.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("response", type=Path)
    parser.add_argument("batches", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    result = score(json.loads(args.response.read_text()), read_batches(args.batches), manifest)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
