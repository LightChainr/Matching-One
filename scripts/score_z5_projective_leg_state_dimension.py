#!/usr/bin/env python3
"""Minimal complex state-dimension score for the P250 projective-leg row."""

from __future__ import annotations

import argparse
import cmath
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from score_z5_charged_threepoint import zero_score
from score_z5_projective_leg_cross_scale import add_covariance, jackknife_covariance
from score_z5_projective_leg_pair_transfer import CHANNELS, means, read_batches, transfer


TARGET_FIT_LAST = 4
TARGET_HOLDOUT = 5
SOURCE_LAST = 3
SOURCE_L = math.sqrt(65.0)
TARGET_L = math.sqrt(101.0)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def series(values: Mapping[str, float], last: int) -> dict[tuple[str, int], list[complex]]:
    return {
        channel: [transfer(values, distance, *channel) for distance in range(1, last + 1)]
        for channel in CHANNELS
    }


def solve_complex(matrix: Sequence[Sequence[complex]], rhs: Sequence[complex]) -> list[complex]:
    size = len(rhs)
    rows = [[complex(value) for value in matrix[i]] + [complex(rhs[i])] for i in range(size)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(rows[row][column]))
        if abs(rows[pivot][column]) < 1e-18:
            raise ValueError("singular complex normal matrix")
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = rows[column][column]
        rows[column] = [value / scale for value in rows[column]]
        for row in range(size):
            if row == column:
                continue
            factor = rows[row][column]
            rows[row] = [value - factor * pivot_value for value, pivot_value in zip(rows[row], rows[column])]
    return [rows[index][-1] for index in range(size)]


def fit_recurrence(rows: Mapping[tuple[str, int], Sequence[complex]], rank: int, fit_last: int) -> list[complex]:
    design = []
    observed = []
    for channel in CHANNELS:
        values = rows[channel]
        for end in range(rank, fit_last):
            observed.append(values[end])
            design.append([values[end - offset] for offset in range(1, rank + 1)])
    normal = [[
        sum(row[i].conjugate() * row[j] for row in design)
        for j in range(rank)
    ] for i in range(rank)]
    rhs = [sum(row[i].conjugate() * value for row, value in zip(design, observed)) for i in range(rank)]
    return solve_complex(normal, rhs)


def recurrence_residual(
    rows: Mapping[tuple[str, int], Sequence[complex]], coefficients: Sequence[complex], end: int
) -> list[float]:
    rank = len(coefficients)
    output = []
    for channel in CHANNELS:
        values = rows[channel]
        value = values[end - 1] - sum(
            coefficients[offset - 1] * values[end - 1 - offset]
            for offset in range(1, rank + 1)
        )
        output.extend((value.real, value.imag))
    return output


def training_residuals(rows, coefficients, fit_last: int) -> list[float]:
    output = []
    rank = len(coefficients)
    for end in range(rank + 1, fit_last + 1):
        output.extend(recurrence_residual(rows, coefficients, end))
    return output


def recurrence_roots(coefficients: Sequence[complex]) -> list[complex]:
    if len(coefficients) == 1:
        return [coefficients[0]]
    if len(coefficients) == 2:
        first, second = coefficients
        discriminant = cmath.sqrt(first * first + 4.0 * second)
        return [(first + discriminant) / 2.0, (first - discriminant) / 2.0]
    return []


def coefficients_from_roots(roots: Sequence[complex]) -> list[complex]:
    if len(roots) == 1:
        return [roots[0]]
    if len(roots) == 2:
        return [roots[0] + roots[1], -roots[0] * roots[1]]
    raise ValueError("root transport is implemented only for ranks one and two")


def transport_coefficients(coefficients: Sequence[complex], scale: float) -> list[complex]:
    roots = recurrence_roots(coefficients)
    transported = [cmath.exp(scale * cmath.log(root)) for root in roots]
    return coefficients_from_roots(transported)


def fit_amplitudes(rows, roots: Sequence[complex], fit_last: int) -> dict[str, list[list[float]]]:
    output = {}
    for channel in CHANNELS:
        values = rows[channel][:fit_last]
        design = [[root**distance for root in roots] for distance in range(fit_last)]
        normal = [[sum(row[i].conjugate() * row[j] for row in design) for j in range(len(roots))] for i in range(len(roots))]
        rhs = [sum(row[i].conjugate() * value for row, value in zip(design, values)) for i in range(len(roots))]
        amplitudes = solve_complex(normal, rhs)
        output[f"{channel[0]}_r{channel[1]}"] = [[value.real, value.imag] for value in amplitudes]
    return output


def complex_payload(values: Sequence[complex]) -> list[dict]:
    return [{
        "re": value.real,
        "im": value.imag,
        "abs": abs(value),
        "phase_radians": math.atan2(value.imag, value.real),
    } for value in values]


def heldout_score(target_batches: Sequence[dict], rank: int) -> tuple[dict, list[list[complex]]]:
    full_rows = series(means(target_batches), TARGET_HOLDOUT)
    coefficients = fit_recurrence(full_rows, rank, TARGET_FIT_LAST)
    residual = recurrence_residual(full_rows, coefficients, TARGET_HOLDOUT)
    deleted_residuals = []
    deleted_coefficients = []
    for omitted in range(len(target_batches)):
        deleted_rows = series(means(target_batches, omitted), TARGET_HOLDOUT)
        deleted = fit_recurrence(deleted_rows, rank, TARGET_FIT_LAST)
        deleted_coefficients.append(deleted)
        deleted_residuals.append(recurrence_residual(deleted_rows, deleted, TARGET_HOLDOUT))
    covariance = jackknife_covariance(deleted_residuals)
    return ({
        "rank": rank,
        "fit_distances": list(range(1, TARGET_FIT_LAST + 1)),
        "heldout_distance": TARGET_HOLDOUT,
        "coefficients": complex_payload(coefficients),
        "eigenvalues": complex_payload(recurrence_roots(coefficients)),
        "channel_amplitudes": fit_amplitudes(full_rows, recurrence_roots(coefficients), TARGET_FIT_LAST),
        "training_residual": training_residuals(full_rows, coefficients, TARGET_FIT_LAST),
        "heldout_residual_order": [f"{hand}_r{charge}_{part}" for hand, charge in CHANNELS for part in ("re", "im")],
        "heldout_residual": residual,
        "heldout_covariance": covariance,
        "heldout_zero_score": zero_score(residual, covariance),
    }, deleted_coefficients)


def source_constraint(
    source_batches: Sequence[dict], target_coefficients: Sequence[complex],
    deleted_target_coefficients: Sequence[Sequence[complex]], scale: float,
) -> dict:
    source_rows = series(means(source_batches), SOURCE_LAST)
    coefficients = transport_coefficients(target_coefficients, scale)
    rank = len(coefficients)
    residual = []
    for end in range(rank + 1, SOURCE_LAST + 1):
        residual.extend(recurrence_residual(source_rows, coefficients, end))
    source_deleted_residuals = []
    for omitted in range(len(source_batches)):
        rows = series(means(source_batches, omitted), SOURCE_LAST)
        row = []
        for end in range(rank + 1, SOURCE_LAST + 1):
            row.extend(recurrence_residual(rows, coefficients, end))
        source_deleted_residuals.append(row)
    source_covariance = jackknife_covariance(source_deleted_residuals)
    target_deleted_residuals = []
    for target in deleted_target_coefficients:
        transported = transport_coefficients(target, scale)
        row = []
        for end in range(rank + 1, SOURCE_LAST + 1):
            row.extend(recurrence_residual(source_rows, transported, end))
        target_deleted_residuals.append(row)
    covariance = add_covariance(source_covariance, jackknife_covariance(target_deleted_residuals))
    return {
        "transport_scale_on_complex_log_eigenvalue": scale,
        "transported_coefficients": complex_payload(coefficients),
        "transported_eigenvalues": complex_payload(recurrence_roots(coefficients)),
        "residual": residual,
        "covariance": covariance,
        "zero_score": zero_score(residual, covariance),
    }


def block_hankel_singular_values(rows) -> list[float]:
    vectors = []
    for channel in CHANNELS:
        values = rows[channel]
        vectors.extend([[values[d], values[d + 1]] for d in range(3)])
    a = sum(abs(row[0]) ** 2 for row in vectors)
    d = sum(abs(row[1]) ** 2 for row in vectors)
    b = sum(row[0].conjugate() * row[1] for row in vectors)
    delta = math.sqrt((a - d) ** 2 + 4.0 * abs(b) ** 2)
    eigenvalues = [(a + d + delta) / 2.0, max((a + d - delta) / 2.0, 0.0)]
    return [math.sqrt(value) for value in eigenvalues]


def image_kernel(distance: int, parent: tuple[int, int], alpha: float) -> float:
    a, b = parent
    period1 = (a, b)
    period2 = (-b, a)
    total = 0.0
    for axis in ((distance, 0), (0, distance)):
        for first in (-1, 0, 1):
            for second in (-1, 0, 1):
                x = axis[0] + first * period1[0] + second * period2[0]
                y = axis[1] + first * period1[1] + second * period2[1]
                total += (x * x + y * y) ** (-alpha / 2.0)
    return total / 2.0


def image_fit(rows, parent: tuple[int, int], fit_last: int) -> dict:
    def at(alpha: float):
        shape = [image_kernel(distance, parent, alpha) for distance in range(1, fit_last + 1)]
        amplitudes = {}
        objective = 0.0
        for channel in CHANNELS:
            observed = rows[channel][:fit_last]
            denominator = sum(value * value for value in shape)
            amplitude = sum(value * target for value, target in zip(shape, observed)) / denominator
            amplitudes[channel] = amplitude
            objective += sum(abs(target - amplitude * value) ** 2 for value, target in zip(shape, observed))
        return objective, amplitudes
    left, right = 0.25, 4.0
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    x1, x2 = right - ratio * (right - left), left + ratio * (right - left)
    f1, f2 = at(x1)[0], at(x2)[0]
    for _ in range(80):
        if f1 > f2:
            left, x1, f1 = x1, x2, f2
            x2 = left + ratio * (right - left)
            f2 = at(x2)[0]
        else:
            right, x2, f2 = x2, x1, f1
            x1 = right - ratio * (right - left)
            f1 = at(x1)[0]
    alpha = (left + right) / 2.0
    objective, amplitudes = at(alpha)
    return {"alpha": alpha, "objective": objective, "amplitudes": amplitudes}


def image_residual(rows, fit, parent: tuple[int, int], heldout: int) -> list[float]:
    shape = image_kernel(heldout, parent, fit["alpha"])
    output = []
    for channel in CHANNELS:
        value = rows[channel][heldout - 1] - fit["amplitudes"][channel] * shape
        output.extend((value.real, value.imag))
    return output


def image_scores(target_batches: Sequence[dict], source_batches: Sequence[dict]) -> dict:
    target_rows = series(means(target_batches), TARGET_HOLDOUT)
    target_fit = image_fit(target_rows, (10, 1), TARGET_FIT_LAST)
    target_residual = image_residual(target_rows, target_fit, (10, 1), TARGET_HOLDOUT)
    target_deleted_fits = []
    target_deleted_residuals = []
    for omitted in range(len(target_batches)):
        rows = series(means(target_batches, omitted), TARGET_HOLDOUT)
        fit = image_fit(rows, (10, 1), TARGET_FIT_LAST)
        target_deleted_fits.append(fit)
        target_deleted_residuals.append(image_residual(rows, fit, (10, 1), TARGET_HOLDOUT))
    target_covariance = jackknife_covariance(target_deleted_residuals)

    source_rows = series(means(source_batches), SOURCE_LAST)
    def source_fit_for_alpha(rows, alpha):
        shape = [image_kernel(distance, (8, 1), alpha) for distance in (1, 2)]
        amplitudes = {}
        for channel in CHANNELS:
            denominator = sum(value * value for value in shape)
            amplitudes[channel] = sum(value * target for value, target in zip(shape, rows[channel][:2])) / denominator
        return {"alpha": alpha, "amplitudes": amplitudes}
    source_fit = source_fit_for_alpha(source_rows, target_fit["alpha"])
    source_residual = image_residual(source_rows, source_fit, (8, 1), 3)
    source_deleted_residuals = []
    for omitted in range(len(source_batches)):
        rows = series(means(source_batches, omitted), SOURCE_LAST)
        fit = source_fit_for_alpha(rows, target_fit["alpha"])
        source_deleted_residuals.append(image_residual(rows, fit, (8, 1), 3))
    source_covariance = jackknife_covariance(source_deleted_residuals)
    target_alpha_residuals = []
    for fit in target_deleted_fits:
        transported = source_fit_for_alpha(source_rows, fit["alpha"])
        target_alpha_residuals.append(image_residual(source_rows, transported, (8, 1), 3))
    source_total_covariance = add_covariance(source_covariance, jackknife_covariance(target_alpha_residuals))
    return {
        "definition": "axis-average nearest 3x3 lattice-image sum of |r+P n|^(-alpha)",
        "target_fit_distances": [1, 2, 3, 4],
        "target_alpha": target_fit["alpha"],
        "target_heldout_d5": {
            "residual": target_residual,
            "covariance": target_covariance,
            "zero_score": zero_score(target_residual, target_covariance),
        },
        "source_amplitudes_fit_distances": [1, 2],
        "source_heldout_d3": {
            "residual": source_residual,
            "covariance": source_total_covariance,
            "zero_score": zero_score(source_residual, source_total_covariance),
        },
        "claim_boundary": "This finite nearest-image kernel is a simple geometry alternative, not the universal torus primary two-point function.",
    }


def score(target_batches, source_batches, manifest) -> dict:
    for label, path in (("target", Path(manifest["inputs"]["target_batches"])), ("source", Path(manifest["inputs"]["source_batches"]))):
        if sha256(path) != manifest["inputs"][f"{label}_sha256"]:
            raise ValueError(f"{label} batch hash changed")
    ranks = {}
    deleted = {}
    for rank in (1, 2, 3):
        ranks[str(rank)], deleted[str(rank)] = heldout_score(target_batches, rank)
        if rank <= 2:
            coefficients = [complex(row["re"], row["im"]) for row in ranks[str(rank)]["coefficients"]]
            ranks[str(rank)]["source_geometry_constraints"] = {
                "lattice_fixed_eigenvalues": source_constraint(source_batches, coefficients, deleted[str(rank)], 1.0),
                "conformal_1_over_L_eigenvalues": source_constraint(
                    source_batches, coefficients, deleted[str(rank)], TARGET_L / SOURCE_L
                ),
            }
        else:
            ranks[str(rank)]["source_geometry_constraints"] = {
                "status": "not_identifiable_from_source_d1_d3_for_rank3"
            }
    alpha = float(manifest["decision_alpha"])
    minimal = None
    for rank in (1, 2):
        row = ranks[str(rank)]
        target_pass = row["heldout_zero_score"]["survival_p"] >= alpha
        source_pass = any(
            constraint["zero_score"]["survival_p"] >= alpha
            for constraint in row["source_geometry_constraints"].values()
        )
        if target_pass and source_pass:
            minimal = rank
            break
    if minimal is None and ranks["3"]["heldout_zero_score"]["survival_p"] >= alpha:
        decision = "rank3_closes_target_but_cross_geometry_dimension_is_not_identifiable"
    elif minimal is None:
        decision = "rank3_also_fails_target_heldout"
    else:
        decision = f"minimal_cross_geometry_recurrence_rank_{minimal}"
    target_rows = series(means(target_batches), TARGET_HOLDOUT)
    singular = block_hankel_singular_values(target_rows)
    return {
        "schema": "matching-one/z5-projective-leg-state-dimension-score/v1",
        "status": "retrospective_no_new_simulation",
        "model_order": ["rank1", "shared_rank2", "rank3_if_needed", "nearest_image_kernel"],
        "target_block_hankel_singular_values": singular,
        "target_block_hankel_s2_over_s1": singular[1] / singular[0],
        "ranks": ranks,
        "nearest_image_kernel": image_scores(target_batches, source_batches),
        "decision_alpha": alpha,
        "decision": decision,
        "claim_boundary": [
            "All rows are existing N325/N505 pair batches; no simulation was run.",
            "Ranks one and two share complex eigenvalues across all four channels and retain channel-specific amplitudes.",
            "The N505 d5 row and the declared N325 recurrence are held out from target fitting.",
            "Rank three lacks enough N325 distances for a cross-geometry constraint.",
        ],
    }


def render(result) -> str:
    lines = ["# P250 minimal projective-leg state dimension", "", "| rank | N505 d5 chi2/df | p | N325 fixed p | N325 1/L p |", "|---:|---:|---:|---:|---:|"]
    for rank in (1, 2, 3):
        row = result["ranks"][str(rank)]
        target = row["heldout_zero_score"]
        if rank <= 2:
            fixed = row["source_geometry_constraints"]["lattice_fixed_eigenvalues"]["zero_score"]["survival_p"]
            scaled = row["source_geometry_constraints"]["conformal_1_over_L_eigenvalues"]["zero_score"]["survival_p"]
            fixed_text, scaled_text = f"{fixed:.6g}", f"{scaled:.6g}"
        else:
            fixed_text = scaled_text = "not identifiable"
        lines.append(f"| {rank} | {target['chi_square']:.6g}/{target['degrees_of_freedom']} | {target['survival_p']:.6g} | {fixed_text} | {scaled_text} |")
    image = result["nearest_image_kernel"]
    lines += [
        "", f"Decision: `{result['decision']}`.",
        f"Nearest-image alpha: `{image['target_alpha']}`; N505 d5 p `{image['target_heldout_d5']['zero_score']['survival_p']}`; N325 d3 p `{image['source_heldout_d3']['zero_score']['survival_p']}`.", "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    result = score(
        read_batches(Path(manifest["inputs"]["target_batches"])),
        read_batches(Path(manifest["inputs"]["source_batches"])), manifest,
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
