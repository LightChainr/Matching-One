#!/usr/bin/env python3
"""Direction versus internal-charge tomography for the P250 rank-two row."""

from __future__ import annotations

import argparse
import cmath
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from score_z5_charged_threepoint import chi_square_survival, zero_score
from score_z5_projective_leg_cross_scale import jackknife_covariance
from score_z5_projective_leg_pair_transfer import CHANNELS, means, read_batches, transfer
from score_z5_projective_leg_state_dimension import (
    fit_amplitudes,
    fit_recurrence,
    recurrence_roots,
    series,
    solve_complex,
)


FIT_LAST = 4
LAST = 5
ZETA5 = cmath.exp(2j * math.pi / 5.0)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def axis_average_difference(values: Mapping[str, float], distance: int, hand: str, charge: int):
    tkey = f"d{distance}_T{charge}_{hand}_"
    akey = f"d{distance}_A{charge}_{hand}_"
    average = complex(values[tkey + "re"], values[tkey + "im"])
    difference = complex(values[akey + "re"], values[akey + "im"])
    return average, difference, average + difference, average - difference


def sorted_rank2_fit(values: Mapping[str, float]):
    rows = series(values, LAST)
    coefficients = fit_recurrence(rows, 2, FIT_LAST)
    roots = sorted(recurrence_roots(coefficients), key=abs, reverse=True)
    payload = fit_amplitudes(rows, roots, FIT_LAST)
    amplitudes = {
        channel: [complex(*value) for value in payload[f"{channel[0]}_r{channel[1]}"]]
        for channel in CHANNELS
    }
    return roots, amplitudes


def direction_factor(q: complex) -> complex:
    """Map the frozen axis relation Ay=q Ax into A/T mode amplitudes.

    With T=(X+Y)/2 and A=(X-Y)/2, a mode satisfying Ay=q Ax has
    A_mode/T_mode=(1-q)/(1+q).  The q=-1 sector is invisible in T and hence
    cannot be diagnosed from roots fitted to the pre-existing T row.
    """
    if abs(1.0 + q) < 1e-12:
        raise ValueError("q=-1 cancels from the axis-average row and is not identifiable here")
    return (1.0 - q) / (1.0 + q)


def candidate_residual(values: Mapping[str, float], q: complex) -> list[float]:
    roots, amplitudes = sorted_rank2_fit(values)
    factor = direction_factor(q)
    output = []
    for hand, charge in CHANNELS:
        second_amplitude = amplitudes[(hand, charge)][1]
        for distance in range(1, LAST + 1):
            _, observed, _, _ = axis_average_difference(values, distance, hand, charge)
            predicted = factor * second_amplitude * roots[1] ** (distance - 1)
            residual = observed - predicted
            output.extend((residual.real, residual.imag))
    return output


def fit_difference_amplitudes(values: Mapping[str, float], roots: Sequence[complex]):
    """Fit both rank-two amplitudes of the archived axis-difference row."""
    design = [[root**distance for root in roots] for distance in range(FIT_LAST)]
    normal = [[
        sum(row[i].conjugate() * row[j] for row in design)
        for j in range(len(roots))
    ] for i in range(len(roots))]
    output = {}
    for channel in CHANNELS:
        observed = [
            axis_average_difference(values, distance, *channel)[1]
            for distance in range(1, FIT_LAST + 1)
        ]
        rhs = [
            sum(row[i].conjugate() * value for row, value in zip(design, observed))
            for i in range(len(roots))
        ]
        output[channel] = solve_complex(normal, rhs)
    return output


def mode_candidate_residual(values: Mapping[str, float], q: complex) -> list[float]:
    """Test only the second mode, leaving the leading A-mode unrestricted."""
    roots, average_amplitudes = sorted_rank2_fit(values)
    difference_amplitudes = fit_difference_amplitudes(values, roots)
    factor = direction_factor(q)
    output = []
    for channel in CHANNELS:
        residual = difference_amplitudes[channel][1] - factor * average_amplitudes[channel][1]
        output.extend((residual.real, residual.imag))
    return output


def heldout_mode_residual(values: Mapping[str, float], q: complex) -> list[float]:
    """Fit a channel-specific leading A amplitude on d1--4 and predict d5."""
    roots, average_amplitudes = sorted_rank2_fit(values)
    factor = direction_factor(q)
    design = [roots[0] ** distance for distance in range(FIT_LAST)]
    denominator = sum(abs(value) ** 2 for value in design)
    output = []
    for channel in CHANNELS:
        second = factor * average_amplitudes[channel][1]
        observed = [
            axis_average_difference(values, distance, *channel)[1]
            for distance in range(1, FIT_LAST + 1)
        ]
        leading = sum(
            basis.conjugate() * (value - second * roots[1] ** distance)
            for distance, (basis, value) in enumerate(zip(design, observed))
        ) / denominator
        prediction = leading * roots[0] ** (LAST - 1) + second * roots[1] ** (LAST - 1)
        residual = axis_average_difference(values, LAST, *channel)[1] - prediction
        output.extend((residual.real, residual.imag))
    return output


def unrestricted_mode_payload(values: Mapping[str, float]) -> dict:
    roots, average_amplitudes = sorted_rank2_fit(values)
    difference_amplitudes = fit_difference_amplitudes(values, roots)
    output = {}
    for channel in CHANNELS:
        average = average_amplitudes[channel][1]
        difference = difference_amplitudes[channel][1]
        q = (average - difference) / (average + difference)
        output[f"{channel[0]}_r{channel[1]}"] = {
            "q_re": q.real,
            "q_im": q.imag,
            "q_abs": abs(q),
            "q_phase": math.atan2(q.imag, q.real),
        }
    return output


def unrestricted_shared_roots_heldout_residual(values: Mapping[str, float]) -> list[float]:
    """Predict A(d5) after fitting both A amplitudes at the frozen T roots."""
    roots, _ = sorted_rank2_fit(values)
    amplitudes = fit_difference_amplitudes(values, roots)
    output = []
    for channel in CHANNELS:
        prediction = sum(amplitudes[channel][index] * roots[index] ** (LAST - 1) for index in range(2))
        residual = axis_average_difference(values, LAST, *channel)[1] - prediction
        output.extend((residual.real, residual.imag))
    return output


def conjugate_residual(values: Mapping[str, float]) -> list[float]:
    output = []
    for hand, charge in CHANNELS:
        for distance in range(1, LAST + 1):
            _, _, xvalue, yvalue = axis_average_difference(values, distance, hand, charge)
            residual = yvalue - xvalue.conjugate()
            output.extend((residual.real, residual.imag))
    return output


def residual_score(point, deleted):
    covariance = jackknife_covariance(deleted)
    heldout_indices = [10 * channel + coordinate for channel in range(len(CHANNELS)) for coordinate in (8, 9)]
    heldout = [point[index] for index in heldout_indices]
    heldout_covariance = [[covariance[i][j] for j in heldout_indices] for i in heldout_indices]
    return {
        "residual_order": [
            f"{hand}_r{charge}_d{distance}_{part}"
            for hand, charge in CHANNELS for distance in range(1, LAST + 1) for part in ("re", "im")
        ],
        "residual": point,
        "covariance": covariance,
        "all_distances_zero_score": zero_score(point, covariance),
        "heldout_d5": {
            "residual": heldout,
            "covariance": heldout_covariance,
            "zero_score": zero_score(heldout, heldout_covariance),
        },
    }


def vector_score(point, deleted):
    covariance = jackknife_covariance(deleted)
    return {
        "residual_order": [f"{hand}_r{charge}_{part}" for hand, charge in CHANNELS for part in ("re", "im")],
        "residual": point,
        "covariance": covariance,
        "zero_score": zero_score(point, covariance),
    }


def wrap_phase(value: float) -> float:
    return math.atan2(math.sin(value), math.cos(value))


def align_roots(reference: Sequence[complex], roots: Sequence[complex]) -> list[complex]:
    direct = abs(reference[0] - roots[0]) + abs(reference[1] - roots[1])
    swapped = abs(reference[0] - roots[1]) + abs(reference[1] - roots[0])
    return list(roots if direct <= swapped else roots[::-1])


def phase_alphabet(values: Mapping[str, float], deleted_values: Sequence[Mapping[str, float]]) -> dict:
    roots, _ = sorted_rank2_fit(values)
    phase = math.atan2(roots[1].imag, roots[1].real)
    deleted_phases = []
    for row in deleted_values:
        candidate, _ = sorted_rank2_fit(row)
        aligned = align_roots(roots, candidate)
        deleted_phases.append(phase + wrap_phase(math.atan2(aligned[1].imag, aligned[1].real) - phase))
    center = sum(deleted_phases) / len(deleted_phases)
    variance = (len(deleted_phases) - 1) / len(deleted_phases) * sum((value - center) ** 2 for value in deleted_phases)
    candidates = {
        "C4_minus_i": -math.pi / 2.0,
        "C4_plus_i": math.pi / 2.0,
        **{f"Z5_j{index}": wrap_phase(2.0 * math.pi * index / 5.0) for index in range(5)},
    }
    scores = {}
    for name, target in candidates.items():
        residual = wrap_phase(phase - target)
        chi_square = residual * residual / variance
        scores[name] = {
            "target_phase": target,
            "wrapped_residual": residual,
            "chi_square": chi_square,
            "degrees_of_freedom": 1,
            "survival_p": chi_square_survival(chi_square, 1),
        }
    return {
        "second_root": {
            "re": roots[1].real, "im": roots[1].imag,
            "abs": abs(roots[1]), "phase": phase,
            "phase_standard_error": math.sqrt(variance),
        },
        "exact_alphabet_scores": scores,
        "deck_generator_statement": "one Z5 deck step has phase 2*pi*j/5; -pi/2 is not in that exact alphabet",
        "scope": "translation-root phase alone may be a compound/non-generator step, so the direction score is the discriminator",
    }


def score(batches: Sequence[dict], manifest: Mapping[str, object]) -> dict:
    path = Path(manifest["input"]["batches"])
    if sha256(path) != manifest["input"]["sha256"]:
        raise ValueError("N505 batch hash changed")
    full_values = means(batches)
    deleted_values = [means(batches, omitted) for omitted in range(len(batches))]
    candidates = {
        "internal_direction_invariant": 1.0 + 0j,
        "spatial_C4_minus_i": -1j,
        "spatial_C4_plus_i": 1j,
        **{f"direction_Z5_j{index}": ZETA5**index for index in range(1, 5)},
    }
    strict_rows = {}
    mode_rows = {}
    for name, q in candidates.items():
        row = residual_score(
            candidate_residual(full_values, q),
            [candidate_residual(values, q) for values in deleted_values],
        )
        row.update({"q_re": q.real, "q_im": q.imag, "q_phase": math.atan2(q.imag, q.real)})
        strict_rows[name] = row
        mode = vector_score(
            mode_candidate_residual(full_values, q),
            [mode_candidate_residual(values, q) for values in deleted_values],
        )
        heldout = vector_score(
            heldout_mode_residual(full_values, q),
            [heldout_mode_residual(values, q) for values in deleted_values],
        )
        mode.update({
            "q_re": q.real,
            "q_im": q.imag,
            "q_phase": math.atan2(q.imag, q.real),
            "heldout_d5": heldout,
        })
        mode_rows[name] = mode
    conjugate = residual_score(
        conjugate_residual(full_values),
        [conjugate_residual(values) for values in deleted_values],
    )
    unrestricted_shared_roots = vector_score(
        unrestricted_shared_roots_heldout_residual(full_values),
        [unrestricted_shared_roots_heldout_residual(values) for values in deleted_values],
    )
    strict_ranking = sorted(strict_rows, key=lambda name: strict_rows[name]["all_distances_zero_score"]["chi_square"])
    ranking = sorted(mode_rows, key=lambda name: mode_rows[name]["zero_score"]["chi_square"])
    alpha = float(manifest["decision_alpha"])
    passing = [
        name for name in ranking
        if mode_rows[name]["zero_score"]["survival_p"] >= alpha
        and mode_rows[name]["heldout_d5"]["zero_score"]["survival_p"] >= alpha
    ]
    if unrestricted_shared_roots["zero_score"]["survival_p"] < alpha:
        decision = "axis_difference_requires_additional_or_different_transfer_state_no_q_identification"
    elif passing == ["internal_direction_invariant"]:
        decision = "direction_invariant_second_state_not_spatial_C4_character"
    elif not passing:
        decision = "no_frozen_direction_character_closes"
    else:
        decision = "direction_character_not_unique"
    return {
        "schema": "matching-one/z5-projective-leg-direction-charge-score/v1",
        "status": "exact_reconstruction_no_new_simulation",
        "reconstruction": {
            "x": "T+A", "y": "T-A",
            "coordinates": "all complex d1-d5 rows, four hand-charge channels",
            "covariance": "delete-one covariance propagated from the original 160 batches",
        },
        "rank2_source_commit": manifest["rank2_source_commit"],
        "phase_alphabet": phase_alphabet(full_values, deleted_values),
        "strict_zero_leading_A_candidates": strict_rows,
        "strict_ranking": strict_ranking,
        "mode_resolved_direction_candidates": mode_rows,
        "unrestricted_second_mode_direction_characters": unrestricted_mode_payload(full_values),
        "unrestricted_A_at_frozen_T_roots_heldout_d5": unrestricted_shared_roots,
        "conjugate_direction_candidate": conjugate,
        "ranking": ranking,
        "decision_alpha": alpha,
        "decision": decision,
        "claim_boundary": [
            "The test distinguishes displacement-direction character from internal charge; it does not rename the second state.",
            "A direction-invariant result does not make the root phase a one-step Z5 generator phase.",
            "If an internal-only account is retained despite a non-Z5 root phase, its transfer step cannot be identified with one deck generator.",
            "The pair observable is C4-even in ensemble, so a spatial character must appear through the frozen mode-resolved x/y relation to survive this score.",
            "No new random stream or rank selection is used.",
            "The mode-resolved score leaves the leading axis-difference amplitude free in each channel and tests only the second root.",
            "Failure of the unrestricted A-row d5 check means the T-row rank-two roots cannot be used as a complete direction-tomography basis.",
        ],
    }


def render(result) -> str:
    lines = ["# P250 direction/internal-charge tomography", "", "| candidate | second-mode chi2/df | p | heldout d5 p |", "|---|---:|---:|---:|"]
    for name in result["ranking"]:
        row = result["mode_resolved_direction_candidates"][name]
        full = row["zero_score"]
        held = row["heldout_d5"]["zero_score"]
        lines.append(f"| {name} | {full['chi_square']:.6g}/{full['degrees_of_freedom']} | {full['survival_p']:.6g} | {held['survival_p']:.6g} |")
    conjugate = result["conjugate_direction_candidate"]
    lines += [
        "", f"Unrestricted A amplitudes at the frozen T roots, d5: p `{result['unrestricted_A_at_frozen_T_roots_heldout_d5']['zero_score']['survival_p']}`.",
        "", f"Conjugate x/y candidate: p `{conjugate['all_distances_zero_score']['survival_p']}`; d5 p `{conjugate['heldout_d5']['zero_score']['survival_p']}`.",
        f"Decision: `{result['decision']}`.", "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    result = score(read_batches(Path(manifest["input"]["batches"])), manifest)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
