#!/usr/bin/env python3
"""Extract and compare the two hand-specific projective-leg annihilator lines."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Callable, Mapping, Sequence

import numpy as np

from norm5_chiral_fixedp_mc import M_MINUS, M_PLUS
from score_z5_projective_leg_bivariate_state import means, pair, read_batches
from score_z5_projective_leg_hankel_rank import MONOMIALS_2, covariance_score
from z5_projective_leg_bivariate_mc import rotate
from z5_projective_leg_cross_scale_mc import CHILD_MINUS, CHILD_PLUS


ALPHA = 0.01
HANDS = ("plus", "minus")
CHARGES = (1, 2)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reflect(point: tuple[int, int]) -> tuple[int, int]:
    return point[0], -point[1]


def rotate_power(point: tuple[int, int], power: int) -> tuple[int, int]:
    for _ in range(power % 4):
        point = rotate(point)
    return point


def transformed_basis(*, alexander_reflection: bool, rotation_power: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        rotate_power(reflect(point) if alexander_reflection else point, rotation_power)
        for point in MONOMIALS_2
    )


def hand_hankel(values: Mapping[str, float], hand: str, basis: Sequence[tuple[int, int]]) -> np.ndarray:
    return np.asarray([
        [pair(values, (left[0] + right[0], left[1] + right[1]), (hand, charge)) for right in basis]
        for charge in CHARGES for left in basis
    ], dtype=complex)


def annihilator_line(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    _, singular_values, right_adjoint = np.linalg.svd(matrix, full_matrices=False)
    line = right_adjoint[-1].conjugate()
    line /= np.linalg.norm(line)
    return line, singular_values


def canonical(line: np.ndarray, pivot: int) -> np.ndarray:
    if abs(line[pivot]) < 1e-14:
        raise ValueError("projective pivot vanished")
    return line / line[pivot]


def comparison_pivot(first: np.ndarray, second: np.ndarray) -> int:
    shared_resolution = np.minimum(np.abs(first), np.abs(second))
    return int(np.argmax(shared_resolution))


def projective_residual(first: np.ndarray, second: np.ndarray, pivot: int) -> list[float]:
    difference = canonical(first, pivot) - canonical(second, pivot)
    output = []
    for index, value in enumerate(difference):
        if index != pivot:
            output.extend((float(value.real), float(value.imag)))
    return output


def line_payload(line: np.ndarray, basis: Sequence[tuple[int, int]]) -> dict:
    pivot = int(np.argmax(np.abs(line)))
    normalized = canonical(line, pivot)
    return {
        "basis": [list(point) for point in basis],
        "normalization_pivot": pivot,
        "coefficients": [
            {"re": float(value.real), "im": float(value.imag), "abs": float(abs(value))}
            for value in normalized
        ],
    }


def bridge_score(
    plus: np.ndarray, minus: np.ndarray, deleted_plus: Sequence[np.ndarray], deleted_minus: Sequence[np.ndarray],
    *, conjugate_plus: bool,
) -> dict:
    mapped = plus.conjugate() if conjugate_plus else plus
    deleted_mapped = [row.conjugate() if conjugate_plus else row for row in deleted_plus]
    pivot = comparison_pivot(mapped, minus)
    point = projective_residual(mapped, minus, pivot)
    deleted = [
        projective_residual(first, second, pivot)
        for first, second in zip(deleted_mapped, deleted_minus)
    ]
    return {
        "coefficient_conjugation": conjugate_plus,
        "comparison_pivot": pivot,
        "score": covariance_score(point, deleted),
    }


def exact_map_gate() -> dict:
    plus = complex(CHILD_PLUS[0][0], CHILD_PLUS[1][0])
    minus = complex(CHILD_MINUS[0][0], CHILD_MINUS[1][0])
    d4_orbit_conjugate_plus = {
        complex(round(value.real), round(value.imag))
        for value in (plus.conjugate(), 1j * plus.conjugate(), -plus.conjugate(), -1j * plus.conjugate())
    }
    return {
        "norm5_multipliers": {"plus": M_PLUS, "minus": M_MINUS},
        "C4_fiber_multipliers_mod5": {"plus": 3, "minus": 2},
        "fiber_multipliers_are_inverse": (3 * 2) % 5 == 1,
        "child_Gaussian_periods": {
            "plus": [int(plus.real), int(plus.imag)],
            "minus": [int(minus.real), int(minus.imag)],
        },
        "same_parent_children_are_exact_D4_reflections": minus in d4_orbit_conjugate_plus,
        "geometry_scope": "The same-parent children are a Hecke hand pair, not exactly D4-isomorphic finite quotients.",
        "deck_row_action": "Any invertible charge permutation/phase acts on Hankel rows and leaves the common right annihilator line invariant.",
        "deck_maps_tested_without_new_parameter": [1, 2, 3, 4],
    }


def score(batches: Sequence[dict], manifest: Mapping[str, object]) -> dict:
    path = Path(manifest["input_batches"])
    if sha256(path) != manifest["input_batches_sha256"]:
        raise ValueError("bivariate batch hash changed")
    values = means(batches)
    deleted_values = [means(batches, index) for index in range(len(batches))]

    basis0 = transformed_basis(alexander_reflection=False, rotation_power=0)
    plus0, plus0_singular = annihilator_line(hand_hankel(values, "plus", basis0))
    deleted_plus0 = [annihilator_line(hand_hankel(row, "plus", basis0))[0] for row in deleted_values]

    candidates = {}
    for alexander in (False, True):
        family = "alexander_reflection" if alexander else "orientation_preserving"
        for rotation_power in range(4):
            basis = transformed_basis(alexander_reflection=alexander, rotation_power=rotation_power)
            minus, singular_values = annihilator_line(hand_hankel(values, "minus", basis))
            deleted_minus = [annihilator_line(hand_hankel(row, "minus", basis))[0] for row in deleted_values]
            for conjugate_plus in (False, True):
                name = f"{family}_R{rotation_power}_{'conjugate' if conjugate_plus else 'linear'}"
                candidates[name] = {
                    "alexander_reflection": alexander,
                    "rotation_power": rotation_power,
                    "target_basis": [list(point) for point in basis],
                    "minus_singular_values": [float(value) for value in singular_values],
                    "minus_line": line_payload(minus, basis),
                    **bridge_score(
                        plus0, minus, deleted_plus0, deleted_minus,
                        conjugate_plus=conjugate_plus,
                    ),
                }

    # Internal C4 covariance is a calibration of the line extraction itself.
    internal = {}
    for hand in HANDS:
        base, _ = annihilator_line(hand_hankel(values, hand, basis0))
        deleted_base = [annihilator_line(hand_hankel(row, hand, basis0))[0] for row in deleted_values]
        for rotation_power in (1, 2, 3):
            basis = transformed_basis(alexander_reflection=False, rotation_power=rotation_power)
            rotated_line, _ = annihilator_line(hand_hankel(values, hand, basis))
            deleted_rotated = [annihilator_line(hand_hankel(row, hand, basis))[0] for row in deleted_values]
            for conjugate_line in (False, True):
                name = f"{hand}_R{rotation_power}_{'conjugate' if conjugate_line else 'linear'}"
                internal[name] = bridge_score(
                    base, rotated_line, deleted_base, deleted_rotated,
                    conjugate_plus=conjugate_line,
                )

    expected = {
        name: row for name, row in candidates.items()
        if row["alexander_reflection"] and row["coefficient_conjugation"]
    }
    surviving_expected = [
        name for name, row in expected.items()
        if row["score"]["finite_batch_survival_p"] >= float(manifest["decision_alpha"])
    ]
    if surviving_expected:
        decision = "Alexander_reflection_conjugation_annihilator_bridge_survives"
    else:
        decision = "all_parameter_free_Alexander_reflection_conjugation_bridges_rejected"
    return {
        "schema": "matching-one/z5-projective-leg-annihilator-bridge/v1",
        "status": "existing_data_two_sector_annihilator_reanalysis",
        "exact_map_gate": exact_map_gate(),
        "source_plus_basis": [list(point) for point in basis0],
        "source_plus_singular_values": [float(value) for value in plus0_singular],
        "source_plus_line": line_payload(plus0, basis0),
        "candidate_maps": candidates,
        "internal_C4_line_calibration": internal,
        "primary_family": "Alexander reflection composed with R^k and complex conjugation, k=0..3",
        "surviving_primary_maps": surviving_expected,
        "decision_alpha": float(manifest["decision_alpha"]),
        "decision": decision,
        "deck_charge_scope": "Deck relabeling and nonzero character phases are invertible Hankel row operations, so all four deck-generator choices induce the identical right-line score.",
        "next_rows_if_rejected": {
            "minimal_degree5_boundary": [[5, 0], [4, 1], [3, 2], [2, 3], [1, 4], [0, 5]],
            "C4_closed_radius5_points": 20,
            "purpose": "shift the degree-two annihilator by every degree-three first-quadrant monomial; only these six degree-five endpoints are new",
            "degree6_requirement": "the full order-two flat-extension matrix requires the degree-six boundary as well",
        },
        "claim_boundary": [
            "A surviving line map is compatibility of truncated sector quotients, not a finite-quotient isomorphism or a completed transfer algebra.",
            "The same-parent plus/minus N505 quotients are not exact D4 reflections, so Alexander is a Hecke-sector hypothesis rather than a microscopic graph map.",
            "A rejected line map excludes only the enumerated parameter-free D4/conjugation bridges; a modulus-dependent intertwiner is not tested.",
        ],
    }


def render(result: Mapping[str, object]) -> str:
    lines = [
        "# P250 two-sector annihilator bridge",
        "",
        "| candidate | Hotelling p | asymptotic p |",
        "|---|---:|---:|",
    ]
    for name, row in result["candidate_maps"].items():
        if row["alexander_reflection"]:
            score_row = row["score"]
            lines.append(
                f"| {name} | {score_row['finite_batch_survival_p']:.6g} | {score_row['asymptotic_survival_p']:.6g} |"
            )
    lines += [
        "",
        f"Surviving primary maps: `{result['surviving_primary_maps']}`.",
        f"Decision: `{result['decision']}`.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    result = score(read_batches(Path(manifest["input_batches"])), manifest)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
