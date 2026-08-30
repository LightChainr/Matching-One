#!/usr/bin/env python3
"""Exact root-translation and anchor-factorization gate for Issue 418.

This is deliberately a semantic certificate, not another archive scorer.  It
uses four deterministic period-loop configurations on each N=505 child and
the full order-101 parent translation orbit.  Charged pair identities are
checked in Z[z]/(1+z+...+z^4), so the pass/fail decision does not depend on
floating-point Fourier transforms.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any, Sequence

from norm5_chiral_fixedp_mc import PHASES
from score_p418_crt_degauging import CRT_MULTIPLIER, DECK_ORDER, exact_section_and_masks, residue
from z5_projective_leg_bivariate_mc import rotation_gauges
from z5_projective_leg_cross_scale_mc import PARENT_GEOMETRY, PARENT_MATRIX, contexts
from z5_projective_leg_multiseparation_mc import ProjectiveLegIndex


SCHEMA = "matching-one/p418-root-translation-semantics/v1"
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/exact/P418-root-translation-semantics/certificate.json"
HANDS = ("plus", "minus")
CHARGES = (1, 2)
GROUP_ORDER = PARENT_GEOMETRY.n
ARCHIVES = {
    "radius4": ROOT / "results/huawei-20260830/P250-projective-leg-bivariate-80k/response_80k.json",
    "radius5": ROOT / "results/huawei-20260830/P250-projective-leg-radius5-morphism-1200k/response_1200k.json",
    "radius6": ROOT / "results/huawei-20260830/P250-projective-leg-radius6-flat-1200k/response_1200k.json",
}
SEMANTIC_FILES = (
    "scripts/z5_projective_leg_bivariate_mc.py",
    "scripts/z5_projective_leg_cross_scale_mc.py",
    "scripts/z5_projective_leg_multiseparation_mc.py",
    "scripts/z5_charged_threepoint_mc.py",
    "scripts/norm5_chiral_fixedp_mc.py",
    "scripts/integer_period_torus.py",
)


def add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return left[0] + right[0], left[1] + right[1]


def scale(multiplier: int, point: tuple[int, int]) -> tuple[int, int]:
    return multiplier * point[0], multiplier * point[1]


def phase(exponent: int) -> complex:
    return complex(*PHASES[exponent % DECK_ORDER])


def canonical_cyclotomic(coefficients: Sequence[int]) -> tuple[int, int, int, int]:
    """Canonical coefficients modulo 1+z+...+z^4."""

    if len(coefficients) != DECK_ORDER:
        raise ValueError("a Z5 cyclotomic row needs five coefficients")
    return tuple(int(value - coefficients[-1]) for value in coefficients[:-1])


def cyclotomic_zero(coefficients: Sequence[int]) -> bool:
    return not any(canonical_cyclotomic(coefficients))


def add_coefficients(left: Sequence[int], right: Sequence[int]) -> list[int]:
    return [int(a + b) for a, b in zip(left, right)]


def subtract_coefficients(left: Sequence[int], right: Sequence[int]) -> list[int]:
    return [int(a - b) for a, b in zip(left, right)]


def scale_coefficients(multiplier: int, row: Sequence[int]) -> list[int]:
    return [int(multiplier * value) for value in row]


def shift_coefficients(row: Sequence[int], exponent: int) -> list[int]:
    output = [0] * DECK_ORDER
    for index, value in enumerate(row):
        output[(index + exponent) % DECK_ORDER] += int(value)
    return output


def convolve_coefficients(left: Sequence[int], right: Sequence[int]) -> list[int]:
    output = [0] * DECK_ORDER
    for first, a in enumerate(left):
        for second, b in enumerate(right):
            output[(first + second) % DECK_ORDER] += int(a) * int(b)
    return output


def pair_coefficients(left: Sequence[int], right: Sequence[int], charge: int) -> list[int]:
    """Numerator of F_r(left) F_-r(right), whose denominator is 25."""

    output = [0] * DECK_ORDER
    for fiber, first in enumerate(left):
        for other_fiber, second in enumerate(right):
            exponent = charge * (other_fiber - fiber)
            output[exponent % DECK_ORDER] += int(first) * int(second)
    return output


def cyclotomic_value(coefficients: Sequence[int], denominator: int) -> complex:
    return sum(int(value) * phase(index) for index, value in enumerate(coefficients)) / denominator


def period_column_loop(context, column: int) -> list[bool]:
    vector = (
        context.geometry.periods.matrix[0][column],
        context.geometry.periods.matrix[1][column],
    )
    steps = [(1 if vector[0] > 0 else -1, 0)] * abs(vector[0])
    steps += [(0, 1 if vector[1] > 0 else -1)] * abs(vector[1])
    active = [False] * context.geometry.n
    point = (0, 0)
    for step in steps:
        active[context.geometry.vertex(point)] = True
        point = add(point, step)
    if context.geometry.vertex(point) != context.geometry.vertex((0, 0)):
        raise AssertionError("period-column witness did not close")
    return active


def configurations(context) -> list[tuple[str, list[bool]]]:
    rows = []
    for column in (0, 1):
        active = period_column_loop(context, column)
        rows.append((f"black_period_column_{column}", active))
        rows.append((f"white_period_column_{column}", [not value for value in active]))
    return rows


def translate_active(context, active: Sequence[bool], translation: tuple[int, int]) -> list[bool]:
    output = [False] * context.geometry.n
    for vertex, point in enumerate(context.geometry.coordinates):
        output[context.geometry.vertex(add(point, translation))] = bool(active[vertex])
    return output


def component_signature(index: ProjectiveLegIndex, vertex: int) -> tuple[Any, ...]:
    component = index._component(vertex)
    # A rank-two generator list may depend on union order; rank, size and the
    # canonical primitive line in rank one are the root-selection invariants.
    line = component.basis[0] if component.rank == 1 else None
    return bool(index.active[vertex]), component.rank, component.size, line


def root_sequences(context, index: ProjectiveLegIndex, section_points: Sequence[tuple[int, int]]) -> list[list[int]]:
    deck_step = (PARENT_MATRIX[0][0], PARENT_MATRIX[1][0])
    return [
        [
            index.scalar(context.geometry.vertex(add(point, scale(fiber, deck_step))))
            for fiber in range(DECK_ORDER)
        ]
        for point in section_points
    ]


def exact_gap_payload(coefficients: Sequence[int], denominator: int) -> dict[str, Any]:
    canonical = canonical_cyclotomic(coefficients)
    value = cyclotomic_value(coefficients, denominator)
    return {
        "canonical_Phi5_coefficients": list(canonical),
        "common_denominator": denominator,
        "complex_abs": abs(value),
        "exact_zero": not any(canonical),
    }


def git_blob_sha256(commit: str, path: str) -> str:
    payload = subprocess.check_output(
        ["git", "show", f"{commit}:{path}"], cwd=ROOT
    )
    return hashlib.sha256(payload).hexdigest()


def archive_provenance_gate() -> dict[str, Any]:
    rows = {}
    for name, path in ARCHIVES.items():
        payload = json.loads(path.read_text())
        commit = payload["manifest_runner_commit"]
        rows[name] = {
            "response": str(path.relative_to(ROOT)),
            "manifest_runner_commit": commit,
            "semantic_blob_sha256": {
                source: git_blob_sha256(commit, source) for source in SEMANTIC_FILES
            },
        }
    common_files = {
        source: len({row["semantic_blob_sha256"][source] for row in rows.values()}) == 1
        for source in SEMANTIC_FILES
    }
    radius4_payload = json.loads(ARCHIVES["radius4"].read_text())
    archived_gauges = radius4_payload["exact_gate"]["rotation_fiber_gate"]["covariant_gauge_mod5"]
    current_gauges = rotation_gauges()
    gauge_matches = {
        hand: archived_gauges[hand] == current_gauges[hand] for hand in HANDS
    }
    passed = all(common_files.values()) and all(gauge_matches.values())
    return {
        "archive_runner_commits": rows,
        "semantic_files_identical_across_radius4_radius5_radius6": common_files,
        "radius4_archived_gauge_equals_reconstructed_gauge": gauge_matches,
        "note": "radius5/radius6 response exact_gate payloads omit the imported gauge array, so equality is certified by the identical historical source blobs; radius4 stores the array and also matches reconstruction exactly",
        "passed": passed,
    }


def translation_gate(context, active: Sequence[bool], section_points: Sequence[tuple[int, int]]) -> dict[str, Any]:
    original = ProjectiveLegIndex(context.geometry, active)
    original_sequences = root_sequences(context, original, section_points)
    translated_sequences = []
    signature_failures = 0
    scalar_failures = 0
    first_failure = None
    checks = 0
    for parent_translation, parent_point in enumerate(PARENT_GEOMETRY.coordinates):
        child_translation = scale(CRT_MULTIPLIER, parent_point)
        moved = ProjectiveLegIndex(
            context.geometry,
            translate_active(context, active, child_translation),
        )
        for vertex, point in enumerate(context.geometry.coordinates):
            target = context.geometry.vertex(add(point, child_translation))
            checks += 1
            old_signature = component_signature(original, vertex)
            new_signature = component_signature(moved, target)
            signature_failures += int(old_signature != new_signature)
            old_scalar = original.scalar(vertex)
            new_scalar = moved.scalar(target)
            scalar_failures += int(old_scalar != new_scalar)
            if first_failure is None and (old_signature != new_signature or old_scalar != new_scalar):
                first_failure = {
                    "parent_translation": parent_translation,
                    "vertex": vertex,
                    "target": target,
                    "old_signature": list(old_signature),
                    "new_signature": list(new_signature),
                    "old_scalar": old_scalar,
                    "new_scalar": new_scalar,
                }
        moved_sequences = root_sequences(context, moved, section_points)
        # s_CRT is a homomorphism: a translated configuration evaluated at
        # x+q must reproduce the old root row at x, fiber by fiber.
        for parent, point in enumerate(PARENT_GEOMETRY.coordinates):
            target_parent = PARENT_GEOMETRY.vertex(add(point, parent_point))
            scalar_failures += sum(
                moved_sequences[target_parent][fiber] != original_sequences[parent][fiber]
                for fiber in range(DECK_ORDER)
            )
        translated_sequences.append(moved_sequences)
    return {
        "component_signature_checks": checks,
        "component_signature_failures": signature_failures,
        "CRT_root_scalar_checks": GROUP_ORDER * GROUP_ORDER * DECK_ORDER,
        "root_scalar_failures": scalar_failures,
        "first_failure": first_failure,
        "original_sequences": original_sequences,
        "translated_sequences": translated_sequences,
    }


def section_gauge_gate(
    context,
    index: ProjectiveLegIndex,
    section_points: Sequence[tuple[int, int]],
    offsets: Sequence[int],
    gauge: Sequence[int],
) -> dict[str, Any]:
    crt = root_sequences(context, index, section_points)
    bfs = [
        [index.scalar(context.field_to_vertex[DECK_ORDER * parent + fiber]) for fiber in range(DECK_ORDER)]
        for parent in range(GROUP_ORDER)
    ]
    offset_failures = 0
    charged_failures = 0
    max_complex_residual = 0.0
    for parent in range(GROUP_ORDER):
        offset = offsets[parent]
        expected = [bfs[parent][(offset + fiber) % DECK_ORDER] for fiber in range(DECK_ORDER)]
        offset_failures += sum(a != b for a, b in zip(crt[parent], expected))
        for charge in range(1, DECK_ORDER):
            raw = sum(bfs[parent][fiber] * phase(-charge * fiber) for fiber in range(DECK_ORDER)) / DECK_ORDER
            crt_value = sum(crt[parent][fiber] * phase(-charge * fiber) for fiber in range(DECK_ORDER)) / DECK_ORDER
            expected_crt = phase(charge * offset) * raw
            residual = abs(crt_value - expected_crt)
            max_complex_residual = max(max_complex_residual, residual)
            charged_failures += int(residual > 1e-12)
            stored = phase(charge * gauge[parent]) * raw
            expected_stored = phase(charge * (gauge[parent] - offset)) * crt_value
            residual = abs(stored - expected_stored)
            max_complex_residual = max(max_complex_residual, residual)
            charged_failures += int(residual > 1e-12)
    return {
        "fiber_scalar_checks": GROUP_ORDER * DECK_ORDER,
        "fiber_offset_failures": offset_failures,
        "charged_DFT_checks": GROUP_ORDER * 4 * 2,
        "charged_DFT_failures": charged_failures,
        "max_complex_roundoff_residual": max_complex_residual,
        "passed": offset_failures == 0 and charged_failures == 0,
        "crt_sequences": crt,
        "bfs_sequences": bfs,
    }


def pair_factorization_gate(
    crt_sequences: Sequence[Sequence[int]],
    bfs_sequences: Sequence[Sequence[int]],
    translated_sequences: Sequence[Sequence[Sequence[int]]],
    residual_phase: Sequence[int],
    gauge: Sequence[int],
    displacement_residue: int,
    charge: int,
) -> dict[str, Any]:
    coordinates_by_residue = {
        residue(point): point for point in PARENT_GEOMETRY.coordinates
    }
    displacement = coordinates_by_residue[displacement_residue]
    phase_counts = [0] * DECK_ORDER
    latent_sum = [0] * DECK_ORDER
    stored_sum = [0] * DECK_ORDER
    anchor_rows = []
    pointwise_failures = 0
    for origin, point in enumerate(PARENT_GEOMETRY.coordinates):
        target = PARENT_GEOMETRY.vertex(add(point, displacement))
        exponent = charge * (residual_phase[origin] - residual_phase[target])
        phase_counts[exponent % DECK_ORDER] += 1
        latent = pair_coefficients(crt_sequences[origin], crt_sequences[target], charge)
        stored = shift_coefficients(latent, exponent)
        direct = shift_coefficients(
            pair_coefficients(bfs_sequences[origin], bfs_sequences[target], charge),
            charge * (gauge[origin] - gauge[target]),
        )
        pointwise_failures += int(not cyclotomic_zero(subtract_coefficients(stored, direct)))
        latent_sum = add_coefficients(latent_sum, latent)
        stored_sum = add_coefficients(stored_sum, stored)
        anchor_rows.append(stored)

    factorized_numerator = convolve_coefficients(phase_counts, latent_sum)
    fixed_configuration_gap = subtract_coefficients(
        scale_coefficients(GROUP_ORDER, stored_sum), factorized_numerator
    )
    anchor_zero_gap = subtract_coefficients(
        scale_coefficients(GROUP_ORDER, anchor_rows[0]), stored_sum
    )

    orbit_stored_sum = [0] * DECK_ORDER
    fixed_anchor_orbit_sum = [0] * DECK_ORDER
    for moved_sequences in translated_sequences:
        for origin, point in enumerate(PARENT_GEOMETRY.coordinates):
            target = PARENT_GEOMETRY.vertex(add(point, displacement))
            exponent = charge * (residual_phase[origin] - residual_phase[target])
            latent = pair_coefficients(moved_sequences[origin], moved_sequences[target], charge)
            stored = shift_coefficients(latent, exponent)
            orbit_stored_sum = add_coefficients(orbit_stored_sum, stored)
            if origin == 0:
                fixed_anchor_orbit_sum = add_coefficients(fixed_anchor_orbit_sum, stored)
    orbit_gap = subtract_coefficients(orbit_stored_sum, factorized_numerator)
    local_phase = charge * (
        residual_phase[0]
        - residual_phase[PARENT_GEOMETRY.vertex(add(PARENT_GEOMETRY.coordinates[0], displacement))]
    )
    fixed_anchor_expected = shift_coefficients(latent_sum, local_phase)
    fixed_anchor_orbit_gap = subtract_coefficients(fixed_anchor_orbit_sum, fixed_anchor_expected)

    distinct_anchors = len({canonical_cyclotomic(row) for row in anchor_rows})
    return {
        "charge": charge,
        "displacement_residue": displacement_residue,
        "displacement_coordinate": list(displacement),
        "phase_counts": phase_counts,
        "pointwise_pair_failures": pointwise_failures,
        "one_anchor": {
            "distinct_exact_values_across_101_anchors": distinct_anchors,
            "anchor0_minus_full_anchor_mean": exact_gap_payload(anchor_zero_gap, 25 * GROUP_ORDER),
        },
        "fixed_configuration_factorization": {
            "formula": "mean_x m(x)m(x+d)^* C_config(x,d) versus A(d) mean_x C_config(x,d)",
            "gap": exact_gap_payload(fixed_configuration_gap, 25 * GROUP_ORDER * GROUP_ORDER),
            "expected_to_hold_per_configuration": False,
        },
        "translation_orbit_factorization": {
            "formula": "mean_config_translation mean_anchor stored pair = A(d) mean_anchor latent pair",
            "gap": exact_gap_payload(orbit_gap, 25 * GROUP_ORDER * GROUP_ORDER),
            "passed": cyclotomic_zero(orbit_gap),
        },
        "fixed_anchor_translation_orbit": {
            "formula": "mean_config_translation at one fixed anchor uses that anchor's phase, not A(d)",
            "local_phase_identity_gap": exact_gap_payload(fixed_anchor_orbit_gap, 25 * GROUP_ORDER),
            "passed": cyclotomic_zero(fixed_anchor_orbit_gap),
        },
    }


def build_certificate() -> dict[str, Any]:
    exact = exact_section_and_masks()
    gauges = rotation_gauges()
    provenance = archive_provenance_gate()
    hand_rows = {}
    total_signature_failures = 0
    total_root_failures = 0
    total_section_failures = 0
    total_orbit_failures = 0
    fixed_configuration_counterexamples = []
    for hand, context in zip(HANDS, contexts()):
        exact_hand = exact["hands"][hand]
        section_points = [scale(CRT_MULTIPLIER, point) for point in PARENT_GEOMETRY.coordinates]
        config_rows = {}
        witness_residues = {
            1: exact["attenuation_witnesses"][f"{hand}_r1"]["residue"],
            2: exact["attenuation_witnesses"][f"{hand}_r2"]["residue"],
        }
        for name, active in configurations(context):
            translation = translation_gate(context, active, section_points)
            index = ProjectiveLegIndex(context.geometry, active)
            section = section_gauge_gate(
                context,
                index,
                section_points,
                exact_hand["fiber_offset_b"],
                gauges[hand],
            )
            pair_rows = {}
            for charge in CHARGES:
                row = pair_factorization_gate(
                    section["crt_sequences"],
                    section["bfs_sequences"],
                    translation["translated_sequences"],
                    exact_hand["residual_phase_u"],
                    gauges[hand],
                    witness_residues[charge],
                    charge,
                )
                pair_rows[f"r{charge}"] = row
                total_orbit_failures += int(not row["translation_orbit_factorization"]["passed"])
                total_section_failures += row["pointwise_pair_failures"]
                if not row["fixed_configuration_factorization"]["gap"]["exact_zero"]:
                    fixed_configuration_counterexamples.append(
                        {
                            "hand": hand,
                            "configuration": name,
                            "charge": charge,
                            "displacement_residue": witness_residues[charge],
                            "gap": row["fixed_configuration_factorization"]["gap"],
                        }
                    )
            total_signature_failures += translation["component_signature_failures"]
            total_root_failures += translation["root_scalar_failures"]
            total_section_failures += section["fiber_offset_failures"] + section["charged_DFT_failures"]
            config_rows[name] = {
                "occupied_vertices": sum(active),
                "translation": {key: value for key, value in translation.items() if not key.endswith("sequences")},
                "section_gauge": {key: value for key, value in section.items() if not key.endswith("sequences")},
                "pair_factorization": pair_rows,
            }
        hand_rows[hand] = {
            "context": context.name,
            "configurations": config_rows,
        }

    passed = provenance["passed"] and not any(
        (total_signature_failures, total_root_failures, total_section_failures, total_orbit_failures)
    )
    if not passed:
        raise AssertionError("P418 root-translation semantic gate failed")
    return {
        "schema": SCHEMA,
        "status": "exact_fixed_configuration_semantic_certificate",
        "issues": [418, 250],
        "base_commit": "8704eee790403e14e5ad75d3465ee1496eaa9c0e",
        "new_monte_carlo": False,
        "geometry": {
            "parent_order": GROUP_ORDER,
            "child_order": contexts()[0].geometry.n,
            "translations": "all 101 parent translations lifted by the unique CRT section s_CRT(x)=405*s_BFS(x)",
            "configurations_per_child": 4,
            "configuration_family": "both period-column black loops and their color complements",
        },
        "exact_arithmetic": {
            "root_field": "integer values {-1,0,+1}",
            "charged_pairs": "integer coefficients in Z[zeta5]/(1+zeta5+...+zeta5^4)",
            "floating_point_role": "reported DFT roundoff residual only; no exact decision uses it",
        },
        "archive_provenance": provenance,
        "hands": hand_rows,
        "summary": {
            "component_signature_failures": total_signature_failures,
            "root_scalar_failures": total_root_failures,
            "section_or_gauge_failures": total_section_failures,
            "translation_orbit_factorization_failures": total_orbit_failures,
            "archive_semantic_provenance_failures": int(not provenance["passed"]),
            "fixed_configuration_factorization_counterexamples": len(fixed_configuration_counterexamples),
            "first_fixed_configuration_counterexample": (
                fixed_configuration_counterexamples[0] if fixed_configuration_counterexamples else None
            ),
            "passed": passed,
        },
        "localization": {
            "root_component_selection": "not the break: active color, rank, size and primitive rank-one line are exactly translation covariant",
            "section_and_gauge": "not the break: BFS and CRT roots plus the stored rotation gauge obey the frozen pointwise relation",
            "single_configuration_warning": "A(d) cannot be pulled outside a full-anchor average for one fixed configuration because C_config(x,d) is not translation constant",
            "sampling_contract": "factorization is restored exactly by the full 101x101 configuration-translation x anchor orbit; one random anchor is unbiased only jointly with the translation-invariant configuration law",
            "remaining_P418_target": "the archive one-anchor/covariance/scorer assembly, not the projective-leg definition or CRT gauge algebra",
        },
        "claim_boundary": [
            "The finite translation-orbit identity is an exact synthetic ensemble, not a Monte Carlo precision claim.",
            "A per-configuration factorization counterexample is expected and does not contradict ensemble stationarity.",
            "This certificate does not by itself identify whether the remaining archive mismatch is sampling variance or covariance/scorer assembly.",
        ],
    }


def markdown(certificate: dict[str, Any]) -> str:
    summary = certificate["summary"]
    first = summary["first_fixed_configuration_counterexample"]
    return "\n".join(
        [
            "# P418 root-translation semantic certificate",
            "",
            "No new Monte Carlo was run. Four deterministic nonzero period-loop witnesses per N505 child were translated through all 101 parent positions.",
            "",
            "| gate | failures |",
            "|---|---:|",
            f"| root component signature | {summary['component_signature_failures']} |",
            f"| real root scalar | {summary['root_scalar_failures']} |",
            f"| BFS/CRT section and stored gauge | {summary['section_or_gauge_failures']} |",
            f"| full translation-orbit factorization | {summary['translation_orbit_factorization_failures']} |",
            f"| historical archive section/gauge provenance | {summary['archive_semantic_provenance_failures']} |",
            "",
            f"There are `{summary['fixed_configuration_factorization_counterexamples']}` exact counterexamples to pulling the mask through a *single fixed configuration*. The first is `{first['hand']}/{first['configuration']}/r{first['charge']}/d={first['displacement_residue']}` with canonical Phi5 gap `{first['gap']['canonical_Phi5_coefficients']}`.",
            "",
            "The radius4/radius5/radius6 archives were produced at different runner commits, but all six files that define the section, gauge, cover, root observable and DFT have identical SHA-256 hashes at those commits; the radius4 stored gauge array also equals the reconstruction exactly. This localizes the Issue 418 failure away from root selection, away from the CRT section/gauge, and away from a historical old4 gauge mismatch. The mask factorization is restored exactly on the full 101 configuration translations times 101 anchors. The remaining target is the archive one-anchor/covariance/scorer assembly.",
            "",
            "Boundary: the fixed-configuration counterexample is expected and is not evidence against stationarity; this gate does not decide sampling variance versus covariance/scorer assembly.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args()
    certificate = build_certificate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(certificate, indent=2) + "\n")
    markdown_path = args.markdown or args.output.with_suffix(".md")
    markdown_path.write_text(markdown(certificate))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
