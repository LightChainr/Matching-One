#!/usr/bin/env python3
"""Decompose P28 composite-tail curvature into K1/K2 birth clocks."""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

import mpmath as mp

from diagnose_p28_threshold_profile_tail_rejection import (
    canonical_modes,
    constant_gls,
    local_effective_exponents,
)
from score_p28_threshold_profile_tails import (
    DEFAULT_MANIFEST,
    ROOT,
    aggregate_counts,
    beta_density_basis,
    binomial_tail_basis,
    correlated_chi_square,
    jackknife_covariance,
    load_histogram,
    profile,
    require,
    sha256,
)


CLOCKS = {"K1": "minus", "K2": "plus"}


def aggregate_clock(
    archive: Mapping[str, Any], orientation: str, clock: str, omitted_batch: int | None = None
) -> tuple[dict[int, int], int]:
    kind = CLOCKS[clock]
    output: dict[int, int] = defaultdict(int)
    for batch in archive["batches"]:
        if batch == omitted_batch:
            continue
        for rank, count in archive["counts"][(orientation, kind, batch)].items():
            output[rank] += count
    return dict(output), sum(output.values())


def cdf_at(n: int, counts: Mapping[int, int], total: int, p: float) -> float:
    basis = binomial_tail_basis(n, p)
    return math.fsum(counts.get(rank, 0) * basis[rank - 1] for rank in range(1, n + 1)) / total


def density_at(n: int, counts: Mapping[int, int], total: int, p: float) -> float:
    basis = beta_density_basis(n, p)
    return math.fsum(counts.get(rank, 0) * basis[rank - 1] for rank in range(1, n + 1)) / total


def enrich_profile(
    n: int, counts: Mapping[int, int], total: int, z_grid: Sequence[float]
) -> dict[str, Any]:
    result = profile(n, counts, total, z_grid)
    result["sample_count"] = total
    result["expected_tail_count"] = {
        side: total * result["tail_probability"][side] for side in ("left", "right")
    }
    result["standardized_density"] = {
        side: [math.exp(value) for value in result["log_density"][side]]
        for side in ("left", "right")
    }
    result["cdf"] = {}
    for side, sign in (("left", -1.0), ("right", 1.0)):
        result["cdf"][side] = [
            cdf_at(n, counts, total, result["center"] + sign * z * result["scale"])
            for z in z_grid
        ]
    return result


def component_profiles(
    archive: Mapping[str, Any], z_grid: Sequence[float], omitted: int | None = None
) -> dict[str, Any]:
    output = {}
    for clock in CLOCKS:
        output[clock] = {}
        for orientation in ("first", "second"):
            counts, total = aggregate_clock(archive, orientation, clock, omitted)
            output[clock][orientation] = enrich_profile(archive["N"], counts, total, z_grid)
    return output


def component_vector(
    profiles: Mapping[str, Any], z_grid: Sequence[float]
) -> tuple[list[float], list[dict[str, Any]]]:
    modes = canonical_modes(z_grid, "stretched_4_over_3")
    output = []
    labels = []
    for clock in CLOCKS:
        for orientation in ("first", "second"):
            for side in ("left", "right"):
                values = profiles[clock][orientation]["log_density"][side]
                for mode, row in enumerate(modes, start=2):
                    output.append(math.fsum(weight * value for weight, value in zip(row, values)))
                    labels.append({
                        "clock": clock,
                        "orientation": orientation,
                        "side": side,
                        "curvature_mode": mode,
                    })
    return output, labels


def log_mixture_decomposition(first_density: float, second_density: float) -> dict[str, float]:
    require(first_density > 0.0 and second_density > 0.0, "component densities must be positive")
    total = first_density + second_density
    first_weight = first_density / total
    second_weight = second_density / total
    component_shape = first_weight * math.log(first_density) + second_weight * math.log(second_density)
    separation_entropy = (
        -first_weight * math.log(first_weight)
        - second_weight * math.log(second_weight)
        - math.log(2.0)
    )
    composite = math.log(total / 2.0)
    return {
        "log_composite": composite,
        "responsibility_K1": first_weight,
        "responsibility_K2": second_weight,
        "component_shape": component_shape,
        "separation_entropy": separation_entropy,
        "reconstruction_error": composite - component_shape - separation_entropy,
    }


def mixture_decomposition(
    archive: Mapping[str, Any], z_grid: Sequence[float], omitted: int | None = None
) -> dict[str, Any]:
    output = {}
    max_density_error = 0.0
    max_cdf_error = 0.0
    max_log_error = 0.0
    for orientation in ("first", "second"):
        combined_counts, combined_total = aggregate_counts(archive, orientation, omitted)
        composite = profile(archive["N"], combined_counts, combined_total, z_grid)
        clock_data = {}
        for clock in CLOCKS:
            clock_data[clock] = aggregate_clock(archive, orientation, clock, omitted)
        output[orientation] = {}
        for side, sign in (("left", -1.0), ("right", 1.0)):
            rows = []
            for index, z in enumerate(z_grid):
                p = composite["center"] + sign * z * composite["scale"]
                densities = {
                    clock: composite["scale"] * density_at(
                        archive["N"], clock_data[clock][0], clock_data[clock][1], p
                    )
                    for clock in CLOCKS
                }
                cdfs = {
                    clock: cdf_at(archive["N"], clock_data[clock][0], clock_data[clock][1], p)
                    for clock in CLOCKS
                }
                decomposition = log_mixture_decomposition(densities["K1"], densities["K2"])
                reconstructed_density = (densities["K1"] + densities["K2"]) / 2
                observed_density = math.exp(composite["log_density"][side][index])
                reconstructed_cdf = (cdfs["K1"] + cdfs["K2"]) / 2
                observed_cdf = cdf_at(archive["N"], combined_counts, combined_total, p)
                max_density_error = max(max_density_error, abs(reconstructed_density - observed_density))
                max_cdf_error = max(max_cdf_error, abs(reconstructed_cdf - observed_cdf))
                max_log_error = max(max_log_error, abs(decomposition["reconstruction_error"]))
                rows.append({
                    "z": z,
                    "p": p,
                    "standardized_density": densities,
                    "cdf": cdfs,
                    **decomposition,
                })
            output[orientation][side] = rows
    return {
        "rows": output,
        "reconstruction_gate": {
            "maximum_density_absolute_error": max_density_error,
            "maximum_cdf_absolute_error": max_cdf_error,
            "maximum_log_identity_error": max_log_error,
            "passed": max(max_density_error, max_cdf_error, max_log_error) < 1e-12,
        },
    }


def decomposition_vectors(
    decomposition: Mapping[str, Any], z_grid: Sequence[float]
) -> tuple[dict[str, list[float]], list[dict[str, Any]]]:
    modes = canonical_modes(z_grid, "stretched_4_over_3")
    vectors = {"composite": [], "component_shape": [], "separation_entropy": []}
    labels = []
    for orientation in ("first", "second"):
        for side in ("left", "right"):
            rows = decomposition["rows"][orientation][side]
            values = {
                "composite": [row["log_composite"] for row in rows],
                "component_shape": [row["component_shape"] for row in rows],
                "separation_entropy": [row["separation_entropy"] for row in rows],
            }
            for mode, contrast in enumerate(modes, start=2):
                for name in vectors:
                    vectors[name].append(
                        math.fsum(weight * value for weight, value in zip(contrast, values[name]))
                    )
                labels.append({"orientation": orientation, "side": side, "curvature_mode": mode})
    return vectors, labels


def inverse_action(
    target: Sequence[float], covariance: Sequence[Sequence[float]], cutoff: float
) -> tuple[list[float], int]:
    diagonal = [math.sqrt(max(covariance[i][i], 0.0)) for i in range(len(target))]
    correlation = mp.matrix([
        [covariance[i][j] / (diagonal[i] * diagonal[j]) for j in range(len(target))]
        for i in range(len(target))
    ])
    eigenvalues, eigenvectors = mp.eigsy(correlation)
    maximum = max(float(eigenvalues[i]) for i in range(len(target)))
    retained = [i for i in range(len(target)) if float(eigenvalues[i]) > cutoff * maximum]
    standardized = mp.matrix([value / scale for value, scale in zip(target, diagonal)])
    action = mp.matrix(len(target), 1)
    for index in retained:
        vector = eigenvectors[:, index]
        action += vector * ((vector.T * standardized)[0] / eigenvalues[index])
    return [float(action[i] / diagonal[i]) for i in range(len(target))], len(retained)


def bilinear_attribution(
    total: Sequence[float], component: Sequence[float], separation: Sequence[float], covariance: Sequence[Sequence[float]], cutoff: float
) -> dict[str, Any]:
    action, rank = inverse_action(total, covariance, cutoff)
    total_chi = math.fsum(left * right for left, right in zip(total, action))
    component_value = math.fsum(left * right for left, right in zip(component, action))
    separation_value = math.fsum(left * right for left, right in zip(separation, action))
    return {
        "composite_chi_square": total_chi,
        "component_shape_contribution": component_value,
        "separation_entropy_contribution": separation_value,
        "additivity_error": total_chi - component_value - separation_value,
        "effective_df": rank,
        "warning": "bilinear GLS contributions are exactly additive but can be signed",
    }


def diagnose_archive(
    archive: Mapping[str, Any],
    z_grid: Sequence[float],
    cutoff: float,
    minimum_total: int,
    minimum_batch: int,
) -> dict[str, Any]:
    full_components = component_profiles(archive, z_grid)
    full_vector, labels = component_vector(full_components, z_grid)
    full_decomposition = mixture_decomposition(archive, z_grid)
    full_mix_vectors, mix_labels = decomposition_vectors(full_decomposition, z_grid)
    deleted_components = []
    deleted_composite = []
    local_deleted: dict[tuple[str, str, str], list[list[float]]] = {
        (clock, orientation, side): []
        for clock in CLOCKS
        for orientation in ("first", "second")
        for side in ("left", "right")
    }
    for omitted in archive["batches"]:
        profiles = component_profiles(archive, z_grid, omitted)
        deleted_components.append(component_vector(profiles, z_grid)[0])
        decomposition = mixture_decomposition(archive, z_grid, omitted)
        deleted_composite.append(decomposition_vectors(decomposition, z_grid)[0]["composite"])
        for key in local_deleted:
            local_deleted[key].append(
                local_effective_exponents(profiles[key[0]][key[1]]["log_density"][key[2]], z_grid)
            )
    component_covariance = jackknife_covariance(deleted_components)
    composite_covariance = jackknife_covariance(deleted_composite)
    component_scores = {}
    component_side_scores = {}
    for clock in CLOCKS:
        indices = [index for index, label in enumerate(labels) if label["clock"] == clock]
        component_scores[clock] = correlated_chi_square(
            [full_vector[index] for index in indices],
            [[component_covariance[left][right] for right in indices] for left in indices],
            cutoff,
        )
        component_side_scores[clock] = {}
        for side in ("left", "right"):
            side_indices = [
                index for index, label in enumerate(labels)
                if label["clock"] == clock and label["side"] == side
            ]
            component_side_scores[clock][side] = correlated_chi_square(
                [full_vector[index] for index in side_indices],
                [[component_covariance[left][right] for right in side_indices] for left in side_indices],
                cutoff,
            )
    count_gates = []
    for clock in CLOCKS:
        for orientation in ("first", "second"):
            for side in ("left", "right"):
                expected = full_components[clock][orientation]["expected_tail_count"][side]
                per_batch = expected / len(archive["batches"])
                count_gates.append({
                    "clock": clock,
                    "orientation": orientation,
                    "side": side,
                    "expected_tail_count": expected,
                    "expected_tail_count_per_batch": per_batch,
                    "passed": expected >= minimum_total and per_batch >= minimum_batch,
                })
    local = {}
    for key, deleted in local_deleted.items():
        clock, orientation, side = key
        values = local_effective_exponents(
            full_components[clock][orientation]["log_density"][side], z_grid
        )
        covariance = jackknife_covariance(deleted)
        local[f"{clock}_{orientation}_{side}"] = {
            "beta_effective": values,
            "jackknife_se": [math.sqrt(covariance[i][i]) for i in range(len(values))],
            "constant_beta_score": constant_gls(values, covariance, cutoff),
        }
    return {
        "N": archive["N"],
        "component_standardized_profiles": full_components,
        "component_count_gates": count_gates,
        "component_pure_4_over_3_scores": component_scores,
        "component_pure_4_over_3_side_scores": component_side_scores,
        "component_local_effective_exponents": local,
        "mixture_on_composite_coordinates": full_decomposition,
        "mixture_residual_attribution": bilinear_attribution(
            full_mix_vectors["composite"],
            full_mix_vectors["component_shape"],
            full_mix_vectors["separation_entropy"],
            composite_covariance,
            cutoff,
        ),
        "mixture_labels": mix_labels,
    }


def combine_component_scores(blocks: Sequence[Mapping[str, Any]], clock: str) -> dict[str, Any]:
    chi_square = math.fsum(block["component_pure_4_over_3_scores"][clock]["chi_square"] for block in blocks)
    df = sum(block["component_pure_4_over_3_scores"][clock]["effective_df"] for block in blocks)
    probability = mp.gammainc(mp.mpf(df) / 2, chi_square / 2, mp.inf, regularized=True)
    return {"chi_square": chi_square, "effective_df": df, "survival_probability": float(probability)}


def run(manifest_path: Path) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    z_grid = [float(value) for value in manifest["tail_design"]["absolute_z_grid"]]
    cutoff = float(manifest["tail_design"]["covariance_relative_eigenvalue_cutoff"])
    minimum_total = int(manifest["tail_design"]["minimum_expected_tail_count_per_curve"])
    minimum_batch = int(manifest["tail_design"]["minimum_expected_tail_count_per_batch"])
    blocks = []
    for row in manifest["heldout_archives"]:
        path = ROOT / row["path"]
        require(sha256(path) == row["sha256"], f"checksum mismatch: {row['path']}")
        blocks.append(
            diagnose_archive(
                load_histogram(path), z_grid, cutoff, minimum_total, minimum_batch
            )
        )
    attribution = {
        field: math.fsum(block["mixture_residual_attribution"][field] for block in blocks)
        for field in (
            "composite_chi_square",
            "component_shape_contribution",
            "separation_entropy_contribution",
            "additivity_error",
        )
    }
    return {
        "schema": "matching-one/p28-birth-clock-mixture-diagnostic/v1",
        "issue": 28,
        "related_issue": 337,
        "status": "post_reveal_mechanism_diagnostic_no_new_model_vote",
        "clock_convention": {"K1": "histogram kind minus; first ambient-H1 birth", "K2": "histogram kind plus; second ambient-H1 birth"},
        "manifest_sha256": sha256(manifest_path),
        "component_global_pure_4_over_3_scores": {
            clock: combine_component_scores(blocks, clock) for clock in CLOCKS
        },
        "global_mixture_residual_attribution": attribution,
        "component_count_gate_summary": {
            "all_passed": all(
                gate["passed"] for block in blocks for gate in block["component_count_gates"]
            ),
            "failed": [
                {"N": block["N"], **gate}
                for block in blocks
                for gate in block["component_count_gates"]
                if not gate["passed"]
            ],
            "boundary": (
                "Component scores are post-reveal diagnostics. A clock/side cell below the frozen "
                "composite count gate cannot independently carry a formal elimination."
            ),
        },
        "blocks": blocks,
        "interpretation": (
            "Self-standardized K1 and K2 scores test intrinsic clock curvature.  The exact log-mixture "
            "identity separates responsibility-weighted component shape from the entropy of clock separation "
            "on the composite coordinates; its GLS contributions use the frozen composite covariance."
        ),
        "new_monte_carlo": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = run(args.manifest.resolve())
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
