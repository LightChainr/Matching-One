#!/usr/bin/env python3
"""Post-reveal decomposition of the frozen P28 pure-tail rejection."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import mpmath as mp

from score_p28_threshold_profile_tails import (
    DEFAULT_MANIFEST,
    ROOT,
    aggregate_counts,
    correlated_chi_square,
    jackknife_covariance,
    load_histogram,
    profile,
    require,
    sha256,
)


def null_contrasts(design: Sequence[Sequence[float]], candidates: Sequence[Sequence[float]]) -> list[list[float]]:
    size = len(design[0])
    columns: list[list[float]] = []
    for raw in design:
        require(len(raw) == size, "design columns are not aligned")
        value = list(raw)
        for column in columns:
            projection = math.fsum(a * b for a, b in zip(value, column))
            value = [a - projection * b for a, b in zip(value, column)]
        norm = math.sqrt(math.fsum(a * a for a in value))
        require(norm > 1e-12, "degenerate diagnostic design")
        columns.append([a / norm for a in value])
    contrasts: list[list[float]] = []
    for raw in candidates:
        value = list(raw)
        for column in columns + contrasts:
            projection = math.fsum(a * b for a, b in zip(value, column))
            value = [a - projection * b for a, b in zip(value, column)]
        norm = math.sqrt(math.fsum(a * a for a in value))
        if norm > 1e-10:
            contrasts.append([a / norm for a in value])
        if len(contrasts) == size - len(columns):
            break
    require(len(contrasts) == size - len(columns), "candidate modes did not span the residual space")
    return contrasts


def design_for_model(z_grid: Sequence[float], model: str) -> list[list[float]]:
    if model == "stretched_4_over_3":
        return [[1.0] * len(z_grid), [z ** (4 / 3) for z in z_grid]]
    if model == "post_reveal_4_over_3_plus_2_over_3":
        return [
            [1.0] * len(z_grid),
            [z ** (4 / 3) for z in z_grid],
            [z ** (2 / 3) for z in z_grid],
        ]
    raise ValueError(f"unknown model: {model}")


def canonical_modes(z_grid: Sequence[float], model: str) -> list[list[float]]:
    center = math.fsum(z_grid) / len(z_grid)
    centered = [value - center for value in z_grid]
    candidates = [[value**power for value in centered] for power in range(2, 7)]
    return null_contrasts(design_for_model(z_grid, model), candidates)


def curve_vector(profiles: Mapping[str, Mapping[str, Any]], z_grid: Sequence[float], model: str) -> tuple[list[float], list[dict[str, Any]]]:
    modes = canonical_modes(z_grid, model)
    output = []
    labels = []
    for orientation in ("first", "second"):
        for side in ("left", "right"):
            values = profiles[orientation]["log_density"][side]
            for mode, row in enumerate(modes, start=2):
                output.append(math.fsum(weight * value for weight, value in zip(row, values)))
                labels.append({"orientation": orientation, "side": side, "curvature_mode": mode})
    return output, labels


def raw_residual_vector(profiles: Mapping[str, Mapping[str, Any]], z_grid: Sequence[float], model: str) -> tuple[list[float], list[dict[str, Any]]]:
    modes = canonical_modes(z_grid, model)
    output = []
    labels = []
    for orientation in ("first", "second"):
        for side in ("left", "right"):
            values = profiles[orientation]["log_density"][side]
            residual = [
                math.fsum(mode[index] * math.fsum(a * b for a, b in zip(mode, values)) for mode in modes)
                for index in range(len(z_grid))
            ]
            output.extend(residual)
            labels.extend(
                {"orientation": orientation, "side": side, "z": z} for z in z_grid
            )
    return output, labels


def archive_profiles(archive: Mapping[str, Any], z_grid: Sequence[float], omitted: int | None = None) -> dict[str, Any]:
    result = {}
    for orientation in ("first", "second"):
        counts, total = aggregate_counts(archive, orientation, omitted)
        result[orientation] = profile(archive["N"], counts, total, z_grid)
    return result


def subset_score(vector: Sequence[float], covariance: Sequence[Sequence[float]], indices: Sequence[int], cutoff: float) -> dict[str, Any]:
    return correlated_chi_square(
        [vector[index] for index in indices],
        [[covariance[left][right] for right in indices] for left in indices],
        cutoff,
    )


def signed_gls_contributions(residual: Sequence[float], covariance: Sequence[Sequence[float]], cutoff: float) -> dict[str, Any]:
    diagonal = [math.sqrt(max(covariance[i][i], 0.0)) for i in range(len(residual))]
    require(all(value > 0.0 for value in diagonal), "zero raw-residual variance")
    correlation = mp.matrix([
        [covariance[i][j] / (diagonal[i] * diagonal[j]) for j in range(len(residual))]
        for i in range(len(residual))
    ])
    eigenvalues, eigenvectors = mp.eigsy(correlation)
    maximum = max(float(eigenvalues[i]) for i in range(len(residual)))
    retained = [i for i in range(len(residual)) if float(eigenvalues[i]) > cutoff * maximum]
    standardized = mp.matrix([value / scale for value, scale in zip(residual, diagonal)])
    inverse_action = mp.matrix(len(residual), 1)
    for index in retained:
        component = (eigenvectors[:, index].T * standardized)[0]
        inverse_action += eigenvectors[:, index] * (component / eigenvalues[index])
    contributions = [float(standardized[i] * inverse_action[i]) for i in range(len(residual))]
    return {"contributions": contributions, "sum": math.fsum(contributions), "effective_df": len(retained)}


def local_effective_exponents(values: Sequence[float], z_grid: Sequence[float]) -> list[float]:
    slopes = [
        -(values[index + 1] - values[index]) / (z_grid[index + 1] - z_grid[index])
        for index in range(len(z_grid) - 1)
    ]
    require(all(value > 0.0 for value in slopes), "tail density is not decreasing")
    midpoints = [(left + right) / 2 for left, right in zip(z_grid, z_grid[1:])]
    return [
        1.0 + math.log(slopes[index + 1] / slopes[index]) / math.log(midpoints[index + 1] / midpoints[index])
        for index in range(len(slopes) - 1)
    ]


def constant_gls(values: Sequence[float], covariance: Sequence[Sequence[float]], cutoff: float) -> dict[str, Any]:
    matrix = mp.matrix(covariance)
    eigenvalues, eigenvectors = mp.eigsy(matrix)
    maximum = max(float(eigenvalues[i]) for i in range(len(values)))
    retained = [i for i in range(len(values)) if float(eigenvalues[i]) > cutoff * maximum]
    inverse = mp.matrix(len(values))
    for index in retained:
        vector = eigenvectors[:, index]
        inverse += (vector * vector.T) / eigenvalues[index]
    ones = mp.matrix([1.0] * len(values))
    observed = mp.matrix(values)
    estimate = (ones.T * inverse * observed)[0] / (ones.T * inverse * ones)[0]
    residual = observed - estimate * ones
    chi_square = (residual.T * inverse * residual)[0]
    df = max(len(retained) - 1, 0)
    probability = mp.gammainc(mp.mpf(df) / 2, chi_square / 2, mp.inf, regularized=True) if df else mp.mpf(1)
    return {
        "constant_beta_hat": float(estimate),
        "chi_square": float(chi_square),
        "effective_df": df,
        "survival_probability": float(probability),
    }


def diagnose_archive(archive: Mapping[str, Any], z_grid: Sequence[float], cutoff: float, model: str) -> dict[str, Any]:
    full_profiles = archive_profiles(archive, z_grid)
    full, labels = curve_vector(full_profiles, z_grid, model)
    raw, raw_labels = raw_residual_vector(full_profiles, z_grid, model)
    deleted = []
    raw_deleted = []
    local_deleted: dict[tuple[str, str], list[list[float]]] = {
        (orientation, side): []
        for orientation in ("first", "second")
        for side in ("left", "right")
    }
    for omitted in archive["batches"]:
        profiles = archive_profiles(archive, z_grid, omitted)
        deleted.append(curve_vector(profiles, z_grid, model)[0])
        raw_deleted.append(raw_residual_vector(profiles, z_grid, model)[0])
        for key in local_deleted:
            local_deleted[key].append(local_effective_exponents(profiles[key[0]]["log_density"][key[1]], z_grid))
    covariance = jackknife_covariance(deleted)
    raw_covariance = jackknife_covariance(raw_deleted)
    overall = correlated_chi_square(full, covariance, cutoff)
    grouped: dict[str, Any] = {}
    for field in ("side", "orientation", "curvature_mode"):
        grouped[field] = {}
        for value in sorted({str(label[field]) for label in labels}):
            indices = [index for index, label in enumerate(labels) if str(label[field]) == value]
            grouped[field][value] = subset_score(full, covariance, indices, cutoff)
    signed = signed_gls_contributions(raw, raw_covariance, cutoff)
    by_z = {str(z): 0.0 for z in z_grid}
    by_side = {"left": 0.0, "right": 0.0}
    for label, contribution in zip(raw_labels, signed["contributions"]):
        by_z[str(label["z"])] += contribution
        by_side[label["side"]] += contribution
    local = {}
    for orientation in ("first", "second"):
        for side in ("left", "right"):
            key = (orientation, side)
            full_beta = local_effective_exponents(full_profiles[orientation]["log_density"][side], z_grid)
            beta_covariance = jackknife_covariance(local_deleted[key])
            local[f"{orientation}_{side}"] = {
                "beta_effective": full_beta,
                "evaluation_z": [
                    math.sqrt(((z_grid[i] + z_grid[i + 1]) / 2) * ((z_grid[i + 1] + z_grid[i + 2]) / 2))
                    for i in range(len(z_grid) - 2)
                ],
                "jackknife_se": [math.sqrt(beta_covariance[i][i]) for i in range(len(full_beta))],
                "constant_beta_score": constant_gls(full_beta, beta_covariance, cutoff),
            }
    return {
        "N": archive["N"],
        "model": model,
        "overall": overall,
        "marginal_nonadditive_decomposition": grouped,
        "signed_additive_gls_attribution": {
            "warning": "point contributions can be negative under correlated GLS; only grouped sums and the total are invariant",
            "by_z": by_z,
            "by_side": by_side,
            "sum": signed["sum"],
        },
        "local_effective_exponent": local,
    }


def combine_blocks(blocks: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    chi_square = math.fsum(block["overall"]["chi_square"] for block in blocks)
    df = sum(block["overall"]["effective_df"] for block in blocks)
    probability = mp.gammainc(mp.mpf(df) / 2, chi_square / 2, mp.inf, regularized=True)
    return {"chi_square": chi_square, "effective_df": df, "survival_probability": float(probability)}


def run(manifest_path: Path) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    z_grid = [float(value) for value in manifest["tail_design"]["absolute_z_grid"]]
    cutoff = float(manifest["tail_design"]["covariance_relative_eigenvalue_cutoff"])
    archives = []
    for row in manifest["heldout_archives"]:
        path = ROOT / row["path"]
        require(sha256(path) == row["sha256"], f"checksum mismatch: {row['path']}")
        archives.append((row, load_histogram(path)))
    pure_blocks = [diagnose_archive(archive, z_grid, cutoff, "stretched_4_over_3") for _, archive in archives]
    corrected_blocks = [
        diagnose_archive(archive, z_grid, cutoff, "post_reveal_4_over_3_plus_2_over_3")
        for _, archive in archives
    ]
    lineage = {
        "P43": combine_blocks([pure_blocks[0]]),
        "P50": combine_blocks([pure_blocks[1]]),
        "P57": combine_blocks(pure_blocks[2:]),
    }
    return {
        "schema": "matching-one/p28-tail-rejection-diagnostic/v1",
        "issue": 28,
        "status": "post_reveal_diagnostic_no_new_vote",
        "manifest_sha256": sha256(manifest_path),
        "frozen_model": {
            "name": "stretched_4_over_3",
            "global": combine_blocks(pure_blocks),
            "by_lineage": lineage,
            "blocks": pure_blocks,
        },
        "descriptive_nested_correction": {
            "name": "a-c*z^(4/3)+d*z^(2/3)",
            "post_reveal": True,
            "not_a_model_vote": True,
            "global": combine_blocks(corrected_blocks),
            "blocks": corrected_blocks,
        },
        "structural_attribution": (
            "The frozen score fits intercept and decay separately for every size, orientation, and side. "
            "Its rejection therefore cannot be caused by left/right amplitude asymmetry, geometry amplitude, "
            "or cross-size coefficient drift; it is within-tail curvature on the five frozen z points."
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
