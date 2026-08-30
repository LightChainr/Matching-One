#!/usr/bin/env python3
"""Score frozen tail families against committed P28 threshold histograms."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Mapping, Sequence

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "analysis/p28_threshold_profile_tail_manifest.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_histogram(path: Path) -> dict[str, Any]:
    counts: dict[tuple[str, str, int], dict[int, int]] = defaultdict(lambda: defaultdict(int))
    samples_by_key: dict[tuple[str, str, int], int] = {}
    n_values: set[int] = set()
    batches: set[int] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            n_values.add(int(row["n"]))
            batch = int(row["batch"])
            batches.add(batch)
            key = (row["orientation"], row["kind"], batch)
            counts[key][int(row["k"])] += int(row["count"])
            samples_by_key[key] = int(row["samples"])
    require(len(n_values) == 1, f"{path}: N is not unique")
    require(batches == set(range(len(batches))), f"{path}: batch ids are not consecutive from zero")
    for key, rank_counts in counts.items():
        require(sum(rank_counts.values()) == samples_by_key[key], f"{path}: sparse counts do not sum to samples for {key}")
    for orientation in ("first", "second"):
        for kind in ("minus", "plus"):
            require(
                all((orientation, kind, batch) in counts for batch in batches),
                f"{path}: missing {orientation}/{kind} batch",
            )
    return {"N": n_values.pop(), "batches": tuple(sorted(batches)), "counts": counts}


def aggregate_counts(
    archive: Mapping[str, Any], orientation: str, omitted_batch: int | None = None
) -> tuple[dict[int, int], int]:
    output: dict[int, int] = defaultdict(int)
    for batch in archive["batches"]:
        if batch == omitted_batch:
            continue
        for kind in ("minus", "plus"):
            for rank, count in archive["counts"][(orientation, kind, batch)].items():
                output[rank] += count
    return dict(output), sum(output.values())


def mixture_location_scale(n: int, counts: Mapping[int, int], total: int) -> tuple[float, float]:
    mean = math.fsum(count * rank for rank, count in counts.items()) / (total * (n + 1))
    raw_second = math.fsum(
        count * rank * (rank + 1) for rank, count in counts.items()
    ) / (total * (n + 1) * (n + 2))
    variance = raw_second - mean * mean
    require(variance > 0.0, "mixture variance is not positive")
    return mean, math.sqrt(variance)


def beta_density_basis(n: int, p: float) -> list[float]:
    require(0.0 < p < 1.0, "tail grid left [0,1]")
    values = []
    for rank in range(1, n + 1):
        log_value = (
            math.log(n)
            + math.lgamma(n)
            - math.lgamma(rank)
            - math.lgamma(n - rank + 1)
            + (rank - 1) * math.log(p)
            + (n - rank) * math.log1p(-p)
        )
        values.append(math.exp(log_value))
    return values


def binomial_tail_basis(n: int, p: float) -> list[float]:
    logs = [
        math.lgamma(n + 1)
        - math.lgamma(occupied + 1)
        - math.lgamma(n - occupied + 1)
        + occupied * math.log(p)
        + (n - occupied) * math.log1p(-p)
        for occupied in range(n + 1)
    ]
    maximum = max(logs)
    masses = [math.exp(value - maximum) for value in logs]
    normalization = math.fsum(masses)
    masses = [value / normalization for value in masses]
    tails = [0.0] * (n + 1)
    running = 0.0
    for occupied in range(n, -1, -1):
        running += masses[occupied]
        tails[occupied] = running
    return tails[1:]


def profile(
    n: int, counts: Mapping[int, int], total: int, z_grid: Sequence[float]
) -> dict[str, Any]:
    center, scale = mixture_location_scale(n, counts, total)
    log_density: dict[str, list[float]] = {"left": [], "right": []}
    for side, sign in (("left", -1.0), ("right", 1.0)):
        for z in z_grid:
            p = center + sign * z * scale
            basis = beta_density_basis(n, p)
            density = math.fsum(counts.get(rank, 0) * basis[rank - 1] for rank in range(1, n + 1)) / total
            standardized = scale * density
            require(standardized > 0.0 and math.isfinite(standardized), "nonpositive tail density")
            log_density[side].append(math.log(standardized))
    outer = z_grid[-1]
    left_p = center - outer * scale
    right_p = center + outer * scale
    left_tails = binomial_tail_basis(n, left_p)
    right_tails = binomial_tail_basis(n, right_p)
    left_probability = math.fsum(
        counts.get(rank, 0) * left_tails[rank - 1] for rank in range(1, n + 1)
    ) / total
    right_cdf = math.fsum(
        counts.get(rank, 0) * right_tails[rank - 1] for rank in range(1, n + 1)
    ) / total
    return {
        "center": center,
        "scale": scale,
        "log_density": log_density,
        "tail_probability": {"left": left_probability, "right": 1.0 - right_cdf},
    }


def orthogonal_contrasts(x: Sequence[float]) -> list[list[float]]:
    size = len(x)
    columns: list[list[float]] = []
    for raw in ([1.0] * size, list(x)):
        value = list(raw)
        for column in columns:
            projection = math.fsum(a * b for a, b in zip(value, column))
            value = [a - projection * b for a, b in zip(value, column)]
        norm = math.sqrt(math.fsum(a * a for a in value))
        require(norm > 1e-12, "degenerate regression design")
        columns.append([a / norm for a in value])
    contrasts: list[list[float]] = []
    for index in range(size):
        value = [1.0 if row == index else 0.0 for row in range(size)]
        for column in columns + contrasts:
            projection = math.fsum(a * b for a, b in zip(value, column))
            value = [a - projection * b for a, b in zip(value, column)]
        norm = math.sqrt(math.fsum(a * a for a in value))
        if norm > 1e-10:
            contrasts.append([a / norm for a in value])
    require(len(contrasts) == size - 2, "wrong contrast dimension")
    return contrasts


def contrast_vector(profiles: Mapping[str, Mapping[str, Any]], z_grid: Sequence[float], alpha: float) -> list[float]:
    basis = orthogonal_contrasts([z**alpha for z in z_grid])
    output = []
    for orientation in ("first", "second"):
        for side in ("left", "right"):
            values = profiles[orientation]["log_density"][side]
            output.extend(math.fsum(weight * value for weight, value in zip(row, values)) for row in basis)
    return output


def jackknife_covariance(deleted: Sequence[Sequence[float]]) -> list[list[float]]:
    batches = len(deleted)
    width = len(deleted[0])
    means = [math.fsum(row[column] for row in deleted) / batches for column in range(width)]
    factor = (batches - 1) / batches
    return [
        [
            factor * math.fsum(
                (row[left] - means[left]) * (row[right] - means[right]) for row in deleted
            )
            for right in range(width)
        ]
        for left in range(width)
    ]


def correlated_chi_square(residual: Sequence[float], covariance: Sequence[Sequence[float]], cutoff: float) -> dict[str, Any]:
    diagonal = [math.sqrt(max(covariance[i][i], 0.0)) for i in range(len(residual))]
    require(all(value > 0.0 for value in diagonal), "zero residual variance")
    correlation = mp.matrix([
        [covariance[i][j] / (diagonal[i] * diagonal[j]) for j in range(len(residual))]
        for i in range(len(residual))
    ])
    eigenvalues, eigenvectors = mp.eigsy(correlation)
    maximum = max(float(eigenvalues[i]) for i in range(len(residual)))
    retained = [i for i in range(len(residual)) if float(eigenvalues[i]) > cutoff * maximum]
    standardized = mp.matrix([value / scale for value, scale in zip(residual, diagonal)])
    chi_square = mp.mpf("0")
    for index in retained:
        component = (eigenvectors[:, index].T * standardized)[0]
        chi_square += component * component / eigenvalues[index]
    probability = mp.gammainc(mp.mpf(len(retained)) / 2, chi_square / 2, mp.inf, regularized=True)
    return {
        "chi_square": float(chi_square),
        "effective_df": len(retained),
        "nominal_df": len(residual),
        "survival_probability": float(probability),
        "correlation_eigenvalues": [float(eigenvalues[i]) for i in range(len(residual))],
    }


def score_archive(archive: Mapping[str, Any], z_grid: Sequence[float], alpha: float, cutoff: float) -> dict[str, Any]:
    full_profiles = {}
    for orientation in ("first", "second"):
        counts, total = aggregate_counts(archive, orientation)
        full_profiles[orientation] = profile(archive["N"], counts, total, z_grid)
    full_vector = contrast_vector(full_profiles, z_grid, alpha)
    deleted_vectors = []
    for omitted in archive["batches"]:
        deleted_profiles = {}
        for orientation in ("first", "second"):
            counts, total = aggregate_counts(archive, orientation, omitted)
            deleted_profiles[orientation] = profile(archive["N"], counts, total, z_grid)
        deleted_vectors.append(contrast_vector(deleted_profiles, z_grid, alpha))
    covariance = jackknife_covariance(deleted_vectors)
    score = correlated_chi_square(full_vector, covariance, cutoff)
    score["N"] = archive["N"]
    score["profiles"] = full_profiles
    return score


def gate_archive(archive: Mapping[str, Any], z_grid: Sequence[float], minimum_total: int, minimum_batch: int) -> dict[str, Any]:
    rows = []
    passed = True
    for orientation in ("first", "second"):
        counts, total = aggregate_counts(archive, orientation)
        observed = profile(archive["N"], counts, total, z_grid)
        batches = len(archive["batches"])
        for side in ("left", "right"):
            expected = total * observed["tail_probability"][side]
            row_passed = expected >= minimum_total and expected / batches >= minimum_batch
            passed = passed and row_passed
            rows.append({
                "orientation": orientation,
                "side": side,
                "expected_tail_count": expected,
                "expected_tail_count_per_batch": expected / batches,
                "passed": row_passed,
            })
    return {"N": archive["N"], "passed": passed, "rows": rows}


def run(manifest_path: Path) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    require(manifest["schema"] == "matching-one/p28-production-tail-elimination/v1", "wrong manifest schema")
    design = manifest["tail_design"]
    z_grid = [float(value) for value in design["absolute_z_grid"]]
    source_gates = []
    heldout_gates = []
    heldout = []
    for role, gate_rows in (("source_archives", source_gates), ("heldout_archives", heldout_gates)):
        for row in manifest[role]:
            path = ROOT / row["path"]
            require(sha256(path) == row["sha256"], f"checksum mismatch: {row['path']}")
            archive = load_histogram(path)
            require(archive["N"] == row["N"], f"N mismatch: {row['path']}")
            gate_rows.append(gate_archive(
                archive,
                z_grid,
                int(design["minimum_expected_tail_count_per_curve"]),
                int(design["minimum_expected_tail_count_per_batch"]),
            ))
            if role == "heldout_archives":
                heldout.append((row, archive))
    model_results = {}
    for model, alpha in design["models"].items():
        blocks = [
            score_archive(archive, z_grid, float(alpha), float(design["covariance_relative_eigenvalue_cutoff"]))
            for _, archive in heldout
        ]
        chi_square = math.fsum(block["chi_square"] for block in blocks)
        effective_df = sum(block["effective_df"] for block in blocks)
        nominal_df = sum(block["nominal_df"] for block in blocks)
        probability = float(mp.gammainc(mp.mpf(effective_df) / 2, chi_square / 2, mp.inf, regularized=True))
        powered = all(row["passed"] for row in heldout_gates) and effective_df >= math.ceil(0.75 * nominal_df)
        if not powered:
            status = "underpowered"
        elif probability < float(design["decision_alpha"]):
            status = "eliminated"
        else:
            status = "survives"
        model_results[model] = {
            "fixed_tail_exponent": alpha,
            "status": status,
            "chi_square": chi_square,
            "effective_df": effective_df,
            "nominal_df": nominal_df,
            "survival_probability": probability,
            "blocks": blocks,
        }
    return {
        "schema": "matching-one/p28-production-tail-elimination-result/v1",
        "issue": 28,
        "manifest": str(manifest_path.relative_to(ROOT)),
        "manifest_sha256": sha256(manifest_path),
        "base_commit": manifest["base_commit"],
        "source_count_gates": source_gates,
        "heldout_count_gates": heldout_gates,
        "models": model_results,
        "covariance": manifest["covariance"],
        "interpretation_level": manifest["interpretation_level"],
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
