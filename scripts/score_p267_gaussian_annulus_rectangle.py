#!/usr/bin/env python3
"""Frozen heldout scorer for the P267 Gaussian x annulus context rectangle.

The Gaussian input is the future general-period P253-compatible batch stream.
The annulus input is the already-revealed P253 contrast vector.  The contexts
are independent; complete covariance is retained inside each context.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import yaml


LAMBDAS = (0.0, 0.5, 1.0)
LAMBDA_LABELS = ("0", "1/2", "1")
CHANNELS = ("A_plus", "A_minus")
RADII = (2, 4, 7, 8)
GAUSSIAN_SEED = 26725360829
GAUSSIAN_COUNTER = (26725300000, 26725500000)
BOOTSTRAP_SEED = 267253255
ANNULUS_SOURCE_COMMIT = "3123b73"
ANNULUS_SOURCE_PATH = "results/server-20260829/P225-norm5-multiradius/analysis.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def jackknife_covariance(replicates: np.ndarray) -> np.ndarray:
    count = replicates.shape[0]
    centered = replicates - replicates.mean(axis=0)
    return (count - 1) / count * centered.T @ centered


def gaussian_row(lam: float) -> np.ndarray:
    # x3-2*x2+x1-lambda*(x2-2*x1+x0)
    return np.asarray((-lam, 1 + 2 * lam, -2 - lam, 1), dtype=float)


def radial_basis(lam: float, coordinate: float) -> np.ndarray:
    if lam == 0:
        third = 1.0 if abs(coordinate) < 1e-14 else 0.0
    elif lam == 1:
        third = coordinate * coordinate
    else:
        third = lam ** coordinate
    return np.asarray((1.0, coordinate, third), dtype=float)


def annulus_row(lam: float) -> np.ndarray:
    coordinate = np.log2(np.asarray(RADII, dtype=float) / RADII[0])
    calibration = np.stack([radial_basis(lam, value) for value in coordinate[:3]])
    interpolation = radial_basis(lam, coordinate[3]) @ np.linalg.inv(calibration)
    return np.r_[-interpolation, 1.0]


def residual_matrix(lambda_gaussian: float, lambda_annulus: float) -> np.ndarray:
    matrix = np.zeros((6, 24), dtype=float)
    grow = gaussian_row(lambda_gaussian)
    arow = annulus_row(lambda_annulus)
    for series in range(4):
        matrix[series, 4 * series:4 * (series + 1)] = grow
    for series in range(2):
        matrix[4 + series, 16 + 4 * series:16 + 4 * (series + 1)] = arow
    return matrix


def psd_inverse(matrix: np.ndarray, rcond: float = 1e-11) -> tuple[np.ndarray, int]:
    symmetric = (matrix + matrix.T) / 2
    values, vectors = np.linalg.eigh(symmetric)
    scale = max(float(np.max(np.abs(values))), 1.0e-300)
    if float(np.min(values)) < -1e-8 * scale:
        raise ValueError("covariance has a material negative eigenvalue")
    keep = values > rcond * scale
    inverse = (vectors[:, keep] / values[keep]) @ vectors[:, keep].T
    return inverse, int(np.count_nonzero(keep))


def chi_square_survival_even(value: float, degrees: int) -> float:
    if degrees <= 0 or degrees % 2:
        raise ValueError("even positive chi-square degrees required")
    half = value / 2
    return math.exp(-half) * sum(
        half ** power / math.factorial(power) for power in range(degrees // 2)
    )


def fixed_score(point: np.ndarray, covariance: np.ndarray,
                lambda_gaussian: float, lambda_annulus: float) -> dict[str, Any]:
    transform = residual_matrix(lambda_gaussian, lambda_annulus)
    residual = transform @ point
    residual_covariance = transform @ covariance @ transform.T
    precision, rank = psd_inverse(residual_covariance)
    total = float(residual @ precision @ residual)
    gaussian_covariance = residual_covariance[:4, :4]
    annulus_covariance = residual_covariance[4:, 4:]
    gaussian_precision, gaussian_rank = psd_inverse(gaussian_covariance)
    annulus_precision, annulus_rank = psd_inverse(annulus_covariance)
    return {
        "lambda_Gaussian": lambda_gaussian,
        "lambda_annulus": lambda_annulus,
        "residual_order": [
            "Gaussian_lineage65_A_plus", "Gaussian_lineage65_A_minus",
            "Gaussian_lineage85_A_plus", "Gaussian_lineage85_A_minus",
            "annulus_N425_A_plus", "annulus_N425_A_minus",
        ],
        "residual": residual.tolist(),
        "residual_covariance": residual_covariance.tolist(),
        "chi_square": total,
        "effective_rank": rank,
        "chi_square_survival_reference": chi_square_survival_even(total, rank),
        "context_components": {
            "Gaussian": {
                "chi_square": float(residual[:4] @ gaussian_precision @ residual[:4]),
                "effective_rank": gaussian_rank,
            },
            "annulus": {
                "chi_square": float(residual[4:] @ annulus_precision @ residual[4:]),
                "effective_rank": annulus_rank,
            },
        },
    }


def all_fixed_scores(point: np.ndarray, covariance: np.ndarray) -> dict[str, Any]:
    grid: dict[str, dict[str, dict[str, Any]]] = {}
    score_values = np.zeros((3, 3), dtype=float)
    for first, (label_g, lambda_g) in enumerate(zip(LAMBDA_LABELS, LAMBDAS)):
        grid[label_g] = {}
        for second, (label_a, lambda_a) in enumerate(zip(LAMBDA_LABELS, LAMBDAS)):
            score = fixed_score(point, covariance, lambda_g, lambda_a)
            grid[label_g][label_a] = score
            score_values[first, second] = score["chi_square"]
    diagonal = np.diag(score_values)
    best_shared = int(np.argmin(diagonal))
    best_pair_flat = int(np.argmin(score_values))
    best_pair = np.unravel_index(best_pair_flat, score_values.shape)
    return {
        "grid": grid,
        "shared_diagonal": [grid[label][label] for label in LAMBDA_LABELS],
        "best_shared_lambda": LAMBDA_LABELS[best_shared],
        "best_context_pair": {
            "Gaussian": LAMBDA_LABELS[best_pair[0]],
            "annulus": LAMBDA_LABELS[best_pair[1]],
        },
        "minimum_shared_chi_square": float(diagonal[best_shared]),
        "minimum_context_pair_chi_square": float(score_values[best_pair]),
        "delta_min_shared_minus_min_context_pair": float(
            diagonal[best_shared] - score_values[best_pair]
        ),
    }


def null_projection(point: np.ndarray, covariance: np.ndarray, lam: float) -> np.ndarray:
    transform = residual_matrix(lam, lam)
    residual_covariance = transform @ covariance @ transform.T
    precision, _ = psd_inverse(residual_covariance)
    return point - covariance @ transform.T @ precision @ (transform @ point)


def psd_factor(covariance: np.ndarray) -> np.ndarray:
    symmetric = (covariance + covariance.T) / 2
    values, vectors = np.linalg.eigh(symmetric)
    scale = max(float(np.max(np.abs(values))), 1.0e-300)
    if float(np.min(values)) < -1e-8 * scale:
        raise ValueError("base covariance has a material negative eigenvalue")
    values = np.clip(values, 0, None)
    return vectors * np.sqrt(values)


def bootstrap_delta(point: np.ndarray, covariance: np.ndarray, observed_delta: float,
                    draws: int, seed: int = BOOTSTRAP_SEED) -> dict[str, Any]:
    if draws <= 0:
        return {
            "draws": 0,
            "seed": seed,
            "status": "not_run",
        }
    factor = psd_factor(covariance)
    rng = np.random.default_rng(seed)
    results: dict[str, Any] = {}
    for label, lam in zip(LAMBDA_LABELS, LAMBDAS):
        mean = null_projection(point, covariance, lam)
        exceed = 0
        generated = 0
        while generated < draws:
            count = min(10000, draws - generated)
            simulated = mean + rng.standard_normal((count, point.size)) @ factor.T
            score_grid = np.zeros((count, 3, 3), dtype=float)
            for first, lambda_g in enumerate(LAMBDAS):
                for second, lambda_a in enumerate(LAMBDAS):
                    transform = residual_matrix(lambda_g, lambda_a)
                    residual_covariance = transform @ covariance @ transform.T
                    precision, _ = psd_inverse(residual_covariance)
                    residual = simulated @ transform.T
                    score_grid[:, first, second] = np.einsum(
                        "bi,ij,bj->b", residual, precision, residual
                    )
            best_shared = np.minimum.reduce(
                (score_grid[:, 0, 0], score_grid[:, 1, 1], score_grid[:, 2, 2])
            )
            delta = best_shared - score_grid.min(axis=(1, 2))
            exceed += int(np.count_nonzero(delta >= observed_delta - 1e-12))
            generated += count
        results[label] = {
            "null": f"shared lambda={label}",
            "exceedances": exceed,
            "p_value_plus_one": (exceed + 1) / (draws + 1),
        }
    worst = max(row["p_value_plus_one"] for row in results.values())
    return {
        "draws": draws,
        "seed": seed,
        "method": (
            "fixed-covariance Gaussian parametric bootstrap after covariance-metric "
            "projection onto each composite shared-lambda null"
        ),
        "per_shared_null": results,
        "worst_case_p_value": worst,
    }


def manifest_designs(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        design
        for lineage in manifest["missing_Gaussian_acquisition"]["cover_lineages"]
        for design in lineage["designs"]
    ]


def read_gaussian_batches(batch_path: Path, metadata_path: Path,
                          manifest: dict[str, Any], require_production: bool
                          ) -> tuple[np.ndarray, np.ndarray, list[str], dict[str, Any]]:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("schema") != "matching-one/general-period-multiradius-pivotal/v1":
        raise ValueError("Gaussian metadata schema mismatch")
    expected_designs = manifest_designs(manifest)
    observed_designs = metadata.get("designs")
    expected_metadata = [
        {
            "label": row["label"], "N": row["N"],
            "period_matrix": row["matrix"], "smith_invariants": row["smith"],
        }
        for row in expected_designs
    ]
    observed_projection = [
        {key: row[key] for key in expected_metadata[0]}
        for row in observed_designs
    ] if observed_designs else []
    if observed_projection != expected_metadata:
        raise ValueError("Gaussian design order/matrix/Smith contract mismatch")
    if metadata.get("seed") != GAUSSIAN_SEED:
        raise ValueError("Gaussian seed mismatch")
    if (metadata.get("replica_counter_first"),
        metadata.get("replica_counter_last_exclusive")) != GAUSSIAN_COUNTER:
        raise ValueError("Gaussian counter domain mismatch")
    if metadata.get("radii") != [2] or metadata.get("cutoff") != "euclidean":
        raise ValueError("Gaussian radius/cutoff contract mismatch")
    if not math.isclose(float(metadata.get("p")), 0.592746050790,
                        rel_tol=0, abs_tol=5e-13):
        raise ValueError("Gaussian p mismatch")
    if require_production:
        if metadata.get("samples_per_design") != 200000 or metadata.get("batches") != 200:
            raise ValueError("production Gaussian stream requires 200k/200 batches")
        binary_hash = metadata.get("binary_sha256", "")
        if len(binary_hash) != 64 or any(ch not in "0123456789abcdef" for ch in binary_hash):
            raise ValueError("production metadata requires frozen lowercase binary SHA256")

    integer_fields = (
        "n", "m00", "m01", "m10", "m11", "smith1", "smith2", "radius",
        "batch", "counter_first", "counter_last_exclusive", "samples",
        "common_field_digest", "primal_pivotal", "matching_pivotal",
        "primal_h4", "matching_h4", "h4_plus", "h4_minus",
    )
    rows: dict[tuple[str, int], dict[str, int]] = {}
    with batch_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not {"label", *integer_fields}.issubset(reader.fieldnames):
            raise ValueError("Gaussian batch schema mismatch")
        for raw in reader:
            label = raw["label"]
            if label not in {row["label"] for row in expected_designs}:
                raise ValueError(f"unexpected Gaussian design {label}")
            row = {field: int(raw[field]) for field in integer_fields}
            if row["radius"] != 2:
                raise ValueError("Gaussian scorer accepts only frozen R2 rows")
            if row["h4_plus"] != row["primal_h4"] + row["matching_h4"]:
                raise ValueError("Gaussian plus identity failed")
            if row["h4_minus"] != row["primal_h4"] - row["matching_h4"]:
                raise ValueError("Gaussian minus identity failed")
            expected = next(item for item in expected_designs if item["label"] == label)
            matrix = expected["matrix"]
            if (row["n"], row["m00"], row["m01"], row["m10"], row["m11"],
                row["smith1"], row["smith2"]) != (
                    expected["N"], matrix[0][0], matrix[0][1], matrix[1][0],
                    matrix[1][1], expected["smith"][0], expected["smith"][1]):
                raise ValueError(f"Gaussian arithmetic mismatch for {label}")
            key = label, row["batch"]
            if key in rows:
                raise ValueError(f"duplicate Gaussian row {key}")
            rows[key] = row

    batches = sorted({batch for _, batch in rows})
    expected_grid = {(row["label"], batch) for row in expected_designs for batch in batches}
    if set(rows) != expected_grid:
        raise ValueError("Gaussian design/batch grid mismatch")
    if len(batches) != metadata.get("batches"):
        raise ValueError("Gaussian batch count mismatch")
    for batch in batches:
        by_size: dict[int, list[dict[str, int]]] = {}
        for design in expected_designs:
            by_size.setdefault(design["N"], []).append(rows[design["label"], batch])
        for size, pair in by_size.items():
            if len(pair) != 2:
                raise ValueError(f"N{size} does not have one orientation pair")
            for field in (
                "samples", "counter_first", "counter_last_exclusive",
                "common_field_digest",
            ):
                if pair[0][field] != pair[1][field]:
                    raise ValueError(f"N{size} common-field coupling failed for {field}")

    sum_fields = (
        "samples", "primal_pivotal", "matching_pivotal", "h4_plus", "h4_minus",
    )

    def aggregate(label: str, omitted: int | None = None) -> dict[str, int]:
        result = {field: 0 for field in sum_fields}
        for batch in batches:
            if batch == omitted:
                continue
            row = rows[label, batch]
            for field in sum_fields:
                result[field] += row[field]
        return result

    lineages = manifest["missing_Gaussian_acquisition"]["cover_lineages"]
    order: list[str] = []

    def vector(omitted: int | None = None) -> np.ndarray:
        values: list[float] = []
        if omitted is None:
            order.clear()
        for lineage_index, lineage in enumerate(lineages):
            designs = lineage["designs"]
            for channel in CHANNELS:
                column = "h4_plus" if channel == "A_plus" else "h4_minus"
                for size in lineage["sizes"]:
                    pair = [row for row in designs if row["N"] == size]
                    if len(pair) != 2:
                        raise ValueError(f"lineage {lineage_index} N{size} pair mismatch")
                    amplitudes = []
                    for design in pair:
                        total = aggregate(design["label"], omitted)
                        pivotal = total["primal_pivotal"] + total["matching_pivotal"]
                        if pivotal == 0:
                            raise ValueError(f"zero pivotal denominator for {design['label']}")
                        amplitudes.append(total[column] / pivotal)
                    values.append(amplitudes[0] - amplitudes[1])
                    if omitted is None:
                        order.append(f"Gaussian_L{lineage_index}_N{size}_{channel}")
        return np.asarray(values, dtype=float)

    point = vector()
    delete = np.stack([vector(batch) for batch in batches])
    covariance = jackknife_covariance(delete)
    audit = {
        "batches": len(batches),
        "samples_per_design": metadata["samples_per_design"],
        "pivotal_normalization": "(first H4/pivotal mass) minus (second H4/pivotal mass)",
        "common_field_digest_verified_within_each_equal_N_pair": True,
    }
    return point, covariance, order, audit


def read_annulus(analysis_path: Path) -> tuple[np.ndarray, np.ndarray, list[str]]:
    payload = json.loads(analysis_path.read_text(encoding="utf-8"))
    if payload.get("schema") != "matching-one/norm5-multiradius-pivotal-score/v1":
        raise ValueError("annulus analysis schema mismatch")
    block = payload["contrast_vector"]
    source_order = block["order"]
    selected = [f"N425_R{radius}_Delta_{channel}"
                for channel in CHANNELS for radius in RADII]
    indices = [source_order.index(label) for label in selected]
    point = np.asarray(block["point"], dtype=float)[indices]
    covariance = np.asarray(block["covariance"], dtype=float)[np.ix_(indices, indices)]
    return point, covariance, selected


def block_diagonal(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    result = np.zeros((first.shape[0] + second.shape[0],
                       first.shape[1] + second.shape[1]), dtype=float)
    result[:first.shape[0], :first.shape[1]] = first
    result[first.shape[0]:, first.shape[1]:] = second
    return result


def synthetic_point(lambda_gaussian: float, lambda_annulus: float) -> np.ndarray:
    values: list[float] = []
    for series in range(4):
        coefficients = np.asarray((0.1 * (series + 1), -0.03, 0.2), dtype=float)
        for coordinate in range(4):
            values.append(float(radial_basis(lambda_gaussian, coordinate) @ coefficients))
    radial_coordinate = np.log2(np.asarray(RADII, dtype=float) / RADII[0])
    for series in range(2):
        coefficients = np.asarray((-0.15 * (series + 1), 0.04, -0.12), dtype=float)
        for coordinate in radial_coordinate:
            values.append(float(radial_basis(lambda_annulus, coordinate) @ coefficients))
    return np.asarray(values)


def synthetic_recovery() -> dict[str, Any]:
    covariance = np.diag(np.linspace(1.0, 2.0, 24)) * 1e-8
    shared = all_fixed_scores(synthetic_point(0.5, 0.5), covariance)
    enriched = all_fixed_scores(synthetic_point(0.5, 1.0), covariance)
    if shared["best_shared_lambda"] != "1/2" or shared["best_context_pair"] != {
            "Gaussian": "1/2", "annulus": "1/2"}:
        raise AssertionError("synthetic shared-generator recovery failed")
    if enriched["best_context_pair"] != {"Gaussian": "1/2", "annulus": "1"}:
        raise AssertionError("synthetic context-enriched recovery failed")
    return {
        "shared_truth": {
            "truth": {"Gaussian": "1/2", "annulus": "1/2"},
            "recovered_shared": shared["best_shared_lambda"],
            "recovered_pair": shared["best_context_pair"],
        },
        "context_enriched_truth": {
            "truth": {"Gaussian": "1/2", "annulus": "1"},
            "recovered_shared": enriched["best_shared_lambda"],
            "recovered_pair": enriched["best_context_pair"],
            "delta": enriched["delta_min_shared_minus_min_context_pair"],
        },
    }


def score(batch_path: Path, metadata_path: Path, annulus_path: Path,
          manifest_path: Path, bootstrap_draws: int = 100000,
          require_production: bool = True) -> dict[str, Any]:
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "matching-one/p267-gaussian-annulus-missing-cells/v1":
        raise ValueError("manifest schema mismatch")
    if require_production and manifest.get("production_authorized") is not True:
        raise ValueError("manifest has not authorized production scoring")
    gaussian_point, gaussian_covariance, gaussian_order, audit = read_gaussian_batches(
        batch_path, metadata_path, manifest, require_production
    )
    annulus_point, annulus_covariance, annulus_order = read_annulus(annulus_path)
    point = np.r_[gaussian_point, annulus_point]
    covariance = block_diagonal(gaussian_covariance, annulus_covariance)
    scores = all_fixed_scores(point, covariance)
    bootstrap = bootstrap_delta(
        point, covariance, scores["delta_min_shared_minus_min_context_pair"],
        bootstrap_draws,
    )
    return {
        "schema": "matching-one/p267-gaussian-annulus-context-score/v1",
        "status": "frozen_heldout_context_rectangle_score",
        "sources": {
            "Gaussian_batches": {"path": str(batch_path), "sha256": sha256(batch_path)},
            "Gaussian_metadata": {"path": str(metadata_path), "sha256": sha256(metadata_path)},
            "annulus_analysis": {
                "path_at_commit": ANNULUS_SOURCE_PATH,
                "commit": ANNULUS_SOURCE_COMMIT,
                "sha256": sha256(annulus_path),
            },
            "manifest": {"path": str(manifest_path), "sha256": sha256(manifest_path)},
        },
        "dependency_groups": {
            "Gaussian": {
                "order": gaussian_order,
                "covariance": "complete 16x16 delete-one batch covariance",
                **audit,
            },
            "annulus": {
                "order": annulus_order,
                "covariance": "existing complete P253 N425 8x8 covariance",
            },
            "between_contexts": "zero; disjoint frozen seed/counter domains",
        },
        "base_point_order": gaussian_order + annulus_order,
        "base_point": point.tolist(),
        "base_covariance": covariance.tolist(),
        "fixed_scores": scores,
        "bootstrap": bootstrap,
        "synthetic_gate": synthetic_recovery(),
        "interpretation_boundary": {
            "exact": (
                "Each fixed score is the full-covariance GLS norm of four Gaussian "
                "and two annulus frozen recurrence residuals."
            ),
            "mechanism_inference": (
                "A calibrated positive Delta supports context-specific effective transfer "
                "generators over one shared candidate lambda."
            ),
            "not_implied": "context gain alone is not path/state memory",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--annulus-analysis", type=Path)
    parser.add_argument(
        "--manifest", type=Path,
        default=Path("analysis/p267_gaussian_annulus_missing_cells_20260829.yaml"),
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bootstrap-draws", type=int, default=100000)
    parser.add_argument("--allow-smoke", action="store_true")
    parser.add_argument("--synthetic-self-test", action="store_true")
    args = parser.parse_args()
    if args.synthetic_self_test:
        print(json.dumps(synthetic_recovery(), indent=2, sort_keys=True))
        return
    required: Iterable[tuple[str, Path | None]] = (
        ("--batches", args.batches), ("--metadata", args.metadata),
        ("--annulus-analysis", args.annulus_analysis), ("--output", args.output),
    )
    missing = [name for name, value in required if value is None]
    if missing:
        parser.error("required unless --synthetic-self-test: " + ", ".join(missing))
    assert args.batches and args.metadata and args.annulus_analysis and args.output
    payload = score(
        args.batches, args.metadata, args.annulus_analysis, args.manifest,
        args.bootstrap_draws, not args.allow_smoke,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
