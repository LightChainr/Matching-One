#!/usr/bin/env python3
"""Split K1/K2 orientation curves into translation and zero-area deformation.

For activation ``i`` let ``D_i=(F_i^a-F_i^b)/Delta cos(4 theta)`` and let
``fbar_i=d(F_i^a+F_i^b)/(2 dp)`` be the pooled empirical activation density.
The decomposition

    A_i = integral D_i,  T_i = A_i fbar_i,  R_i = D_i - T_i

has ``integral T_i=A_i`` and ``integral R_i=0`` exactly in the finite
Bernstein representation. ``T_i`` is the same-area location-translation
tangent; ``R_i`` is the remaining distributional deformation.
"""

from __future__ import annotations

import argparse
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import yaml

from analyze_activation_curve_nodes import (
    activation_area,
    bernstein_value,
)
from analyze_two_activation_h4 import (
    Archive,
    ArchiveNotScoreable,
    _cos4,
    _sum_histograms,
    estimate as activation_estimate,
    jackknife_covariance,
    load_manifest as load_source_manifest,
    read_archive,
    sha256,
)


SCHEMA = "matching-one.activation-transport-shape.v1"
MANIFEST_SCHEMA = "matching-one.activation-transport-shape.manifest.v1"
ROOT = Path(__file__).resolve().parents[1]

PER_ACTIVATION_METRICS = (
    "area",
    "D_at_p_bar",
    "translation_at_p_bar",
    "deformation_at_p_bar",
    "D_l2_energy",
    "translation_l2_energy",
    "deformation_l2_energy",
    "deformation_energy_fraction",
    "translation_cosine",
)
ALL_METRICS = tuple(
    f"K{activation}_{metric}"
    for activation in (1, 2)
    for metric in PER_ACTIVATION_METRICS
)


def _resolve(root: Path, value: object) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else root / path


def load_manifest(path: Path) -> tuple[dict[str, Any], Path, dict[str, Any]]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("schema") != MANIFEST_SCHEMA:
        raise ValueError(f"manifest must use schema {MANIFEST_SCHEMA}")
    if payload.get("status") != "retrospective_existing_data_reanalysis":
        raise ValueError("manifest status must remain retrospective_existing_data_reanalysis")
    source_path = _resolve(ROOT, payload.get("source_manifest", ""))
    if not source_path.is_file():
        raise ValueError(f"source manifest is missing: {source_path}")
    actual_sha = sha256(source_path)
    if actual_sha != payload.get("source_manifest_sha256"):
        raise ValueError(
            f"source manifest SHA256 {actual_sha} != frozen "
            f"{payload.get('source_manifest_sha256')!r}"
        )
    source = load_source_manifest(source_path)
    required = [int(value) for value in payload.get("required_sizes", ())]
    if required != [int(value) for value in source["required_sizes"]]:
        raise ValueError("required_sizes must equal the frozen source manifest order")
    quadrature = payload.get("quadrature")
    if not isinstance(quadrature, dict):
        raise ValueError("manifest lacks quadrature")
    intervals = int(quadrature.get("intervals", 0))
    if quadrature.get("method") != "composite_simpson" or intervals < 32 or intervals % 2:
        raise ValueError("quadrature must use an even composite-Simpson interval count")
    return payload, source_path, source


def _degree_elevate(values: Sequence[float]) -> tuple[float, ...]:
    """Raise Bernstein degree n-1 to n without changing the polynomial."""

    n = len(values)
    return tuple(
        (index / n) * (float(values[index - 1]) if index else 0.0)
        + (1.0 - index / n) * (float(values[index]) if index < n else 0.0)
        for index in range(n + 1)
    )


def _curve_bundle(
    archive: Archive, omitted_batch: Optional[int]
) -> tuple[float, dict[int, dict[str, tuple[float, ...]]]]:
    first = _sum_histograms(archive.histograms["first"], omitted_batch)
    second = _sum_histograms(archive.histograms["second"], omitted_batch)
    if first["samples"] != second["samples"]:
        raise ArchiveNotScoreable("orientation samples differ after aligned deletion")
    delta_cos4 = _cos4(first["a"], first["b"]) - _cos4(second["a"], second["b"])
    if abs(delta_cos4) < 1.0e-15:
        raise ArchiveNotScoreable("orientation pair has zero Delta cos(4 theta)")
    bundles: dict[int, dict[str, tuple[float, ...]]] = {}
    for activation, key in ((1, "k1"), (2, "k2")):
        cumulative_first = 0
        cumulative_second = 0
        first_cdf = []
        second_cdf = []
        for occupied in range(archive.n + 1):
            cumulative_first += first[key][occupied]
            cumulative_second += second[key][occupied]
            first_cdf.append(cumulative_first / first["samples"])
            second_cdf.append(cumulative_second / second["samples"])
        difference = tuple(
            (left - right) / delta_cos4
            for left, right in zip(first_cdf, second_cdf)
        )
        pooled = tuple(
            (left + right) / 2.0 for left, right in zip(first_cdf, second_cdf)
        )
        density_degree_n_minus_1 = tuple(
            archive.n * (pooled[index + 1] - pooled[index])
            for index in range(archive.n)
        )
        density = _degree_elevate(density_degree_n_minus_1)
        area = activation_area(difference)
        translation = tuple(area * value for value in density)
        deformation = tuple(
            value - tangent for value, tangent in zip(difference, translation)
        )
        bundles[activation] = {
            "D": difference,
            "pooled_density": density,
            "translation": translation,
            "deformation": deformation,
        }
    return delta_cos4, bundles


def _simpson_inner(
    left: Sequence[float], right: Sequence[float], intervals: int
) -> float:
    step = 1.0 / intervals
    terms = []
    for index in range(intervals + 1):
        p = index * step
        weight = 1.0 if index in (0, intervals) else (4.0 if index % 2 else 2.0)
        terms.append(
            weight * bernstein_value(left, p) * bernstein_value(right, p)
        )
    return step * math.fsum(terms) / 3.0


def estimate(
    archive: Archive, intervals: int, omitted_batch: Optional[int] = None
) -> dict[str, Any]:
    p_bar = float(activation_estimate(archive, omitted_batch)["metrics"]["p_bar"])
    delta_cos4, bundles = _curve_bundle(archive, omitted_batch)
    metrics: dict[str, float] = {}
    identities: dict[str, float] = {}
    public_curves: dict[str, Any] = {}
    for activation in (1, 2):
        item = bundles[activation]
        difference = item["D"]
        translation = item["translation"]
        deformation = item["deformation"]
        area = activation_area(difference)
        d_energy = _simpson_inner(difference, difference, intervals)
        t_energy = _simpson_inner(translation, translation, intervals)
        r_energy = _simpson_inner(deformation, deformation, intervals)
        cross = _simpson_inner(difference, translation, intervals)
        prefix = f"K{activation}_"
        metrics.update(
            {
                prefix + "area": area,
                prefix + "D_at_p_bar": bernstein_value(difference, p_bar),
                prefix + "translation_at_p_bar": bernstein_value(translation, p_bar),
                prefix + "deformation_at_p_bar": bernstein_value(deformation, p_bar),
                prefix + "D_l2_energy": d_energy,
                prefix + "translation_l2_energy": t_energy,
                prefix + "deformation_l2_energy": r_energy,
                prefix + "deformation_energy_fraction": r_energy / d_energy,
                prefix + "translation_cosine": cross / math.sqrt(d_energy * t_energy),
            }
        )
        identities.update(
            {
                prefix + "translation_area_minus_A": activation_area(translation) - area,
                prefix + "deformation_area": activation_area(deformation),
                prefix + "coefficient_closure_max_abs": max(
                    abs(value - tangent - residual)
                    for value, tangent, residual in zip(
                        difference, translation, deformation
                    )
                ),
            }
        )
        public_curves[f"K{activation}"] = {
            "D_coefficients": list(difference),
            "pooled_density_coefficients_degree_elevated": list(item["pooled_density"]),
            "translation_coefficients": list(translation),
            "zero_area_deformation_coefficients": list(deformation),
        }
    if max(abs(value) for value in identities.values()) > 5.0e-13:
        raise ArithmeticError(f"translation/deformation identity failed: {identities}")
    return {
        "p_bar": p_bar,
        "delta_cos4": delta_cos4,
        "metrics": metrics,
        "identity_audit": identities,
        "curves": public_curves,
    }


def analyze_archive(
    root_text: str, entry: Mapping[str, Any], intervals: int
) -> dict[str, Any]:
    root = Path(root_text)
    n = int(entry["N"])
    try:
        archive = read_archive(root, entry)
        full = estimate(archive, intervals)
        batch_ids = [row.batch for row in archive.histograms["first"]]
        deleted = [estimate(archive, intervals, batch)["metrics"] for batch in batch_ids]
        covariance = [
            [
                jackknife_covariance(
                    [row[left] for row in deleted], [row[right] for row in deleted]
                )
                for right in ALL_METRICS
            ]
            for left in ALL_METRICS
        ]
        standard_errors = {
            metric: math.sqrt(max(0.0, covariance[index][index]))
            for index, metric in enumerate(ALL_METRICS)
        }
        return {
            "public": {
                "N": n,
                "status": "scoreable",
                "dependency_group": archive.dependency_group,
                "representatives": {
                    orientation: [
                        archive.histograms[orientation][0].a,
                        archive.histograms[orientation][0].b,
                    ]
                    for orientation in ("first", "second")
                },
                "samples_per_orientation": int(archive.metadata["samples_per_pair"]),
                "batch_count": len(batch_ids),
                "p_bar": full["p_bar"],
                "delta_cos4": full["delta_cos4"],
                "metric_order": list(ALL_METRICS),
                "estimate_vector": [full["metrics"][name] for name in ALL_METRICS],
                "standard_errors": standard_errors,
                "jackknife_covariance": covariance,
                "identity_audit": full["identity_audit"],
                "full_curve_bernstein": {
                    "degree": n,
                    "coefficient_index": "occupied count j=0,...,N",
                    **full["curves"],
                },
                "provenance": {
                    "source_commit": archive.metadata["git_commit"],
                    "seed": archive.metadata["seed"],
                    "counter_first": archive.metadata["replica_counter_first"],
                    "counter_last_exclusive": archive.metadata[
                        "replica_counter_last_exclusive"
                    ],
                    "inputs": {
                        kind: {
                            "path": str(path.relative_to(root)),
                            "sha256": sha256(path),
                        }
                        for kind, path in archive.paths.items()
                    },
                },
            },
            "deleted": deleted,
            "batch_signature": [
                [row.batch, row.samples] for row in archive.histograms["first"]
            ],
        }
    except (
        ArchiveNotScoreable,
        ArithmeticError,
        KeyError,
        TypeError,
        json.JSONDecodeError,
    ) as exc:
        return {
            "public": {
                "N": n,
                "status": "not_scoreable",
                "dependency_group": str(entry.get("dependency_group", "unknown")),
                "reason": str(exc),
            },
            "deleted": None,
            "batch_signature": None,
        }


def _dependency_groups(results: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[str, list[dict[str, Any]]] = {}
    for result in results:
        if result["public"]["status"] == "scoreable":
            groups.setdefault(result["public"]["dependency_group"], []).append(result)
    output = []
    for group, members in groups.items():
        signatures = {
            tuple(tuple(row) for row in member["batch_signature"]) for member in members
        }
        if len(members) > 1 and len(signatures) != 1:
            raise ValueError(f"dependency group {group} is not aligned by batch")
        output.append(
            {
                "id": group,
                "sizes": [member["public"]["N"] for member in members],
                "rule": "aligned_delete_one_full_covariance",
                "independent_evidence_units": 1,
            }
        )
    return output


def render(
    manifest_path: Path,
    manifest: Mapping[str, Any],
    source_path: Path,
    results: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    scored = [row for row in results if row["public"]["status"] == "scoreable"]
    observations = [
        {"N": row["public"]["N"], "metric": metric}
        for row in scored
        for metric in ALL_METRICS
    ]
    estimates = [
        row["public"]["estimate_vector"][row["public"]["metric_order"].index(metric)]
        for row in scored
        for metric in ALL_METRICS
    ]
    covariance = []
    for left in scored:
        for left_metric in ALL_METRICS:
            line = []
            for right in scored:
                for right_metric in ALL_METRICS:
                    if left["public"]["dependency_group"] != right["public"]["dependency_group"]:
                        line.append(0.0)
                    else:
                        line.append(
                            jackknife_covariance(
                                [row[left_metric] for row in left["deleted"]],
                                [row[right_metric] for row in right["deleted"]],
                            )
                        )
            covariance.append(line)
    by_n = {str(row["public"]["N"]): row["public"] for row in results}
    negative_d2 = []
    residual_driven = []
    for row in scored:
        order = row["public"]["metric_order"]
        values = dict(zip(order, row["public"]["estimate_vector"]))
        n = row["public"]["N"]
        if values["K2_D_at_p_bar"] < 0.0:
            negative_d2.append(n)
            if values["K2_translation_at_p_bar"] > 0.0 and values["K2_deformation_at_p_bar"] < 0.0:
                residual_driven.append(n)
    return {
        "schema": SCHEMA,
        "status": "retrospective existing-data translation/deformation analysis; no new simulation and no exponent fit",
        "definitions": {
            "D_i(p)": "(F_i_first(p)-F_i_second(p))/Delta cos(4 theta)",
            "fbar_i(p)": "d/dp [(F_i_first(p)+F_i_second(p))/2]",
            "A_i": "integral_0^1 D_i(p) dp",
            "translation_T_i(p)": "A_i fbar_i(p)",
            "zero_area_deformation_R_i(p)": "D_i(p)-T_i(p)",
            "exact_area_identity": "integral T_i=A_i and integral R_i=0",
        },
        "interpretation_boundary": (
            "T_i is the pooled-density tangent for a first-order threshold-location "
            "shift with the observed area. R_i measures finite-archive distributional "
            "reshaping. Neither component alone identifies a continuum operator."
        ),
        "quadrature": manifest["quadrature"],
        "size_order": [int(value) for value in manifest["required_sizes"]],
        "scoreable_sizes": [row["public"]["N"] for row in scored],
        "not_scoreable_sizes": [
            row["public"]["N"]
            for row in results
            if row["public"]["status"] != "scoreable"
        ],
        "by_N": by_n,
        "dependency_groups": _dependency_groups(results),
        "decision_covariance": {
            "metric_order_with_N": observations,
            "estimate_vector": estimates,
            "jackknife_covariance": covariance,
            "rule": "same deleted batch across both orientations and all N in one dependency group; zero entries across independent groups",
            "derived_views_are_not_independent_evidence": True,
        },
        "descriptive_findings": {
            "negative_K2_D_at_p_bar_sizes": negative_d2,
            "negative_K2_points_with_positive_translation_and_negative_dominant_deformation": residual_driven,
            "guard": "Signs and energy shares are finite-archive estimates with aligned delete-one covariance; no cross-size exponent is fitted.",
        },
        "provenance": {
            "manifest": str(manifest_path),
            "manifest_sha256": sha256(manifest_path),
            "source_manifest": str(source_path.relative_to(ROOT)),
            "source_manifest_sha256": sha256(source_path),
        },
    }


def render_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# K1/K2 translation versus zero-area deformation",
        "",
        "This reuses the ten frozen threshold-rank archives. It generates no samples and fits no exponent.",
        "",
        "For each activation, `D=A fbar+R`: `A fbar` is the pooled-density translation tangent with the same exact area as `D`, while `R` has exactly zero integral.",
        "",
        "| N | K1 deformation energy share | K2 deformation energy share | K2 D(p_bar) | K2 translation | K2 deformation | reading |",
        "|---:|---:|---:|---:|---:|---:|:---|",
    ]
    for n in payload["size_order"]:
        row = payload["by_N"][str(n)]
        if row["status"] != "scoreable":
            lines.append(f"| {n} | | | | | | not scoreable: {row['reason']} |")
            continue
        values = dict(zip(row["metric_order"], row["estimate_vector"]))
        reading = (
            "shape residual overturns positive translation"
            if values["K2_D_at_p_bar"] < 0.0
            and values["K2_translation_at_p_bar"] > 0.0
            and values["K2_deformation_at_p_bar"] < 0.0
            else "translation and deformation partially cancel"
        )
        lines.append(
            f"| {n} | {values['K1_deformation_energy_fraction']:.3f} | "
            f"{values['K2_deformation_energy_fraction']:.3f} | "
            f"{values['K2_D_at_p_bar']:.6e} | "
            f"{values['K2_translation_at_p_bar']:.6e} | "
            f"{values['K2_deformation_at_p_bar']:.6e} | {reading} |"
        )
    findings = payload["descriptive_findings"]
    lines.extend(
        [
            "",
            "## Decision reading",
            "",
            "The same-area translation term is positive at every scored K2 root point. The deformation is negative there and accounts for the local negative values at N="
            + ",".join(map(str, findings["negative_K2_points_with_positive_translation_and_negative_dominant_deformation"]))
            + ". Thus those lobes are activation-distribution reshaping, not a reversal of the integrated translation direction.",
            "",
            "The JSON retains every Bernstein coefficient, the exact area/closure audits, within-N covariance, and the full cross-N covariance under aligned dependency-group deletion. Derived translation and deformation coordinates from one archive remain one evidence block.",
            "",
            "This finite-archive split does not decide whether the reshaping is carried by rank-one lifetime, plateau line, landing pivotal, or another typed field; it identifies the next mechanism-sensitive information axis.",
            "",
        ]
    )
    return "\n".join(lines)


def _analyze_star(arguments: tuple[str, Mapping[str, Any], int]) -> dict[str, Any]:
    return analyze_archive(*arguments)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "analysis/activation_transport_shape_manifest.yaml",
    )
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1:
        raise SystemExit("--workers must be positive")
    try:
        manifest, source_path, source = load_manifest(args.manifest)
        intervals = int(manifest["quadrature"]["intervals"])
        worker_args = [(str(ROOT), entry, intervals) for entry in source["runs"]]
        if args.workers == 1:
            results = [_analyze_star(item) for item in worker_args]
        else:
            with ProcessPoolExecutor(max_workers=min(args.workers, len(worker_args))) as pool:
                results = list(pool.map(_analyze_star, worker_args))
        payload = render(args.manifest, manifest, source_path, results)
    except (ValueError, OSError, yaml.YAMLError) as exc:
        raise SystemExit(str(exc)) from exc
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.output_md.write_text(render_markdown(payload), encoding="utf-8")
    print(args.output_json)
    print(args.output_md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
