#!/usr/bin/env python3
"""Resolve the full K1/K2 orientation curves and their stable node branches.

This is a zero-production companion to ``analyze_two_activation_h4.py``.  It
uses the same ten frozen histogram archives and the same aligned delete-one
batches.  The complete empirical curves are serialized in the Bernstein basis,
so their integrals are exact finite sums and their values remain numerically
stable away from the matching root.

For orientation ``first - second`` and the exact angular contrast ``Delta c4``:

    Di(p) = (Fi_first(p) - Fi_second(p)) / Delta c4,
    Ai    = integral_0^1 Di(p) dp
          = -Delta E[Ki] / ((N + 1) Delta c4).

No Monte Carlo samples are generated and no scaling exponent is fitted.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import yaml

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


SCHEMA = "matching-one.activation-curve-nodes.v1"
MANIFEST_SCHEMA = "matching-one.activation-curve-nodes.manifest.v1"
ESTIMATE_METRICS = (
    "A1",
    "A2",
    "D1_at_p_bar",
    "D2_at_p_bar",
    "D1_prime_at_p_bar",
    "D2_prime_at_p_bar",
)
DECISION_METRICS = ("A2", "D2_at_p_bar", "D2_prime_at_p_bar")


def _resolve(root: Path, value: object) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else root / path


def load_manifest(path: Path, root: Path) -> tuple[dict[str, Any], Path, dict[str, Any]]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("schema") != MANIFEST_SCHEMA:
        raise ValueError(f"manifest must use schema {MANIFEST_SCHEMA}")
    if payload.get("status") != "retrospective_existing_data_reanalysis":
        raise ValueError("manifest status must remain retrospective_existing_data_reanalysis")
    source_path = _resolve(root, payload.get("source_manifest", ""))
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
    policy = payload.get("node_policy")
    if not isinstance(policy, dict):
        raise ValueError("manifest lacks node_policy")
    full_domain = [float(value) for value in policy.get("full_domain", ())]
    critical = [float(value) for value in policy.get("critical_window", ())]
    if (
        len(full_domain) != 2
        or len(critical) != 2
        or not 0.0 < full_domain[0] < critical[0] < critical[1] < full_domain[1] < 1.0
    ):
        raise ValueError("node domains must be nested strictly inside (0,1)")
    for name in ("point_grid_intervals", "jackknife_grid_intervals"):
        value = int(policy.get(name, 0))
        if value < 64:
            raise ValueError(f"{name} must be at least 64")
    if policy.get("node_kind") != "simple_sign_changing":
        raise ValueError("only simple_sign_changing nodes are supported")
    return payload, source_path, source


def bernstein_value(coefficients: Sequence[float], p: float) -> float:
    """Evaluate a Bernstein polynomial using a mode-centred binomial recurrence."""

    degree = len(coefficients) - 1
    if degree < 0:
        raise ValueError("a Bernstein coefficient vector cannot be empty")
    if p <= 0.0:
        return float(coefficients[0])
    if p >= 1.0:
        return float(coefficients[-1])
    mode = min(degree, max(0, int((degree + 1) * p)))
    log_probability = (
        math.lgamma(degree + 1)
        - math.lgamma(mode + 1)
        - math.lgamma(degree - mode + 1)
        + mode * math.log(p)
        + (degree - mode) * math.log1p(-p)
    )
    probability_at_mode = math.exp(log_probability)
    terms = [float(coefficients[mode]) * probability_at_mode]
    probability = probability_at_mode
    odds = p / (1.0 - p)
    for occupied in range(mode, degree):
        probability *= (degree - occupied) * odds / (occupied + 1)
        terms.append(float(coefficients[occupied + 1]) * probability)
    probability = probability_at_mode
    inverse_odds = (1.0 - p) / p
    for occupied in range(mode, 0, -1):
        probability *= occupied * inverse_odds / (degree - occupied + 1)
        terms.append(float(coefficients[occupied - 1]) * probability)
    return math.fsum(terms)


def bernstein_derivative(coefficients: Sequence[float], p: float) -> float:
    degree = len(coefficients) - 1
    if degree <= 0:
        return 0.0
    differences = [
        degree * (float(right) - float(left))
        for left, right in zip(coefficients, coefficients[1:])
    ]
    return bernstein_value(differences, p)


def activation_area(coefficients: Sequence[float]) -> float:
    """Integrate a Bernstein polynomial exactly from its coefficients."""

    return math.fsum(float(value) for value in coefficients) / len(coefficients)


def _mean_rank(counts: Sequence[int], samples: int) -> float:
    return math.fsum(rank * count for rank, count in enumerate(counts)) / samples


def _curve_coefficients(
    archive: Archive, omitted_batch: Optional[int]
) -> tuple[dict[str, Any], dict[str, Any], float, dict[str, tuple[float, ...]]]:
    first = _sum_histograms(archive.histograms["first"], omitted_batch)
    second = _sum_histograms(archive.histograms["second"], omitted_batch)
    if first["samples"] != second["samples"]:
        raise ArchiveNotScoreable("orientation samples differ after aligned deletion")
    delta_cos4 = _cos4(first["a"], first["b"]) - _cos4(second["a"], second["b"])
    if abs(delta_cos4) < 1.0e-15:
        raise ArchiveNotScoreable("orientation pair has zero Delta cos(4 theta)")
    output: dict[str, tuple[float, ...]] = {}
    for label, key in (("D1", "k1"), ("D2", "k2")):
        cumulative_first = 0
        cumulative_second = 0
        coefficients = []
        for occupied in range(archive.n + 1):
            cumulative_first += first[key][occupied]
            cumulative_second += second[key][occupied]
            coefficients.append(
                (
                    cumulative_first / first["samples"]
                    - cumulative_second / second["samples"]
                )
                / delta_cos4
            )
        output[label] = tuple(coefficients)
    return first, second, delta_cos4, output


def _center_width(
    archive: Archive,
    first: Mapping[str, Any],
    second: Mapping[str, Any],
    delta_cos4: float,
) -> dict[str, Any]:
    denominator = archive.n + 1.0

    def coordinates(row: Mapping[str, Any]) -> dict[str, float]:
        k1 = _mean_rank(row["k1"], int(row["samples"]))
        k2 = _mean_rank(row["k2"], int(row["samples"]))
        center = (k1 + k2) / (2.0 * denominator)
        width = (k2 - k1) / denominator
        return {
            "E_K1": k1,
            "E_K2": k2,
            "E_C": center,
            "E_W": width,
            "E_K2_over_N_plus_1": k2 / denominator,
            "identity_residual": k2 / denominator - center - width / 2.0,
        }

    left = coordinates(first)
    right = coordinates(second)
    delta_center = left["E_C"] - right["E_C"]
    delta_width = left["E_W"] - right["E_W"]
    center_contribution = -delta_center / delta_cos4
    width_contribution = -delta_width / (2.0 * delta_cos4)
    return {
        "definitions": {
            "C": "(K1+K2)/(2(N+1))",
            "W": "(K2-K1)/(N+1)",
            "identity": "K2/(N+1)=C+W/2",
        },
        "first": left,
        "second": right,
        "delta_first_minus_second": {
            "delta_E_C": delta_center,
            "delta_E_W": delta_width,
        },
        "A2_decomposition": {
            "from_center_C": center_contribution,
            "from_width_W_over_2": width_contribution,
            "sum": center_contribution + width_contribution,
        },
    }


def curve_estimate(archive: Archive, omitted_batch: Optional[int] = None) -> dict[str, Any]:
    root_estimate = activation_estimate(archive, omitted_batch)
    p_bar = float(root_estimate["metrics"]["p_bar"])
    first, second, delta_cos4, coefficients = _curve_coefficients(
        archive, omitted_batch
    )
    center_width = _center_width(archive, first, second, delta_cos4)
    mean_deltas = {
        label: _mean_rank(first[key], first["samples"])
        - _mean_rank(second[key], second["samples"])
        for label, key in (("A1", "k1"), ("A2", "k2"))
    }
    metrics = {
        "A1": activation_area(coefficients["D1"]),
        "A2": activation_area(coefficients["D2"]),
        "D1_at_p_bar": bernstein_value(coefficients["D1"], p_bar),
        "D2_at_p_bar": bernstein_value(coefficients["D2"], p_bar),
        "D1_prime_at_p_bar": bernstein_derivative(coefficients["D1"], p_bar),
        "D2_prime_at_p_bar": bernstein_derivative(coefficients["D2"], p_bar),
    }
    expected_areas = {
        name: -mean_deltas[name] / ((archive.n + 1.0) * delta_cos4)
        for name in ("A1", "A2")
    }
    identity_audit = {
        "A1_minus_rank_identity": metrics["A1"] - expected_areas["A1"],
        "A2_minus_rank_identity": metrics["A2"] - expected_areas["A2"],
        "A2_minus_C_W_decomposition": metrics["A2"]
        - center_width["A2_decomposition"]["sum"],
        "first_K2_over_N_plus_1_minus_C_plus_W_over_2": center_width["first"][
            "identity_residual"
        ],
        "second_K2_over_N_plus_1_minus_C_plus_W_over_2": center_width["second"][
            "identity_residual"
        ],
        "D1_at_p_bar_minus_two_activation_value": metrics["D1_at_p_bar"]
        - root_estimate["metrics"]["angular_delta_F1"],
        "D2_at_p_bar_minus_two_activation_value": metrics["D2_at_p_bar"]
        - root_estimate["metrics"]["angular_delta_F2"],
    }
    if max(abs(value) for value in identity_audit.values()) > 3.0e-13:
        raise ArithmeticError("activation-curve algebraic identity failed")
    return {
        "p_bar": p_bar,
        "delta_cos4": delta_cos4,
        "metrics": metrics,
        "coefficients": coefficients,
        "center_width": center_width,
        "rank_area_formula": expected_areas,
        "identity_audit": identity_audit,
    }


def _bisect_root(
    coefficients: Sequence[float], lower: float, upper: float, tolerance: float
) -> float:
    left_value = bernstein_value(coefficients, lower)
    right_value = bernstein_value(coefficients, upper)
    if left_value == 0.0:
        return lower
    if right_value == 0.0:
        return upper
    if left_value * right_value > 0.0:
        raise ArithmeticError("root bracket does not change sign")
    for _ in range(80):
        middle = (lower + upper) / 2.0
        middle_value = bernstein_value(coefficients, middle)
        if middle_value == 0.0 or upper - lower <= tolerance:
            return middle
        if left_value * middle_value < 0.0:
            upper = middle
            right_value = middle_value
        else:
            lower = middle
            left_value = middle_value
    return (lower + upper) / 2.0


def find_sign_changing_nodes(
    coefficients: Sequence[float],
    domain: Sequence[float],
    grid_intervals: int,
    tolerance: float,
) -> list[float]:
    """Find simple sign-changing nodes; the endpoint zeros are out of scope."""

    lower, upper = (float(domain[0]), float(domain[1]))
    step = (upper - lower) / grid_intervals
    roots: list[float] = []
    left = lower
    left_value = bernstein_value(coefficients, left)
    for index in range(1, grid_intervals + 1):
        right = lower + index * step
        right_value = bernstein_value(coefficients, right)
        if left_value == 0.0:
            roots.append(left)
        elif right_value == 0.0 or left_value * right_value < 0.0:
            roots.append(_bisect_root(coefficients, left, right, tolerance))
        left = right
        left_value = right_value
    unique = []
    for root in roots:
        if not unique or abs(root - unique[-1]) > 4.0 * tolerance:
            unique.append(root)
    return unique


def _node_spectrum(
    point_nodes: Sequence[float],
    deleted_nodes: Sequence[Sequence[float]],
    domain: Sequence[float],
) -> dict[str, Any]:
    point = [value for value in point_nodes if domain[0] <= value <= domain[1]]
    deleted = [
        [value for value in roots if domain[0] <= value <= domain[1]]
        for roots in deleted_nodes
    ]
    branch_counts = Counter(len(roots) for roots in deleted)
    stable = all(len(roots) == len(point) for roots in deleted)
    payload: dict[str, Any] = {
        "status": "scoreable" if stable else "not_scoreable",
        "domain": [float(domain[0]), float(domain[1])],
        "point_estimate_nodes": point,
        "point_estimate_branch_count": len(point),
        "delete_one_branch_count_histogram": {
            str(count): frequency for count, frequency in sorted(branch_counts.items())
        },
    }
    if not stable:
        payload["reason"] = (
            "delete-one branch count changes; node positions are not aligned or scored"
        )
        payload["node_standard_errors"] = None
        payload["jackknife_node_covariance"] = None
        return payload
    covariance = [
        [
            jackknife_covariance(
                [roots[left] for roots in deleted],
                [roots[right] for roots in deleted],
            )
            for right in range(len(point))
        ]
        for left in range(len(point))
    ]
    payload["node_standard_errors"] = [
        math.sqrt(max(0.0, covariance[index][index]))
        for index in range(len(point))
    ]
    payload["jackknife_node_covariance"] = covariance
    return payload


def analyze_archive(
    root_text: str, entry: Mapping[str, Any], policy: Mapping[str, Any]
) -> dict[str, Any]:
    root = Path(root_text)
    n = int(entry["N"])
    try:
        archive = read_archive(root, entry)
        full = curve_estimate(archive)
        full_nodes = {
            label: find_sign_changing_nodes(
                full["coefficients"][label],
                policy["full_domain"],
                int(policy["point_grid_intervals"]),
                float(policy["root_tolerance"]),
            )
            for label in ("D1", "D2")
        }
        batch_ids = [row.batch for row in archive.histograms["first"]]
        deleted = []
        deleted_nodes: dict[str, list[list[float]]] = {"D1": [], "D2": []}
        for batch in batch_ids:
            item = curve_estimate(archive, batch)
            deleted.append(item["metrics"])
            for label in ("D1", "D2"):
                deleted_nodes[label].append(
                    find_sign_changing_nodes(
                        item["coefficients"][label],
                        policy["full_domain"],
                        int(policy["jackknife_grid_intervals"]),
                        float(policy["root_tolerance"]),
                    )
                )
        covariance = [
            [
                jackknife_covariance(
                    [item[left] for item in deleted],
                    [item[right] for item in deleted],
                )
                for right in ESTIMATE_METRICS
            ]
            for left in ESTIMATE_METRICS
        ]
        standard_errors = {
            metric: math.sqrt(max(0.0, covariance[index][index]))
            for index, metric in enumerate(ESTIMATE_METRICS)
        }
        spectra = {}
        for label in ("D1", "D2"):
            spectra[label] = {
                "full_domain": _node_spectrum(
                    full_nodes[label], deleted_nodes[label], policy["full_domain"]
                ),
                "critical_window": _node_spectrum(
                    full_nodes[label], deleted_nodes[label], policy["critical_window"]
                ),
            }
        critical_k2 = spectra["D2"]["critical_window"]
        nearest = None
        if critical_k2["status"] == "scoreable" and critical_k2["point_estimate_nodes"]:
            nearest = min(
                critical_k2["point_estimate_nodes"],
                key=lambda value: abs(value - full["p_bar"]),
            )
        local_context = {
            "nearest_scoreable_K2_node": nearest,
            "p_bar_minus_node": full["p_bar"] - nearest if nearest is not None else None,
            "positive_integrated_A2": full["metrics"]["A2"] > 0.0,
            "negative_D2_at_p_bar": full["metrics"]["D2_at_p_bar"] < 0.0,
            "nearby_node_above_p_bar": (
                nearest is not None
                and nearest > full["p_bar"]
                and nearest - full["p_bar"] <= float(policy["nearby_node_radius"])
            ),
        }
        local_context["classification"] = (
            "local_negative_lobe_before_nearby_zero_not_integrated_reversal"
            if local_context["positive_integrated_A2"]
            and local_context["negative_D2_at_p_bar"]
            and local_context["nearby_node_above_p_bar"]
            else "no_local_negative_with_resolved_nearby_upper_node"
        )
        public = {
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
            "curve_values_at_p_bar": {
                "D1": full["metrics"]["D1_at_p_bar"],
                "D2": full["metrics"]["D2_at_p_bar"],
                "D1_prime": full["metrics"]["D1_prime_at_p_bar"],
                "D2_prime": full["metrics"]["D2_prime_at_p_bar"],
            },
            "integrated_areas": {
                "A1": full["metrics"]["A1"],
                "A2": full["metrics"]["A2"],
                "rank_formula": full["rank_area_formula"],
            },
            "center_width_identity": full["center_width"],
            "full_curve_bernstein": {
                "degree": n,
                "coefficient_index": "occupied edge count j=0,...,N",
                "D1_coefficients": list(full["coefficients"]["D1"]),
                "D2_coefficients": list(full["coefficients"]["D2"]),
            },
            "node_spectra": spectra,
            "local_K2_context": local_context,
            "estimate_vector_order": list(ESTIMATE_METRICS),
            "estimate_vector": [full["metrics"][name] for name in ESTIMATE_METRICS],
            "standard_errors": standard_errors,
            "jackknife_covariance": covariance,
            "identity_audit": full["identity_audit"],
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
        }
        return {
            "public": public,
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
                "rule": (
                    "aligned_delete_one_full_covariance"
                    if len(members) > 1
                    else "independent_archive"
                ),
                "independent_evidence_units": 1,
            }
        )
    return output


def render(
    manifest_path: Path,
    manifest: Mapping[str, Any],
    source_manifest_path: Path,
    results: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    scored = [row for row in results if row["public"]["status"] == "scoreable"]
    observations = [
        {"N": row["public"]["N"], "metric": metric}
        for row in scored
        for metric in DECISION_METRICS
    ]
    estimates = [
        row["public"]["estimate_vector"][
            row["public"]["estimate_vector_order"].index(metric)
        ]
        for row in scored
        for metric in DECISION_METRICS
    ]
    covariance = []
    for left in scored:
        for left_metric in DECISION_METRICS:
            line = []
            for right in scored:
                for right_metric in DECISION_METRICS:
                    if left["public"]["dependency_group"] != right["public"][
                        "dependency_group"
                    ]:
                        line.append(0.0)
                    else:
                        line.append(
                            jackknife_covariance(
                                [item[left_metric] for item in left["deleted"]],
                                [item[right_metric] for item in right["deleted"]],
                            )
                        )
            covariance.append(line)
    by_n = {str(row["public"]["N"]): row["public"] for row in results}
    negative_sizes = [
        row["public"]["N"]
        for row in scored
        if row["public"]["curve_values_at_p_bar"]["D2"] < 0.0
    ]
    local_explanations = [
        row["public"]["N"]
        for row in scored
        if row["public"]["local_K2_context"]["classification"]
        == "local_negative_lobe_before_nearby_zero_not_integrated_reversal"
    ]
    signs = {
        1 if row["public"]["integrated_areas"]["A2"] > 0.0 else -1
        if row["public"]["integrated_areas"]["A2"] < 0.0
        else 0
        for row in scored
    }
    return {
        "schema": SCHEMA,
        "status": "retrospective existing-data full-curve analysis; no new simulation and no exponent fit",
        "definitions": {
            "K1": "K_minus",
            "K2": "K_plus",
            "D1(p)": "(F1_first(p)-F1_second(p))/Delta cos(4 theta)",
            "D2(p)": "(F2_first(p)-F2_second(p))/Delta cos(4 theta)",
            "A_i": "integral_0^1 D_i(p) dp = -Delta E[Ki]/((N+1) Delta cos(4 theta))",
            "C": "(K1+K2)/(2(N+1))",
            "W": "(K2-K1)/(N+1)",
            "K2_identity": "K2/(N+1)=C+W/2",
        },
        "node_policy": manifest["node_policy"],
        "size_order": [int(value) for value in manifest["required_sizes"]],
        "scoreable_sizes": [row["public"]["N"] for row in scored],
        "not_scoreable_sizes": [
            row["public"]["N"]
            for row in results
            if row["public"]["status"] == "not_scoreable"
        ],
        "by_N": by_n,
        "dependency_groups": _dependency_groups(results),
        "decision_covariance": {
            "metric_order_with_N": observations,
            "estimate_vector": estimates,
            "jackknife_covariance": covariance,
            "shared_stream_rule": (
                "N=65,85,130,170 use the same aligned deleted batch and form one "
                "dependency block; distinct counter intervals have zero cross-block entries"
            ),
            "views_must_not_be_added_as_independent_evidence": True,
        },
        "descriptive_findings": {
            "A2_has_one_nonzero_sign_across_all_scoreable_sizes": len(signs) == 1
            and 0 not in signs,
            "A2_sign": "positive" if signs == {1} else "mixed_or_zero",
            "negative_D2_at_p_bar_sizes": negative_sizes,
            "negative_points_explained_by_scoreable_nearby_upper_node": local_explanations,
            "upper_node_statement_uses_point_estimate_not_significant_ordering": True,
            "guard": (
                "The area sign and finite-archive node topology are descriptive. "
                "The node-above-p_bar statement is about the point estimate, not a "
                "significant ordering after node-position uncertainty. No scaling "
                "exponent or continuum operator is identified."
            ),
        },
        "provenance": {
            "manifest": str(manifest_path),
            "manifest_sha256": sha256(manifest_path),
            "source_manifest": str(
                source_manifest_path.relative_to(Path(__file__).resolve().parents[1])
            ),
            "source_manifest_sha256": sha256(source_manifest_path),
        },
    }


def render_markdown(payload: Mapping[str, Any]) -> str:
    lines = [
        "# K1/K2 activation-curve node spectrum",
        "",
        "This is a retrospective full-curve analysis of existing threshold-rank "
        "histograms. It generates no Monte Carlo samples and fits no exponent.",
        "",
        "The exact orientation-normalized curves are "
        "`Di(p)=(Fi_first-Fi_second)/Delta cos(4 theta)`. Their areas obey "
        "`Ai=-Delta E[Ki]/((N+1) Delta cos(4 theta))`.",
        "",
        "| N | A1 | A2 | D2(p_bar) | D2'(p_bar) | nearest stable K2 node | K2 critical branches | local reading |",
        "|---:|---:|---:|---:|---:|---:|:---|:---|",
    ]
    for n in payload["size_order"]:
        row = payload["by_N"][str(n)]
        if row["status"] != "scoreable":
            lines.append(f"| {n} | | | | | | not scoreable | {row['reason']} |")
            continue
        areas = row["integrated_areas"]
        values = row["curve_values_at_p_bar"]
        spectrum = row["node_spectra"]["D2"]["critical_window"]
        nearest = row["local_K2_context"]["nearest_scoreable_K2_node"]
        nearest_text = f"{nearest:.9f}" if nearest is not None else "none"
        lines.append(
            f"| {n} | {areas['A1']:.6e} | {areas['A2']:.6e} | "
            f"{values['D2']:.6e} | {values['D2_prime']:.6e} | {nearest_text} | "
            f"{spectrum['status']} ({spectrum['point_estimate_branch_count']}) | "
            f"{row['local_K2_context']['classification']} |"
        )
    findings = payload["descriptive_findings"]
    if findings["negative_points_explained_by_scoreable_nearby_upper_node"] == findings[
        "negative_D2_at_p_bar_sizes"
    ]:
        local_reading = (
            "At each of those sizes the stable matching-root-window branch has a "
            "point estimate just above `p_bar`; the negative point is therefore a "
            "local lobe next to a zero crossing, not a reversal of the integrated K2 "
            "response. Node-position standard errors remain in the JSON, so this does "
            "not claim significant ordering of the node and `p_bar`."
        )
    else:
        local_reading = (
            "The full-curve point estimates place a nearby zero above `p_bar`, but "
            "at least one delete-one node spectrum is not scoreable; the local-node "
            "reading therefore remains descriptive rather than branch-stable."
        )
    lines.extend(
        [
            "",
            "## What the full curves add",
            "",
            f"All {len(payload['scoreable_sizes'])} scoreable archives have the same "
            f"integrated K2 direction, `{findings['A2_sign']}`. The point value "
            f"`D2(p_bar)` is negative only at N={','.join(map(str, findings['negative_D2_at_p_bar_sizes']))}. "
            + local_reading,
            "",
            "The JSON stores every Bernstein coefficient, so `D1(p)` and `D2(p)` can "
            "be reconstructed over the entire unit interval. Endpoint zeros are "
            "structural and excluded from the node score. Full-domain and "
            "critical-window branch stability are reported separately; whenever a "
            "delete-one replicate changes the branch count, that node spectrum is "
            "marked `not_scoreable` instead of silently aligning different roots.",
            "",
            "## Center-width identity",
            "",
            "With `C=(K1+K2)/(2(N+1))` and `W=(K2-K1)/(N+1)`, the exact identity "
            "`K2/(N+1)=C+W/2` splits every A2 area into center and width contributions. "
            "The per-size residuals of this identity and of the rank-area formula are "
            "serialized under `identity_audit`.",
            "",
            "## Dependence boundary",
            "",
            "Both directions are deleted as one aligned batch. N=65,85,130,170 also "
            "share a counter stream, so the JSON contains their cross-N covariance for "
            "`A2`, `D2(p_bar)` and `D2'(p_bar)`. Those four views are one dependency "
            "block and must not be added as independent evidence.",
            "",
            "This finite-archive node map fits no free exponent and does not by itself "
            "identify a continuum operator.",
            "",
        ]
    )
    return "\n".join(lines)


def _analyze_archive_star(
    arguments: tuple[str, Mapping[str, Any], Mapping[str, Any]]
) -> dict[str, Any]:
    return analyze_archive(*arguments)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "analysis/activation_curve_nodes_manifest.yaml",
    )
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1:
        raise SystemExit("--workers must be positive")
    try:
        manifest, source_path, source = load_manifest(args.manifest, root)
        worker_args = [
            (str(root), entry, manifest["node_policy"]) for entry in source["runs"]
        ]
        if args.workers == 1:
            results = [_analyze_archive_star(item) for item in worker_args]
        else:
            with ProcessPoolExecutor(
                max_workers=min(args.workers, len(worker_args))
            ) as pool:
                results = list(pool.map(_analyze_archive_star, worker_args))
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
