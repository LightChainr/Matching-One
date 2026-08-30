#!/usr/bin/env python3
"""Retrospective fixed-K density-clock score for the two P267 production blocks."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import mpmath as mp

import score_external_observer_transfer as transfer
import score_marked_birth_path as base


MANIFEST_SCHEMA = "matching-one/p267-density-clock-orthogonal-retrospective/v1"
SCORE_SCHEMA = "matching-one/p267-density-clock-orthogonal-score/v1"
METRICS = (
    "P4_raw_re", "P4_raw_im",
    "P4_clock_re", "P4_clock_im",
    "P4_clock_D_re", "P4_clock_D_im",
    "P4_clock_S_re", "P4_clock_S_im",
    "retained_fraction",
    "beta_raw_first", "beta_raw_second",
    "beta0_first", "beta0_second",
    "O0_mean_first", "O0_mean_second",
    "mu_cancellation_first", "mu_cancellation_second",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fixed_k_mu(n: int, k: int) -> mp.mpf:
    """Exact E[V-E+F0|K=k] on a translation-regular N-site square quotient."""
    if not 0 <= k <= n:
        raise ValueError("K outside [0,N]")

    def ratio(order: int) -> mp.mpf:
        if k < order:
            return mp.mpf(0)
        numerator = math.prod(range(k - order + 1, k + 1))
        denominator = math.prod(range(n - order + 1, n + 1))
        return mp.mpf(numerator) / denominator

    return mp.mpf(k) - 2 * n * ratio(2) + n * ratio(4)


def complex_mean(row: base.PathRow, stem: str) -> mp.mpc:
    return mp.mpc(
        row.values[f"{stem}_re"] / row.samples,
        row.values[f"{stem}_im"] / row.samples,
    )


def conditional_orientation(
    rows: Sequence[base.PathRow], p: mp.mpf
) -> dict[str, mp.mpf | mp.mpc]:
    """Remove the exact density clock and both source conditional means."""
    n = len(rows)
    weights = base.binomial_weights(n - 1, p)
    clock_d_terms = []
    clock_s_terms = []
    gram0_terms = []
    norm0_terms = []
    o0_terms = []
    cancellation = []
    for k, (row, weight) in enumerate(zip(rows, weights)):
        samples = mp.mpf(row.samples)
        mu = fixed_k_mu(n, k)
        o = row.values["sum_O_ext"] / samples
        jd = complex_mean(row, "sum_J_D")
        js = complex_mean(row, "sum_J_S")
        ojd = complex_mean(row, "sum_O_ext_J_D")
        ojs = complex_mean(row, "sum_O_ext_J_S")
        gram = mp.mpc(
            row.values["sum_J_D_conj_J_S_re"] / samples,
            row.values["sum_J_D_conj_J_S_im"] / samples,
        )
        norm_s = row.values["sum_abs_J_S2"] / samples

        # Expanded in the frozen order: subtract exact mu from O, then the
        # fixed-K source mean.  The simplified form is recorded as an audit.
        expanded_d = ojd - mu * jd - jd * (o - mu)
        expanded_s = ojs - mu * js - js * (o - mu)
        simplified_d = ojd - o * jd
        simplified_s = ojs - o * js
        cancellation.extend((abs(expanded_d - simplified_d), abs(expanded_s - simplified_s)))
        canonical = weight * mp.mpf(n) / (n - k)
        clock_d_terms.append(canonical * expanded_d)
        clock_s_terms.append(canonical * expanded_s)
        gram0_terms.append(weight * (gram - jd * mp.conj(js)))
        norm0_terms.append(weight * (norm_s - abs(js) ** 2))
        o0_terms.append(weight * (o - mu))

    clock_d = mp.fsum(clock_d_terms)
    clock_s = mp.fsum(clock_s_terms)
    gram0 = mp.fsum(gram0_terms)
    norm0 = mp.fsum(norm0_terms)
    if norm0 <= 0:
        raise ValueError("non-positive fixed-K residual S norm")
    beta0 = mp.re(gram0) / norm0

    raw = base.orientation_observables(rows, p)
    raw_norm = raw["Gram_abs_J_S2"]
    if raw_norm <= 0:
        raise ValueError("non-positive raw S norm")
    beta_raw = mp.re(raw["Gram_J_D_conj_J_S_re"]) / raw_norm
    raw_d = mp.mpc(
        raw["connected_O_ext_J_D_re"], raw["connected_O_ext_J_D_im"]
    )
    raw_s = mp.mpc(
        raw["connected_O_ext_J_S_re"], raw["connected_O_ext_J_S_im"]
    )
    return {
        "raw_D": raw_d,
        "raw_S": raw_s,
        "raw_perp": raw_d - beta_raw * raw_s,
        "clock_D": clock_d,
        "clock_S": clock_s,
        "clock_perp": clock_d - beta0 * clock_s,
        "beta_raw": beta_raw,
        "beta0": beta0,
        "Gram0_imag": mp.im(gram0),
        "O0_mean": mp.fsum(o0_terms),
        "mu_cancellation": max(cancellation),
    }


def projected(
    first: Sequence[base.PathRow], second: Sequence[base.PathRow]
) -> tuple[mp.mpf, dict[str, mp.mpf]]:
    center = base.intrinsic_center(first, second)
    left = conditional_orientation(first, center)
    right = conditional_orientation(second, center)
    leverage = base.cos4(first[0].a, first[0].b) - base.cos4(second[0].a, second[0].b)
    if leverage == 0:
        raise ValueError("zero H4 leverage")

    def p4(name: str) -> mp.mpc:
        return (left[name] - right[name]) / leverage

    raw = p4("raw_perp")
    clock = p4("clock_perp")
    clock_d = p4("clock_D")
    clock_s = p4("clock_S")
    if abs(raw) == 0:
        raise ValueError("zero raw P4 response")
    return center, {
        "P4_raw_re": mp.re(raw), "P4_raw_im": mp.im(raw),
        "P4_clock_re": mp.re(clock), "P4_clock_im": mp.im(clock),
        "P4_clock_D_re": mp.re(clock_d), "P4_clock_D_im": mp.im(clock_d),
        "P4_clock_S_re": mp.re(clock_s), "P4_clock_S_im": mp.im(clock_s),
        "retained_fraction": abs(clock) / abs(raw),
        "beta_raw_first": left["beta_raw"], "beta_raw_second": right["beta_raw"],
        "beta0_first": left["beta0"], "beta0_second": right["beta0"],
        "O0_mean_first": left["O0_mean"], "O0_mean_second": right["O0_mean"],
        "mu_cancellation_first": left["mu_cancellation"],
        "mu_cancellation_second": right["mu_cancellation"],
    }


def subtract_group(total: Sequence[base.PathRow], omitted: Sequence[base.PathRow]) -> list[base.PathRow]:
    if len(total) != len(omitted):
        raise ValueError("delete-one path lengths differ")
    output = []
    for all_row, one_row in zip(total, omitted):
        if (all_row.n, all_row.a, all_row.b, all_row.orientation, all_row.k) != (
            one_row.n, one_row.a, one_row.b, one_row.orientation, one_row.k
        ):
            raise ValueError("delete-one path identities differ")
        output.append(base.PathRow(
            n=all_row.n,
            a=all_row.a,
            b=all_row.b,
            orientation=all_row.orientation,
            batch=-1,
            samples=all_row.samples - one_row.samples,
            k=all_row.k,
            values={name: all_row.values[name] - one_row.values[name] for name in base.VALUE_COLUMNS},
        ))
    return output


def delete_one_covariance(rows: Sequence[Mapping[str, mp.mpf]]) -> list[list[mp.mpf]]:
    count = len(rows)
    means = {name: mp.fsum(row[name] for row in rows) / count for name in METRICS}
    factor = mp.mpf(count - 1) / count
    return [[
        factor * mp.fsum(
            (row[left] - means[left]) * (row[right] - means[right]) for row in rows
        )
        for right in METRICS
    ] for left in METRICS]


def score_size(groups: Mapping[tuple[int, str, int], list[base.PathRow]]) -> dict[str, object]:
    sizes = {key[0] for key in groups}
    if len(sizes) != 1:
        raise ValueError("path input must contain one size")
    n = sizes.pop()
    batches = sorted(
        {key[2] for key in groups if key[1] == "first"}
        & {key[2] for key in groups if key[1] == "second"}
    )
    if batches != list(range(len(batches))) or len(batches) < 2:
        raise ValueError("aligned zero-based batches required")
    totals = {
        side: base.combine([groups[(n, side, batch)] for batch in batches])
        for side in ("first", "second")
    }
    center, point = projected(totals["first"], totals["second"])
    delete_centers = []
    delete_rows = []
    for omitted in batches:
        leave = {
            side: subtract_group(totals[side], groups[(n, side, omitted)])
            for side in ("first", "second")
        }
        one_center, one_point = projected(leave["first"], leave["second"])
        delete_centers.append(one_center)
        delete_rows.append(one_point)
    covariance = delete_one_covariance(delete_rows)
    standard_error = {
        name: mp.sqrt(max(mp.mpf(0), covariance[index][index]))
        for index, name in enumerate(METRICS)
    }
    center_mean = mp.fsum(delete_centers) / len(delete_centers)
    center_se = mp.sqrt(
        mp.mpf(len(delete_centers) - 1) / len(delete_centers)
        * mp.fsum((value - center_mean) ** 2 for value in delete_centers)
    )
    return {
        "N": n,
        "batches": len(batches),
        "samples_per_orientation": totals["first"][0].samples,
        "intrinsic_center": base._text(center),
        "intrinsic_center_delete_one_se": base._text(center_se),
        "metric_order": list(METRICS),
        "point": {name: base._text(point[name]) for name in METRICS},
        "standard_error": {name: base._text(standard_error[name]) for name in METRICS},
        "delete_one_covariance": [[base._text(value) for value in row] for row in covariance],
    }


def selected_covariance(report: Mapping[str, object], names: Sequence[str]) -> list[list[float]]:
    order = report["metric_order"]
    source = report["delete_one_covariance"]
    indices = [order.index(name) for name in names]
    return [[float(source[i][j]) for j in indices] for i in indices]


def complex_point(report: Mapping[str, object], prefix: str) -> complex:
    return complex(float(report["point"][prefix + "_re"]), float(report["point"][prefix + "_im"]))


def chi2_survival_2d(value: float) -> float:
    return math.exp(-0.5 * value)


def log10_chi2_survival_2d(value: float) -> float:
    return -0.5 * value / math.log(10.0)


def nonzero_summary(report: Mapping[str, object], prefix: str, alpha: float) -> dict[str, object]:
    names = [prefix + "_re", prefix + "_im"]
    value = complex_point(report, prefix)
    covariance = selected_covariance(report, names)
    chi2 = transfer.mahalanobis(value, covariance)
    p = chi2_survival_2d(chi2)
    return {
        "complex": [value.real, value.imag],
        "covariance_re_im": covariance,
        "mahalanobis_chi2_2d": chi2,
        "p_value": p,
        "log10_p_value": log10_chi2_survival_2d(chi2),
        "decision": "nonzero at frozen alpha" if p < alpha else "zero null survives",
    }


def sum_covariance(first: Sequence[Sequence[float]], second: Sequence[Sequence[float]]) -> list[list[float]]:
    return [[first[i][j] + second[i][j] for j in range(len(first[i]))] for i in range(len(first))]


def block_summary(by_size: Mapping[str, Mapping[str, object]], alpha: float) -> dict[str, object]:
    small, large = by_size["N325"], by_size["N425"]
    z1, z2 = complex_point(small, "P4_clock"), complex_point(large, "P4_clock")
    c1 = selected_covariance(small, ("P4_clock_re", "P4_clock_im"))
    c2 = selected_covariance(large, ("P4_clock_re", "P4_clock_im"))
    return {
        "single_size": {
            "N325": nonzero_summary(small, "P4_clock", alpha),
            "N425": nonzero_summary(large, "P4_clock", alpha),
        },
        "transfer_N425_over_N325": transfer.ratio_summary(
            z1, z2, transfer.block_diagonal(c1, c2)
        ),
        "retained_fraction": {
            name: {
                "point": float(report["point"]["retained_fraction"]),
                "standard_error": float(report["standard_error"]["retained_fraction"]),
            }
            for name, report in by_size.items()
        },
    }


def compatibility(
    first: Mapping[str, object], second: Mapping[str, object], prefix: str, alpha: float
) -> dict[str, object]:
    names = (prefix + "_re", prefix + "_im")
    z1, z2 = complex_point(first, prefix), complex_point(second, prefix)
    covariance = sum_covariance(
        selected_covariance(first, names), selected_covariance(second, names)
    )
    residual = z2 - z1
    chi2 = transfer.mahalanobis(residual, covariance)
    p = chi2_survival_2d(chi2)
    return {
        "residual_second_minus_first": [residual.real, residual.imag],
        "covariance_re_im": covariance,
        "mahalanobis_chi2_2d": chi2,
        "p_value": p,
        "decision": "blocks differ at frozen alpha" if p < alpha else "independent blocks compatible",
    }


def validate_and_load(root: Path, manifest_path: Path) -> tuple[dict, dict[str, dict[str, dict]]]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != MANIFEST_SCHEMA or manifest.get("status") != "retrospective_protocol_locked_before_density_clock_score":
        raise ValueError("manifest is not the frozen P267 density-clock contract")
    intervals = []
    loaded: dict[str, dict[str, dict]] = {}
    required = set(manifest["required_path_columns"])
    for block, by_size in manifest["blocks"].items():
        loaded[block] = {}
        for size, item in by_size.items():
            path = root / item["path"]
            metadata_path = root / item["metadata"]
            if sha256(path) != item["path_sha256"] or sha256(metadata_path) != item["metadata_sha256"]:
                raise ValueError(f"frozen input digest mismatch for {block}/{size}")
            with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
                reader = csv.reader(handle)
                header = next(reader)
                rows = sum(1 for _ in reader)
            if not required.issubset(header) or rows != int(item["path_rows"]):
                raise ValueError(f"path schema/row count mismatch for {block}/{size}")
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            observed = (
                int(metadata["seed"]), int(metadata["replica_counter_first"]),
                int(metadata["replica_counter_last_exclusive"]), str(metadata["git_commit"]),
            )
            expected = (int(item["seed"]), *map(int, item["counter"]), str(item["producer_commit"]))
            if observed != expected:
                raise ValueError(f"metadata identity mismatch for {block}/{size}")
            intervals.append((block, size, expected[1], expected[2]))
            loaded[block][size] = base.read_path(path)
    for i, left in enumerate(intervals):
        for right in intervals[i + 1:]:
            if not (left[3] <= right[2] or right[3] <= left[2]):
                raise ValueError(f"counter domains overlap: {left}, {right}")
    return manifest, loaded


def build(root: Path, manifest_path: Path) -> dict[str, object]:
    manifest, loaded = validate_and_load(root, manifest_path)
    alpha = float(manifest["decision_alpha"])
    scores = {
        block: {size: score_size(groups) for size, groups in by_size.items()}
        for block, by_size in loaded.items()
    }
    summaries = {block: block_summary(by_size, alpha) for block, by_size in scores.items()}
    blocks = list(scores)
    comparisons = {
        size: {
            "clock": compatibility(scores[blocks[0]][size], scores[blocks[1]][size], "P4_clock", alpha),
            "raw": compatibility(scores[blocks[0]][size], scores[blocks[1]][size], "P4_raw", alpha),
        }
        for size in ("N325", "N425")
    }
    return {
        "schema": SCORE_SCHEMA,
        "status": "retrospective frozen density-clock score complete",
        "manifest": str(manifest_path.relative_to(root)),
        "manifest_sha256": sha256(manifest_path),
        "decision_alpha": alpha,
        "definitions": manifest["fixed_K_contract"],
        "block_scores": scores,
        "block_summaries": summaries,
        "independent_block_compatibility": comparisons,
        "sufficiency": manifest["sufficiency"],
        "forbidden_observer": "O_far conditional centering was not defined or scored",
        "interpretation_boundary": "A surviving fixed-K residual is a lattice mechanism result, not a continuum field identity or exponent fit.",
    }


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path,
        default=root / "analysis/p267_density_clock_orthogonal_20260830.json",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dps", type=int, default=50)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    result = build(root, args.manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
