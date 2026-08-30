#!/usr/bin/env python3
"""Score aligned P321 equal-area rectangle roots and fixed N^-2/N^-3 flow."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import mpmath as mp

from analyze_threshold_rank_orientation import add_histograms, read_histograms, validate_moments
from analyze_threshold_ranks import matching_root


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DESIGN = ROOT / "analysis/p321_equal_area_rectangle_design.json"
DEFAULT_ORACLE = ROOT / "predictions/p321_thermal_q4_aspect_ratio_20260830.json"
RHO_ORDER = ("1", "16/9", "9/4", "4", "9")


def _jackknife_covariance(left: Sequence[float], right: Sequence[float]) -> float:
    if len(left) != len(right) or len(left) < 2:
        raise ValueError("aligned nontrivial jackknife vectors required")
    mean_left = math.fsum(left) / len(left)
    mean_right = math.fsum(right) / len(right)
    return (len(left) - 1) / len(left) * math.fsum(
        (x - mean_left) * (y - mean_right) for x, y in zip(left, right)
    )


def _orientation_bytes(path: Path, orientation: str) -> bytes:
    lines = path.read_bytes().splitlines(keepends=True)
    if not lines:
        raise ValueError(f"empty CSV: {path}")
    header = next(csv.reader([lines[0].decode("utf-8")]))
    try:
        index = header.index("orientation")
    except ValueError as error:
        raise ValueError("CSV lacks orientation column") from error
    selected = [lines[0]]
    for line in lines[1:]:
        row = next(csv.reader([line.decode("utf-8")]))
        if row[index] == orientation:
            selected.append(line)
    if len(selected) == 1:
        raise ValueError(f"no {orientation} rows in {path}")
    return b"".join(selected)


def _aggregate_root(
    records: Mapping[tuple[int, str, int], Mapping[str, Any]],
    n: int,
    orientation: str,
    drop_batch: int | None = None,
) -> float:
    selected = [
        dict(records[key]) for key in sorted(records)
        if key[0] == n and key[1] == orientation and key[2] != drop_batch
    ]
    if len(selected) < 2:
        raise ValueError("at least two retained batches required")
    samples = sum(int(row["samples"]) for row in selected)
    minus = add_histograms(selected, "minus")
    plus = add_histograms(selected, "plus")
    return float(matching_root(n, samples, minus, plus))


def score_campaign(campaign_dir: Path, design: Mapping[str, Any]) -> dict[str, Any]:
    manifest = json.loads((campaign_dir / "campaign.json").read_text(encoding="utf-8"))
    if manifest.get("schema") != "matching-one/p321-equal-area-campaign/v1":
        raise ValueError("wrong campaign schema")
    n = int(manifest["N"])
    design_rows = {
        row["aspect_ratio"]: row for row in design["rows"] if int(row["N"]) == n
    }
    if set(design_rows) != {"1/1", "16/9", "9/4", "4/1", "9/1"}:
        raise ValueError("design aspect-ratio set drift")
    aliases = {"16/9": "16/9", "9/4": "9/4", "4/1": "4", "9/1": "9"}
    expected_common = {
        "samples_per_pair": int(manifest["samples_per_shape"]),
        "batches": int(manifest["batches"]),
        "seed": int(manifest["seed"]),
        "replica_counter_first": int(manifest["replica_counter_first"]),
        "replica_counter_last_exclusive": int(manifest["replica_counter_last_exclusive"]),
    }
    loaded: dict[str, Mapping[tuple[int, str, int], Mapping[str, Any]]] = {}
    square_hist_bytes = None
    square_moment_bytes = None
    elapsed_seconds = 0.0
    for run in manifest["runs"]:
        raw_rho = run["rho"]
        rho = aliases.get(raw_rho)
        if rho is None or rho in loaded:
            raise ValueError("unknown or duplicate rectangle rho")
        histogram = campaign_dir / run["histogram"]
        moments = campaign_dir / run["moments"]
        metadata = json.loads((campaign_dir / run["metadata"]).read_text(encoding="utf-8"))
        for key, expected in expected_common.items():
            if metadata.get(key) != expected:
                raise ValueError(f"metadata drift for {rho}: {key}")
        design_metadata = metadata["designs"][0]
        if design_metadata["first_period_matrix"] != design_rows["1/1"]["period_matrix_row_major"]:
            raise ValueError("square period matrix drift")
        if design_metadata["second_period_matrix"] != design_rows[raw_rho]["period_matrix_row_major"]:
            raise ValueError("rectangle period matrix drift")
        current_square_hist = _orientation_bytes(histogram, "first")
        current_square_moments = _orientation_bytes(moments, "first")
        if square_hist_bytes is None:
            square_hist_bytes = current_square_hist
            square_moment_bytes = current_square_moments
        elif current_square_hist != square_hist_bytes or current_square_moments != square_moment_bytes:
            raise ValueError("repeated square histogram/moments are not byte-identical")
        records = read_histograms(histogram)
        validate_moments(moments, records)
        loaded[rho] = records
        elapsed_seconds += float(metadata["elapsed_seconds"])
    if set(loaded) != {"16/9", "9/4", "4", "9"}:
        raise ValueError("campaign lacks a rectangle")

    baseline = loaded["16/9"]
    batch_ids = sorted({key[2] for key in baseline if key[1] == "first"})
    if batch_ids != list(range(int(manifest["batches"]))):
        raise ValueError("batch ids are not complete and zero-based")

    def root_vector(drop_batch: int | None = None) -> list[float]:
        values = [_aggregate_root(baseline, n, "first", drop_batch)]
        values.extend(_aggregate_root(loaded[rho], n, "second", drop_batch) for rho in RHO_ORDER[1:])
        return values

    roots = root_vector()
    delete_one = [root_vector(batch) for batch in batch_ids]
    covariance = [
        [
            _jackknife_covariance(
                [row[i] for row in delete_one], [row[j] for row in delete_one]
            )
            for j in range(len(RHO_ORDER))
        ]
        for i in range(len(RHO_ORDER))
    ]
    contrast_matrix = [[-1.0 if j == 0 else (1.0 if j == i else 0.0) for j in range(5)] for i in range(1, 5)]
    contrasts = [roots[i] - roots[0] for i in range(1, 5)]
    contrast_covariance = [
        [
            math.fsum(contrast_matrix[i][a] * covariance[a][b] * contrast_matrix[j][b] for a in range(5) for b in range(5))
            for j in range(4)
        ]
        for i in range(4)
    ]
    return {
        "N": n,
        "root_order": list(RHO_ORDER),
        "roots": roots,
        "root_covariance": covariance,
        "root_standard_errors": [math.sqrt(max(0.0, covariance[i][i])) for i in range(5)],
        "contrast_order": list(RHO_ORDER[1:]),
        "contrasts_to_square": contrasts,
        "contrast_covariance": contrast_covariance,
        "contrast_standard_errors": [math.sqrt(max(0.0, contrast_covariance[i][i])) for i in range(4)],
        "batches": len(batch_ids),
        "samples_per_shape": int(manifest["samples_per_shape"]),
        "square_histograms_byte_identical": True,
        "square_moments_byte_identical": True,
        "elapsed_seconds_all_four_pairs": elapsed_seconds,
    }


def _block_inverse(covariance: Sequence[Sequence[float]]) -> mp.matrix:
    return mp.inverse(mp.matrix(covariance))


def fit_fixed_model(
    campaigns: Sequence[Mapping[str, Any]], oracle: Mapping[str, Any]
) -> dict[str, Any]:
    ordered = sorted(campaigns, key=lambda row: int(row["N"]))
    if len(ordered) < 3:
        return {
            "status": "insufficient_scales_for_fixed_N^-2_N^-3_fit",
            "required_N": [144, 576, 1296],
            "available_N": [int(row["N"]) for row in ordered],
        }
    if [int(row["N"]) for row in ordered] != [144, 576, 1296]:
        raise ValueError("fixed P321 fit requires exactly N=144,576,1296")
    parameter_count = 1 + 2 * len(RHO_ORDER)
    normal = mp.matrix(parameter_count, parameter_count)
    rhs = mp.matrix(parameter_count, 1)
    for campaign in ordered:
        n = int(campaign["N"])
        inverse = _block_inverse(campaign["root_covariance"])
        design_rows = []
        for rho_index in range(len(RHO_ORDER)):
            row = [mp.mpf(0)] * parameter_count
            row[0] = 1
            row[1 + rho_index] = mp.mpf(n) ** -2
            row[1 + len(RHO_ORDER) + rho_index] = mp.mpf(n) ** -3
            design_rows.append(row)
        x = mp.matrix(design_rows)
        y = mp.matrix([float(value) for value in campaign["roots"]])
        normal += x.T * inverse * x
        rhs += x.T * inverse * y
    parameter_covariance = mp.inverse(normal)
    beta = parameter_covariance * rhs
    chi2 = mp.mpf(0)
    for campaign in ordered:
        inverse = _block_inverse(campaign["root_covariance"])
        residual = mp.matrix(5, 1)
        n = mp.mpf(campaign["N"])
        for rho_index in range(5):
            prediction = beta[0] + beta[1 + rho_index] * n**-2 + beta[6 + rho_index] * n**-3
            residual[rho_index] = mp.mpf(campaign["roots"][rho_index]) - prediction
        chi2 += (residual.T * inverse * residual)[0]
    c_n = [float(beta[1 + i]) for i in range(5)]
    d_n = [float(beta[6 + i]) for i in range(5)]
    c_covariance = [[float(parameter_covariance[1 + i, 1 + j]) for j in range(5)] for i in range(5)]
    rho_values = [1.0, 16 / 9, 9 / 4, 4.0, 9.0]
    c_width = [value / rho**2 for value, rho in zip(c_n, rho_values)]
    c_width_covariance = [
        [c_covariance[i][j] / (rho_values[i] ** 2 * rho_values[j] ** 2) for j in range(5)]
        for i in range(5)
    ]

    oracle_by_rho = {row["rho"]: row for row in oracle["records"]}
    primary_indices = (1, 2, 3)
    factors = [float(oracle_by_rho[RHO_ORDER[i]]["width_C_over_square_C"]) for i in primary_indices]
    residuals = [c_width[i] - factor * c_width[0] for i, factor in zip(primary_indices, factors)]
    residual_covariance = []
    for i, factor_i in zip(primary_indices, factors):
        row = []
        for j, factor_j in zip(primary_indices, factors):
            row.append(
                c_width_covariance[i][j]
                - factor_i * c_width_covariance[0][j]
                - factor_j * c_width_covariance[i][0]
                + factor_i * factor_j * c_width_covariance[0][0]
            )
        residual_covariance.append(row)
    e4_chi2 = float((mp.matrix(residuals).T * mp.inverse(mp.matrix(residual_covariance)) * mp.matrix(residuals))[0])
    return {
        "status": "fixed_N^-2_N^-3_gls_fit",
        "model": "p(N,rho)=pc+C_N(rho) N^-2+D_N(rho) N^-3",
        "free_exponent_fit": False,
        "pc": float(beta[0]),
        "C_N_order": list(RHO_ORDER),
        "C_N": c_n,
        "D_N": d_n,
        "C_N_covariance": c_covariance,
        "C_width_equals_C_N_over_rho_squared": c_width,
        "C_width_covariance": c_width_covariance,
        "fit_chi_square": float(chi2),
        "fit_degrees_of_freedom": 15 - parameter_count,
        "conditional_thermal_Q4_E4_score": {
            "status": "conditional_on_Virasoro_transparent_homology_projector",
            "rho_order": [RHO_ORDER[i] for i in primary_indices],
            "frozen_width_C_over_square_C": factors,
            "residuals": residuals,
            "residual_covariance": residual_covariance,
            "chi_square": e4_chi2,
            "degrees_of_freedom": len(primary_indices),
            "endpoint_rho_9_is_diagnostic_only": True,
        },
    }


def score(
    campaign_dirs: Sequence[Path], design_path: Path = DEFAULT_DESIGN,
    oracle_path: Path = DEFAULT_ORACLE,
) -> dict[str, Any]:
    design = json.loads(design_path.read_text(encoding="utf-8"))
    oracle = json.loads(oracle_path.read_text(encoding="utf-8"))
    campaigns = [score_campaign(path, design) for path in campaign_dirs]
    if len({row["N"] for row in campaigns}) != len(campaigns):
        raise ValueError("duplicate N campaign")
    return {
        "schema": "matching-one/p321-equal-area-covariance-score/v1",
        "status": "variance_smoke" if all(row["samples_per_shape"] <= 100_000 for row in campaigns) else "scored",
        "root_semantics": "individual root of P2-P0; contrasts are rectangle minus square",
        "campaigns": campaigns,
        "scale_fit": fit_fixed_model(campaigns, oracle),
        "model_selection_rule": "fixed N^-2/N^-3 only; E4 conditional prediction frozen before target data and never tuned on a smoke",
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-dir", action="append", required=True, type=Path)
    parser.add_argument("--design", type=Path, default=DEFAULT_DESIGN)
    parser.add_argument("--oracle", type=Path, default=DEFAULT_ORACLE)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    mp.mp.dps = 60
    result = score(args.campaign_dir, args.design, args.oracle)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, allow_nan=False, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
