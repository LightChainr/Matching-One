#!/usr/bin/env python3
"""Re-express threshold-rank moments in essential H1 birth coordinates.

For every permutation, C=(K1+K2)/2 is the birth midpoint and
W=K2-K1 is the lifetime of the rank-one ambient-homology state.  Existing
per-batch sufficient statistics determine their first two moments and their
covariance exactly, without replaying random fields.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from statistics import mean, stdev


DEFAULT_INPUTS = (
    "results/server-20260828/P43-heldout-fullcurve-500m/raw/n185.moments.csv",
    "results/server-20260828/P43-heldout-fullcurve-500m/raw/n265.moments.csv",
    "results/server-20260828/P49-fullcurve-doubling-100m/raw/n130.moments.csv",
    "results/server-20260828/P49-fullcurve-doubling-100m/raw/n170.moments.csv",
    "results/server-20260829/P50-n145-n290-fullcurve/raw/n145_100m.moments.csv",
    "results/server-20260829/P50-n145-n290-fullcurve/raw/n290_100m.moments.csv",
    "results/server-20260829/P57-norm5-500m/raw/n325_500m.moments.csv",
    "results/server-20260829/P57-norm5-500m/raw/n425_500m.moments.csv",
    "results/server-20260829/P154-norm4-variance-pilot/raw/n260_10m.moments.csv",
    "results/server-20260829/P154-norm4-variance-pilot/raw/n340_10m.moments.csv",
)


def cos4(a: int, b: int) -> float:
    n = a * a + b * b
    return (a**4 - 6 * a * a * b * b + b**4) / (n * n)


def row_coordinates(row: dict[str, str]) -> dict[str, float]:
    samples = int(row["samples"])
    n = int(row["n"])
    scale = n + 1
    km = int(row["sum_kminus"]) / samples
    kp = int(row["sum_kplus"]) / samples
    km2 = int(row["sum_kminus2"]) / samples
    kp2 = int(row["sum_kplus2"]) / samples
    product = int(row["sum_product"]) / samples
    gap2 = int(row["sum_gap2"]) / samples

    c = (km + kp) / (2 * scale)
    w = (kp - km) / scale
    c2 = (km2 + 2 * product + kp2) / (4 * scale * scale)
    w2 = gap2 / (scale * scale)
    cw = (kp2 - km2) / (2 * scale * scale)
    var_c = max(0.0, c2 - c * c)
    var_w = max(0.0, w2 - w * w)
    cov_cw = cw - c * w
    correlation = (
        cov_cw / math.sqrt(var_c * var_w) if var_c > 0 and var_w > 0 else 0.0
    )
    return {
        "C": c,
        "W": w,
        "var_C": var_c,
        "var_W": var_w,
        "cov_CW": cov_cw,
        "corr_CW": correlation,
    }


def paired_mean_se(values: list[float]) -> tuple[float, float]:
    center = mean(values)
    standard_error = stdev(values) / math.sqrt(len(values)) if len(values) > 1 else 0.0
    return center, standard_error


def analyze_file(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"empty moments file: {path}")

    by_orientation: dict[str, dict[int, tuple[dict[str, str], dict[str, float]]]] = {}
    for row in rows:
        orientation = row["orientation"]
        batch = int(row["batch"])
        if batch in by_orientation.setdefault(orientation, {}):
            raise ValueError(f"duplicate batch {batch} for {orientation} in {path}")
        by_orientation[orientation][batch] = (row, row_coordinates(row))
    if set(by_orientation) != {"first", "second"}:
        raise ValueError(f"expected first/second orientations in {path}")

    batches = sorted(set(by_orientation["first"]) & set(by_orientation["second"]))
    if len(batches) != len(by_orientation["first"]) or len(batches) != len(
        by_orientation["second"]
    ):
        raise ValueError(f"unaligned common-field batches in {path}")

    first_row = by_orientation["first"][batches[0]][0]
    second_row = by_orientation["second"][batches[0]][0]
    n = int(first_row["n"])
    if int(second_row["n"]) != n:
        raise ValueError("same-file orientations must have the same N")
    c4_first = cos4(int(first_row["a"]), int(first_row["b"]))
    c4_second = cos4(int(second_row["a"]), int(second_row["b"]))
    delta_c4 = c4_first - c4_second
    if delta_c4 == 0:
        raise ValueError(f"zero H4 contrast in {path}")

    delta_c_batches = [
        by_orientation["first"][batch][1]["C"]
        - by_orientation["second"][batch][1]["C"]
        for batch in batches
    ]
    delta_w_batches = [
        by_orientation["first"][batch][1]["W"]
        - by_orientation["second"][batch][1]["W"]
        for batch in batches
    ]
    delta_c, se_c = paired_mean_se(delta_c_batches)
    delta_w, se_w = paired_mean_se(delta_w_batches)

    orientation_summary = {}
    for orientation in ("first", "second"):
        coordinate_rows = [
            by_orientation[orientation][batch][1] for batch in batches
        ]
        row0 = by_orientation[orientation][batches[0]][0]
        orientation_summary[orientation] = {
            "a": int(row0["a"]),
            "b": int(row0["b"]),
            "cos4": cos4(int(row0["a"]), int(row0["b"])),
            "mean_C": mean(row["C"] for row in coordinate_rows),
            "mean_W": mean(row["W"] for row in coordinate_rows),
            "mean_var_C": mean(row["var_C"] for row in coordinate_rows),
            "mean_var_W": mean(row["var_W"] for row in coordinate_rows),
            "mean_corr_CW": mean(row["corr_CW"] for row in coordinate_rows),
        }

    return {
        "path": str(path),
        "N": n,
        "batches": len(batches),
        "samples_per_orientation": sum(
            int(by_orientation["first"][batch][0]["samples"]) for batch in batches
        ),
        "orientations": orientation_summary,
        "H4_contrast": {
            "delta_cos4_first_minus_second": delta_c4,
            "delta_C": delta_c,
            "se_delta_C": se_c,
            "z_delta_C": delta_c / se_c if se_c else None,
            "delta_W": delta_w,
            "se_delta_W": se_w,
            "z_delta_W": delta_w / se_w if se_w else None,
            "C_per_delta_cos4": delta_c / delta_c4,
            "W_per_delta_cos4": delta_w / delta_c4,
            "N_13_over_8_C_per_delta_cos4": n ** (13 / 8) * delta_c / delta_c4,
            "N_1_W_per_delta_cos4": n * delta_w / delta_c4,
        },
    }


def loglog_fit(records: list[dict], coordinate: str) -> dict:
    points = []
    for record in records:
        value = record["H4_contrast"][f"{coordinate}_per_delta_cos4"]
        if value:
            points.append((math.log(record["N"]), math.log(abs(value)), math.copysign(1, value)))
    xbar = mean(point[0] for point in points)
    ybar = mean(point[1] for point in points)
    denominator = sum((point[0] - xbar) ** 2 for point in points)
    slope = sum((x - xbar) * (y - ybar) for x, y, _ in points) / denominator
    intercept = ybar - slope * xbar
    residuals = [y - (intercept + slope * x) for x, y, _ in points]
    return {
        "points": len(points),
        "signs": [int(point[2]) for point in points],
        "log_amplitude_slope": slope,
        "effective_decay_exponent": -slope,
        "intercept": intercept,
        "rms_log_residual": math.sqrt(mean(value * value for value in residuals)),
        "boundary": "unweighted retrospective cross-lineage diagnostic, not a frozen exponent score",
    }


def fixed_power_score(records: list[dict], coordinate: str, exponent: float) -> dict:
    amplitudes = []
    for record in records:
        h4 = record["H4_contrast"]
        delta_cos4 = h4["delta_cos4_first_minus_second"]
        value = record["N"] ** exponent * h4[f"delta_{coordinate}"] / delta_cos4
        standard_error = (
            record["N"] ** exponent
            * h4[f"se_delta_{coordinate}"]
            / abs(delta_cos4)
        )
        if standard_error > 0:
            amplitudes.append((value, standard_error, record["N"]))
    precision = sum(1 / standard_error**2 for _, standard_error, _ in amplitudes)
    common = sum(
        value / standard_error**2 for value, standard_error, _ in amplitudes
    ) / precision
    chi2 = sum(
        ((value - common) / standard_error) ** 2
        for value, standard_error, _ in amplitudes
    )
    return {
        "coordinate": coordinate,
        "fixed_decay_exponent": exponent,
        "common_scaled_amplitude": common,
        "common_scaled_amplitude_se": 1 / math.sqrt(precision),
        "chi2": chi2,
        "dof": len(amplitudes) - 1,
        "sizes": [n for _, _, n in amplitudes],
        "scaled_amplitudes": [value for value, _, _ in amplitudes],
        "scaled_standard_errors": [standard_error for _, standard_error, _ in amplitudes],
        "boundary": "retrospective common-amplitude diagnostic across distinct archived lineages",
    }


def build_report(inputs: list[Path]) -> dict:
    records = [analyze_file(path) for path in inputs]
    records.sort(key=lambda record: record["N"])
    high_statistics = [
        record for record in records if "P154-norm4-variance-pilot" not in record["path"]
    ]
    return {
        "schema": "matching-one.essential-birth-clock.v1",
        "issues": [28, 200, 215, 269, 276],
        "coordinates": {
            "K1": "first occupation rank with ambient H1 rank at least one",
            "K2": "first occupation rank with ambient H1 rank two",
            "C": "(K1+K2)/(2(N+1)); complement-odd centered birth clock",
            "W": "(K2-K1)/(N+1); complement-even rank-one lifetime",
            "exact_moment_bridge": "existing sums of K1,K2,K1^2,K2^2,K1*K2 determine mean/variance/covariance of C,W exactly",
        },
        "datasets": records,
        "retrospective_scaling_diagnostics": {
            "C": loglog_fit(records, "C"),
            "W": loglog_fit(records, "W"),
            "high_statistics_fixed_power_scores": {
                "C_N_minus_13_over_8": fixed_power_score(
                    high_statistics, "C", 13 / 8
                ),
                "W_N_minus_1": fixed_power_score(high_statistics, "W", 1),
                "W_N_minus_5_over_4": fixed_power_score(
                    high_statistics, "W", 5 / 4
                ),
                "W_N_minus_11_over_8": fixed_power_score(
                    high_statistics, "W", 11 / 8
                ),
                "W_N_minus_3_over_2": fixed_power_score(
                    high_statistics, "W", 3 / 2
                ),
            },
            "high_statistics_rule": "exclude only the named P154 10M variance-pilot files; retain P43/P49/P50/P57 archives",
        },
        "research_consequence": {
            "state_split": "C and W are topology-defined coordinates, not a fitted latent basis: C records clock translation while W records rank-one persistence/shear.",
            "next_archive": "future streams should add the rank-one line ell and its Smith index iota at activation, producing the marked state (C,W,ell,iota).",
            "claim_boundary": "all coordinate identities are exact; cross-lineage scaling fits are exploratory reuse of correlated archives and do not create new evidence blocks.",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    inputs = args.inputs or [Path(path) for path in DEFAULT_INPUTS]
    rendered = json.dumps(build_report(inputs), indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
