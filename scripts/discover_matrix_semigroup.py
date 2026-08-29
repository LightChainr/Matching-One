#!/usr/bin/env python3
"""Discover a minimal matrix Gaussian-semigroup model for the P48 channels.

This is deliberately a discovery scorer, not a preregistered test.  It rebuilds
the intrinsic-center slope from the committed full-curve histograms and uses
aligned delete-one-batch covariance throughout.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import analyze_matching_parity_derivatives_fast as p48  # noqa: E402


STATE_METRICS = ("I_S", "I_Du", "T_D", "T_Su")
ORIGINAL_METRICS = ("S_scaled", "Dp_scaled", "D_scaled", "Sp_scaled")
ALPHAS = {"P4_S": 1.0, "P4_D": 13.0 / 8.0,
          "P4_S_prime": 5.0 / 4.0, "P4_D_prime": 5.0 / 8.0}


def _point(first: p48.H, second: p48.H) -> Dict[str, float]:
    center, _delta_cos4, projected = p48.project(first, second)
    first_obs, second_obs = p48.obs(first, center), p48.obs(second, center)
    # M = 2D for each orientation, hence d[(M1+M2)/2]/dp = Dp1+Dp2.
    mbar_prime = first_obs["Dp"] + second_obs["Dp"]
    n = first.n
    slope_scale = mbar_prime / mp.power(n, mp.mpf(3) / 8)
    state = {
        "I_S": n * projected["P4_S"],
        "I_Du": n * projected["P4_D_prime"] / mbar_prime,
        "T_D": mp.power(n, mp.mpf(13) / 8) * projected["P4_D"],
        "T_Su": mp.power(n, mp.mpf(13) / 8) * projected["P4_S_prime"] / mbar_prime,
    }
    original = {
        "S_scaled": n * projected["P4_S"],
        "Dp_scaled": mp.power(n, mp.mpf(5) / 8) * projected["P4_D_prime"],
        "D_scaled": mp.power(n, mp.mpf(13) / 8) * projected["P4_D"],
        "Sp_scaled": mp.power(n, mp.mpf(5) / 4) * projected["P4_S_prime"],
    }
    return {
        **{key: float(value) for key, value in state.items()},
        **{key: float(value) for key, value in original.items()},
        "mbar_prime": float(mbar_prime),
        "slope_scale": float(slope_scale),
        "rho_I": float(projected["P4_D_prime"] / (mbar_prime * projected["P4_S"])),
        "rho_T": float(projected["P4_S_prime"] / (mbar_prime * projected["P4_D"])),
    }


def _covariance(rows: Sequence[Sequence[float]]) -> List[List[float]]:
    count, width = len(rows), len(rows[0])
    means = [math.fsum(row[j] for row in rows) / count for j in range(width)]
    factor = (count - 1.0) / count
    return [[factor * math.fsum((row[i] - means[i]) * (row[j] - means[j]) for row in rows)
             for j in range(width)] for i in range(width)]


def _analyze_aligned(paths: Sequence[Path]) -> Tuple[List[Dict[str, object]], List[List[float]]]:
    datasets = [p48.read(path) for path in paths]
    totals: Dict[int, Tuple[p48.H, p48.H]] = {}
    batches = None
    for data in datasets:
        sizes = {key[0] for key in data}
        if len(sizes) != 1:
            raise ValueError("each input file must contain exactly one N")
        n = sizes.pop()
        grouped = {orientation: [data[key] for key in sorted(data)
                                 if key[1] == orientation]
                   for orientation in ("first", "second")}
        ids = [row.batch for row in grouped["first"]]
        if ids != [row.batch for row in grouped["second"]]:
            raise ValueError(f"N={n}: first/second batch ids differ")
        if batches is None:
            batches = ids
        elif ids != batches:
            raise ValueError("aligned files do not share batch ids")
        totals[n] = (p48.combine(grouped["first"]), p48.combine(grouped["second"]))
    sizes = sorted(totals)
    full = []
    for n in sizes:
        row = _point(*totals[n])
        row.update({"N": n, "source": "+".join(path.name for path in paths)})
        full.append(row)
    metric_order = STATE_METRICS + ORIGINAL_METRICS + ("rho_I", "rho_T", "slope_scale")
    deleted_vectors = []
    assert batches is not None
    for batch in batches:
        vector = []
        for data, n in zip(datasets, sizes):
            first_total, second_total = totals[n]
            first_row = data[(n, "first", batch)]
            second_row = data[(n, "second", batch)]
            deleted = _point(p48.remove(first_total, first_row),
                             p48.remove(second_total, second_row))
            vector.extend(deleted[key] for key in metric_order)
        deleted_vectors.append(vector)
    return full, _covariance(deleted_vectors)


def _block_diag(blocks: Sequence[Sequence[Sequence[float]]]) -> List[List[float]]:
    size = sum(len(block) for block in blocks)
    out = [[0.0] * size for _ in range(size)]
    offset = 0
    for block in blocks:
        for i, row in enumerate(block):
            out[offset + i][offset:offset + len(block)] = row
        offset += len(block)
    return out


def _inverse(matrix: Sequence[Sequence[float]]) -> List[List[float]]:
    return [[float(value) for value in row] for row in mp.inverse(mp.matrix(matrix)).tolist()]


def _matmul(left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]) -> List[List[float]]:
    return [[math.fsum(a * b for a, b in zip(row, column))
             for column in zip(*right)] for row in left]


def _transpose(matrix: Sequence[Sequence[float]]) -> List[List[float]]:
    return [list(column) for column in zip(*matrix)]


def _gls(values: Sequence[float], covariance: Sequence[Sequence[float]],
         design: Sequence[Sequence[float]], names: Sequence[str]) -> Dict[str, object]:
    cov_inv = _inverse(covariance)
    xt_ci = _matmul(_transpose(design), cov_inv)
    parameter_cov = _inverse(_matmul(xt_ci, design))
    influence = _matmul(parameter_cov, xt_ci)
    beta = [row[0] for row in _matmul(influence, [[value] for value in values])]
    fitted = [row[0] for row in _matmul(design, [[value] for value in beta])]
    residual = [value - fit for value, fit in zip(values, fitted)]
    solved = _matmul(cov_inv, [[value] for value in residual])
    chi_square = math.fsum(value * solved[i][0] for i, value in enumerate(residual))
    observations, parameters = len(values), len(beta)
    aic = chi_square + 2 * parameters
    aicc = (aic + 2 * parameters * (parameters + 1) / (observations - parameters - 1)
            if observations > parameters + 1 else math.inf)
    bic = chi_square + parameters * math.log(observations)
    return {
        "parameters": dict(zip(names, beta)),
        "parameter_covariance": parameter_cov,
        "chi_square": chi_square,
        "dof": observations - parameters,
        "parameter_count": parameters,
        "aic": aic,
        "aicc": aicc,
        "bic": bic,
    }


def _extract(points: Sequence[Mapping[str, object]], covariance: Sequence[Sequence[float]],
             metrics: Sequence[str]) -> Tuple[List[float], List[List[float]], List[Tuple[int, str]]]:
    all_order = STATE_METRICS + ORIGINAL_METRICS + ("rho_I", "rho_T", "slope_scale")
    width = len(all_order)
    indices, labels, values = [], [], []
    for point_index, point in enumerate(points):
        for metric in metrics:
            indices.append(point_index * width + all_order.index(metric))
            labels.append((int(point["N"]), metric))
            values.append(float(point[metric]))
    selected = [[covariance[i][j] for j in indices] for i in indices]
    return values, selected, labels


def _design(labels: Sequence[Tuple[int, str]], mode: str) -> Tuple[List[List[float]], List[str]]:
    metrics = sorted({metric for _n, metric in labels})
    terms: List[Tuple[str, str]] = [(metric, "constant") for metric in metrics]
    if mode in ("ordinary", "jordan"):
        feature = "inverse_N" if mode == "ordinary" else "log_N"
        terms += [(metric, feature) for metric in metrics if metric in ("T_D", "T_Su")]
    elif mode in ("ordinary_Su", "jordan_Su"):
        feature = "inverse_N" if mode == "ordinary_Su" else "log_N"
        terms.append(("T_Su", feature))
    rows = []
    for n, metric in labels:
        row = []
        for target, feature in terms:
            if metric != target:
                row.append(0.0)
            elif feature == "constant":
                row.append(1.0)
            elif feature == "inverse_N":
                row.append(1.0 / n)
            else:
                row.append(math.log(n))
        rows.append(row)
    return rows, [f"{metric}:{feature}" for metric, feature in terms]


def _prediction(fit: Mapping[str, object], metric: str, n: int, mode: str) -> Dict[str, float]:
    params = fit["parameters"]
    cov = fit["parameter_covariance"]
    names = list(params)
    vector = [0.0] * len(names)
    vector[names.index(f"{metric}:constant")] = 1.0
    feature_name = None
    if mode.startswith("ordinary") and f"{metric}:inverse_N" in names:
        feature_name, feature = f"{metric}:inverse_N", 1.0 / n
    elif mode.startswith("jordan") and f"{metric}:log_N" in names:
        feature_name, feature = f"{metric}:log_N", math.log(n)
    if feature_name:
        vector[names.index(feature_name)] = feature
    mean = math.fsum(vector[i] * float(params[name]) for i, name in enumerate(names))
    variance = math.fsum(vector[i] * cov[i][j] * vector[j]
                         for i in range(len(vector)) for j in range(len(vector)))
    return {"N": n, "mean": mean, "se": math.sqrt(max(0.0, variance))}


def analyze(paths: Sequence[Sequence[Path]]) -> Dict[str, object]:
    groups = [_analyze_aligned(group) for group in paths]
    points = [point for group_points, _cov in groups for point in group_points]
    covariance = _block_diag([cov for _points, cov in groups])

    values, state_cov, labels = _extract(points, covariance, STATE_METRICS)
    fits = {}
    for mode in ("pure", "ordinary_Su", "jordan_Su", "ordinary", "jordan"):
        design, names = _design(labels, mode)
        fits[mode] = _gls(values, state_cov, design, names)
    baseline = fits["pure"]
    for fit in fits.values():
        fit["delta_aicc_vs_pure"] = fit["aicc"] - baseline["aicc"]
        fit["delta_bic_vs_pure"] = fit["bic"] - baseline["bic"]

    original_values, original_cov, original_labels = _extract(points, covariance, ORIGINAL_METRICS)
    original_design, original_names = _design(original_labels, "pure")
    original_pure = _gls(original_values, original_cov, original_design, original_names)

    rho_values, rho_cov, rho_labels = _extract(points, covariance, ("rho_I", "rho_T"))
    rho_design, rho_names = _design(rho_labels, "pure")
    rho_score = _gls(rho_values, rho_cov, rho_design, rho_names)

    predictions = {}
    for mode in ("ordinary_Su", "jordan_Su", "ordinary", "jordan"):
        predictions[mode] = {
            "T_Su": [_prediction(fits[mode], "T_Su", n, mode) for n in (370, 530)],
            "T_D": [_prediction(fits[mode], "T_D", n, mode) for n in (370, 530)],
        }

    point_by_n = {int(point["N"]): point for point in points}
    discriminators = []
    for n in (65, 85):
        y0 = float(point_by_n[n]["T_Su"])
        y2 = float(point_by_n[2 * n]["T_Su"])
        ordinary_y5 = (8 * y2 - 3 * y0) / 5
        jordan_y5 = y0 + (y2 - y0) * math.log(5) / math.log(2)
        discriminators.append({
            "base_N": n,
            "new_N": 5 * n,
            "observable": "U_T2=N^(13/8)*P4[S_prime]/Mbar_prime",
            "ordinary_completion_U_T2": ordinary_y5,
            "jordan_completion_U_T2": jordan_y5,
            "model_gap_U_T2": jordan_y5 - ordinary_y5,
            "rank1_null_combination": "3*U(N)-8*U(2N)+5*U(5N)",
            "rank1_value": 0.0,
            "jordan_value_from_existing_anchors": 5 * (jordan_y5 - ordinary_y5),
        })

    return {
        "classification": "post-reveal discovery; not preregistered",
        "coordinate_change": {
            "u": "Mbar(p) at the intrinsic center",
            "Mbar_prime_definition": "d[(M_first+M_second)/2]/dp = Dp_first+Dp_second",
            "I_state": ["N*P4[S]", "N*P4[D_prime]/Mbar_prime"],
            "T_state": ["N^(13/8)*P4[D]", "N^(13/8)*P4[S_prime]/Mbar_prime"],
        },
        "points": points,
        "scores": fits,
        "original_coordinate_pure_score": original_pure,
        "intrinsic_coordinate_pure_score": fits["pure"],
        "rho_constant_score": rho_score,
        "predictions": predictions,
        "single_new_experiment_discriminators": discriminators,
        "exact_matrices": {
            "cross_either_affine_on_1_S": [[1, 0], [1, -1]],
            "cross_either_centered_on_1_SminusHalf": [[1, 0], [0, -1]],
            "ordinary_feature_transfer": "R_q=diag(1,q^-1) on (1,N^-1)",
            "jordan_feature_transfer": "J_q=[[1,0],[log(q),1]] on (1,log(N))",
            "composition": "R_q1 R_q2=R_(q1*q2), J_q1 J_q2=J_(q1*q2)",
        },
        "minimal_hidden_dimension": {
            "value": 3,
            "decomposition": "one matching-even rank-1 state plus one matching-odd rank-2 state",
            "readout": "(I_S,I_Du)=h_I z_I; (T_D,T_Su)=H_T z_T",
        },
    }


def _write_csv(path: Path, points: Iterable[Mapping[str, object]]) -> None:
    import csv
    rows = list(points)
    fields = ["N", "source", "mbar_prime", "slope_scale", "rho_I", "rho_T"] + list(STATE_METRICS) + list(ORIGINAL_METRICS)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row[field] for field in fields} for row in rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    mp.mp.dps = 40
    groups = [
        [ROOT / "results/server-20260828/P45-root-amplitude/n65.hist.csv",
         ROOT / "results/server-20260828/P45-root-amplitude/n85.hist.csv",
         ROOT / "results/server-20260828/P49-fullcurve-doubling-100m/raw/n130.hist.csv",
         ROOT / "results/server-20260828/P49-fullcurve-doubling-100m/raw/n170.hist.csv"],
        [ROOT / "results/server-20260828/P43-heldout-fullcurve-500m/raw/n185.hist.csv"],
        [ROOT / "results/server-20260828/P43-heldout-fullcurve-500m/raw/n265.hist.csv"],
    ]
    result = analyze(groups)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "discovery.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    _write_csv(args.output / "state_points.csv", result["points"])
    print(args.output / "discovery.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
