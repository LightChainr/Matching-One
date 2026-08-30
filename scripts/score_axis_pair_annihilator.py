#!/usr/bin/env python3
"""Score adjacent-axis Mertens-Ziff annihilator shapes from coupled rank histograms.

For one upper size L and lower size L-1 define

    G_L(p) = L^(13/4) M_L(p),
    F_L(p) = G_L(p) - G_(L-1)(p).

If

    M_L(pc) = A L^-13/4 [1 + C L^-q + ...],
    M'_L(pc) = B L^3/4 [1 + ...],

then at a fixed reconstruction coordinate p_ref=pc+delta,

    F_L(p_ref)
      = C0 [L^-q-(L-1)^-q]
      + T  [L^4-(L-1)^4]
      + ...,

where T=B*delta is a nuisance.  Thus q can be challenged without externally
supplying pc and without first solving a high-precision accelerated root.

The same histograms also reconstruct the accelerated root F_L(p)=0.  The
coupled adjacent-size engine uses an exact restriction coupling; this scorer
reports its empirical variance gain relative to an independence counterfactual.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Dict, Mapping, Sequence, Tuple

from analyze_p48_retrospective import (
    covariance_of_mean,
    inverse,
    matmul,
    pseudovalues,
    quadratic,
    transpose,
)
from score_p49_fullcurve_doubling import tail_and_derivative


@dataclass
class PairHistogram:
    pair_L: int
    n: int
    L: int
    role: str
    batch: int
    samples: int
    minus: list[int]
    plus: list[int]


def read_pair_histograms(paths: Sequence[Path]) -> Dict[Tuple[int, str, int], PairHistogram]:
    required = {"pair_L", "n", "L", "role", "batch", "samples", "kind", "k", "count"}
    records: Dict[Tuple[int, str, int], PairHistogram] = {}
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            missing = required - set(reader.fieldnames or ())
            if missing:
                raise ValueError(f"{path}: missing fields {sorted(missing)}")
            for raw in reader:
                pair_L = int(raw["pair_L"])
                L = int(raw["L"])
                n = int(raw["n"])
                role = raw["role"]
                batch = int(raw["batch"])
                samples = int(raw["samples"])
                kind = raw["kind"]
                rank = int(raw["k"])
                count = int(raw["count"])
                if role not in ("upper", "lower") or kind not in ("minus", "plus"):
                    raise ValueError("unknown role/kind")
                if pair_L < 3 or L not in (pair_L, pair_L - 1) or n != L * L:
                    raise ValueError("inconsistent adjacent-axis geometry")
                if (role == "upper") != (L == pair_L):
                    raise ValueError("role does not match L")
                key = (pair_L, role, batch)
                if key not in records:
                    records[key] = PairHistogram(
                        pair_L, n, L, role, batch, samples, [0] * (n + 1), [0] * (n + 1)
                    )
                row = records[key]
                if (row.n, row.L, row.samples) != (n, L, samples):
                    raise ValueError("metadata changed within batch")
                getattr(row, kind)[rank] += count
    if not records:
        raise ValueError("no histogram rows")

    for row in records.values():
        if sum(row.minus) != row.samples or sum(row.plus) != row.samples:
            raise ValueError("histogram total differs from samples")

    pair_sizes = sorted({key[0] for key in records})
    for pair_L in pair_sizes:
        signatures = []
        for role in ("upper", "lower"):
            selected = sorted(row for key, row in records.items() if key[:2] == (pair_L, role))
            if not selected:
                raise ValueError(f"missing role {role} at pair L={pair_L}")
            batches = [row.batch for row in selected]
            if batches != list(range(len(batches))) or len(batches) < 2:
                raise ValueError("batches must be contiguous and at least two")
            signatures.append((batches, [row.samples for row in selected]))
        if signatures[0] != signatures[1]:
            raise ValueError("upper/lower batch alignment is absent")
    return records


def grouped(records: Mapping[Tuple[int, str, int], PairHistogram]):
    return {
        pair_L: {
            role: sorted(
                (row for key, row in records.items() if key[:2] == (pair_L, role)),
                key=lambda row: row.batch,
            )
            for role in ("upper", "lower")
        }
        for pair_L in sorted({key[0] for key in records})
    }


def aggregate(rows: Sequence[PairHistogram], omitted: int = -1) -> dict:
    n = rows[0].n
    minus = [0] * (n + 1)
    plus = [0] * (n + 1)
    samples = 0
    for row in rows:
        if row.batch == omitted:
            continue
        samples += row.samples
        for k, count in enumerate(row.minus):
            minus[k] += count
        for k, count in enumerate(row.plus):
            plus[k] += count
    if samples <= 0:
        raise ValueError("cannot omit all batches")
    return {"L": rows[0].L, "n": n, "samples": samples, "minus": minus, "plus": plus}


def matching_value(row: Mapping[str, object], p: float) -> Tuple[float, float]:
    minus, d_minus = tail_and_derivative(row["minus"], int(row["samples"]), p)
    plus, d_plus = tail_and_derivative(row["plus"], int(row["samples"]), p)
    return minus + plus - 1.0, d_minus + d_plus


def solve_root(function) -> float:
    lower, upper = 0.55, 0.63
    f_lower, f_upper = function(lower), function(upper)
    if not f_lower <= 0.0 <= f_upper:
        raise ValueError(f"root is not bracketed: {f_lower}, {f_upper}")
    for _ in range(58):
        midpoint = (lower + upper) / 2.0
        if function(midpoint) < 0.0:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def pair_stat(by_role, p_ref: float, omitted: int = -1) -> dict:
    upper = aggregate(by_role["upper"], omitted)
    lower = aggregate(by_role["lower"], omitted)
    L = int(upper["L"])
    if int(lower["L"]) != L - 1:
        raise ValueError("not an adjacent pair")
    wu = L ** (13.0 / 4.0)
    wl = (L - 1) ** (13.0 / 4.0)

    mu, du = matching_value(upper, p_ref)
    ml, dl = matching_value(lower, p_ref)
    f_ref = wu * mu - wl * ml

    def f_curve(p: float) -> float:
        return wu * matching_value(upper, p)[0] - wl * matching_value(lower, p)[0]

    ann_root = solve_root(f_curve)
    ordinary_upper = solve_root(lambda p: matching_value(upper, p)[0])
    ordinary_lower = solve_root(lambda p: matching_value(lower, p)[0])
    return {
        "pair_L": L,
        "upper_N": int(upper["n"]),
        "lower_N": int(lower["n"]),
        "M_upper_p_ref": mu,
        "M_lower_p_ref": ml,
        "Mprime_upper_p_ref": du,
        "Mprime_lower_p_ref": dl,
        "F_p_ref": f_ref,
        "annihilator_root": ann_root,
        "ordinary_root_upper": ordinary_upper,
        "ordinary_root_lower": ordinary_lower,
    }


def jackknife_pseudovalue_cov(full: Sequence[float], deleted: Sequence[Sequence[float]]):
    columns = [
        pseudovalues(full[j], [row[j] for row in deleted])
        for j in range(len(full))
    ]
    rows = [[columns[j][i] for j in range(len(full))] for i in range(len(deleted))]
    return covariance_of_mean(rows), rows


def pair_jackknife(by_role, p_ref: float) -> dict:
    full = pair_stat(by_role, p_ref)
    batches = len(by_role["upper"])
    deleted = [pair_stat(by_role, p_ref, omitted=i) for i in range(batches)]
    fields = ["M_upper_p_ref", "M_lower_p_ref", "F_p_ref", "annihilator_root"]
    vector = [full[name] for name in fields]
    deleted_vectors = [[row[name] for name in fields] for row in deleted]
    covariance, _ = jackknife_pseudovalue_cov(vector, deleted_vectors)
    se = {fields[i]: math.sqrt(max(covariance[i][i], 0.0)) for i in range(len(fields))}

    L = full["pair_L"]
    wu = L ** (13.0 / 4.0)
    wl = (L - 1) ** (13.0 / 4.0)
    var_mu = covariance[0][0]
    var_ml = covariance[1][1]
    cov_ul = covariance[0][1]
    independent_var_f = wu * wu * var_mu + wl * wl * var_ml
    coupled_var_f = covariance[2][2]
    gain = independent_var_f / coupled_var_f if coupled_var_f > 0 else None
    corr = cov_ul / math.sqrt(var_mu * var_ml) if var_mu > 0 and var_ml > 0 else None

    full.update({
        "batches": batches,
        "se": se,
        "within_pair_M_correlation": corr,
        "F_independence_counterfactual_se": math.sqrt(max(independent_var_f, 0.0)),
        "F_coupled_se": math.sqrt(max(coupled_var_f, 0.0)),
        "F_variance_gain_independent_over_coupled": gain,
    })
    return full


def weighted_linear_fit(
    features: Sequence[Sequence[float]],
    values: Sequence[float],
    variances: Sequence[float],
):
    """Weighted linear fit with explicit design-column conditioning.

    Finite-size basis columns can differ by many orders of magnitude, e.g.
    `[1, L^-10]` or `[Delta L^-q, Delta L^4]`.  Solving the raw normal
    equations can therefore trip a scale-based singularity guard even when the
    columns are linearly independent.  Normalize each design column before the
    solve, then map coefficients/covariance back to the original basis.
    """
    if not features or len(features) != len(values) or len(values) != len(variances):
        raise ValueError("features, values and variances must be nonempty and aligned")
    if any(v <= 0.0 or not math.isfinite(v) for v in variances):
        raise ValueError("variances must be finite and positive")
    p = len(features[0])
    if p == 0 or any(len(row) != p for row in features):
        raise ValueError("design matrix is empty or ragged")

    scales = []
    for column in range(p):
        scale = max(abs(row[column]) for row in features)
        if scale == 0.0 or not math.isfinite(scale):
            raise ValueError("design column is identically zero or nonfinite")
        scales.append(scale)
    scaled_features = [
        [row[column] / scales[column] for column in range(p)]
        for row in features
    ]

    weights = [1.0 / v for v in variances]
    normal = [[0.0] * p for _ in range(p)]
    rhs = [[0.0] for _ in range(p)]
    for x, y, weight in zip(scaled_features, values, weights):
        for i in range(p):
            rhs[i][0] += weight * x[i] * y
            for j in range(p):
                normal[i][j] += weight * x[i] * x[j]

    scaled_covariance = inverse(normal)
    scaled_parameters = matmul(scaled_covariance, rhs)
    beta = [scaled_parameters[j][0] / scales[j] for j in range(p)]
    covariance = [
        [scaled_covariance[i][j] / (scales[i] * scales[j]) for j in range(p)]
        for i in range(p)
    ]
    fitted = [sum(a * b for a, b in zip(x, beta)) for x in features]
    residual = [y - yhat for y, yhat in zip(values, fitted)]
    chi2 = sum(r * r / v for r, v in zip(residual, variances))
    return beta, covariance, residual, chi2


def fit_f_shape(rows: Sequence[dict], q: float, train_max_L: int) -> dict:
    train = [row for row in rows if row["pair_L"] <= train_max_L]
    held = [row for row in rows if row["pair_L"] > train_max_L]
    if len(train) < 3 or not held:
        raise ValueError("need >=3 train pairs and >=1 heldout pair")

    def feature(L: int):
        return [L ** (-q) - (L - 1) ** (-q), L ** 4 - (L - 1) ** 4]

    x_train = [feature(row["pair_L"]) for row in train]
    y_train = [row["F_p_ref"] for row in train]
    v_train = [row["F_coupled_se"] ** 2 for row in train]
    beta, beta_cov, _, train_chi2 = weighted_linear_fit(x_train, y_train, v_train)
    x_held = [feature(row["pair_L"]) for row in held]
    prediction = [sum(a * b for a, b in zip(x, beta)) for x in x_held]
    held_residual = [row["F_p_ref"] - pred for row, pred in zip(held, prediction)]
    # Target pairs are required to use distinct RNG domains, so target covariance
    # is diagonal. Add source-parameter prediction covariance.
    residual_cov = []
    for i, (row_i, x_i) in enumerate(zip(held, x_held)):
        current = []
        for j, (row_j, x_j) in enumerate(zip(held, x_held)):
            source = sum(
                x_i[a] * beta_cov[a][b] * x_j[b]
                for a in range(2) for b in range(2)
            )
            target = row_i["F_coupled_se"] ** 2 if i == j else 0.0
            current.append(source + target)
        residual_cov.append(current)
    return {
        "q": q,
        "w_ann": 4.0 + q,
        "train_L": [row["pair_L"] for row in train],
        "heldout_L": [row["pair_L"] for row in held],
        "parameters": {
            "correction_amplitude": beta[0],
            "thermal_mistuning_nuisance": beta[1],
        },
        "parameter_covariance": beta_cov,
        "train_chi_square": train_chi2,
        "train_df": len(train) - 2,
        "heldout_prediction": prediction,
        "heldout_residual": held_residual,
        "heldout_residual_covariance": residual_cov,
        "heldout_chi_square": quadratic(held_residual, residual_cov),
        "heldout_df": len(held),
    }


def fit_root_power(rows: Sequence[dict], w: float, train_max_L: int) -> dict:
    train = [row for row in rows if row["pair_L"] <= train_max_L]
    held = [row for row in rows if row["pair_L"] > train_max_L]
    if len(train) < 3 or not held:
        raise ValueError("need >=3 train pairs and >=1 heldout pair")

    def feature(L: int):
        return [1.0, L ** (-w)]

    x_train = [feature(row["pair_L"]) for row in train]
    y_train = [row["annihilator_root"] for row in train]
    v_train = [row["se"]["annihilator_root"] ** 2 for row in train]
    beta, beta_cov, _, train_chi2 = weighted_linear_fit(x_train, y_train, v_train)
    x_held = [feature(row["pair_L"]) for row in held]
    pred = [sum(a * b for a, b in zip(x, beta)) for x in x_held]
    residual = [row["annihilator_root"] - yhat for row, yhat in zip(held, pred)]
    residual_cov = []
    for i, (row_i, x_i) in enumerate(zip(held, x_held)):
        current = []
        for j, (row_j, x_j) in enumerate(zip(held, x_held)):
            source = sum(
                x_i[a] * beta_cov[a][b] * x_j[b]
                for a in range(2) for b in range(2)
            )
            target = row_i["se"]["annihilator_root"] ** 2 if i == j else 0.0
            current.append(source + target)
        residual_cov.append(current)
    return {
        "w": w,
        "q": w - 4.0,
        "train_L": [row["pair_L"] for row in train],
        "heldout_L": [row["pair_L"] for row in held],
        "parameters": {"pc": beta[0], "amplitude": beta[1]},
        "parameter_covariance": beta_cov,
        "train_chi_square": train_chi2,
        "train_df": len(train) - 2,
        "heldout_prediction": pred,
        "heldout_residual": residual,
        "heldout_residual_covariance": residual_cov,
        "heldout_chi_square": quadratic(residual, residual_cov),
        "heldout_df": len(held),
    }


def calculate(
    paths: Sequence[Path],
    p_ref: float,
    train_max_L: int,
    q_candidates: Sequence[float],
) -> dict:
    records = read_pair_histograms(paths)
    groups = grouped(records)
    rows = [pair_jackknife(groups[L], p_ref) for L in sorted(groups)]
    f_models = [fit_f_shape(rows, q, train_max_L) for q in q_candidates]
    root_models = [fit_root_power(rows, 4.0 + q, train_max_L) for q in q_candidates]
    return {
        "format_version": 1,
        "p_ref": p_ref,
        "train_max_L": train_max_L,
        "pairs": rows,
        "F_shape_models": f_models,
        "accelerated_root_models": root_models,
        "primary_candidate_order": list(q_candidates),
        "interpretation": {
            "q=2": "relative L^-2 correction / w_ann=6",
            "q=3": "V_<1,4> scalar candidate / w_ann=7",
            "q=4": "nonlinear H4/H12 sideband candidate / w_ann=8",
            "q=6": "next ordinary thermal spin-4 quasiprimary / w_ann=10",
        },
    }


def report(payload: dict) -> str:
    lines = [
        "# Coupled adjacent-axis annihilator score",
        "",
        f"Reconstruction coordinate: `p_ref={payload['p_ref']}`.",
        "",
        "## Pair diagnostics",
        "",
        "| upper L | F(p_ref) | coupled SE | CRN gain | ann root | root SE |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in payload["pairs"]:
        gain = row["F_variance_gain_independent_over_coupled"]
        lines.append(
            f"| {row['pair_L']} | {row['F_p_ref']:.7g} | {row['F_coupled_se']:.3g} | "
            f"{gain:.3g} | {row['annihilator_root']:.15g} | "
            f"{row['se']['annihilator_root']:.3g} |"
        )
    lines.extend([
        "",
        "## Frozen F-shape held-out challenge",
        "",
        "| q | w_ann | heldout chi-square / df |",
        "|---:|---:|---:|",
    ])
    for model in payload["F_shape_models"]:
        lines.append(
            f"| {model['q']:g} | {model['w_ann']:g} | "
            f"{model['heldout_chi_square']:.6g} / {model['heldout_df']} |"
        )
    lines.extend([
        "",
        "## Accelerated-root cross-check",
        "",
        "| w | heldout chi-square / df | fitted pc |",
        "|---:|---:|---:|",
    ])
    for model in payload["accelerated_root_models"]:
        lines.append(
            f"| {model['w']:g} | {model['heldout_chi_square']:.6g} / "
            f"{model['heldout_df']} | {model['parameters']['pc']:.15g} |"
        )
    lines.extend([
        "",
        "Primary scientific score is the fixed-p F-shape challenge because it "
        "removes the leading H4 amplitude algebraically and treats threshold "
        "mistuning as an explicit nuisance. The accelerated root is a cross-check, "
        "not the model-selection source.",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("histograms", nargs="+", type=Path)
    parser.add_argument("--p-ref", type=float, default=0.592746050790)
    parser.add_argument("--train-max-L", type=int, required=True)
    parser.add_argument("--q", nargs="+", type=float, default=[2.0, 3.0, 4.0, 6.0])
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()
    payload = calculate(args.histograms, args.p_ref, args.train_max_L, args.q)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.report.write_text(report(payload), encoding="utf-8")
    print(args.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
