#!/usr/bin/env python3
"""Issue #8: L^{-4} / L^{-6} cancellation audit of matched estimators.

Post-processing of Issue #7 exact microcanonical totals and Issue #9
Newman-Ziff roots. No Monte Carlo is run.

Ordinary correction models (on finite-size roots p_L):

    A: p_L = pc + a4 L^{-4}
    B: p_L = pc + a4 L^{-4} + a6 L^{-6}
    C: p_L = pc + a4 L^{-4} + a6 L^{-6} + a8 L^{-8}

Richardson cancellation uses weights that depend only on the lattice sizes.
Uncertainty of a cancelled combination is the batch SD/SE of
p_super[b] = sum_i w_i p_{L_i}[b], which keeps the matched-root covariance
already present in each p_L[b]. Training residuals are diagnostics, not
confidence intervals. Min-max intercepts are model spread, not a CI.
Ranking of estimators uses withheld-size prediction error only.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import platform
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Iterable, Sequence

import mpmath as mp

SCRIPT_DIR = Path(__file__).resolve().parent
REPO = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import analyze_newman_ziff as nz  # noqa: E402
import matching_annihilator as ma  # noqa: E402

SIZES_SHARED = (16, 24, 32, 48, 64, 96, 128, 192, 256)
SIZES_INDEPENDENT = (32, 64, 128, 256)
OBSERVABLES = ("H", "V", "either", "both", "M")
WRAP_OBS = ("H", "V", "either", "both")
MODELS = {
    "A": (4,),
    "B": (4, 6),
    "C": (4, 6, 8),
}
ROLLING_FOLDS = (
    (64, 96),
    (96, 128),
    (128, 192),
    (192, 256),
)
N_BATCHES = 20
N_DIGITS = 18
AMP_DIGITS = 18
WEIGHT_DIGITS = 24


def fmt(value, digits: int = N_DIGITS) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    try:
        if isinstance(value, float) and math.isnan(value):
            return "nan"
    except TypeError:
        pass
    try:
        mpv = mp.mpf(value)
    except (ValueError, TypeError):
        return "nan"
    if mpv == mp.mpf("nan") or mp.isinf(mpv):
        if mp.isinf(mpv):
            return "inf" if mpv > 0 else "-inf"
        return "nan"
    return mp.nstr(mpv, n=digits, strip_zeros=False)


def to_mpf(value) -> mp.mpf:
    if isinstance(value, mp.mpf):
        return value
    if value is None:
        return mp.mpf("nan")
    if isinstance(value, float) and math.isnan(value):
        return mp.mpf("nan")
    try:
        return mp.mpf(value)
    except (ValueError, TypeError):
        return mp.mpf("nan")


def is_finite(value) -> bool:
    mpv = to_mpf(value)
    return mpv == mpv and not mp.isinf(mpv)


def sign_label(value) -> str:
    mpv = to_mpf(value)
    if not is_finite(mpv):
        return "nan"
    floor = mp.mpf("1e-40")
    if abs(mpv) < floor:
        return "0"
    return "+" if mpv > 0 else "-"


def mean_sd_se(values: Sequence[mp.mpf]) -> tuple[mp.mpf, mp.mpf, mp.mpf, int]:
    finite = [to_mpf(v) for v in values if is_finite(v)]
    m = len(finite)
    if m == 0:
        nan = mp.mpf("nan")
        return nan, nan, nan, 0
    mu = mp.fsum(finite) / m
    if m == 1:
        return mu, mp.mpf(0), mp.mpf(0), 1
    var = mp.fsum((v - mu) ** 2 for v in finite) / (m - 1)
    sd = mp.sqrt(var)
    se = sd / mp.sqrt(m)
    return mu, sd, se, m


def rms(values: Iterable[mp.mpf]) -> mp.mpf:
    seq = [to_mpf(v) for v in values]
    if not seq:
        return mp.mpf("nan")
    return mp.sqrt(mp.fsum(v * v for v in seq) / len(seq))


def design_row(size: int, powers: Sequence[int]) -> list[mp.mpf]:
    L = mp.mpf(size)
    return [mp.mpf(1), *(L ** (-power) for power in powers)]


def matrix_condition(matrix: mp.matrix) -> mp.mpf:
    if matrix.rows == 0 or matrix.cols == 0:
        return mp.mpf("nan")
    gram = matrix.T * matrix if matrix.rows >= matrix.cols else matrix * matrix.T
    try:
        eigenvalues = mp.eigsy(gram, eigvals_only=True)
    except (ValueError, ZeroDivisionError):
        return mp.inf
    positive = [value for value in eigenvalues if value > 0]
    if len(positive) < gram.rows:
        return mp.inf
    return mp.sqrt(max(positive) / min(positive))


def fit_correction(
    sizes: Sequence[int],
    values: Sequence[mp.mpf],
    powers: Sequence[int],
) -> dict:
    n_params = 1 + len(powers)
    if len(sizes) < n_params:
        return {"status": "too_few_points"}
    matrix = mp.matrix([design_row(size, powers) for size in sizes])
    target = mp.matrix(list(values))
    try:
        coefficients, _resid = mp.qr_solve(matrix, target)
    except (ValueError, ZeroDivisionError):
        return {"status": "singular"}
    predicted = [mp.fsum(term * coefficients[i] for i, term in enumerate(design_row(size, powers))) for size in sizes]
    residuals = [pred - obs for pred, obs in zip(predicted, values)]
    coeff = {
        "pc": coefficients[0],
        "a4": mp.mpf("nan"),
        "a6": mp.mpf("nan"),
        "a8": mp.mpf("nan"),
        "condition_number": matrix_condition(matrix),
        "training_rmse": rms(residuals),
        "training_max_abs": max(abs(r) for r in residuals),
        "n_train": len(sizes),
        "status": "ok",
        "coefficients": coefficients,
        "powers": powers,
    }
    power_to_idx = {power: 1 + i for i, power in enumerate(powers)}
    for power, name in ((4, "a4"), (6, "a6"), (8, "a8")):
        if power in power_to_idx:
            coeff[name] = coefficients[power_to_idx[power]]
    return coeff


def predict_p(size: int, fit: dict) -> mp.mpf:
    if fit.get("status") != "ok":
        return mp.mpf("nan")
    row = design_row(size, fit["powers"])
    return mp.fsum(term * fit["coefficients"][i] for i, term in enumerate(row))


def solve_cancel_weights(sizes: Sequence[int], exponents: Sequence[mp.mpf]) -> dict:
    matrix = ma.constraint_matrix(sizes, exponents)
    weights = ma.solve_weights(matrix)
    weight_list = [weights[i] for i in range(len(sizes))]
    l1 = mp.fsum(abs(w) for w in weight_list)
    l2 = mp.sqrt(mp.fsum(w * w for w in weight_list))
    return {
        "weights": weight_list,
        "sum_abs_weights": l1,
        "l2_norm": l2,
        "condition_number": ma.matrix_condition(matrix),
        "status": "ok",
    }


@dataclass
class RootTables:
    pooled: dict  # (L, obs, mode) -> dict
    batches: dict  # (L, obs, mode, batch) -> mpf
    issue7: list


def load_issue7(path: Path) -> list[dict]:
    rows = []
    for csv_path in sorted(path.glob("L*_microcanonical.csv")):
        with csv_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            n_rows = 0
            k_max = -1
            fields = reader.fieldnames or []
            for row in reader:
                n_rows += 1
                k_max = max(k_max, int(row["k"]))
            label = csv_path.name
            L = int(label[1:3]) if label[1].isdigit() else -1
            rows.append(
                {
                    "file": label,
                    "L": L,
                    "n_k": n_rows,
                    "k_max": k_max,
                    "N": k_max,
                    "columns": ",".join(fields),
                }
            )
    return rows


def load_roots(issue9: Path) -> RootTables:
    pooled = {}
    with (issue9 / "root_sequence.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (int(row["L"]), row["observable"], row["mode"])
            pooled[key] = {
                "L": int(row["L"]),
                "observable": row["observable"],
                "mode": row["mode"],
                "pooled_root": to_mpf(row["pooled_root"]),
                "batch_sd": to_mpf(row["batch_sd"]),
                "batch_se": to_mpf(row["batch_se"]),
                "replicas": int(row["replicas"]),
                "runtime": row["runtime"],
                "status": row["status"],
            }
    batches = {}
    with (issue9 / "roots_by_batch.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (int(row["L"]), row["observable"], row["mode"], int(row["batch"]))
            batches[key] = to_mpf(row["root"])
    return RootTables(pooled=pooled, batches=batches, issue7=[])


def sizes_for_mode(mode: str) -> tuple[int, ...]:
    return SIZES_SHARED if mode == "shared" else SIZES_INDEPENDENT


def available_sizes(tables: RootTables, observable: str, mode: str) -> list[int]:
    out = []
    for L in sizes_for_mode(mode):
        rec = tables.pooled.get((L, observable, mode))
        if rec is not None and is_finite(rec["pooled_root"]):
            out.append(L)
    return out


def next_size(sizes: Sequence[int], current: int) -> int | None:
    larger = [s for s in sizes if s > current]
    return larger[0] if larger else None


def topology_of(observable: str, mode: str) -> str:
    family = "cluster_matching" if observable == "M" else f"wrapping_{observable}"
    return f"torus_{mode}:{family}"


class Microcanonical:
    def __init__(self, issue9: Path):
        self.issue9 = issue9
        self.cache: dict[tuple[int, str], dict] = {}

    def load(self, L: int, mode: str) -> dict | None:
        key = (L, mode)
        if key in self.cache:
            return self.cache[key]
        name = f"microcanonical_L{L:03d}.csv"
        if mode == "independent":
            name = f"microcanonical_L{L:03d}_independent.csv"
        path = self.issue9 / name
        if not path.is_file():
            self.cache[key] = None  # type: ignore[assignment]
            return None
        n, replicas, cl_g, cl_gs, wrap_g, wrap_gs, stats = nz.read_pooled_csv(str(path))
        rec = {
            "n": n,
            "replicas": replicas,
            "cl_g": cl_g,
            "cl_gs": cl_gs,
            "wrap_g": wrap_g,
            "wrap_gs": wrap_gs,
            "M": nz.make_M(cl_g, cl_gs, n),
            "D": {
                "H": nz.make_D_wrap(wrap_g["H"], wrap_gs["H"], n),
                "V": nz.make_D_wrap(wrap_g["V"], wrap_gs["V"], n),
                "either": nz.make_D_wrap(wrap_g["E"], wrap_gs["E"], n),
                "both": nz.make_D_wrap(wrap_g["B"], wrap_gs["B"], n),
            },
        }
        self.cache[key] = rec
        return rec

    def matching_fn(self, L: int, mode: str, observable: str):
        rec = self.load(L, mode)
        if rec is None:
            return None
        if observable == "M":
            return rec["M"]
        return rec["D"][observable]


def combine_batches(
    tables: RootTables,
    observable: str,
    mode: str,
    sizes: Sequence[int],
    weights: Sequence[mp.mpf],
) -> dict:
    supers = []
    for batch in range(N_BATCHES):
        acc = mp.mpf(0)
        ok = True
        for size, weight in zip(sizes, weights):
            value = tables.batches.get((size, observable, mode, batch))
            if value is None or not is_finite(value):
                ok = False
                break
            acc += weight * value
        if ok:
            supers.append(acc)
    mu, sd, se, n = mean_sd_se(supers)
    largest = sizes[-1]
    ref = tables.pooled.get((largest, observable, mode), {})
    ref_sd = ref.get("batch_sd", mp.mpf("nan"))
    if n < 2:
        sd = mp.mpf("nan")
        se = mp.mpf("nan")
        emp_amp = mp.mpf("nan")
    else:
        emp_amp = sd / ref_sd if is_finite(sd) and is_finite(ref_sd) and ref_sd != 0 else mp.mpf("nan")
    pooled_combo = mp.mpf(0)
    pooled_ok = True
    for size, weight in zip(sizes, weights):
        rec = tables.pooled.get((size, observable, mode))
        if rec is None or not is_finite(rec["pooled_root"]):
            pooled_ok = False
            break
        pooled_combo += weight * rec["pooled_root"]
    if not pooled_ok:
        pooled_combo = mp.mpf("nan")
    return {
        "combined_pooled": pooled_combo,
        "combined_batch_mean": mu,
        "batch_sd": sd,
        "batch_se": se,
        "n_batches": n,
        "empirical_noise_amp": emp_amp,
        "ref_L": largest,
        "ref_batch_sd": ref_sd,
    }


def ordinary_fits(tables: RootTables) -> list[dict]:
    rows = []
    for mode in ("shared", "independent"):
        size_grid = sizes_for_mode(mode)
        train_max_candidates = size_grid
        for observable in OBSERVABLES:
            avail = available_sizes(tables, observable, mode)
            for train_max in train_max_candidates:
                train_pool = [L for L in avail if L <= train_max]
                withheld_L = next_size(size_grid, train_max)
                withheld_rec = (
                    tables.pooled.get((withheld_L, observable, mode)) if withheld_L is not None else None
                )
                for L_min in train_pool:
                    train_sizes = [L for L in train_pool if L >= L_min]
                    train_vals = [tables.pooled[(L, observable, mode)]["pooled_root"] for L in train_sizes]
                    for model_name, powers in MODELS.items():
                        fit = fit_correction(train_sizes, train_vals, powers)
                        row = {
                            "observable": observable,
                            "mode": mode,
                            "topology": topology_of(observable, mode),
                            "model": model_name,
                            "powers": ",".join(map(str, powers)),
                            "L_min": L_min,
                            "train_max": train_max,
                            "n_train": len(train_sizes),
                            "train_sizes": " ".join(map(str, train_sizes)),
                            "pc": "",
                            "a4": "",
                            "a6": "",
                            "a8": "",
                            "condition_number": "",
                            "training_rmse": "",
                            "training_max_abs": "",
                            "withheld_L": withheld_L if withheld_L is not None else "",
                            "withheld_true": "",
                            "predicted": "",
                            "signed_error": "",
                            "abs_error": "",
                            "withheld_batch_se": "",
                            "standardized_error": "",
                            "status": fit.get("status", "ok"),
                        }
                        if fit.get("status") == "ok":
                            row.update(
                                {
                                    "pc": fmt(fit["pc"]),
                                    "a4": fmt(fit["a4"], AMP_DIGITS),
                                    "a6": fmt(fit["a6"], AMP_DIGITS),
                                    "a8": fmt(fit["a8"], AMP_DIGITS),
                                    "condition_number": fmt(fit["condition_number"], 12),
                                    "training_rmse": fmt(fit["training_rmse"]),
                                    "training_max_abs": fmt(fit["training_max_abs"]),
                                }
                            )
                            if withheld_rec is not None and is_finite(withheld_rec["pooled_root"]):
                                pred = predict_p(withheld_L, fit)
                                true = withheld_rec["pooled_root"]
                                err = pred - true
                                se = withheld_rec["batch_se"]
                                row["withheld_true"] = fmt(true)
                                row["predicted"] = fmt(pred)
                                row["signed_error"] = fmt(err)
                                row["abs_error"] = fmt(abs(err))
                                row["withheld_batch_se"] = fmt(se)
                                row["standardized_error"] = fmt(err / se) if is_finite(se) and se != 0 else "nan"
                        rows.append(row)
    return rows


def select_on_training(
    tables: RootTables,
    observable: str,
    mode: str,
    train_sizes: Sequence[int],
) -> dict:
    """Choose model and L_min using only training sizes (inner one-step holdout)."""
    if len(train_sizes) < 3:
        return {"status": "too_few_points"}
    inner_holdout = train_sizes[-1]
    inner_train_pool = list(train_sizes[:-1])
    candidates = []
    for L_min in inner_train_pool:
        inner_sizes = [L for L in inner_train_pool if L >= L_min]
        inner_vals = [tables.pooled[(L, observable, mode)]["pooled_root"] for L in inner_sizes]
        for model_name, powers in MODELS.items():
            if len(inner_sizes) < 1 + len(powers):
                continue
            fit = fit_correction(inner_sizes, inner_vals, powers)
            if fit.get("status") != "ok":
                continue
            pred = predict_p(inner_holdout, fit)
            true = tables.pooled[(inner_holdout, observable, mode)]["pooled_root"]
            if not is_finite(pred) or not is_finite(true):
                continue
            abs_err = abs(pred - true)
            complexity = len(powers)
            candidates.append(
                {
                    "model": model_name,
                    "powers": powers,
                    "L_min": L_min,
                    "inner_abs_error": abs_err,
                    "complexity": complexity,
                    "inner_n_train": len(inner_sizes),
                }
            )
    if not candidates:
        return {"status": "no_inner_candidate"}
    candidates.sort(key=lambda c: (c["inner_abs_error"], c["complexity"], -c["L_min"]))
    best = candidates[0]
    final_sizes = [L for L in train_sizes if L >= best["L_min"]]
    final_vals = [tables.pooled[(L, observable, mode)]["pooled_root"] for L in final_sizes]
    fit = fit_correction(final_sizes, final_vals, best["powers"])
    return {
        "status": fit.get("status", "ok"),
        "model": best["model"],
        "powers": best["powers"],
        "L_min": best["L_min"],
        "inner_abs_error": best["inner_abs_error"],
        "inner_holdout": inner_holdout,
        "n_candidates": len(candidates),
        "fit": fit,
        "train_sizes": final_sizes,
    }


def rolling_predictions(tables: RootTables) -> list[dict]:
    rows = []
    mode = "shared"
    size_grid = SIZES_SHARED
    for observable in OBSERVABLES:
        avail = available_sizes(tables, observable, mode)
        for train_max, withheld_L in ROLLING_FOLDS:
            train_sizes = [L for L in avail if L <= train_max]
            withheld = tables.pooled.get((withheld_L, observable, mode))
            selection = select_on_training(tables, observable, mode, train_sizes)
            row = {
                "observable": observable,
                "mode": mode,
                "train_max": train_max,
                "train_sizes": " ".join(map(str, train_sizes)),
                "withheld_L": withheld_L,
                "selected_model": selection.get("model", ""),
                "selected_L_min": selection.get("L_min", ""),
                "inner_holdout_L": selection.get("inner_holdout", ""),
                "inner_abs_error": fmt(selection["inner_abs_error"]) if "inner_abs_error" in selection else "",
                "n_inner_candidates": selection.get("n_candidates", ""),
                "pc": "",
                "a4": "",
                "a6": "",
                "a8": "",
                "condition_number": "",
                "training_rmse": "",
                "predicted": "",
                "true_value": "",
                "signed_error": "",
                "abs_error": "",
                "withheld_batch_se": "",
                "standardized_error": "",
                "status": selection.get("status", "fail"),
            }
            fit = selection.get("fit")
            if fit and fit.get("status") == "ok":
                row["pc"] = fmt(fit["pc"])
                row["a4"] = fmt(fit["a4"], AMP_DIGITS)
                row["a6"] = fmt(fit["a6"], AMP_DIGITS)
                row["a8"] = fmt(fit["a8"], AMP_DIGITS)
                row["condition_number"] = fmt(fit["condition_number"], 12)
                row["training_rmse"] = fmt(fit["training_rmse"])
                if withheld is not None and is_finite(withheld["pooled_root"]):
                    pred = predict_p(withheld_L, fit)
                    true = withheld["pooled_root"]
                    err = pred - true
                    se = withheld["batch_se"]
                    row["predicted"] = fmt(pred)
                    row["true_value"] = fmt(true)
                    row["signed_error"] = fmt(err)
                    row["abs_error"] = fmt(abs(err))
                    row["withheld_batch_se"] = fmt(se)
                    row["standardized_error"] = fmt(err / se) if is_finite(se) and se != 0 else "nan"
                    row["status"] = "ok"
            rows.append(row)
    return rows


def adjacent_pairs(sizes: Sequence[int]) -> list[tuple[int, int]]:
    return [(sizes[i], sizes[i + 1]) for i in range(len(sizes) - 1)]


def consecutive_triples(sizes: Sequence[int]) -> list[tuple[int, int, int]]:
    return [(sizes[i], sizes[i + 1], sizes[i + 2]) for i in range(len(sizes) - 2)]


def cancel_L4_rows(tables: RootTables) -> list[dict]:
    rows = []
    exp4 = [mp.mpf(4)]
    for mode in ("shared", "independent"):
        size_grid = sizes_for_mode(mode)
        for observable in OBSERVABLES:
            avail = available_sizes(tables, observable, mode)
            for L1, L2 in adjacent_pairs(avail):
                solved = solve_cancel_weights((L1, L2), exp4)
                weights = solved["weights"]
                combo = combine_batches(tables, observable, mode, (L1, L2), weights)
                withheld_L = next_size(size_grid, L2)
                withheld = tables.pooled.get((withheld_L, observable, mode)) if withheld_L else None
                row = {
                    "observable": observable,
                    "mode": mode,
                    "L1": L1,
                    "L2": L2,
                    "w1": fmt(weights[0], WEIGHT_DIGITS),
                    "w2": fmt(weights[1], WEIGHT_DIGITS),
                    "sum_abs_weights": fmt(solved["sum_abs_weights"], 18),
                    "l2_norm": fmt(solved["l2_norm"], 18),
                    "condition_number": fmt(solved["condition_number"], 12),
                    "combined_pooled": fmt(combo["combined_pooled"]),
                    "combined_batch_mean": fmt(combo["combined_batch_mean"]),
                    "batch_sd": fmt(combo["batch_sd"]),
                    "batch_se": fmt(combo["batch_se"]),
                    "n_batches": combo["n_batches"],
                    "noise_amp_l1": fmt(solved["sum_abs_weights"], 18),
                    "noise_amp_l2": fmt(solved["l2_norm"], 18),
                    "noise_amp_empirical": fmt(combo["empirical_noise_amp"], 18),
                    "ref_L": combo["ref_L"],
                    "ref_batch_sd": fmt(combo["ref_batch_sd"]),
                    "withheld_L": withheld_L if withheld_L is not None else "",
                    "withheld_true": "",
                    "signed_error": "",
                    "abs_error": "",
                    "standardized_error": "",
                    "status": "ok",
                }
                if withheld is not None and is_finite(withheld["pooled_root"]) and is_finite(combo["combined_pooled"]):
                    true = withheld["pooled_root"]
                    err = combo["combined_pooled"] - true
                    se = withheld["batch_se"]
                    row["withheld_true"] = fmt(true)
                    row["signed_error"] = fmt(err)
                    row["abs_error"] = fmt(abs(err))
                    row["standardized_error"] = fmt(err / se) if is_finite(se) and se != 0 else "nan"
                rows.append(row)
    return rows


def cancel_L4_L6_rows(tables: RootTables) -> list[dict]:
    rows = []
    exponents = [mp.mpf(4), mp.mpf(6)]
    for mode in ("shared", "independent"):
        size_grid = sizes_for_mode(mode)
        for observable in OBSERVABLES:
            avail = available_sizes(tables, observable, mode)
            for L1, L2, L3 in consecutive_triples(avail):
                solved = solve_cancel_weights((L1, L2, L3), exponents)
                weights = solved["weights"]
                combo = combine_batches(tables, observable, mode, (L1, L2, L3), weights)
                withheld_L = next_size(size_grid, L3)
                withheld = tables.pooled.get((withheld_L, observable, mode)) if withheld_L else None
                row = {
                    "observable": observable,
                    "mode": mode,
                    "L1": L1,
                    "L2": L2,
                    "L3": L3,
                    "w1": fmt(weights[0], WEIGHT_DIGITS),
                    "w2": fmt(weights[1], WEIGHT_DIGITS),
                    "w3": fmt(weights[2], WEIGHT_DIGITS),
                    "sum_abs_weights": fmt(solved["sum_abs_weights"], 18),
                    "l2_norm": fmt(solved["l2_norm"], 18),
                    "condition_number": fmt(solved["condition_number"], 12),
                    "combined_pooled": fmt(combo["combined_pooled"]),
                    "combined_batch_mean": fmt(combo["combined_batch_mean"]),
                    "batch_sd": fmt(combo["batch_sd"]),
                    "batch_se": fmt(combo["batch_se"]),
                    "n_batches": combo["n_batches"],
                    "noise_amp_l1": fmt(solved["sum_abs_weights"], 18),
                    "noise_amp_l2": fmt(solved["l2_norm"], 18),
                    "noise_amp_empirical": fmt(combo["empirical_noise_amp"], 18),
                    "ref_L": combo["ref_L"],
                    "ref_batch_sd": fmt(combo["ref_batch_sd"]),
                    "withheld_L": withheld_L if withheld_L is not None else "",
                    "withheld_true": "",
                    "signed_error": "",
                    "abs_error": "",
                    "standardized_error": "",
                    "status": "ok",
                }
                if withheld is not None and is_finite(withheld["pooled_root"]) and is_finite(combo["combined_pooled"]):
                    true = withheld["pooled_root"]
                    err = combo["combined_pooled"] - true
                    se = withheld["batch_se"]
                    row["withheld_true"] = fmt(true)
                    row["signed_error"] = fmt(err)
                    row["abs_error"] = fmt(abs(err))
                    row["standardized_error"] = fmt(err / se) if is_finite(se) and se != 0 else "nan"
                rows.append(row)
    return rows


def matching_annihilator_rows(tables: RootTables, micro: Microcanonical) -> list[dict]:
    rows = []
    exponent = mp.mpf(13) / mp.mpf(4)
    for mode in ("shared", "independent"):
        size_grid = sizes_for_mode(mode)
        for observable in OBSERVABLES:
            avail = available_sizes(tables, observable, mode)
            for L1, L2 in adjacent_pairs(avail):
                f1 = micro.matching_fn(L1, mode, observable)
                f2 = micro.matching_fn(L2, mode, observable)
                row = {
                    "observable": observable,
                    "mode": mode,
                    "L1": L1,
                    "L2": L2,
                    "prefactor_L1": fmt(mp.power(L1, exponent), 18),
                    "prefactor_L2": fmt(mp.power(L2, exponent), 18),
                    "annihilator_root": "",
                    "annihilator_status": "",
                    "ordinary_root_L1": fmt(tables.pooled[(L1, observable, mode)]["pooled_root"]),
                    "ordinary_root_L2": fmt(tables.pooled[(L2, observable, mode)]["pooled_root"]),
                    "ordinary_L1_status": tables.pooled[(L1, observable, mode)]["status"],
                    "ordinary_L2_status": tables.pooled[(L2, observable, mode)]["status"],
                    "withheld_L": "",
                    "withheld_ordinary_root": "",
                    "signed_error": "",
                    "abs_error": "",
                    "standardized_error": "",
                    "weight_l1": "",
                    "weight_l2": "",
                    "sum_abs_weights": "",
                    "condition_number": "",
                    "status": "no_microcanonical",
                }
                wsol = solve_cancel_weights((L1, L2), [exponent])
                row["weight_l1"] = fmt(wsol["weights"][0], WEIGHT_DIGITS)
                row["weight_l2"] = fmt(wsol["weights"][1], WEIGHT_DIGITS)
                row["sum_abs_weights"] = fmt(wsol["sum_abs_weights"], 18)
                row["condition_number"] = fmt(wsol["condition_number"], 12)
                withheld_L = next_size(size_grid, L2)
                if withheld_L is not None:
                    row["withheld_L"] = withheld_L
                    wrec = tables.pooled.get((withheld_L, observable, mode))
                    if wrec is not None:
                        row["withheld_ordinary_root"] = fmt(wrec["pooled_root"])
                if f1 is None or f2 is None:
                    rows.append(row)
                    continue
                pref1 = float(L1 ** (13.0 / 4.0))
                pref2 = float(L2 ** (13.0 / 4.0))

                def A(p: float, a=f1, b=f2, c=pref1, d=pref2) -> float:
                    return c * a(p) - d * b(p)

                root, status, _ival = nz.find_root(A)
                row["annihilator_root"] = nz.fmt18(root)
                row["annihilator_status"] = status
                row["status"] = status
                if withheld_L is not None:
                    wrec = tables.pooled.get((withheld_L, observable, mode))
                    if wrec is not None and is_finite(wrec["pooled_root"]) and root == root:
                        true = wrec["pooled_root"]
                        err = to_mpf(root) - true
                        se = wrec["batch_se"]
                        row["signed_error"] = fmt(err)
                        row["abs_error"] = fmt(abs(err))
                        row["standardized_error"] = fmt(err / se) if is_finite(se) and se != 0 else "nan"
                rows.append(row)
    return rows


def amplitude_sign_rows(ordinary: list[dict]) -> list[dict]:
    rows = []
    for row in ordinary:
        if row["status"] != "ok":
            continue
        if row["model"] not in ("B", "C"):
            continue
        a4 = to_mpf(row["a4"])
        a6 = to_mpf(row["a6"])
        a8 = to_mpf(row["a8"]) if row["model"] == "C" else mp.mpf("nan")
        ratio = a4 / a6 if is_finite(a6) and a6 != 0 else mp.mpf("nan")
        rows.append(
            {
                "observable": row["observable"],
                "mode": row["mode"],
                "topology": row["topology"],
                "model": row["model"],
                "L_min": row["L_min"],
                "train_max": row["train_max"],
                "n_train": row["n_train"],
                "a4": row["a4"],
                "a6": row["a6"],
                "a8": row["a8"],
                "a4_sign": sign_label(a4),
                "a6_sign": sign_label(a6),
                "a8_sign": sign_label(a8) if row["model"] == "C" else "",
                "a4_over_a6": fmt(ratio, 18),
                "condition_number": row["condition_number"],
                "training_rmse": row["training_rmse"],
            }
        )
    return rows


def noise_amplification_rows(c4: list[dict], c46: list[dict], ann: list[dict]) -> list[dict]:
    rows = []
    for row in c4:
        rows.append(
            {
                "method": "L4_cancel",
                "observable": row["observable"],
                "mode": row["mode"],
                "sizes": f"{row['L1']} {row['L2']}",
                "sum_abs_weights": row["sum_abs_weights"],
                "l2_norm": row["l2_norm"],
                "condition_number": row["condition_number"],
                "batch_sd_combined": row["batch_sd"],
                "batch_se_combined": row["batch_se"],
                "ref_L": row["ref_L"],
                "ref_batch_sd": row["ref_batch_sd"],
                "empirical_noise_amp": row["noise_amp_empirical"],
                "n_batches": row["n_batches"],
                "note": "empirical amp = batch_sd(sum w_i p_Li[b]) / batch_sd(p_Lmax); sizes independent across L, matching covariance is inside each p_L",
            }
        )
    for row in c46:
        rows.append(
            {
                "method": "L4_L6_cancel",
                "observable": row["observable"],
                "mode": row["mode"],
                "sizes": f"{row['L1']} {row['L2']} {row['L3']}",
                "sum_abs_weights": row["sum_abs_weights"],
                "l2_norm": row["l2_norm"],
                "condition_number": row["condition_number"],
                "batch_sd_combined": row["batch_sd"],
                "batch_se_combined": row["batch_se"],
                "ref_L": row["ref_L"],
                "ref_batch_sd": row["ref_batch_sd"],
                "empirical_noise_amp": row["noise_amp_empirical"],
                "n_batches": row["n_batches"],
                "note": "empirical amp = batch_sd(sum w_i p_Li[b]) / batch_sd(p_Lmax)",
            }
        )
    for row in ann:
        rows.append(
            {
                "method": "matching_annihilator_13_4",
                "observable": row["observable"],
                "mode": row["mode"],
                "sizes": f"{row['L1']} {row['L2']}",
                "sum_abs_weights": row["sum_abs_weights"],
                "l2_norm": "",
                "condition_number": row["condition_number"],
                "batch_sd_combined": "",
                "batch_se_combined": "",
                "ref_L": row["L2"],
                "ref_batch_sd": "",
                "empirical_noise_amp": "",
                "n_batches": "",
                "note": "weights annihilate L^{-13/4} in M_L or D_L; per-batch annihilator SD unavailable (pooled microcanonical only)",
            }
        )
    return rows


def parse_abs(row: dict, key: str = "abs_error") -> mp.mpf | None:
    text = row.get(key, "")
    if text in ("", None, "nan"):
        return None
    value = to_mpf(text)
    return value if is_finite(value) else None


def median_worst(rows: list[dict], key: str = "abs_error") -> tuple[int, mp.mpf, mp.mpf]:
    vals = [parse_abs(r, key) for r in rows]
    vals = [v for v in vals if v is not None]
    if not vals:
        nan = mp.mpf("nan")
        return 0, nan, nan
    return len(vals), mp.mpf(str(median(vals))), max(vals)


def verdict(ratio: mp.mpf) -> str:
    if not is_finite(ratio):
        return "NA"
    if ratio < 1:
        return "BETTER"
    if ratio > 1:
        return "WORSE"
    return "SAME"


def rolling_aligned_cancel(rows: list[dict], L_fields: tuple[str, ...]) -> list[dict]:
    """Keep the cancellation that uses the largest training sizes for each rolling fold."""
    wanted = []
    last_field = L_fields[-1]
    for train_max, withheld_L in ROLLING_FOLDS:
        wanted.append((train_max, withheld_L))
    out = []
    for row in rows:
        if row.get("mode") != "shared":
            continue
        withheld = row.get("withheld_L")
        last = int(row[last_field])
        if withheld in ("", None):
            continue
        withheld_i = int(withheld)
        # aligned iff last participating size equals the rolling train_max
        if (last, withheld_i) in [(t, w) for t, w in ROLLING_FOLDS]:
            out.append(row)
    return out


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def run_cmd(args: list[str]) -> str:
    try:
        proc = subprocess.run(args, check=False, capture_output=True, text=True)
        return (proc.stdout or proc.stderr or "").rstrip()
    except OSError as exc:
        return f"(unavailable: {exc})"


def write_environment(path: Path, repo: Path, issue7_meta: list[dict]) -> None:
    pip_freeze = run_cmd([sys.executable, "-m", "pip", "freeze"])
    git_head = run_cmd(["git", "-C", str(repo), "rev-parse", "HEAD"])
    git_branch = run_cmd(["git", "-C", str(repo), "rev-parse", "--abbrev-ref", "HEAD"])
    git_status = run_cmd(["git", "-C", str(repo), "status", "--short"])
    lines = [
        run_cmd(["date", "-u"]),
        run_cmd(["uname", "-a"]),
        "",
        run_cmd(["lscpu"]),
        "",
        run_cmd(["free", "-h"]),
        "",
        f"python: {sys.version.replace(chr(10), ' ')}",
        f"platform: {platform.platform()}",
        f"mpmath: {mp.__version__}",
        f"mp.dps: {mp.mp.dps}",
        f"git_branch: {git_branch}",
        f"git_head: {git_head}",
        f"git_status:",
        git_status,
        "",
        "pip freeze:",
        pip_freeze,
        "",
        "issue-7 microcanonical files read:",
    ]
    for rec in issue7_meta:
        lines.append(f"  {rec['file']}: L={rec['L']} n_k={rec['n_k']} k_max={rec['k_max']} columns={rec['columns']}")
    lines.append("")
    lines.append("issue-9 inputs: root_sequence.csv, roots_by_batch.csv, microcanonical_L*.csv")
    lines.append("no Monte Carlo was run for Issue #8")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def md_fmt(value, digits: int = 8) -> str:
    if not is_finite(value):
        return "n/a"
    return mp.nstr(to_mpf(value), n=digits, strip_zeros=False)


def write_report(
    path: Path,
    tables: RootTables,
    issue7_meta: list[dict],
    ordinary: list[dict],
    rolling: list[dict],
    c4: list[dict],
    c46: list[dict],
    ann: list[dict],
    amps: list[dict],
) -> None:
    a: list[str] = []
    a.append("# Issue #8 — L^{-4} / L^{-6} cancellation in matched estimators")
    a.append("")
    a.append("Post-processing of Issue #7 exact torus microcanonical totals and Issue #9")
    a.append("Newman–Ziff matched roots. **No Monte Carlo was run.** Weights for")
    a.append("Richardson cancellation depend only on lattice sizes. Training residuals")
    a.append("are fit diagnostics and are **not** uncertainties. Intercept min–max across")
    a.append("models is **model spread**, not a statistical confidence interval.")
    a.append("Estimator ranking uses withheld-size prediction error only, never closeness")
    a.append("to a published threshold.")
    a.append("")
    a.append("## Inputs")
    a.append("")
    a.append("Issue #7 exact microcanonical files:")
    a.append("")
    for rec in issue7_meta:
        a.append(f"- `{rec['file']}`: L={rec['L']}, occupancy k=0..{rec['k_max']}")
    a.append("")
    a.append("Issue #9 shared roots on L = 16, 24, 32, 48, 64, 96, 128, 192, 256")
    a.append("(independent baseline on L = 32, 64, 128, 256). Observables: wrapping")
    a.append("H, V, either, both, and cluster matching M. Shared-mode wrapping roots")
    a.append("coincide to printing precision because of the Issue #7 wrapping identity.")
    a.append("Cluster matching M is noisier; some large-L M roots left `[0.590, 0.595]`")
    a.append("and are used as numeric roots, not discarded.")
    a.append("")
    a.append("## Protocol")
    a.append("")
    a.append("Rolling one-step folds (shared mode):")
    a.append("")
    a.append("- train through 64 → predict 96")
    a.append("- train through 96 → predict 128")
    a.append("- train through 128 → predict 192")
    a.append("- train through 192 → predict 256")
    a.append("")
    a.append("On each fold the ordinary correction **model and L_min are chosen only on")
    a.append("training sizes** (inner one-step holdout of the last training size; ties")
    a.append("break toward the simpler model, then larger L_min). The chosen pair is")
    a.append("frozen, refit on all training sizes, then scored on the withheld size.")
    a.append("")
    a.append("L^{-4} cancellation uses adjacent sizes with `w1+w2=1` and")
    a.append("`w1 L1^{-4}+w2 L2^{-4}=0`. L^{-4}/L^{-6} cancellation uses three consecutive")
    a.append("sizes and also annihilates L^{-6}. Combined uncertainty is the SD/SE of")
    a.append("the 20 batch combinations `p_super[b]=sum_i w_i p_{L_i}[b]`.")
    a.append("")
    a.append("Matching annihilator (when microcanonical statistics exist):")
    a.append("")
    a.append("```")
    a.append("A(p) = L1^{13/4} F_{L1}(p) - L2^{13/4} F_{L2}(p) = 0")
    a.append("```")
    a.append("")
    a.append("with `F = M_L` for cluster matching and `F = D_L^x` for wrapping class x.")
    a.append("The annihilator root is scored against the ordinary withheld-size root.")
    a.append("Per-batch annihilator SD is not available: Issue #9 stored pooled")
    a.append("microcanonical sufficient statistics, not per-batch histograms.")
    a.append("")
    a.append("Marks BETTER / SAME / WORSE compare **median absolute withheld error**")
    a.append("to the ordinary rolling estimator on the same observable. They describe")
    a.append("out-of-sample numbers only.")
    a.append("")

    c4_aligned = rolling_aligned_cancel(c4, ("L1", "L2"))
    c46_aligned = rolling_aligned_cancel(c46, ("L1", "L2", "L3"))
    ann_aligned = rolling_aligned_cancel(ann, ("L1", "L2"))

    a.append("## Per-observable withheld-size comparison (shared, rolling-aligned)")
    a.append("")
    a.append("Ordinary = nested model/L_min selection predicting `p_{L_next}`.")
    a.append("L^{-4}-cancel = last adjacent training pair, scored vs `p_{L_next}`.")
    a.append("L^{-4}/L^{-6}-cancel = last consecutive training triple, scored vs `p_{L_next}`.")
    a.append("Matching annihilator = `A(p)=0` on last adjacent training pair, scored vs ordinary `p_{L_next}`.")
    a.append("")

    for observable in OBSERVABLES:
        a.append(f"### Observable `{observable}`")
        a.append("")
        ord_rows = [r for r in rolling if r["observable"] == observable]
        c4_rows = [r for r in c4_aligned if r["observable"] == observable]
        c46_rows = [r for r in c46_aligned if r["observable"] == observable]
        ann_rows = [r for r in ann_aligned if r["observable"] == observable]

        n_o, med_o, worst_o = median_worst(ord_rows)
        n4, med4, worst4 = median_worst(c4_rows)
        n46, med46, worst46 = median_worst(c46_rows)
        na, meda, worsta = median_worst(ann_rows)

        ratio4 = med4 / med_o if is_finite(med4) and is_finite(med_o) and med_o != 0 else mp.mpf("nan")
        ratio46 = med46 / med_o if is_finite(med46) and is_finite(med_o) and med_o != 0 else mp.mpf("nan")
        ratioa = meda / med_o if is_finite(meda) and is_finite(med_o) and med_o != 0 else mp.mpf("nan")

        def amp_stats(rows, key="noise_amp_empirical"):
            vals = [to_mpf(r.get(key, "nan")) for r in rows]
            vals = [v for v in vals if is_finite(v)]
            if not vals:
                return mp.mpf("nan"), mp.mpf("nan")
            return mp.mpf(str(median(vals))), max(vals)

        emp4_med, emp4_max = amp_stats(c4_rows)
        l14_med, l14_max = amp_stats(c4_rows, "noise_amp_l1")
        emp46_med, emp46_max = amp_stats(c46_rows)
        l146_med, l146_max = amp_stats(c46_rows, "sum_abs_weights")
        l1a_med, l1a_max = amp_stats(ann_rows, "sum_abs_weights")

        a.append("| estimator | n_success | median_abs_err | worst_abs_err | ratio vs ordinary | mark | noise amp (emp median / L1 median) |")
        a.append("|---|---:|---:|---:|---:|---|---|")
        a.append(
            f"| ordinary rolling | {n_o} | {md_fmt(med_o)} | {md_fmt(worst_o)} | 1 | — | n/a |"
        )
        a.append(
            f"| L^{-4} cancel | {n4} | {md_fmt(med4)} | {md_fmt(worst4)} | {md_fmt(ratio4, 6)} | {verdict(ratio4)} | emp {md_fmt(emp4_med, 6)} / L1 {md_fmt(l14_med, 6)} |"
        )
        a.append(
            f"| L^{-4}/L^{-6} cancel | {n46} | {md_fmt(med46)} | {md_fmt(worst46)} | {md_fmt(ratio46, 6)} | {verdict(ratio46)} | emp {md_fmt(emp46_med, 6)} / L1 {md_fmt(l146_med, 6)} |"
        )
        a.append(
            f"| matching annihilator | {na} | {md_fmt(meda)} | {md_fmt(worsta)} | {md_fmt(ratioa, 6)} | {verdict(ratioa)} | L1 {md_fmt(l1a_med, 6)} (no batch SD) |"
        )
        a.append("")
        a.append("Fold-level ordinary rolling:")
        a.append("")
        a.append("| train_max | withheld_L | model | L_min | predicted | true | signed err | abs err | std err (err/SE) |")
        a.append("|---:|---:|---|---:|---|---|---|---|---|")
        for r in ord_rows:
            a.append(
                f"| {r['train_max']} | {r['withheld_L']} | {r['selected_model']} | {r['selected_L_min']} | "
                f"{r['predicted'] or 'n/a'} | {r['true_value'] or 'n/a'} | {r['signed_error'] or 'n/a'} | "
                f"{r.get('abs_error') or 'n/a'} | {r.get('standardized_error') or 'n/a'} |"
            )
        a.append("")
        a.append("Fold-level L^{-4} cancel (last training pair):")
        a.append("")
        a.append("| L1 | L2 | withheld_L | combined | true | abs err | emp noise amp | L1 amp |")
        a.append("|---:|---:|---:|---|---|---|---|---|")
        for r in c4_rows:
            a.append(
                f"| {r['L1']} | {r['L2']} | {r['withheld_L']} | {r['combined_pooled']} | "
                f"{r.get('withheld_true') or 'n/a'} | {r.get('abs_error') or 'n/a'} | "
                f"{r['noise_amp_empirical']} | {r['noise_amp_l1']} |"
            )
        a.append("")
        a.append("Fold-level L^{-4}/L^{-6} cancel (last training triple):")
        a.append("")
        a.append("| L1 | L2 | L3 | withheld_L | combined | true | abs err | emp noise amp | sum\\|w\\| | cond |")
        a.append("|---:|---:|---:|---:|---|---|---|---|---|---|")
        for r in c46_rows:
            a.append(
                f"| {r['L1']} | {r['L2']} | {r['L3']} | {r['withheld_L']} | {r['combined_pooled']} | "
                f"{r.get('withheld_true') or 'n/a'} | {r.get('abs_error') or 'n/a'} | "
                f"{r['noise_amp_empirical']} | {r['sum_abs_weights']} | {r['condition_number']} |"
            )
        a.append("")
        a.append("Fold-level matching annihilator:")
        a.append("")
        a.append("| L1 | L2 | withheld_L | annihilator root | ordinary withheld | abs err | status | sum\\|w\\| |")
        a.append("|---:|---:|---:|---|---|---|---|---|")
        for r in ann_rows:
            a.append(
                f"| {r['L1']} | {r['L2']} | {r['withheld_L']} | {r['annihilator_root'] or 'n/a'} | "
                f"{r.get('withheld_ordinary_root') or 'n/a'} | {r.get('abs_error') or 'n/a'} | "
                f"{r['status']} | {r['sum_abs_weights']} |"
            )
        a.append("")

    a.append("## All adjacent-pair / triple constructions (shared), not just rolling-aligned")
    a.append("")
    a.append("These rows use every adjacent pair or consecutive triple and score against")
    a.append("the next larger size. They are extra diagnostics; the marks above use only")
    a.append("the rolling-aligned subset.")
    a.append("")
    a.append("| observable | method | n_success | median_abs_err | worst_abs_err | median emp noise amp | median L1 amp |")
    a.append("|---|---|---:|---:|---:|---:|---:|")
    for observable in OBSERVABLES:
        for label, subset, emp_key, l1_key in (
            ("L4 all pairs", [r for r in c4 if r["mode"] == "shared" and r["observable"] == observable], "noise_amp_empirical", "noise_amp_l1"),
            ("L4L6 all triples", [r for r in c46 if r["mode"] == "shared" and r["observable"] == observable], "noise_amp_empirical", "sum_abs_weights"),
            ("annihilator all pairs", [r for r in ann if r["mode"] == "shared" and r["observable"] == observable], "", "sum_abs_weights"),
        ):
            n, med, worst = median_worst(subset)
            emp_vals = [to_mpf(r.get(emp_key, "nan")) for r in subset] if emp_key else []
            emp_vals = [v for v in emp_vals if is_finite(v)]
            l1_vals = [to_mpf(r.get(l1_key, "nan")) for r in subset]
            l1_vals = [v for v in l1_vals if is_finite(v)]
            emp_med = mp.mpf(str(median(emp_vals))) if emp_vals else mp.mpf("nan")
            l1_med = mp.mpf(str(median(l1_vals))) if l1_vals else mp.mpf("nan")
            a.append(
                f"| {observable} | {label} | {n} | {md_fmt(med)} | {md_fmt(worst)} | {md_fmt(emp_med, 6)} | {md_fmt(l1_med, 6)} |"
            )
    a.append("")

    a.append("## Ordinary full-window intercepts (model spread, not a CI)")
    a.append("")
    a.append("Fits on all shared sizes with L_min = 16. Training RMSE is a residual,")
    a.append("not an uncertainty.")
    a.append("")
    a.append("| observable | model | pc | a4 | a6 | a8 | cond | train RMSE | withheld L=none |")
    a.append("|---|---|---|---|---|---|---|---|---|")
    intercepts: dict[str, list[mp.mpf]] = defaultdict(list)
    for row in ordinary:
        if row["mode"] != "shared":
            continue
        if int(row["train_max"]) != 256 or int(row["L_min"]) != 16:
            continue
        if row["status"] != "ok":
            continue
        intercepts[row["observable"]].append(to_mpf(row["pc"]))
        a.append(
            f"| {row['observable']} | {row['model']} | {row['pc']} | {row['a4']} | {row['a6']} | {row['a8']} | "
            f"{row['condition_number']} | {row['training_rmse']} | full data |"
        )
    a.append("")
    a.append("Model spread of those intercepts (max − min over A/B/C at L_min=16, train_max=256):")
    a.append("")
    a.append("| observable | n_models | min pc | max pc | spread (max-min) |")
    a.append("|---|---:|---|---|---|")
    for observable in OBSERVABLES:
        vals = intercepts.get(observable, [])
        if not vals:
            a.append(f"| {observable} | 0 | n/a | n/a | n/a |")
            continue
        a.append(
            f"| {observable} | {len(vals)} | {fmt(min(vals))} | {fmt(max(vals))} | {fmt(max(vals) - min(vals))} |"
        )
    a.append("")
    a.append("That spread is **not** a statistical confidence interval.")
    a.append("")

    a.append("## Amplitude sign test (Model B/C)")
    a.append("")
    a.append("Signs of fitted `a4` and `a6` on every successful training window.")
    a.append("Cancellation is only meaningful if the signed amplitudes are stable")
    a.append("across observable / window / topology rather than accidental averages.")
    a.append("")

    def sign_summary(subset: list[dict], model: str) -> str:
        rel = [r for r in subset if r["model"] == model]
        if not rel:
            return "no fits"
        a4s = [r["a4_sign"] for r in rel]
        a6s = [r["a6_sign"] for r in rel]
        ratios = [to_mpf(r["a4_over_a6"]) for r in rel if is_finite(to_mpf(r["a4_over_a6"]))]
        return (
            f"n={len(rel)}; a4 signs { {s: a4s.count(s) for s in sorted(set(a4s))} }; "
            f"a6 signs { {s: a6s.count(s) for s in sorted(set(a6s))} }; "
            f"median a4/a6={md_fmt(mp.mpf(str(median(ratios))), 6) if ratios else 'n/a'}"
        )

    a.append("| observable | mode | Model B | Model C |")
    a.append("|---|---|---|---|")
    for observable in OBSERVABLES:
        for mode in ("shared", "independent"):
            sub = [r for r in amps if r["observable"] == observable and r["mode"] == mode]
            if not sub:
                continue
            a.append(f"| {observable} | {mode} | {sign_summary(sub, 'B')} | {sign_summary(sub, 'C')} |")
    a.append("")

    a.append("## Noise amplification")
    a.append("")
    a.append("L1 = sum |w_i| (fully correlated bound). L2 = sqrt(sum w_i^2) (independent")
    a.append("Gaussian bound). Empirical = batch SD of the weighted combination divided")
    a.append("by batch SD of the largest participating size. Because different L are")
    a.append("independent campaigns, matching covariance lives inside each p_L and the")
    a.append("empirical factor should track L2 more closely than L1.")
    a.append("")
    a.append("| method | observable | median L1 | median L2 | median empirical | max empirical | median cond |")
    a.append("|---|---|---:|---:|---:|---:|---:|")
    for method, subset, l1k, l2k, empk, condk in (
        ("L4_cancel", [r for r in c4 if r["mode"] == "shared"], "sum_abs_weights", "l2_norm", "noise_amp_empirical", "condition_number"),
        ("L4_L6_cancel", [r for r in c46 if r["mode"] == "shared"], "sum_abs_weights", "l2_norm", "noise_amp_empirical", "condition_number"),
        ("annihilator_13/4", [r for r in ann if r["mode"] == "shared"], "sum_abs_weights", "", "", "condition_number"),
    ):
        for observable in OBSERVABLES:
            rel = [r for r in subset if r["observable"] == observable]
            def medkey(k):
                if not k:
                    return mp.mpf("nan")
                vals = [to_mpf(r.get(k, "nan")) for r in rel]
                vals = [v for v in vals if is_finite(v)]
                return mp.mpf(str(median(vals))) if vals else mp.mpf("nan")
            def maxkey(k):
                if not k:
                    return mp.mpf("nan")
                vals = [to_mpf(r.get(k, "nan")) for r in rel]
                vals = [v for v in vals if is_finite(v)]
                return max(vals) if vals else mp.mpf("nan")
            a.append(
                f"| {method} | {observable} | {md_fmt(medkey(l1k), 6)} | {md_fmt(medkey(l2k), 6)} | "
                f"{md_fmt(medkey(empk), 6)} | {md_fmt(maxkey(empk), 6)} | {md_fmt(medkey(condk), 6)} |"
            )
    a.append("")
    a.append("## Notes")
    a.append("")
    a.append("- Shared wrapping H/V/either/both roots are identical at printing precision,")
    a.append("  so their ordinary/cancellation tables repeat. Cluster matching M does not.")
    a.append("- Independent-mode sizes are 32, 64, 128, 256 only; they appear in the CSVs")
    a.append("  but not in the specified rolling folds.")
    a.append("- Issue #7 L=2..5 totals were read and are the calibration source of Issue #9;")
    a.append("  they are not used as rolling-fit points.")
    a.append("- Fitted `a4` and `a6` signs flip across training windows for wrapping")
    a.append("  observables. The implied correction `a4 L^{-4}` is smaller than the")
    a.append("  batch SE of `p_L` (~1e-4), so the signed amplitudes are noise-dominated")
    a.append("  rather than a stable cancellation pattern. Cluster matching M is noisier")
    a.append("  still; several annihilator roots left `[0.590, 0.595]`.")
    a.append("- Batch combinations skip a batch if any participating size has a non-finite")
    a.append("  root. Empirical noise amp is undefined when fewer than two complete batches")
    a.append("  remain (reported as nan, not 0).")
    a.append("- No hardware or resource recommendation is made.")
    a.append("")
    path.write_text("\n".join(a) + "\n", encoding="utf-8")


ORDINARY_FIELDS = [
    "observable", "mode", "topology", "model", "powers", "L_min", "train_max", "n_train",
    "train_sizes", "pc", "a4", "a6", "a8", "condition_number", "training_rmse",
    "training_max_abs", "withheld_L", "withheld_true", "predicted", "signed_error",
    "abs_error", "withheld_batch_se", "standardized_error", "status",
]
ROLLING_FIELDS = [
    "observable", "mode", "train_max", "train_sizes", "withheld_L", "selected_model",
    "selected_L_min", "inner_holdout_L", "inner_abs_error", "n_inner_candidates",
    "pc", "a4", "a6", "a8", "condition_number", "training_rmse", "predicted",
    "true_value", "signed_error", "abs_error", "withheld_batch_se", "standardized_error",
    "status",
]
C4_FIELDS = [
    "observable", "mode", "L1", "L2", "w1", "w2", "sum_abs_weights", "l2_norm",
    "condition_number", "combined_pooled", "combined_batch_mean", "batch_sd", "batch_se",
    "n_batches", "noise_amp_l1", "noise_amp_l2", "noise_amp_empirical", "ref_L",
    "ref_batch_sd", "withheld_L", "withheld_true", "signed_error", "abs_error",
    "standardized_error", "status",
]
C46_FIELDS = [
    "observable", "mode", "L1", "L2", "L3", "w1", "w2", "w3", "sum_abs_weights",
    "l2_norm", "condition_number", "combined_pooled", "combined_batch_mean", "batch_sd",
    "batch_se", "n_batches", "noise_amp_l1", "noise_amp_l2", "noise_amp_empirical",
    "ref_L", "ref_batch_sd", "withheld_L", "withheld_true", "signed_error", "abs_error",
    "standardized_error", "status",
]
ANN_FIELDS = [
    "observable", "mode", "L1", "L2", "prefactor_L1", "prefactor_L2", "annihilator_root",
    "annihilator_status", "ordinary_root_L1", "ordinary_root_L2", "ordinary_L1_status",
    "ordinary_L2_status", "withheld_L", "withheld_ordinary_root", "signed_error",
    "abs_error", "standardized_error", "weight_l1", "weight_l2", "sum_abs_weights",
    "condition_number", "status",
]
AMP_FIELDS = [
    "observable", "mode", "topology", "model", "L_min", "train_max", "n_train",
    "a4", "a6", "a8", "a4_sign", "a6_sign", "a8_sign", "a4_over_a6",
    "condition_number", "training_rmse",
]
NOISE_FIELDS = [
    "method", "observable", "mode", "sizes", "sum_abs_weights", "l2_norm",
    "condition_number", "batch_sd_combined", "batch_se_combined", "ref_L",
    "ref_batch_sd", "empirical_noise_amp", "n_batches", "note",
]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--issue7", type=Path, default=REPO / "results" / "issue-7")
    parser.add_argument("--issue9", type=Path, default=REPO / "results" / "issue-9")
    parser.add_argument("--outdir", type=Path, default=REPO / "results" / "issue-8")
    parser.add_argument("--dps", type=int, default=80)
    args = parser.parse_args()
    if args.dps < 40:
        raise SystemExit("dps must be at least 40")
    mp.mp.dps = args.dps
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    issue7_meta = load_issue7(args.issue7)
    if not issue7_meta:
        raise SystemExit(f"no Issue #7 microcanonical CSVs in {args.issue7}")
    tables = load_roots(args.issue9)
    tables.issue7 = issue7_meta
    micro = Microcanonical(args.issue9)

    ordinary = ordinary_fits(tables)
    rolling = rolling_predictions(tables)
    c4 = cancel_L4_rows(tables)
    c46 = cancel_L4_L6_rows(tables)
    print("computing matching annihilators from microcanonical statistics...", file=sys.stderr)
    ann = matching_annihilator_rows(tables, micro)
    amps = amplitude_sign_rows(ordinary)
    noise = noise_amplification_rows(c4, c46, ann)

    write_csv(outdir / "ordinary_fit.csv", ordinary, ORDINARY_FIELDS)
    write_csv(outdir / "rolling_predictions.csv", rolling, ROLLING_FIELDS)
    write_csv(outdir / "cancel_L4.csv", c4, C4_FIELDS)
    write_csv(outdir / "cancel_L4_L6.csv", c46, C46_FIELDS)
    write_csv(outdir / "matching_annihilator.csv", ann, ANN_FIELDS)
    write_csv(outdir / "amplitude_signs.csv", amps, AMP_FIELDS)
    write_csv(outdir / "noise_amplification.csv", noise, NOISE_FIELDS)
    write_environment(outdir / "environment.txt", REPO, issue7_meta)
    write_report(
        outdir / "REPORT.md",
        tables,
        issue7_meta,
        ordinary,
        rolling,
        c4,
        c46,
        ann,
        amps,
    )
    print(f"wrote {outdir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
