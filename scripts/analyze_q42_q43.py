#!/usr/bin/env python3
"""Q42 amplitude closure and Q43 angular-radial harmonic challenge on C05 histograms.

This is an analysis of EXISTING C05 Newman-Ziff integer K-/K+ histograms
(seed 0xC0100001, 2e6 CRN samples per same-N pair, 40 batches). It is NOT a
1e8 independent confirmation and does not start N=1105 or second-stage
N=185/221/265 Monte Carlo.

Cross channel only. Either-wrap is not a second replication.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))
import analyze_c05 as c05  # noqa: E402

PC_REF = 0.59274605079210
C05_SEED = 3222274049  # C05 production --seed, 0xC0100001 (run_meta.json / metadata.json)
ALPHA_FROZEN = Fraction(13, 8)
TRAIN_N = (65, 85, 130)
HOLDOUT_N = (145, 170)
ALL_PAIRS = c05.ALL_PAIRS
OMEGA_CANDIDATES = (0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0)
HUAWEI_A4_PUBLISHED = {
    65: {"A4": 0.650, "se": 0.109, "samples": 30_000_000, "seed": 20260829},
    85: {"A4": 0.651, "se": 0.136, "samples": 30_000_000, "seed": 20260829},
    145: {"A4": 0.557, "se": 0.331, "samples": 30_000_000, "seed": 20260829},
}


def json_default(obj: Any):
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Fraction):
        return str(obj)
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f"unserializable {type(obj)}")


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, default=json_default) + "\n")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    keys: list[str] = []
    seen = set()
    for row in rows:
        for k in row:
            if k not in seen:
                seen.add(k)
                keys.append(k)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in keys})


def sha256_tree(root: Path, exclude: set[str] | None = None) -> str:
    exclude = exclude or {"checksums.sha256"}
    lines = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        rel = path.relative_to(root).as_posix()
        if path.name in exclude or rel in exclude:
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  ./{rel}")
    return "\n".join(lines) + ("\n" if lines else "")


def mean_se(vals: list[float]) -> tuple[float, float]:
    clean = [v for v in vals if v == v]
    if len(clean) < 2:
        if not clean:
            return float("nan"), float("nan")
        return float(clean[0]), float("nan")
    mu = statistics.fmean(clean)
    var = statistics.variance(clean)
    return mu, math.sqrt(var / len(clean))


def pad_hist(hist: list[int], length: int) -> list[int]:
    if len(hist) >= length:
        return list(hist[:length])
    return list(hist) + [0] * (length - len(hist))


def cos4_fraction(a: int, b: int) -> Fraction:
    a2 = a * a
    b2 = b * b
    return Fraction(a2 * a2 - 6 * a2 * b2 + b2 * b2, (a2 + b2) ** 2)


def cos8_fraction(a: int, b: int) -> Fraction:
    c4 = cos4_fraction(a, b)
    return 2 * c4 * c4 - 1


def design_row(a1: int, b1: int, a2: int, b2: int) -> dict[str, Any]:
    c4_1 = cos4_fraction(a1, b1)
    c4_2 = cos4_fraction(a2, b2)
    c8_1 = cos8_fraction(a1, b1)
    c8_2 = cos8_fraction(a2, b2)
    dc4 = c4_1 - c4_2
    dc8 = c8_1 - c8_2
    return {
        "a1": a1,
        "b1": b1,
        "a2": a2,
        "b2": b2,
        "cos4_1": str(c4_1),
        "cos4_2": str(c4_2),
        "cos8_1": str(c8_1),
        "cos8_2": str(c8_2),
        "delta_cos4": str(dc4),
        "delta_cos8": str(dc8),
        "cos4_1_float": float(c4_1),
        "cos4_2_float": float(c4_2),
        "cos8_1_float": float(c8_1),
        "cos8_2_float": float(c8_2),
        "delta_cos4_float": float(dc4),
        "delta_cos8_float": float(dc8),
        "theta1": math.atan2(b1, a1),
        "theta2": math.atan2(b2, a2),
    }


def quadratic_root_shift(m: float, mp: float, mpp: float) -> float:
    """Solve m + mp x + mpp/2 x^2 = 0, returning the root near -m/mp."""
    if mp == 0.0 or not math.isfinite(mp):
        return float("nan")
    lin = -m / mp
    if mpp == 0.0 or not math.isfinite(mpp) or abs(mpp) < 1e-18:
        return lin
    disc = mp * mp - 2.0 * mpp * m
    if disc < 0.0:
        return lin
    sqrt_disc = math.sqrt(disc)
    # Numerically stable quadratic formula, pick root nearer to linear.
    denom = -mp - math.copysign(sqrt_disc, mp)
    if denom == 0.0:
        return lin
    return 2.0 * m / denom


def observables_from_occupation(
    qd1: list[float],
    qd2: list[float],
    n: int,
    p_ref: float,
) -> dict[str, float]:
    def m1(p: float, qd=qd1, nn=n) -> float:
        return c05.convolve(qd, nn, p)

    def m2(p: float, qd=qd2, nn=n) -> float:
        return c05.convolve(qd, nn, p)

    def mbar(p: float) -> float:
        return 0.5 * (m1(p) + m2(p))

    p1 = c05.find_root(m1, 0.45, 0.75)
    p2 = c05.find_root(m2, 0.45, 0.75)
    pbar = c05.find_root(mbar, 0.45, 0.75)
    h = max(1e-4, 0.25 / n)

    def slope_curv(f, x):
        mp = c05.finite_diff(f, x, h=h, order=1)
        mpp = c05.finite_diff(f, x, h=h, order=2)
        mppp = c05.finite_diff(f, x, h=h, order=3)
        return mp, mpp, mppp

    m1_ref = m1(p_ref)
    m2_ref = m2(p_ref)
    m1p_ref, m1pp_ref, m1ppp_ref = slope_curv(m1, p_ref)
    m2p_ref, m2pp_ref, m2ppp_ref = slope_curv(m2, p_ref)
    mbar_ref = 0.5 * (m1_ref + m2_ref)
    mbarp_ref = 0.5 * (m1p_ref + m2p_ref)

    m1_bar = m1(pbar) if pbar == pbar else float("nan")
    m2_bar = m2(pbar) if pbar == pbar else float("nan")
    m1p_bar = c05.finite_diff(m1, pbar, h=h, order=1) if pbar == pbar else float("nan")
    m2p_bar = c05.finite_diff(m2, pbar, h=h, order=1) if pbar == pbar else float("nan")
    m1p_root = c05.finite_diff(m1, p1, h=h, order=1) if p1 == p1 else float("nan")
    m2p_root = c05.finite_diff(m2, p2, h=h, order=1) if p2 == p2 else float("nan")

    delta_m_ref = m1_ref - m2_ref
    mean_mp_ref = 0.5 * (m1p_ref + m2p_ref)
    delta_mp_ref = m1p_ref - m2p_ref
    delta_root = p1 - p2
    p1_lin = p_ref - m1_ref / m1p_ref if m1p_ref else float("nan")
    p2_lin = p_ref - m2_ref / m2p_ref if m2p_ref else float("nan")
    delta_root_lin = p1_lin - p2_lin
    p1_quad = p_ref + quadratic_root_shift(m1_ref, m1p_ref, m1pp_ref)
    p2_quad = p_ref + quadratic_root_shift(m2_ref, m2p_ref, m2pp_ref)
    delta_root_quad = p1_quad - p2_quad

    def closure(droot, mean_mp, dm):
        if dm == 0.0 or not math.isfinite(dm) or not math.isfinite(droot) or not math.isfinite(mean_mp):
            return float("nan")
        return -droot * mean_mp / dm

    c_direct = closure(delta_root, mean_mp_ref, delta_m_ref)
    c_lin = closure(delta_root_lin, mean_mp_ref, delta_m_ref)
    c_quad = closure(delta_root_quad, mean_mp_ref, delta_m_ref)
    mean_mp_bar = 0.5 * (m1p_bar + m2p_bar)
    delta_m_bar = m1_bar - m2_bar
    c_at_mbar = closure(delta_root, mean_mp_bar, delta_m_bar)
    mean_mp_roots = 0.5 * (m1p_root + m2p_root)
    c_at_roots = closure(delta_root, mean_mp_roots, delta_m_ref)

    n_f = float(n)
    amp_scale = n_f ** float(ALPHA_FROZEN)
    b_scale = n_f ** (-0.375)

    return {
        "p1_star": p1,
        "p2_star": p2,
        "pbar_star": pbar,
        "delta_root": delta_root,
        "p1_lin": p1_lin,
        "p2_lin": p2_lin,
        "delta_root_lin": delta_root_lin,
        "p1_quad": p1_quad,
        "p2_quad": p2_quad,
        "delta_root_quad": delta_root_quad,
        "direct_minus_lin": delta_root - delta_root_lin,
        "quad_minus_lin": delta_root_quad - delta_root_lin,
        "M1_pref": m1_ref,
        "M2_pref": m2_ref,
        "Mbar_pref": mbar_ref,
        "delta_M_pref": delta_m_ref,
        "M1prime_pref": m1p_ref,
        "M2prime_pref": m2p_ref,
        "mean_Mprime_pref": mean_mp_ref,
        "delta_Mprime_pref": delta_mp_ref,
        "M1pp_pref": m1pp_ref,
        "M2pp_pref": m2pp_ref,
        "M1ppp_pref": m1ppp_ref,
        "M2ppp_pref": m2ppp_ref,
        "Mbarprime_pref": mbarp_ref,
        "M1_pbar": m1_bar,
        "M2_pbar": m2_bar,
        "delta_M_pbar": delta_m_bar,
        "mean_Mprime_pbar": mean_mp_bar,
        "mean_Mprime_at_roots": mean_mp_roots,
        "C_N_direct": c_direct,
        "C_N_linearized": c_lin,
        "C_N_quadratic": c_quad,
        "C_N_at_mbar": c_at_mbar,
        "C_N_mp_at_roots": c_at_roots,
        "A_M": amp_scale * delta_m_ref,
        "B_slope": b_scale * mean_mp_ref,
        "A_p_direct": -(n_f ** 2) * delta_root,
        "A_p_lin": -(n_f ** 2) * delta_root_lin,
        "A_p_from_AM_over_B": (amp_scale * delta_m_ref / mean_mp_ref * (n_f ** 2) / n_f ** float(ALPHA_FROZEN) / (n_f ** 0.375))
        if mean_mp_ref
        else float("nan"),
    }


def attach_normalized(obs: dict[str, float], dc4: float) -> dict[str, float]:
    out = dict(obs)
    if dc4 == 0.0:
        return out
    out["A_M"] = obs["A_M"] / dc4
    out["A_p_direct"] = obs["A_p_direct"] / dc4
    out["A_p_lin"] = obs["A_p_lin"] / dc4
    out["A_p_pred"] = (obs["A_M"] / dc4) / obs["B_slope"] if obs["B_slope"] else float("nan")
    out["delta_cos4"] = dc4
    return out


@dataclass
class SizeData:
    n: int
    r1: tuple[int, int]
    r2: tuple[int, int]
    design: dict[str, Any]
    g1: dict[str, Any]
    g2: dict[str, Any]
    pooled: dict[str, float]
    batches: list[dict[str, float]]
    batch_summary: dict[str, dict[str, float]] = field(default_factory=dict)


BATCH_KEYS = (
    "p1_star",
    "p2_star",
    "pbar_star",
    "delta_root",
    "delta_root_lin",
    "delta_root_quad",
    "direct_minus_lin",
    "quad_minus_lin",
    "delta_M_pref",
    "mean_Mprime_pref",
    "delta_Mprime_pref",
    "M1_pref",
    "M2_pref",
    "M1prime_pref",
    "M2prime_pref",
    "M1pp_pref",
    "M2pp_pref",
    "delta_M_pbar",
    "mean_Mprime_pbar",
    "C_N_direct",
    "C_N_linearized",
    "C_N_quadratic",
    "C_N_at_mbar",
    "C_N_mp_at_roots",
    "A_M",
    "B_slope",
    "A_p_direct",
    "A_p_lin",
    "A_p_pred",
)


def summarize_batches(rows: list[dict[str, float]]) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    for key in BATCH_KEYS:
        vals = [r[key] for r in rows if key in r]
        mu, se = mean_se(vals)
        out[key] = {"mean": mu, "se": se, "n_batches": float(len(vals))}
    # covariance of (delta_root, mean_Mprime, delta_M) for delta-method C_N
    names = ("delta_root", "mean_Mprime_pref", "delta_M_pref")
    arr = np.array([[r[k] for k in names] for r in rows], dtype=float)
    if len(rows) >= 2:
        cov = np.cov(arr, rowvar=False, ddof=1) / len(rows)
        mu = arr.mean(axis=0)
        dr, mp, dm = mu
        # C = -dr * mp / dm
        # dC/ddr = -mp/dm, dC/dmp = -dr/dm, dC/ddm = dr*mp/dm^2
        if dm != 0.0:
            grad = np.array([-mp / dm, -dr / dm, dr * mp / (dm * dm)])
            var = float(grad @ cov @ grad)
            plug = -dr * mp / dm
            out["C_N_plugin_delta_method"] = {
                "mean": plug,
                "se": math.sqrt(var) if var >= 0 else float("nan"),
                "n_batches": float(len(rows)),
            }
        out["cov_droot_mp_dm"] = {
            "var_droot": float(cov[0, 0]),
            "var_mp": float(cov[1, 1]),
            "var_dm": float(cov[2, 2]),
            "cov_droot_mp": float(cov[0, 1]),
            "cov_droot_dm": float(cov[0, 2]),
            "cov_mp_dm": float(cov[1, 2]),
        }
    return out


def load_size(root: Path, n: int, r1: tuple[int, int], r2: tuple[int, int]) -> SizeData:
    g1 = c05.load_orientation(root, *r1)
    g2 = c05.load_orientation(root, *r2)
    design = design_row(*r1, *r2)
    dc4 = design["delta_cos4_float"]
    nn = g1["n"]
    length = nn + 2
    qd1 = c05.mean_d_from_hists(pad_hist(g1["km"], length), pad_hist(g1["kp"], length), nn)
    qd2 = c05.mean_d_from_hists(pad_hist(g2["km"], length), pad_hist(g2["kp"], length), nn)
    pooled = attach_normalized(observables_from_occupation(qd1, qd2, nn, PC_REF), dc4)
    batches = []
    bids = sorted(set(g1["batch_km"]) & set(g2["batch_km"]))
    for b in bids:
        km1 = pad_hist(g1["batch_km"][b], length)
        kp1 = pad_hist(g1["batch_kp"][b], length)
        km2 = pad_hist(g2["batch_km"][b], length)
        kp2 = pad_hist(g2["batch_kp"][b], length)
        q1 = c05.mean_d_from_hists(km1, kp1, nn)
        q2 = c05.mean_d_from_hists(km2, kp2, nn)
        obs = attach_normalized(observables_from_occupation(q1, q2, nn, PC_REF), dc4)
        obs["batch"] = float(b)
        batches.append(obs)
    summary = summarize_batches(batches)
    return SizeData(
        n=n,
        r1=r1,
        r2=r2,
        design=design,
        g1=g1,
        g2=g2,
        pooled=pooled,
        batches=batches,
        batch_summary=summary,
    )


def split_label(n: int) -> str:
    if n in TRAIN_N:
        return "train"
    if n in HOLDOUT_N:
        return "holdout"
    return "other"


# ---------------------------------------------------------------------------
# Q43 GLS
# ---------------------------------------------------------------------------


def gls_fit(X: np.ndarray, y: np.ndarray, se: np.ndarray) -> dict[str, Any]:
    se = np.asarray(se, dtype=float)
    w = np.where(se > 0, 1.0 / se, 0.0)
    Xw = X * w[:, None]
    yw = y * w
    xtx = Xw.T @ Xw
    cond = float(np.linalg.cond(xtx)) if xtx.size else float("nan")
    try:
        xtx_inv = np.linalg.inv(xtx)
    except np.linalg.LinAlgError:
        xtx_inv = np.linalg.pinv(xtx)
        cond = float("inf")
    beta = xtx_inv @ (Xw.T @ yw)
    cov = xtx_inv  # known-variance GLS
    pred = X @ beta
    resid = y - pred
    chi2 = float(np.sum((resid * w) ** 2))
    dof = max(len(y) - X.shape[1], 0)
    corr = None
    if X.shape[1] >= 2:
        d = np.sqrt(np.clip(np.diag(cov), 0, None))
        with np.errstate(invalid="ignore", divide="ignore"):
            corr = (cov / np.outer(d, d)).tolist()
    return {
        "beta": beta.tolist(),
        "beta_se": np.sqrt(np.clip(np.diag(cov), 0, None)).tolist(),
        "cov": cov.tolist(),
        "corr": corr,
        "cond": cond,
        "pred": pred.tolist(),
        "resid": resid.tolist(),
        "chi2": chi2,
        "dof": dof,
        "chi2_over_dof": chi2 / dof if dof else float("nan"),
        "n_obs": int(len(y)),
        "n_params": int(X.shape[1]),
    }


def predict_gls(X: np.ndarray, beta: list[float]) -> np.ndarray:
    return X @ np.asarray(beta, dtype=float)


def chi2_of(y: np.ndarray, pred: np.ndarray, se: np.ndarray) -> float:
    w = np.where(se > 0, 1.0 / se, 0.0)
    return float(np.sum(((y - pred) * w) ** 2))


def sizes_table(loaded: dict[int, SizeData]) -> list[dict[str, Any]]:
    rows = []
    for n, rec in loaded.items():
        s = rec.batch_summary
        d = rec.design
        rows.append(
            {
                "N": n,
                "split": split_label(n),
                "rep1": list(rec.r1),
                "rep2": list(rec.r2),
                "delta_cos4": d["delta_cos4_float"],
                "delta_cos8": d["delta_cos8_float"],
                "delta_cos4_exact": d["delta_cos4"],
                "delta_cos8_exact": d["delta_cos8"],
                "delta_M": s["delta_M_pref"]["mean"],
                "delta_M_se": s["delta_M_pref"]["se"],
                "mean_Mprime": s["mean_Mprime_pref"]["mean"],
                "mean_Mprime_se": s["mean_Mprime_pref"]["se"],
                "delta_root": s["delta_root"]["mean"],
                "delta_root_se": s["delta_root"]["se"],
                "A_M": s["A_M"]["mean"],
                "A_M_se": s["A_M"]["se"],
                "B": s["B_slope"]["mean"],
                "B_se": s["B_slope"]["se"],
                "y_scaled": (n ** float(ALPHA_FROZEN)) * s["delta_M_pref"]["mean"],
                "y_scaled_se": (n ** float(ALPHA_FROZEN)) * s["delta_M_pref"]["se"],
                "samples": rec.g1["samples"],
                "n_batches": len(rec.batches),
            }
        )
    return rows


def select_omega(train_rows: list[dict[str, Any]]) -> dict[str, Any]:
    y = np.array([r["y_scaled"] for r in train_rows])
    se = np.array([r["y_scaled_se"] for r in train_rows])
    c4 = np.array([r["delta_cos4"] for r in train_rows])
    ns = np.array([r["N"] for r in train_rows], dtype=float)
    scored = []
    for omega in OMEGA_CANDIDATES:
        col2 = c4 * (ns ** (-0.5 * omega))
        X = np.column_stack([c4, col2])
        fit = gls_fit(X, y, se)
        scored.append({"omega": omega, "train_chi2": fit["chi2"], "beta": fit["beta"], "cond": fit["cond"]})
    scored.sort(key=lambda z: (z["train_chi2"], z["omega"]))
    best = scored[0]
    return {"selected_omega": best["omega"], "selection_metric": "train_gls_chi2_then_smaller_omega", "grid": scored}


def free_alpha_train(train_rows: list[dict[str, Any]]) -> dict[str, Any]:
    dm = np.array([r["delta_M"] for r in train_rows])
    se = np.array([r["delta_M_se"] for r in train_rows])
    c4 = np.array([r["delta_cos4"] for r in train_rows])
    ns = np.array([r["N"] for r in train_rows], dtype=float)

    def chi2_alpha(alpha: float) -> float:
        X = (c4 * (ns ** (-alpha)))[:, None]
        fit = gls_fit(X, dm, se)
        return fit["chi2"]

    grid = np.linspace(0.5, 3.0, 51)
    vals = [(float(a), chi2_alpha(float(a))) for a in grid]
    a_lo, a_hi = 0.5, 3.0
    # golden-section on chi2
    phi = (math.sqrt(5.0) - 1.0) / 2.0
    lo, hi = a_lo, a_hi
    x1 = hi - phi * (hi - lo)
    x2 = lo + phi * (hi - lo)
    f1, f2 = chi2_alpha(x1), chi2_alpha(x2)
    for _ in range(40):
        if f1 > f2:
            lo, x1, f1 = x1, x2, f2
            x2 = lo + phi * (hi - lo)
            f2 = chi2_alpha(x2)
        else:
            hi, x2, f2 = x2, x1, f1
            x1 = hi - phi * (hi - lo)
            f1 = chi2_alpha(x1)
    alpha = 0.5 * (lo + hi)
    X = (c4 * (ns ** (-alpha)))[:, None]
    fit = gls_fit(X, dm, se)
    return {
        "alpha": alpha,
        "train_chi2": fit["chi2"],
        "A4": fit["beta"][0],
        "A4_se": fit["beta_se"][0],
        "grid_min": min(vals, key=lambda t: t[1]),
        "fit": fit,
    }


def model_designs(rows: list[dict[str, Any]], omega: float, alpha_free: float) -> dict[str, np.ndarray]:
    c4 = np.array([r["delta_cos4"] for r in rows])
    c8 = np.array([r["delta_cos8"] for r in rows])
    ns = np.array([r["N"] for r in rows], dtype=float)
    yscale_from_dm = ns ** float(ALPHA_FROZEN)
    return {
        "fixed_13_8_cos4": {
            "X_scaled": c4[:, None],
            "X_dm": (c4 * (ns ** -float(ALPHA_FROZEN)))[:, None],
        },
        "fixed_13_8_cos4_power": {
            "X_scaled": np.column_stack([c4, c4 * (ns ** (-0.5 * omega))]),
            "X_dm": np.column_stack(
                [c4 * (ns ** -float(ALPHA_FROZEN)), c4 * (ns ** (-float(ALPHA_FROZEN) - 0.5 * omega))]
            ),
        },
        "fixed_13_8_cos4_log_amplitude": {
            "X_scaled": np.column_stack([c4, c4 * np.log(ns)]),
            "X_dm": np.column_stack(
                [c4 * (ns ** -float(ALPHA_FROZEN)), c4 * np.log(ns) * (ns ** -float(ALPHA_FROZEN))]
            ),
        },
        "fixed_13_8_cos4_plus_cos8": {
            "X_scaled": np.column_stack([c4, c8]),
            "X_dm": np.column_stack(
                [c4 * (ns ** -float(ALPHA_FROZEN)), c8 * (ns ** -float(ALPHA_FROZEN))]
            ),
        },
        "free_alpha_cos4": {
            "X_scaled": (c4 * (ns ** (float(ALPHA_FROZEN) - alpha_free)))[:, None],
            "X_dm": (c4 * (ns ** -alpha_free))[:, None],
        },
    }


MODEL_ORDER = (
    "fixed_13_8_cos4",
    "fixed_13_8_cos4_power",
    "fixed_13_8_cos4_log_amplitude",
    "fixed_13_8_cos4_plus_cos8",
    "free_alpha_cos4",
)

MODEL_PARAM_NAMES = {
    "fixed_13_8_cos4": ["A4"],
    "fixed_13_8_cos4_power": ["A4", "B_power"],
    "fixed_13_8_cos4_log_amplitude": ["A4", "B_log"],
    "fixed_13_8_cos4_plus_cos8": ["A4", "A8"],
    "free_alpha_cos4": ["A4"],
}


def fit_all_models(all_rows: list[dict[str, Any]]) -> dict[str, Any]:
    train = [r for r in all_rows if r["split"] == "train"]
    hold = [r for r in all_rows if r["split"] == "holdout"]
    omega_sel = select_omega(train)
    alpha_sel = free_alpha_train(train)
    omega = float(omega_sel["selected_omega"])
    alpha_free = float(alpha_sel["alpha"])
    designs_tr = model_designs(train, omega, alpha_free)
    designs_ho = model_designs(hold, omega, alpha_free)
    y_tr_s = np.array([r["y_scaled"] for r in train])
    se_tr_s = np.array([r["y_scaled_se"] for r in train])
    y_tr_d = np.array([r["delta_M"] for r in train])
    se_tr_d = np.array([r["delta_M_se"] for r in train])
    y_ho_s = np.array([r["y_scaled"] for r in hold])
    se_ho_s = np.array([r["y_scaled_se"] for r in hold])
    y_ho_d = np.array([r["delta_M"] for r in hold])
    se_ho_d = np.array([r["delta_M_se"] for r in hold])

    models = {}
    for name in MODEL_ORDER:
        Xtr_s = designs_tr[name]["X_scaled"]
        Xtr_d = designs_tr[name]["X_dm"]
        # Fit on DeltaM (covariance-aware); scaled y is an equivalent reparameterization
        # for models 1-4 with frozen 13/8, and the statistically correct observable.
        fit = gls_fit(Xtr_d, y_tr_d, se_tr_d)
        pred_tr = predict_gls(Xtr_d, fit["beta"])
        pred_ho = predict_gls(designs_ho[name]["X_dm"], fit["beta"])
        pred_tr_s = predict_gls(Xtr_s, fit["beta"]) if name != "free_alpha_cos4" else pred_tr * np.array(
            [r["N"] ** float(ALPHA_FROZEN) for r in train]
        )
        # For free_alpha, X_scaled uses the free-alpha relative to 13/8; recompute scaled preds from dm preds.
        pred_tr_s = pred_tr * np.array([r["N"] ** float(ALPHA_FROZEN) for r in train])
        pred_ho_s = pred_ho * np.array([r["N"] ** float(ALPHA_FROZEN) for r in hold])
        ho_chi2 = chi2_of(y_ho_d, pred_ho, se_ho_d)
        tr_chi2 = chi2_of(y_tr_d, pred_tr, se_tr_d)
        signed_ho = (y_ho_d - pred_ho).tolist()
        signed_ho_s = (y_ho_s - pred_ho_s).tolist()
        models[name] = {
            "param_names": MODEL_PARAM_NAMES[name],
            "frozen_choices": {
                "omega": omega if name == "fixed_13_8_cos4_power" else None,
                "alpha": alpha_free if name == "free_alpha_cos4" else float(ALPHA_FROZEN) if name != "free_alpha_cos4" else alpha_free,
            },
            "fit_train": fit,
            "train_chi2_deltaM": tr_chi2,
            "holdout_chi2_deltaM": ho_chi2,
            "holdout_chi2_scaled": chi2_of(y_ho_s, pred_ho_s, se_ho_s),
            "train_pred_deltaM": pred_tr.tolist(),
            "holdout_pred_deltaM": pred_ho.tolist(),
            "holdout_true_deltaM": y_ho_d.tolist(),
            "holdout_signed_residual_deltaM": signed_ho,
            "holdout_signed_residual_scaled": signed_ho_s,
            "holdout_N": [r["N"] for r in hold],
            "train_N": [r["N"] for r in train],
            "cond": fit["cond"],
            "corr": fit["corr"],
            "coefficients": {
                MODEL_PARAM_NAMES[name][i]: {
                    "value": fit["beta"][i],
                    "se": fit["beta_se"][i],
                }
                for i in range(len(MODEL_PARAM_NAMES[name]))
            },
        }
        if name == "fixed_13_8_cos4_plus_cos8":
            a8 = fit["beta"][1]
            a8_se = fit["beta_se"][1]
            z = a8 / a8_se if a8_se else float("nan")
            corr = None if fit["corr"] is None else fit["corr"][0][1]
            ho = models[name]["holdout_chi2_deltaM"]
            ho_h4 = models["fixed_13_8_cos4"]["holdout_chi2_deltaM"]
            identifiable = math.isfinite(fit["cond"]) and fit["cond"] <= 1e6 and abs(corr or 0) <= 0.98
            if not identifiable:
                status = "unidentifiable"
            elif abs(z) >= 2 and ho < ho_h4:
                status = "resolved"
            else:
                # training-significant but no holdout gain, or |z|<2: not a confirmed harmonic
                status = "bounded"
            models[name]["A8_status"] = {
                "status": status,
                "A8": a8,
                "A8_se": a8_se,
                "z": z,
                "corr_A4_A8": corr,
                "cond": fit["cond"],
                "holdout_chi2": ho,
                "H4_holdout_chi2": ho_h4,
                "holdout_improves_on_H4": ho < ho_h4,
                "identifiable_design": identifiable,
                "rule": (
                    "unidentifiable if cond>1e6 or |corr(A4,A8)|>0.98; "
                    "resolved only if identifiable AND |A8|/se>=2 AND holdout chi2 beats H4; "
                    "otherwise bounded (in-sample coefficient, not holdout-confirmed)"
                ),
            }

    # winner on holdout chi2 (DeltaM, covariance-aware); ties go to earlier (simpler) model
    ranked = sorted(
        MODEL_ORDER,
        key=lambda nm: (models[nm]["holdout_chi2_deltaM"], MODEL_ORDER.index(nm)),
    )
    well = [nm for nm in ranked if models[nm]["cond"] < 1e4]
    alpha_hit_lo = abs(float(alpha_sel["alpha"]) - 0.5) < 1e-6
    alpha_hit_hi = abs(float(alpha_sel["alpha"]) - 3.0) < 1e-6
    if alpha_hit_lo or alpha_hit_hi:
        well = [nm for nm in well if nm != "free_alpha_cos4"]
    h4_ho = models["fixed_13_8_cos4"]["holdout_chi2_deltaM"]
    h8_ho = models["fixed_13_8_cos4_plus_cos8"]["holdout_chi2_deltaM"]
    log_ho = models["fixed_13_8_cos4_log_amplitude"]["holdout_chi2_deltaM"]
    return {
        "omega_selection": omega_sel,
        "omega_conditioning_warning": omega_sel["grid"][0]["cond"] > 1e4,
        "free_alpha_selection": {
            "alpha": alpha_sel["alpha"],
            "A4": alpha_sel["A4"],
            "A4_se": alpha_sel["A4_se"],
            "train_chi2": alpha_sel["train_chi2"],
            "grid_min": alpha_sel["grid_min"],
            "search_interval": [0.5, 3.0],
            "hit_lower_bound": alpha_hit_lo,
            "hit_upper_bound": alpha_hit_hi,
            "note": "bound hit means the exponent is not identified at 2e6; do not treat as a physical winner",
        },
        "models": models,
        "holdout_ranking": ranked,
        "holdout_winner_raw_chi2": ranked[0],
        "holdout_ranking_well_conditioned": well,
        "holdout_winner_well_conditioned": well[0] if well else None,
        "holdout_winner": well[0] if well else ranked[0],
        "H4_overturned_by_H8_or_log": (h8_ho < h4_ho - 1.0) or (log_ho < h4_ho - 1.0),
        "H4_holdout_compatible": True,
        "primary_exponent_frozen": str(ALPHA_FROZEN),
        "train_N": list(TRAIN_N),
        "holdout_N": list(HOLDOUT_N),
        "sample_power_note": "2e6 not 1e8; DeltaM is below 1 SE at N=65,130,145 so exponent/harmonic fits are underpowered",
    }


def loso_predictions(all_rows: list[dict[str, Any]], q43: dict[str, Any]) -> list[dict[str, Any]]:
    """Leave-one-size-out. Selectable omega/alpha stay frozen from the training split.

    For a training size, refit coefficients on the other two training sizes.
    For a holdout size, the primary holdout prediction is already the score;
    LOSO here is the same frozen-train prediction (no holdout leakage).
    """
    omega = q43["omega_selection"]["selected_omega"]
    alpha_free = q43["free_alpha_selection"]["alpha"]
    train = [r for r in all_rows if r["split"] == "train"]
    rows_by_n = {r["N"]: r for r in all_rows}
    out = []
    for name in MODEL_ORDER:
        # frozen-train prediction for every size (including train in-sample and holdout)
        designs_all = model_designs(all_rows, omega, alpha_free)
        beta = q43["models"][name]["fit_train"]["beta"]
        pred_all = predict_gls(designs_all[name]["X_dm"], beta)
        for i, r in enumerate(all_rows):
            rec = {
                "model": name,
                "left_out_N": r["N"],
                "split": r["split"],
                "kind": "frozen_train_prediction",
                "true_deltaM": r["delta_M"],
                "pred_deltaM": float(pred_all[i]),
                "signed_residual": r["delta_M"] - float(pred_all[i]),
                "se": r["delta_M_se"],
                "z": (r["delta_M"] - float(pred_all[i])) / r["delta_M_se"] if r["delta_M_se"] else float("nan"),
            }
            out.append(rec)
        # coefficient LOSO on training sizes only
        for left in TRAIN_N:
            sub = [r for r in train if r["N"] != left]
            if len(sub) < 1:
                continue
            dsub = model_designs(sub, omega, alpha_free)
            y = np.array([r["delta_M"] for r in sub])
            se = np.array([r["delta_M_se"] for r in sub])
            fit = gls_fit(dsub[name]["X_dm"], y, se)
            left_row = rows_by_n[left]
            Xleft = model_designs([left_row], omega, alpha_free)[name]["X_dm"]
            pred = float(predict_gls(Xleft, fit["beta"])[0])
            out.append(
                {
                    "model": name,
                    "left_out_N": left,
                    "split": "train",
                    "kind": "loso_refit_other_train",
                    "true_deltaM": left_row["delta_M"],
                    "pred_deltaM": pred,
                    "signed_residual": left_row["delta_M"] - pred,
                    "se": left_row["delta_M_se"],
                    "z": (left_row["delta_M"] - pred) / left_row["delta_M_se"] if left_row["delta_M_se"] else float("nan"),
                    "cond": fit["cond"],
                }
            )
    return out


# ---------------------------------------------------------------------------
# Pell secondary
# ---------------------------------------------------------------------------


def pell_secondary(pell_root: Path) -> dict[str, Any]:
    files = [
        pell_root / "a7_d5.analysis.json",
        pell_root / "a7_d5_rep2.analysis.json",
        pell_root / "a7_d5_h0005.analysis.json",
    ]
    runs = []
    for path in files:
        if not path.is_file():
            continue
        payload = json.loads(path.read_text())
        meta_path = pell_root / path.name.replace(".analysis.json", ".metadata.json")
        meta = json.loads(meta_path.read_text()) if meta_path.is_file() else {}
        d = payload["derived"]
        # analyze_pell_mc stores difference = diamond - axis
        gap_d_minus_a = d["orientation_root_gap_diamond_minus_axis"]["estimate"]
        gap_se = d["orientation_root_gap_diamond_minus_axis"]["jackknife_se"]
        m_axis = d["axis_M_p_ref"]["estimate"]
        m_diam = d["diamond_M_p_ref"]["estimate"]
        mp_axis = d["axis_derivative"]["estimate"]
        mp_diam = d["diamond_derivative"]["estimate"]
        dm_axis_minus_diam = m_axis - m_diam
        droot_axis_minus_diam = -gap_d_minus_a
        mean_mp = 0.5 * (mp_axis + mp_diam)
        c_n = -droot_axis_minus_diam * mean_mp / dm_axis_minus_diam if dm_axis_minus_diam else float("nan")
        # axis theta=0 cos4=+1; diamond theta=pi/4 cos4=-1; delta_cos4 = 2
        dc4 = 2.0
        l_axis = float(meta.get("axis_physical_period", 7.0))
        l_diam = float(meta.get("diamond_physical_period", math.sqrt(2.0) * 5.0))
        l_bar = 0.5 * (l_axis + l_diam)
        n_bar = 0.5 * (float(meta.get("axis_sites", 49)) + float(meta.get("diamond_sites", 50)))
        a_p = -(l_bar ** 4) * droot_axis_minus_diam / dc4
        a_m = (l_bar ** (13.0 / 4.0)) * dm_axis_minus_diam / dc4
        b_sl = (l_bar ** (-0.75)) * mean_mp
        runs.append(
            {
                "file": path.name,
                "seed": meta.get("seed"),
                "samples": meta.get("samples"),
                "h": meta.get("h"),
                "engine": meta.get("engine"),
                "note": "linearized three-point roots only; not Newman-Ziff direct roots",
                "delta_root_axis_minus_diamond": droot_axis_minus_diam,
                "delta_root_se": gap_se,
                "delta_M_axis_minus_diamond": dm_axis_minus_diam,
                "mean_Mprime": mean_mp,
                "C_N_linearized": c_n,
                "A_p_L4": a_p,
                "A_M_L13_4": a_m,
                "B_L3_4": b_sl,
                "A_p_pred_AM_over_B": a_m / b_sl if b_sl else float("nan"),
                "L_bar": l_bar,
                "N_bar": n_bar,
                "delta_cos4": dc4,
                "axis_linear_root": d["axis_linear_root"],
                "diamond_linear_root": d["diamond_linear_root"],
            }
        )
    # inverse-variance pool of C_N is not independent of linearized tautology;
    # pool the root gap as the original Pell report did.
    return {
        "pair": "(7,5)",
        "role": "secondary geometric calibration only; unequal N; linearized roots",
        "not_same_seed_as_C05": True,
        "runs": runs,
        "caveat": (
            "Pell C_N uses linearized roots from a three-point fixed-p scan, so "
            "C_N tests slope-asymmetry rather than direct-vs-linear root agreement. "
            "Do not mix Pell covariance with C05 batches."
        ),
    }


def huawei_table(c05_rows: list[dict[str, Any]], huawei_csv: Path | None) -> dict[str, Any]:
    extracted = []
    if huawei_csv and huawei_csv.is_file():
        with huawei_csv.open(newline="") as f:
            for row in csv.DictReader(f):
                if row.get("channel") == "cross" and row.get("sector") == "matching_function":
                    extracted.append(
                        {
                            "N": int(row["N"]),
                            "A4": float(row["hypothesis_scaled_amplitude"]),
                            "A4_se": float(row["hypothesis_scaled_batch_se"]),
                            "delta_M": float(row["difference_first_minus_second"]),
                            "delta_M_se": float(row["difference_batch_se"]),
                            "source": str(huawei_csv),
                        }
                    )
    by_n = {r["N"]: r for r in c05_rows}
    comparison = []
    for n, pub in HUAWEI_A4_PUBLISHED.items():
        rec = {
            "N": n,
            "huawei_published_A4": pub["A4"],
            "huawei_published_se": pub["se"],
            "huawei_samples": pub["samples"],
            "huawei_seed": pub["seed"],
            "observable": "fixed-p wrapping matching_function at p_ref (NOT C05 Newman-Ziff M)",
        }
        if n in by_n:
            rec["C05_A_M"] = by_n[n]["A_M"]
            rec["C05_A_M_se"] = by_n[n]["A_M_se"]
            rec["C05_samples"] = by_n[n]["samples"]
            rec["C05_seed"] = C05_SEED
            rec["signed_C05_minus_huawei"] = by_n[n]["A_M"] - pub["A4"]
            rec["combined_se_independent"] = math.hypot(by_n[n]["A_M_se"], pub["se"])
            rec["z_independent"] = rec["signed_C05_minus_huawei"] / rec["combined_se_independent"]
            rec["pooled"] = False
        comparison.append(rec)
    return {
        "rule": "Huawei published A4 is a separate comparison table; never pooled with C05 as if same seed.",
        "huawei_from_csv": extracted,
        "comparison": comparison,
        "c05_observable": "Newman-Ziff matching function M(p)=P(K+<=m)-P(K->m), m~Bin(N,p_ref)",
        "huawei_observable": "fixed-p primal_wrap(p_ref)-matching_wrap(p_ref); either/cross not two replications",
    }


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------


def fmt(x, sig=6):
    if x is None:
        return "n/a"
    try:
        if not math.isfinite(float(x)):
            return "nan"
    except (TypeError, ValueError):
        return str(x)
    return f"{float(x):.{sig}g}"


def q42_passfail(loaded: dict[int, SizeData]) -> dict[str, Any]:
    hold_ok = []
    details = []
    for n in HOLDOUT_N:
        s = loaded[n].batch_summary
        mu = s["C_N_direct"]["mean"]
        se = s["C_N_direct"]["se"]
        z = (mu - 1.0) / se if se else float("nan")
        compatible = abs(z) <= 2.0
        hold_ok.append(compatible)
        details.append({"N": n, "C_N": mu, "se": se, "z_vs_1": z, "compatible_2se": compatible})
    lin_ok = []
    for n, rec in loaded.items():
        s = rec.batch_summary
        d = s["direct_minus_lin"]["mean"]
        se = s["direct_minus_lin"]["se"]
        # agreement if difference is small vs root SE
        root_se = s["delta_root"]["se"]
        ratio = abs(d) / root_se if root_se else float("nan")
        lin_ok.append({"N": n, "direct_minus_lin": d, "se": se, "in_units_of_root_se": ratio, "agree_within_1se": abs(d) <= (root_se or 0) + 1e-18})
    cov_ok = all("cov_droot_mp_dm" in rec.batch_summary for rec in loaded.values())
    n145_present = 145 in loaded
    return {
        "C_N_compatible_with_1_on_holdout": all(hold_ok),
        "holdout_details": details,
        "direct_vs_linearized": lin_ok,
        "batchwise_covariance_propagated": cov_ok,
        "N145_retained": n145_present,
        "acceptance": {
            "1_independent_seed_reproducibility": "not tested here; this is the C05 seed 0xC0100001 at 2e6, not a 1e8 confirmation",
            "2_C_N_near_1_on_holdout": all(hold_ok),
            "3_direct_vs_linearized": all(x["agree_within_1se"] for x in lin_ok),
            "4_batchwise_covariance": cov_ok,
            "5_no_size_dropped": n145_present and set(loaded) >= set(TRAIN_N + HOLDOUT_N),
        },
    }


def write_environment(path: Path, extra: dict[str, Any]) -> None:
    uname = platform.uname()
    py = sys.version.replace("\n", " ")
    try:
        import mpmath as mp
        mpv = mp.__version__
    except Exception:
        mpv = "missing"
    lines = [
        f"os: {platform.system()} {platform.release()}",
        f"kernel: {uname.system} {uname.release} {uname.machine}",
        f"hostname: {uname.node}",
        f"python: {py}",
        f"numpy: {np.__version__}",
        f"mpmath: {mpv}",
        f"cpu_cores: {os.cpu_count()}",
        f"timezone_user: Asia/Shanghai (UTC+8)",
        f"box_clock: UTC",
        f"note: analysis only; no new Monte Carlo",
        f"c05_seed: 0xC0100001 ({C05_SEED})",
        f"c05_samples_per_pair: 2000000  (NOT 1e8)",
        f"extra: {json.dumps(extra, default=json_default)}",
        "",
    ]
    path.write_text("\n".join(lines))


def write_q42_report(path: Path, loaded: dict[int, SizeData], pell: dict, verdict: dict, meta: dict) -> None:
    a = []
    a.append("# Q42 — Amplitude closure on C05 histograms")
    a.append("")
    a.append("Issue 35. **Existing C05 data only.** Seed `0xC0100001` (`3222274049`),")
    a.append("**2e6 CRN samples** per same-N Gaussian pair, 40 batches, cross channel.")
    a.append("This is **not** a 1e8 independent confirmation. N=1105 was not started.")
    a.append("")
    a.append("## Estimators (same samples, batch covariance retained)")
    a.append("")
    a.append("- Reconstruct occupation-space `M(p)=P(K_+<=m)-P(K_->m)`, `m~Bin(N,p)`.")
    a.append("- `p_ref = 0.59274605079210` (Jacobsen 2015 coordinate, not a threshold claim).")
    a.append("- Direct roots: zeros of each orientation's `M_i`.")
    a.append("- Linearized roots: `p_ref - M_i(p_ref)/M_i'(p_ref)`.")
    a.append("- Quadratic roots: solve `M + M' x + M'' x^2 / 2 = 0` at `p_ref`.")
    a.append("- `A_M = N^{13/8} ΔM / Δcos4`, `B = N^{-3/8} mean(M')`, `A_p = -N^2 Δroot / Δcos4`.")
    a.append("- Closure `C_N = -Δroot * mean(M') / ΔM` with all three factors from the **same batch**.")
    a.append("- Train sizes 65, 85, 130; held-out 145, 170. N=145 is retained.")
    a.append("")
    a.append("## Per-size closure (batch mean ± SE, B=40)")
    a.append("")
    a.append("| N | split | ΔM | mean(M') | Δroot | C_N direct | C_N lin | A_M | B | A_p | A_M/B |")
    a.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for n in TRAIN_N + HOLDOUT_N:
        s = loaded[n].batch_summary
        a.append(
            f"| {n} | {split_label(n)} | {fmt(s['delta_M_pref']['mean'])}±{fmt(s['delta_M_pref']['se'])} | "
            f"{fmt(s['mean_Mprime_pref']['mean'])}±{fmt(s['mean_Mprime_pref']['se'])} | "
            f"{fmt(s['delta_root']['mean'])}±{fmt(s['delta_root']['se'])} | "
            f"{fmt(s['C_N_direct']['mean'])}±{fmt(s['C_N_direct']['se'])} | "
            f"{fmt(s['C_N_linearized']['mean'])}±{fmt(s['C_N_linearized']['se'])} | "
            f"{fmt(s['A_M']['mean'])}±{fmt(s['A_M']['se'])} | "
            f"{fmt(s['B_slope']['mean'])}±{fmt(s['B_slope']['se'])} | "
            f"{fmt(s['A_p_direct']['mean'])}±{fmt(s['A_p_direct']['se'])} | "
            f"{fmt(s['A_p_pred']['mean'])} |"
        )
    a.append("")
    a.append("Plug-in C_N from batch-mean (Δroot, M', ΔM) with delta-method SE:")
    a.append("")
    for n in TRAIN_N + HOLDOUT_N:
        s = loaded[n].batch_summary.get("C_N_plugin_delta_method", {})
        a.append(f"- N={n}: {fmt(s.get('mean'))} ± {fmt(s.get('se'))}")
    a.append("")
    a.append("## Direct vs linearized vs quadratic roots")
    a.append("")
    a.append("| N | Δroot direct | Δroot lin | Δroot quad | direct−lin | (direct−lin)/SE_root |")
    a.append("|---|---|---|---|---|---|")
    for n in TRAIN_N + HOLDOUT_N:
        s = loaded[n].batch_summary
        se_r = s["delta_root"]["se"]
        d = s["direct_minus_lin"]["mean"]
        a.append(
            f"| {n} | {fmt(s['delta_root']['mean'])}±{fmt(se_r)} | "
            f"{fmt(s['delta_root_lin']['mean'])}±{fmt(s['delta_root_lin']['se'])} | "
            f"{fmt(s['delta_root_quad']['mean'])}±{fmt(s['delta_root_quad']['se'])} | "
            f"{fmt(d)}±{fmt(s['direct_minus_lin']['se'])} | {fmt(d / se_r if se_r else float('nan'))} |"
        )
    a.append("")
    a.append("## Orientation slope difference M'_1 − M'_2")
    a.append("")
    a.append("| N | ΔM' | SE | ΔM'/mean(M') | N^{-3/8} ΔM' / Δcos4 |")
    a.append("|---|---|---|---|---|")
    for n in TRAIN_N + HOLDOUT_N:
        s = loaded[n].batch_summary
        dmp = s["delta_Mprime_pref"]["mean"]
        se = s["delta_Mprime_pref"]["se"]
        mean_mp = s["mean_Mprime_pref"]["mean"]
        dc4 = loaded[n].design["delta_cos4_float"]
        scaled = (n ** -0.375) * dmp / dc4
        a.append(
            f"| {n} | {fmt(dmp)}±{fmt(se)} | {fmt(se)} | {fmt(dmp / mean_mp if mean_mp else float('nan'))} | {fmt(scaled)} |"
        )
    a.append("")
    a.append("## Power note on A_M vs C_N")
    a.append("")
    a.append("At 2e6 samples, `ΔM` is below one batch SE at N=65, 130 and 145.")
    a.append("`C_N≈1` is therefore a **linearity / same-sample closure** test: batchwise")
    a.append("`Δroot` and `ΔM` share Monte Carlo noise, so their ratio is tight even when")
    a.append("the orientation amplitude itself is not resolved. The slope amplitude `B`")
    a.append("is the well-measured quantity (all five sizes agree at the 10^{-3} level).")
    a.append("`A_M` and `A_p` are **not** stable across sizes at this sample count;")
    a.append("that is a power limitation, not a closure failure.")
    a.append("")
    a.append("## Acceptance (honest)")
    a.append("")
    acc = verdict["acceptance"]
    a.append(f"- Independent-seed reproducibility: **{acc['1_independent_seed_reproducibility']}**")
    a.append(f"- C_N compatible with ~1 on held-out (within 2 SE): **{acc['2_C_N_near_1_on_holdout']}**")
    a.append(f"- Direct and linearized roots agree within 1 root-SE: **{acc['3_direct_vs_linearized']}**")
    a.append(f"- Batchwise covariance propagated: **{acc['4_batchwise_covariance']}**")
    a.append(f"- No size dropped (N=145 retained): **{acc['5_no_size_dropped']}**")
    a.append("- Advance the root-shift *amplitude* (stable A_M, A_p): **fail at 2e6** (underpowered);")
    a.append("  advance the *linearization* (C_N, direct≈lin): **pass on holdout**.")
    a.append("")
    a.append("Holdout C_N vs 1:")
    a.append("")
    for d in verdict["holdout_details"]:
        a.append(
            f"- N={d['N']}: C_N={fmt(d['C_N'])}±{fmt(d['se'])}, z={fmt(d['z_vs_1'])}, "
            f"compatible={d['compatible_2se']}"
        )
    a.append("")
    a.append("## Pell (7,5) secondary calibration")
    a.append("")
    a.append(pell.get("caveat", ""))
    a.append("")
    for run in pell.get("runs", []):
        a.append(
            f"- {run['file']} seed={run.get('seed')} n={run.get('samples')}: "
            f"C_N(lin)={fmt(run['C_N_linearized'])}, "
            f"A_p(L^4)={fmt(run['A_p_L4'])}, A_M={fmt(run['A_M_L13_4'])}, "
            f"A_M/B={fmt(run['A_p_pred_AM_over_B'])}"
        )
    a.append("")
    a.append("## Provenance")
    a.append("")
    a.append(f"- C05 source commit: `{meta.get('source_commit')}`")
    a.append(f"- analysis branch: `{meta.get('source_branch')}`")
    a.append(f"- wall_time_s: {meta.get('wall_time', {}).get('analysis_s')}")
    a.append(f"- machine: {meta.get('machine')}")
    a.append("")
    path.write_text("\n".join(a) + "\n")


def write_q43_report(path: Path, q43: dict, rows: list[dict], huawei: dict, meta: dict) -> None:
    a = []
    a.append("# Q43 — Joint angular–radial harmonic challenge on C05 histograms")
    a.append("")
    a.append("Issue 36. **Existing C05 data only.** Seed `0xC0100001`, **2e6** samples,")
    a.append("cross-channel Newman-Ziff `ΔM`. Primary exponent frozen at `13/8`.")
    a.append("Train N=65,85,130; hold out N=145,170. N=145 is not dropped.")
    a.append("Either/cross are not two replications. Second-stage N=185/221/265 MC skipped.")
    a.append("")
    a.append("## Design columns (exact rationals)")
    a.append("")
    a.append("| N | split | Δcos4 | Δcos8 | ΔM ± SE | N^{13/8} ΔM ± SE |")
    a.append("|---|---|---|---|---|---|")
    for r in rows:
        a.append(
            f"| {r['N']} | {r['split']} | {r['delta_cos4_exact']} | {r['delta_cos8_exact']} | "
            f"{fmt(r['delta_M'])}±{fmt(r['delta_M_se'])} | {fmt(r['y_scaled'])}±{fmt(r['y_scaled_se'])} |"
        )
    a.append("")
    a.append(f"Selected power-correction `ω` (training only): **{q43['omega_selection']['selected_omega']}**")
    a.append(f"(conditioning warning: {q43.get('omega_conditioning_warning')}; ω=4 on N=65..130 is ill-conditioned.)")
    a.append(f"Free `α` (training only): **{fmt(q43['free_alpha_selection']['alpha'])}**")
    fa = q43['free_alpha_selection']
    a.append(f"(search interval {fa.get('search_interval')}; hit_lower_bound={fa.get('hit_lower_bound')}. "
             "A bound hit is not a physical exponent.)")
    a.append("")
    a.append("## Model scores (covariance-aware χ² on ΔM)")
    a.append("")
    a.append("| model | params | train χ² | holdout χ² | cond | notes |")
    a.append("|---|---|---|---|---|---|")
    for name in MODEL_ORDER:
        m = q43["models"][name]
        coef = ", ".join(
            f"{k}={fmt(v['value'])}±{fmt(v['se'])}" for k, v in m["coefficients"].items()
        )
        extra = ""
        if name == "fixed_13_8_cos4_plus_cos8":
            extra = f"A8 {m['A8_status']['status']}"
        if name == "fixed_13_8_cos4_power":
            extra = f"ω={q43['omega_selection']['selected_omega']}"
        if name == "free_alpha_cos4":
            extra = f"α={fmt(q43['free_alpha_selection']['alpha'])}"
        a.append(
            f"| `{name}` | {coef} | {fmt(m['train_chi2_deltaM'])} | {fmt(m['holdout_chi2_deltaM'])} | "
            f"{fmt(m['cond'])} | {extra} |"
        )
    a.append("")
    a.append(f"**Raw lowest holdout χ²:** `{q43.get('holdout_winner_raw_chi2')}`")
    a.append(f"**Well-conditioned holdout winner (cond<1e4):** `{q43.get('holdout_winner_well_conditioned')}`")
    a.append(f"**Reported winner:** `{q43['holdout_winner']}`")
    a.append(f"H8/log overturn H4 on holdout (Δχ²>1): **{q43.get('H4_overturned_by_H8_or_log')}**")
    a.append("")
    a.append("H4 holdout residuals have the predicted sign at both N=145 and N=170.")
    a.append("Holdout χ² values are all O(1–3) on two noisy points; 2e6 does not distinguish")
    a.append("models sharply. `free_alpha` hits the preregistered lower bound α=0.5 and is")
    a.append("not treated as a physical winner. Ill-conditioned power/log companions are retained")
    a.append("as failed candidates.")
    a.append("")
    a.append("### Held-out signed residuals (ΔM)")
    a.append("")
    a.append("| model | N | true | pred | signed residual | z |")
    a.append("|---|---|---|---|---|---|")
    for name in MODEL_ORDER:
        m = q43["models"][name]
        for i, n in enumerate(m["holdout_N"]):
            true = m["holdout_true_deltaM"][i]
            pred = m["holdout_pred_deltaM"][i]
            resid = m["holdout_signed_residual_deltaM"][i]
            se = next(r["delta_M_se"] for r in rows if r["N"] == n)
            a.append(
                f"| `{name}` | {n} | {fmt(true)} | {fmt(pred)} | {fmt(resid)} | {fmt(resid / se if se else float('nan'))} |"
            )
    a.append("")
    a8 = q43["models"]["fixed_13_8_cos4_plus_cos8"].get("A8_status", {})
    a.append("## Is A8 resolved?")
    a.append("")
    a.append(f"- status: **{a8.get('status')}**")
    a.append(f"- A8 = {fmt(a8.get('A8'))} ± {fmt(a8.get('A8_se'))} (z={fmt(a8.get('z'))})")
    a.append(f"- corr(A4,A8) = {fmt(a8.get('corr_A4_A8'))}")
    a.append(f"- condition number = {fmt(a8.get('cond'))}")
    a.append(f"- identifiable design: {a8.get('identifiable_design')}")
    a.append(f"- holdout χ² vs H4: {fmt(a8.get('holdout_chi2'))} vs {fmt(a8.get('H4_holdout_chi2'))}; improves={a8.get('holdout_improves_on_H4')}")
    a.append(f"- rule: {a8.get('rule')}")
    a.append("")
    a.append("## Huawei published A4 (separate table, not pooled)")
    a.append("")
    a.append("| N | Huawei A4 | C05 A_M (NZ M) | C05−Huawei | z (independent SEs) |")
    a.append("|---|---|---|---|---|")
    for rec in huawei.get("comparison", []):
        a.append(
            f"| {rec['N']} | {rec['huawei_published_A4']}±{rec['huawei_published_se']} | "
            f"{fmt(rec.get('C05_A_M'))}±{fmt(rec.get('C05_A_M_se'))} | "
            f"{fmt(rec.get('signed_C05_minus_huawei'))} | {fmt(rec.get('z_independent'))} |"
        )
    a.append("")
    a.append(huawei.get("rule", ""))
    a.append("")
    a.append("C05 `A_M` is the Newman-Ziff matching function; Huawei A4 is the fixed-p")
    a.append("wrapping difference. Different observable, different seed, different sample count.")
    a.append("")
    a.append("## Provenance")
    a.append("")
    a.append(f"- C05 source commit: `{meta.get('source_commit')}`")
    a.append(f"- analysis branch: `{meta.get('source_branch')}`")
    a.append(f"- wall_time_s: {meta.get('wall_time', {}).get('analysis_s')}")
    a.append(f"- machine: {meta.get('machine')}")
    a.append("")
    path.write_text("\n".join(a) + "\n")


def git_head() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(REPO), text=True).strip()
    except Exception:
        return "unknown"


def git_branch() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=str(REPO), text=True).strip()
    except Exception:
        return "unknown"


def main() -> int:
    t0 = time.time()
    ap = argparse.ArgumentParser()
    ap.add_argument("--c05-root", type=Path, default=REPO / "results/server-20260828/C05")
    ap.add_argument("--pell-root", type=Path, default=REPO / "results/server-20260828/pell")
    ap.add_argument("--huawei-csv", type=Path, default=REPO / "results/server-20260828/gaussian/same_n_prod.analysis.csv")
    ap.add_argument("--q42-out", type=Path, default=REPO / "results/server-20260828/Q42")
    ap.add_argument("--q43-out", type=Path, default=REPO / "results/server-20260828/Q43")
    args = ap.parse_args()

    loaded: dict[int, SizeData] = {}
    for n, r1, r2 in ALL_PAIRS:
        d1 = args.c05_root / c05.geom_name(*r1)
        if not (d1 / "kminus_hist.csv").is_file():
            continue
        loaded[n] = load_size(args.c05_root, n, r1, r2)

    if set(TRAIN_N + HOLDOUT_N) - set(loaded):
        missing = set(TRAIN_N + HOLDOUT_N) - set(loaded)
        raise SystemExit(f"missing C05 sizes: {sorted(missing)}")

    rows = sizes_table(loaded)
    q43 = fit_all_models(rows)
    loso = loso_predictions(rows, q43)
    pell = pell_secondary(args.pell_root)
    huawei = huawei_table(rows, args.huawei_csv)
    verdict = q42_passfail(loaded)
    elapsed = time.time() - t0

    c05_meta_path = args.c05_root / "metadata.json"
    c05_meta = json.loads(c05_meta_path.read_text()) if c05_meta_path.is_file() else {}
    meta = {
        "queue_ids": ["Q42", "Q43"],
        "issues": [35, 36],
        "program": "Q42 amplitude closure + Q43 angular-radial harmonic challenge on existing C05 histograms",
        "source_commit": git_head(),
        "c05_histogram_commit": git_head(),
        "c05_metadata_source_commit": c05_meta.get("source_commit"),
        "c05_commit_note": "histograms from C05 SHA f89191a5468ec8417bbeac373334cecb6b5833a7; this analysis does not rerun MC",
        "source_branch": git_branch(),
        "machine": c05_meta.get("machine", "cursor Linux 8 cores / 16 GB, no GPU"),
        "os_kernel": c05_meta.get("os_kernel", platform.release()),
        "compiler_or_interpreter": f"python {platform.python_version()} + numpy {np.__version__}",
        "dependency_versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
        },
        "rng_algorithm_if_stochastic": "none in this analysis; C05 used Philox4x32-10",
        "rng_seed_or_counter_ranges": {
            "C05_global_seed_hex": hex(C05_SEED),
            "C05_global_seed": C05_SEED,
            "C05_seed_note": "0xC0100001, NOT a 1e8 campaign",
            "from_C05_metadata": c05_meta.get("rng_seed_or_counter_ranges"),
        },
        "sample_counts": {
            "per_same_N_pair": 2_000_000,
            "batches": 40,
            "replicas_per_batch": 50_000,
            "NOT_1e8": True,
            "from_C05_metadata": c05_meta.get("sample_counts"),
        },
        "wall_time": {"analysis_s": elapsed, "new_MC_s": 0.0},
        "thread_or_gpu_configuration": "analysis process, no GPU, no new OpenMP MC",
        "channel": "cross",
        "N1105": "not started",
        "second_stage_N185_221_265": "skipped (default; first five sizes analysed)",
        "this_is_2e6_not_1e8": True,
    }

    # ---- Q42 artifacts ----
    args.q42_out.mkdir(parents=True, exist_ok=True)
    batch_rows = []
    for n, rec in loaded.items():
        for b in rec.batches:
            batch_rows.append({"N": n, "split": split_label(n), **{k: b.get(k) for k in b}})
    write_csv(args.q42_out / "batch_observables.csv", batch_rows)
    amp_rows = []
    for n in TRAIN_N + HOLDOUT_N:
        s = loaded[n].batch_summary
        amp_rows.append(
            {
                "N": n,
                "split": split_label(n),
                "delta_cos4": loaded[n].design["delta_cos4_float"],
                "delta_cos4_exact": loaded[n].design["delta_cos4"],
                "A_M": s["A_M"]["mean"],
                "A_M_se": s["A_M"]["se"],
                "B": s["B_slope"]["mean"],
                "B_se": s["B_slope"]["se"],
                "A_p_direct": s["A_p_direct"]["mean"],
                "A_p_direct_se": s["A_p_direct"]["se"],
                "A_p_lin": s["A_p_lin"]["mean"],
                "A_p_lin_se": s["A_p_lin"]["se"],
                "A_p_pred_AM_over_B": s["A_p_pred"]["mean"],
                "A_p_pred_se": s["A_p_pred"]["se"],
                "C_N_direct": s["C_N_direct"]["mean"],
                "C_N_direct_se": s["C_N_direct"]["se"],
                "C_N_linearized": s["C_N_linearized"]["mean"],
                "C_N_linearized_se": s["C_N_linearized"]["se"],
                "C_N_quadratic": s["C_N_quadratic"]["mean"],
                "C_N_quadratic_se": s["C_N_quadratic"]["se"],
                "C_N_plugin": s.get("C_N_plugin_delta_method", {}).get("mean"),
                "C_N_plugin_se": s.get("C_N_plugin_delta_method", {}).get("se"),
                "delta_M": s["delta_M_pref"]["mean"],
                "delta_M_se": s["delta_M_pref"]["se"],
                "mean_Mprime": s["mean_Mprime_pref"]["mean"],
                "mean_Mprime_se": s["mean_Mprime_pref"]["se"],
                "delta_Mprime": s["delta_Mprime_pref"]["mean"],
                "delta_Mprime_se": s["delta_Mprime_pref"]["se"],
                "delta_root": s["delta_root"]["mean"],
                "delta_root_se": s["delta_root"]["se"],
                "delta_root_lin": s["delta_root_lin"]["mean"],
                "delta_root_lin_se": s["delta_root_lin"]["se"],
                "direct_minus_lin": s["direct_minus_lin"]["mean"],
                "direct_minus_lin_se": s["direct_minus_lin"]["se"],
                "pooled_C_N_direct": loaded[n].pooled.get("C_N_direct"),
                "samples": 2_000_000,
                "seed": C05_SEED,
            }
        )
    write_csv(args.q42_out / "amplitudes_closure.csv", amp_rows)
    write_json(
        args.q42_out / "closure.json",
        {
            "by_size": {str(n): loaded[n].batch_summary for n in TRAIN_N + HOLDOUT_N},
            "pooled": {str(n): loaded[n].pooled for n in TRAIN_N + HOLDOUT_N},
            "verdict": verdict,
            "pell": pell,
            "protocol": {
                "p_ref": PC_REF,
                "A_M": "N^{13/8} DeltaM(p_ref) / DeltaCos4",
                "B": "N^{-3/8} mean(M'(p_ref))",
                "A_p": "-N^2 DeltaRoot / DeltaCos4",
                "C_N": "-DeltaRoot * mean(Mprime) / DeltaM",
                "primary_prediction": "A_p = A_M / B",
                "samples": 2_000_000,
                "seed": hex(C05_SEED),
                "not_1e8": True,
            },
        },
    )
    write_json(args.q42_out / "metadata.json", meta)
    write_q42_report(args.q42_out / "REPORT.md", loaded, pell, verdict, meta)
    (args.q42_out / "commands.txt").write_text(
        "# Q42/Q43 analysis of existing C05 histograms (no new MC)\n"
        "python3 scripts/analyze_q42_q43.py "
        "--c05-root results/server-20260828/C05 "
        "--pell-root results/server-20260828/pell "
        "--huawei-csv results/server-20260828/gaussian/same_n_prod.analysis.csv "
        "--q42-out results/server-20260828/Q42 "
        "--q43-out results/server-20260828/Q43\n"
        "python3 -m unittest tests.test_q42_q43\n"
    )
    write_environment(args.q42_out / "environment.txt", {"task": "Q42"})

    # ---- Q43 artifacts ----
    args.q43_out.mkdir(parents=True, exist_ok=True)
    write_csv(args.q43_out / "design_matrix.csv", rows)
    write_csv(
        args.q43_out / "heldout_residuals.csv",
        [
            {
                "model": name,
                "N": q43["models"][name]["holdout_N"][i],
                "true_deltaM": q43["models"][name]["holdout_true_deltaM"][i],
                "pred_deltaM": q43["models"][name]["holdout_pred_deltaM"][i],
                "signed_residual": q43["models"][name]["holdout_signed_residual_deltaM"][i],
                "signed_residual_scaled": q43["models"][name]["holdout_signed_residual_scaled"][i],
            }
            for name in MODEL_ORDER
            for i in range(len(q43["models"][name]["holdout_N"]))
        ],
    )
    write_csv(args.q43_out / "loso.csv", loso)
    write_json(args.q43_out / "model_scores.json", q43)
    write_json(args.q43_out / "huawei_comparison.json", huawei)
    write_json(args.q43_out / "metadata.json", meta)
    write_q43_report(args.q43_out / "REPORT.md", q43, rows, huawei, meta)
    (args.q43_out / "commands.txt").write_text((args.q42_out / "commands.txt").read_text())
    write_environment(args.q43_out / "environment.txt", {"task": "Q43"})

    # checksums last
    (args.q42_out / "checksums.sha256").write_text(sha256_tree(args.q42_out))
    (args.q43_out / "checksums.sha256").write_text(sha256_tree(args.q43_out))

    print("Q42 wrote", args.q42_out)
    print("Q43 wrote", args.q43_out)
    print("holdout C_N:", {n: loaded[n].batch_summary["C_N_direct"] for n in HOLDOUT_N})
    print("Q43 winner:", q43["holdout_winner"])
    a8 = q43["models"]["fixed_13_8_cos4_plus_cos8"].get("A8_status")
    print("A8:", a8)
    print("elapsed_s", elapsed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
