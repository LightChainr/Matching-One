#!/usr/bin/env python3
"""Reconstruct M, S, D, P4 and thermal even/odd projectors from C05 rank histograms.

Training protocol (frozen before held-out inspection in this script's control flow):
  - u-grid is chosen from N=65 orientation-averaged Mbar only
  - scaling models are fit on N in TRAIN_N
  - N in HOLDOUT_N are scored after the fit is frozen

CROSS channel only. Either-wrap rows are not a second replication of M.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
import sys
from pathlib import Path

import mpmath as mp

PC_REF = 0.59274605079210
TRAIN_N = (65, 85, 130)
HOLDOUT_N = (145, 170)
ALL_PAIRS = (
    (65, (8, 1), (7, 4)),
    (85, (9, 2), (7, 6)),
    (130, (11, 3), (9, 7)),
    (145, (12, 1), (9, 8)),
    (170, (13, 1), (11, 7)),
)


def geom_name(a: int, b: int) -> str:
    return f"g_{a}_{b}"


def cos4_exact(a: int, b: int) -> float:
    a2 = a * a
    b2 = b * b
    return (a2 * a2 - 6 * a2 * b2 + b2 * b2) / (a2 + b2) ** 2


def load_hist(path: Path) -> list[int]:
    rows = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            rows.append((int(row["k"]), int(row["count"])))
    rows.sort()
    n_plus_1 = rows[-1][0]
    hist = [0] * (n_plus_1 + 1)
    for k, c in rows:
        hist[k] = c
    return hist


def load_batch_hist(path: Path) -> dict[int, list[int]]:
    by_batch: dict[int, list[tuple[int, int]]] = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            b = int(row["batch"])
            by_batch.setdefault(b, []).append((int(row["k"]), int(row["count"])))
    out = {}
    for b, pairs in by_batch.items():
        pairs.sort()
        hist = [0] * (pairs[-1][0] + 1)
        for k, c in pairs:
            hist[k] = c
        out[b] = hist
    return out


def binomial_weights(n: int, p: float, cutoff: float = 1e-20):
    if p <= 0.0:
        return 0, [1.0], 1.0
    if p >= 1.0:
        return n, [1.0], 1.0
    mode = int((n + 1) * p)
    mode = min(max(mode, 0), n)
    ratio_up = p / (1.0 - p)
    ratio_dn = (1.0 - p) / p
    up = [1.0]
    acc = 1.0
    k = mode
    while k < n:
        k += 1
        acc *= (n - k + 1) / k * ratio_up
        if acc < cutoff and k - mode > 64:
            break
        if acc == 0.0:
            break
        up.append(acc)
    dn = []
    acc = 1.0
    k = mode
    while k > 0:
        acc *= k / (n - k + 1) * ratio_dn
        k -= 1
        if acc < cutoff and mode - k > 64:
            break
        if acc == 0.0:
            break
        dn.append(acc)
    kmin = mode - len(dn)
    w = list(reversed(dn)) + up
    tot = math.fsum(w)
    return kmin, w, tot


def convolve(qk: list[float], n: int, p: float) -> float:
    kmin, w, tot = binomial_weights(n, p)
    s = 0.0
    for i, wt in enumerate(w):
        idx = kmin + i
        if 0 <= idx <= n:
            s += wt * qk[idx]
    return s / tot


def cdf_from_hist(hist: list[int]) -> list[float]:
    total = float(sum(hist))
    out = [0.0] * len(hist)
    acc = 0.0
    for k, c in enumerate(hist):
        acc += c
        out[k] = acc / total if total else 0.0
    return out


def wrap_from_kplus(kp_hist: list[int], n: int) -> list[float]:
    """P(K_plus <= k) for k=0..n (occupation of primal)."""
    cdf = cdf_from_hist(kp_hist)
    return [cdf[k] for k in range(n + 1)]


def wrap_hat_from_kminus(km_hist: list[int], n: int) -> list[float]:
    """P(matching wraps | matching occupation m) = P(K_minus > n-m)."""
    total = float(sum(km_hist))
    # sf[k] = P(K_minus > k)
    sf = [0.0] * (n + 2)
    acc = total
    for k in range(n + 2):
        acc -= km_hist[k] if k < len(km_hist) else 0
        sf[k] = acc / total if total else 0.0
    q = [0.0] * (n + 1)
    for m in range(n + 1):
        q[m] = sf[n - m]
    return q


def mean_d_from_hists(km_hist: list[int], kp_hist: list[int], n: int) -> list[float]:
    rg = wrap_from_kplus(kp_hist, n)
    # M occupation-k: P(K+<=k) - P(K- > k)
    total = float(sum(km_hist))
    sf = [0.0] * (n + 2)
    acc = total
    for k in range(n + 2):
        acc -= km_hist[k] if k < len(km_hist) else 0
        sf[k] = acc / total if total else 0.0
    return [rg[k] - sf[k] for k in range(n + 1)]


def reconstruct_curves(km_hist: list[int], kp_hist: list[int], n: int, p_grid: list[float]):
    qg = wrap_from_kplus(kp_hist, n)
    qh = wrap_hat_from_kminus(km_hist, n)
    qd = mean_d_from_hists(km_hist, kp_hist, n)
    rows = []
    for p in p_grid:
        rg = convolve(qg, n, p)
        rh = convolve(qh, n, p)
        rh_comp = convolve(qh, n, 1.0 - p)
        m = convolve(qd, n, p)
        rows.append(
            {
                "p": p,
                "R_G": rg,
                "R_hat": rh,
                "R_hat_comp": rh_comp,
                "M": m,
                "S": 0.5 * (rg + rh),
                "D": 0.5 * (rg - rh),
            }
        )
    return rows


def find_root(f, lo=0.4, hi=0.8):
    flo, fhi = f(lo), f(hi)
    if flo == 0.0:
        return lo
    if fhi == 0.0:
        return hi
    if flo * fhi > 0.0:
        grid = [lo + i * (hi - lo) / 400.0 for i in range(401)]
        vals = [f(x) for x in grid]
        found = None
        for i in range(len(grid) - 1):
            if vals[i] == 0.0:
                return grid[i]
            if vals[i] * vals[i + 1] <= 0.0:
                found = (grid[i], grid[i + 1])
                break
        if found is None:
            return float("nan")
        lo, hi = found
        flo, fhi = f(lo), f(hi)
    a, b, fa, fb = lo, hi, flo, fhi
    for _ in range(80):
        m = 0.5 * (a + b)
        fm = f(m)
        if fm == 0.0 or (b - a) < 1e-15:
            return m
        if fa * fm <= 0.0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


def solve_mbar_equals(mbar, target, lo=0.05, hi=0.95):
    return find_root(lambda p: mbar(p) - target, lo, hi)


def finite_diff(f, x, h=1e-4, order=1):
    if order == 1:
        return (f(x + h) - f(x - h)) / (2.0 * h)
    if order == 2:
        return (f(x + h) - 2.0 * f(x) + f(x - h)) / (h * h)
    if order == 3:
        return (f(x + 2 * h) - 2.0 * f(x + h) + 2.0 * f(x - h) - f(x - 2 * h)) / (2.0 * h**3)
    if order == 5:
        return (
            f(x + 3 * h)
            - 4.0 * f(x + 2 * h)
            + 5.0 * f(x + h)
            - 5.0 * f(x - h)
            + 4.0 * f(x - 2 * h)
            - f(x - 3 * h)
        ) / (2.0 * h**5)
    raise ValueError(order)


def beta_pdf(p: float, k: int, n: int) -> float:
    if k < 1 or k > n or p <= 0.0 or p >= 1.0:
        return 0.0
    # Beta(k, n+1-k)
    logB = math.lgamma(k) + math.lgamma(n + 1 - k) - math.lgamma(n + 1)
    return math.exp((k - 1) * math.log(p) + (n - k) * math.log(1.0 - p) - logB)


def mixture_density(km_hist, kp_hist, n, p):
    total = float(sum(km_hist))
    if total == 0:
        return 0.0
    s = 0.0
    for k in range(1, n + 1):
        w = (km_hist[k] + kp_hist[k]) / (2.0 * total)
        if w:
            s += w * beta_pdf(p, k, n)
    return s


def hist_moments(hist: list[int]):
    total = float(sum(hist))
    mean = sum(k * c for k, c in enumerate(hist)) / total
    var = sum((k - mean) ** 2 * c for k, c in enumerate(hist)) / total
    m3 = sum((k - mean) ** 3 * c for k, c in enumerate(hist)) / total
    m4 = sum((k - mean) ** 4 * c for k, c in enumerate(hist)) / total
    sd = math.sqrt(var) if var > 0 else float("nan")
    skew = m3 / sd**3 if sd == sd and sd > 0 else float("nan")
    kurt = m4 / sd**4 if sd == sd and sd > 0 else float("nan")
    return {"mean": mean, "var": var, "sd": sd, "skew": skew, "kurtosis": kurt, "n": total}


def quantile_from_hist(hist: list[int], q: float) -> float:
    total = float(sum(hist))
    acc = 0.0
    for k, c in enumerate(hist):
        acc += c
        if acc / total >= q:
            return float(k)
    return float(len(hist) - 1)


def p4(x1: float, x2: float, c1: float, c2: float) -> float:
    den = c1 - c2
    if abs(den) < 1e-18:
        return float("nan")
    return (x1 - x2) / den


def load_orientation(root: Path, a: int, b: int):
    name = geom_name(a, b)
    d = root / name
    km = load_hist(d / "kminus_hist.csv")
    kp = load_hist(d / "kplus_hist.csv")
    n = len(km) - 2  # bins 0..N+1
    meta = {}
    mp_ = d / "run_meta.json"
    if mp_.is_file():
        meta = json.loads(mp_.read_text())
        n = int(meta.get("N", n))
    samples = int(meta.get("samples", sum(km)))
    return {
        "name": name,
        "a": a,
        "b": b,
        "n": n,
        "theta": math.atan2(b, a),
        "cos4": cos4_exact(a, b),
        "km": km,
        "kp": kp,
        "samples": samples,
        "meta": meta,
        "batch_km": load_batch_hist(d / "batch_kminus_hist.csv") if (d / "batch_kminus_hist.csv").is_file() else {},
        "batch_kp": load_batch_hist(d / "batch_kplus_hist.csv") if (d / "batch_kplus_hist.csv").is_file() else {},
    }


def jackknife_p4(g1, g2, p: float, field: str):
    batches = sorted(set(g1["batch_km"]) & set(g2["batch_km"]))
    if len(batches) < 2:
        return float("nan"), float("nan")
    vals = []
    n = g1["n"]
    c1, c2 = g1["cos4"], g2["cos4"]
    for b in batches:
        km1, kp1 = g1["batch_km"][b], g1["batch_kp"][b]
        km2, kp2 = g2["batch_km"][b], g2["batch_kp"][b]
        curve1 = reconstruct_curves(km1, kp1, n, [p])[0]
        curve2 = reconstruct_curves(km2, kp2, n, [p])[0]
        vals.append(p4(curve1[field], curve2[field], c1, c2))
    mu = statistics.fmean(vals)
    if len(vals) < 2:
        return mu, 0.0
    var = statistics.variance(vals)
    se = math.sqrt(var / len(vals))
    return mu, se


def lin_log_fit(xs: list[float], ys: list[float]):
    """Fit log|y| = a + b log x; return (amp_sign, exponent, predictions)."""
    lx = [math.log(x) for x in xs]
    signs = [1.0 if y >= 0 else -1.0 for y in ys]
    ly = [math.log(abs(y)) if y != 0 else float("-inf") for y in ys]
    if any(not math.isfinite(v) for v in ly):
        return {"ok": False}
    n = len(xs)
    mx = sum(lx) / n
    my = sum(ly) / n
    sxx = sum((x - mx) ** 2 for x in lx)
    sxy = sum((lx[i] - mx) * (ly[i] - my) for i in range(n))
    b = sxy / sxx if sxx else float("nan")
    a = my - b * mx
    pred = [signs[i] * math.exp(a + b * math.log(xs[i])) for i in range(n)]
    return {"ok": True, "log_amp": a, "exponent": b, "pred": pred, "sign": signs[0]}


def rmse(y, yhat):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(y, yhat)) / len(y))


def analyze(root: Path, outdir: Path):
    mp.mp.dps = 40
    p_grid = [i / 200.0 for i in range(1, 200)]
    available = []
    for n, r1, r2 in ALL_PAIRS:
        d1 = root / geom_name(*r1)
        d2 = root / geom_name(*r2)
        if d1.is_dir() and d2.is_dir() and (d1 / "kminus_hist.csv").is_file():
            available.append((n, r1, r2))

    geoms = {}
    for n, r1, r2 in available:
        geoms[n] = (load_orientation(root, *r1), load_orientation(root, *r2))

    # --- freeze u from N=65 only ---
    frozen = {"u_grid": [0.0], "u_source": None, "p_star_train": None}
    if 65 in geoms:
        g1, g2 = geoms[65]
        n = g1["n"]
        qd1 = mean_d_from_hists(g1["km"], g1["kp"], n)
        qd2 = mean_d_from_hists(g2["km"], g2["kp"], n)

        def mbar65(p):
            return 0.5 * (convolve(qd1, n, p) + convolve(qd2, n, p))

        p_star = find_root(mbar65, 0.45, 0.75)
        frozen["p_star_train"] = p_star
        frozen["u_source"] = "N=65 Mbar_N root and |Mbar| coverage"
        # probe reachable |Mbar|
        span = min(abs(mbar65(0.2)), abs(mbar65(0.85)))
        candidates = [0.0, 0.05, 0.1, 0.2, 0.4]
        u_grid = [u for u in candidates if u == 0.0 or u < 0.8 * span]
        if 0.0 not in u_grid:
            u_grid.insert(0, 0.0)
        frozen["u_grid"] = u_grid
        frozen["mbar_span_at_edges"] = span
    else:
        frozen["u_grid"] = [0.0]
        frozen["note"] = "N=65 missing; u-grid fallback {0}"

    derived = {
        "convention": {
            "K_plus": "smallest black occupation with primal CROSS wrap; N+1 if never",
            "K_minus": "N-m_star+1, m_star first matching CROSS wrap on reverse perm; 0 if never",
            "beta": "T|K=k ~ Beta(k,N+1-k) for k=1..N",
            "channel": "cross only; either is not a second M replication",
        },
        "frozen": frozen,
        "train_N": [n for n in TRAIN_N if n in geoms],
        "holdout_N": [n for n in HOLDOUT_N if n in geoms],
        "sizes": {},
        "thermal": {},
        "scaling": {},
    }

    curves_out = []
    thermal_rows = []

    for n, r1, r2 in available:
        g1, g2 = geoms[n]
        c1, c2 = g1["cos4"], g2["cos4"]
        L = math.sqrt(n)
        qd1 = mean_d_from_hists(g1["km"], g1["kp"], g1["n"])
        qd2 = mean_d_from_hists(g2["km"], g2["kp"], g2["n"])
        qg1 = wrap_from_kplus(g1["kp"], g1["n"])
        qg2 = wrap_from_kplus(g2["kp"], g2["n"])
        qh1 = wrap_hat_from_kminus(g1["km"], g1["n"])
        qh2 = wrap_hat_from_kminus(g2["km"], g2["n"])

        def m1(p, qd=qd1, nn=g1["n"]):
            return convolve(qd, nn, p)

        def m2(p, qd=qd2, nn=g2["n"]):
            return convolve(qd, nn, p)

        def mbar(p):
            return 0.5 * (m1(p) + m2(p))

        p_star = find_root(mbar, 0.45, 0.75)
        h = max(1e-4, 0.25 / n)
        mprime = finite_diff(mbar, p_star, h=h, order=1)
        m3 = finite_diff(mbar, p_star, h=h, order=3)
        m5 = finite_diff(mbar, p_star, h=max(h, 2e-4), order=5)
        kappa3 = m3 / mprime**3 if mprime else float("nan")
        kappa5 = m5 / mprime**5 if mprime else float("nan")
        dens = mixture_density(g1["km"], g1["kp"], g1["n"], p_star)

        def field_at(g, qg, qh, qd, p):
            rg = convolve(qg, g["n"], p)
            rh = convolve(qh, g["n"], p)
            return {
                "R_G": rg,
                "R_hat": rh,
                "M": convolve(qd, g["n"], p),
                "S": 0.5 * (rg + rh),
                "D": 0.5 * (rg - rh),
            }

        # pooled curves at p_star
        f1s = field_at(g1, qg1, qh1, qd1, p_star)
        f2s = field_at(g2, qg2, qh2, qd2, p_star)
        p4_at_star = {k: p4(f1s[k], f2s[k], c1, c2) for k in ("R_G", "R_hat", "M", "S", "D")}
        p4_se = {}
        for k in ("S", "D", "M"):
            mu, se = jackknife_p4(g1, g2, p_star, k)
            p4_se[k] = {"batch_mean": mu, "batch_se": se}

        km_mom = hist_moments(g1["km"])
        # mixture of two orientations' ranks for gap
        gap_summary = {
            "mean_Kminus_1": hist_moments(g1["km"])["mean"],
            "mean_Kplus_1": hist_moments(g1["kp"])["mean"],
            "mean_Kminus_2": hist_moments(g2["km"])["mean"],
            "mean_Kplus_2": hist_moments(g2["kp"])["mean"],
            "mean_gap_1": hist_moments(g1["kp"])["mean"] - hist_moments(g1["km"])["mean"],
            "mean_gap_2": hist_moments(g2["kp"])["mean"] - hist_moments(g2["km"])["mean"],
        }

        size_row = {
            "N": n,
            "L_phys": L,
            "reps": [list(r1), list(r2)],
            "cos4": [c1, c2],
            "delta_cos4": c1 - c2,
            "samples": [g1["samples"], g2["samples"]],
            "p_star_Mbar": p_star,
            "Mbar_at_pref": mbar(PC_REF),
            "Mprime_at_pstar": mprime,
            "density_Mprime_over_2": 0.5 * mprime,
            "beta_mixture_density_at_pstar": dens,
            "kappa3_at_pstar": kappa3,
            "kappa5_at_pstar": kappa5,
            "P4_at_pstar": p4_at_star,
            "P4_batch": p4_se,
            "rank_moments": gap_summary,
            "quantiles_Kplus_1": {
                "q05": quantile_from_hist(g1["kp"], 0.05),
                "q50": quantile_from_hist(g1["kp"], 0.5),
                "q95": quantile_from_hist(g1["kp"], 0.95),
            },
            "median_root_orientation": [find_root(m1, 0.45, 0.75), find_root(m2, 0.45, 0.75)],
        }
        derived["sizes"][str(n)] = size_row

        # thermal even/odd for frozen u
        trow = {"N": n, "L_phys": L, "u": []}
        for u in frozen["u_grid"]:
            if u == 0.0:
                pp = pm = p_star
            else:
                pm = solve_mbar_equals(mbar, -u, 0.05, p_star if p_star == p_star else 0.6)
                pp = solve_mbar_equals(mbar, +u, p_star if p_star == p_star else 0.6, 0.95)
            entry = {"u": u, "p_minus": pm, "p_plus": pp, "fields": {}}
            for key in ("R_G", "R_hat", "M", "S", "D"):
                if u == 0.0:
                    x4p = p4_at_star[key]
                    x4m = x4p
                else:
                    f1p = field_at(g1, qg1, qh1, qd1, pp)
                    f2p = field_at(g2, qg2, qh2, qd2, pp)
                    f1m = field_at(g1, qg1, qh1, qd1, pm)
                    f2m = field_at(g2, qg2, qh2, qd2, pm)
                    x4p = p4(f1p[key], f2p[key], c1, c2)
                    x4m = p4(f1m[key], f2m[key], c1, c2)
                entry["fields"][key] = {
                    "P4_plus": x4p,
                    "P4_minus": x4m,
                    "thermal_even": 0.5 * (x4p + x4m),
                    "thermal_odd": 0.5 * (x4p - x4m),
                }
            trow["u"].append(entry)
        thermal_rows.append(trow)
        derived["thermal"][str(n)] = trow

        for p in (0.5, PC_REF, p_star):
            f1 = field_at(g1, qg1, qh1, qd1, p)
            f2 = field_at(g2, qg2, qh2, qd2, p)
            curves_out.append(
                {
                    "N": n,
                    "p": p,
                    "tag": "pstar" if p == p_star else ("pref" if p == PC_REF else "p50"),
                    **{f"1_{k}": f1[k] for k in f1},
                    **{f"2_{k}": f2[k] for k in f2},
                    **{f"P4_{k}": p4(f1[k], f2[k], c1, c2) for k in f1},
                }
            )

    # --- scaling, fit on TRAIN_N, score HOLDOUT_N ---
    def collect_even(key, u=0.0):
        xs, ys, ns = [], [], []
        for rec in thermal_rows:
            n = rec["N"]
            for ent in rec["u"]:
                if abs(ent["u"] - u) < 1e-15:
                    xs.append(math.sqrt(n))
                    ys.append(ent["fields"][key]["thermal_even"])
                    ns.append(n)
        return ns, xs, ys

    scaling = {}
    for key in ("D", "S", "M"):
        ns, xs, ys = collect_even(key, 0.0)
        train_mask = [n in TRAIN_N for n in ns]
        hold_mask = [n in HOLDOUT_N for n in ns]
        n_tr = [n for n, m in zip(ns, train_mask) if m]
        x_tr = [x for x, m in zip(xs, train_mask) if m]
        y_tr = [y for y, m in zip(ys, train_mask) if m]
        n_ho = [n for n, m in zip(ns, hold_mask) if m]
        x_ho = [x for x, m in zip(xs, hold_mask) if m]
        y_ho = [y for y, m in zip(ys, hold_mask) if m]

        def model_pred(exponent, x_fit, y_fit, x_pred):
            # y ~ A * L^exponent, A from mean y L^{-exponent} on train
            if not x_fit:
                return []
            amps = [y_fit[i] / (x_fit[i] ** exponent) for i in range(len(x_fit))]
            A = sum(amps) / len(amps)
            return A, [A * (x ** exponent) for x in x_pred]

        entry = {
            "N": ns,
            "L": xs,
            "P4_thermal_even_u0": ys,
            "signed_train": y_tr,
            "signed_holdout": y_ho,
        }
        if len(x_tr) >= 1:
            for name, exp in (("L^-13/4", -13.0 / 4.0), ("L^-2", -2.0), ("L^-1", -1.0)):
                A, pred_tr = model_pred(exp, x_tr, y_tr, x_tr)
                _, pred_ho = model_pred(exp, x_tr, y_tr, x_ho) if x_ho else (A, [])
                entry[name] = {
                    "A_train": A,
                    "train_rmse": rmse(y_tr, pred_tr) if pred_tr else None,
                    "holdout_pred": pred_ho,
                    "holdout_true": y_ho,
                    "holdout_signed_error": [y_ho[i] - pred_ho[i] for i in range(len(y_ho))]
                    if pred_ho
                    else [],
                    "holdout_rmse": rmse(y_ho, pred_ho) if pred_ho else None,
                }
            # log companion: y ~ A L^e (1 + B log L) with e=-13/4, fit B on train
            e0 = -13.0 / 4.0
            if len(x_tr) >= 2:
                # y / L^e = A + A B log L  => linear in log L
                z = [y_tr[i] / (x_tr[i] ** e0) for i in range(len(x_tr))]
                lx = [math.log(x) for x in x_tr]
                mx, mz = sum(lx) / len(lx), sum(z) / len(z)
                sxx = sum((x - mx) ** 2 for x in lx)
                sxy = sum((lx[i] - mx) * (z[i] - mz) for i in range(len(lx)))
                slope = sxy / sxx if sxx else 0.0
                intercept = mz - slope * mx  # A
                B = slope / intercept if intercept else float("nan")

                def pred_log(xv):
                    return [intercept * (1.0 + B * math.log(x)) * (x ** e0) for x in xv]

                ptr = pred_log(x_tr)
                pho = pred_log(x_ho) if x_ho else []
                entry["L^-13/4_(1+B_logL)"] = {
                    "A_train": intercept,
                    "B_train": B,
                    "train_rmse": rmse(y_tr, ptr),
                    "holdout_pred": pho,
                    "holdout_true": y_ho,
                    "holdout_signed_error": [y_ho[i] - pho[i] for i in range(len(y_ho))] if pho else [],
                    "holdout_rmse": rmse(y_ho, pho) if pho else None,
                }
            free = lin_log_fit(x_tr, y_tr) if len(x_tr) >= 2 else {"ok": False}
            if free.get("ok"):
                Afree = math.exp(free["log_amp"]) * (1.0 if y_tr[0] >= 0 else -1.0)
                # actually sign is in data; reconstruct A from first point
                efree = free["exponent"]
                A, pred_tr = model_pred(efree, x_tr, y_tr, x_tr)
                _, pred_ho = model_pred(efree, x_tr, y_tr, x_ho) if x_ho else (A, [])
                entry["free_exponent"] = {
                    "exponent_train": efree,
                    "A_train": A,
                    "train_rmse": rmse(y_tr, pred_tr),
                    "holdout_pred": pred_ho,
                    "holdout_true": y_ho,
                    "holdout_signed_error": [y_ho[i] - pred_ho[i] for i in range(len(y_ho))]
                    if pred_ho
                    else [],
                    "holdout_rmse": rmse(y_ho, pred_ho) if pred_ho else None,
                }
        # pairwise effective exponent between consecutive sizes
        eff = []
        for i in range(len(xs) - 1):
            if ys[i] == 0 or ys[i + 1] == 0 or ys[i] * ys[i + 1] <= 0:
                weff = float("nan")
            else:
                weff = math.log(abs(ys[i + 1] / ys[i])) / math.log(xs[i + 1] / xs[i])
            eff.append({"N_pair": [ns[i], ns[i + 1]], "w_eff": weff, "y": [ys[i], ys[i + 1]]})
        entry["pairwise_effective_exponent"] = eff
        scaling[key] = entry

    derived["scaling"] = scaling
    # target-test summaries
    d0 = scaling.get("D", {})
    s0 = scaling.get("S", {})
    derived["target_tests"] = {
        "P4D_thermal_even_u0_like_L_m13_4": {
            "pairwise_w_eff": d0.get("pairwise_effective_exponent"),
            "target": -13.0 / 4.0,
            "train_rmse_fixed": (d0.get("L^-13/4") or {}).get("train_rmse"),
            "holdout_rmse_fixed": (d0.get("L^-13/4") or {}).get("holdout_rmse"),
            "holdout_rmse_log": (d0.get("L^-13/4_(1+B_logL)") or {}).get("holdout_rmse"),
            "holdout_rmse_free": (d0.get("free_exponent") or {}).get("holdout_rmse"),
            "free_exponent_train": (d0.get("free_exponent") or {}).get("exponent_train"),
            "note": "Underpowered if fewer than three training sizes; negative results preserved.",
        },
        "P4S_thermal_even_u0_like_L_m2": {
            "pairwise_w_eff": s0.get("pairwise_effective_exponent"),
            "target": -2.0,
            "train_rmse_L_m2": (s0.get("L^-2") or {}).get("train_rmse"),
            "holdout_rmse_L_m2": (s0.get("L^-2") or {}).get("holdout_rmse"),
            "holdout_rmse_L_m13_4": (s0.get("L^-13/4") or {}).get("holdout_rmse"),
        },
        "log_alternative_improves_holdout_P4D": None,
    }
    h_fixed = (d0.get("L^-13/4") or {}).get("holdout_rmse")
    h_log = (d0.get("L^-13/4_(1+B_logL)") or {}).get("holdout_rmse")
    if h_fixed is not None and h_log is not None:
        derived["target_tests"]["log_alternative_improves_holdout_P4D"] = h_log < h_fixed
        derived["target_tests"]["holdout_rmse_delta_log_minus_fixed"] = h_log - h_fixed

    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "derived_summary.json").write_text(json.dumps(derived, indent=2, default=str) + "\n")
    with (outdir / "curves_at_anchors.csv").open("w", newline="") as f:
        if curves_out:
            w = csv.DictWriter(f, fieldnames=list(curves_out[0].keys()))
            w.writeheader()
            w.writerows(curves_out)

    # combined pooled hists copy pointers
    combined_km = []
    combined_kp = []
    for n, r1, r2 in available:
        g1, g2 = geoms[n]
        for g in (g1, g2):
            for k, c in enumerate(g["km"]):
                combined_km.append({"geom": g["name"], "N": n, "k": k, "count": c})
            for k, c in enumerate(g["kp"]):
                combined_kp.append({"geom": g["name"], "N": n, "k": k, "count": c})
    with (outdir / "kminus_hist.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["geom", "N", "k", "count"])
        w.writeheader()
        w.writerows(combined_km)
    with (outdir / "kplus_hist.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["geom", "N", "k", "count"])
        w.writeheader()
        w.writerows(combined_kp)

    return derived


def write_report(derived: dict, outdir: Path, extra: dict):
    lines = []
    a = lines.append
    a("# C05 / P33 — Thermal-coordinate tomography from threshold ranks")
    a("")
    a("Engine: issue-9 Philox Fisher–Yates Newman–Ziff on the C00 `HomologyUnionFind`")
    a("(exact `adj(P)/det(P)` windings). Channel: **cross** only. Either-wrap is an")
    a("exact-test diagnostic, not a second matching-function replication.")
    a("")
    a("## Off-by-one convention")
    a("")
    conv = derived["convention"]
    for k, v in conv.items():
        a(f"- `{k}`: {v}")
    a("")
    a("## Exact tests")
    a("")
    a(f"- overall: **{extra.get('exact_status', 'n/a')}**")
    a("")
    a("## Frozen analysis choices (from N=65 only)")
    a("")
    a("```json")
    a(json.dumps(derived.get("frozen", {}), indent=2, default=str))
    a("```")
    a("")
    a("Training sizes for scaling: "
      + ", ".join(str(n) for n in derived.get("train_N", [])))
    a("Held-out sizes: " + ", ".join(str(n) for n in derived.get("holdout_N", [])))
    a("")
    a("## Per-size summary")
    a("")
    a("| N | L | samples | p*_Mbar | P4[D](p*) | P4[S](p*) | kappa3 | mean gap |")
    a("|---|---|---|---|---|---|---|---|")
    for n in sorted((int(k) for k in derived.get("sizes", {})), key=int):
        s = derived["sizes"][str(n)]
        samp = s["samples"][0]
        p4d = s["P4_at_pstar"]["D"]
        p4s = s["P4_at_pstar"]["S"]
        gap = 0.5 * (s["rank_moments"]["mean_gap_1"] + s["rank_moments"]["mean_gap_2"])
        a(
            f"| {n} | {s['L_phys']:.4f} | {samp} | {s['p_star_Mbar']:.8f} | "
            f"{p4d:.6e} | {p4s:.6e} | {s['kappa3_at_pstar']:.4f} | {gap:.3f} |"
        )
    a("")
    a("Batch jackknife (signed P4, not p-values):")
    a("")
    for n in sorted((int(k) for k in derived.get("sizes", {})), key=int):
        s = derived["sizes"][str(n)]
        bd = s["P4_batch"].get("D", {})
        bs = s["P4_batch"].get("S", {})
        a(
            f"- N={n}: P4[D] batch mean {bd.get('batch_mean')} se {bd.get('batch_se')}; "
            f"P4[S] batch mean {bs.get('batch_mean')} se {bs.get('batch_se')}"
        )
    a("")
    a("## Target tests")
    a("")
    a("P4[D] thermal-even at u=0 vs L^{-13/4}≈L^{-3.25}:")
    a("")
    a("```json")
    a(json.dumps(derived.get("target_tests", {}).get("P4D_thermal_even_u0_like_L_m13_4"), indent=2, default=str))
    a("```")
    a("")
    a("P4[S] thermal-even at u=0 vs L^{-2}:")
    a("")
    a("```json")
    a(json.dumps(derived.get("target_tests", {}).get("P4S_thermal_even_u0_like_L_m2"), indent=2, default=str))
    a("```")
    a("")
    a(
        "Log alternative improves held-out P4[D]: "
        + str(derived.get("target_tests", {}).get("log_alternative_improves_holdout_P4D"))
    )
    a("")
    a("## Remaining")
    a("")
    for line in extra.get("remaining", []):
        a(f"- {line}")
    a("")
    (outdir / "REPORT.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("results/server-20260828/C05"))
    ap.add_argument("--exact-json", type=Path, default=None)
    args = ap.parse_args()
    extra = {"remaining": []}
    ej = args.exact_json or (args.root / "exact" / "exact_tests.json")
    if ej.is_file():
        payload = json.loads(ej.read_text())
        extra["exact_status"] = payload.get("overall")
    derived = analyze(args.root, args.root)
    write_report(derived, args.root, extra)
    print("wrote", args.root / "derived_summary.json")
    print("wrote", args.root / "REPORT.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
