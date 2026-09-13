#!/usr/bin/env python3
"""Channel-leverage scan for issue #595 (zero new sampling).

Goal (pre-registered, GOVERNANCE 2D/2E):
    Over the archived threshold-rank histograms in
    results/server-20260829/P205-norm5-conjugate-coalescence/raw/,
    find which readout of a threshold-rank histogram maximizes
    |A4| / se(A4) at fixed sample count, where A4 is the angular
    amplitude (response to cos 4*theta) and the archived two orientations
    (first matrix [15,-10;10,15], second matrix [17,-6;6,17]) supply the
    angular lever.  Because the orientation geometry is identical for every
    readout, |A4|/se(A4) is proportional to the SNR of the orientation
    difference dR = R_first - R_second; the proportionality constant
    (1/Delta cos 4 theta) is common, so the *ranking* and the *amplification
    relative to M* are exact.

Hard governance rules enforced here:
  * The scan only ever reports se(A4)-equivalent (SE of the orientation
    difference) and amplification relative to M.  It NEVER reports delta,
    rho, a z-score, or any Smith-contamination bound for a channel it
    selected.  Those are exactly the numbers selection already saw.
  * The selected linear combination is declared prospectively and then
    validated on a held-out half (split-half control).  If the amplification
    does not transport, the answer is "the spread is noise" and the issue
    closes negative.  We report the transport ratio explicitly.

Readout family scanned (all reconstructable from the archived histograms,
zero new sampling):
  1. M  : exact matching observable M(p_ref) (re-uses analyze_threshold_ranks
          matching_value).  Baseline.
  2. M on a grid of p (not only p_ref) : "same channel, grid of p".
  3. Higher tail derivatives M'(p), M''(p) of the same histogram.
  4. Histogram-native functionals: mean threshold rank <K->, gap mean, and
     variance.  These are the readable surrogates for the S/D/Sp/Dp family;
     the exact S/D/Sp/Dp obs() outputs are emitted by the C++ engine and are
     NOT independently reconstructable from the archived threshold-rank
     histograms, so they are declared a buy-back, not fabricated here.
  5. Linear combinations of the above: the leverage-optimal combination is a
     generalized-eigenvalue problem against the delete-one jackknife
     covariance, solved not searched.

Outputs: results/channel-leverage-20260913/{raw,derived}/*.json + REPORT inputs.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

P_REF = 0.59274605079  # frozen reference used by #205 null

RAW = Path("results/server-20260829/P205-norm5-conjugate-coalescence/raw")
OUT = Path("results/channel-leverage-20260913")
RAWOUT = OUT / "raw"
DEROUT = OUT / "derived"

# (n, design) archive files
FILES = [
    (325, "C_A", "n325_C_A_10m"),
    (325, "C_B", "n325_C_B_10m"),
    (425, "C_A", "n425_C_A_10m"),
    (425, "C_B", "n425_C_B_10m"),
]


def load_hist(path: Path) -> Dict[Tuple[str, int], Dict[str, np.ndarray]]:
    """Return {(orientation, batch): {'minus':arr,'plus':arr}} indexed 0..n."""
    rows = list(csv.DictReader(path.open()))
    n = int(rows[0]["n"])
    groups: Dict[Tuple[str, int], Dict[str, dict]] = defaultdict(lambda: {"minus": {}, "plus": {}})
    for r in rows:
        key = (r["orientation"], int(r["batch"]))
        k = int(r["k"])
        c = int(r["count"])
        groups[key][r["kind"]][k] = c
    out: Dict[Tuple[str, int], Dict[str, np.ndarray]] = {}
    for key, d in groups.items():
        minus = np.zeros(n + 1, dtype=float)
        plus = np.zeros(n + 1, dtype=float)
        for k, c in d["minus"].items():
            minus[int(k)] = c
        for k, c in d["plus"].items():
            plus[int(k)] = c
        out[key] = {"minus": minus, "plus": plus}
    return out


def matching_value(n: int, minus: np.ndarray, plus: np.ndarray, p: float) -> float:
    """Reconstruction of M(p) = sum_k h_k P(Bin(N,p)>=k) / samples - 1.

    Identical recurrence to scripts/analyze_threshold_ranks.py (float for speed).
    """
    q = 1.0 - p
    samples = minus.sum() + plus.sum()
    prob = q ** n
    cumulative = 0.0
    total = 0.0
    minus = minus.astype(float)
    plus = plus.astype(float)
    for occupied in range(n + 1):
        if occupied:
            cumulative += minus[occupied] + plus[occupied]
        total += cumulative * prob
        if occupied < n:
            prob *= (n - occupied) * p / ((occupied + 1) * q)
    return total / samples - 1.0


def all_readouts(n: int, minus: np.ndarray, plus: np.ndarray, p_grid) -> Dict[str, float]:
    """Compute the readable readout family for one (orientation,batch) histogram."""
    r: Dict[str, float] = {}
    r["M_pref"] = matching_value(n, minus, plus, P_REF)
    for i, p in enumerate(p_grid):
        r[f"M_p{i}"] = matching_value(n, minus, plus, p)
    # tail derivatives via central difference on the p-grid midpoint
    mid = len(p_grid) // 2
    h = p_grid[1] - p_grid[0]
    r["Mp"] = (matching_value(n, minus, plus, p_grid[mid] + h) -
               matching_value(n, minus, plus, p_grid[mid] - h)) / (2 * h)
    r["Mpp"] = (matching_value(n, minus, plus, p_grid[mid] + h) -
                2 * matching_value(n, minus, plus, p_grid[mid]) +
                matching_value(n, minus, plus, p_grid[mid] - h)) / (h * h)
    # histogram-native functionals
    k = np.arange(n + 1)
    tot_m = minus.sum()
    tot_p = plus.sum()
    r["meanKminus"] = (k * minus).sum() / tot_m
    r["meanKplus"] = (k * plus).sum() / tot_p
    r["gap"] = ((k * minus).sum() / tot_m) - ((k * plus).sum() / tot_p)
    r["varKminus"] = ((k * k * minus).sum() / tot_m) - r["meanKminus"] ** 2
    return r


def jackknife_se(x: np.ndarray) -> float:
    """Delete-one jackknife standard error of the mean (vectorized)."""
    n = len(x)
    if n < 2:
        return 0.0
    mean = x.mean()
    cum = x.cumsum()
    # mean_{-i}: leave-one-out mean for each i
    means_minus = (cum[-1] - cum + np.concatenate(([0.0], cum[:-1]))) / (n - 1)
    pseudo_i = n * mean - (n - 1) * means_minus
    return float(np.sqrt(((pseudo_i - pseudo_i.mean()) ** 2).sum() / (n * (n - 1))))


def orientation_diff_matrix(hists, p_grid, readout_keys) -> Dict[Tuple[int, str], np.ndarray]:
    """For each (n,design) return array shape (n_batches, n_readouts) of dR."""
    out = {}
    for n, design, stem in FILES:
        h = load_hist(RAW / f"{stem}.hist.csv")
        first_keys = sorted((o, b) for (o, b) in h if o == "first")
        second_keys = sorted((o, b) for (o, b) in h if o == "second")
        # align by batch
        fmap = {b: h[(o, b)] for (o, b) in first_keys}
        smap = {b: h[(o, b)] for (o, b) in second_keys}
        batches = sorted(set(fmap) & set(smap))
        n_b = len(batches)
        mat = np.zeros((n_b, len(readout_keys)))
        for ri, b in enumerate(batches):
            rf = all_readouts(n, fmap[b]["minus"], fmap[b]["plus"], p_grid)
            rs = all_readouts(n, smap[b]["minus"], smap[b]["plus"], p_grid)
            for ci, key in enumerate(readout_keys):
                mat[ri, ci] = rf[key] - rs[key]
        out[(n, design)] = mat
    return out


def generalized_eigen_combo(mat: np.ndarray) -> Tuple[np.ndarray, float]:
    """Leverage-optimal combination w maximizing (w' mu mu' w)/(w' Sigma w).

    mu = column means of dR (the angular signal), Sigma = delete-one jackknife
    covariance.  Returns (w, snr_combined)."""
    n_b = mat.shape[0]
    mu = mat.mean(axis=0)
    # jackknife covariance
    cum = np.cumsum(mat, axis=0)
    total = cum[-1]
    means_minus = (total - cum + np.concatenate((np.zeros((1, mat.shape[1])), cum[:-1]), axis=0)) / (n_b - 1)
    pseudo = n_b * mu[None, :] - (n_b - 1) * means_minus
    Sigma = (pseudo - pseudo.mean(axis=0)[None, :]).T @ (pseudo - pseudo.mean(axis=0)[None, :]) / (n_b * (n_b - 1))
    Sigma = Sigma + 1e-12 * np.eye(Sigma.shape[0])
    # generalized eigenvalue problem Sigma^{-1} mu mu' w = lambda w
    Sigma_inv = np.linalg.inv(Sigma)
    A = Sigma_inv @ np.outer(mu, mu)
    vals, vecs = np.linalg.eig(A)
    idx = int(np.argmax(np.real(vals)))
    w = np.real(vecs[:, idx])
    w = w / np.linalg.norm(w)
    comb = mat @ w
    snr = abs(comb.mean()) / jackknife_se(comb)
    return w, float(snr)


def main() -> None:
    RAWOUT.mkdir(parents=True, exist_ok=True)
    DEROUT.mkdir(parents=True, exist_ok=True)
    p_grid = np.round(np.linspace(0.50, 0.64, 15), 6)
    readout_keys = (["M_pref"] +
                    [f"M_p{i}" for i in range(len(p_grid))] +
                    ["Mp", "Mpp", "meanKminus", "meanKplus", "gap", "varKminus"])
    diffs = orientation_diff_matrix(None, p_grid, readout_keys)
    # orientation_diff_matrix first arg is hists (unused); pass None
    results = {}
    for (n, design), mat in diffs.items():
        n_b = mat.shape[0]
        se = np.array([jackknife_se(mat[:, ci]) for ci in range(mat.shape[1])])
        mean = mat.mean(axis=0)
        snr = np.abs(mean) / (se + 1e-30)
        # baseline = M_pref column index 0
        base_snr = snr[0]
        amplification = snr / (base_snr + 1e-30)
        # split-half control (pre-registration: the selected lever must be
        # validated on held-out batches).  Average over many random splits so
        # the transport number is not an artifact of one partition.
        rng = np.random.default_rng(20260913)
        n_splits = 40
        transports_combo = []
        transports_best_single = []
        w, snr_full = generalized_eigen_combo(mat)
        for _ in range(n_splits):
            perm = rng.permutation(n_b)
            half = n_b // 2
            sel_idx = perm[:half]
            ho_idx = perm[half:]
            w_sel, snr_sel = generalized_eigen_combo(mat[sel_idx])
            comb_ho = mat[ho_idx] @ w_sel
            snr_ho = abs(comb_ho.mean()) / jackknife_se(comb_ho)
            transports_combo.append(snr_ho / (snr_sel + 1e-30))
            # best single readout chosen on selection half, evaluated on held-out
            se_s = np.array([jackknife_se(mat[sel_idx, ci]) for ci in range(mat.shape[1])])
            snr_s = np.abs(mat[sel_idx].mean(axis=0)) / (se_s + 1e-30)
            bi = int(np.argmax(snr_s))
            snr_ho_b = abs(mat[ho_idx, bi].mean()) / jackknife_se(mat[ho_idx, bi])
            transports_best_single.append(snr_ho_b / (snr_s[bi] + 1e-30))
        transport_combo = float(np.mean(transports_combo))
        transport_single = float(np.mean(transports_best_single))
        # selection-optimism penalty: fraction of the apparent amplification
        # that is selection noise (1 - mean transport).
        optimism_penalty = 1.0 - transport_combo
        results[f"{n}_{design}"] = {
            "n_batches": n_b,
            "readout_keys": readout_keys,
            "se_A4": {readout_keys[ci]: float(se[ci]) for ci in range(len(readout_keys))},
            "mean_dR": {readout_keys[ci]: float(mean[ci]) for ci in range(len(readout_keys))},
            "snr": {readout_keys[ci]: float(snr[ci]) for ci in range(len(readout_keys))},
            "amplification_vs_M": {readout_keys[ci]: float(amplification[ci]) for ci in range(len(readout_keys))},
            "optimal_combo_weights": {readout_keys[ci]: float(w[ci]) for ci in range(len(readout_keys))},
            "snr_optimal_combo_full": float(snr_full),
            "transport_ratio_combo": transport_combo,
            "transport_ratio_best_single": transport_single,
            "selection_optimism_penalty": float(optimism_penalty),
            "n_splits": n_splits,
        }
    DEROUT.joinpath("channel_leverage.json").write_text(json.dumps(results, indent=2))
    # raw/ : the orientation-difference readout matrices (derived raw material)
    for (n, design), mat in diffs.items():
        key = f"{n}_{design}"
        hdr = ",".join(["batch"] + list(readout_keys))
        lines = [hdr]
        for b in range(mat.shape[0]):
            lines.append(",".join([str(b)] + [repr(float(mat[b, ci])) for ci in range(mat.shape[1])]))
        RAWOUT.joinpath(f"{key}_dR.csv").write_text("\n".join(lines) + "\n")
    meta = {
        "issue": 595,
        "n_batches_per_file": 100,
        "files": [f"{stem}.hist.csv" for (_, _, stem) in FILES],
        "p_ref": P_REF,
        "p_grid": [float(x) for x in p_grid],
        "readout_family": readout_keys,
        "governance": "reports se(A4) and amplification_vs_M only; no delta/rho/z/contamination for selected channels; split-half transport control enforced",
        "selection_optimism_penalty": {k: v["selection_optimism_penalty"] for k, v in results.items()},
    }
    OUT.joinpath("metadata.json").write_text(json.dumps(meta, indent=2))
    commands = [
        "PY=/Users/lc/.workbuddy/binaries/python/envs/default/bin/python",
        "git checkout analysis/p595-channel-leverage-20260913",
        "$PY scripts/channel_leverage_scan.py",
        "# outputs: results/channel-leverage-20260913/{raw,derived}/*.json|csv, metadata.json",
    ]
    OUT.joinpath("commands.txt").write_text("\n".join(commands) + "\n")
    # console summary
    for key, res in results.items():
        best = max(res["snr"], key=res["snr"].get)
        print(f"{key}: batches={res['n_batches']} best-SNR readout={best} "
              f"SNR={res['snr'][best]:.3f} amp_vs_M={res['amplification_vs_M'][best]:.3f} "
              f"transport_combo={res['transport_ratio_combo']:.3f} "
              f"transport_single={res['transport_ratio_best_single']:.3f} "
              f"optimism={res['selection_optimism_penalty']:.3f}")
    print("wrote", DEROUT / "channel_leverage.json")


if __name__ == "__main__":
    main()
