#!/usr/bin/env python3
"""Analyze C01 same-N Gaussian orientation batch aggregates.

Consumes *.t*.moments.json / *.batches.csv produced by
``src/gaussian_orientation_mc.cpp``.  Replica-level first and second moments
are recovered from the stored 20-channel sum and Gram matrix; batch jackknife
standard errors are the reported uncertainties.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


CHANNEL_NAMES = ("cross", "both", "either", "direction_0", "direction_1")
VAR_NAMES = [
    f"o{ori}_{graph}_{ch}"
    for ori in (1, 2)
    for graph in ("primal", "matching")
    for ch in CHANNEL_NAMES
]
NVAR = 20


def cos4(a: int, b: int) -> float:
    a2 = a * a
    b2 = b * b
    return (a2 * a2 - 6 * a2 * b2 + b2 * b2) / (a2 + b2) ** 2


def theta_deg(a: int, b: int) -> float:
    return math.degrees(math.atan2(b, a))


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def combine_moments(moments: dict[str, Any]) -> dict[str, Any]:
    batches = moments["batches"]
    n_batch = len(batches)
    if n_batch < 2:
        raise ValueError("need at least two batches")
    total_n = 0
    total_sum = [0] * NVAR
    total_gram = [0] * (NVAR * NVAR)
    batch_means: list[list[float]] = []
    batch_n: list[int] = []
    for row in batches:
        n = int(row["samples"])
        s = [int(v) for v in row["sum"]]
        g = [int(v) for v in row["gram"]]
        if len(s) != NVAR or len(g) != NVAR * NVAR:
            raise ValueError("moment vector width mismatch")
        total_n += n
        for i in range(NVAR):
            total_sum[i] += s[i]
        for i in range(NVAR * NVAR):
            total_gram[i] += g[i]
        batch_n.append(n)
        batch_means.append([s[i] / n for i in range(NVAR)])
    mean = [total_sum[i] / total_n for i in range(NVAR)]
    cov = [[0.0] * NVAR for _ in range(NVAR)]
    for i in range(NVAR):
        for j in range(NVAR):
            exy = total_gram[i * NVAR + j] / total_n
            cov[i][j] = exy - mean[i] * mean[j]
    return {
        "n": total_n,
        "n_batch": n_batch,
        "mean": mean,
        "cov": cov,
        "batch_means": batch_means,
        "batch_n": batch_n,
        "t": int(moments["t"]),
        "variables": moments.get("variables", VAR_NAMES),
    }


def batch_se(values: list[float]) -> float:
    n = len(values)
    if n < 2:
        return float("nan")
    mean = math.fsum(values) / n
    var = math.fsum((v - mean) ** 2 for v in values) / (n - 1)
    return math.sqrt(var / n)


def jackknife_se(values: list[float]) -> float:
    n = len(values)
    total = math.fsum(values)
    loo = [(total - v) / (n - 1) for v in values]
    center = math.fsum(loo) / n
    return math.sqrt((n - 1) / n * math.fsum((v - center) ** 2 for v in loo))


def lincomb_batch_values(batch_means: list[list[float]], weights: list[float]) -> list[float]:
    out = []
    for row in batch_means:
        out.append(math.fsum(w * row[i] for i, w in enumerate(weights) if w))
    return out


def replica_var(cov: list[list[float]], weights: list[float]) -> float:
    acc = 0.0
    for i, wi in enumerate(weights):
        if wi == 0:
            continue
        for j, wj in enumerate(weights):
            if wj == 0:
                continue
            acc += wi * wj * cov[i][j]
    return acc


def weights_for(kind: str, channel: str) -> list[float]:
    c = CHANNEL_NAMES.index(channel)
    w = [0.0] * NVAR
    # layout: o1 primal 0-4, o1 matching 5-9, o2 primal 10-14, o2 matching 15-19
    i1p, i1m, i2p, i2m = c, 5 + c, 10 + c, 15 + c
    if kind == "primal_delta":
        w[i1p] = 1.0
        w[i2p] = -1.0
    elif kind == "matching_delta":
        w[i1m] = 1.0
        w[i2m] = -1.0
    elif kind == "M_delta":
        # D_N = M(theta1)-M(theta2), M = primal - matching
        w[i1p] = 1.0
        w[i1m] = -1.0
        w[i2p] = -1.0
        w[i2m] = 1.0
    elif kind == "S_delta":
        # matching-even orientation difference
        w[i1p] = 0.5
        w[i1m] = 0.5
        w[i2p] = -0.5
        w[i2m] = -0.5
    elif kind == "o1_primal":
        w[i1p] = 1.0
    elif kind == "o1_matching":
        w[i1m] = 1.0
    elif kind == "o2_primal":
        w[i2p] = 1.0
    elif kind == "o2_matching":
        w[i2m] = 1.0
    elif kind == "o1_M":
        w[i1p] = 1.0
        w[i1m] = -1.0
    elif kind == "o2_M":
        w[i2p] = 1.0
        w[i2m] = -1.0
    elif kind == "o1_S":
        w[i1p] = 0.5
        w[i1m] = 0.5
    elif kind == "o2_S":
        w[i2p] = 0.5
        w[i2m] = 0.5
    else:
        raise ValueError(kind)
    return w


def estimate(combined: dict[str, Any], kind: str, channel: str) -> dict[str, float]:
    w = weights_for(kind, channel)
    batch_vals = lincomb_batch_values(combined["batch_means"], w)
    mean = math.fsum(m * n for m, n in zip(
        [math.fsum(ww * combined["mean"][i] for i, ww in enumerate(w) if ww)],
        [1],
    ))
    # Direct population mean from combined.mean
    mean = math.fsum(ww * combined["mean"][i] for i, ww in enumerate(w))
    se = batch_se(batch_vals)
    var = replica_var(combined["cov"], w)
    n = combined["n"]
    replica_se = math.sqrt(max(var, 0.0) / n) if n else float("nan")
    z = mean / se if se else float("nan")
    return {
        "mean": mean,
        "batch_se": se,
        "jackknife_se": jackknife_se(batch_vals),
        "replica_var": var,
        "replica_se": replica_se,
        "z_batch": z,
        "n": float(n),
        "n_batch": float(combined["n_batch"]),
    }


def correlation(combined: dict[str, Any], i: int, j: int) -> float:
    cii = combined["cov"][i][i]
    cjj = combined["cov"][j][j]
    if cii <= 0 or cjj <= 0:
        return float("nan")
    return combined["cov"][i][j] / math.sqrt(cii * cjj)


def analyze_run(moments_path: Path, metadata_path: Path) -> dict[str, Any]:
    moments = load_json(moments_path)
    metadata = load_json(metadata_path)
    combined = combine_moments(moments)
    a1, b1 = metadata["rep1"]
    a2, b2 = metadata["rep2"]
    n_sites = int(metadata["N"])
    dcos = cos4(a1, b1) - cos4(a2, b2)
    n13 = n_sites ** (13 / 8)
    payload: dict[str, Any] = {
        "moments_path": str(moments_path),
        "metadata_path": str(metadata_path),
        "mode": metadata["mode"],
        "N": n_sites,
        "rep1": [a1, b1],
        "rep2": [a2, b2],
        "theta1_deg": theta_deg(a1, b1),
        "theta2_deg": theta_deg(a2, b2),
        "cos4_1": cos4(a1, b1),
        "cos4_2": cos4(a2, b2),
        "delta_cos4": dcos,
        "t": combined["t"],
        "p": metadata["p"],
        "samples": combined["n"],
        "batches": combined["n_batch"],
        "seed": metadata["seed"],
        "replica_begin": metadata["replica_begin"],
        "replica_end": metadata["replica_end"],
        "elapsed_seconds": metadata.get("elapsed_seconds"),
        "threads": metadata.get("threads"),
        "channels": {},
        "means": {},
        "correlations": {},
    }
    for ch in CHANNEL_NAMES:
        ch_block = {}
        for kind in (
            "o1_primal", "o1_matching", "o2_primal", "o2_matching",
            "o1_M", "o2_M", "o1_S", "o2_S",
            "primal_delta", "matching_delta", "M_delta", "S_delta",
        ):
            ch_block[kind] = estimate(combined, kind, ch)
        mdelta = ch_block["M_delta"]
        sdelta = ch_block["S_delta"]
        if abs(dcos) < 1e-18:
            ch_block["A4_M"] = {"mean": None, "batch_se": None}
            ch_block["A4_S_times_N"] = {"mean": None, "batch_se": None}
        else:
            ch_block["A4_M"] = {
                "mean": n13 * mdelta["mean"] / dcos,
                "batch_se": abs(n13 / dcos) * mdelta["batch_se"],
            }
            ch_block["A4_S_times_N"] = {
                "mean": n_sites * sdelta["mean"] / dcos,
                "batch_se": abs(n_sites / dcos) * sdelta["batch_se"],
            }
        payload["channels"][ch] = ch_block
        payload["means"][ch] = {
            "o1_primal": ch_block["o1_primal"]["mean"],
            "o1_matching": ch_block["o1_matching"]["mean"],
            "o2_primal": ch_block["o2_primal"]["mean"],
            "o2_matching": ch_block["o2_matching"]["mean"],
        }
        # CRN correlation of primal wrapping and of matching-even S
        i1p = CHANNEL_NAMES.index(ch)
        payload["correlations"][ch] = {
            "primal_o1_o2": correlation(combined, i1p, 10 + i1p),
            "matching_o1_o2": correlation(combined, 5 + i1p, 15 + i1p),
        }
    payload["covariance_20"] = combined["cov"]
    payload["variable_means_20"] = combined["mean"]
    payload["variable_names"] = VAR_NAMES
    return payload


def write_long_form_rows(analysis: dict[str, Any], stage: str) -> list[dict[str, Any]]:
    rows = []
    for ch, block in analysis["channels"].items():
        for kind in (
            "o1_primal", "o1_matching", "o2_primal", "o2_matching",
            "o1_M", "o2_M", "o1_S", "o2_S",
            "primal_delta", "matching_delta", "M_delta", "S_delta",
        ):
            est = block[kind]
            rows.append({
                "stage": stage,
                "mode": analysis["mode"],
                "N": analysis["N"],
                "a1": analysis["rep1"][0],
                "b1": analysis["rep1"][1],
                "a2": analysis["rep2"][0],
                "b2": analysis["rep2"][1],
                "t": analysis["t"],
                "p": analysis["p"],
                "channel": ch,
                "observable": kind,
                "mean": est["mean"],
                "batch_se": est["batch_se"],
                "z_batch": est["z_batch"],
                "samples": analysis["samples"],
                "batches": analysis["batches"],
            })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--moments", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    payload = analyze_run(args.moments, args.metadata)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    either = payload["channels"]["either"]
    a4 = either["A4_M"]["mean"]
    a4se = either["A4_M"]["batch_se"]
    a4txt = "n/a" if a4 is None else f"{a4:+.6g} +/- {a4se:.3g}"
    print(
        f"N={payload['N']} t={payload['t']} mode={payload['mode']} "
        f"D_M_either={either['M_delta']['mean']:+.6g} +/- {either['M_delta']['batch_se']:.3g}  "
        f"D_S_either={either['S_delta']['mean']:+.6g} +/- {either['S_delta']['batch_se']:.3g}  "
        f"A4_M={a4txt}"
    )
    print(f"wrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
