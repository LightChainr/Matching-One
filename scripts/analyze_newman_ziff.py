#!/usr/bin/env python3
"""Canonical binomial convolution and matched-root analysis for Issue #9.

Q(p) = sum_k Binomial(N,k) p^k (1-p)^{N-k} Q_k
computed with a mode-centered recurrence (no huge binomial coefficients).

D_L^x(p) = R_G^x(p) - R_G*^x(1-p)
M_L(p)   = N_L(p) - Nhat_L(1-p) - N (p - 2 p^2 + p^4)

Roots are reported to >=18 significant digits. Uncertainty is batch SD / SE,
not print width. Reverse permutation is matching complement, not antithetic.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import struct
import sys
from collections import defaultdict

import mpmath as mp

P_LO = 0.590
P_HI = 0.595
PC_REF = 0.59274605079210


def read_batch_bin(path: str):
    with open(path, "rb") as f:
        magic = f.read(4)
        if magic != b"NZB1":
            raise ValueError(f"bad magic in {path}: {magic!r}")
        L, n_batches, n_k = struct.unpack("<III", f.read(12))
        (rpb,) = struct.unpack("<Q", f.read(8))
        rec = struct.Struct("<10Q")
        batches = []
        for _ in range(n_batches):
            cl_g = [0.0] * n_k
            cl_gs = [0.0] * n_k
            wraps = {name: [0.0] * n_k for name in ("H", "V", "E", "B")}
            wraps_s = {name: [0.0] * n_k for name in ("H", "V", "E", "B")}
            for k in range(n_k):
                row = rec.unpack(f.read(rec.size))
                cl_g[k] = row[0] / rpb
                cl_gs[k] = row[1] / rpb
                wraps["H"][k] = row[2] / rpb
                wraps["V"][k] = row[3] / rpb
                wraps["E"][k] = row[4] / rpb
                wraps["B"][k] = row[5] / rpb
                wraps_s["H"][k] = row[6] / rpb
                wraps_s["V"][k] = row[7] / rpb
                wraps_s["E"][k] = row[8] / rpb
                wraps_s["B"][k] = row[9] / rpb
            batches.append(
                {
                    "replicas": rpb,
                    "cl_g": cl_g,
                    "cl_gs": cl_gs,
                    "wrap_g": wraps,
                    "wrap_gs": wraps_s,
                }
            )
        return L, n_k - 1, rpb, batches


def read_pooled_csv(path: str):
    rows = []
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append(row)
    n = int(rows[-1]["k"])
    replicas = int(rows[0]["replicas"])
    cl_g = [0.0] * (n + 1)
    cl_gs = [0.0] * (n + 1)
    wrap_g = {name: [0.0] * (n + 1) for name in ("H", "V", "E", "B")}
    wrap_gs = {name: [0.0] * (n + 1) for name in ("H", "V", "E", "B")}
    stats = {
        "cl_g2": [0.0] * (n + 1),
        "cl_gs2": [0.0] * (n + 1),
        "cl_xy": [0.0] * (n + 1),
        "h_xy": [0.0] * (n + 1),
        "v_xy": [0.0] * (n + 1),
        "e_xy": [0.0] * (n + 1),
        "b_xy": [0.0] * (n + 1),
        "h_g": [0.0] * (n + 1),
        "h_gs": [0.0] * (n + 1),
        "v_g": [0.0] * (n + 1),
        "v_gs": [0.0] * (n + 1),
        "e_g": [0.0] * (n + 1),
        "e_gs": [0.0] * (n + 1),
        "b_g": [0.0] * (n + 1),
        "b_gs": [0.0] * (n + 1),
    }
    for row in rows:
        k = int(row["k"])
        cl_g[k] = float(row["sum_clusters_G"]) / replicas
        cl_gs[k] = float(row["sum_clusters_Gstar"]) / replicas
        wrap_g["H"][k] = float(row["sum_wrap_H_G"]) / replicas
        wrap_g["V"][k] = float(row["sum_wrap_V_G"]) / replicas
        wrap_g["E"][k] = float(row["sum_wrap_E_G"]) / replicas
        wrap_g["B"][k] = float(row["sum_wrap_B_G"]) / replicas
        wrap_gs["H"][k] = float(row["sum_wrap_H_Gstar"]) / replicas
        wrap_gs["V"][k] = float(row["sum_wrap_V_Gstar"]) / replicas
        wrap_gs["E"][k] = float(row["sum_wrap_E_Gstar"]) / replicas
        wrap_gs["B"][k] = float(row["sum_wrap_B_Gstar"]) / replicas
        stats["cl_g2"][k] = float(row["sum_clusters_G_sq"]) / replicas
        stats["cl_gs2"][k] = float(row["sum_clusters_Gstar_sq"]) / replicas
        stats["cl_xy"][k] = float(row["sum_clusters_G_Gstar"]) / replicas
        stats["h_xy"][k] = float(row["sum_H_G_H_Gstar"]) / replicas
        stats["v_xy"][k] = float(row["sum_V_G_V_Gstar"]) / replicas
        stats["e_xy"][k] = float(row["sum_E_G_E_Gstar"]) / replicas
        stats["b_xy"][k] = float(row["sum_B_G_B_Gstar"]) / replicas
        stats["h_g"][k] = wrap_g["H"][k]
        stats["h_gs"][k] = wrap_gs["H"][k]
        stats["v_g"][k] = wrap_g["V"][k]
        stats["v_gs"][k] = wrap_gs["V"][k]
        stats["e_g"][k] = wrap_g["E"][k]
        stats["e_gs"][k] = wrap_gs["E"][k]
        stats["b_g"][k] = wrap_g["B"][k]
        stats["b_gs"][k] = wrap_gs["B"][k]
    return n, replicas, cl_g, cl_gs, wrap_g, wrap_gs, stats


def binomial_weights(n: int, p: float, cutoff: float = 1e-20):
    if p <= 0.0:
        return 0, [1.0], 1.0
    if p >= 1.0:
        return n, [1.0], 1.0
    mode = int((n + 1) * p)
    if mode < 0:
        mode = 0
    if mode > n:
        mode = n
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


def convolve(qk, n, p):
    kmin, w, tot = binomial_weights(n, p)
    s = 0.0
    for i, wt in enumerate(w):
        s += wt * qk[kmin + i]
    return s / tot


def chi_square(p: float) -> float:
    return p - 2.0 * p * p + p * p * p * p


def make_D_wrap(qg, qgs, n):
    def D(p: float) -> float:
        return convolve(qg, n, p) - convolve(qgs, n, p)

    return D


def make_M(cl_g, cl_gs, n):
    def M(p: float) -> float:
        return convolve(cl_g, n, p) - convolve(cl_gs, n, p) - n * chi_square(p)

    return M


def find_root(f, lo=P_LO, hi=P_HI):
    flo = f(lo)
    fhi = f(hi)
    used_lo, used_hi = lo, hi
    if flo == 0.0:
        return lo, "endpoint", (lo, hi)
    if fhi == 0.0:
        return hi, "endpoint", (lo, hi)
    if flo * fhi > 0.0:
        expand_lo, expand_hi = 0.58, 0.61
        flo2 = f(expand_lo)
        fhi2 = f(expand_hi)
        if flo2 * fhi2 > 0.0:
            grid = [0.58 + i * 0.001 for i in range(31)]
            vals = [f(x) for x in grid]
            found = None
            for i in range(len(grid) - 1):
                if vals[i] == 0.0:
                    found = (grid[i], grid[i])
                    break
                if vals[i] * vals[i + 1] <= 0.0:
                    found = (grid[i], grid[i + 1])
                    break
            if found is None:
                return float("nan"), "no_sign_change", (lo, hi)
            used_lo, used_hi = found
            flo = f(used_lo)
            fhi = f(used_hi)
        else:
            used_lo, used_hi = expand_lo, expand_hi
            flo, fhi = flo2, fhi2
    a, b = used_lo, used_hi
    fa, fb = flo, fhi
    if fa * fb > 0.0:
        return float("nan"), "no_sign_change", (a, b)
    for _ in range(80):
        m = 0.5 * (a + b)
        fm = f(m)
        if fm == 0.0 or (b - a) < 1e-18:
            a = b = m
            break
        if fa * fm <= 0.0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    root = 0.5 * (a + b)
    status = "ok"
    if root < P_LO or root > P_HI:
        status = "outside_requested_interval"
    return root, status, (used_lo, used_hi)


def fmt18(x: float) -> str:
    if x != x:
        return "nan"
    return mp.nstr(mp.mpf(x), 18, strip_zeros=False)


def mean_sd_se(xs):
    xs = [x for x in xs if x == x]
    m = len(xs)
    if m == 0:
        return float("nan"), float("nan"), float("nan")
    mu = math.fsum(xs) / m
    if m == 1:
        return mu, 0.0, 0.0
    var = math.fsum((x - mu) ** 2 for x in xs) / (m - 1)
    sd = math.sqrt(var)
    se = sd / math.sqrt(m)
    return mu, sd, se


def var_diff(ex, ey, ex2, ey2, exy):
    return ex2 + ey2 - 2.0 * exy - (ex - ey) ** 2


def load_perf(path):
    if not os.path.isfile(path):
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def write_report(indir, seq_rows, var_rows, perf_rows):
    rng = {}
    cal = {}
    rp = os.path.join(indir, "rng_validation.json")
    cp = os.path.join(indir, "exact_calibration.json")
    if os.path.isfile(rp):
        with open(rp) as f:
            rng = json.load(f)
    if os.path.isfile(cp):
        with open(cp) as f:
            cal = json.load(f)
    lines = []
    a = lines.append
    a("# Issue #9 — Coupled Newman–Ziff CPU campaign")
    a("")
    a("Estimator: one Fisher–Yates permutation of `N=L^2` sites feeds **G** (square NN,")
    a("forward) and **G\\*** (NN+NNN matching, reverse / complement). Reverse permutation")
    a("is the matching complement, **not** an antithetic coupling. Independent mode uses")
    a("a second Philox stream (`stream_id=1`) with the same replica counts.")
    a("")
    a("RNG: Philox4x32-10 (Random123). Counter layout: `key=(global_seed, batch_id)`,")
    a("`ctr=(draw_counter, replica_id, stream_id, 0)`. Bounded integers: Lemire unbiased.")
    a("Production does not use `std::mt19937`.")
    a("")
    a("## RNG known-answer tests")
    a("")
    a(f"- source: `{rng.get('source', 'n/a')}`")
    a(f"- all_pass: **{rng.get('all_pass', 'n/a')}**")
    for kat in rng.get("kats", []):
        a(f"  - ctr={kat.get('ctr')} expected={kat.get('expected')} pass={kat.get('pass')}")
    a("")
    a("## Exact calibration (before production)")
    a("")
    a(f"- overall: **{cal.get('overall', 'n/a')}**")
    a(f"- exact source: `{cal.get('exact_source', 'n/a')}`")
    for section in ("subset_enumeration", "hand_crafted_and_prefix", "exhaustive_permutations", "mc_convergence"):
        block = cal.get(section, {})
        a(f"- {section}:")
        if isinstance(block, dict):
            for k, v in block.items():
                a(f"  - {k}: {v}")
    a("")
    a("## Root sequence")
    a("")
    a("Search window requested: `[0.590, 0.595]`. Printed digits are 18 significant;")
    a("uncertainty is the batch standard error, not print width.")
    a("")
    a("| L | observable | mode | pooled_root | batch_sd | batch_se | replicas | runtime_s |")
    a("|---|---|---|---|---|---|---|---|")
    for r in seq_rows:
        a(
            f"| {r['L']} | {r['observable']} | {r['mode']} | {r['pooled_root']} | "
            f"{r['batch_sd']} | {r['batch_se']} | {r['replicas']} | {r['runtime']} |"
        )
    a("")
    a("## Shared vs independent variance")
    a("")
    a("| L | observable | var_ratio (indep/shared roots) | runtime_ratio | var_ratio (X-Y at k~pc N) |")
    a("|---|---|---|---|---|")
    for r in var_rows:
        a(
            f"| {r['L']} | {r['observable']} | {r['variance_independent_over_shared']} | "
            f"{r['runtime_independent_over_shared']} | {r['var_diff_independent_over_shared']} |"
        )
    a("")
    a("## Performance")
    a("")
    a("| L | mode | replicas | threads | wall_s | site_updates/s | replicas/s |")
    a("|---|---|---|---|---|---|---|")
    for r in perf_rows:
        a(
            f"| {r.get('L')} | {r.get('mode')} | {r.get('replicas')} | {r.get('threads')} | "
            f"{r.get('wall_seconds')} | {r.get('site_updates_per_sec')} | {r.get('replicas_per_sec')} |"
        )
    a("")
    a("## Notes")
    a("")
    a("- Wrapping engine is Issue #7 DisplacementDSU in src/torus_connectivity.hpp;")
    a("  this campaign does not fork a second winding implementation.")
    a("- Microcanonical CSVs store aggregate sufficient statistics only (no permutations).")
    a("- Antithetic mode is not implemented.")
    a("")
    path = os.path.join(indir, "REPORT.md")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="results/issue-9")
    args = ap.parse_args()
    indir = args.indir
    mp.mp.dps = 50

    perf_rows = load_perf(os.path.join(indir, "performance.csv"))
    perf_map = {}
    for r in perf_rows:
        perf_map[(int(r["L"]), r["mode"])] = r

    batch_csv_rows = []
    seq_rows = []
    var_scratch = defaultdict(dict)

    bin_dir = os.path.join(indir, "batch_bin")
    if not os.path.isdir(bin_dir):
        print(f"no batch_bin in {indir}", file=sys.stderr)
        return 1

    jobs = []
    for fn in sorted(os.listdir(bin_dir)):
        if not fn.endswith(".bin"):
            continue
        body = fn[:-4]
        parts = body.split("_", 1)
        if len(parts) != 2:
            continue
        mode = parts[1]
        jobs.append((os.path.join(bin_dir, fn), mode))

    for path, mode in jobs:
        L, n, rpb, batches = read_batch_bin(path)
        n_batches = len(batches)
        replicas_total = rpb * n_batches
        runtime = 0.0
        pr = perf_map.get((L, mode))
        if pr:
            runtime = float(pr["wall_seconds"])

        pooled_cl_g = [0.0] * (n + 1)
        pooled_cl_gs = [0.0] * (n + 1)
        pooled_g = {name: [0.0] * (n + 1) for name in ("H", "V", "E", "B")}
        pooled_gs = {name: [0.0] * (n + 1) for name in ("H", "V", "E", "B")}
        for b in batches:
            for k in range(n + 1):
                pooled_cl_g[k] += b["cl_g"][k]
                pooled_cl_gs[k] += b["cl_gs"][k]
                for name in ("H", "V", "E", "B"):
                    pooled_g[name][k] += b["wrap_g"][name][k]
                    pooled_gs[name][k] += b["wrap_gs"][name][k]
        inv = 1.0 / n_batches
        for k in range(n + 1):
            pooled_cl_g[k] *= inv
            pooled_cl_gs[k] *= inv
            for name in ("H", "V", "E", "B"):
                pooled_g[name][k] *= inv
                pooled_gs[name][k] *= inv

        observables = [
            ("H", make_D_wrap(pooled_g["H"], pooled_gs["H"], n), "wrap"),
            ("V", make_D_wrap(pooled_g["V"], pooled_gs["V"], n), "wrap"),
            ("either", make_D_wrap(pooled_g["E"], pooled_gs["E"], n), "wrap"),
            ("both", make_D_wrap(pooled_g["B"], pooled_gs["B"], n), "wrap"),
            ("M", make_M(pooled_cl_g, pooled_cl_gs, n), "cluster"),
        ]
        batch_funcs_wrap = {
            "H": ("H", "H"),
            "V": ("V", "V"),
            "either": ("E", "E"),
            "both": ("B", "B"),
        }

        for obs_name, f_pooled, kind in observables:
            pooled_root, status, _ival = find_root(f_pooled)
            batch_roots = []
            for bi, b in enumerate(batches):
                if kind == "wrap":
                    gkey, skey = batch_funcs_wrap[obs_name]
                    fb = make_D_wrap(b["wrap_g"][gkey], b["wrap_gs"][skey], n)
                else:
                    fb = make_M(b["cl_g"], b["cl_gs"], n)
                rt, st, _ = find_root(fb)
                batch_roots.append(rt)
                batch_csv_rows.append(
                    {
                        "L": L,
                        "observable": obs_name,
                        "mode": mode,
                        "batch": bi,
                        "root": fmt18(rt),
                        "status": st,
                        "replicas": rpb,
                    }
                )
            _mu, sd, se = mean_sd_se(batch_roots)
            seq_rows.append(
                {
                    "L": L,
                    "observable": obs_name,
                    "mode": mode,
                    "pooled_root": fmt18(pooled_root),
                    "batch_sd": fmt18(sd),
                    "batch_se": fmt18(se),
                    "replicas": replicas_total,
                    "runtime": f"{runtime:.10g}",
                    "status": status,
                }
            )
            var_scratch[(L, obs_name)][mode] = {
                "root_var": sd * sd if sd == sd else float("nan"),
                "runtime": runtime,
                "pooled_root": pooled_root,
            }

        csv_name = f"microcanonical_L{L:03d}.csv"
        if mode == "independent":
            csv_name = f"microcanonical_L{L:03d}_independent.csv"
        csv_path = os.path.join(indir, csv_name)
        if os.path.isfile(csv_path):
            _n, _rep, _cg, _cgs, _wg, _wgs, stats = read_pooled_csv(csv_path)
            kstar = int(round(PC_REF * n))
            if kstar < 0:
                kstar = 0
            if kstar > n:
                kstar = n
            for obs_name, xk, yk, xy in (
                ("H", stats["h_g"], stats["h_gs"], stats["h_xy"]),
                ("V", stats["v_g"], stats["v_gs"], stats["v_xy"]),
                ("either", stats["e_g"], stats["e_gs"], stats["e_xy"]),
                ("both", stats["b_g"], stats["b_gs"], stats["b_xy"]),
                ("M", None, None, None),
            ):
                if obs_name == "M":
                    vd = var_diff(
                        _cg[kstar],
                        _cgs[kstar],
                        stats["cl_g2"][kstar],
                        stats["cl_gs2"][kstar],
                        stats["cl_xy"][kstar],
                    )
                else:
                    vd = var_diff(xk[kstar], yk[kstar], xk[kstar], yk[kstar], xy[kstar])
                var_scratch[(L, obs_name)].setdefault(mode, {})
                var_scratch[(L, obs_name)][mode]["var_diff_kstar"] = vd
                var_scratch[(L, obs_name)][mode]["kstar"] = kstar

        print(f"analyzed L={L} mode={mode} batches={n_batches} rpb={rpb}", file=sys.stderr)

    seq_path = os.path.join(indir, "root_sequence.csv")
    with open(seq_path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "L", "observable", "mode", "pooled_root", "batch_sd", "batch_se",
                "replicas", "runtime", "status",
            ],
        )
        w.writeheader()
        for r in sorted(seq_rows, key=lambda z: (int(z["L"]), z["mode"], z["observable"])):
            w.writerow(r)

    batch_path = os.path.join(indir, "roots_by_batch.csv")
    with open(batch_path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["L", "observable", "mode", "batch", "root", "status", "replicas"],
        )
        w.writeheader()
        for r in sorted(batch_csv_rows, key=lambda z: (int(z["L"]), z["mode"], z["observable"], int(z["batch"]))):
            w.writerow(r)

    var_rows = []
    for (L, obs), modes in sorted(var_scratch.items()):
        if "shared" in modes and "independent" in modes:
            vs = modes["shared"].get("root_var", float("nan"))
            vi = modes["independent"].get("root_var", float("nan"))
            ratio = vi / vs if (vs == vs and vs > 0.0 and vi == vi) else float("nan")
            rs = modes["shared"].get("runtime", float("nan"))
            ri = modes["independent"].get("runtime", float("nan"))
            rratio = ri / rs if (rs == rs and rs > 0.0 and ri == ri) else float("nan")
            vds = modes["shared"].get("var_diff_kstar", float("nan"))
            vdi = modes["independent"].get("var_diff_kstar", float("nan"))
            dratio = vdi / vds if (vds == vds and vds > 0.0 and vdi == vdi) else float("nan")
            var_rows.append(
                {
                    "L": L,
                    "observable": obs,
                    "variance_independent_over_shared": fmt18(ratio),
                    "runtime_independent_over_shared": fmt18(rratio),
                    "var_diff_independent_over_shared": fmt18(dratio),
                    "var_diff_shared": fmt18(vds),
                    "var_diff_independent": fmt18(vdi),
                    "root_var_shared": fmt18(vs),
                    "root_var_independent": fmt18(vi),
                    "kstar": modes["shared"].get("kstar", ""),
                }
            )
    var_path = os.path.join(indir, "variance_comparison.csv")
    with open(var_path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "L", "observable", "variance_independent_over_shared",
                "runtime_independent_over_shared", "var_diff_independent_over_shared",
                "var_diff_shared", "var_diff_independent", "root_var_shared",
                "root_var_independent", "kstar",
            ],
        )
        w.writeheader()
        for r in var_rows:
            w.writerow(r)

    report = write_report(indir, seq_rows, var_rows, perf_rows)
    print(f"wrote {seq_path}")
    print(f"wrote {batch_path}")
    print(f"wrote {var_path}")
    print(f"wrote {report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
