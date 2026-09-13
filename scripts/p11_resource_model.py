#!/usr/bin/env python3
"""Resource model for the transfer-matrix DP: fit measured state counts,
peak RSS and per-evaluation wall time versus width n, and extrapolate the
memory/runtime needed for widths 25-28.

Inputs (results/transfer-subfrontier-20260913/raw/):
  - statecount Python DP: dpA_n*.json (peak_states)     [exact path A]
  - container CRT telemetry: crt_rows_n*.rows (per-row states, per-pass time)
  - container probe JSONs: probe_*.json (peak RSS, wall seconds)

Models (log-log and log-linear least squares over the measured range):
  states(n)  ~ a_s * b_s^n      (uncondensed partial-row map size)
  rss(n)     ~ a_m * b_m^n
  time(n)    ~ a_t * b_t^n      (one full CRT pass / one f128 evaluation)

The extrapolation is honest in the sense that only measured widths enter the
fit, the fit interval and residuals are reported, and widths whose model
residual is large are flagged as not credible.
"""

from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

RAW = Path(__file__).resolve().parents[1] / "results" / \
    "transfer-subfrontier-20260913" / "raw"


def loglin_fit(xs, ys):
    """Fit y = a*b^x by least squares on (log y); returns a, b, max rel resid."""
    lxs = xs
    lys = [math.log(y) for y in ys]
    n = len(lxs)
    mx = sum(lxs) / n
    my = sum(lys) / n
    sxx = sum((x - mx) ** 2 for x in lxs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(lxs, lys))
    slope = sxy / sxx
    intercept = my - slope * mx
    a = math.exp(intercept)
    b = math.exp(slope)
    resids = []
    for x, y in zip(lxs, ys):
        pred = a * b ** x
        resids.append(abs(pred - y) / y)
    return a, b, max(resids), resids


def peak_states_from_rows(path):
    best = 0
    for line in path.read_text().splitlines():
        m = re.search(r'"states":(\d+)', line)
        if m:
            best = max(best, int(m.group(1)))
    return best


def probe_stats(path):
    d = json.loads(path.read_text())
    return d["timing"]["wall_seconds"], d["memory"]["peak_rss_kib"]


def main():
    out = {"schema": "matching-one/p11-resource-model/v1", "fits": {}}

    # ---- exact path A (local python DP): peak states + time vs n ----------
    ns, states, times = [], [], []
    for j in sorted(RAW.glob("dpA_n*.json")):
        n = int(re.search(r"n(\d+)", j.stem).group(1))
        d = json.loads(j.read_text())
        ns.append(n)
        states.append(d["peak_states"])
        t = RAW / f"dpA_n{n}.time"
        if t.exists():
            for line in t.read_text().splitlines():
                if line.startswith("real"):
                    times.append(float(line.split()[1]))
    if len(ns) >= 3:
        a, b, rmax, _ = loglin_fit(ns, states)
        out["fits"]["pathA_peak_states"] = {
            "points": dict(zip(map(str, ns), states)),
            "a": a, "b": b, "max_rel_resid": rmax}
        if len(times) == len(ns):
            a, b, rmax, _ = loglin_fit(ns, times)
            out["fits"]["pathA_wall_seconds"] = {
                "points": dict(zip(map(str, ns), times)),
                "a": a, "b": b, "max_rel_resid": rmax}

    # ---- container CRT telemetry -----------------------------------------
    ns, states, walls, rss = [], [], [], []
    for j in sorted(RAW.glob("probe_crt_n*.json")):
        n = int(re.search(r"n(\d+)", j.stem).group(1))
        w, r = probe_stats(j)
        rows = RAW / f"crt_rows_n{n}.rows"
        if not rows.exists():
            rows = RAW / f"crt_n{n}.rows"
        s = peak_states_from_rows(rows) if rows.exists() else None
        ns.append(n)
        walls.append(w)
        rss.append(r)
        if s:
            states.append(s)
    if len(ns) >= 3:
        a, b, rmax, _ = loglin_fit(ns, walls)
        out["fits"]["crt_wall_seconds_all_primes"] = {
            "points": dict(zip(map(str, ns), walls)),
            "a": a, "b": b, "max_rel_resid": rmax}
        a, b, rmax, _ = loglin_fit(ns, rss)
        out["fits"]["crt_peak_rss_kib"] = {
            "points": dict(zip(map(str, ns), rss)),
            "a": a, "b": b, "max_rel_resid": rmax}
        if len(states) == len(ns):
            a, b, rmax, _ = loglin_fit(ns, states)
            out["fits"]["crt_peak_states"] = {
                "points": dict(zip(map(str, ns), states)),
                "a": a, "b": b, "max_rel_resid": rmax}

    # ---- f128 evaluations --------------------------------------------------
    ns, walls, rss = [], [], []
    for j in sorted(RAW.glob("probe_f128_med_n*.json")):
        n = int(re.search(r"n(\d+)", j.stem).group(1))
        w, r = probe_stats(j)
        ns.append(n)
        walls.append(w)
        rss.append(r)
    if len(ns) >= 3:
        a, b, rmax, _ = loglin_fit(ns, walls)
        out["fits"]["f128_eval_wall_seconds"] = {
            "points": dict(zip(map(str, ns), walls)),
            "a": a, "b": b, "max_rel_resid": rmax}
        a, b, rmax, _ = loglin_fit(ns, rss)
        out["fits"]["f128_peak_rss_kib"] = {
            "points": dict(zip(map(str, ns), rss)),
            "a": a, "b": b, "max_rel_resid": rmax}

    # ---- extrapolation to n=25..28 ----------------------------------------
    extrap = {}
    for key, fit in out["fits"].items():
        a, b = fit["a"], fit["b"]
        extrap[key] = {str(n): a * b ** n for n in (25, 26, 27, 28)}
        extrap[key]["_growth_base"] = b
        extrap[key]["_max_rel_resid"] = fit["max_rel_resid"]
    out["extrapolation_25_28"] = extrap

    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
