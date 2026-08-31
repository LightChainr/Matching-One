#!/usr/bin/env python3
"""One fixed four-coupling score from exact (K,g,q) counts; no enumeration."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
import math
from pathlib import Path
import subprocess
import sys
import time

from p337_closed_source_score import Interval as I, interval_json, middle, sha

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/p337_closed_source_finite_coupling_contract.json"
BASELINE = ROOT / "results/p337-closed-source-n25/latest.json"
FREEZE = "b70dc4bd2fddd7676e9536b42bf912ee00ad302f"
N, DELTA = 25, F(1152, 625)
FIELDS = ("z", "q", "e", "s", "qs", "es")


def read_counts(path):
    with path.open(newline="") as handle:
        rows = [{key: int(value) for key, value in row.items()}
                for row in csv.DictReader(handle)]
    totals = [0] * (N+1)
    for row in rows:
        k, g, q, count = (row[key] for key in ("k", "g", "q", "count"))
        if not (0 <= k <= N and g >= 0 and q in (-1, 0, 1) and count > 0):
            raise ValueError("invalid exact histogram row")
        totals[k] += count
    if totals != [math.comb(N, k) for k in range(N+1)]:
        raise ValueError("incomplete exact configuration population")
    return rows


def coefficients(rows, multiplier):
    result = {key: [F(0)] * (N+1) for key in FIELDS}
    for row in rows:
        k, g, q, count = (row[key] for key in ("k", "g", "q", "count"))
        weight = F(count, multiplier**g)
        for key, value in zip(FIELDS, (1, q, q*q, -g, -q*g, -q*q*g)):
            result[key][k] += weight * value
    return result


def polynomial(coeffs, h):
    value = 0
    for coefficient in reversed(coeffs):
        value = value*h + coefficient
    return value


def convolve(a, b):
    result = [F(0)] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def root_interval(pair, steps):
    # Z_a,Z_b are positive: this numerator has the sign of pooled Q.
    one = convolve(pair[0]["q"], pair[1]["z"])
    two = convolve(pair[1]["q"], pair[0]["z"])
    numerator = [a+b for a, b in zip(one, two)]
    lo, hi = F(0), F(1)
    if polynomial(numerator, lo) >= 0:
        raise ValueError("empty configuration is not rank zero")
    while polynomial(numerator, hi) < 0:
        hi *= 2
    if polynomial(numerator, hi) == 0:
        return I.of(hi)
    for _ in range(steps):
        mid = (lo+hi)/2
        value = polynomial(numerator, mid)
        if value == 0:
            return I.of(mid)
        if value < 0:
            lo = mid
        else:
            hi = mid
    return I(lo, hi)


def derivative(coeffs):
    return [k*c for k, c in enumerate(coeffs)][1:]


def normalized_moments(data, h):
    raw = {}
    for key, coeffs in data.items():
        first = derivative(coeffs)
        raw[key] = [polynomial(coeffs, h), polynomial(first, h),
                    polynomial(derivative(first), h)]
    z, zp, zpp = raw["z"]
    normalized = {}
    for key in FIELDS[1:]:
        value, first, second = raw[key]
        mean = value/z
        mean_p = (first-mean*zp)/z
        mean_pp = (second-mean*zpp-2*mean_p*zp)/z
        normalized[key] = [mean, mean_p, mean_pp]
    return normalized


def score(pair, h):
    packet = [normalized_moments(data, h) for data in pair]
    D = (packet[0]["q"][1]+packet[1]["q"][1])/2
    B = (packet[0]["e"][1]-packet[1]["e"][1])/DELTA
    T = (packet[0]["q"][2]+packet[1]["q"][2])/2
    H = (packet[0]["e"][2]-packet[1]["e"][2])/DELTA
    if D.lo <= 0:
        raise ArithmeticError("matching slope enclosure is not strictly positive")
    jq, jqp, jep = [], [], []
    for row in packet:
        q, qp, _ = row["q"]
        e, ep, _ = row["e"]
        s, sp, _ = row["s"]
        qs, qsp, _ = row["qs"]
        _, esp, _ = row["es"]
        jq.append(qs-q*s)
        jqp.append(qsp-qp*s-q*sp)
        jep.append(esp-ep*s-e*sp)
    jQ, jQp = sum(jq)/2, sum(jqp)/2
    jYp = (jep[0]-jep[1])/DELTA
    terms = {
        "direct": jYp/D,
        "root_motion": -H*jQ/(D**2),
        "slope_source": -B*jQp/(D**2),
        "slope_root": B*T*jQ/(D**3),
    }
    return {"U_over_A": B/D, "V_over_A": sum(terms.values()), "terms": terms,
            "D_h": D, "critical_h_tangent": -jQ/D}


def from_json(packet):
    return I(F(packet["lower_fraction"]), F(packet["upper_fraction"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    start_utc = datetime.now(timezone.utc).isoformat()
    contract = json.loads(CONTRACT.read_text())
    if contract["multipliers_m"] != [2, 4, 8, 16] or contract["N"] != N:
        raise ValueError("this scorer only implements the fixed decision contract")
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    paths = [args.counts_dir.resolve()/name for name in ("axis.csv", "tilted.csv")]
    rows = [read_counts(path) for path in paths]
    baseline = json.loads(BASELINE.read_text())
    initial_u = from_json(baseline["U25_over_A_enclosure"])
    initial_v = from_json(baseline["source_enclosures"]["sstar"]["V_over_A"])
    if initial_u.lo <= 0 or initial_v.lo <= 0:
        raise ValueError("the recorded positive baseline premise is not resolved")
    records, downturn_at = [], []
    with localcontext() as context:
        context.prec = 70
        area = Decimal(N)**(Decimal(13)/Decimal(8))/2

        def full_value(interval):
            value = middle(interval)
            return float(area*Decimal(value.numerator)/Decimal(value.denominator))

        for m in contract["multipliers_m"]:
            pair = [coefficients(data, m) for data in rows]
            h = root_interval(pair, contract["root_bisection_steps"])
            scored = score(pair, h)
            p = m*h/(1+m*h)
            logit_velocity = 1+scored["critical_h_tangent"]/h
            down = scored["V_over_A"].hi < 0
            if down:
                downturn_at.append(m)
            record = {
                "m": m, "t_approx": math.log(m),
                "critical_h": interval_json(h), "critical_p": interval_json(p),
                "U_over_A": interval_json(scored["U_over_A"]),
                "dU_dt_over_A": interval_json(scored["V_over_A"]),
                "D_h": interval_json(scored["D_h"]),
                "critical_h_tangent": interval_json(scored["critical_h_tangent"]),
                "critical_logit_velocity": interval_json(logit_velocity),
                "terms_over_A": {key: interval_json(value) for key, value in scored["terms"].items()},
                "numerical_values": {"h": float(middle(h)), "p": float(middle(p)),
                                     "U": full_value(scored["U_over_A"]),
                                     "dU_dt": full_value(scored["V_over_A"])},
                "downturn_resolved": down,
            }
            records.append(record)
            print(json.dumps({"m": m, **record["numerical_values"], "downturn_resolved": down}), flush=True)
    decision = ("positive_source_global_U_turnover_resolved" if downturn_at else
                "fixed_grid_does_not_resolve_downturn_stop_without_extension")
    result = {
        "schema": "matching-one.p337-closed-source-finite-coupling.score.v1",
        "status": "completed_fixed_four_coupling_score", "decision": decision,
        "contract_freeze_commit": FREEZE, "contract": contract,
        "code_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "input_counts": [{"path": str(path.relative_to(ROOT)), "sha256": sha(path)} for path in paths],
        "baseline": {"path": str(BASELINE.relative_to(ROOT)), "sha256": sha(BASELINE),
                     "used": "published positive U(0),U_t(0) signs only; no baseline rescoring",
                     "numerical_values": {key: baseline["numerical_values"][key] for key in ("U25", "V25_Sstar")}},
        "downturn_resolved_at_m": downturn_at, "scores": records,
        "interval_scope": "rational outward enclosures of finite count polynomials, not statistical confidence intervals; positive area factor displayed numerically",
        "source_chart": "-g at fixed h produces the same full U(t) as Sstar at fixed p, by a common homogeneous thermal reparametrization",
        "dependency_group": "p337-N25-axis-tilted-exhaustive-configuration-populations; same finite graphs as baseline, no independent random evidence",
        "new_random_samples": 0, "cloud_jobs": 0, "tests_run": 0,
    }
    (out/"score.json").write_text(json.dumps(result, indent=2)+"\n")
    table = "\n".join(f"| {r['m']} | {r['numerical_values']['p']:.12g} | {r['numerical_values']['U']:.12g} | {r['numerical_values']['dU_dt']:.12g} |" for r in records)
    report = f"""# Fixed finite-coupling closed-source turnover

Decision: `{decision}`. Negative exact U_t/A enclosures occur at
`m={downturn_at}`. The already published U(0)>0 and U_t(0)>0, together
with the unique analytic positive-coupling root, imply a positive-coupling
maximum before any resolved negative-derivative point.

| m=exp(t) | critical p | U | dU/dt |
|---:|---:|---:|---:|
{table}

These are four preselected finite-source laws, not a fitted scan.
The freeze is `{FREEZE}`. Every prescribed value is reported; no grid
extension or new source is allowed by this decision. The baseline values
are imported from the prior exact packet without rescoring.

## Exact finite law and root-complete score

The integer histogram records (K,g,q), where
`g=2N+1-K-Sstar=2K-(beta1+beta_null)`.
With m=exp(t) and h=p/((1-p)m), its normalized weights are h^K/m^g.
Each geometry is normalized separately before the pooled matching mean
Q and P4(E). The pooled-root numerator has degree at most50. The root
uses the frozen 128 rational bisections; all reported enclosures are
outward rational arithmetic bounds, not sampling confidence intervals.

Let D=Q_h, B=Y_h, T=Q_hh, H=Y_hh, and jO=Cov(O,-g) within each geometry.
The four-term response is
`U_t/A=jY_h/D-H*jQ/D^2-B*jQ_h/D^2+B*T*jQ/D^3`.
The moving root is h_t=-jQ/D. Its common p-to-h Jacobian cancels in U;
using -g in this chart exactly differentiates the same homogeneous
Sstar curve. It is not a source substitution on an inhomogeneous
checkerboard saturation chart.

## Scope and next consequence

The result concerns the fixed N25 axis/tilted pair and the named source.
It resolves the monotone-amplification alternative if a downturn is
present; it does not identify a thermodynamic transition or continuum
field. The new full histogram enumerates the same finite configuration
populations as the old first-moment packet; it is not an independent
statistical vote. There were no new random samples, cloud jobs or test
campaign. The finite-volume strong-coupling argument predicts U->0;
the table is the frozen finite-coupling consequence, not a fit of its
asymptotic rate.
"""
    (out/"REPORT.md").write_text(report)
    receipt = {"started_utc": start_utc, "finished_utc": datetime.now(timezone.utc).isoformat(),
               "elapsed_seconds": time.perf_counter()-started, "command": sys.argv,
               "python": sys.version, "hashes": {"contract": sha(CONTRACT), "scorer": sha(__file__),
               "interval_backend": sha(ROOT/"scripts/p337_closed_source_score.py"), "score": sha(out/"score.json")},
               "new_random_samples": 0, "cloud_jobs": 0, "tests_run": 0}
    (out/"run.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"decision": decision, "elapsed_seconds": receipt["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
