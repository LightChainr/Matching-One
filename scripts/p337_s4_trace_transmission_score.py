#!/usr/bin/env python3
"""One fixed Q4 central-trace response at the saved original pooled root."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import time

from p337_closed_source_score import Interval as I, interval_json, middle, sha
from p337_closed_source_finite_score import polynomial, derivative

ROOT = Path(__file__).resolve().parents[1]
FREEZE = "55fdba789a576d8d4c507372b7834f92cf506c80"
CONTRACT = ROOT / "analysis/p337_s4_trace_transmission_contract.json"
OLD = ROOT / "results/p337-closed-source-finite-coupling"
N, DELTA = 25, F(1152, 625)


def read_rows(path):
    with path.open(newline="") as handle:
        return [{k: int(v) for k, v in row.items()} for row in csv.DictReader(handle)]


def from_json(packet):
    return I(F(packet["lower_fraction"]), F(packet["upper_fraction"]))


def old_coefficients(rows):
    out = {key: [F(0)]*(N+1) for key in ("z", "q", "e")}
    for row in rows:
        k, q = row["k"], row["q"]
        weight = F(row["count"], 2**row["g"])
        for key, multiplier in (("z", 1), ("q", q), ("e", q*q)):
            out[key][k] += weight*multiplier
    return out


def trace_coefficients(rows, old_rows):
    out = [F(0)]*(N+1)
    projected_q = [F(0)]*(N+1)
    projected_e = [F(0)]*(N+1)
    marginal = defaultdict(int)
    for row in rows:
        k, g, q, count = (row[key] for key in ("k", "g", "q", "count"))
        bad2, n3 = row["bad2"], row["n_bad3"]
        if bad2 not in (0, 1) or not 0 <= n3 <= N or count <= 0:
            raise ValueError("invalid incoming seam-constraint count")
        beta = F(1, 6) + F(int(bad2 == 0), 2) - F(2, 3*4**n3)
        if q != 0 and beta:
            raise ArithmeticError("rank0/rank2 character annihilation failed")
        weight = F(count, 2**g)*beta
        out[k] += weight
        projected_q[k] += q*weight
        projected_e[k] += q*q*weight
        marginal[k, g, q] += count
    expected = {(r["k"], r["g"], r["q"]): r["count"] for r in old_rows}
    if dict(marginal) != expected:
        raise ValueError("new seam counts do not retain the locked occupation population")
    if any(projected_q) or any(projected_e):
        raise ArithmeticError("filtered q/E numerator is not identically zero")
    return out


def raw_jet(coeffs, h):
    first = derivative(coeffs)
    return [polynomial(coeffs, h), polynomial(first, h),
            polynomial(derivative(first), h)]


def geometry_packet(base, trace, h):
    z, zh, zhh = raw_jet(base["z"], h)
    if z.lo <= 0:
        raise ArithmeticError("original partition bound is not positive")
    result = {}
    for key in ("q", "e"):
        v, vh, vhh = raw_jet(base[key], h)
        mean = v/z
        first = (vh-mean*zh)/z
        second = (vhh-mean*zhh-2*first*zh)/z
        result[key] = [mean, first, second]
    raw_f, raw_fh, _ = raw_jet(trace, h)
    f = raw_f/z
    result["f"] = [f, (raw_fh-f*zh)/z]
    return result


def score(packets, saved):
    # D and R are the already published original values, not a newly found root
    # or a re-score of the old four-coupling experiment.
    D, R = from_json(saved["D_h"]), from_json(saved["U_over_A"])
    if D.lo <= 0:
        raise ArithmeticError("saved original thermal denominator is not positive")
    one, two = packets
    a = (one["q"][0]-two["q"][0])/2
    ah = (one["q"][1]-two["q"][1])/2
    e = (one["e"][0]+two["e"][0])/2
    eh = (one["e"][1]+two["e"][1])/2
    Y = (one["e"][0]-two["e"][0])/DELTA
    Yhh = (one["e"][2]-two["e"][2])/DELTA
    Mhh = (one["q"][2]+two["q"][2])/2
    H = Yhh-R*Mhh
    coefficients = {"common_thermal": -Y/D,
                    "geometric_thermal": (R*a-2*e/DELTA)/D,
                    "geometric_value": (R*ah-2*eh/DELTA+(a/D)*H)/D}
    fc_h = (one["f"][1]+two["f"][1])/2
    fd_h = (one["f"][1]-two["f"][1])/2
    fd = (one["f"][0]-two["f"][0])/2
    arguments = {"common_thermal": fc_h, "geometric_thermal": fd_h,
                 "geometric_value": fd}
    terms = {key: coefficients[key]*arguments[key] for key in coefficients}
    value = sum(terms.values())
    return {"value": value, "coefficients": coefficients, "arguments": arguments,
            "terms": terms, "root_h_derivative": a*fd/D}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    contract = json.loads(CONTRACT.read_text())
    if contract["N"] != N or contract["m"] != 2 or contract["geometries"] != [[5, 0], [4, 3]]:
        raise ValueError("this score implements only the fixed Q4 character question")
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    saved_path = OLD / "score/score.json"
    saved = next(row for row in json.loads(saved_path.read_text())["scores"] if row["m"] == 2)
    h = from_json(saved["critical_h"])
    packets, inputs = [], []
    for label in ("axis", "tilted"):
        old_path, new_path = OLD/f"{label}.csv", args.counts_dir/f"{label}.csv"
        old_rows, rows = read_rows(old_path), read_rows(new_path)
        base = old_coefficients(old_rows)
        trace = trace_coefficients(rows, old_rows)
        packets.append(geometry_packet(base, trace, h))
        inputs.append({"geometry": label,
                       "new_counts": str(new_path.resolve()), "new_sha256": sha(new_path),
                       "old_counts": str(old_path.relative_to(ROOT)), "old_sha256": sha(old_path),
                       "trace_coefficients_h": [str(x) for x in trace]})
    result = score(packets, saved)
    value = result["value"]
    decision = "nonzero_normalization_transmission" if value.lo > 0 or value.hi < 0 else "zero_not_excluded_stop_fixed_score"
    with localcontext() as ctx:
        ctx.prec = 60
        A = Decimal(N)**(Decimal(13)/Decimal(8))/2
        def scaled(interval):
            v = middle(interval)
            return float(A*Decimal(v.numerator)/Decimal(v.denominator))
        full_value = scaled(value)
        full_terms = {key: scaled(v) for key, v in result["terms"].items()}
    record = {"schema": "p337.s4-trace-transmission.score.v1",
              "decision": decision, "freeze_commit": FREEZE,
              "code_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "created_utc": datetime.now(timezone.utc).isoformat(), "contract": contract,
              "root": {"source": str(saved_path.relative_to(ROOT)), "sha256": sha(saved_path),
                       "m": 2, "critical_h": saved["critical_h"], "new_root_search": False},
              "inputs": inputs, "V_beta_over_A": interval_json(value), "V_beta_approx": full_value,
              "terms_over_A": {key: interval_json(v) for key, v in result["terms"].items()},
              "terms_approx": full_terms,
              "transmission_coefficients": {key: interval_json(v) for key, v in result["coefficients"].items()},
              "trace_arguments": {key: interval_json(v) for key, v in result["arguments"].items()},
              "root_h_derivative": interval_json(result["root_h_derivative"]),
              "geometry_trace_fraction": [{"f": interval_json(p["f"][0]), "f_h": interval_json(p["f"][1])} for p in packets],
              "direct_q_E_numerators": "identically_zero_as_exact_polynomials",
              "old_population_marginals": "retained_exactly",
              "source_coordinate": "epsilon coefficient of canonical Q4 [2,2] trace; not t or logQ",
              "new_random_samples": 0, "Q4_trace_responses_scored": 1,
              "new_coupling_points": 0,
              "elapsed_seconds": time.perf_counter()-started}
    (out/"score.json").write_text(json.dumps(record, indent=2)+"\n")
    (out/"REPORT.md").write_text(
        "# Fixed Q4 canonical trace transmission\n\n"
        f"Decision: **{decision}**.\n\n"
        f"Original pooled-root `V_beta={full_value:.16g}`; the exact reduced interval is in score.json.\n\n"
        "The direct q/E numerators vanish as exact polynomials. This is the response through the original normalizers, root and thermal denominator. "
        "The fixed source coordinate is the central-trace coefficient epsilon, not log Q. No Q1 field activation or continuum identity is inferred.\n\n"
        "The old m2 root and D/U bounds were imported unchanged; only missing seam-component constraints and this one score were new.\n")
    print(json.dumps({"decision": decision, "V_beta": full_value,
                      "V_beta_over_A": interval_json(value), "terms": full_terms,
                      "elapsed_seconds": record["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
