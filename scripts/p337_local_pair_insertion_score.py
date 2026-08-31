#!/usr/bin/env python3
"""One fixed local four-port insertion response; no enumeration or root search."""
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
import time

from p337_closed_source_score import Interval as I, interval_json, middle, sha
from p337_closed_source_finite_score import normalized_moments, from_json

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/p337_local_pair_insertion_contract.json"
BASELINE = ROOT / "results/p337-closed-source-n25/latest.json"
OLD = ROOT / "results/p337-closed-source-finite-coupling"
FREEZE = "d7f15e68"
N, DELTA = 25, F(1152, 625)


def commit(ref):
    return subprocess.check_output(["git", "rev-parse", ref+"^{commit}"], cwd=ROOT, text=True).strip()


def rows(path):
    with path.open(newline="") as handle:
        return [{key: int(value) for key, value in row.items()} for row in csv.DictReader(handle)]


def coefficients(path, old_path):
    old = {key: [F(0)]*(N+1) for key in ("z", "q", "e")}
    for row in rows(old_path):
        k, q, count = (row[key] for key in ("k", "q", "count"))
        for key, value in (("z", 1), ("q", q), ("e", q*q)):
            old[key][k] += count*value
    data = rows(path)
    if len(data) != N+1:
        raise ValueError("one complete K=0..25 profile is required")
    result = {**old, **{key: [F(0)]*(N+1) for key in ("s", "qs", "es")}}
    for k, row in enumerate(data):
        if row["k"] != k or row["count"] != math.comb(N, k):
            raise ValueError("the full occupation population is not retained")
        if [row["count"], row["sum_q"], row["sum_e"]] != [old[key][k] for key in ("z", "q", "e")]:
            raise ValueError("the original q/E population differs from the locked archive")
        if not -2*row["count"] <= row["sum_s2"] <= 0:
            raise ValueError("the predeclared local source sign or scale differs")
        for key, field in (("s", "sum_s2"), ("qs", "sum_qs2"), ("es", "sum_es2")):
            result[key][k] = F(row[field], 2)
    return result


def response(pair, h, D, ratio):
    packets = [normalized_moments(data, h) for data in pair]
    T = (packets[0]["q"][2]+packets[1]["q"][2])/2
    H = (packets[0]["e"][2]-packets[1]["e"][2])/DELTA
    sources = []
    for row in packets:
        q, qh, _ = row["q"]
        e, eh, _ = row["e"]
        s, sh, _ = row["s"]
        qs, qsh, _ = row["qs"]
        es, esh, _ = row["es"]
        sources.append({"j_q": qs-q*s, "j_q_h": qsh-qh*s-q*sh,
                        "j_e": es-e*s, "j_e_h": esh-eh*s-e*sh,
                        "mean_s": s, "mean_s_h": sh,
                        "raw_qs": qs, "raw_es": es})
    jq = (sources[0]["j_q"]+sources[1]["j_q"])/2
    jqh = (sources[0]["j_q_h"]+sources[1]["j_q_h"])/2
    jyh = (sources[0]["j_e_h"]-sources[1]["j_e_h"])/DELTA
    terms = {"direct_centered": jyh/D, "root_motion": -H*jq/(D**2),
             "slope_source": -ratio*jqh/D, "slope_root": ratio*T*jq/(D**2)}
    return {"V_over_A": sum(terms.values()), "terms": terms,
            "root_h_tangent": -jq/D, "geometry_source": sources,
            "jQ": jq, "jQ_h": jqh, "jY_h": jyh}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--kernel-gate-commit", required=True)
    parser.add_argument("--topology-gate-commit", required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    contract = json.loads(CONTRACT.read_text())
    if contract["N"] != N or contract["geometries"] != [[5, 0], [4, 3]]:
        raise ValueError("only the frozen N25 local pair decision is implemented")
    gates = {"kernel": commit(args.kernel_gate_commit), "topology": commit(args.topology_gate_commit)}
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    saved = json.loads(BASELINE.read_text())
    p = from_json(saved["root_enclosure"])
    h = p/(1-p)
    Dp = from_json(saved["matching_slope_enclosure"])
    D = Dp/((1+h)**2)
    ratio = from_json(saved["U25_over_A_enclosure"])
    if D.lo <= 0:
        raise ArithmeticError("the imported original denominator is not positive")
    pair, inputs = [], []
    for label in ("axis", "tilted"):
        path, old_path = args.counts_dir.resolve()/f"{label}.csv", OLD/f"{label}.csv"
        data = coefficients(path, old_path)
        pair.append(data)
        inputs.append({"geometry": label, "new_source_path": str(path), "new_source_sha256": sha(path),
                       "old_population_path": str(old_path.relative_to(ROOT)), "old_population_sha256": sha(old_path),
                       "source_coefficients_y": {key: [str(value) for value in data[key]] for key in ("s", "qs", "es")},
                       "old_q_E_marginals_retained": True, "configurations": 2**N})
    scored = response(pair, h, D, ratio)
    value = scored["V_over_A"]
    decision = "local_pair_original_U_zero_response_rejected" if value.lo > 0 or value.hi < 0 else "unresolved_stop_fixed_insertion"
    with localcontext() as context:
        context.prec = 60
        area = Decimal(N)**(Decimal(13)/Decimal(8))/2
        def full(interval):
            mid = middle(interval)
            return float(area*Decimal(mid.numerator)/Decimal(mid.denominator))
        numerical = {"V": full(value), "terms": {key: full(v) for key, v in scored["terms"].items()}}
    result = {"schema": "p337.local-pair-insertion.score.v1", "status": "completed_one_fixed_response",
              "decision": decision, "contract": contract, "freeze_commit": commit(FREEZE),
              "code_commit": commit("HEAD"), "accepted_theory_gates": gates,
              "created_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs,
              "V_over_A": interval_json(value), "numerical_values": numerical,
              "terms_over_A": {key: interval_json(v) for key, v in scored["terms"].items()},
              "root_h_tangent": interval_json(scored["root_h_tangent"]),
              "root_p_tangent": interval_json(scored["root_h_tangent"]/((1+h)**2)),
              "geometry_source": [{key: interval_json(v) for key, v in packet.items()} for packet in scored["geometry_source"]],
              "source_covariances": {key: interval_json(scored[key]) for key in ("jQ", "jQ_h", "jY_h")},
              "baseline": {"path": str(BASELINE.relative_to(ROOT)), "sha256": sha(BASELINE),
                           "imported_root_p": saved["root_enclosure"], "converted_root_h": interval_json(h),
                           "imported_Dp": saved["matching_slope_enclosure"], "converted_Dh": interval_json(D),
                           "imported_U_over_A": saved["U25_over_A_enclosure"], "new_root_search": False},
              "normalization": "S=-site-average(t); fixed-origin first-source global moments are exactly equal by translation symmetry",
              "interval_scope": "exact rational outward bounds for finite count polynomials, not statistical confidence intervals",
              "new_random_samples": 0, "old_sources_rescored": False, "response_scores": 1,
              "boundary": "local C4-averaged End(pair) tensor; not a pure fixed-cut central projector, a free-index primary, or the full seam trace",
              "elapsed_seconds": time.perf_counter()-started}
    (out/"score.json").write_text(json.dumps(result, indent=2)+"\n")
    (out/"REPORT.md").write_text(
        "# Fixed local four-port pair insertion into original U\n\n"
        f"Decision: **{decision}**. Original-U response **{numerical['V']:.16g}**.\n\n"
        "The prescribed local tensor is C4 averaged and site averaged. Its source is -t at Q1; fixed-origin counting supplies all first-source global moments by exact translation symmetry. "
        "The original Q1 root, D and U/A were imported. Direct q/E source moments, covariance centering, thermal derivatives, pooled-root movement and slope response are all retained.\n\n"
        "This is one finite microscopic tensor perturbation, not a continuum field assignment or an equality to the full seam trace. The full rational bounds, four terms and source moments are in score.json. "
        "No root search, source fitting, Q/seam/support-radius scan or old-source rescore was run.\n")
    print(json.dumps({"decision": decision, **numerical, "elapsed_seconds": result["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
