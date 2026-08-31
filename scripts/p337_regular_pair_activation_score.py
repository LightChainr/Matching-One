#!/usr/bin/env python3
"""One canonical regular local tensor's mixed logQ/epsilon response."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
import math
from pathlib import Path
import time

from p337_closed_source_score import interval_json, middle, sha
from p337_closed_source_finite_score import from_json
from p337_local_pair_insertion_score import commit, response, rows

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/p337_regular_pair_activation_contract.json"
BASELINE = ROOT / "results/p337-closed-source-n25/latest.json"
OLD = ROOT / "results/p337-closed-source-finite-coupling"
FREEZE = "25f70f68"
N = 25


def coefficients(path, old_path):
    result = {key: [F(0)]*(N+1) for key in ("z", "q", "e", "s", "qs", "es")}
    for row in rows(old_path):
        k, q, count = (row[key] for key in ("k", "q", "count"))
        for key, value in (("z", 1), ("q", q), ("e", q*q)):
            result[key][k] += count*value
    data = rows(path)
    if len(data) != N+1:
        raise ValueError("a complete K=0..25 profile is required")
    for k, row in enumerate(data):
        if row["k"] != k or row["count"] != math.comb(N, k):
            raise ValueError("the full occupation population must be retained")
        if [row["count"], row["sum_q"], row["sum_e"]] != [result[key][k] for key in ("z", "q", "e")]:
            raise ValueError("the original q/E population differs from the locked archive")
        if not -2*row["count"] <= row["sum_a4"] <= 4*row["count"]:
            raise ValueError("the prescribed integer activation range differs")
        for key, field in (("s", "sum_a4"), ("qs", "sum_qa4"), ("es", "sum_ea4")):
            result[key][k] = F(row[field], 4)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--algebra-commit", required=True)
    parser.add_argument("--geometry-commit", required=True)
    parser.add_argument("--producer-commit", required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    contract = json.loads(CONTRACT.read_text())
    if contract["N"] != N or contract["geometries"] != [[5, 0], [4, 3]]:
        raise ValueError("only the frozen canonical N25 activation is implemented")
    gates = {key: commit(getattr(args, key+"_commit")) for key in ("algebra", "geometry", "producer")}
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
        inputs.append({"geometry": label, "new_source_file": path.name, "new_source_sha256": sha(path),
                       "old_population_path": str(old_path.relative_to(ROOT)), "old_population_sha256": sha(old_path),
                       "activation_coefficients_y": {key: [str(value) for value in data[key]] for key in ("s", "qs", "es")},
                       "original_q_E_marginals_retained": True, "configurations": 2**N})
    scored = response(pair, h, D, ratio)
    value = scored["V_over_A"]
    decision = "canonical_mixed_response_zero_rejected" if value.lo > 0 or value.hi < 0 else "unresolved_stop_fixed_completion"
    with localcontext() as context:
        context.prec = 60
        area = Decimal(N)**(Decimal(13)/Decimal(8))/2
        def full(interval):
            mid = middle(interval)
            return float(area*Decimal(mid.numerator)/Decimal(mid.denominator))
        numerical = {"W": full(value), "terms": {key: full(v) for key, v in scored["terms"].items()}}
    result = {"schema": "p337.regular-pair-activation.score.v1", "status": "completed_one_fixed_mixed_response",
              "decision": decision, "contract": contract, "freeze_commit": commit(FREEZE),
              "code_commit": commit("HEAD"), "accepted_gates": gates,
              "created_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs,
              "W_over_A": interval_json(value), "numerical_values": numerical,
              "direct_epsilon_response_at_Q1": "0 exactly for every finite graph and thermal activity",
              "terms_over_A": {key: interval_json(v) for key, v in scored["terms"].items()},
              "root_h_mixed_tangent": interval_json(scored["root_h_tangent"]),
              "root_p_mixed_tangent": interval_json(scored["root_h_tangent"]/((1+h)**2)),
              "geometry_source": [{key: interval_json(v) for key, v in packet.items()} for packet in scored["geometry_source"]],
              "source_covariances": {key: interval_json(scored[key]) for key in ("jQ", "jQ_h", "jY_h")},
              "baseline": {"path": str(BASELINE.relative_to(ROOT)), "sha256": sha(BASELINE),
                           "root_p": saved["root_enclosure"], "D_p": saved["matching_slope_enclosure"],
                           "U_over_A": saved["U25_over_A_enclosure"], "new_root_search": False},
              "normalization": "site-average a; coefficientwise beta_reg(1)=0 makes this the full mixed logQ/epsilon derivative",
              "interval_scope": "exact finite rational bounds, not statistical confidence intervals",
              "new_random_samples": 0, "old_sources_rescored": False, "response_scores": 1,
              "boundary": "canonical c(Q)=1 local completion only; no continuum or completion-independence claim",
              "elapsed_seconds": time.perf_counter()-started}
    (out/"score.json").write_text(json.dumps(result, indent=2)+"\n")
    (out/"REPORT.md").write_text(
        "# Canonical regular local pair Q activation\n\n"
        f"Decision: **{decision}**. Full mixed original-U response **W={numerical['W']:.16g}**.\n\n"
        "Here W=d_logQ d_epsilon U at Q1,epsilon0. The direct epsilon response is exactly zero. "
        "The completion is the fixed C4-averaged i(I-P1)i^dagger, with site-average normalization. "
        "All source covariances, thermal derivatives, induced root motion and slope terms are retained.\n\n"
        "The root, denominator and U/A were imported; only the fixed activation crossmoments are new. "
        "Same exact N25 populations, no independent stochastic evidence, no root search or fitted completion. "
        "Full rational bounds, coefficients and source packets are in score.json.\n")
    print(json.dumps({"decision": decision, **numerical, "elapsed_seconds": result["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
