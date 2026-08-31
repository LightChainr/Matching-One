#!/usr/bin/env python3
"""One frozen canonical joint-Q original-U readout, with a fixed contact split."""
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
CONTRACT = ROOT / "analysis/p337_regular_pair_joint_contract.json"
BASELINE = ROOT / "results/p337-closed-source-n25/latest.json"
OLD = ROOT / "results/p337-closed-source-finite-coupling"
FREEZE = "4ce4dfe894c9fe96f268c61cf21eb6585dba5418"
N = 25
PARTS = {"total": "", "adjacent": "_adj", "nonadjacent": "_far"}


def coefficients(path, old_path):
    base = {key: [F(0)]*(N+1) for key in ("z", "q", "e")}
    for row in rows(old_path):
        k, q, count = (row[key] for key in ("k", "q", "count"))
        for key, value in (("z", 1), ("q", q), ("e", q*q)):
            base[key][k] += count*value
    data = rows(path)
    if len(data) != N+1:
        raise ValueError("one complete K=0..25 profile is required")
    parts = {name: {**base, **{key: [F(0)]*(N+1) for key in ("s", "qs", "es")}}
             for name in PARTS}
    for k, row in enumerate(data):
        if row["k"] != k or row["count"] != math.comb(N, k):
            raise ValueError("the full occupation population must be retained")
        if [row["count"], row["sum_q"], row["sum_e"]] != [base[key][k] for key in ("z", "q", "e")]:
            raise ValueError("the original q/E population differs from the locked archive")
        for key, field in (("s", "sum_b16"), ("qs", "sum_qb16"), ("es", "sum_eb16")):
            if row[field] != row[field+"_adj"] + row[field+"_far"]:
                raise ValueError("the frozen contact split does not sum to the total")
            for name, suffix in PARTS.items():
                parts[name][key][k] = F(row[field+suffix], 16*N)
        if abs(row["sum_b16"]) > 43*(N-1)*row["count"]:
            raise ValueError("joint kernel integer units differ from the frozen table")
    return parts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--response-proof-commit", required=True)
    parser.add_argument("--adjacent-proof-commit", required=True)
    parser.add_argument("--producer-commit", required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    contract = json.loads(CONTRACT.read_text())
    if contract["N"] != N or contract["geometries"] != [[5, 0], [4, 3]]:
        raise ValueError("only the fixed canonical N25 joint response is implemented")
    gates = {key: commit(getattr(args, key+"_commit"))
             for key in ("response_proof", "adjacent_proof", "producer")}
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    saved = json.loads(BASELINE.read_text())
    p = from_json(saved["root_enclosure"])
    h = p/(1-p)
    D = from_json(saved["matching_slope_enclosure"])/((1+h)**2)
    ratio = from_json(saved["U25_over_A_enclosure"])
    if D.lo <= 0:
        raise ArithmeticError("the imported original denominator is not positive")

    data, inputs = [], []
    for label in ("axis", "tilted"):
        path = args.counts_dir.resolve()/f"{label}.csv"
        old_path = OLD/f"{label}.csv"
        parts = coefficients(path, old_path)
        data.append(parts)
        inputs.append({"geometry": label, "source_file": path.name, "source_sha256": sha(path),
                       "old_population": str(old_path.relative_to(ROOT)), "old_population_sha256": sha(old_path),
                       "configurations": 2**N,
                       "source_coefficients_y": {name: {key: [str(v) for v in parts[name][key]]
                                                       for key in ("s", "qs", "es")} for name in PARTS}})

    scored = {name: response([row[name] for row in data], h, D, ratio) for name in PARTS}
    packets, numbers = {}, {}
    with localcontext() as context:
        context.prec = 60
        area = Decimal(N)**(Decimal(13)/Decimal(8))/2

        def full(interval):
            value = middle(interval)
            return float(area*Decimal(value.numerator)/Decimal(value.denominator))

        for name, packet in scored.items():
            value = packet["V_over_A"]
            packets[name] = {
                "J2_over_A": interval_json(value),
                "zero_excluded": value.lo > 0 or value.hi < 0,
                "terms_over_A": {key: interval_json(v) for key, v in packet["terms"].items()},
                "root_h_mixed_tangent": interval_json(packet["root_h_tangent"]),
                "root_p_mixed_tangent": interval_json(packet["root_h_tangent"]/((1+h)**2)),
                "geometry_source": [{key: interval_json(v) for key, v in row.items()}
                                    for row in packet["geometry_source"]],
                "source_covariances": {key: interval_json(packet[key]) for key in ("jQ", "jQ_h", "jY_h")}}
            numbers[name] = {"J2": full(value), "terms": {key: full(v) for key, v in packet["terms"].items()}}

    decisions = {
        "primary": "additive_linear_first_Q_global_closure_rejected" if packets["total"]["zero_excluded"]
                   else "total_joint_response_unresolved_stop_fixed_model",
        "secondary": "adjacent_contact_only_global_closure_rejected" if packets["nonadjacent"]["zero_excluded"]
                     else "nonadjacent_global_response_unresolved_stop_fixed_split"}
    result = {
        "schema": "p337.regular-pair-joint.score.v1", "status": "completed_fixed_joint_response",
        "decisions": decisions, "contract": contract, "freeze_commit": commit(FREEZE),
        "code_commit": commit("HEAD"), "accepted_gates": gates,
        "created_utc": datetime.now(timezone.utc).isoformat(), "inputs": inputs,
        "parts": packets, "numerical_values": numbers,
        "baseline": {"path": str(BASELINE.relative_to(ROOT)), "sha256": sha(BASELINE),
                     "root_p": saved["root_enclosure"], "D_p": saved["matching_slope_enclosure"],
                     "U_over_A": saved["U25_over_A_enclosure"], "new_root_search": False},
        "normalization": "s2=ordered-pair-sum(g)/N^2; translation-reduced g16 sum/(16*N), with no further factor two",
        "direct_Q1_epsilon_derivatives": "all zero on the finite regular root branch",
        "interval_scope": "exact finite rational enclosures, not statistical confidence intervals",
        "dependency": "same full exact N25 populations; adjacent and nonadjacent are additive components, not independent evidence",
        "new_random_samples": 0, "old_sources_rescored": False, "fixed_readouts": 1,
        "field_assignment": None, "counterterm_fitted": False,
        "elapsed_seconds": time.perf_counter()-started}
    (out/"score.json").write_text(json.dumps(result, indent=2)+"\n")
    report = ["# Canonical joint first-Q transmission into original U", "",
              f"Primary: **{decisions['primary']}**.", f"Secondary: **{decisions['secondary']}**.", "",
              "| Fixed pair class | J2 |", "|---|---:|"]
    report += [f"| {name} | {values['J2']:.16g} |" for name, values in numbers.items()]
    report += ["", "J2=d_logQ d_epsilon^2 U for the canonical site-average regular interaction. "
               "The exact original root, slope and U/A are imported; all centering, thermal, root and slope terms remain. "
               "The nonadjacent part includes every distinct non-NN pair on N25, not a continuum long-distance limit. "
               "No counterterm, geometry or contact class was selected after readout. Full rational intervals and provenance are in score.json."]
    (out/"REPORT.md").write_text("\n".join(report)+"\n")
    print(json.dumps({"decisions": decisions, "numerical_values": numbers,
                      "elapsed_seconds": result["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
