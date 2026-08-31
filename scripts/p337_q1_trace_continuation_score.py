#!/usr/bin/env python3
"""Fixed Q1 trace packets; run only after the two declared theory gates pass."""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import time

from p337_closed_source_score import Interval as I, interval_json, middle, sha
from p337_s4_trace_transmission_score import (
    read_rows, from_json, geometry_packet, score as normalization_map,
)

ROOT = Path(__file__).resolve().parents[1]
FREEZE = "964ef2032effbe59f9158c158cf06a2c0844d7ee"
CONTRACT = ROOT / "analysis/p337_q1_trace_continuation_contract.json"
BASELINE = ROOT / "results/p337-closed-source-n25/latest.json"
OLD = ROOT / "results/p337-closed-source-finite-coupling"
N = 25
PACKETS = ("primary_beta1", "secondary_fixed_gauge_H")
ALLOWED = {-1: {(0, 0)}, 0: {(0, 0), (1, 1), (1, 2), (0, 1)}, 1: {(1, 1)}}
HASHES = {
    "axis": {"old": "2d23fecc98d276d9ad15ad1867199cd308f0570cb5040ef94eb6b923b4c53458",
             "seam": "d589cf29ec770292674e510018d849bf44f10fb08b882d38ad30ba468db06127"},
    "tilted": {"old": "225031e612929ed922ba75c55e76703d59990f5283e7ac39b94f022841798da5",
               "seam": "107f5081734d9f5d3259303b28a6f71083f5822d297ceed96f73348ac1ec1d5d"},
}


def git_commit(value):
    return subprocess.check_output(["git", "rev-parse", value+"^{commit}"], cwd=ROOT, text=True).strip()


def iid_coefficients(rows):
    result = {key: [F(0)]*(N+1) for key in ("z", "q", "e")}
    for row in rows:
        k, q, count = row["k"], row["q"], row["count"]
        for key, multiplier in (("z", 1), ("q", q), ("e", q*q)):
            result[key][k] += count*multiplier
    return result


def trace_coefficients(rows, old_rows):
    result = {key: [F(0)]*(N+1) for key in PACKETS}
    direct = {key: {o: [F(0)]*(N+1) for o in ("q", "e")} for key in PACKETS}
    marginal, support = defaultdict(int), defaultdict(int)
    rejected = []
    for row in rows:
        k, g, q, count = (row[key] for key in ("k", "g", "q", "count"))
        pattern = (row["bad2"], row["n_bad3"])
        marginal[k, g, q] += count
        if q not in ALLOWED or pattern not in ALLOWED[q]:
            rejected.append(row)
            continue
        kind = "A" if q == 0 and pattern == (1, 2) else "B" if q == 0 and pattern == (0, 1) else "zero"
        support[kind] += count
        beta = F(-1) if kind in ("A", "B") else F(0)
        H = F(k+g+(3 if kind == "A" else 1), 2) if kind in ("A", "B") else F(0)
        for key, value in zip(PACKETS, (beta, H)):
            result[key][k] += count*value
            direct[key]["q"][k] += count*q*value
            direct[key]["e"][k] += count*q*q*value
    expected = {(r["k"], r["g"], r["q"]): r["count"] for r in old_rows}
    if dict(marginal) != expected:
        raise ValueError("seam archive is not the locked complete occupation population")
    if any(value for packet in direct.values() for coeffs in packet.values() for value in coeffs):
        raise ArithmeticError("direct q/E trace numerator is not identically zero")
    return result, {"unsupported_patterns": rejected, "support_configuration_counts": dict(support),
                    "occupation_marginal_retained": True, "direct_q_E_numerators": "identically_zero"}


def pack_score(result, geometry, scale):
    value = result["value"]
    return {
        "decision": "nonzero_normalization_transmission" if value.lo > 0 or value.hi < 0 else "zero_not_excluded_stop_fixed_score",
        "response_over_A": interval_json(value), "response_approx": scale(value),
        "three_terms_over_A": {k: interval_json(v) for k, v in result["terms"].items()},
        "three_terms_approx": {k: scale(v) for k, v in result["terms"].items()},
        "transmission_coefficients": {k: interval_json(v) for k, v in result["coefficients"].items()},
        "trace_arguments": {k: interval_json(v) for k, v in result["arguments"].items()},
        "root_h_tangent": interval_json(result["root_h_derivative"]),
        "geometry_trace_fraction": [{"f": interval_json(p["f"][0]), "f_h": interval_json(p["f"][1])} for p in geometry],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-dir", type=Path, default=ROOT/"results/p337-s4-trace-transmission")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--packing-gate-commit", required=True)
    parser.add_argument("--character-gate-commit", required=True)
    args = parser.parse_args()
    start = time.perf_counter()
    contract = json.loads(CONTRACT.read_text())
    if contract["N"] != N or contract["geometries"] != [[5, 0], [4, 3]]:
        raise ValueError("only the predeclared Q1 continuation is implemented")
    theory = {"packing_commit": git_commit(args.packing_gate_commit),
              "character_commit": git_commit(args.character_gate_commit),
              "authorization": "run only after coordinator confirms both proofs; commit IDs record the completed gates"}
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    inputs, bases, traces, gates = [], [], [], []
    for label in ("axis", "tilted"):
        old_path, path = OLD/f"{label}.csv", args.counts_dir/f"{label}.csv"
        if sha(old_path) != HASHES[label]["old"] or sha(path) != HASHES[label]["seam"]:
            raise ValueError("the fixed old/seam population hash differs")
        old_rows, rows = read_rows(old_path), read_rows(path)
        trace, gate = trace_coefficients(rows, old_rows)
        bases.append(iid_coefficients(old_rows)); traces.append(trace)
        gates.append({"geometry": label, **gate})
        inputs.append({"geometry": label, "old_path": str(old_path.relative_to(ROOT)),
                       "old_sha256": sha(old_path), "seam_path": str(path.resolve()), "seam_sha256": sha(path),
                       "coefficients_y": {key: [str(v) for v in vals] for key, vals in trace.items()}})
    shared = {"schema": "p337.q1-trace-continuation.score.v1", "freeze_commit": FREEZE,
              "code_commit": git_commit("HEAD"), "contract": contract, "contract_sha256": sha(CONTRACT),
              "theory_gates": theory, "finite_pattern_gates": gates, "inputs": inputs,
              "created_utc": datetime.now(timezone.utc).isoformat(),
              "new_enumerations": 0, "new_samples": 0, "root_searches": 0, "new_Q_or_seam_scan": False}
    if any(gate["unsupported_patterns"] for gate in gates):
        shared["status"] = "not_scored_unsupported_mod6_pattern"
        (out/"score.json").write_text(json.dumps(shared, indent=2)+"\n")
        print(shared["status"])
        return

    saved = json.loads(BASELINE.read_text())
    p = from_json(saved["root_enclosure"])
    h = p/(1-p)
    Dp = from_json(saved["matching_slope_enclosure"])
    Dh = Dp/((1+h)**2)
    adapted = {"D_h": interval_json(Dh), "U_over_A": saved["U25_over_A_enclosure"]}
    outputs = {}
    with localcontext() as context:
        context.prec = 60
        area = Decimal(N)**(Decimal(13)/Decimal(8))/2
        def scale(value):
            v = middle(value)
            return float(area*Decimal(v.numerator)/Decimal(v.denominator))
        for key in PACKETS:
            geometry = [geometry_packet(base, trace[key], h) for base, trace in zip(bases, traces)]
            outputs[key] = pack_score(normalization_map(geometry, adapted), geometry, scale)
    shared.update({"status": "completed_fixed_two_packet_score", "outputs": outputs,
        "primary_scope": "epsilon insertion beta1=-1_A-1_B into the actual iid Q1 occupation law",
        "secondary_scope": "raw-trace Q-derivative packet in y^K Q^(-(K+g)/2) gauge; additive attribution, not a gauge-invariant share or mixed epsilon-Q derivative",
        "root": {"source": str(BASELINE.relative_to(ROOT)), "sha256": sha(BASELINE),
                 "imported_p": saved["root_enclosure"], "converted_h_y": interval_json(h),
                 "imported_Dp": saved["matching_slope_enclosure"], "converted_Dh": adapted["D_h"],
                 "imported_U_over_A": saved["U25_over_A_enclosure"],
                 "coordinate_conversion": "h=p/(1-p); D_h=D_p/(1+h)^2; U/A unchanged"},
        "reused_map_sha256": sha(ROOT/"scripts/p337_s4_trace_transmission_score.py"),
        "source_score_or_baseline_root_recomputed": False,
        "boundary": "two prescribed views of one exact population; no regular endpoint activation, continuum field identity or independent second evidence block",
        "elapsed_seconds": time.perf_counter()-start})
    (out/"score.json").write_text(json.dumps(shared, indent=2)+"\n")
    primary, secondary = (outputs[key] for key in PACKETS)
    (out/"REPORT.md").write_text(
        "# Fixed Q1 closed-trace continuation\n\n"
        f"Primary: **{primary['decision']}**, response {primary['response_approx']:.16g}.\n\n"
        f"Secondary fixed-gauge raw-trace packet: {secondary['response_approx']:.16g}. "
        "This is additive attribution in the declared reduced partition; it is neither a gauge-invariant share nor a mixed epsilon-Q derivative.\n\n"
        "Both direct q/E numerators are exactly zero. The imported iid root, original D and U/A remain fixed; only the prescribed normalization/root/slope transmission is scored. "
        "The interval results and separate common-thermal, geometric-thermal and geometric-value terms are in score.json. "
        "No enumeration, root search, Q scan, old-source rescore or new random sample was used.\n")
    print(json.dumps({"primary": primary["response_approx"], "primary_decision": primary["decision"],
                      "secondary_fixed_gauge": secondary["response_approx"],
                      "elapsed_seconds": shared["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
