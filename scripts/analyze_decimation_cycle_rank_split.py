#!/usr/bin/env python3
"""Remove the closed cycle-gas source's explicit rank bias at original global U.

One exact reduction of saved N25 coefficients. No enumeration or sampling.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import io
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

import analyze_decimation_plaquette_u as base

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/decimation_cycle_rank_contract.json"


def saved_interval(item):
    return base.Interval(F(item["lower_fraction"]), F(item["upper_fraction"]))


def git_blob(commit, path):
    return subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    started_at = datetime.now(timezone.utc).isoformat()
    contract = json.loads(CONTRACT.read_text())
    counts_dir = ROOT / contract["baseline_directory"]
    baseline_path = counts_dir / "score/latest.json"
    baseline = json.loads(baseline_path.read_text())
    external_commit = contract["closed_source_commit"]
    external_dir = contract["closed_source_directory"]
    external_bytes = git_blob(external_commit, external_dir + "/latest.json")
    external = json.loads(external_bytes)
    assert external["contract"]["geometries"] == contract["geometries"]
    assert external["contract"]["delta_cos4"] == contract["delta_cos4"]
    pair, input_profiles = [], []
    for name in ("axis.csv", "tilted.csv"):
        data = base.profiles(counts_dir/name)
        other_bytes = git_blob(external_commit, external_dir + "/" + name)
        other = list(csv.DictReader(io.StringIO(other_bytes.decode())))
        # This coefficient identity establishes a common ensemble/root for the
        # linear response subtraction; it is not an independent replication.
        if len(other) != 26:
            raise ValueError("closed-source baseline lacks a full N25 profile")
        for k, row in enumerate(other):
            if int(row["k"]) != k:
                raise ValueError("external K order differs")
            for field in ("count", "q", "e"):
                if int(row[field if field == "count" else "sum_"+field]) != data[field][k]:
                    raise ValueError("the two baseline coefficient profiles differ")
        pair.append(data)
        input_profiles.append({"name": name, "local_sha256": base.sha(counts_dir/name),
                               "external_sha256": hashlib.sha256(other_bytes).hexdigest(),
                               "q_E_and_count_coefficients_identical": True})
    a, b = saved_interval(baseline["root_enclosure"]), saved_interval(external["root_enclosure"])
    root = base.Interval(max(a.lo, b.lo), min(a.hi, b.hi))
    if root.lo > root.hi:
        raise ValueError("saved exact root enclosures do not overlap")
    packet = [base.moments(data, root) for data in pair]
    D = sum(row["q"][1] for row in packet)/2
    Qpp = sum(row["q"][2] for row in packet)/2
    Yp = (packet[0]["e"][1]-packet[1]["e"][1])/base.DELTA
    Ypp = (packet[0]["e"][2]-packet[1]["e"][2])/base.DELTA
    jq, jqp, jep = [], [], []
    for row in packet:
        q, qp, _ = row["q"]
        e, ep, _ = row["e"]
        jq.append(e-q*q)
        jqp.append(ep-2*q*qp)
        jep.append(qp*(1-e)-q*ep)
    jQ, jQp = sum(jq)/2, sum(jqp)/2
    jYp = (jep[0]-jep[1])/base.DELTA
    terms = {"direct": jYp/D, "root_motion": -Ypp*jQ/(D**2),
             "slope_source": -Yp*jQp/(D**2), "slope_root": Yp*Qpp*jQ/(D**3)}
    vq = sum(terms.values())
    vs = saved_interval(external["source_enclosures"]["sstar"]["V_over_A"])
    residual = vs-vq
    quantities = {"V_q": vq, "V_Sstar": vs, "two_V_beta_null": residual,
                  "V_beta_null": residual/2, "two_V_beta1": vs+vq}
    bounds = {name: base.interval_json(x) for name, x in quantities.items()}
    decision = ("fixed_rank_bias_only_U_alias_rejected" if bounds["two_V_beta_null"]["excludes_zero"]
                else "rank_bias_only_U_comparison_unresolved")
    with localcontext() as ctx:
        ctx.prec = 55
        area = Decimal(25)**(Decimal(13)/8)/2
        values = {name: float(area*base.decimal_fraction(base.middle(x))) for name, x in quantities.items()}
        values["U25"] = float(area*base.decimal_fraction(base.middle(Yp/D)))
    result = {"schema": "matching-one.decimation-cycle-rank.v1",
              "status": "completed_saved_exact_coefficient_reduction", "contract": contract,
              "code_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "source_profile_identity": input_profiles, "root_enclosure": base.interval_json(root),
              "source_enclosures_over_A": bounds,
              "q_source_terms_over_A": {name: base.interval_json(x) for name, x in terms.items()},
              "q_source_root_tangent": base.interval_json(-jQ/D),
              "numerical_values": values, "decision": decision,
              "interpretation": "the residual is the fixed twice-ambient-null-cycle response, not regression remainder or independent evidence",
              "new_configurations": 0, "new_random_samples": 0, "root_searches": 0}
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out/"latest.json").write_text(json.dumps(result, indent=2)+"\n")
    report = f"""# Is the closed source only an explicit rank bias?

## The fixed rank-bias-only prediction has an exact answer

Decision: **`{decision}`** on the fixed N25 Gaussian pair.
The closed source has original-U response {values['V_Sstar']:+.12g}; its
explicit q source contributes {values['V_q']:+.12g}. Their difference,
{values['two_V_beta_null']:+.12g}, is exactly twice the ambient-null graph-cycle
source response. The rational enclosure of that difference/A
{'excludes' if bounds['two_V_beta_null']['excludes_zero'] else 'contains'} zero.
This answers the fixed prediction V_Sstar=V_q, without choosing a new rank
coefficient or altering the source after observing the split.

## A coefficient-fixed cycle source, not a free residual

The configuration identity is
`Sstar=2*beta_null+q-3*K+2*N+2`, with q=ambient rank−1 and
`beta_null=dim ker[H1(occupied NN graph;R) -> H1(torus;R)]`.
The common K source is a Bernoulli-odds reparameterization: its derivative
of root/slope-normalized U is zero. Consequently
`V_Sstar=2*V_beta_null+V_q` in the same bulk exp(tS) units.
The computed V_beta_null is {values['V_beta_null']:+.12g}.
The alternate graph-cycle basis gives 2*V_beta1={values['two_V_beta1']:+.12g};
this is one algebraic change of basis, not a second independent experiment.

Ambient-null graph cycles include zero-winding combinations of nontrivial
cycles. Beta_null is not a count of elementary full faces or the cellular
hole count of a filled-cell complex. Its source is a specified finite graph
statistic; no claim of a local continuum operator follows.

## The original observer and root normalization are retained

The two geometries are (5,0) and (4,3), DeltaCos4=1152/625,
`U=A*Y_p/Q_p` at the pooled Q=0 root, A=25^(13/8)/2.
Per geometry the exact rank algebra q^3=q gives
`j_q=E-q_mean^2` and `j_E=q_mean*(1-E)`.
Their p derivatives, source root motion and both denominator corrections
give the same four-term original-U response formula as the parent reader.
All terms and rational enclosures are saved in latest.json.

Only the saved integer coefficient profiles and the committed Sstar response
were consumed. Their q/E/count coefficient arrays agree exactly, establishing
the same finite ensembles and root. The saved root bracket was reused;
there was no enumeration, sampling, production replay or new root search.

## What this separates and what remains open

The claim is one exact finite-observer source separation. These deterministic
N25 calculations are coordinates of the same complete configuration sets,
not independent statistical votes. The axis Z5xZ5 and tilted Z25 quotients
have different Smith classes; no large-N amplitude or continuum mechanism is
identified. The completed independent F4 block's inconclusive decision is
unchanged, as are the P154/P334 source-specific decisions.

This excludes only the explicit unit-coefficient q alias. One scalar response
cannot exclude a post-hoc fitted c*q source. A nonzero beta_null contribution
can still act through its occupancy/rank conditional means; it does not imply
that a source centered within every K/rank becomes visible to global U.

This equilibrium decomposition does not itself compute interior Xi. The later
execution delivery f5c4a74a20bad8589c39e1034cfb209462110dbe,
`results/p337-endpoint-defect/score/REPORT.md`, now completes that separate
calculation: Xi=-10.755718407564073 and R=U*U_st−U_s*U_t=27.766563581230237
have nonzero rational enclosures. The fixed source-independent gain model
and mixed thermal-only null fail. That delivered result supersedes the earlier
scorer-only snapshot; no duplicate defect run is the next assignment.

## Sources and reproduction

The fixed action comes from {contract['action_commit']}:
{contract['action_path']}. The complete source response is pinned at
{external_commit}:{external_dir}/latest.json. The baseline is
{contract['baseline_commit']}:{contract['baseline_directory']}.
The exact term values are reported as a short list rather than a trend chart:
there is one finite-pair comparison and no estimated scaling series.

Run `python scripts/analyze_decimation_cycle_rank_split.py --output-dir NEW_DIRECTORY`.
"""
    (out/"REPORT.md").write_text(report)
    receipt = {"started_utc": started_at, "finished_utc": datetime.now(timezone.utc).isoformat(),
               "elapsed_seconds": time.perf_counter()-started, "command": sys.argv,
               "python": sys.version, "machine": platform.machine(),
               "hashes": {"contract": base.sha(CONTRACT), "script": base.sha(__file__),
                          "helper": base.sha(base.__file__), "baseline_result": base.sha(baseline_path),
                          "external_result": hashlib.sha256(external_bytes).hexdigest(),
                          "result": base.sha(out/"latest.json")},
               "new_enumerations": 0, "new_random_samples": 0, "cloud_jobs": 0,
               "root_searches": 0, "science_tests_run": 0}
    (out/"run.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"decision": decision, "values": values,
                      "elapsed_seconds": receipt["elapsed_seconds"]}))


if __name__ == "__main__":
    main()
