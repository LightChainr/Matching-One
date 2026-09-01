#!/usr/bin/env python3
"""Compare two named Q tangents by exact linear reduction of saved bounds."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

from analyze_decimation_plaquette_u import Interval, interval_json, middle

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/weak_q_path_comparison_contract.json"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def saved_bound(item):
    return Interval(F(item["lower_fraction"]), F(item["upper_fraction"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    started_at = datetime.now(timezone.utc).isoformat()
    contract = json.loads(CONTRACT.read_text())
    raw = subprocess.check_output(["git", "show", contract["source_commit"]+":"+contract["source_path"]], cwd=ROOT)
    source = json.loads(raw)
    if source["contract"]["N"] != contract["N"] or source["contract"]["geometries"] != contract["geometries"]:
        raise ValueError("the two responses must belong to the fixed common finite ensemble")
    inputs = [saved_bound(source["source_enclosures"][name]["V_over_A"]) for name in ("sstar", "bv")]
    output = {}
    with localcontext() as ctx:
        ctx.prec = 65
        area = Decimal(contract["N"])**(Decimal(13)/8)/2
        for name, row in zip(contract["output_order"], contract["fixed_matrix_rows"]):
            bound = sum(value*F(coefficient) for value, coefficient in zip(inputs, row))
            mid = middle(bound)
            output[name] = {"response_over_A": interval_json(bound),
                            "response_approx": float(area*Decimal(mid.numerator)/Decimal(mid.denominator)),
                            "sign": "positive" if bound.lo > 0 else "negative" if bound.hi < 0 else "unresolved"}
    opposite = output["tied_Q_path"]["sign"] == "positive" and output["rank_projected_site_RC_Q_path"]["sign"] == "negative"
    code_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    result = {"schema": "matching-one.weak-q-path-comparison.v1",
              "status": "completed_exact_saved_interval_linear_reduction",
              "contract": contract, "code_commit": code_commit,
              "input": {"commit": contract["source_commit"], "path": contract["source_path"], "sha256": digest(raw)},
              "shared_root_enclosure": source["root_enclosure"],
              "responses": output,
              "decision": "named_Q_paths_have_opposite_original_U_tangents" if opposite else "opposite_sign_prediction_unresolved",
              "dependence": "same exact N25 populations and root-complete source responses; interval subtraction is enclosing, not statistical independence",
              "boundary": "rank_projected_ordinary_site_RC_not_unprojected_site_RC_no_cross_size_or_field_identity_claim",
              "new_random_samples": 0, "new_enumerations": 0, "root_searches": 0}
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out/"latest.json").write_text(json.dumps(result, indent=2)+"\n")
    rows = "\n".join(f"| {name} | {v['response_approx']:+.15g} | {v['sign']} |" for name, v in output.items())
    report = f"""# Two named Q paths already have opposite original-U tangents

Decision: **{result['decision']}** on the same N25 iid baseline.

| Root-adjusted derivative at Q=1 | Response | Exact sign |
|---|---:|---|
{rows}

The identity is `V_Sstar/2 = V_(C_B-r/2) + V_B`, with
`V_B=V_Bvac` on the thermal quotient. The first term is the Q tangent of
the rank-projected ordinary site-RC law; the left side is the fixed tied-edge
family. Both use the same rank projection and original q/E observer.
The factor1/2 converts the pre-existing t derivative to log-Q, since Q=exp(2t).
At Q=1 the log-Q and Q derivatives coincide. Each path follows its perturbed
pooled root, retaining its individual geometry normalizers and slope motion.

The local edge control reverses the response sign; it cannot be removed by
a common density reparameterization. This excludes equality of these two
named finite-observer tangents, not Potts universality or a continuum field.
It does not compare the tied law with *unprojected* ordinary site-RC.

Only the saved rational source enclosures were transformed. No profile,
root or scorer was rerun. Both source values were known before this reduction;
this is an exact algebraic consequence, not a prospective independent test.
Full reduced-response bounds, common-root provenance and the fixed matrix
are stored in latest.json. The area factor is displayed numerically only.

The ordinary endpoint's regular Q-activation exclusion and the remaining
trace/confluent interface are explained in
[the mechanism note](../../notes/weak-q-paths-and-regular-selection.md).
Neither this N25 sign comparison nor its source control establishes a
cross-size logarithmic velocity. P154/P334/F4 production decisions stand.

Reproduce: python scripts/analyze_weak_q_path_comparison.py --output-dir NEW_DIRECTORY.
"""
    (out/"REPORT.md").write_text(report)
    receipt = {"started_utc": started_at, "finished_utc": datetime.now(timezone.utc).isoformat(),
               "elapsed_seconds": time.perf_counter()-started, "command": sys.argv,
               "python": sys.version, "machine": platform.machine(), "code_commit": code_commit,
               "hashes": {"script": digest(Path(__file__).read_bytes()), "contract": digest(CONTRACT.read_bytes()),
                          "interval_helper": digest((ROOT/"scripts/analyze_decimation_plaquette_u.py").read_bytes()),
                          "result": digest((out/"latest.json").read_bytes())},
               "new_enumerations": 0, "new_random_samples": 0, "root_searches": 0,
               "science_tests": 0, "cloud_jobs": 0}
    (out/"run.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"decision": result["decision"], "responses": {k:v["response_approx"] for k,v in output.items()},
                      "elapsed_seconds": receipt["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
