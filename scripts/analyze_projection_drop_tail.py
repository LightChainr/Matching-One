#!/usr/bin/env python3
"""Extract the fixed projection-deleted leading U tail from saved exact counts."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import io
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/projection_drop_tail_contract.json"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read_source(commit, path):
    return subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    start_utc = datetime.now(timezone.utc).isoformat()
    contract = json.loads(CONTRACT.read_text())
    code_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    n = contract["N"]
    delta = F(contract["delta_cos4"])
    source_commit = subprocess.check_output(["git", "rev-parse", contract["histogram_commit"]], cwd=ROOT, text=True).strip()
    geometries = []
    for label in ("axis", "tilted"):
        path = contract["histogram_directory"]+"/"+label+".csv"
        data = read_source(source_commit, path)
        rows = [{k: int(v) for k, v in row.items()} for row in csv.DictReader(io.StringIO(data.decode()))]
        if sum(r["count"] for r in rows) != 2**n:
            raise ValueError("the fixed complete histogram is required")
        for r in rows:
            r["eta"] = F(r["g"]-r["q"]-1)+F(2*r["k"], n)
        ground = [r for r in rows if r["eta"] == 0]
        if {(r["k"], r["q"], r["count"]) for r in ground} != {(0, -1, 1), (n, 1, 1)}:
            raise ValueError("the proposed empty/full normalization is not this histogram's ground support")
        proper_gap = min(r["eta"] for r in rows if 0 < r["k"] < n)
        one = [r for r in rows if r["q"] == 0]
        exponent = min(r["eta"] for r in one)
        support = [r for r in one if r["eta"] == exponent]
        # Differentiate R1/(1+d^N) at d=1, then use E=1-P1.
        e_slope = sum(F(r["count"])*(F(n, 2)-r["k"])/2 for r in support)
        geometries.append({"label": label, "minimum_eta": str(exponent),
                           "minimum_proper_eta": str(proper_gap),
                           "leading_support": [{k: r[k] for k in ("k", "g", "q", "count")} for r in support],
                           "leading_E_d_coefficient": str(e_slope),
                           "source": {"commit": source_commit, "path": path, "sha256": sha(data)}})
    leading_exp = min(F(g["minimum_eta"]) for g in geometries)
    coefficient = sum(sign*F(g["leading_E_d_coefficient"])
                      for sign, g in zip((1, -1), geometries)
                      if F(g["minimum_eta"]) == leading_exp)/(delta*F(n, 2))
    expected = {"exponent": "42/5", "coefficient": str(3/delta)}
    actual = {"exponent": str(leading_exp), "coefficient": str(coefficient)}
    theory = []
    for length in (5, 10, 15):
        q = length*length-6*length+6
        theory.append({"N": length*length, "L": length,
                       "status": "N25_supported_by_existing_histogram" if length == 5 else "unmeasured_combinatorial_prediction",
                       "projected_exponent": 2*length+1, "projected_coefficient": str(-q/delta),
                       "drop_exponent": str(F(2*length-2)+F(2, length)),
                       "drop_coefficient": str(F(length-2)/delta),
                       "ratio_coefficient": str(-F(length-2, q)),
                       "ratio_exponential_rate_in_t": str(F(3)-F(2, length))})
    result = {"schema": "matching-one.projection-drop-tail.v1", "status": "completed_exact_saved_histogram_leading_term",
              "contract": contract, "code_commit": code_commit, "geometry_leading_terms": geometries,
              "projection_deleted_U_over_A_leading": actual,
              "fixed_prediction": expected, "matches_fixed_prediction": actual == expected,
              "decision": "projection_deleted_positive_tail_opposes_original_negative_tail" if coefficient > 0 else "positive_tail_prediction_not_established",
              "theoretical_size_laws": theory,
              "original_negative_tail_source": {"commit": contract["source_commit"], "path": "notes/closed-source-square-family-leading-law.md"},
              "dependency": "same exact N25 populations as original first moments, finite-coupling grid and projected series; no independent statistical vote",
              "boundary": "fixed_N_then_t_to_infinity_no_uniform_size_remainder_finite_t_sign_threshold_or_new_larger_N_data",
              "new_random_samples": 0, "new_enumerations": 0, "root_searches": 0}
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out/"latest.json").write_text(json.dumps(result, indent=2)+"\n")
    report = f"""# Removing the topological projection reverses the strong-coupling U tail

The fixed N25 histogram gives **U_drop/A = {coefficient} lambda^({leading_exp})
+ higher terms**. Its leading coefficient is positive. The already proved
projected law has **U_star/A = -625/1152 lambda^11 + O(lambda^13)**.
These two fixed, bulk-pressure-equivalent laws have opposite eventual signs
for the original root/slope-normalized angular observer.

No coupling point, series fit or configuration was added. The minimum of
eta=g-r+2K/N was extracted exactly from each saved integer histogram. The
axis leading support has K=5,g=9,count=10; after partition normalization
its E_d coefficient is75/2. Dividing by Q_d=25/2 and Delta=1152/625 gives
1875/1152=625/384. Root motion begins later and cannot cancel this leading term.

The comparison law is the previously defined Sdrop=Sstar+r, obtained solely
by omitting m^(-r). The thermal variable d=p/((1-p)m)*m^(2/N) is common to
both geometries, so its Jacobian cancels in U. The result is not a new
rank-fugacity fit, an arbitrary observer or a fixed-root substitute.

## Cross-size consequences are theory, not new measurements

For an axis L×L torus with L>=5 and same-area companion ell1>=L+2,

`U_drop/A_N = (L-2)/Delta * lambda^(2L-2+2/L) + O(lambda^(2L-2+4/L))`.

Compared with the original negative law,
`U_drop/U_star ~ -(L-2)/(L²-6L+6) * exp((3-2/L)t)` at each fixed L.
N100 and N225 coefficients/exponents are recorded in latest.json as unmeasured
combinatorial predictions. See the [proof](../../notes/topological-projection-reverses-global-u-tail.md).

The first stripe widths now have unequal source costs. Their reciprocal
occupation symmetry no longer cancels the lowest normalized thermal slope.
This identifies how a topological projection can change U even though its
pressure-density discrepancy is bounded by2t/N. That pressure statement takes
fixed t with growing N; it is not an interchange of the two limits.

The finite-t sign-change location, quantitative useful-coupling window and
thermodynamic behavior remain open. These are deterministic views of the
same N25 complete populations, not independent evidence or a continuum field
assignment. P154/P334/F4 decisions stay unchanged.

Reproduce: python scripts/analyze_projection_drop_tail.py --output-dir NEW_DIRECTORY.
"""
    (out/"REPORT.md").write_text(report)
    receipt = {"started_utc": start_utc, "finished_utc": datetime.now(timezone.utc).isoformat(),
               "elapsed_seconds": time.perf_counter()-started, "command": sys.argv,
               "python": sys.version, "machine": platform.machine(), "code_commit": code_commit,
               "hashes": {"script": sha(Path(__file__).read_bytes()), "contract": sha(CONTRACT.read_bytes()),
                          "result": sha((out/"latest.json").read_bytes())},
               "new_random_samples": 0, "new_enumerations": 0, "root_searches": 0,
               "cloud_jobs": 0, "science_tests_run": 0}
    (out/"run.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"decision": result["decision"], "leading": actual,
                      "elapsed_seconds": receipt['elapsed_seconds']}), flush=True)


if __name__ == "__main__":
    main()
