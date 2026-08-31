#!/usr/bin/env python3
"""Complete the fixed weighted-jump-only test, using an exact 1/8 subset.

Only two new cross-moments; saved endpoint/defect marginals and root reused.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time

import analyze_decimation_plaquette_u as base

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/defect_reweight_contract.json"
CPP = ROOT / "scripts/exact_defect_reweight_cross_moments.cpp"
I = base.Interval
DEGREE = 25


def blob(commit, path):
    return subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)


def unpack(data):
    return list(csv.DictReader(io.StringIO(data.decode())))


def interval(item):
    return I(F(item["lower_fraction"]), F(item["upper_fraction"]))


def jet(values, p):
    # Sums contain multiplicities. The four pinned neighbors do not change
    # the Bernoulli degree: all weights remain p^K*(1-p)^(25-K).
    poly = [0]*(DEGREE+1)
    for k, value in enumerate(values):
        for j in range(DEGREE-k+1):
            poly[k+j] += value*math.comb(DEGREE-k, j)*(-1)**j
    result = []
    for d in range(3):
        value = I.of(0)
        for k in range(DEGREE, d-1, -1):
            value = value*p+poly[k]*math.comb(k, d)*math.factorial(d)
        result.append(value)
    return result


def column(rows, key):
    return [int(r[key]) for r in rows]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    started_utc = datetime.now(timezone.utc).isoformat()
    contract = json.loads(CONTRACT.read_text())
    code_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    command = ["clang++", "-O3", "-std=c++17", str(CPP)]
    with tempfile.TemporaryDirectory(prefix="matching-defect-cross-") as tmp:
        binary = Path(tmp)/"defect-cross"
        subprocess.run(command+["-o", str(binary)], check=True, capture_output=True, text=True)

        def run(item):
            name, (a, b) = item
            path = output/(name+".csv")
            proc = subprocess.run([str(binary), str(a), str(b), str(path)],
                                  check=True, capture_output=True, text=True)
            return json.loads(proc.stdout)

        with ThreadPoolExecutor(max_workers=2) as pool:
            runs = list(pool.map(run, zip(("first", "second"), contract["parent_geometries"])))
    enumeration_seconds = time.perf_counter()-started
    source_receipts = []

    def source(commit, path):
        data = blob(commit, path)
        source_receipts.append({"commit": commit, "path": path,
                                "sha256": hashlib.sha256(data).hexdigest()})
        return data

    child_result = json.loads(source(contract["baseline_commit"], contract["baseline_directory"]+"/latest.json"))
    root = interval(child_result["root_enclosure"])
    root = I(1-root.hi, 1-root.lo)
    full_result = json.loads(source(contract["defect_result_commit"], contract["defect_directory"]+"/score/score.json"))
    xi_full = interval(full_result["rational_enclosures"]["Xi_U_t_epsilon_over_A"])
    packet = []
    for child_name, defect_name in zip(("axis", "tilted"), ("first", "second")):
        child = unpack(source(contract["baseline_commit"], contract["baseline_directory"]+"/"+child_name+".csv"))
        defect = unpack(source(contract["defect_result_commit"], contract["defect_directory"]+"/"+defect_name+".csv"))
        subset = unpack((output/(defect_name+".csv")).read_bytes())
        if any(len(x) != 26 for x in (child, defect, subset)):
            raise ValueError("all K=0..25 rows are required")
        for k in range(26):
            expected = 2*math.comb(21, k-2) if 2 <= k <= 23 else 0
            if (int(subset[k]["k"]) != k or int(subset[k]["count"]) != expected
                or int(defect[k]["count"]) != math.comb(25, k)
                or int(child[k]["count"]) != math.comb(25, k)):
                raise ValueError("incomplete conditional or full configuration population")
        plus = {f: [(-1 if f in ("q", "qsstar") else 1)*v
                    for v in reversed(column(child, "sum_"+f))]
                for f in ("q", "e", "sstar", "qsstar", "esstar")}
        minus_s = column(defect, "sum_sstar")
        sj = jet([a-b for a, b in zip(minus_s, plus["sstar"])], root)
        row = {"q": jet(plus["q"], root), "e": jet(plus["e"], root)}
        for obs in ("q", "e"):
            cross = [int(d["sum_"+obs+"sstar"])
                     +int(a["sum_sminus_"+obs+"plus"])-int(a["sum_"+obs+"sstar"])
                     for d, a in zip(defect, subset)]
            joint = jet([v-w for v, w in zip(cross, plus[obs+"sstar"])], root)
            mu = row[obs]
            cov = joint[0]-sj[0]*mu[0]
            cov_p = joint[1]-sj[1]*mu[0]-sj[0]*mu[1]
            row["h"+obs] = [25*(1-root)*cov, 25*((1-root)*cov_p-cov)]
        packet.append(row)
    delta = F(contract["delta_cos4"])
    D = sum(r["q"][1] for r in packet)/2
    Qpp = sum(r["q"][2] for r in packet)/2
    Yp = (packet[0]["e"][1]-packet[1]["e"][1])/delta
    Ypp = (packet[0]["e"][2]-packet[1]["e"][2])/delta
    hQ = sum(r["hq"][0] for r in packet)/2
    hQp = sum(r["hq"][1] for r in packet)/2
    hYp = (packet[0]["he"][1]-packet[1]["he"][1])/delta
    terms = {"direct": hYp/D, "root_motion": -Ypp*hQ/D**2,
             "slope_source": -Yp*hQp/D**2, "slope_root": Yp*Qpp*hQ/D**3}
    xi_reweight = sum(terms.values())
    quantities = {"Xi_total": xi_full, "Xi_reweight": xi_reweight,
                  "Xi_weighted_jump": xi_full-xi_reweight}
    bounds = {k: base.interval_json(v) for k, v in quantities.items()}
    decision = ("weighted_rank_jump_only_rejected" if bounds["Xi_reweight"]["excludes_zero"]
                else "weighted_rank_jump_only_not_excluded_by_this_projection")
    with localcontext() as context:
        context.prec = 65
        area = Decimal(50)**(Decimal(13)/8)/2
        values = {k: float(area*base.decimal_fraction(base.middle(v))) for k, v in quantities.items()}
    result = {"schema": "matching-one.defect-reweight.v1", "status": "completed_exact_subset_reconstruction",
              "decision": decision, "contract": contract, "code_commit": code_commit,
              "numerical_values": values, "enclosures_over_A": bounds,
              "reweight_terms_over_A": {k: base.interval_json(v) for k, v in terms.items()},
              "root_enclosure": base.interval_json(root), "source_files": source_receipts,
              "population": {"new_subset_configurations_per_geometry": 2**22,
                             "full_population_per_geometry": 2**25,
                             "complement_cross_correction": "exactly_zero_by_alternating_face_rank_surgery",
                             "new_random_samples": 0, "new_root_searches": 0},
              "subset_files": [{"path": name+".csv", "sha256": base.sha(output/(name+".csv"))}
                               for name in ("first", "second")],
              "dependency": "same paired exact finite populations as the prior full Xi, no independent evidence vote"}
    (output/"latest.json").write_text(json.dumps(result, indent=2)+"\n")
    report = f"""# Baseline reweighting is a measured part of the one-hole U response

Decision: **{decision}** on the fixed N50 parent pair.

| Mixed original-U contribution | Exact coefficients evaluated numerically |
|---|---:|
| Baseline reweighting | {values['Xi_reweight']:+.15g} |
| Weighted observable jump | {values['Xi_weighted_jump']:+.15g} |
| Total, imported from the completed one-hole result | {values['Xi_total']:+.15g} |

The two prescribed terms oppose each other: baseline reweighting partially
offsets the negative weighted jump. The one-term model misses a real positive
contribution; the full normalized defect operator already contains it, without
an added source or an adjustable mixing coefficient.

The reweighting contribution has a rational enclosure excluding zero:
{bounds['Xi_reweight']['excludes_zero']}. This is the prescribed contribution
of Cov(w,O_intact), not a newly fitted residual or a rank-preserving population
share. Both rank-changing and rank-preserving configurations can contribute.
The jump contribution is total minus reweighting, so the two add exactly by
definition and do not count as independent evidence.

## How the missing information was obtained

The prior full defect packet lacks Sminus*qplus and Sminus*Eplus. The exact
single-defect topology restricts q/E changes to alternating four-neighbor
patterns. For every other configuration Oplus=Ominus; hence

`full sum(Sminus*Oplus) = old full sum(Sminus*Ominus)`
`+ alternating sum[Sminus*(Oplus-Ominus)]`.

Only 2 alternating patterns times 2^21 remaining bits were enumerated per
geometry: one eighth of the full population. Intact and defective observers
use identical B configurations. The Bernoulli degree stays25, including the
two forced occupied and two forced vacant neighbors. No baseline enumeration,
random sampling, root search or test suite was performed. Raw subset sums,
input commit/hash pointers and the exact arithmetic result are included.

## Fixed observer, source and chart

Source Sstar=C+F4+Bvac; pA=s+(1-s)p, pB=p, epsilon=1-s. At zero source the
baseline-reweighting insertion vanishes identically for every p. Its mixed
jet is hO=25(1-p)Cov(Sminus-Splus,Oplus). Each geometry is normalized before
pooling; the saved complementary root is reused. All four terms in
Xi_reweight/A = hY_p/D - Y_pp*hQ/D^2 - Y_p*hQ_p/D^2 + Y_p*Q_pp*hQ/D^3
are included, with D=Q_p and A=50^(13/8)/2. The derivative of25(1-p) is retained.

This decides the fixed weighted-rank-jump-only response model; it does not
identify a continuum field, finite-interior law or an independent production
effect. The earlier source-independent gain rejection, larger-N F4 unresolved
result and P154/P334 fixed decisions remain unchanged. The two contributions
are not two new adjustable sources. Do not fit another relative coefficient
to restore a failed one-term model.

Definition and mechanism algebra: notes/decimation-closed-source-and-global-u.md
at7132f0c2. Operator proof: bc17b81d:notes/checkerboard-single-defect-source.md.
Prior total Xi: f5c4a74a:results/p337-endpoint-defect/score/score.json.
Run: python scripts/analyze_defect_reweight.py --output-dir NEW_DIRECTORY.
"""
    (output/"REPORT.md").write_text(report)
    receipt = {"started_utc": started_utc, "finished_utc": datetime.now(timezone.utc).isoformat(),
               "elapsed_seconds": time.perf_counter()-started, "enumeration_and_compile_seconds": enumeration_seconds,
               "command": sys.argv, "compile_command": command+["-o", "TEMP_BINARY"],
               "enumerations": runs, "python": sys.version, "machine": platform.machine(),
               "code_commit": code_commit, "hashes": {"script": base.sha(__file__), "cpp": base.sha(CPP),
                   "helper": base.sha(base.__file__), "contract": base.sha(CONTRACT),
                   "result": base.sha(output/"latest.json")},
               "cloud_jobs": 0, "new_random_samples": 0, "science_tests_run": 0}
    (output/"run.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"decision": decision, "values": values, "elapsed_seconds": receipt['elapsed_seconds']}), flush=True)


if __name__ == "__main__":
    main()
