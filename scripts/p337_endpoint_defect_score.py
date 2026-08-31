#!/usr/bin/env python3
"""Fixed single-hole test of source-independent gain in the original global U.

Only consumes complete integer profiles. No sampling, endpoint re-enumeration,
source fitting or finite differences. All reported sign decisions use rational
enclosures, with the positive irrational area factor applied for display only.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

from p337_closed_source_score import Interval as I, interval_json, middle

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/p337-closed-source-n25"
DEGREE = 25
PARENT_N = 50
DELTA = F(-1152, 625)
ORDER = 3
FIELDS = ("count", "q", "e", "sstar", "qsstar", "esstar")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_counts(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != DEGREE + 1:
        raise ValueError("requires every free-B occupancy K=0..25")
    data = {f: [] for f in FIELDS}
    for k, row in enumerate(rows):
        if int(row["k"]) != k or int(row["count"]) != math.comb(DEGREE, k):
            raise ValueError("incomplete fixed finite configuration population")
        for field in FIELDS:
            data[field].append(int(row[field if field == "count" else "sum_"+field]))
    return data


def parent_baseline(child):
    # K_B(parent)=25-K(child); source is fixed, q complements, E does not.
    return {field: [(-1 if field in ("q", "qsstar") else 1)*v
                    for v in reversed(values)] for field, values in child.items()}


def power_coefficients(values):
    # Enumerated sums already contain their binomial multiplicities.
    out = [0]*(DEGREE+1)
    for k, value in enumerate(values):
        for j in range(DEGREE-k+1):
            out[k+j] += value*math.comb(DEGREE-k, j)*(-1)**j
    return out


def evaluate_jet(values, root):
    poly = power_coefficients(values)
    jet = []
    for derivative in range(ORDER+1):
        value = I.of(0)
        for k in range(DEGREE, derivative-1, -1):
            value = value*root + poly[k]*math.comb(k, derivative)
        jet.append(value)
    return jet  # Taylor coefficients: derivative/factorial(derivative)


def add(a, b):
    return [x+y for x, y in zip(a, b)]


def scale(a, c):
    return [c*x for x in a]


def sub(a, b):
    return add(a, scale(b, -1))


def multiply(a, b):
    return [sum((a[j]*b[k-j] for j in range(k+1)), I.of(0))
            for k in range(ORDER+1)]


def normalized_partials(base, defect, root):
    h0 = {f: evaluate_jet(v, root) for f, v in base.items()}
    hd = {f: evaluate_jet(v, root) for f, v in defect.items()}
    # x=s-1: derivative of raw mixture is -25*(1-p)*(Hd-H0).
    dose_s = [-DEGREE*(1-root), I.of(DEGREE), I.of(0), I.of(0)]
    z_t = h0["sstar"]
    z_st = multiply(dose_s, sub(hd["sstar"], h0["sstar"]))
    result = {}
    for obs in ("q", "e"):
        baseline = h0[obs]
        obs_s = multiply(dose_s, sub(hd[obs], baseline))
        obs_t = sub(h0[obs+"sstar"], multiply(baseline, z_t))
        raw_st = multiply(dose_s, sub(hd[obs+"sstar"], h0[obs+"sstar"]))
        obs_st = sub(sub(raw_st, multiply(obs_s, z_t)), multiply(baseline, z_st))
        result[obs] = {(0, 0): baseline, (1, 0): obs_s,
                       (0, 1): obs_t, (1, 1): obs_st}
    return result


def project(pair, obs, weights):
    return {key: [weights[0]*pair[0][obs][key][p]+weights[1]*pair[1][obs][key][p]
                  for p in range(ORDER+1)]
            for key in pair[0][obs]}


def partial(packet, p=0, s=0, t=0):
    return math.factorial(p)*packet[(s, t)][p]


def comoving_score(q, y):
    D = partial(q, p=1)
    B = partial(y, p=1)
    if D.lo <= 0:
        raise ArithmeticError("the endpoint pooled root must have a positive slope")
    ps = -partial(q, s=1)/D
    pt = -partial(q, t=1)/D
    pst = -(partial(q, s=1, t=1)+partial(q, p=1, s=1)*pt
            +partial(q, p=1, t=1)*ps+partial(q, p=2)*ps*pt)/D

    def moving_slope(packet):
        value = partial(packet, p=1)
        ds = partial(packet, p=1, s=1)+partial(packet, p=2)*ps
        dt = partial(packet, p=1, t=1)+partial(packet, p=2)*pt
        dst = (partial(packet, p=1, s=1, t=1)
               +partial(packet, p=2, s=1)*pt
               +partial(packet, p=2, t=1)*ps
               +partial(packet, p=3)*ps*pt+partial(packet, p=2)*pst)
        return value, ds, dt, dst

    _, Bs, Bt, Bst = moving_slope(y)
    _, Ds, Dt, Dst = moving_slope(q)
    f = B/D
    fs = Bs/D-B*Ds/(D**2)
    ft = Bt/D-B*Dt/(D**2)
    fst_terms = {
        "numerator_mixed": Bst/D,
        "numerator_s_denominator_t": -Bs*Dt/(D**2),
        "numerator_t_denominator_s": -Bt*Ds/(D**2),
        "denominator_mixed": -B*Dst/(D**2),
        "denominator_product": 2*B*Ds*Dt/(D**3),
    }
    fst = sum(fst_terms.values(), I.of(0))
    gain_residual = f*fst-fs*ft
    return {
        "U_over_A": f, "U_s_over_A": fs, "U_t_over_A": ft,
        "U_st_over_A": fst, "Xi_U_t_epsilon_over_A": -fst,
        "R_over_A_squared": gain_residual,
        "gain_residual_over_A": gain_residual/f,
        "d_s_Ut_over_U": gain_residual/(f**2),
        "root_s": ps, "root_t": pt, "root_st": pst,
        "matching_slope": D, "mixed_ratio_terms": fst_terms,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--defect-axis", type=Path, required=True)
    parser.add_argument("--defect-tilted", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    utc = datetime.now(timezone.utc).isoformat()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    child_paths = [BASE/"axis.csv", BASE/"tilted.csv"]
    defect_paths = [args.defect_axis.resolve(), args.defect_tilted.resolve()]
    base_score = json.loads((BASE/"latest.json").read_text())
    saved_root = base_score["root_enclosure"]
    root = I(1-F(saved_root["upper_fraction"]), 1-F(saved_root["lower_fraction"]))
    child = [read_counts(p) for p in child_paths]
    defect = [read_counts(p) for p in defect_paths]
    pair = [normalized_partials(parent_baseline(c), d, root) for c, d in zip(child, defect)]
    q = project(pair, "q", (F(1, 2), F(1, 2)))
    y = project(pair, "e", (1/DELTA, -1/DELTA))
    score = comoving_score(q, y)
    enclosures = {k: interval_json(v) for k, v in score.items() if k != "mixed_ratio_terms"}
    enclosures["mixed_ratio_terms"] = {k: interval_json(v) for k, v in score["mixed_ratio_terms"].items()}
    with localcontext() as context:
        context.prec = 70
        A = Decimal(PARENT_N)**(Decimal(13)/Decimal(8))/2

        def approx(value, factor=Decimal(1)):
            mid = middle(value)
            return float(factor*Decimal(mid.numerator)/Decimal(mid.denominator))

        numerical = {"p0_endpoint": approx(root)}
        for name in ("U", "U_s", "U_t", "U_st", "Xi_U_t_epsilon"):
            numerical[name] = approx(score[name+"_over_A"], A)
        numerical["R"] = approx(score["R_over_A_squared"], A*A)
        numerical["d_s_Ut_over_U"] = approx(score["d_s_Ut_over_U"])
        numerical["source_free_gain_slope"] = numerical["U_s"]/numerical["U"]
        numerical["gain_predicted_U_st"] = numerical["U_s"]*numerical["U_t"]/numerical["U"]
        numerical["U_st_minus_gain_prediction"] = approx(score["gain_residual_over_A"], A)
    primary = ("source_independent_endpoint_gain_rejected" if
               enclosures["R_over_A_squared"]["excludes_zero"] else
               "source_independent_endpoint_gain_unresolved")
    secondary = ("interior_thermal_only_mixed_null_rejected" if
                 enclosures["Xi_U_t_epsilon_over_A"]["excludes_zero"] else
                 "interior_thermal_only_mixed_null_unresolved")
    result = {
        "schema": "matching-one.p337-endpoint-defect.score.v1",
        "status": "completed_exact_finite_one_defect",
        "decision": primary, "secondary_decision": secondary,
        "freeze_commit": subprocess.check_output(["git", "rev-parse", "9024fdbf"], cwd=ROOT, text=True).strip(),
        "execution_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "model": "bulk exp(t*(C+F+Bvac)), pA=s+(1-s)p, pB=p",
        "parent_geometries": [[5, 5], [1, 7]], "parent_N": PARENT_N,
        "free_B_sites": DEGREE, "delta_cos4": str(DELTA),
        "root_source": "complement of the prior complete N25 root enclosure",
        "root_enclosure": interval_json(root), "numerical_values": numerical,
        "rational_enclosures": enclosures,
        "inputs": [{"path": str(p), "sha256": sha(p)} for p in child_paths+defect_paths+[BASE/"latest.json"]],
        "score_sha256": sha(__file__),
        "normalization": "each geometry first; Z_t and Z_st retained; root and thermal slope comove",
        "dependence": "one exact N50 single-defect population paired with the fixed N25 endpoint; no independent evidence votes",
        "boundary": "finite endpoint tangent only, no full interior curve or continuum field; original F4/P154/P334 stops unchanged",
        "new_random_samples": 0, "tests_run": 0,
    }
    (out/"score.json").write_text(json.dumps(result, indent=2)+"\n")
    rows = "\n".join(f"| {key} | {value:.15g} |" for key, value in numerical.items())
    report = f"""# One saturation defect: the closed source and original global U

Primary decision: **{primary}**.
Secondary fixed question: **{secondary}**.

| Endpoint quantity | Numerical evaluation of exact coefficients |
|---|---:|
{rows}

The source, graph pair and original global-U normalization are fixed.
U_s and U_st use s increasing toward full saturation; epsilon=1-s, so
Xi=U_t,epsilon=-U_st. The primary R=U*U_st-U_s*U_t eliminates an unknown
source-independent geometric gain. Decisions use the rational enclosures
of R/A50^2 and Xi/A50, not the rounded displayed values. Every ratio term,
root displacement and slope displacement is included in `score.json`.

One removed A site represents all25 positions by translation. Each geometry
enumerates all2^25 free-B configurations exactly. Its Bernstein degree is25.
The intact endpoint uses the prior N25 coefficients by complement, with no
new baseline enumeration or resampling. Source normalizers are retained
separately for each geometry before the pooled root and angular projection.
The p-dependent defect dose25(1-p) is differentiated, not frozen at the root.

This is a finite mechanistic equality test, not a confidence interval or
a continuum field identification. A nonzero R excludes the local scalar-gain
extension, not the exact saturated identity. No extra source, fitted curve,
defect class, sample extension or old-experiment rescue enters this result.
The original F4/P154/P334 stop decisions remain unchanged.

Frozen specification: `notes/checkerboard-endpoint-defect-decision-freeze.md`
at9024fdbf. See the producer contract and receipts in the parent result folder.
"""
    (out/"REPORT.md").write_text(report)
    receipt = {"started_utc": utc, "finished_utc": datetime.now(timezone.utc).isoformat(),
               "elapsed_seconds": time.perf_counter()-started, "command": sys.argv,
               "python": sys.version, "score_sha256": sha(out/"score.json"),
               "new_random_samples": 0, "tests_run": 0, "cloud_jobs": 0}
    (out/"run.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"decision": primary, "secondary_decision": secondary,
                      "values": numerical, "elapsed_seconds": receipt["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
