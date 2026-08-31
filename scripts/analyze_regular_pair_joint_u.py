#!/usr/bin/env python3
"""Compute one fixed canonical two-insertion original-U response, exactly.

New origin-vacant crossmoments only; existing full q/E population and root.
The first-Q two-insertion tensor source is not a mark covariance.
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

from analyze_decimation_plaquette_u import Interval as I, interval_json, middle

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/regular_pair_joint_u_contract.json"
CPP = ROOT / "scripts/regular_pair_joint_u_exact.cpp"
N = 25
DELTA = F(1152, 625)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def interval(item):
    return I(F(item["lower_fraction"]), F(item["upper_fraction"]))


def rows(data):
    return [{key: int(value) for key, value in row.items()}
            for row in csv.DictReader(io.StringIO(data.decode()))]


def raw_jet(coefficients, h):
    values = []
    for order in range(3):
        value = I.of(0)
        for k in range(N, order-1, -1):
            value = value*h + coefficients[k]*math.comb(k, order)*math.factorial(order)
        values.append(value)
    return values


def normalized_jet(raw, z):
    value = raw[0]/z[0]
    first = (raw[1]-value*z[1])/z[0]
    second = (raw[2]-value*z[2]-2*first*z[1])/z[0]
    return [value, first, second]


def coefficients(old_rows, new_rows):
    data = {key: [F(0)]*(N+1) for key in ("z", "q", "e", "s", "qs", "es")}
    for row in old_rows:
        k, q, count = (row[key] for key in ("k", "q", "count"))
        for key, value in (("z", 1), ("q", q), ("e", q*q)):
            data[key][k] += count*value
    if data["z"] != [math.comb(N, k) for k in range(N+1)]:
        raise ValueError("the imported full baseline population is incomplete")
    if len(new_rows) != N+1:
        raise ValueError("the source profile must contain every K=0..25")
    for k, row in enumerate(new_rows):
        expected = math.comb(N-1, k) if k < N else 0
        if row["k"] != k or row["count"] != expected:
            raise ValueError("the prescribed origin-vacant exact population is incomplete")
        # Deterministic accounting in the same traversal, not a second test run:
        # translation invariance fixes the one-vacancy marginal of q and E.
        for key, field in (("q", "sum_q"), ("e", "sum_e")):
            if N*row[field] != (N-k)*data[key][k]:
                raise ValueError("origin-vacant q/E marginal differs from the existing baseline")
        for key, field in (("s", "sum_G16"), ("qs", "sum_G16_q"), ("es", "sum_G16_E")):
            data[key][k] = F(row[field], 16*N)
    return data


def response(pair, h, D, ratio):
    packets = []
    for data in pair:
        z = raw_jet(data["z"], h)
        packets.append({key: normalized_jet(raw_jet(values, h), z)
                        for key, values in data.items() if key != "z"})
    source = []
    for row in packets:
        q, qh, _ = row["q"]
        e, eh, _ = row["e"]
        s, sh, _ = row["s"]
        qs, qsh, _ = row["qs"]
        es, esh, _ = row["es"]
        source.append({"mean_S2": s, "mean_S2_h": sh,
                       "j_q": qs-q*s, "j_q_h": qsh-qh*s-q*sh,
                       "j_e": es-e*s, "j_e_h": esh-eh*s-e*sh})
    jq = sum(row["j_q"] for row in source)/2
    jqh = sum(row["j_q_h"] for row in source)/2
    jyh = (source[0]["j_e_h"]-source[1]["j_e_h"])/DELTA
    mh2 = sum(row["q"][2] for row in packets)/2
    yh2 = (packets[0]["e"][2]-packets[1]["e"][2])/DELTA
    terms = {"direct_centered": jyh/D, "root_motion": -yh2*jq/D**2,
             "slope_source": -ratio*jqh/D, "slope_root": ratio*mh2*jq/D**2}
    return {"J2_over_A": sum(terms.values()), "terms": terms,
            "root_h_joint_tangent": -jq/D, "source_packets": source}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    started_utc = datetime.now(timezone.utc).isoformat()
    contract_bytes = CONTRACT.read_bytes()
    contract = json.loads(contract_bytes)
    if contract["N"] != N or contract["geometries"] != [[5, 0], [4, 3]]:
        raise ValueError("only the fixed N25 comparison is implemented")
    code_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    source_files = []

    def source(commit, path):
        data = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        source_files.append({"commit": commit, "path": path, "sha256": sha(data)})
        return data

    kernel = source(contract["kernel"]["commit"], contract["kernel"]["path"])
    if sha(kernel) != contract["kernel"]["sha256"]:
        raise ValueError("the pinned signed joint kernel differs")
    baseline = contract["baseline"]
    baseline_bytes = source(baseline["commit"], baseline["root_path"])
    if sha(baseline_bytes) != baseline["root_sha256"]:
        raise ValueError("the imported original root packet differs")
    saved = json.loads(baseline_bytes)
    p = interval(saved["root_enclosure"])
    h = p/(1-p)
    D = interval(saved["matching_slope_enclosure"])/(1+h)**2
    ratio = interval(saved["U25_over_A_enclosure"])
    if D.lo <= 0:
        raise ArithmeticError("the saved matching slope is not strictly positive")
    old = {name: rows(source(baseline["commit"], baseline["counts_directory"]+"/"+name+".csv"))
           for name in ("axis", "tilted")}
    with tempfile.TemporaryDirectory(prefix="matching-regular-joint-") as tmp:
        temporary = Path(tmp)
        kernel_path = temporary/"kernel.tsv"
        kernel_path.write_bytes(kernel)
        binary = temporary/"joint-exact"
        compile_command = ["clang++", "-O3", "-std=c++17", str(CPP), "-o", str(binary)]
        subprocess.run(compile_command, check=True, capture_output=True, text=True)
        compiler = subprocess.check_output(["clang++", "--version"], text=True).splitlines()[0]

        def run(item):
            name, (a, b) = item
            command = [str(binary), str(a), str(b), str(kernel_path), str(out/(name+".csv"))]
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            receipt = json.loads(result.stdout)
            receipt.update({"geometry": name, "command": command, "output": name+".csv",
                            "output_sha256": sha((out/(name+".csv")).read_bytes())})
            print(json.dumps({"phase": "new_joint_moments_completed", **receipt}), flush=True)
            return receipt

        with ThreadPoolExecutor(max_workers=2) as pool:
            runs = list(pool.map(run, zip(("axis", "tilted"), contract["geometries"])))
    scored_start = time.perf_counter()
    pair = [coefficients(old[name], rows((out/(name+".csv")).read_bytes()))
            for name in ("axis", "tilted")]
    scored = response(pair, h, D, ratio)
    score_seconds = time.perf_counter()-scored_start
    bound = interval_json(scored["J2_over_A"])
    decision = "global_additive_first_Q_closure_rejected" if bound["excludes_zero"] else "unresolved_at_saved_root_precision"
    with localcontext() as context:
        context.prec = 65
        area = Decimal(N)**(Decimal(13)/8)/2

        def full(x):
            value = middle(x)
            return float(area*Decimal(value.numerator)/Decimal(value.denominator))

        values = {"J2": full(scored["J2_over_A"]),
                  "terms": {key: full(value) for key, value in scored["terms"].items()}}
    result = {"schema": "matching-one.regular-pair-joint-u.v1", "status": "completed_fixed_exact_joint_response",
              "definition_commit": code_commit, "contract": contract, "decision": decision,
              "J2_over_A": bound, "numerical_values": values,
              "terms_over_A": {key: interval_json(value) for key, value in scored["terms"].items()},
              "root_p_joint_tangent": interval_json(scored["root_h_joint_tangent"]/(1+h)**2),
              "source_packets": [{key: interval_json(value) for key, value in row.items()} for row in scored["source_packets"]],
              "coefficient_packets": [{"geometry": name, "S2_coefficients_h": {key: [str(v) for v in data[key]] for key in ("s", "qs", "es")}}
                                      for name, data in zip(("axis", "tilted"), pair)],
              "source_files": source_files, "saved_root_p": saved["root_enclosure"],
              "source_unit_divisor": 16*N, "normalization_population": 2**N,
              "new_origin_vacant_configurations_per_geometry": 2**(N-1),
              "new_random_samples": 0, "new_root_searches": 0, "old_sources_rescored": False,
              "response_scores": 1, "tests_run": 0, "cloud_jobs": 0,
              "score_seconds": score_seconds, "elapsed_seconds": time.perf_counter()-started,
              "boundary": contract["boundary"]}
    report = ["# The canonical joint Q activation in original U", "",
              f"Decision: **{decision}**. **J2={values['J2']:+.16g}**.", "",
              "J2=partial_logQ partial_epsilon^2 U at Q1,epsilon0 for fixed Kreg=K2+K0.",
              "Each vacant-site vertex uses epsilon/N; the original separately normalized q/E, pooled root and thermal-slope quotient are retained.", "",
              "| Complete original-U term | Value |", "|---|---:|"]
    report += [f"| {key} | {value:+.16g} |" for key, value in values["terms"].items()]
    report += ["", "The first-Q effective logweight of a linear additive model predicts J2=0 exactly.",
               "This is a joint tensor closure, not Cov(a_x,a_y), the conditional 13/8 contraction or the spatial C at fixed separation.",
               "Both adjacent and nonadjacent vertices enter. Adjacent vacant marks share one physical isolated edge ID.", "",
               "One exact 2^24 origin-vacant traversal per geometry supplies only the missing joint source moments.",
               "Translation-invariant moments use sum_y g16_0y/(16N) and the full 2^25 baseline, not a conditional probability law.",
               "Old root and q/E populations were imported; no old source score, new root search, Monte Carlo or test campaign.",
               "The finite N25 populations are shared with earlier exact work, not independent evidence or a continuum field assignment.", "",
               "See [latest.json](latest.json) for exact rational enclosures and [run.json](run.json) for commands and hashes.", ""]
    contents = {"latest.json": json.dumps(result, indent=2)+"\n", "REPORT.md": "\n".join(report)}
    for name, text in contents.items():
        (out/name).write_text(text)
    receipt = {"schema": "matching-one.regular-pair-joint-u.run.v1", "definition_commit": code_commit,
               "started_utc": started_utc, "finished_utc": datetime.now(timezone.utc).isoformat(),
               "command": [sys.executable, *sys.argv], "python": sys.version, "machine": platform.machine(),
               "compiler": compiler, "compile_command": compile_command, "compile_count": 1,
               "geometry_runs": runs, "source_files": source_files,
               "contract_sha256": sha(contract_bytes), "script_sha256": sha(Path(__file__).read_bytes()),
               "producer_sha256": sha(CPP.read_bytes()),
               "interval_helper_sha256": sha((ROOT/"scripts/analyze_decimation_plaquette_u.py").read_bytes()),
               "output_sha256": {name: sha(text.encode()) for name, text in contents.items()},
               "elapsed_seconds": result["elapsed_seconds"], "score_seconds": score_seconds,
               "new_random_samples": 0, "old_scores_rerun": 0, "new_root_searches": 0, "tests_run": 0, "cloud_jobs": 0}
    (out/"run.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"decision": decision, **values, "elapsed_seconds": result["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
