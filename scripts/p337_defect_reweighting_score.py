#!/usr/bin/env python3
"""One fixed exact test of the defect baseline-reweighting mixed-U contribution.

Only the two absent cross moments are newly enumerated. The endpoint root,
intact/defect marginals and complete U_st enclosure are committed inputs.
This never re-evaluates the previous Xi or gain residual R.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time

from p337_endpoint_defect_score import (
    DEGREE, PARENT_N, DELTA, read_counts, parent_baseline, evaluate_jet,
    sub, scale, multiply,
)
from p337_closed_source_score import Interval as I, interval_json, middle

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT/"analysis/p337_defect_reweighting_contract.json"
CPP = ROOT/"scripts/p337_defect_reweighting_exact.cpp"
BASE = ROOT/"results/p337-closed-source-n25"
DEFECT = ROOT/"results/p337-endpoint-defect"


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read_cross(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != DEGREE+1:
        raise ValueError("cross profile must retain all free-B counts")
    result = {"q": [], "e": []}
    for k, row in enumerate(rows):
        if int(row["k"]) != k or int(row["count"]) != math.comb(DEGREE, k):
            raise ValueError("incomplete paired configuration population")
        result["q"].append(int(row["sum_Sdef_qintact"]))
        result["e"].append(int(row["sum_Sdef_Eintact"]))
    return result


def cross_covariances(base, defect, cross, root):
    h0 = {f: evaluate_jet(v, root) for f, v in base.items()}
    sd = evaluate_jet(defect["sstar"], root)
    delta_s = sub(sd, h0["sstar"])
    packet = {"covariance": {}, "uncentered": {}, "centering": {}}
    for obs in ("q", "e"):
        raw = sub(evaluate_jet(cross[obs], root), h0[obs+"sstar"])
        centering = scale(multiply(delta_s, h0[obs]), -1)
        packet["uncentered"][obs] = raw
        packet["centering"][obs] = centering
        packet["covariance"][obs] = sub(raw, scale(centering, -1))
    return h0, packet


def project_jets(pair, obs):
    weights = (F(1, 2), F(1, 2)) if obs == "q" else (1/DELTA, -1/DELTA)
    return [weights[0]*pair[0][obs][k]+weights[1]*pair[1][obs][k]
            for k in range(4)]


def linear_mixed_contribution(covariances, ordinary, root):
    q0, y0 = project_jets(ordinary, "q"), project_jets(ordinary, "e")
    cq, cy = project_jets(covariances, "q"), project_jets(covariances, "e")
    D, B, T, H = q0[1], y0[1], 2*q0[2], 2*y0[2]
    if D.lo <= 0:
        raise ArithmeticError("the inherited endpoint root has no positive slope")
    dose, dose_p = -DEGREE*(1-root), DEGREE  # s increases toward saturation
    gq, gy = dose*cq[0], dose*cy[0]
    gqp, gyp = dose*cq[1]+dose_p*cq[0], dose*cy[1]+dose_p*cy[0]
    terms = {
        "fixed_p_E_mixed_jet": gyp/D,
        "root_motion_E_curvature": -H*gq/(D**2),
        "slope_mixed_jet": -B*gqp/(D**2),
        "root_motion_slope": B*T*gq/(D**3),
    }
    dose_split = {
        "thermal_derivative_of_covariance": dose*(cy[1]/D-B*cq[1]/(D**2)),
        "explicit_dose_derivative": dose_p*(cy[0]/D-B*cq[0]/(D**2)),
        "root_motion": -H*gq/(D**2)+B*T*gq/(D**3),
    }
    return {"U_st_over_A": sum(terms.values(), I.of(0)),
            "root_st_contribution": -gq/D, "terms": terms, "dose_split": dose_split,
            "gQ": gq, "gQp": gqp, "gY": gy, "gYp": gyp, "D": D}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    started, utc = time.perf_counter(), datetime.now(timezone.utc).isoformat()
    contract = json.loads(CONTRACT.read_text())
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    binary = Path(tempfile.mkdtemp(prefix="p337-defect-reweighting-"))/"enumerate"
    command = ["/usr/bin/clang++", "-O3", "-std=c++17", str(CPP), "-o", str(binary)]
    subprocess.run(command, check=True)
    binary_hash = sha(binary)

    def enumerate_one(item):
        geometry, name = item
        path = out/f"{name}.csv"
        cmd = [str(binary), *(str(x) for x in geometry), str(path)]
        completed = subprocess.run(cmd, check=True, capture_output=True, text=True)
        receipt = json.loads(completed.stdout)
        receipt.update(geometry=geometry, file=path.name, sha256=sha(path),
                       command=cmd, exit_code=completed.returncode,
                       paired_intact_and_defect=True)
        print(json.dumps(receipt), flush=True)
        return receipt

    with ThreadPoolExecutor(max_workers=contract["workers"]) as pool:
        receipts = list(pool.map(enumerate_one, zip(contract["geometries"], ("first", "second"))))
    child_paths = [BASE/"axis.csv", BASE/"tilted.csv"]
    defect_paths = [DEFECT/"first.csv", DEFECT/"second.csv"]
    cross_paths = [out/"first.csv", out/"second.csv"]
    base_score_path = BASE/"latest.json"
    saved_root = json.loads(base_score_path.read_text())["root_enclosure"]
    root = I(1-F(saved_root["upper_fraction"]), 1-F(saved_root["lower_fraction"]))
    ordinary, packets = [], []
    for c, d, x in zip(child_paths, defect_paths, cross_paths):
        h0, packet = cross_covariances(parent_baseline(read_counts(c)), read_counts(d), read_cross(x), root)
        ordinary.append(h0)
        packets.append(packet)
    scored = {name: linear_mixed_contribution([x[name] for x in packets], ordinary, root)
              for name in ("covariance", "uncentered", "centering")}
    main_score = scored["covariance"]
    imported_path = DEFECT/"score/score.json"
    imported = json.loads(imported_path.read_text())["rational_enclosures"]["U_st_over_A"]
    previous_full = I(F(imported["lower_fraction"]), F(imported["upper_fraction"]))
    weighted_rank = previous_full-main_score["U_st_over_A"]
    with localcontext() as ctx:
        ctx.prec = 70
        amplitude = Decimal(PARENT_N)**(Decimal(13)/Decimal(8))/2

        def numeric(value, factor=amplitude):
            mid = middle(value)
            return float(factor*Decimal(mid.numerator)/Decimal(mid.denominator))

        numerical = {
            "endpoint_p0": numeric(root, Decimal(1)),
            "U_st_reweighting": numeric(main_score["U_st_over_A"]),
            "Xi_reweighting_epsilon": -numeric(main_score["U_st_over_A"]),
            "weighted_rank_jump_component_from_previous_total": numeric(weighted_rank),
            "previous_full_U_st_imported_not_recomputed": numeric(previous_full),
            "uncentered_cross_difference_U_st": numeric(scored["uncentered"]["U_st_over_A"]),
            "normalization_centering_U_st": numeric(scored["centering"]["U_st_over_A"]),
            "root_st_reweighting_contribution": numeric(main_score["root_st_contribution"], Decimal(1)),
        }
        numerical["dose_split_U_st"] = {key: numeric(value) for key, value in main_score["dose_split"].items()}
    enclosure = interval_json(main_score["U_st_over_A"])
    result = {
        "schema": "matching-one.p337-defect-reweighting.score.v1",
        "status": "completed_exact_fixed_operator_split",
        "decision": "weighted_rank_jump_only_model_rejected" if enclosure["excludes_zero"] else
                    "weighted_rank_jump_only_model_not_resolved",
        "contract": contract, "contract_freeze_commit": "e6a900d9d644b26278f01c17bdfb6f27f3903b75",
        "code_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "root_enclosure": interval_json(root), "U_st_reweighting_over_A": enclosure,
        "weighted_rank_jump_component_over_A": interval_json(weighted_rank),
        "imported_previous_full_U_st_over_A": imported,
        "complete_root_slope_terms_over_A": {k: interval_json(v) for k, v in main_score["terms"].items()},
        "alternative_dose_split_over_A": {k: interval_json(v) for k, v in main_score["dose_split"].items()},
        "alternative_covariance_centering_split_over_A": {
            name: interval_json(scored[name]["U_st_over_A"]) for name in ("uncentered", "centering")},
        "raw_mixed_mean_jets": {k: interval_json(main_score[k]) for k in ("gQ", "gQp", "gY", "gYp")},
        "per_geometry_deltaS_covariance_jets": [
            {obs: [interval_json(jet) for jet in packet["covariance"][obs][:2]]
             for obs in ("q", "e")} for packet in packets],
        "numerical_values": numerical,
        "inputs": [{"path": str(path.relative_to(ROOT)), "sha256": sha(path)}
                   for path in child_paths+defect_paths+cross_paths+[base_score_path, imported_path]],
        "scope": "signed operator contribution with exact paired configurations; neither a population share nor a fitted percentage",
        "uncertainty": "rational computational bounds, no sampling error; conservative arithmetic over the saved root enclosure",
        "old_Xi_or_R_recomputed": False, "new_random_samples": 0, "cloud_jobs": 0, "tests_run": 0,
    }
    (out/"score.json").write_text(json.dumps(result, indent=2)+"\n")
    run = {
        "started_utc": utc, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter()-started, "enumerations": receipts,
        "compile_command": command, "binary_sha256": binary_hash,
        "source_sha256": {"contract": sha(CONTRACT), "cpp": sha(CPP), "scorer": sha(Path(__file__)),
                          "jet_backend": sha(ROOT/"scripts/p337_endpoint_defect_score.py"),
                          "interval_backend": sha(ROOT/"scripts/p337_closed_source_score.py")},
        "python": sys.version, "machine": platform.machine(), "command": sys.argv,
        "result_sha256": sha(out/"score.json"), "new_random_samples": 0, "cloud_jobs": 0,
        "tests_run": 0,
    }
    (out/"run.json").write_text(json.dumps(run, indent=2)+"\n")
    print(json.dumps({"decision": result["decision"], "values": numerical,
                      "elapsed_seconds": run["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()
