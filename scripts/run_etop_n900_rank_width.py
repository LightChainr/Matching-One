#!/usr/bin/env python3
"""Collect and score the third-scale intrinsic rank-width experiment."""
import argparse
import csv
from fractions import Fraction
import hashlib
import json
import subprocess

import numpy as np
from scipy.stats import chi2

import run_etop_n100_three_modulus as engine
from etop_thermal_transport import describe

ROOT = engine.ROOT
OUT = ROOT/"results/etop-n900-rank-width"
CONTRACT = ROOT/"experiments/etop_n900_rank_width_20260831.json"


def rank_moments(contract):
    n, batches = contract["area"], contract["batches"]
    per_batch = contract["samples_per_shape_pair"]//batches
    coefficients = np.zeros((batches, 2, n+1))
    for s, shape in enumerate(contract["shapes"]):
        delta = float(Fraction(shape["delta_cos4"]))
        with (OUT/"raw"/(shape["name"]+".hist.csv")).open() as stream:
            for row in csv.DictReader(stream):
                sign = 1 if row["orientation"] == "first" else -1
                coefficients[int(row["batch"]), s, int(row["k"])] += sign*int(row["count"])/(per_batch*delta)
    steps = np.cumsum(coefficients[:, 1]-coefficients[:, 0], axis=-1)[:, :n]
    low, high = np.arange(n)/n, (np.arange(n)+1)/n
    kernels = np.array([(high**(j+1)-low**(j+1))/(j+1) for j in range(3)]).T
    batch_moments = steps@kernels
    mean = batch_moments.mean(axis=0)
    loo = (batches*mean-batch_moments)/(batches-1)
    def value(m):
        mu = m[1]/m[0]
        return np.array([m[0], mu, (m[2]/m[0]-mu*mu)*n**.75])
    return {"labels": ["signed_step_area", "mean_p", "centered_second_z"],
            **describe(value(mean), np.array([value(v) for v in loo])),
            "batch_raw_moments": batch_moments.tolist()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--freeze", required=True)
    args = parser.parse_args()
    contract = json.loads(CONTRACT.read_text())
    engine.OUT, engine.FREEZE = OUT, args.freeze
    if args.run:
        engine.launch(contract)
    batches = np.column_stack([engine.read_batch_vectors(contract, s) for s in contract["shapes"]])
    mean = batches.mean(axis=0)
    covariance = np.cov(batches, rowvar=False, ddof=1)/len(batches)
    ranks = rank_moments(contract)
    source_bytes = subprocess.check_output(["git", "show", contract["source_prediction_commit"]+":"+contract["source_prediction_path"]], cwd=ROOT)
    source = json.loads(source_bytes)["post_reveal_working_fingerprint"]
    predictions = np.asarray(source["estimate"])[1:]
    prediction_cov = np.asarray(source["covariance"])[1:, 1:]
    observed, target_var = ranks["mean"][2], ranks["covariance"][2][2]
    residual = observed-predictions
    comparison_cov = prediction_cov+target_var
    comparisons = []
    for i, name in enumerate(("quarter_power_width", "fixed_critical_width_profile")):
        variance = comparison_cov[i, i]
        comparisons.append({"name": name, "prediction": float(predictions[i]),
            "prediction_se": float(np.sqrt(prediction_cov[i, i])), "target_minus_prediction": float(residual[i]),
            "total_se": float(np.sqrt(variance)), "z": float(residual[i]/np.sqrt(variance)),
            "nominal_p_value": float(chi2.sf(residual[i]**2/variance, 1))})
    anchor, anchor_var = predictions[1], prediction_cov[1, 1]
    gamma = None
    if observed > 0:
        denom = 2*np.log(900/400)
        gamma = {"estimate": float(.375-np.log(observed/anchor)/denom),
                 "se": float(np.sqrt(target_var/observed**2+anchor_var/anchor**2)/denom),
                 "scope": "Finite N400-to-N900 signed-profile effective width; not an asymptotic exponent."}
    result = {"contract": contract, "prediction_freeze_commit": args.freeze,
        "shape_order": [s["name"] for s in contract["shapes"]], "field_order": engine.FIELDS,
        "mean": mean.tolist(), "covariance": covariance.tolist(), "batch_vectors": batches.tolist(),
        "rank_profile": ranks, "source_prediction_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "prediction_comparisons": comparisons, "comparison_covariance": comparison_cov.tolist(),
        "effective_width_400_900": gamma,
        "raw_sha256": {str(p.relative_to(ROOT)): engine.digest(p) for p in sorted((OUT/"raw").iterdir())},
        "boundary": contract["interpretation"]}
    engine.dump(OUT/"score.json", result)
    lines = ["# N900: third-scale intrinsic rank-clock width", "",
        "32M new shared counters, 800 aligned batches, two modulus pairs. This block is independent of N100/N400; predictions retain their common N400 anchor covariance.", "",
        f"Measured rank-step centered z variance: {observed:.10g} +/- {np.sqrt(target_var):.6g}.", "",
        "| conditional prediction | expected Vz | observed minus expected | total SE | z | nominal p |",
        "|---|---:|---:|---:|---:|---:|"]
    for row in comparisons:
        lines.append(f"| {row['name']} | {row['prediction']:.9g} | {row['target_minus_prediction']:.9g} | {row['total_se']:.6g} | {row['z']:.6g} | {row['nominal_p_value']:.6g} |")
    lines += ["", "The two prediction comparisons share both the target and N400 anchor; they are not independent tests. The complete comparison covariance is saved, and no forced model ranking is reported.", "",
              f"Finite N400-to-N900 effective width: {gamma}.", "", contract["interpretation"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
