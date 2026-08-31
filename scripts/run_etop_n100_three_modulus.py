#!/usr/bin/env python3
"""Acquire and score the frozen three-shape N100 experiment locally."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.linalg import block_diag
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import binom, chi2

from etop_modulus_survivors import FIELDS, interpolation_weight

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "experiments/etop_n100_three_modulus_20260831.json"
OUT = ROOT / "results/etop-n100-three-modulus"
FREEZE = "4c1ec50"


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def launch(contract):
    (OUT/"raw").mkdir(parents=True, exist_ok=True)
    (OUT/"logs").mkdir(exist_ok=True)
    build = ROOT/"build/etop-n100"
    build.mkdir(parents=True, exist_ok=True)
    binary = build/"threshold_rank_integer_period_mc"
    source = ROOT/"src/threshold_rank_integer_period_mc.cpp"
    subprocess.run(["clang++", "-O3", "-std=c++17", str(source), "-o", str(binary)], check=True)
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()

    def one(shape):
        prefix = OUT/"raw"/shape["name"]
        if prefix.with_suffix(".metadata.json").exists():
            raise RuntimeError("A completed target already exists; do not overwrite it")
        cmd = [str(binary), "--samples", str(contract["samples_per_shape_pair"]),
               "--batches", str(contract["batches"]), "--seed", str(contract["seed"]),
               "--replica-offset", str(contract["replica_offset"]), "--threads", "1",
               "--first-matrix", *map(str, shape["first"]),
               "--second-matrix", *map(str, shape["second"]), "--git-commit", head,
               "--output-prefix", str(prefix)]
        started = datetime.now(timezone.utc).isoformat()
        tic = time.monotonic()
        with (OUT/"logs"/(shape["name"]+".log")).open("w") as log:
            result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        receipt = {"shape": shape, "command": cmd, "started_utc": started,
                   "elapsed_seconds": time.monotonic()-tic, "exit_code": result.returncode,
                   "engine_sha256": digest(source), "binary_sha256": digest(binary),
                   "runner_commit": head, "prediction_freeze_commit": FREEZE}
        dump(OUT/"logs"/(shape["name"]+".receipt.json"), receipt)
        result.check_returncode()
        print(f"Completed {shape['name']} in {receipt['elapsed_seconds']:.2f}s", flush=True)
        return receipt
    with ThreadPoolExecutor(max_workers=3) as executor:
        list(executor.map(one, contract["shapes"]))


def read_batch_vectors(contract, shape):
    n, batches = contract["area"], contract["batches"]
    per_batch = contract["samples_per_shape_pair"]//batches
    hist = np.zeros((batches, 2, 2, n+1))
    prefix = OUT/"raw"/shape["name"]
    with prefix.with_suffix(".hist.csv").open() as stream:
        for row in csv.DictReader(stream):
            b, o, k = int(row["batch"]), int(row["orientation"] == "second"), int(row["k"])
            kind = int(row["kind"] == "plus")
            hist[b, o, kind, k] += int(row["count"])
    if not np.all(hist.sum(axis=-1) == per_batch):
        raise ValueError("Incomplete threshold histograms")
    tail = binom.sf(np.arange(n+1)-1, n, contract["fixed_p"])
    f = hist@tail/per_batch
    values = np.zeros((batches, 2, 4))
    values[:, :, 0] = f[:, :, 0]+f[:, :, 1]-1
    values[:, :, 1] = 1+f[:, :, 1]-f[:, :, 0]
    with prefix.with_suffix(".moments.csv").open() as stream:
        for row in csv.DictReader(stream):
            b, o = int(row["batch"]), int(row["orientation"] == "second")
            km, kp = int(row["sum_kminus"]), int(row["sum_kplus"])
            values[b, o, 2] = (km+kp)/(2*(n+1)*per_batch)
            values[b, o, 3] = (kp-km)/((n+1)*per_batch)
    return (values[:, 0]-values[:, 1])/float(Fraction(shape["delta_cos4"]))


def wald(mean, covariance):
    statistic = float(mean@np.linalg.solve(covariance, mean))
    return {"estimate": mean.tolist(), "covariance": covariance.tolist(),
            "standard_error": np.sqrt(np.diag(covariance)).tolist(),
            "z": (mean/np.sqrt(np.diag(covariance))).tolist(),
            "chi_square": statistic, "df": len(mean), "p_value": float(chi2.sf(statistic, len(mean)))}


def weight_score(weights, mean, covariance):
    projection = np.kron(np.array(weights).reshape(1, 3), np.eye(4))
    return wald(projection@mean, projection@covariance@projection.T)


def same_area_models(mean, covariance):
    models = {}
    for model in ("affine_E4", "affine_height_E4", "affine_height_squared", "affine_log_height"):
        t4, ts = interpolation_weight(model, 4j), interpolation_weight(model, .5+1j)
        secant = (ts-1)/(t4-1)
        weights = [secant-1, -secant, 1.]
        models[model] = {"weights": weights, "secant_ratio": secant,
                         **weight_score(weights, mean, covariance)}
    def objective(theta):
        c, s = np.cos(theta), np.sin(theta)
        return weight_score([s-c, -s, c], mean, covariance)["chi_square"]
    grid = np.linspace(-np.pi/2, np.pi/2, 2049)
    q = np.array([objective(t) for t in grid])
    candidates = [(q[0], grid[0])]
    for i in range(1, len(grid)-1):
        if q[i] <= q[i-1] and q[i] <= q[i+1]:
            fit = minimize_scalar(objective, bounds=(grid[i-1], grid[i+1]), method="bounded")
            candidates.append((fit.fun, fit.x))
    best, theta = min(candidates)
    models["free_common_secant"] = {"secant_ratio": float(np.tan(theta)),
        "chi_square": float(best), "df": 3, "p_value": float(chi2.sf(best, 3)),
        "scope": "Any one shared scalar shape coordinate with independent four-vector offset/slope; retrospective flexible comparison"}
    flat = np.kron(np.array([[-1., 1., 0.], [-1., 0., 1.]]), np.eye(4))
    models["no_shape_response"] = wald(flat@mean, flat@covariance@flat.T)
    return models


def source_transfer(mean, covariance):
    old = json.loads(subprocess.check_output(["git", "show", f"{FREEZE}:results/etop-modulus-survivors/latest.json"], cwd=ROOT))
    rows = old["input"]["rows"]
    x, y = (np.array(rows[k]["estimate"]) for k in ("tau_i", "tau_2i"))
    sx, sy = (np.array(rows[k]["covariance"]) for k in ("tau_i", "tau_2i"))
    observed = np.r_[x, y, mean]
    cov = block_diag(sx, sy, covariance)
    whitening = np.linalg.inv(np.linalg.cholesky(cov))
    whitened = whitening@observed
    result = {}
    for model in ("affine_E4", "affine_height_E4", "affine_height_squared", "affine_log_height"):
        t4, ts = interpolation_weight(model, 4j), interpolation_weight(model, .5+1j)
        def profile(gains, detail=False):
            d = np.diag(gains)
            matrix = np.zeros((20, 8))
            matrix[:8] = np.eye(8)
            matrix[8:12, 4:] = d
            matrix[12:16, :4], matrix[12:16, 4:] = (1-t4)*d, t4*d
            matrix[16:20, :4], matrix[16:20, 4:] = (1-ts)*d, ts*d
            matrix = whitening@matrix
            latent = np.linalg.lstsq(matrix, whitened, rcond=None)[0]
            residual = whitened-matrix@latent
            return (float(residual@residual), latent) if detail else float(residual@residual)
        starts = [np.full(4, q) for q in (.1, .3, .6, 1.)]
        starts.append(np.clip(mean[:4]/y, -3, 3))
        fits = [minimize(profile, start, method="Nelder-Mead",
                        options={"maxiter": 6000, "xatol": 1e-9, "fatol": 1e-9}) for start in starts]
        fit = min(fits, key=lambda r: r.fun)
        result[model] = {"chi_square": float(fit.fun), "df": 8,
                        "p_value": float(chi2.sf(fit.fun, 8)), "area_gains": fit.x.tolist(),
                        "optimizer_success": bool(fit.success),
                        "latent_N50_i_2i": profile(fit.x, True)[1].tolist()}
    return result


def score(contract):
    blocks = [read_batch_vectors(contract, shape) for shape in contract["shapes"]]
    batches = np.column_stack(blocks)
    mean = batches.mean(axis=0)
    covariance = np.cov(batches, rowvar=False, ddof=1)/len(batches)
    result = {"contract": contract, "prediction_freeze_commit": FREEZE,
        "shape_order": [s["name"] for s in contract["shapes"]], "field_order": FIELDS,
        "mean": mean.tolist(), "covariance": covariance.tolist(),
        "batch_vectors": batches.tolist(), "same_area_models": same_area_models(mean, covariance),
        "source_shape_transfer": source_transfer(mean, covariance),
        "raw_sha256": {str(p.relative_to(ROOT)): digest(p) for p in sorted((OUT/"raw").iterdir())},
        "inference": "Gaussian covariance/profile comparisons; same-area primary is a fixed linear null. Models/readouts reuse one new common-random block; scores are not independent evidence."}
    dump(OUT/"score.json", result)
    lines = ["# N100: a new three-modulus response experiment", "",
             "Each pair has two million new shared-counter permutations in 200 aligned batches. All three shapes were frozen before acquisition; the full 12x12 covariance is retained.", "",
             "| shape | A_top | E_top | C | W |", "|---|---:|---:|---:|---:|"]
    for i, shape in enumerate(contract["shapes"]):
        fields = [f"{mean[4*i+j]:.8g} +/- {np.sqrt(covariance[4*i+j,4*i+j]):.3g}" for j in range(4)]
        lines.append("| "+shape["name"]+" | "+" | ".join(fields)+" |")
    lines += ["", "## Same-area shape comparison: no N50 calibration or area exponent", "",
              "| model | chi-square / df | p |", "|---|---:|---:|"]
    for name, row in result["same_area_models"].items():
        lines.append(f"| {name} | {row['chi_square']:.7g} / {row['df']} | {row['p_value']:.6g} |")
    lines += ["", "The first row is the frozen primary. Other fixed rows are declared comparators; free-common-secant is an exploratory one-parameter relaxation. None is a continuum field count or E4 identity proof.", "",
              "## Source-informed N50-to-N100 shape transfer", "",
              "Independent gains are profiled for A/E/C/W, with both uncertain source vectors and all same-stream target covariance. This tests the additional cross-area shape-separability hypothesis without dividing by the weak E denominator.", "",
              "| model | chi-square / df | p |", "|---|---:|---:|"]
    for name, row in result["source_shape_transfer"].items():
        lines.append(f"| {name} | {row['chi_square']:.7g} / {row['df']} | {row['p_value']:.6g} |")
    lines += ["", result["inference"], "", f"Prediction freeze: {FREEZE}. Exact geometry: b9e4ea1; three-modulus null: 964d770. No old source p-values are added to the new scores.", ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    contract = json.loads(CONTRACT.read_text())
    if args.run:
        launch(contract)
    score(contract)
