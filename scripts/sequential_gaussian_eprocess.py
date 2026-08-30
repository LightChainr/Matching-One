#!/usr/bin/env python3
"""Predeclared Gaussian likelihood-ratio e-processes.

The implementation is dependency-free and intentionally limited to frozen,
simple Gaussian hypotheses.  It is not a license to estimate means or
covariances from interim production data.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path


def _solve(matrix, vector):
    n = len(vector)
    a = [list(map(float, row)) + [float(vector[i])] for i, row in enumerate(matrix)]
    if len(a) != n or any(len(row) != n + 1 for row in a):
        raise ValueError("covariance must be square and match the mean dimension")
    for col in range(n):
        pivot = max(range(col, n), key=lambda row: abs(a[row][col]))
        if abs(a[pivot][col]) < 1e-15:
            raise ValueError("covariance is singular")
        a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        a[col] = [value / scale for value in a[col]]
        for row in range(n):
            if row == col:
                continue
            scale = a[row][col]
            a[row] = [x - scale * y for x, y in zip(a[row], a[col])]
    return [a[i][-1] for i in range(n)]


def _cholesky(matrix):
    n = len(matrix)
    if any(len(row) != n for row in matrix):
        raise ValueError("covariance must be square")
    out = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            value = float(matrix[i][j]) - sum(out[i][k] * out[j][k] for k in range(j))
            if i == j:
                if value <= 0.0:
                    raise ValueError("covariance must be positive definite")
                out[i][j] = math.sqrt(value)
            else:
                out[i][j] = value / out[j][j]
    return out


def log_lr_increment(observation, null_mean, alternative_mean, covariance):
    """Return log p_alt(x)/p_null(x) for equal-covariance Gaussians."""
    if not (len(observation) == len(null_mean) == len(alternative_mean)):
        raise ValueError("observation and means must have the same dimension")
    delta = [a - n for n, a in zip(null_mean, alternative_mean)]
    weight = _solve(covariance, delta)
    midpoint = [(n + a) / 2.0 for n, a in zip(null_mean, alternative_mean)]
    return sum(w * (x - m) for w, x, m in zip(weight, observation, midpoint))


def run_path(observations, null_mean, alternative_mean, covariance, alpha,
             minimum_batches, maximum_batches):
    """Apply symmetric e-value boundaries and return a frozen decision record."""
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0, 1)")
    if not 1 <= minimum_batches <= maximum_batches:
        raise ValueError("invalid batch limits")
    boundary = math.log(1.0 / alpha)
    log_e = 0.0
    trace = []
    for batch, observation in enumerate(observations, 1):
        if batch > maximum_batches:
            break
        log_e += log_lr_increment(observation, null_mean, alternative_mean, covariance)
        trace.append({"batch": batch, "log_e_alternative_vs_null": log_e})
        if batch >= minimum_batches and log_e >= boundary:
            return {"decision": "alternative", "batches": batch, "log_e": log_e, "trace": trace}
        if batch >= minimum_batches and log_e <= -boundary:
            return {"decision": "null", "batches": batch, "log_e": log_e, "trace": trace}
    return {"decision": "inconclusive", "batches": len(trace), "log_e": log_e, "trace": trace}


def _draw(mean, chol, rng):
    z = [rng.gauss(0.0, 1.0) for _ in mean]
    return [mean[i] + sum(chol[i][j] * z[j] for j in range(i + 1)) for i in range(len(mean))]


def calibrate(config, simulations, seed):
    rng = random.Random(seed)
    models = config["models"]
    covariance = config["batch_covariance"]
    chol = _cholesky(covariance)
    alpha = float(config["alpha"])
    minimum = int(config["minimum_batches"])
    maximum = int(config["maximum_batches"])
    result = {
        "schema_version": 1,
        "seed": seed,
        "simulations_per_model_test": simulations,
        "guarantee_scope": config["guarantee_scope"],
        "tests": {},
    }
    for test in config["ordered_tests"]:
        null_name, alternative_name = test["null"], test["alternative"]
        null_mean, alt_mean = models[null_name], models[alternative_name]
        delta = [a - n for n, a in zip(null_mean, alt_mean)]
        information = sum(d * w for d, w in zip(delta, _solve(covariance, delta)))
        rows = {}
        for generating_name, mean in models.items():
            counts = {"null": 0, "alternative": 0, "inconclusive": 0}
            total_batches = 0
            for _ in range(simulations):
                observations = (_draw(mean, chol, rng) for _ in range(maximum))
                path = run_path(observations, null_mean, alt_mean, covariance,
                                alpha, minimum, maximum)
                counts[path["decision"]] += 1
                total_batches += path["batches"]
            rows[generating_name] = {
                "decision_rates": {key: value / simulations for key, value in counts.items()},
                "mean_batches": total_batches / simulations,
                "compute_fraction_vs_fixed": total_batches / (simulations * maximum),
            }
        result["tests"][test["name"]] = {
            "null": null_name,
            "alternative": alternative_name,
            "mahalanobis_information_per_batch": information,
            "expected_log_e_drift_at_endpoint": information / 2.0,
            "ville_error_bound_each_direction": alpha,
            "models": rows,
        }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path)
    parser.add_argument("--simulations", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=126)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    result = calibrate(config, args.simulations, args.seed)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
