#!/usr/bin/env python3
"""Fixed block-fugacity prediction against saved width8 stationary response."""
from __future__ import annotations

import hashlib
import io
import json
import platform
import subprocess
import sys
import time
import types
from pathlib import Path

import numpy as np
import scipy
from scipy import linalg

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/p398_block_count_measure_contract.json"
OUTPUT = ROOT / "results/p398-block-count-measure"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def decode(value):
    array = np.asarray(value)
    return array[..., 0] + 1j * array[..., 1]


def encode(value):
    array = np.asarray(value)
    return np.stack((array.real, array.imag), axis=-1).tolist()


def maxabs(value):
    return float(np.max(np.abs(value)))


def main():
    start = time.perf_counter()
    output = OUTPUT / "latest.json"
    if output.exists():
        raise ValueError("Existing scientific result is preserved; reproduce in a separate checkout")
    contract = json.loads(CONTRACT.read_text())
    commit, directory = contract["source_commit"], contract["source_directory"]
    records = []

    def read(name):
        path = f"{directory}/{name}"
        data = subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)
        records.append({"commit": commit, "path": path, "sha256": sha(data)})
        return data

    # Reuse the pinned geometry definitions in isolated in-process modules.
    # Their command-line analysis entrypoints are never run.
    for name in ("frozen_model", "frozen_rate_model"):
        data = read(name + ".py")
        module = types.ModuleType(name)
        module.__file__ = f"{commit}/{directory}/{name}.py"
        sys.modules[name] = module
        exec(compile(data, module.__file__, "exec"), module.__dict__)
    model = sys.modules["frozen_rate_model"]
    old = json.loads(read("results/latest.json"))
    archived = np.load(io.BytesIO(read("previous_character_i_zero.npz")))
    states, parts, q, source, _, _, complement = model.construct(8)
    pi = np.asarray(old["stationary_probability"])
    dpi = np.asarray(old["stationary_derivative"])
    score = dpi / pi
    blocks = np.asarray([len(set(s)) for s in states], dtype=float)
    centered = blocks - pi @ blocks
    candidate_score = -2 * centered
    candidate_dpi = pi * candidate_score
    h, metric = archived["mass"], archived["metric"]
    if maxabs(source - archived["source"]) > 1e-12:
        raise ValueError("changed source/basis ordering")
    candidate_metric = (q.conj().T @ q.multiply(candidate_dpi[:, None])).toarray()

    def cov(m, propagator):
        return (source.conj().T @ m @ propagator @ source).conj()

    c0 = decode(old["integrals"]["C0"])
    candidate_dc0 = cov(candidate_metric, np.eye(len(h)))
    true_dc0 = decode(old["integrals"]["C0prime"])
    inverse = linalg.inv(h)
    integrated_u = decode(old["integrals"]["integrated_U"])
    candidate_integral_pi = linalg.solve(c0, cov(candidate_metric, inverse) - candidate_dc0 @ integrated_u)
    true_integral_pi = decode(old["integrals"]["integrated_Uprime_pi"])
    dynamic_integral = decode(old["integrals"]["integrated_Uprime_generator"])
    candidate_integral = candidate_integral_pi + dynamic_integral
    true_integral = decode(old["integrals"]["integrated_Uprime"])
    lag_rows = []
    for row in old["lag_results"]:
        lag = float(row["lag"])
        propagator = linalg.expm(-h * lag)
        candidate_pi_term = linalg.solve(c0, cov(candidate_metric, propagator)
                                         - candidate_dc0 @ decode(row["U"]))
        prediction = candidate_pi_term + decode(row["Uprime_generator"])
        lag_rows.append({"lag": lag, "candidate_Uprime": encode(prediction),
                         "candidate_Uprime_pi": encode(candidate_pi_term),
                         "true_Uprime": row["Uprime"], "true_Uprime_pi": row["Uprime_pi"],
                         "reused_true_generator_term": row["Uprime_generator"]})

    conditional_score = np.zeros_like(score)
    conditional_rows = []
    witness = None
    for b in sorted(set(blocks)):
        indices = np.flatnonzero(blocks == b)
        mass = float(pi[indices].sum())
        mean = float(pi[indices] @ score[indices] / mass)
        conditional_score[indices] = mean
        low, high = indices[np.argmin(score[indices])], indices[np.argmax(score[indices])]
        span = float(score[high] - score[low])
        conditional_rows.append({"blocks": int(b), "states": len(indices), "pi_mass": mass,
                                 "mean_score": mean, "min_score": float(score[low]),
                                 "max_score": float(score[high]), "score_span": span})
        if witness is None or span > witness["score_span"]:
            witness = {"blocks": int(b), "score_span": span,
                       "states": [{"state_index": int(i), "partition": list(states[i]),
                                   "pi": float(pi[i]), "true_score": float(score[i]),
                                   "candidate_score": float(candidate_score[i])} for i in (low, high)]}
    score_variance = float(pi @ score**2)
    candidate_error = float(pi @ (score - candidate_score)**2)
    within_block_variance = float(pi @ (score - conditional_score)**2)
    forward, perturbation = parts[0] + parts[1], parts[0] - parts[1]
    true_equation = forward @ dpi + perturbation @ pi
    candidate_equation = forward @ candidate_dpi + perturbation @ pi
    # These are deterministic discrepancies, not statistical hypothesis tests.
    result = {
        "schema": "matching-one.p398-block-count-measure.v1",
        "status": "computed_fixed_candidate_and_block_sufficiency_diagnostic",
        "contract": contract,
        "execution_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "inputs": records, "code_sha256": sha(Path(__file__).read_bytes()),
        "contract_sha256": sha(CONTRACT.read_bytes()),
        "state_count": len(states), "character_dimension": len(h),
        "stationary_measure": {
            "mean_blocks": float(pi @ blocks), "true_score_mean": float(pi @ score),
            "candidate_score_mean": float(pi @ candidate_score),
            "true_score_variance": score_variance,
            "candidate_weighted_MSE": candidate_error,
            "candidate_relative_RMS_error": float(np.sqrt(candidate_error / score_variance)),
            "within_equal_block_score_variance": within_block_variance,
            "within_equal_block_fraction": within_block_variance / score_variance,
            "true_stationary_equation_max_residual": maxabs(true_equation),
            "candidate_stationary_equation_max_residual": maxabs(candidate_equation),
            "candidate_K_odd_residual": maxabs(candidate_score[complement] + candidate_score),
            "groups": conditional_rows, "equal_block_witness": witness,
            "pi": pi.tolist(), "true_score": score.tolist(), "candidate_score": candidate_score.tolist(),
            "conditional_score_given_blocks": conditional_score.tolist()},
        "static": {"C0": encode(c0), "true_C0prime": encode(true_dc0),
                   "candidate_C0prime": encode(candidate_dc0)},
        "integrals": {"true_Uprime": encode(true_integral), "candidate_Uprime": encode(candidate_integral),
                      "true_pi_term": encode(true_integral_pi), "candidate_pi_term": encode(candidate_integral_pi),
                      "reused_true_generator_term": encode(dynamic_integral)},
        "lag_results": lag_rows,
        "environment": {"python": platform.python_version(), "platform": platform.platform(),
                        "numpy": np.__version__, "scipy": scipy.__version__},
        "elapsed_seconds": time.perf_counter() - start,
        "new_eta_points": 0, "stationary_solves": 0, "Frechet_recomputations": 0,
        "new_MC_samples": 0, "server_operations": 0, "test_suites": 0,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"seconds": result["elapsed_seconds"],
                      "measure": {k: v for k, v in result["stationary_measure"].items()
                                  if k not in ("groups", "pi", "true_score", "candidate_score", "conditional_score_given_blocks")},
                      "static_cross": {"true": float(true_dc0[0, 1].real), "candidate": float(candidate_dc0[0, 1].real)},
                      "integral_cross_true": [float(true_integral[i, j].real) for i, j in ((0, 1), (1, 0))],
                      "integral_cross_candidate": [float(candidate_integral[i, j].real) for i, j in ((0, 1), (1, 0))]}, indent=2))


if __name__ == "__main__":
    main()
