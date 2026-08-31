#!/usr/bin/env python3
"""Occupancy/rank source response decomposition on saved norm-4 configurations."""
from __future__ import annotations

import hashlib
import json
import math
import platform
import subprocess
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy

from analyze_norm4_source_endpoint_1m import load_profile
from analyze_norm4_source_thermal import binomial_moments, direction_values

ROOT = Path(__file__).resolve().parents[1]
NS = (65, 85, 130, 170, 260, 340)
CONTRACT = ROOT / "analysis/norm4_global_source_projection_contract.json"
OUTPUT = ROOT / "results/norm4-global-source-projection"
FIELDS = ("v", "rootdot", "rank1_rootdot", "v_direct", "v_rootmotion",
          "v_slope_source", "v_slope_root")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_response(sums, samples, n, p0):
    rows = [direction_values(binomial_moments(sums[g], samples, p0, n)) for g in range(2)]
    d = np.mean([r["q_p"] for r in rows])
    delta = float(Fraction("1152/845" if n in (65, 130, 260) else "2304/1445"))
    b = (rows[0]["E_p"] - rows[1]["E_p"]) / delta
    h = (rows[0]["E_pp"] - rows[1]["E_pp"]) / delta
    t = np.mean([r["q_pp"] for r in rows])
    jq = np.mean([r["Jq"] for r in rows])
    jqp = np.mean([r["Jq_p"] for r in rows])
    je = np.mean([r["JE"] for r in rows])
    jep = (rows[0]["JE_p"] - rows[1]["JE_p"]) / delta
    ep = np.mean([r["E_p"] for r in rows])
    rootdot = -jq / d
    scale = n ** (13 / 8) / 2
    pieces = (scale * jep / d, scale * h * rootdot / d,
              -scale * b * jqp / d**2, -scale * b * t * rootdot / d**2)
    values = (math.fsum(pieces), rootdot, -je - ep * rootdot, *pieces)
    return dict(zip(FIELDS, map(float, values))), float(scale * b / d)


def point(sums, samples, n, p0):
    # load_profile uses density units; the entire analysis here uses bulk units.
    bulk = sums.copy()
    bulk[..., 2:] *= n
    occupancy = bulk.copy()
    occupancy[..., 3] = bulk[..., 0] * bulk[..., 2] / samples
    occupancy[..., 4] = bulk[..., 1] * bulk[..., 2] / samples
    rank = bulk.copy()
    rank[..., 2] = 0
    rank[..., 3:] -= occupancy[..., 3:]
    out = {"p0": float(p0)}
    for component, profile in (("total", bulk), ("occupancy", occupancy), ("rank", rank)):
        response, u = source_response(profile, samples, n, p0)
        out.update({f"{component}.{key}": value for key, value in response.items()})
        out["U"] = u
    return out


def vectorize(points):
    values = {f"N{n}.{key}": value for n in NS for key, value in points[n].items()}
    for start, middle, end in ((65, 130, 260), (85, 170, 340)):
        for model, coeff in (("q2", (1, -3, 2)), ("Jordan", (1, -2, 1))):
            for component in ("total", "occupancy", "rank"):
                values[f"{model}.{start}.{component}.v"] = math.fsum(
                    c * points[n][f"{component}.v"] for n, c in zip((start, middle, end), coeff))
    return values


def addback_error(points):
    return max(abs(row[f"total.{field}"] - row[f"occupancy.{field}"] - row[f"rank.{field}"])
               for row in points.values() for field in FIELDS)


def main():
    started = time.perf_counter()
    destination = OUTPUT / "latest.json"
    if destination.exists():
        raise ValueError("Saved result exists; use a separate checkout for reproduction")
    contract = json.loads(CONTRACT.read_text())
    source_path = ROOT / contract["source_result"]
    if digest(source_path) != contract["source_result_sha256"]:
        raise ValueError("source result changed from the declared input")
    source = json.loads(source_path.read_text())
    source_index = {label: i for i, label in enumerate(source["labels"])}
    geometry_path = ROOT / "analysis/p40_source_thermal_chain_candidates.json"
    runs = {run["N"]: run for run in json.loads(geometry_path.read_text())["runs"]}
    profiles, samples, inputs = {}, {}, []
    for n in NS:
        base_path = ROOT / f"results/norm4-source-thermal/raw/n{n}.csv"
        profiles[n] = load_profile(base_path, n, 1000, runs[n])
        inputs.append({"path": str(base_path.relative_to(ROOT)), "sha256": digest(base_path)})
        samples[n] = 100000
        if n in (260, 340):
            increment = ROOT / f"results/norm4-source-endpoint-1m/increment/raw/n{n}.csv"
            profiles[n] += load_profile(increment, n, 9000, runs[n])
            inputs.append({"path": str(increment.relative_to(ROOT)), "sha256": digest(increment)})
            samples[n] = 1000000
    totals = {n: profiles[n].sum(axis=0) for n in NS}
    points = {n: point(totals[n], samples[n], n, source["by_N"][str(n)]["points"]["p0"])
              for n in NS}
    central_map = vectorize(points)
    labels = list(central_map)
    central = np.array(list(central_map.values()))
    covariance = np.zeros((len(labels), len(labels)))
    groups = {}
    addback = addback_error(points)
    saved_difference = {field: 0.0 for field in ("v", "rootdot", "rank1_rootdot", "U")}
    old_fields = {"v": "Udot_fugacity", "rootdot": "rootdot_fugacity",
                  "rank1_rootdot": "root_comoving_rank1_fugacity", "U": "U"}

    def compare_saved(rows, saved):
        for n, row in rows.items():
            for field, old in old_fields.items():
                actual = row["U" if field == "U" else f"total.{field}"]
                saved_difference[field] = max(saved_difference[field], abs(actual - saved(n, old)))

    compare_saved(points, lambda n, old: source["by_N"][str(n)]["points"][old])
    for name, previous in source["covariance_contributions"].items():
        if not name.startswith("source:"):
            continue
        if previous["delete_one_batch_ids"] != list(range(100)):
            raise ValueError("unaligned source omissions")
        saved_vectors = np.asarray(previous["delete_one_vectors"])
        vectors = []
        for batch in range(100):
            changed = dict(points)
            for n in previous["Ns"]:
                p0 = saved_vectors[batch, source_index[f"N{n}.p0"]]
                changed[n] = point(totals[n] - profiles[n][batch], samples[n] * .99, n, p0)
            addback = max(addback, addback_error(changed))
            compare_saved({n: changed[n] for n in previous["Ns"]},
                          lambda n, old: saved_vectors[batch, source_index[f"N{n}.{old}"]])
            vectors.append(list(vectorize(changed).values()))
        vectors = np.asarray(vectors)
        deviation = vectors - vectors.mean(axis=0)
        contribution = .99 * deviation.T @ deviation
        covariance += contribution
        groups[name] = {"Ns": previous["Ns"], "delete_one_batch_ids": list(range(100)),
                        "batch_counts": previous["batch_counts"],
                        "delete_one_vectors": vectors.tolist(), "covariance": contribution.tolist(),
                        "operation": "same original paired omission; recompute m_g(K) and rank source on retained samples"}
    errors = np.sqrt(np.maximum(0, covariance.diagonal()))
    estimates = {label: {"value": float(value), "se": float(error),
                         "z": float(value / error) if error else None}
                 for label, value, error in zip(labels, central, errors)}
    if not np.isfinite(central).all() or not np.isfinite(covariance).all():
        raise ValueError("nonfinite output")
    if max(saved_difference.values()) > 1e-6 or addback > 1e-6:
        raise ValueError(f"changed response semantics: {saved_difference}, addback={addback}")
    result = {
        "schema": "matching-one.norm4-global-source-projection.v1",
        "status": "computed_existing_data_source_decomposition",
        "execution_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "contract": contract, "labels": labels, "estimates": estimates,
        "covariance": covariance.tolist(), "covariance_contributions": groups,
        "samples_per_N": samples, "inputs": inputs,
        "source_result_sha256": digest(source_path), "geometry_sha256": digest(geometry_path),
        "code_sha256": digest(Path(__file__)), "contract_sha256": digest(CONTRACT),
        "identities": {"max_total_minus_components_central_and_LOO": addback,
                       "max_difference_from_saved_source_central_and_LOO": saved_difference,
                       "spatial_residual_global_response": "exactly_zero_by_conditional_expectation_not_a_test"},
        "environment": {"python": platform.python_version(), "platform": platform.platform(),
                        "numpy": np.__version__, "scipy": scipy.__version__},
        "elapsed_seconds": time.perf_counter() - started,
        "new_random_samples": 0, "configuration_replays": 0, "root_finders": 0,
        "server_operations": 0, "test_suites": 0,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(destination.relative_to(ROOT)), "seconds": result["elapsed_seconds"],
                      "identities": result["identities"], "covariance_dimension": len(labels)}))
    for n in NS:
        print(n, {component: {field: estimates[f"N{n}.{component}.{field}"]
                             for field in ("v", "rootdot", "rank1_rootdot")}
                  for component in ("total", "occupancy", "rank")})


if __name__ == "__main__":
    main()
