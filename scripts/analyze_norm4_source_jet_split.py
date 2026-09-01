#!/usr/bin/env python3
"""Split one archived cluster-source tangent into topology and within-sector jets.

The raw source is primary. W=L[R] is the unique compensation defined in
notes/p40-thermal-clock-source-quotient.md, not a new independent evidence row.
This consumes complete Newman--Ziff K profiles, not fixed-p importance samples.
"""
from __future__ import annotations

import json
import math
import platform
import subprocess
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy

import analyze_norm4_source_thermal as parent


ROOT = parent.ROOT
NS = parent.NS
DIRECTIONS = ("first", "second")
STATES = np.array([-1.0, 0.0, 1.0])
FIELDS = ("raw_bulk", "within_sector_bulk", "topology_only_bulk")
DEFINITION = ROOT / "notes/p40-thermal-clock-source-quotient.md"
DESTINATION = parent.OUTPUT / "jet-split.json"


def conditional_table(sums, samples, p, n):
    """Recover conditional (S,K) moments without requiring S squared.

    parent.read_raw has already divided s, qs and es by N. Each K has
    exactly `samples` prefixes, so its rank-sector counts partition samples.
    Binomial integration includes its normalizing denominator, as in parent.
    """
    q, e, s, qs, es = sums.T
    counts = np.column_stack(((e - q) / 2, samples - e, (e + q) / 2))
    source_sums = np.column_stack(((es - qs) / 2, s - es, (es + qs) / 2))
    k = np.arange(n + 1, dtype=float)[:, None]
    joint = np.concatenate((counts, counts * k, source_sums), axis=1)
    integrated = parent.binomial_moments(joint, samples, p, n)[0]
    probabilities = integrated[:3]
    if not np.isfinite(integrated).all() or np.any(probabilities <= 0):
        raise ValueError(f"N{n}: conditional source needs all three observed sectors")
    mean_k = integrated[3:6] / probabilities
    mean_s = integrated[6:9] / probabilities

    # Center on each sector before integration to avoid subtracting large K^2
    # moments. This uses source first moments only, never S^2.
    dk = k - mean_k[None, :]
    centered = np.concatenate(
        (counts * dk * dk,
         (source_sums - counts * mean_s[None, :]) * dk), axis=1)
    centered_moments = parent.binomial_moments(centered, samples, p, n)[0]
    var_k = centered_moments[:3] / probabilities
    cov_sk = centered_moments[3:6] / probabilities
    if not np.isfinite(centered_moments).all() or np.any(var_k < 0):
        raise ValueError(f"N{n}: invalid within-sector conditional moments")
    return {
        "probabilities": probabilities,
        "mean_K": mean_k,
        "mean_S_density": mean_s,
        "var_K": var_k,
        "cov_S_density_K": cov_sk,
    }


def split_at_root(sums, samples, n, delta, bracket):
    raw, raw_diagnostic = parent.at_root(sums, samples, n, delta, bracket)
    p0, d = raw["p0"], raw["D"]
    tables = [conditional_table(sums[g], samples, p0, n) for g in range(2)]
    numerator = math.fsum(float(t["probabilities"] @ t["cov_S_density_K"])
                         for t in tables)
    denominator = math.fsum(float(t["probabilities"] @ t["var_K"])
                           for t in tables)
    if denominator < 0 or not math.isfinite(denominator):
        raise ValueError(f"N{n}: invalid common within-sector clock variance")
    # When K is exactly determined by q, this clock direction is redundant.
    beta = numerator / denominator if denominator > 0 else 0.0
    if denominator == 0 and numerator != 0:
        raise ValueError(f"N{n}: zero within-sector clock variance with nonzero covariance")

    jq_p, je_p, diagnostics = [], [], {}
    for g, table in zip(DIRECTIONS, tables):
        probabilities = table["probabilities"]
        gamma = table["cov_S_density_K"] - beta * table["var_K"]
        row = raw_diagnostic["direction"][g]
        scale = p0 * (1 - p0)
        jq_p.append(float(probabilities @ ((STATES - row["q"]) * gamma) / scale))
        je_p.append(float(probabilities @ ((STATES**2 - row["E"]) * gamma) / scale))
        diagnostics[g] = {
            "sector_order": [-1, 0, 1],
            **{key: value.tolist() for key, value in table.items()},
            "gamma_density": gamma.tolist(),
            "topology_counterterm_density":
                (table["mean_S_density"] - beta * table["mean_K"]).tolist(),
            "Jq_prime_R_density": jq_p[-1],
            "JE_prime_R_density": je_p[-1],
            "cov_R_K_density": float(probabilities @ gamma),
        }

    # R_g=S_g-E[S_g|q]-beta*(K-E[K|q]); the conditional functions and beta
    # above are FROZEN while differentiating p. Jq(R)=JE(R)=0 at this root,
    # so only the direct and source-slope pieces remain; root motion is zero.
    p4_ep = (raw_diagnostic["direction"]["first"]["E_p"]
             - raw_diagnostic["direction"]["second"]["E_p"]) / delta
    prefactor = n ** (13 / 8) / 2
    direct_density = prefactor * (je_p[0] - je_p[1]) / (delta * d)
    slope_density = -prefactor * p4_ep * math.fsum(jq_p) / (2 * d * d)
    within_bulk = n * math.fsum((direct_density, slope_density))
    raw_bulk = raw["Udot_fugacity"]
    values = dict(zip(FIELDS, (raw_bulk, within_bulk, raw_bulk - within_bulk)))
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError(f"N{n}: nonfinite source split")
    diagnostic = {
        "p0": p0, "D": d, "permutations": samples,
        "raw_U_state": raw,
        "b_star_density": beta, "b_star_bulk": n * beta,
        "pooled_clock_projection_numerator": numerator,
        "pooled_clock_projection_denominator": denominator,
        "clock_direction_redundant": denominator == 0,
        "directions": diagnostics,
        "within_sector_pieces_bulk": {
            "direct": n * direct_density,
            "slope_source": n * slope_density,
            "rootmotion": 0.0, "slope_root": 0.0,
        },
        "pooled_cov_R_K_density": math.fsum(
            diagnostics[g]["cov_R_K_density"] for g in DIRECTIONS) / 2,
        "source_identity": "S_g=R_g+b_star*K+f_g(q), f_g=E[S_g|q]-b_star*E[K|q]",
        "response_identity": "L[S]=L[R]+L[f_g(q)]; the common b_star*K has zero response",
    }
    return values, diagnostic


def vectorize(by_n):
    return {f"N{n}.{field}": by_n[n][field] for n in NS for field in FIELDS}


def main():
    started = time.perf_counter()
    if DESTINATION.exists():
        raise ValueError("saved jet split exists; do not overwrite or repeat")
    manifest = json.loads(parent.MANIFEST.read_text())
    hypotheses = json.loads(parent.HYPOTHESES.read_text())
    runs = {run["N"]: run for run in manifest["runs"]}
    bracket = hypotheses["root_bracket"]
    selection = manifest["selection"]
    samples = int(selection["selected_samples_per_N"])
    batches = int(selection["selected_batches_per_N"])
    batch_samples = int(selection["selected_samples_per_batch"])
    if batches != 100 or samples != batches * batch_samples:
        raise ValueError("source split requires the parent 100-batch sampling contract")
    raw = {n: parent.read_raw(parent.OUTPUT / "raw" / f"n{n}.csv", n) for n in NS}
    totals = {n: raw[n].sum(axis=0) for n in NS}
    deltas = {n: float(Fraction(runs[n]["delta_cos4"])) for n in NS}

    central_by_n, diagnostics = {}, {}
    for n in NS:
        central_by_n[n], diagnostics[n] = split_at_root(
            totals[n], samples, n, deltas[n], bracket)
    central_map = vectorize(central_by_n)
    labels = list(central_map)
    central = np.array(list(central_map.values()))
    covariance = np.zeros((len(labels), len(labels)))
    groups = {}
    for n in NS:
        groups.setdefault(runs[n]["dependency_group"], []).append(n)
    group_results = {}
    for group, ns in groups.items():
        vectors = []
        roots, betas = {n: [] for n in ns}, {n: [] for n in ns}
        for batch in range(batches):
            omitted = dict(central_by_n)
            for n in ns:
                omitted[n], diagnostic = split_at_root(
                    totals[n] - raw[n][batch], samples - batch_samples,
                    n, deltas[n], bracket)
                roots[n].append(diagnostic["p0"])
                betas[n].append(diagnostic["b_star_density"])
            vectors.append(list(vectorize(omitted).values()))
        vectors = np.asarray(vectors)
        deviations = vectors - vectors.mean(axis=0)
        contribution = (batches - 1) / batches * deviations.T @ deviations
        covariance += contribution
        group_results[group] = {
            "Ns": ns, "delete_one_batch_ids": list(range(batches)),
            "delete_one_vectors": vectors.tolist(),
            "covariance_contribution": contribution.tolist(),
            "delete_one_roots": roots, "delete_one_b_star_density": betas,
            "operation": "Omit the same permutation batch in every N of this group; other groups remain central. Refind each root and refit conditional moments/b_star, then freeze those coefficients for thermal derivatives.",
        }
    errors = np.sqrt(np.maximum(0, np.diag(covariance)))
    estimates = {
        label: {"value": float(value), "se": float(error),
                "z": float(value / error) if error > 0 else None}
        for label, value, error in zip(labels, central, errors)
    }
    code_paths = (Path(__file__), Path(parent.__file__),
                  ROOT / "scripts/analyze_p40_source_thermal.py",
                  parent.MANIFEST, parent.HYPOTHESES, DEFINITION)
    result = {
        "schema": "matching-one.norm4-source-jet-split.v1",
        "status": "computed_existing_production_source_decomposition",
        "execution_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "labels": labels, "estimates": estimates,
        "covariance": covariance.tolist(), "by_N": diagnostics,
        "dependency_groups": group_results,
        "replay_manifest": manifest,
        "definition": {
            "primary_source": "S=(CB+CW)/N; no replacement of the primary raw source",
            "bulk_coordinate": "All reported v,W,v-W are derivatives for the common fugacity exp(t*(CB+CW)); multiply density-source responses by N",
            "within_source": "R_g=S_g-E_g[S_g|q]-b_star*(K-E_g[K|q])",
            "b_star": "sum_g E_g[Cov(S,K|q)] / sum_g E_g[Var(K|q)], equal geometry weights and one common coefficient",
            "thermal_derivative": "Jf_prime(R)=sum_r P(q=r)*(f(r)-mean(f))*[Cov(S,K|r)-b_star*Var(K|r)]/[p0*(1-p0)]",
            "within_response": "W=L[R]; Jq(R)=JE(R)=0 at the baseline, so only direct and source-slope contributions remain",
            "topology_only_remainder": "v-W=L[f_g(q)], f_g(q)=E_g[S_g|q]-b_star*E_g[K|q]; the common b_star*K is annihilated by L",
            "coefficient_convention": "Refit at each full/leave-one-out root, then hold b_star and the conditional functions fixed when differentiating p and lambda",
        },
        "interpretation": {
            "question": "Does the same raw source act on the intrinsic angular thermal jet through three-state reweighting, within-sector cluster/occupancy coupling, or opposing contributions?",
            "evidence": "One exact decomposition of the same source and archived samples, not independent evidence votes; do not infer dominance from unstable ratios or add component significances",
            "within_scope": "A resolved W detects this thermal-jet coupling beyond one common K clock and arbitrary fixed per-geometry functions of q; it does not escape all functions of (q,K), count fields or identify a continuum energy operator",
            "topology_scope": "The remainder is a source projection at the root, with geometry-dependent three-state coefficients; it is not a universal microscopic field or an independently sampled source",
            "covariance": "The complete 18-coordinate covariance is redundant because raw=within+topology. It is saved without inversion or an omnibus component test",
            "missing_statistic": "No S^2 mark is available or used; no residual-source variance, efficiency, or variance-normalized attribution is reported",
        },
        "estimator": "Reuse the parent complete-prefix Binomial integration and original U/root/slope definition. K prefixes of one permutation are correlated readouts, not independent samples. New compensated marks are algebraic projections of the same stored q,E,S,qS,ES profiles.",
        "uncertainty": "Sum the same three group-specific aligned-delete-one covariance contributions as the parent; 100 full output vectors per group are retained",
        "finite_subset_scope": "100000 original production permutations per N; these errors do not inherit the full 1.9B/1B production precision",
        "inputs": [{"path": str((parent.OUTPUT / "raw" / f"n{n}.csv").relative_to(ROOT)),
                    "sha256": parent.sha(parent.OUTPUT / "raw" / f"n{n}.csv")} for n in NS],
        "code": [{"path": str(path.relative_to(ROOT)), "sha256": parent.sha(path)}
                 for path in code_paths],
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "machine": platform.machine()},
        "elapsed_seconds": time.perf_counter() - started,
        "new_samples": 0, "new_configuration_replays": 0,
        "server_actions": 0, "test_suites": [],
    }
    with DESTINATION.open("x") as handle:
        handle.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(DESTINATION), "elapsed_seconds": result["elapsed_seconds"],
                      "estimates": estimates}, indent=2))


if __name__ == "__main__":
    main()
