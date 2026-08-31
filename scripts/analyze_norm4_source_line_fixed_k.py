#!/usr/bin/env python3
"""Decompose existing conditional-line source covariance at saved roots only."""
from __future__ import annotations

import json
import math
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import scipy

import analyze_norm4_source_line as primary

ROOT = primary.ROOT
NS = primary.NS
CONTRACT = ROOT / "analysis/norm4_source_line_fixed_k_contract.json"
COMPONENTS = primary.COMPONENTS


def decompose(sums, samples, n, p0, rootdot, primary_values, tail):
    """All-K plug-in identity and a separately labelled Bessel sensitivity."""
    k = np.arange(n + 1, dtype=float)
    logw = (primary.old.gammaln(n + 1) - primary.old.gammaln(k + 1)
            - primary.old.gammaln(n - k + 1)
            + k * math.log(p0) + (n - k) * math.log1p(-p0))
    w = np.exp(logw)
    w /= w.sum()
    cdf = np.cumsum(w)
    lower = int(np.searchsorted(cdf, tail))
    upper = min(n, int(np.searchsorted(cdf, 1 - tail)))
    diagnostic_window = (k >= lower) & (k <= upper)
    values = {"p0": float(p0), "rootdot_s": float(rootdot)}
    diagnostics = {}
    for g, direction in enumerate(("first", "second")):
        a, c = sums[g, :, 0], sums[g, :, 1]
        b = sums[g, :, 2] + 1j * sums[g, :, 3]
        d = sums[g, :, 4] + 1j * sums[g, :, 5]
        if np.any(a < 0) or np.any(a > samples):
            raise ValueError(f"N{n} {direction}: invalid rank-one counts")
        zero = a == 0
        if np.any(b[zero] != 0) or np.any(c[zero] != 0) or np.any(d[zero] != 0):
            raise ValueError(f"N{n} {direction}: zero support has nonzero marked products")
        z = float(w @ a)
        if z <= 0:
            raise ValueError(f"N{n} {direction}: no empirical rank-one mass")
        product = np.zeros(n + 1, dtype=complex)
        np.divide(b * c, a, out=product, where=a > 0)
        numerator = d - product
        # A singleton has zero sample covariance; this is not a population null.
        numerator[a <= 1] = 0
        mu = (w @ b) / z
        mean_s = float(w @ c) / z
        within = (w @ numerator) / z
        between = (w @ product) / z - mu * mean_s
        total = (w @ d) / z - mu * mean_s
        correction = np.zeros(n + 1)
        np.divide(a, a - 1, out=correction, where=a > 1)
        bessel = (w @ (numerator * correction)) / z

        def saved(quantity):
            return complex(primary_values[f"N{n}.{direction}.{quantity}_re"],
                           primary_values[f"N{n}.{direction}.{quantity}_im"])

        slope = saved("mu_p")
        rootclock = rootdot * slope
        moving = saved("nu_s")
        quantities = {"within": within, "between": between, "Ctotal": total,
                      "moving_nu": moving, "rootclock": rootclock,
                      "within_bessel_sensitivity": bessel,
                      "bessel_minus_plugin": bessel - within,
                      "mu": mu, "mu_p": slope,
                      "primary_pred_E": saved("pred_E"),
                      "primary_resid_E": saved("resid_E")}
        for quantity, value in quantities.items():
            values[f"{direction}.{quantity}_re"] = float(value.real)
            values[f"{direction}.{quantity}_im"] = float(value.imag)
        values[f"{direction}.rank1_probability"] = z / samples
        small = a <= 1
        positive = a[a > 0]
        diagnostics[direction] = {
            "samples": int(samples), "rank1_empirical_mass": z / samples,
            "binomial_weighted_rank1_count_not_independent_ESS": z,
            "rank1_counts_by_K": a.astype(np.int64).tolist(),
            "binomial_mass_a_zero": float(w[zero].sum()),
            "binomial_mass_a_le_one": float(w[small].sum()),
            "rank1_empirical_mass_a_le_one": float(w[small] @ a[small]) / samples,
            "rank1_mass_fraction_a_le_one": float(w[small] @ a[small]) / z,
            "rank1_mass_fraction_bessel_supported": float(w[a > 1] @ a[a > 1]) / z,
            "minimum_positive_rank1_count": int(positive.min()),
            "diagnostic_binomial_quantile_range": [lower, upper],
            "minimum_rank1_count_in_diagnostic_quantile_range": int(a[diagnostic_window].min()),
            "diagnostic_quantile_range_is_not_an_analysis_cut": True,
            "conditional_source_mean_bulk": mean_s,
            "algebraic_total_minus_within_between_abs": float(abs(total - within - between)),
            "same_primary_conditional_covariance_difference_abs": float(abs(total - saved("conditional_cov_s"))),
            "same_primary_moving_minus_total_rootclock_abs": float(abs(moving - total - rootclock)),
            "singleton_boundary": "a=1 gives an empirical zero within, not a population-zero claim",
        }
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError("nonfinite fixed-K decomposition")
    return values, {"N": n, "p0": p0, "binomial_weights": w.tolist(), "directions": diagnostics}


def main():
    started = time.perf_counter()
    contract = json.loads(CONTRACT.read_text())
    destination = ROOT / contract["output"]
    if destination.exists():
        raise ValueError("fixed-K result exists; refusing overwrite or repeated aggregation")
    if tuple(contract["Ns"]) != NS or contract["batches"] != 100:
        raise ValueError("fixed-K contract must preserve six N and 100 aligned batches")
    primary_path, source_path = ROOT / contract["primary_result"], ROOT / contract["source_result"]
    prior = json.loads(primary_path.read_text())
    source = json.loads(source_path.read_text())
    if prior["execution_commit"] != contract["primary_execution_commit"]:
        raise ValueError("primary result differs from the fixed-K contract")
    if primary.digest(source_path) != contract["source_result_sha256"]:
        raise ValueError("saved root/source result differs from the fixed-K contract")
    if prior["source_result"]["sha256"] != contract["source_result_sha256"]:
        raise ValueError("primary and fixed-K analysis must share the identical root/source archive")
    line_contract = prior["contract"]
    manifest_path = ROOT / prior["source_manifest"]["path"]
    manifest = json.loads(manifest_path.read_text())
    runs = {run["N"]: run for run in manifest["runs"]}
    prior_values = {label: float(prior["estimates"][label]["value"]) for label in prior["labels"]}
    source_index = {label: i for i, label in enumerate(source["labels"])}
    prior_inputs = {row["path"]: row["sha256"] for row in prior["inputs"]}
    per_batch = {n: int(line_contract["marked_permutations_by_N"][str(n)]) // 100 for n in NS}
    line_profiles, totals, inputs = {}, {}, []
    for n in NS:
        old_path = primary.old.OUTPUT / "raw" / f"n{n}.csv"
        source_profiles = primary.load_profile(old_path, n, 1000, runs[n])
        paths = [(old_path, "original_source_input_for_existing_loader")]
        if n in (260, 340):
            increment = ROOT / "results/norm4-source-endpoint-1m/increment/raw" / f"n{n}.csv"
            source_profiles += primary.load_profile(increment, n, 9000, runs[n])
            paths.append((increment, "same_nested_source_union_for_existing_loader"))
        line_path = primary_path.parent / "raw" / f"n{n}.csv"
        line_profiles[n], _ = primary.load_line_profile(line_path, n, per_batch[n], runs[n], source_profiles)
        paths.append((line_path, "same_existing_conditional_line_marks"))
        for path, role in paths:
            record = primary.input_record(path, role, N=n)
            if prior_inputs.get(record["path"]) != record["sha256"]:
                raise ValueError(f"input differs from primary result: {record['path']}")
            inputs.append(record)
        totals[n] = line_profiles[n].sum(axis=0)
    tail = float(contract["small_count_sensitivity"]["diagnostic_binomial_tail_probability"])
    points, diagnostics = {}, {}
    for n in NS:
        scalar = source["by_N"][str(n)]["points"]
        if float(scalar["p0"]) != prior_values[f"N{n}.p0"]:
            raise ValueError("primary and source central roots differ")
        points[n], diagnostics[n] = decompose(totals[n], 100 * per_batch[n], n,
            float(scalar["p0"]), float(scalar["rootdot_fugacity"]), prior_values, tail)
    central_map = primary.vectorize(points)
    labels = list(central_map)
    central = np.asarray(list(central_map.values()))
    covariance = np.zeros((len(labels), len(labels)))
    contributions, covered = {}, []
    for group in line_contract["dependency_groups"]:
        group_id, sizes = group["id"], group["Ns"]
        covered.extend(sizes)
        prior_group = prior["covariance_contributions"][group_id]
        source_group = source["covariance_contributions"]["source:" + group_id]
        if (prior_group["Ns"] != sizes or source_group["Ns"] != sizes
                or prior_group["delete_one_batch_ids"] != list(range(100))
                or source_group["delete_one_batch_ids"] != list(range(100))):
            raise ValueError("source/primary delete-one dependency alignment differs")
        for n in sizes:
            if (prior_group["batch_counts"] != [per_batch[n]] * 100
                    or source_group["batch_counts"] != [per_batch[n]] * 100):
                raise ValueError("source/line nested batch sample counts differ")
        prior_loo = np.asarray(prior_group["delete_one_vectors"])
        source_loo = np.asarray(source_group["delete_one_vectors"])
        if (prior_loo.shape != (100, len(prior["labels"]))
                or source_loo.shape != (100, len(source["labels"]))):
            raise ValueError("incomplete saved source or primary leave-one-out vectors")
        vectors, loo_diagnostics = [], []
        for batch in range(100):
            inherited = dict(zip(prior["labels"], map(float, prior_loo[batch])))
            changed = dict(points)
            batch_diagnostics = {}
            for n in sizes:
                p0 = float(source_loo[batch, source_index[f"N{n}.p0"]])
                rootdot = float(source_loo[batch, source_index[f"N{n}.rootdot_fugacity"]])
                if p0 != inherited[f"N{n}.p0"]:
                    raise ValueError("source/primary omitted root differs")
                changed[n], diag = decompose(totals[n] - line_profiles[n][batch],
                    99 * per_batch[n], n, p0, rootdot, inherited, tail)
                batch_diagnostics[n] = {g: {key: value for key, value in fields.items()
                    if key != "rank1_counts_by_K"}
                    for g, fields in diag["directions"].items()}
            vectors.append(list(primary.vectorize(changed).values()))
            loo_diagnostics.append(batch_diagnostics)
        vectors = np.asarray(vectors)
        centered = vectors - vectors.mean(axis=0)
        component = 99 / 100 * centered.T @ centered
        covariance += component
        contributions[group_id] = {
            "Ns": sizes, "batch_counts": prior_group["batch_counts"],
            "delete_one_batch_ids": list(range(100)), "delete_one_vectors": vectors.tolist(),
            "covariance": component.tolist(), "small_count_and_identity_diagnostics": loo_diagnostics,
            "operation": "Same existing source/line batch omitted; source p0/rootdot and primary conditional vectors reused, other dependency groups held central",
        }
    if sorted(covered) != sorted(NS):
        raise ValueError("the three dependency groups must cover each N once")
    errors = np.sqrt(np.maximum(0, np.diag(covariance)))
    estimates = {label: {"value": float(value), "se": float(error),
        "z": float(value / error) if error > 0 else None}
        for label, value, error in zip(labels, central, errors)}
    selections = {n: [f"N{n}.{g}.within_{part}" for g, part in COMPONENTS] for n in NS}
    by_n = {n: primary.selected_zero(labels, central, covariance, selections[n]) for n in NS}
    selected_all = [label for n in NS for label in selections[n]]
    joint = primary.selected_zero(labels, central, covariance, selected_all)
    result = {
        "schema": "matching-one.norm4-source-line-fixed-k.v1",
        "status": "computed_same_archive_fixed_K_covariance_decomposition",
        "execution_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "command": [sys.executable, *sys.argv], "contract": contract,
        "labels": labels, "estimates": estimates, "covariance": covariance.tolist(),
        "by_N": diagnostics, "covariance_contributions": contributions,
        "within_zero": {"by_N_four_real": by_n, "joint_24_real": joint,
                        "scope": contract["definitions"]["null"]},
        "summary": {"within_joint_nominal_p": joint.get("nominal_p"),
            "within_joint_supported_rank": joint.get("supported_rank"),
            "within_per_N_nominal_p": {n: by_n[n].get("nominal_p") for n in NS},
            "maximum_abs_within_component_z": max(abs(estimates[label]["z"] or 0) for label in selected_all)},
        "inputs": inputs,
        "primary_result": primary.input_record(primary_path, "same_primary_central_and_full_leave_one_out_vectors"),
        "source_result": primary.input_record(source_path, "same_saved_matching_roots_and_root_motion"),
        "source_manifest": primary.input_record(manifest_path, "unchanged_original_geometry_and_counter_contract"),
        "code": [primary.input_record(path, "analysis_code_or_contract") for path in
            (Path(__file__).resolve(), CONTRACT, Path(primary.__file__), ROOT / "scripts/norm4_source_line_core.py",
             ROOT / "scripts/analyze_norm4_source_thermal.py", ROOT / "scripts/analyze_norm4_source_endpoint_1m.py")],
        "uncertainty": "Full same-three-group aligned covariance includes within, between, total, copied moving response/primary E residual, rootclock and Bessel sensitivity; K rows and derived views are not independent evidence",
        "small_count_boundary": contract["small_count_sensitivity"],
        "interpretation": contract["interpretation"],
        "root_finder_calls": 0, "configuration_replays_by_this_script": 0,
        "new_samples": 0, "server_actions": 0, "test_suites": [],
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "machine": platform.machine()},
        "elapsed_seconds": time.perf_counter() - started,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("x") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"output": primary.display(destination), "elapsed_seconds": result["elapsed_seconds"],
                      **result["summary"]}, indent=2))


if __name__ == "__main__":
    main()
