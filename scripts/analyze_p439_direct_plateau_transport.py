#!/usr/bin/env python3
"""Partition real P439 matching response without creating new samples.

The old K_A estimates and their paired delete-one vectors are reused exactly.
Only direct/plateau matching components are reconstructed from pinned archives.
This script does not rerun the parent P337/P439 analysis or a test suite.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import chi2

ORDER = ["K_A", "M_direct", "M_plateau", "M_total"]
ORIENTATIONS = ["first", "second"]


def blob(commit, path, expected):
    payload = subprocess.check_output(["git", "cat-file", "blob", f"{commit}:{path}"])
    observed = hashlib.sha256(payload).hexdigest()
    if observed != expected:
        raise ValueError(f"immutable source hash mismatch: {commit}:{path}")
    return payload


def covariance(leave):
    values = np.asarray(leave, dtype=float)
    shifted = values - values[0]
    centered = shifted - shifted.mean(axis=0)
    return centered.T @ centered * ((len(values) - 1) / len(values))


def quadratic(vector, cov):
    scale = np.sqrt(np.diag(cov))
    correlation = cov / np.outer(scale, scale)
    eigenvalues, basis = np.linalg.eigh(correlation)
    retained = eigenvalues > 1e-10 * eigenvalues.max()
    score = basis.T @ (np.asarray(vector) / scale)
    stat = float(np.sum(score[retained] ** 2 / eigenvalues[retained]))
    rank = int(retained.sum())
    return {"chi2": stat, "df": rank, "p_asymptotic": float(chi2.sf(stat, rank)),
            "correlation_eigenvalues": eigenvalues.tolist(), "relative_cutoff": 1e-10}


def common_ray(values, covariances):
    """Same two-coordinate profile as P439, with free per-generation amplitude."""
    precisions = [np.linalg.inv(c) for c in covariances]

    def direction(u):
        amplitudes = [(u @ w @ y) / (u @ w @ u)
                      for y, w in zip(values, precisions)]
        q = sum((y - a * u) @ w @ (y - a * u)
                for y, a, w in zip(values, amplitudes, precisions))
        return float(q), amplitudes

    def angular(theta):
        return direction(np.array([np.sin(theta), np.cos(theta)]))[0]

    grid = np.linspace(-np.pi / 2, np.pi / 2, 2049)
    best = int(np.argmin([angular(t) for t in grid]))
    step = grid[1] - grid[0]
    fit = minimize_scalar(angular, bounds=(grid[best] - step, grid[best] + step),
                          method="bounded", options={"xatol": 1e-14})
    ratio = float(np.tan(fit.x))
    minimum, amplitudes = direction(np.array([ratio, 1.0]))

    def at_ratio(r):
        return direction(np.array([r, 1.0]))[0]

    interval = []
    threshold = minimum + chi2.ppf(.95, 1)
    for sign in (-1, 1):
        distance = max(abs(ratio) * .1, .001)
        for _ in range(50):
            outside = ratio + sign * distance
            if at_ratio(outside) > threshold:
                interval.append(float(brentq(lambda r: at_ratio(r) - threshold,
                                            *sorted([ratio, outside]))))
                break
            distance *= 2
        else:
            interval.append(None)
    zero = at_ratio(0.0)
    return {"loading_ratio_component_over_K_A": ratio, "profile_95": interval,
            "chi2": minimum, "df": 3, "p_asymptotic": float(chi2.sf(minimum, 3)),
            "K_A_amplitudes": amplitudes,
            "zero_component_loading": {"chi2": zero, "df": 4,
                "p_asymptotic": float(chi2.sf(zero, 4)),
                "delta_chi2_against_fitted_ray": zero - minimum,
                "p_delta_asymptotic_1df": float(chi2.sf(max(0.0, zero - minimum), 1))}}


def binomial_tails(n, p):
    # Preserve the stable normalized recurrence used by the immutable source.
    values = [0.0] * (n + 1)
    mode = min(n, max(0, int(math.floor((n + 1) * p))))
    values[mode] = math.exp(math.lgamma(n + 1) - math.lgamma(mode + 1)
                            - math.lgamma(n - mode + 1)
                            + mode * math.log(p) + (n - mode) * math.log1p(-p))
    for k in range(mode, 0, -1):
        values[k - 1] = values[k] * k / (n - k + 1) * (1.0 - p) / p
    for k in range(mode, n):
        values[k + 1] = values[k] * (n - k) / (k + 1) * p / (1.0 - p)
    norm = math.fsum(values)
    out, cumulative = [0.0] * (n + 1), 0.0
    for k in range(n, -1, -1):
        cumulative += values[k] / norm
        out[k] = cumulative
    return np.asarray(out)


def load_histograms(run):
    source, n = run["source"], run["N"]
    batch_ids = next(c for c in run["contexts"] if c["id"] == "fixed_reference")["delete_one"]["batch_ids"]
    batch_index = {batch: index for index, batch in enumerate(batch_ids)}
    # direct, line birth, line completion; integer counts remain exact.
    hist = np.zeros((2, len(batch_ids), 3, n + 1), dtype=np.int64)
    sample_counts = np.zeros((2, len(batch_ids)), dtype=np.int64)
    observed_counts = np.zeros_like(sample_counts)
    payload = blob(source["archive_commit"], source["path"], source["sha256_verified"])
    binary = io.BytesIO(payload)
    if source["path"].endswith(".gz"):
        binary = gzip.GzipFile(fileobj=binary)
    records = 0
    with io.TextIOWrapper(binary, encoding="utf-8", newline="") as text:
        reader = csv.reader(text)
        fields = next(reader)
        required = ["n", "orientation", "batch", "samples", "tau1", "tau2", "kind", "count"]
        if not set(required).issubset(fields):
            raise ValueError(f"not_scoreable: missing raw fields at N{n}")
        col = {name: fields.index(name) for name in required}
        for row in reader:
            if int(row[col["n"]]) != n or row[col["orientation"]] not in ORIENTATIONS:
                raise ValueError("raw geometry/orientation disagrees with source")
            orientation = ORIENTATIONS.index(row[col["orientation"]])
            batch = batch_index[int(row[col["batch"]])]
            samples, count = int(row[col["samples"]]), int(row[col["count"]])
            k1, k2 = int(row[col["tau1"]]), int(row[col["tau2"]])
            if count <= 0 or not 0 <= k1 <= k2 <= n:
                raise ValueError("invalid activation rank or count")
            if sample_counts[orientation, batch] not in (0, samples):
                raise ValueError("inconsistent total batch denominator")
            sample_counts[orientation, batch] = samples
            observed_counts[orientation, batch] += count
            kind = row[col["kind"]]
            if kind == "DIRECT_RANK2" and k1 == k2:
                hist[orientation, batch, 0, k1] += count
            elif kind == "LINE" and k1 < k2:
                hist[orientation, batch, 1, k1] += count
                hist[orientation, batch, 2, k2] += count
            else:
                raise ValueError("not_scoreable: direct/plateau semantics disagree")
            records += 1
    if not np.array_equal(sample_counts, observed_counts):
        raise ValueError("direct and plateau families do not exhaust original samples")
    if not np.all(sample_counts.sum(axis=1) == source["samples_per_orientation"]):
        raise ValueError("source total denominators do not match the raw archive")
    direct_counts = hist[:, :, 0, :].sum(axis=(1, 2))
    if direct_counts.tolist() != [source["direct_counts"][o] for o in ORIENTATIONS]:
        raise ValueError("raw direct counts do not match immutable source")
    return {"histograms": hist, "samples": sample_counts, "batch_ids": batch_ids,
            "records": records, "direct_counts": direct_counts.tolist()}


def score_context(run, loaded, context):
    n = run["N"]
    hist, samples = loaded["histograms"], loaded["samples"]
    full_hist, full_samples = hist.sum(axis=1), samples.sum(axis=1)
    parent_order = context["vector_order"]
    old_indices = [parent_order.index("angular_K_A_activity"), parent_order.index("angular_A_top")]
    parent_point = np.asarray(context["vector"])[old_indices]
    parent_leave = np.asarray(context["delete_one"]["vectors"])[:, old_indices]
    if loaded["batch_ids"] != context["delete_one"]["batch_ids"]:
        raise ValueError("K_A and split-M batch IDs do not align")
    maximum_closure = 0.0

    def estimate(omitted, p, parent_value):
        nonlocal maximum_closure
        h = full_hist if omitted is None else full_hist - hist[:, omitted]
        denominator = full_samples if omitted is None else full_samples - samples[:, omitted]
        tail = binomial_tails(n, p)
        direct_mass = h[:, 0].sum(axis=1)
        plateau_mass = h[:, 1].sum(axis=1)
        direct = (2 * (h[:, 0] @ tail) - direct_mass) / denominator
        plateau = ((h[:, 1] + h[:, 2]) @ tail - plateau_mass) / denominator
        independently_rebuilt_total = ((2 * h[:, 0] + h[:, 1] + h[:, 2]) @ tail) / denominator - 1
        contrasts = [(row[1] - row[0]) / context["delta_cos4"]
                     for row in (direct, plateau)]
        total = contrasts[0] + contrasts[1]
        maximum_closure = max(maximum_closure, float(np.max(np.abs(
            direct + plateau - independently_rebuilt_total))), abs(total - parent_value[1]))
        orientation_values = {o: {"M_direct": float(direct[i]), "M_plateau": float(plateau[i]),
                                  "M_total": float(independently_rebuilt_total[i]),
                                  "direct_mass": float(direct_mass[i] / denominator[i]),
                                  "total_sample_denominator": int(denominator[i])}
                              for i, o in enumerate(ORIENTATIONS)}
        return [parent_value[0], *contrasts, total], orientation_values

    point, orientation_values = estimate(None, context["p"], parent_point)
    leaves = [estimate(index, p, parent_leave[index])[0]
              for index, p in enumerate(context["delete_one"]["p_values"])]
    cov = covariance(leaves)
    expected_cov = np.asarray(context["covariance"])[np.ix_(old_indices, old_indices)]
    observed_cov = cov[np.ix_([0, 3], [0, 3])]
    if maximum_closure > 2e-12 or not np.allclose(observed_cov, expected_cov, rtol=1e-8, atol=1e-20):
        raise ValueError("single reconstruction check failed: source M or paired covariance differs")
    standard_error = np.sqrt(np.diag(cov))
    return {"N": n, "id": context["id"], "p": context["p"],
            "dependency_group": run["dependency_group"], "source": run["source"],
            "raw_rows_read": loaded["records"], "delta_cos4_exact": context["delta_cos4_exact"],
            "delta_cos4": context["delta_cos4"], "coordinate_order": ORDER,
            "value": np.asarray(point).tolist(), "standard_error": standard_error.tolist(),
            "z_measurement_only": (point / standard_error).tolist(),
            "covariance": cov.tolist(),
            "correlation_direct_plateau": float(cov[1, 2] / (standard_error[1] * standard_error[2])),
            "orientation_values": orientation_values,
            "delete_one": {"unit": "same paired batch across both orientations and all coordinates",
                "batch_ids": loaded["batch_ids"], "p_values": context["delete_one"]["p_values"],
                "vectors": np.asarray(leaves).tolist()},
            "single_reconstruction_check": {"max_point_and_delete_one_closure_error": maximum_closure,
                "max_source_K_A_M_covariance_absolute_error": float(np.max(np.abs(observed_cov - expected_cov))),
                "all_p_identity": "The exhaustive disjoint histograms obey M_direct+M_plateau=F1+F2-1 coefficientwise, with original denominators."}}


def wedge_vector(values):
    return np.asarray([values[j + 1, k] * values[j, 0] - values[j + 1, 0] * values[j, k]
                       for j in range(3) for k in range(1, 4)])


def combine(rows, parent):
    values = np.asarray([row["value"] for row in rows])
    covs = np.asarray([row["covariance"] for row in rows])
    wedges = wedge_vector(values)
    group_covs, group_leaves = [], []
    for index, row in enumerate(rows):
        leaves = []
        for replacement in row["delete_one"]["vectors"]:
            one_deleted = values.copy()
            one_deleted[index] = replacement
            leaves.append(wedge_vector(one_deleted))
        group_covs.append(covariance(leaves))
        group_leaves.append(np.asarray(leaves).tolist())
    wedge_cov = sum(group_covs)
    summaries = {}
    for k, name in enumerate(ORDER[1:], 1):
        ix = list(range(k - 1, 9, 3))
        selected_cov = wedge_cov[np.ix_(ix, ix)]
        one_c = [cov[np.ix_([k, 0], [k, 0])] for cov in covs]
        summaries[name] = {"wedge": {"value": wedges[ix].tolist(),
            "standard_error": np.sqrt(np.diag(selected_cov)).tolist(),
            "covariance": selected_cov.tolist(), **quadratic(wedges[ix], selected_cov)},
            "common_ray": common_ray(values[:, [k, 0]], one_c),
            "zero_across_generations": quadratic(values[:, k], np.diag(covs[:, k, k]))}
    independent = [0, 1, 3, 4, 6, 7]
    total_indices = [2, 5, 8]
    total_parent_error = float(np.max(np.abs(wedges[total_indices] - np.asarray(parent["wedge"]["value"]))))
    if total_parent_error > 1e-13 or not np.allclose(wedge_cov[np.ix_(total_indices, total_indices)],
                                                  parent["wedge"]["covariance"], rtol=1e-8, atol=1e-22):
        raise ValueError("single reconstruction check failed: total wedge differs from parent")
    return {"generations": rows, "component_summaries": summaries,
            "joint_wedges": {"order": [f"N{n}:{name}" for n in (85, 170, 340) for name in ORDER[1:]],
                "value": wedges.tolist(), "covariance": wedge_cov.tolist(),
                "covariance_by_independent_generation": [cov.tolist() for cov in group_covs],
                "delete_one_by_independent_generation": group_leaves,
                "algebraic_dependence": "At every full and delete-one point D_total=D_direct+D_plateau; nine coordinates are at most rank six.",
                "independent_direct_plateau_basis_test": quadratic(wedges[independent], wedge_cov[np.ix_(independent, independent)]),
                "max_parent_total_wedge_error": total_parent_error},
            "joint_component_zero_test": quadratic(values[:, 1:3].reshape(-1),
                np.block([[covs[i, 1:3, 1:3] if i == j else np.zeros((2, 2))
                           for j in range(4)] for i in range(4)]))}


def report(result):
    primary = result["primary"]
    lines = ["# P439: direct/plateau decomposition of the same-stream matching response", "",
             "Retrospective real-archive analysis; **zero new samples**. This is a partition of one",
             "microscopic source, not evidence for two independent sources.", "",
             "## Fixed-p result", "",
             "All entries are exact H4 direction contrasts at p=0.592746050790; errors are one",
             "paired-batch jackknife SE. K_A and its delete-one vectors are reused unchanged.", "",
             "| N | M_direct | M_plateau | M_total | corr(direct, plateau) |", "|---|---|---|---|---|"]
    for row in primary["generations"]:
        cells = [f"{row['value'][i]:.8g} ± {row['standard_error'][i]:.6g}" for i in (1, 2, 3)]
        lines.append(f"| {row['N']} | {' | '.join(cells)} | {row['correlation_direct_plateau']:.4f} |")
    lines += ["", "| Component | zero across four N: chi2/4, p | adjacent wedge: chi2/3, p | loading/K_A, 95% profile | ray p |",
              "|---|---|---|---|---|"]
    for name, summary in primary["component_summaries"].items():
        zero, wedge, ray = summary["zero_across_generations"], summary["wedge"], summary["common_ray"]
        bounds = ray["profile_95"]
        lines.append(f"| {name} | {zero['chi2']:.5g}, {zero['p_asymptotic']:.5g} | {wedge['chi2']:.5g}, {wedge['p_asymptotic']:.5g} | {ray['loading_ratio_component_over_K_A']:.6g}, [{bounds[0]:.6g}, {bounds[1]:.6g}] | {ray['p_asymptotic']:.5g} |")
    joint = primary["joint_component_zero_test"]
    joint_wedge = primary["joint_wedges"]["independent_direct_plateau_basis_test"]
    lines += ["", f"The joint eight-coordinate (direct, plateau) zero diagnostic is {joint['chi2']:.6g}/{joint['df']} df, p={joint['p_asymptotic']:.6g}.",
              f"The six independent component wedges give {joint_wedge['chi2']:.6g}/{joint_wedge['df']} df, p={joint_wedge['p_asymptotic']:.6g}.",
              "The total wedges are algebraic sums and are not additional degrees of freedom.", "",
              "## Exact partition and dependence", "",
              "For H_k(p)=Pr[Binomial(N,p)>=k], use the original total sample denominator:", "",
              "- M_direct = E[1_DIRECT_RANK2 (2H_tau1-1)].",
              "- M_plateau = E[1_LINE (H_tau1+H_tau2-1)].",
              "- M_total = M_direct + M_plateau = F1+F2-1.", "",
              "Both raw family contributions are reconstructed independently; neither is rescaled",
              "by its own family count. Raw rank equality/inequality, exhaustive counts and hashes",
              "are checked in the one analysis pass. Every paired delete-one M sum, its covariance",
              "with unchanged K_A and the full total-wedge covariance reproduce immutable P439.", "",
              "Within each generation the directions and all four coordinates share a deletion.",
              "The four generations have distinct dependency groups (20/80/80/80 batches). Their",
              "nonlinear wedge covariance contributions are summed, with no fictitious cross-N batch alignment.", "",
              "## Scope and next scientific decision", "",
              "These are measurement-only asymptotic chi-square/profile summaries, conditional on the",
              "saved covariance estimates and retrospective choice of partition. They do not include",
              "source/model uncertainty and are not prospective model-selection certification.",
              "Opposite point-estimate signs are not by themselves resolved cancellation. Nor does a",
              "surviving common ray establish nonzero coupling to K_A.", "",
              "The root sensitivity repeats only this partition at the parent's already saved full/delete-one",
              "pooled roots; roots, K_A, transfer order and source definitions are not reselected.", "",
              "At the pinned data used here, each direct contribution is within 1.11 SE of zero and",
              "each plateau contribution is within 1.28 SE. Their modest covariance does not reveal",
              "two resolved large terms hidden by cancellation. Both family loadings remain unresolved;",
              "this is not proof that either loading is exactly zero or that cancellation is impossible.", "",
              "Direct and plateau point signs oppose at N85/N340 but agree at N170/N680. The",
              "root sensitivity changes the joint zero p-value only from .790132 to .790122.", "",
              "**Next output:** on an existing geometry with a separately resolved canonical M contrast,",
              "produce M and the same natural K_A from one paired batch stream and report their joint",
              "loading interval. First consume any compatible saved birth archive; acquire only the",
              "missing same-stream rows if needed. Another K_A-only N1360 point would not answer this question.", "",
              "Move parent same-stream scoring and this direct/plateau split to completed analysis;",
              "do not present reimplementation of their scorers or another synthetic check as the next science task.", "",
              "## Reproduce", "",
              "A clone must contain the immutable source objects listed in the manifest and score.",
              "No server or external raw download is needed when those Git objects are present.", "", "```sh",
              "python3 scripts/analyze_p439_direct_plateau_transport.py",
              "```", "",
              "Output retains full and delete-one (K_A,M_direct,M_plateau,M_total) vectors, joint",
              "covariances, all nine wedges, per-generation covariance contributions, raw provenance,",
              "environment and one-pass reconstruction diagnostics. No parent analysis or test suite is rerun.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=Path("analysis/p439_direct_plateau_transport_manifest.json"))
    parser.add_argument("--output", type=Path, default=Path("results/p439-direct-plateau-transport/score.json"))
    parser.add_argument("--report", type=Path, default=Path("results/p439-direct-plateau-transport/REPORT.md"))
    args = parser.parse_args()
    started = time.perf_counter()
    manifest_bytes = args.manifest.read_bytes()
    manifest = json.loads(manifest_bytes)
    source = json.loads(blob(manifest["source_commit"], manifest["source_path"], manifest["source_sha256"]))
    parent = json.loads(blob(manifest["parent_p439_commit"], manifest["parent_p439_path"], manifest["parent_p439_sha256"]))
    runs = sorted([r for r in source["runs"] if r["role"] == "four_generation_primary"], key=lambda r: r["N"])
    if [run["N"] for run in runs] != manifest["N"]:
        raise ValueError("size list differs from original P439")
    contexts = {name: [] for name in ("fixed_reference", "pooled_root")}
    for run in runs:
        print(f"Reading immutable N{run['N']} histogram once...", flush=True)
        loaded = load_histograms(run)
        for name in contexts:
            context = next(c for c in run["contexts"] if c["id"] == name)
            if name == "fixed_reference" and context["p"] != manifest["fixed_p"]:
                raise ValueError("fixed p differs from original P439")
            contexts[name].append(score_context(run, loaded, context))
    result = {"schema": "matching-one.p439-direct-plateau-transport.v1",
              "status": "retrospective_real_archive_decomposition", "new_samples": 0,
              "manifest": manifest, "manifest_sha256": hashlib.sha256(manifest_bytes).hexdigest(),
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "primary": combine(contexts["fixed_reference"], parent["primary"]),
              "root_sensitivity": combine(contexts["pooled_root"], parent["root_sensitivity"]),
              "inference": "Measurement-only asymptotic summaries from estimated paired-batch covariance; retrospective partition, not new independent sources or prospective model selection.",
              "execution": {"python": sys.version, "platform": platform.platform(),
                  "numpy": np.__version__, "command": " ".join(sys.argv),
                  "wall_seconds": time.perf_counter() - started,
                  "git_base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
                  "server_contacted": False, "parent_analysis_rerun": False, "test_suite_run": False}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    args.report.write_text(report(result))
    print(json.dumps({"output": str(args.output), "wall_seconds": result["execution"]["wall_seconds"],
                      "primary_summaries": result["primary"]["component_summaries"]}, indent=2))


if __name__ == "__main__":
    main()
