#!/usr/bin/env python3
"""Score the real four-generation same-stream crosswalk, without new sampling.

The immutable P337 artifact already contains the exact canonical reconstruction
and aligned delete-one vectors. Generations are independent evidence blocks;
batch labels are aligned across orientations, NOT across generations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import chi2

SOURCE_COMMIT = "6123955fa1850506d568753ca2d178ab8aeccac8"
SOURCE_PATH = "results/p337-birth-state-current/latest.json"
SOURCE_SHA256 = "7bb3a048b16c0436cf2e28ae6a87eff3c4f388d74f099ea3d7cda579aab8793d"
SIZES = [85, 170, 340, 680]


def covariance(leave):
    leave = np.asarray(leave, dtype=float)
    shifted = leave - leave[0]
    centered = shifted - shifted.mean(axis=0)
    return (centered.T @ centered) * ((len(leave) - 1) / len(leave))


def quadratic(vector, cov):
    """Rank decision in correlation coordinates, not dimensional units."""
    scale = np.sqrt(np.diag(cov))
    corr = cov / np.outer(scale, scale)
    eigenvalues, basis = np.linalg.eigh(corr)
    keep = eigenvalues > 1e-10 * eigenvalues.max()
    z = basis.T @ (vector / scale)
    stat = float(np.sum(z[keep] ** 2 / eigenvalues[keep]))
    rank = int(keep.sum())
    return {"chi2": stat, "df": rank, "p_asymptotic": float(chi2.sf(stat, rank)),
            "correlation_eigenvalues": eigenvalues.tolist(), "relative_cutoff": 1e-10}


def wedges(y):
    return y[1:, 0] * y[:-1, 1] - y[1:, 1] * y[:-1, 0]


def common_ray(y, covs):
    precisions = [np.linalg.inv(c) for c in covs]

    def fit_direction(u):
        amplitudes = np.array([(u @ w @ point) / (u @ w @ u)
                               for point, w in zip(y, precisions)])
        q = sum((point - a*u) @ w @ (point - a*u)
                for point, a, w in zip(y, amplitudes, precisions))
        return float(q), amplitudes

    def angle_score(theta):
        return fit_direction(np.array([np.sin(theta), np.cos(theta)]))[0]

    grid = np.linspace(-np.pi / 2, np.pi / 2, 2049)
    best = int(np.argmin([angle_score(theta) for theta in grid]))
    step = grid[1] - grid[0]
    fit = minimize_scalar(angle_score, bounds=(grid[best]-step, grid[best]+step),
                          method="bounded", options={"xatol": 1e-14})
    ratio = float(np.tan(fit.x))

    def ratio_score(r):
        return fit_direction(np.array([r, 1.0]))[0]

    minimum, amplitudes = fit_direction(np.array([ratio, 1.0]))

    def interval(level):
        threshold = minimum + chi2.ppf(level, 1)
        bounds = []
        for direction in (-1, 1):
            distance = max(abs(ratio) * .1, .001)
            for _ in range(50):
                outside = ratio + direction * distance
                if ratio_score(outside) > threshold:
                    bounds.append(float(brentq(lambda r: ratio_score(r)-threshold,
                                               *sorted([ratio, outside]))))
                    break
                distance *= 2
            else:
                bounds.append(None)
        return bounds

    zero = ratio_score(0.0)
    return {"loading_ratio_first_over_K": ratio, "K_amplitudes": amplitudes.tolist(),
            "chi2": minimum, "df": 3, "p_asymptotic": float(chi2.sf(minimum, 3)),
            "ratio_profile_95": interval(.95), "ratio_profile_99": interval(.99),
            "zero_first_loading": {"chi2": zero, "df": 4,
                "p_asymptotic": float(chi2.sf(zero, 4)),
                "delta_chi2_against_fitted_ray": zero-minimum,
                "p_delta_asymptotic_1df": float(chi2.sf(max(0.0, zero-minimum), 1))}}


def score(source, context_id, observable, save_replicates=False):
    runs = sorted([r for r in source["runs"] if r["role"] == "four_generation_primary"],
                  key=lambda r: r["N"])
    assert [r["N"] for r in runs] == SIZES
    assert len({r["dependency_group"] for r in runs}) == 4
    rows, means, covs, leaves = [], [], [], []
    for run in runs:
        context = next(c for c in run["contexts"] if c["id"] == context_id)
        names = [observable, "angular_K_A_activity"]
        ix = [context["vector_order"].index(name) for name in names]
        mean = np.asarray(context["vector"])[ix]
        cov = np.asarray(context["covariance"])[np.ix_(ix, ix)]
        leave = np.asarray(context["delete_one"]["vectors"])[:, ix]
        # One proportionality check of reused sufficient statistics, not raw replay.
        np.testing.assert_allclose(covariance(leave), cov, rtol=1e-9, atol=1e-20)
        for orientation in context["orientation_values"].values():
            assert abs(orientation["A_top"] - (orientation["F1"]+orientation["F2"]-1)) < 1e-13
        means.append(mean); covs.append(cov); leaves.append(leave)
        row = {"N": run["N"], "p": context["p"], "dependency_group": run["dependency_group"],
               "source": run["source"], "delta_cos4_exact": context["delta_cos4_exact"],
               "coordinates": names, "value": mean.tolist(), "covariance": cov.tolist(),
               "standard_error": np.sqrt(np.diag(cov)).tolist(),
               "z": (mean/np.sqrt(np.diag(cov))).tolist(),
               "correlation": float(cov[0, 1]/np.sqrt(cov[0, 0]*cov[1, 1]))}
        if save_replicates:
            row["batch_ids"] = context["delete_one"]["batch_ids"]
            row["delete_one_vectors"] = leave.tolist()
        rows.append(row)
    means = np.asarray(means)
    group_covs = []
    for i, leave in enumerate(leaves):
        replicate_wedges = []
        for replacement in leave:
            replicate = means.copy()
            replicate[i] = replacement
            replicate_wedges.append(wedges(replicate))
        group_covs.append(covariance(replicate_wedges))
    wedge_cov = sum(group_covs)
    value = wedges(means)
    return {"context": context_id, "generations": rows,
            "wedge": {"N": SIZES[:-1], "value": value.tolist(),
                "standard_error": np.sqrt(np.diag(wedge_cov)).tolist(),
                "covariance": wedge_cov.tolist(),
                "covariance_by_independent_generation": [c.tolist() for c in group_covs],
                **quadratic(value, wedge_cov)},
            "common_ray": common_ray(means, covs)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="optional downloaded immutable source JSON")
    parser.add_argument("--output", type=Path, default=Path("results/p439-same-stream-production/latest.json"))
    args = parser.parse_args()
    payload = (args.input.read_bytes() if args.input else subprocess.check_output(
        ["git", "cat-file", "blob", f"{SOURCE_COMMIT}:{SOURCE_PATH}"]))
    assert hashlib.sha256(payload).hexdigest() == SOURCE_SHA256
    source = json.loads(payload)
    result = {"schema": "matching-one.p439-same-stream-production.v1",
              "status": "retrospective_real_archive_score", "new_samples": 0,
              "source": {"commit": SOURCE_COMMIT, "path": SOURCE_PATH, "sha256": SOURCE_SHA256},
              "covariance_contract": "Aligned orientation deletion within each generation; independent generation covariance contributions summed after full nonlinear wedge recomputation.",
              "inference": "Measurement-only asymptotic chi-square/profile diagnostics; estimated batch covariances are not exact known covariances. No source/model uncertainty was added.",
              "primary": score(source, "fixed_reference", "angular_A_top", True),
              "root_sensitivity": score(source, "pooled_root", "angular_A_top", True),
              "exploratory_not_primary": {
                  name: score(source, "fixed_reference", name) for name in
                  ("angular_A_log_derivative", "angular_E_top")},
              "conclusion": "Common ray not rejected, but a zero M loading also survives. These data do not establish that unmarked M shares the marked-current radial state.",
              "not_done": ["prospective common transfer model selection",
                           "direct02/plateau conditional M decomposition", "N1360 acquisition"]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"output": str(args.output), "wedge": result["primary"]["wedge"],
                      "ray": result["primary"]["common_ray"]}, indent=2))


if __name__ == "__main__":
    main()
