#!/usr/bin/env python3
"""Frozen D-only shape readouts beyond a free center and width, no new MC."""
from fractions import Fraction
import json
from math import comb
from pathlib import Path

import numpy as np
from scipy.stats import chi2

from p267_scalar_clock_transport import load_source, landmarks, exact_empirical_root_intervals

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"results/p267-standardized-rank-shape"
PROTOCOL = ROOT/"experiments/p267_standardized_rank_shape_20260831.json"
ORDERS = (3, 4, 5, 6)
LABELS = ([f"rank_step.standardized_mu{r}" for r in ORDERS]
          +[f"canonical.standardized_mu{r}" for r in ORDERS]
          +[f"canonical.standardized_p.{s}" for s in ("first_peak", "valley", "second_peak")]
          +[f"canonical.unit_area_scaled_height.{s}" for s in ("first_peak", "valley", "second_peak")]
          +["canonical.valley_over_first_peak", "canonical.second_over_first_peak"])


def beta_raw_moments(n):
    j = np.arange(n+1)
    rows = [np.ones(n+1)]
    for r in range(1, 7):
        rows.append(rows[-1]*(j+r)/(n+1+r))
    return np.array(rows)


def shape_moments(f, raw_beta):
    n, j = len(f)-1, np.arange(len(f))
    weights = f/f.sum()
    jm = weights @ j
    jv = weights @ ((j-jm)**2)
    mu_s, var_s = (jm+.5)/n, (jv+1/12)/n**2
    mu_c = (jm+1)/(n+2)
    var_c = jv/(n+2)**2+(weights@((j+1)*(n-j+1)))/((n+2)**2*(n+3))
    if min(var_s, var_c) <= 0:
        raise ValueError("Nonpositive signed central variance; frozen normalization is undefined")
    left, right = j/n-mu_s, (j+1)/n-mu_s
    step = [float(weights @ (n*(right**(r+1)-left**(r+1))/(r+1))/var_s**(r/2)) for r in ORDERS]
    raw = raw_beta @ weights
    canonical = [sum(comb(r, k)*(-mu_c)**(r-k)*raw[k] for k in range(r+1))/var_c**(r/2)
                 for r in ORDERS]
    return np.array(step), np.array(canonical), mu_c, var_c


def features(f, raw_beta, brackets):
    n = len(f)-1
    step, canonical, mu_c, var_c = shape_moments(f, raw_beta)
    positions, heights, curvatures = landmarks(f, brackets)
    standardized_positions = (positions-mu_c)/np.sqrt(var_c)
    # Canonical integral is sum(f)/(n+1).
    standardized_heights = heights*np.sqrt(var_c)*(n+1)/f.sum()
    ratios = heights[1:]/heights[0]
    return np.r_[step, canonical, standardized_positions, standardized_heights, ratios], curvatures


def jackknife_covariance(loo):
    b, centered = len(loo), loo-loo.mean(axis=0)
    return (b-1)/b*(centered.T@centered)


def main():
    protocol = json.loads(PROTOCOL.read_text())
    sources, means, covariances = {}, {}, {}
    for area, commit in protocol["sources"].items():
        n = int(area)
        contract, hashes, bernstein, integers = load_source({"source_commit": commit, "source_directory": f"results/etop-n{n}-three-modulus"})
        batches, int_coefficients = bernstein[0, :, 0], integers[0, :, 0].sum(axis=0)
        intervals = exact_empirical_root_intervals(np.diff(int_coefficients))
        if len(intervals) != 3:
            raise ValueError("Ordered critical pattern changed; do not force the frozen three-landmark auxiliary")
        brackets = [[float(Fraction(a)), float(Fraction(b))] for a, b in intervals]
        b, mean, raw_beta = len(batches), batches.mean(axis=0), beta_raw_moments(n)
        point, curvature = features(mean, raw_beta, brackets)
        loo, curvature_signs = [], []
        for row in batches:
            v, cv = features((b*mean-row)/(b-1), raw_beta, brackets)
            loo.append(v)
            curvature_signs.append(np.sign(cv))
        loo = np.array(loo)
        cov = jackknife_covariance(loo)
        means[n], covariances[n] = point, cov
        sources[area] = {"source_commit": commit, "source_sha256": hashes, "source_contract": contract,
                         "labels": LABELS, "estimate": point.tolist(), "se": np.sqrt(cov.diagonal()).tolist(),
                         "full_covariance": cov.tolist(), "loo_vectors": loo.tolist(),
                         "jackknife_bias_estimate": ((b-1)*(loo.mean(axis=0)-point)).tolist(),
                         "exact_empirical_critical_intervals": intervals,
                         "all_LOO_curvature_types_stable": bool(np.all(np.array(curvature_signs) == np.sign(curvature)))}
    difference = means[400]-means[100]
    covariance = covariances[100]+covariances[400]
    blocks = {}
    for name, indices in (("rank_step_mu3_to_mu6_primary", np.arange(4)),
                          ("canonical_mu3_to_mu6_smoothing_control", np.arange(4, 8)),
                          ("canonical_ordered_peak_coordinates_auxiliary", np.arange(8, 14))):
        delta, cov = difference[indices], covariance[np.ix_(indices, indices)]
        q = float(delta @ np.linalg.solve(cov, delta))
        blocks[name] = {"labels": [LABELS[k] for k in indices], "chi2": q,
                        "df": len(indices), "nominal_p": float(chi2.sf(q, len(indices))),
                        "covariance_condition_number": float(np.linalg.cond(cov))}
    result = {"schema": protocol["schema"], "analysis_freeze": "0e805abb", "new_MC": 0,
              "sources": sources, "N400_minus_N100": {"labels": LABELS,
                  "estimate": difference.tolist(), "se": np.sqrt(covariance.diagonal()).tolist(),
                  "full_covariance": covariance.tolist(), "marginal_z": (difference/np.sqrt(covariance.diagonal())).tolist()},
              "frozen_block_comparisons": blocks,
              "boundary": "Fixed orders and ordinal landmarks, not selected windows. N100/N400 inputs had previously informed the width conjecture; this is a mechanism decomposition, not independent evidence from new samples. Signed-profile moments, not probability-law cumulants. The six peak coordinates exclude redundant height ratios from their omnibus score. Canonical peaks retain finite-N smoothing and are not literal extrema of the noisy rank-step function."}
    OUT.mkdir(exist_ok=True)
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Standardized rank shape beyond center and width", "",
             "| Frozen readout | N100 | N400 | N400-N100 | SE |", "|---|---:|---:|---:|---:|"]
    for i, name in enumerate(LABELS):
        lines.append(f"| {name} | {means[100][i]:.8g} | {means[400][i]:.8g} | {difference[i]:.8g} | {np.sqrt(covariance[i,i]):.5g} |")
    lines += ["", "## Fixed blocks", ""]
    for name, row in blocks.items():
        lines.append(f"- {name}: chi2={row['chi2']:.8g}/{row['df']}, nominal p={row['nominal_p']:.6g}.")
    lines += ["", result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
