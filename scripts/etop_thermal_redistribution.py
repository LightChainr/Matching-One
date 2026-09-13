#!/usr/bin/env python3
"""Resolve the new N100 shape response along p, without further sampling.

The clock-calibrated odd residual has exactly zero integral.  Its nonzero
higher moments describe thermal redistribution, not another integral clock.
"""
from __future__ import annotations

import csv
from fractions import Fraction
import json
from math import comb
from pathlib import Path

import numpy as np
from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/etop-n100-three-modulus"
OUT = ROOT / "results/etop-thermal-redistribution"


def summary(values, influences):
    covariance = np.cov(influences, rowvar=False, ddof=1) / len(influences)
    se = np.sqrt(np.maximum(np.diag(np.atleast_2d(covariance)), 0))
    return {"mean": values.tolist(), "se": se.tolist(),
            "covariance": np.atleast_2d(covariance).tolist()}


def main():
    source = json.loads((SOURCE / "score.json").read_text())
    contract = source["contract"]
    n, b = contract["area"], contract["batches"]
    count = contract["samples_per_shape_pair"] // b
    # Shape x batch x event(K1,K2) x threshold-rank, already P4 normalized.
    counts = np.zeros((3, b, 2, n+1))
    for s, shape in enumerate(contract["shapes"]):
        delta = float(Fraction(shape["delta_cos4"]))
        with (SOURCE / "raw" / (shape["name"] + ".hist.csv")).open() as stream:
            for row in csv.DictReader(stream):
                sign = 1 if row["orientation"] == "first" else -1
                counts[s, int(row["batch"]), int(row["kind"] == "plus"), int(row["k"])] += sign * int(row["count"]) / (count * delta)
    # Constants in A=F1+F2-1 and E=1+F2-F1 cancel in the orientation contrast.
    fields = np.stack((counts[:, :, 0] + counts[:, :, 1],
                       counts[:, :, 1] - counts[:, :, 0]), axis=2)
    d, u = fields[1]-fields[0], fields[2]-fields[0]
    dm, um = d.mean(axis=0), u.mean(axis=0)
    clocks = np.array(source["batch_vectors"])[:, [2, 6, 10]]
    dc, uc = clocks[:, 1]-clocks[:, 0], clocks[:, 2]-clocks[:, 0]
    r = uc.mean()/dc.mean()
    # First-order influence of the ratio, including all same-stream covariance.
    dr = (uc-r*dc)/dc.mean()
    cm = um-r*dm
    ci = (u-um)-r*(d-dm)-dr[:, None, None]*dm
    re4 = source["same_area_models"]["affine_E4"]["secant_ratio"]
    em = um-re4*dm
    ei = (u-um)-re4*(d-dm)

    # Exact integration of each Bernstein tail against z^j, z=N^(3/8)(p-p0).
    # int p^j F_k(p) dp = (1-E[Beta(k,N+1-k)^(j+1)])/(j+1).
    # The 1/(j+1) term cancels since every contrast has zero coefficient sum.
    pc, scale = contract["fixed_p"], n**(3/8)
    k = np.arange(n+1, dtype=float)
    rising = np.ones(n+1)
    beta_moments = []
    for j in range(1, 6):
        rising = rising*(k+j-1)/(n+j)
        beta_moments.append(rising.copy())
    kernels = []
    for j in range(5):
        kernel = sum(comb(j, m)*(-pc)**(j-m)*(-beta_moments[m]/(m+1))
                     for m in range(j+1)) * scale**j
        kernels.append(kernel)
    kernels = np.array(kernels).T
    moment_mean, moment_influence = cm @ kernels, ci @ kernels
    # Zero odd area is an algebraic constraint, not a numerically measured z.
    moment_mean[0, 0] = 0.
    moment_influence[:, 0, 0] = 0.

    grid = np.linspace(0., 1., 801)
    tails = binom.sf(k[:, None]-1, n, grid[None, :])
    profiles = cm @ tails
    influence_profiles = ci @ tails
    fixed_profiles = em @ tails
    fixed_influences = ei @ tails
    # Full covariance at a compact grid; dense points for descriptive lobe roots.
    keep = np.arange(0, len(grid), 10)
    curve_result = summary(profiles[:, keep].reshape(-1),
                           influence_profiles[:, :, keep].reshape(b, -1))
    curve_result.update({"p": grid[keep].tolist(), "layout": "A(p-grid), E(p-grid)"})
    fixed_result = summary(fixed_profiles[:, keep].reshape(-1),
                           fixed_influences[:, :, keep].reshape(b, -1))
    fixed_result.update({"p": grid[keep].tolist(), "layout": "A(p-grid), E(p-grid)"})
    descriptions = []
    for j, name in enumerate(("A_top", "E_top")):
        se = influence_profiles[:, j].std(axis=0, ddof=1)/np.sqrt(b)
        peak = int(np.argmax(np.abs(profiles[j])))
        crosses = np.flatnonzero(profiles[j, :-1]*profiles[j, 1:] < 0)
        # Retain core crossings; tiny-tail oscillations have no physical weight.
        roots = [float(grid[t]-profiles[j, t]*(grid[t+1]-grid[t]) /
                       (profiles[j, t+1]-profiles[j, t]))
                 for t in crosses if .25 <= grid[t] <= .85]
        descriptions.append({"field": name, "core_empirical_zeroes": roots,
            "largest_mean_lobe_p": float(grid[peak]), "largest_mean_lobe": float(profiles[j, peak]),
            "pointwise_se_at_selected_peak": float(se[peak]),
            "positive_area_quadrature": float(np.trapezoid(np.maximum(profiles[j], 0), grid)),
            "negative_area_quadrature": float(np.trapezoid(np.minimum(profiles[j], 0), grid)),
            "scope": "Post-reveal descriptive roots and selected peak; not simultaneous significance."})
    result = {
        "source": "results/etop-n100-three-modulus/score.json",
        "dependency_group": "N100 seed 20260831125401 offset 267100000000; all 3 shapes and all p reuse one block",
        "status": "post-reveal exploratory thermal decomposition; zero additional MC",
        "definition": "D=Y(4i)-Y(2i), U=Y(1/2+i)-Y(2i), R(p)=U(p)-r_C D(p)",
        "clock_secant": float(r), "clock_secant_se": float(dr.std(ddof=1)/np.sqrt(b)),
        "identities": ["integral P4[A_top(p)] dp = -2 C",
                       "integral P4[E_top(p)] dp = -W",
                       "integral R_A(p) dp = 0 exactly by clock calibration",
                       "integral R_E(p) dp = -R_W"],
        "moments": {"definition": "integral [N^(3/8)*(p-p_ref)]^j R(p) dp, j=0..4",
                    "layout": "A(j0..4), E(j0..4)",
                    **summary(moment_mean.reshape(-1), moment_influence.reshape(b, -1))},
        "clock_quotient_profiles": curve_result, "fixed_E4_profiles": fixed_result,
        "empirical_lobes": descriptions,
        "boundary": "Clock calibration is data-adaptive; uncertainties use the correlated ratio influence. No new field count, exact E4 identity, independent p-grid evidence, or continuum exponent is inferred."}
    OUT.mkdir(exist_ok=True)
    (OUT / "thermal_profiles.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Clock-compatible integral, nontrivial thermal redistribution", "",
        "The new N100 stream contains more information than its fixed-p readout. The exact integral of the odd orientation contrast is -2C; that of the even contrast is -W. Therefore using C as the shape coordinate forces the remaining odd response to have zero area, but does not force its p-profile to vanish.", "",
        f"The fitted clock secant is {r:.9g} +/- {result['clock_secant_se']:.4g}. This is a post-reveal decomposition of the existing correlated block, with no new sampling.", "",
        "## Thermal moments of the clock-quotient residual", "",
        "Here z=N^(3/8)(p-p_ref). Integrals are evaluated algebraically from threshold-rank histograms, not by numerical quadrature. The zero odd area is imposed by an exact identity, not a failed-to-reject test.", "",
        "| field | power of z | signed integral | SE |", "|---|---:|---:|---:|"]
    errors = np.array(result["moments"]["se"]).reshape(2, 5)
    for j, name in enumerate(("A_top", "E_top")):
        for power in range(5):
            lines.append(f"| {name} | {power} | {moment_mean[j,power]:.9g} | {errors[j,power]:.5g} |")
    lines += ["", "## Empirical profile structure", ""]
    for row in descriptions:
        lines.append(f"- {row['field']}: core zeroes {row['core_empirical_zeroes']}; largest mean lobe {row['largest_mean_lobe']:.6g} near p={row['largest_mean_lobe_p']:.5g}; positive/negative areas {row['positive_area_quadrature']:.6g}/{row['negative_area_quadrature']:.6g}.")
    lines += ["", "Roots, quadrature lobe areas and peak locations are descriptive. All pointwise readouts share one dependency group; the JSON retains their full compact-grid covariance. A separately propagated fixed-E4 profile is included for comparison.", "", result["boundary"], ""]
    (OUT / "REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
