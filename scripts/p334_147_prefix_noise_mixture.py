#!/usr/bin/env python3
"""Exact conditional-noise weight in the fixed uniform 147-prefix mixture.

Only two already-declared readouts: a canonical binomial tail and its integral.
No new p grid, prefix selection, network solve, continuation sampling, or tests.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = "results/p334-all147-prefix-clocks/full_clocks.json"
OUT = ROOT/"results/p334-147-prefix-noise-mixture"
P_REF, N, K0 = .59274605079, 425, 252


def frac(value):
    return {"exact": str(value), "float": float(value)}


def conditional(row, tails, clock):
    survival = [F(x["exact"]) for x in row["true_survival"]]
    weights_exact = [a-b for a, b in zip(survival, survival[1:])]
    assert sum(weights_exact) == 1 and all(x >= 0 for x in weights_exact)
    weights = np.array([float(x) for x in weights_exact])
    weights /= weights.sum()
    values = np.column_stack((tails, clock))
    shifted = values-values[0]
    mean_shift = weights@shifted
    mean = values[0]+mean_shift
    centered = shifted-mean_shift
    cov = (centered.T*weights)@centered
    exact_clock = [F(N+1-K0-t, N+1) for t in range(1, len(survival))]
    clock_mean = sum(w*x for w, x in zip(weights_exact, exact_clock))
    clock_second = sum(w*x*x for w, x in zip(weights_exact, exact_clock))
    return mean, cov, clock_mean, clock_second-clock_mean**2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, help="Completed local all147 full_clocks.json")
    parser.add_argument("--source-commit", required=True, help="Public commit containing the completed exact clocks")
    args = parser.parse_args()
    frozen = subprocess.check_output(["git", "show", args.source_commit+":"+SOURCE_PATH], cwd=ROOT)
    raw = args.input.read_bytes() if args.input else frozen
    if sha256(raw).digest() != sha256(frozen).digest():
        raise ValueError("Local full-clock file differs from declared commit")
    source = json.loads(raw)
    rows = sorted(source["records"], key=lambda row: row["counter"])
    assert len(rows) == 147 and len({row["counter"] for row in rows}) == 147
    assert all(row["status"] == "solved_full_physical" and int(row["original_row"]["k0"]) == K0 for row in rows)
    assert all(len(row["true_survival"]) == 174 for row in rows)
    t = np.arange(1, 174)
    tails = binom.sf(K0+t-1, N, P_REF)
    clock = (N+1-K0-t)/(N+1)
    means, covariances, exact_means, exact_variances, details = [], [], [], [], []
    for row in rows:
        mean, cov, emean, evar = conditional(row, tails, clock)
        means.append(mean)
        covariances.append(cov)
        exact_means.append(emean)
        exact_variances.append(evar)
        details.append({"counter": row["counter"], "conditional_mean": mean.tolist(),
                        "conditional_covariance": cov.tolist(),
                        "integrated_clock_mean": frac(emean), "integrated_clock_variance": frac(evar)})
    means, covariances = np.array(means), np.array(covariances)
    mean = means.mean(axis=0)
    within = covariances.mean(axis=0)
    centered = means-mean
    between = centered.T@centered/len(rows)  # Uniform finite mixture, ddof=0.
    total = within+between
    removed = within.diagonal()/total.diagonal()
    exact_mean = sum(exact_means)/len(rows)
    exact_within = sum(exact_variances)/len(rows)
    exact_between = sum((m-exact_mean)**2 for m in exact_means)/len(rows)
    exact_total = exact_within+exact_between
    result = {
        "schema": "matching-one/p334-uniform147-conditional-noise-mixture/v1",
        "source_commit": args.source_commit, "source_path": SOURCE_PATH, "source_sha256": sha256(raw).hexdigest(),
        "fixed_mixture": {"prefix_count": 147, "prefix_weights": "uniform 1/147", "variance_denominator": 147,
                          "N": N, "k0": K0, "orientation": "second", "age": 10, "ell": [12, -19],
                          "scope": "The complete supplied 147-member eligible old-source set, with the source exclusions retained; not a population law or all MC checkpoints."},
        "readouts": ["g_T(p_ref)=Pr[Binomial(425,p_ref)>=252+T]", "integral g_T=(174-T)/426"],
        "p_ref": P_REF, "mixture_mean": mean.tolist(),
        "within_prefix_suffix_covariance": within.tolist(),
        "between_prefix_conditional_mean_covariance": between.tolist(),
        "total_one_fresh_suffix_covariance": total.tolist(),
        "exact_conditional_averaging_removable_fraction": removed.tolist(),
        "after_exact_conditional_averaging_covariance": between.tolist(),
        "integrated_clock_exact": {"mean": frac(exact_mean), "within": frac(exact_within),
                                   "between": frac(exact_between), "total": frac(exact_total),
                                   "removable_fraction": frac(exact_within/exact_total)},
        "per_prefix": details,
        "sampling": {"new_MC": 0, "network_reruns": 0, "new_prefixes": 0,
                     "SE_or_population_CI": "not assigned; these are deterministic functionals of the fixed empirical mixture"},
        "boundary": "This is the conditional-stratum noise weight for one fresh uniform suffix per uniformly drawn member of these 147 fixed prefixes. It is not global production variance reduction, CPU speedup, an independent new random block, or a distribution over coupled orientations. Between-prefix variance remains after averaging."
    }
    OUT.mkdir(exist_ok=True)
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Conditional noise weight in the fixed 147-prefix mixture", "",
             "| Readout | Mean | Within-prefix suffix variance | Between-prefix variance | Total | Removable fraction |",
             "|---|---:|---:|---:|---:|---:|"]
    for i, name in enumerate(("canonical g(p_ref)", "integrated g")):
        lines.append(f"| {name} | {mean[i]:.10g} | {within[i,i]:.10g} | {between[i,i]:.10g} | {total[i,i]:.10g} | {100*removed[i]:.6g}% |")
    lines += ["", "The uniform finite-mixture variance uses denominator147, not a sample ddof correction. Baseline readouts are binomial tails, not Bernoulli events. All two-readout covariance is retained.", "", result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
