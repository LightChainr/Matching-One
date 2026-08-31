#!/usr/bin/env python3
"""Turn the two exact physical T laws into conditional canonical K2 readouts.

Reads the solved archive; does not rerun reliability DP or sample prefixes.
"""
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "6358ba49ef390c10a3f501b589ba7ba1d4e05b09"
SOURCE_PATH = "results/p334-contracted-full-clock/full_physical_birth_clock.json"
OUT = ROOT/"results/p334-conditional-thermal-averaging"
P_REF = .59274605079
GRID = np.array([0., .4, .5, .525, .55, .57, .585, P_REF, .6,
                 .615, .63, .65, .675, .7, .75, .8, .9, 1.])


def rational(value):
    return {"exact": str(value), "float": float(value)}


def weighted_covariance(values, probabilities):
    probabilities = np.asarray(probabilities, dtype=float)
    probabilities /= probabilities.sum()
    # Anchor subtraction makes constant p=0/1 columns exactly zero variance.
    shifted = values-values[0]
    mean_shift = probabilities @ shifted
    centered = shifted-mean_shift
    mean = values[0]+mean_shift
    return mean, (centered.T*probabilities) @ centered


def score_record(row, grid=GRID):
    survival = [Fraction(s["exact"]) for s in row["true_survival"]]
    t_values = list(range(1, len(survival)))
    probabilities = [survival[t-1]-survival[t] for t in t_values]
    assert sum(probabilities) == 1 and all(p >= 0 for p in probabilities)
    mean_t = sum(Fraction(t)*p for t, p in zip(t_values, probabilities))
    variance_t = sum((Fraction(t)-mean_t)**2*p for t, p in zip(t_values, probabilities))
    assert mean_t == Fraction(row["mean_true_birth_step"]["exact"])
    n, k0 = row["N"], row["k0"]
    ranks = k0+np.array(t_values)
    weights = np.array(list(map(float, probabilities)))
    g = binom.sf(ranks[:, None]-1, n, grid[None, :])
    area = (n+1-ranks)/(n+1)
    values = np.column_stack((g, area))
    mean, covariance = weighted_covariance(values, weights)
    p_index = int(np.flatnonzero(grid == P_REF)[0])
    exact_area = (Fraction(n+1-k0)-mean_t)/(n+1)
    exact_area_var = variance_t/(n+1)**2
    assert abs(mean[-1]-float(exact_area)) < 1e-14
    assert abs(covariance[-1, -1]-float(exact_area_var)) < 1e-14
    law = [{"T": t, "probability": str(p)} for t, p in zip(t_values, probabilities) if p]
    return {
        "counter": row["counter"], "seed": row["seed"], "N": n, "k0": k0,
        "conditional_T_law": law, "mean_T": rational(mean_t), "variance_T": rational(variance_t),
        "source_K1": "not reconstructed; source record gives k0 and K2 waiting law, not an explicit K1",
        "layout": [f"F2(p={p:.11g})" for p in grid]+["integral_0^1_F2_dp"],
        "p_grid": grid.tolist(), "conditional_mean": mean.tolist(),
        "suffix_conditional_covariance_removed": covariance.tolist(),
        "p_ref": {"p": P_REF, "conditional_F2": float(mean[p_index]),
                  "continuation_E_g_squared": float(covariance[p_index, p_index]+mean[p_index]**2),
                  "continuation_variance": float(covariance[p_index, p_index]),
                  "continuation_sd": float(np.sqrt(covariance[p_index, p_index])),
                  "incorrect_Bernoulli_variance_not_used": float(mean[p_index]*(1-mean[p_index])),
                  "covariance_with_integrated_clock": float(covariance[p_index, -1])},
        "integrated_clock": {"conditional_mean": rational(exact_area),
                             "continuation_variance": rational(exact_area_var),
                             "continuation_sd": float(np.sqrt(float(exact_area_var)))},
        "after_exact_conditional_averaging": {"remaining_suffix_variance_p_ref": 0,
                                               "remaining_suffix_variance_integrated_clock": 0,
                                               "remaining_between_prefix_variance": "not identified by two selected prefixes"},
        "A_E_second_birth_contribution": {"loading": [1, 1],
            "isolated_K2_covariance_at_p_ref": (covariance[p_index, p_index]*np.ones((2, 2))).tolist(),
            "scope": "Only the +F2 term in A=F1+F2-1 and E=1+F2-F1. Full A/E means or covariance are not reconstructed without K1 semantics."},
        "numerics": "T probabilities and clock moments are exact rationals; binomial tails and their weighted thermal covariance are evaluated in double precision. No Monte Carlo error is attached to this deterministic conditional-law calculation."
    }


def main():
    raw = subprocess.check_output(["git", "show", SOURCE_COMMIT+":"+SOURCE_PATH], cwd=ROOT)
    source = json.loads(raw)
    rows = [score_record(row) for row in source["records"]]
    a, b = rows
    area_difference = (Fraction(b["integrated_clock"]["conditional_mean"]["exact"])
                       -Fraction(a["integrated_clock"]["conditional_mean"]["exact"]))
    result = {
        "schema": "matching-one/p334-exact-conditional-K2-thermal-readout/v1",
        "source_commit": SOURCE_COMMIT, "source_path": SOURCE_PATH, "source_sha256": sha256(raw).hexdigest(),
        "new_prefixes": 0, "new_continuations": 0, "reliability_DP_rerun": False,
        "definition": "K2=k0+T; g(K2,p)=Pr[Binomial(N,p)>=K2]; F2_X(p)=E[g|prefix X]",
        "records": rows,
        "between_prefix_mean_comparison": {
            "B_minus_A_F2_at_p_ref": b["p_ref"]["conditional_F2"]-a["p_ref"]["conditional_F2"],
            "B_minus_A_integrated_clock": rational(area_difference),
            "ordering": "The source proves T_B stochastically later than T_A, hence F2_B(p)<=F2_A(p) for every p. This is inherited, not inferred from grid points.",
            "cross_prefix_covariance": "not specified: marginal T laws do not define a coupled A/B suffix stream"},
        "noise_decomposition": "Cov(G)=Cov(E[G|prefix])+E[Cov(G|prefix)]. Exact conditional averaging removes the second term for each averaged prefix; these two selected prefixes do not estimate its ensemble average.",
        "boundary": "A conditional noise-removal calculation, not total production variance reduction or CPU speedup. Prefix frequencies, between-prefix variance, coupled orientations, and algorithmic cost across the production ensemble are not determined here."
    }
    OUT.mkdir(exist_ok=True)
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Exact conditional averaging of the true second-birth clock", "",
             "| Prefix | F2(p_ref) | continuation Var[g(p_ref)] | integrated F2 | continuation Var[integrated F2] |",
             "|---|---:|---:|---:|---:|"]
    for row in rows:
        lines.append(f"| {row['counter']} | {row['p_ref']['conditional_F2']:.10g} | {row['p_ref']['continuation_variance']:.10g} | {row['integrated_clock']['conditional_mean']['float']:.10g} | {row['integrated_clock']['continuation_variance']['float']:.10g} |")
    lines += ["", "These are variances of the already-canonical binomial-tail readout, not F2*(1-F2). Exact conditional averaging makes the suffix contribution deterministic for each fixed prefix. It does not remove between-prefix variation.", "", result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
