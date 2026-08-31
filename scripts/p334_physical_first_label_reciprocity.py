#!/usr/bin/env python3
"""Read two solved site-time tables by first/last reciprocity; no network calls."""
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/p334-physical-first-label-reciprocity"
SOURCE = "1c06230b"
RECIPROCITY = "31c17d48"
COUNTERS = (43042514269, 43042505280)


def packed(x):
    return {"exact": str(x), "value": float(x)}


def git_json(path, hashes):
    data = subprocess.check_output(["git", "show", f"{SOURCE}:{path}"], cwd=ROOT)
    hashes[path] = hashlib.sha256(data).hexdigest()
    return json.loads(data)


def score(clock, marked):
    counts = clock["true_safe_counts"]
    d = len(counts) - 1
    survival = [Fraction(x, comb(d, k)) for k, x in enumerate(counts)]
    mean = sum(survival)
    variance = sum((2*k+1)*x for k, x in enumerate(survival)) - mean**2
    site_rows, profiles, count_rows, means = [], [], [], []
    for site in marked["site_records"]:
        c = [0]
        for k in range(1, d+1):
            c.append(counts[k-1] - c[-1] - site["pivotal_count_by_prior_size"][k-1])
        if any(x < 0 or x > comb(d-1, k-1) for k, x in enumerate(c[1:], 1)):
            raise ValueError("Saved coefficients do not satisfy the safe-count convention")
        profile = [Fraction(1)] + [Fraction(c[k], comb(d-1, k-1)) for k in range(1, d+1)]
        m = sum(profile)
        means.append(m)
        profiles.append([float(x) for x in profile])
        count_rows.append(c)
        site_rows.append({"first_site": site["site"], "role": site["type"],
                          "final_winner_probability": site["birth_probability"],
                          "first_label_mean_T": packed(m),
                          "first_label_mean_child_wait": packed(m-1)})
    if len(means) != d or sum(means)/d != mean:
        raise ValueError("First-label means fail the uniform-label identity")
    if any(sum(c[k] for c in count_rows) != k*counts[k] for k in range(d+1)):
        raise ValueError("First-label curves do not recover the original clock")
    innovation = sum((x-mean)**2 for x in means)/d
    matrix = np.asarray(profiles)
    centered = matrix - np.asarray([float(x) for x in survival])[None, :]
    gram = centered.T @ centered / d
    singular = np.linalg.svd(centered/np.sqrt(d), compute_uv=False)
    eigenvalues = singular**2
    trace = float(eigenvalues.sum())
    energy = np.cumsum(eigenvalues)/trace
    distinct = len({tuple(row) for row in count_rows})
    role_rows = []
    role_variance = Fraction(0)
    for role in ("port_0", "port_1", "interior", "outside_core"):
        indices = [j for j, row in enumerate(site_rows) if row["role"] == role]
        if indices:
            role_mean = sum(means[j] for j in indices)/len(indices)
            role_variance += Fraction(len(indices), d)*(role_mean-mean)**2
            role_rows.append({"role": role, "sites": len(indices),
                              "mean_T_given_first_role": packed(role_mean),
                              "first_label_variance_contribution": packed(sum((means[j]-mean)**2 for j in indices)/d)})
    zero_winner = [j for j, row in enumerate(site_rows)
                   if Fraction(row["final_winner_probability"]["exact"]) == 0]
    h = sum(site["pivotal_count_by_prior_size"][0] for site in marked["site_records"])
    counter = marked["counter"]
    np.savez_compressed(OUT / f"first_label_arrays_{counter}.npz",
                        sites=np.array([r["first_site"] for r in site_rows]),
                        horizons=np.arange(d+1), conditional_survival=matrix,
                        mean_survival=np.asarray([float(x) for x in survival]),
                        covariance_Gram=gram, Gram_eigenvalues=eigenvalues)
    return {"counter": counter, "N": marked["N"], "k0": marked["k0"],
            "seed": marked["seed"], "remaining_labels": d,
            "clock_mean_T": packed(mean), "clock_variance_T": packed(variance),
            "original_direct_h": h, "binary_direct_safe_floor": 0 if h == 0 else None,
            "first_label_innovation_variance": packed(innovation),
            "first_label_fraction_of_clock_variance": packed(innovation/variance),
            "mean_remaining_variance_after_first_label": packed(variance-innovation),
            "first_role_fraction_of_first_label_innovation": packed(role_variance/innovation),
            "roles": role_rows, "zero_final_winner_sites": len(zero_winner),
            "mean_T_after_zero_winner_first": [packed(x) for x in sorted({means[j] for j in zero_winner})],
            "zero_winner_contribution_to_first_label_innovation": packed(sum((means[j]-mean)**2 for j in zero_winner)/d),
            "distinct_first_label_profiles": distinct,
            "Gram_numerical_singular_value_relative_cutoff": 1e-12,
            "Gram_numerically_resolved_rank": int(np.sum(singular > singular[0]*1e-12)),
            "Gram_trace": trace,
            "Gram_eigenvalues_descending": eigenvalues.tolist(),
            "Gram_leading_trace_fractions": (eigenvalues[:5]/trace).tolist(),
            "Gram_modes_for_99pct_trace": int(np.searchsorted(energy, .99)+1),
            "Gram_modes_for_999pct_trace": int(np.searchsorted(energy, .999)+1),
            "first_sites": site_rows}


def main():
    hashes = {}
    clocks = git_json("results/p334-contracted-full-clock/full_physical_birth_clock.json", hashes)
    indexed = {r["counter"]: r for r in clocks["records"]}
    OUT.mkdir(parents=True, exist_ok=True)
    records = [score(indexed[c], git_json(f"results/p334-exact-marked-birth/marked_birth_{c}.json", hashes))
               for c in COUNTERS]
    full_source = subprocess.check_output(["git", "rev-parse", SOURCE], cwd=ROOT, text=True).strip()
    result = {"schema": "matching-one/p334-physical-first-label-reciprocity/v1",
              "source_commit": full_source, "reciprocity_note_commit": RECIPROCITY,
              "source_sha256": hashes, "new_samples": 0, "new_network_or_DP_calls": 0,
              "dependence": "Two previously selected real N425 second-orientation prefixes; algebraic readout of the same complete final-site joint laws",
              "spectrum": "Numerical Gram spectrum only; no exact-rank claim",
              "records": records}
    result["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (OUT / "score.json").write_text(json.dumps(result, indent=2)+"\n")
    for row in records:
        print(row["counter"], "mean", row["clock_mean_T"]["value"],
              "variance", row["clock_variance_T"]["value"],
              "first_variance", row["first_label_innovation_variance"]["value"],
              "fraction", row["first_label_fraction_of_clock_variance"]["value"],
              "numerical_rank", row["Gram_numerically_resolved_rank"],
              "eigen_trace", row["Gram_leading_trace_fractions"])


if __name__ == "__main__":
    main()
