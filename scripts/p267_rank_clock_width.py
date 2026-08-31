#!/usr/bin/env python3
"""Exact finite-N moment decomposition, evaluated on two existing rank archives.

The microcanonical signed step profile has value f_j on [j/N,(j+1)/N).
Its canonical readout is sum_j f_j B_(N,j)(p). No inversion or new MC.
"""
import json
from pathlib import Path

import numpy as np

from p267_scalar_clock_transport import load_source

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"results/p267-rank-clock-width"
SOURCES = {100: "7b30648be558df0652a7ff22143cc87ed399d042",
           400: "3e01b495b5b637b0070705e37b4137a9a0ef0d8b"}
LABELS = ["canonical_variance_z", "rank_step_variance_z", "canonical_minus_rank_step_z",
          "rank_grid_variance_z", "beta_kernel_term_z", "contracted_rank_grid_term_z",
          "rank_step_variance_x_N_quarter"]


def width_values(f):
    n, j = len(f)-1, np.arange(len(f))
    weights = f/f.sum()
    mean_j = weights @ j
    variance_j = weights @ ((j-mean_j)**2)
    scale_squared = n**.75
    beta_term = (weights @ ((j+1)*(n-j+1)))/((n+2)**2*(n+3))*scale_squared
    contracted = variance_j/(n+2)**2*scale_squared
    canonical = contracted+beta_term
    rank_step = (variance_j+1/12)/n**2*scale_squared
    rank_grid = variance_j/n**2*scale_squared
    return np.array([canonical, rank_step, canonical-rank_step, rank_grid,
                     beta_term, contracted, rank_step/n**.25])


def jackknife_covariance(values):
    b = len(values)
    centered = values-values.mean(axis=0)
    return (b-1)/b*(centered.T @ centered)


def main():
    base = json.loads((ROOT/"experiments/p267_scalar_clock_transport_20260831.json").read_text())
    rows, estimates, covariances, replicas = {}, {}, {}, {}
    for n, commit in SOURCES.items():
        protocol = dict(base, source_commit=commit, source_directory=f"results/etop-n{n}-three-modulus")
        contract, hashes, bernstein, _ = load_source(protocol)
        batches = bernstein[0, :, 0]  # D_A only; no new selection among channels.
        b, mean = len(batches), batches.mean(axis=0)
        point = width_values(mean)
        loo = np.array([width_values((b*mean-f)/(b-1)) for f in batches])
        cov = jackknife_covariance(loo)
        estimates[n], covariances[n], replicas[n] = point, cov, loo
        rows[str(n)] = {"source_commit": commit, "source_sha256": hashes,
                        "source_contract": contract, "labels": LABELS,
                        "estimate": point.tolist(), "se": np.sqrt(cov.diagonal()).tolist(),
                        "full_covariance": cov.tolist(), "loo_vectors": loo.tolist(),
                        "signed_weight_boundary": "Weights are the normalized signed microcanonical profile, not the probability law of K1 or K2. Moment identities do not require all weights to be positive."}
    delta = estimates[400]-estimates[100]
    delta_cov = covariances[400]+covariances[100]
    def fraction(change):
        return np.array([change[1]/change[0], change[2]/change[0]])
    fraction_point = fraction(delta)
    fraction_cov = (jackknife_covariance(np.array([fraction(estimates[400]-x) for x in replicas[100]]))
                    +jackknife_covariance(np.array([fraction(x-estimates[100]) for x in replicas[400]])))
    def working_fingerprint(old, new):
        gamma = .375-np.log(new[1]/old[1])/(2*np.log(4))
        return np.array([gamma, new[1]*(900/400)**.25, new[1]])
    fingerprint = working_fingerprint(estimates[100], estimates[400])
    fingerprint_cov = (
        jackknife_covariance(np.array([working_fingerprint(x, estimates[400]) for x in replicas[100]]))
        +jackknife_covariance(np.array([working_fingerprint(estimates[100], x) for x in replicas[400]])))
    result = {"schema": "matching-one/p267-canonical-rank-clock-width/v1",
              "status": "post-reveal physical decomposition of existing independent N100/N400 archives; no new MC",
              "observable": "D_A=P4[A_top](4i)-P4[A_top](2i)",
              "coordinates": {"z": "N^(3/8)*(p-p_ref)", "x": "N^(1/4)*(p-p_ref), post-reveal working fingerprint only"},
              "sources": rows,
              "N400_minus_N100": {"labels": LABELS, "estimate": delta.tolist(),
                  "se": np.sqrt(delta_cov.diagonal()).tolist(), "covariance": delta_cov.tolist()},
              "fraction_of_canonical_width_increase": {"labels": ["retained_in_rank_step", "net_canonicalization_correction"],
                  "estimate": fraction_point.tolist(), "se": np.sqrt(fraction_cov.diagonal()).tolist(),
                  "covariance": fraction_cov.tolist()},
              "post_reveal_working_fingerprint": {
                  "labels": ["two_area_gamma_effective", "N900_rank_variance_z_if_width_N_minus_quarter", "N900_rank_variance_z_if_fixed_critical_width_profile"],
                  "estimate": fingerprint.tolist(), "se": np.sqrt(fingerprint_cov.diagonal()).tolist(),
                  "covariance": fingerprint_cov.tolist(),
                  "definition": "gamma_eff=3/8-log(Vz400/Vz100)/(2 log 4); predictions use the same N400 anchor",
                  "boundary": "A two-size signed-profile effective width, selected after reveal; not a new critical exponent or proof of full-profile collapse. N900 predictions are conditional working targets, not measurements."
              },
              "exact_identities": [
                  "C(p)=sum_j f_j BinomialPMF(j;N,p), S(p)=f_floor(Np), f_N=0",
                  "int C=sum f/(N+1), int S=sum f/N",
                  "mean_C=(mean_J+1)/(N+2), mean_S=(mean_J+1/2)/N",
                  "var_C=var_J/(N+2)^2+E_w[(J+1)(N-J+1)]/[(N+2)^2(N+3)]",
                  "var_S=(var_J+1/12)/N^2"
              ],
              "boundary": "Central moments of signed profiles, not a new threshold estimator or probability distribution. The N^(1/4) coordinate was suggested after both scales were seen; two areas do not prove an exponent or full profile collapse. All reuse remains inside the original two dependency groups."}
    OUT.mkdir(exist_ok=True)
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Width broadening remains in the rank-clock profile", "",
             "| Central z variance/component | N100 | N400 |", "|---|---:|---:|"]
    for i, label in enumerate(LABELS):
        lines.append(f"| {label} | {estimates[100][i]:.9g} +/- {np.sqrt(covariances[100][i,i]):.6g} | {estimates[400][i]:.9g} +/- {np.sqrt(covariances[400][i,i]):.6g} |")
    lines += ["", f"Rank-step width increase: {delta[1]:.9g} +/- {np.sqrt(delta_cov[1,1]):.6g}.",
              f"Fraction of canonical width increase retained in rank-step profile: {fraction_point[0]:.9g} +/- {np.sqrt(fraction_cov[0,0]):.6g}.", "",
              "The conditional-Beta term decreases with N, while the deterministic grid contraction also weakens. Canonical minus rank-step is a net correction, not a pure positive extra-noise variance.", "", result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
