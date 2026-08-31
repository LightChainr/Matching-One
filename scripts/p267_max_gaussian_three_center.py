#!/usr/bin/env python3
"""One frozen positive moment model, then unused moments seven and eight."""
import argparse
import json
from math import comb, factorial
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.stats import chi2

from p267_scalar_clock_transport import load_source

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT/"experiments/p267_max_gaussian_three_center_20260831.json"
OUT = ROOT/"results/p267-max-gaussian-three-center"
LABELS = (["gaussian_variance_fraction"]
          +[f"center_u_{i}" for i in range(3)]
          +[f"weight_{i}" for i in range(3)]
          +["profile_mean_p", "profile_variance_p", "gaussian_variance_p"]
          +[f"center_p_{i}" for i in range(3)]
          +["gaussian_variance_z", "gaussian_variance_x_quarter"]
          +[f"observed_m{r}" for r in (7, 8)]
          +[f"predicted_m{r}" for r in (7, 8)]
          +[f"observed_minus_predicted_m{r}" for r in (7, 8)])


def rank_step_moments(f):
    """Exact bin integration, evaluated in centered coordinates for stability."""
    n, j = len(f)-1, np.arange(len(f), dtype=float)
    w = f/f.sum()
    jm = w@j
    variance_j = w@((j-jm)**2)+1/12
    if variance_j <= 0:
        raise ValueError("Nonpositive rank-step signed variance")
    sigma_j = np.sqrt(variance_j)
    u = (j-jm)/sigma_j
    moments = []
    for r in range(9):
        # Integrate (u+Uniform[-1/(2*sigma_j),1/(2*sigma_j)])**r.
        integrated = sum(comb(r, k)*u**(r-k)/(2**k*(k+1)*sigma_j**k)
                         for k in range(0, r+1, 2))
        moments.append(float(w@integrated))
    moments[:3] = [1., 0., 1.]
    return np.array(moments), (jm+.5)/n, variance_j/n**2, {
        "normalized_negative_bin_mass": float(-w[w < 0].sum()),
        "signed_area": float(f.sum()/n), "nonzero_rank_bins": int(np.count_nonzero(f))}


def heat_moments(moments, variance):
    """Positive variance convolves, negative variance formally deconvolves."""
    return np.array([sum(variance**j*factorial(r)/(2**j*factorial(j)*factorial(r-2*j))*moments[r-2*j]
                         for j in range(r//2+1)) for r in range(len(moments))])


def hankel(moments, size):
    return moments[np.add.outer(np.arange(size), np.arange(size))]


def realization(moments):
    """First H3 PSD boundary; no second family or component count is tried."""
    initial = np.linalg.eigvalsh(hankel(moments[:7], 4))
    if initial[0] <= 1e-10:
        raise ValueError("H3(0) is not positive definite; frozen interior construction not available")
    def minimum(t):
        return np.linalg.eigvalsh(hankel(heat_moments(moments[:7], -t), 4))[0]
    if minimum(1.) >= 0:
        raise ValueError("No nondegenerate H3 boundary before zero residual variance")
    t = brentq(minimum, 0., 1., xtol=2e-14)
    q = heat_moments(moments[:7], -t)
    h2, h3 = hankel(q, 3), hankel(q, 4)
    ev2, ev3 = np.linalg.eigvalsh(h2), np.linalg.eigvalsh(h3)
    if ev2[0] <= 1e-10 or abs(ev3[0]) >= 1e-9:
        raise ValueError("H2 positivity / flat H3 boundary gate failed")
    shifted = q[np.add.outer(np.arange(3), np.arange(3))+1]
    centers = eigh(shifted, h2, eigvals_only=True)
    weights = np.linalg.solve(np.vander(centers, 3, increasing=True).T, q[:3])
    if np.min(weights) <= 0 or np.min(np.diff(centers)) <= 1e-10:
        raise ValueError("No three distinct positive-weight centers")
    atom_moments = np.array([weights@centers**r for r in range(9)])
    predictions = heat_moments(atom_moments, t)
    error = float(np.max(np.abs(predictions[:7]-moments[:7])))
    if error >= 1e-8:
        raise ValueError("Three-center model does not reproduce all construction moments")
    cubic = np.r_[np.linalg.solve(h2, -q[3:6]), 1.]
    certificate = {"H3_at_zero_eigenvalues": initial.tolist(),
                   "deconvolved_moments_0_to_6": q.tolist(),
                   "H2_boundary_eigenvalues": ev2.tolist(),
                   "H3_boundary_eigenvalues": ev3.tolist(),
                   "H3_boundary_numerical_rank": int(np.count_nonzero(ev3 > 1e-9)),
                   "orthogonal_cubic_ascending_coefficients": cubic.tolist(),
                   "orthogonal_cubic_at_centers": np.polynomial.polynomial.polyval(centers, cubic).tolist(),
                   "max_construction_moment_error": error,
                   "interpretation": "Floating-point realization of empirical moments, not an exact population positivity certificate"}
    return t, centers, weights, predictions, certificate


def features_from_moments(moments, mean, variance, n):
    """Direct same-definition entry: standardized m0..m8, mean_p, var_p, N."""
    moments = np.asarray(moments, dtype=float)
    if moments.shape != (9,) or not np.allclose(moments[:3], [1, 0, 1], atol=1e-10):
        raise ValueError("Require the standardized rank-step moments m0..m8")
    t, centers, weights, predicted, certificate = realization(moments)
    vector = np.r_[t, centers, weights, mean, variance, t*variance,
                   mean+np.sqrt(variance)*centers, t*variance*n**.75,
                   t*variance*n**.5, moments[7:9], predicted[7:9],
                   moments[7:9]-predicted[7:9]]
    return vector, {"moments_0_to_8": moments.tolist(),
                    "realization": certificate}


def features(f):
    moments, mean, variance, source_shape = rank_step_moments(f)
    vector, details = features_from_moments(moments, mean, variance, len(f)-1)
    details["source_shape"] = source_shape
    return vector, details


def affine_shape_chart(vector):
    """Four nonredundant shape coordinates on the positive flat-rank3 branch."""
    t, left, middle, right = vector[:4]
    return np.array([t, vector[4], vector[6], (middle-left)/(right-left)])


def jackknife_covariance(loo):
    b = len(loo)
    centered = loo-loo.mean(axis=0)
    return (b-1)/b*centered.T@centered


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-commit", help="Score exactly one new completed committed archive")
    parser.add_argument("--source-directory", help="Git-relative archive path, with the original score.json and raw/*.hist.csv schema")
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    if bool(args.source_commit) != bool(args.source_directory):
        parser.error("--source-commit and --source-directory must be given together")
    protocol = json.loads(PROTOCOL.read_text())
    inputs = ([(args.source_commit, args.source_directory)] if args.source_commit else
              [(commit, f"results/etop-n{area}-three-modulus")
               for area, commit in protocol["source_commits"].items()])
    sources, estimates, covariances = {}, {}, {}
    for commit, source_directory in inputs:
        contract, hashes, bernstein, _ = load_source({
            "source_commit": commit, "source_directory": source_directory})
        n = contract["area"]
        area = str(n)
        batches = bernstein[0, :, 0]  # The same D_A as the rank-width and shape analyses.
        b, mean = len(batches), batches.mean(axis=0)
        row = {"source_commit": commit, "source_directory": source_directory,
               "source_sha256": hashes, "source_contract": contract,
               "common_batches": b, "labels": LABELS}
        try:
            point, details = features(mean)
        except ValueError as error:
            row.update(status="construction_not_available", reason=str(error))
            sources[area] = row
            continue
        loo, failures, max_error, min_weight, min_h2 = [], [], 0., 1., 1.
        for k, f in enumerate(batches):
            try:
                vector, detail = features((b*mean-f)/(b-1))
            except ValueError as error:
                failures.append({"batch": k, "reason": str(error)})
                continue
            loo.append(vector)
            max_error = max(max_error, detail["realization"]["max_construction_moment_error"])
            min_weight = min(min_weight, np.min(vector[4:7]))
            min_h2 = min(min_h2, detail["realization"]["H2_boundary_eigenvalues"][0])
        row.update(status="positive_three_center_realization", estimate=point.tolist(), **details)
        row["LOO_construction"] = {"successful": len(loo), "failed": failures,
                                   "max_m0_to_m6_error": max_error,
                                   "minimum_weight": float(min_weight), "minimum_H2_eigenvalue": min_h2}
        if failures:
            row["uncertainty_status"] = "not_reported: some LOO constructions failed; no family changed"
            sources[area] = row
            continue
        loo = np.array(loo)
        cov = jackknife_covariance(loo)
        chart_loo = np.array([affine_shape_chart(v) for v in loo])
        chart_cov = jackknife_covariance(chart_loo)
        row["affine_shape_chart"] = {
            "labels": ["gaussian_variance_fraction", "left_weight", "right_weight", "relative_middle_gap"],
            "estimate": affine_shape_chart(point).tolist(),
            "se": np.sqrt(np.diag(chart_cov)).tolist(), "full_covariance": chart_cov.tolist(),
            "scope": "A source-defined change of coordinates; not an added omnibus comparison"}
        estimates[n], covariances[n] = point, cov
        delta, residual_cov = point[-2:], cov[-2:, -2:]
        q = float(delta@np.linalg.solve(residual_cov, delta))
        # Pure algebraic decomposition of the already reported quarter-coordinate
        # width, not another fitted model or added comparison family.
        transform = np.zeros((3, len(point)))
        transform[0, 8] = np.sqrt(n)
        transform[1, 14] = 1.
        transform[2] = transform[0]-transform[1]
        width_cov = transform@cov@transform.T
        row["quarter_coordinate_width_decomposition"] = {
            "labels": ["total", "common_gaussian", "between_centers"],
            "estimate": (transform@point).tolist(),
            "se": np.sqrt(np.diag(width_cov)).tolist(),
            "full_covariance": width_cov.tolist(),
            "scope": "An identity within the fitted moment representation, not an independent shape test"}
        row.update(se=np.sqrt(np.diag(cov)).tolist(), full_covariance=cov.tolist(),
                   leave_one_common_batch_out_vectors=loo.tolist(),
                   jackknife_bias_estimate=((b-1)*(loo.mean(axis=0)-point)).tolist(),
                   unused_moment_readout={"orders": [7, 8], "residual": delta.tolist(),
                       "covariance": residual_cov.tolist(), "marginal_z": (delta/np.sqrt(np.diag(residual_cov))).tolist(),
                       "chi2": q, "df": 2, "nominal_p": float(chi2.sf(q, 2)),
                       "scope": "Same-source post-reveal algebraic prediction check, not independent data"})
        sources[area] = row
    result = {"schema": protocol["schema"], "construction_freeze": "191c20e2", "new_MC": 0,
              "protocol": protocol, "sources": sources,
              "boundary": protocol["boundary"]}
    if args.source_commit:
        result["application_stage"] = "One completed-source auxiliary application of the N100/N400-defined construction; no model reselection"
        result["boundary"] = (protocol["boundary"].replace(" N900 is not read.", "")
                              + " The historical protocol records the construction-stage boundary; only the explicitly supplied completed archive is read in this application. Its original primary target is not changed.")
    if set(estimates) == {100, 400}:
        delta = estimates[400]-estimates[100]
        cov = covariances[400]+covariances[100]
        result["independent_source_scale_difference"] = {
            "labels": LABELS, "N400_minus_N100": delta.tolist(),
            "se": np.sqrt(np.diag(cov)).tolist(), "full_covariance": cov.tolist(),
            "scope": "Parameter changes are descriptive; no separate omnibus scan or field counting"}
        old = sources["100"]["quarter_coordinate_width_decomposition"]
        new = sources["400"]["quarter_coordinate_width_decomposition"]
        width_cov = np.array(old["full_covariance"])+np.array(new["full_covariance"])
        result["quarter_coordinate_width_decomposition_difference"] = {
            "labels": old["labels"], "N400_minus_N100": (np.array(new["estimate"])-old["estimate"]).tolist(),
            "se": np.sqrt(np.diag(width_cov)).tolist(), "full_covariance": width_cov.tolist()}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Maximal common Gaussian variance and three positive centers", "",
             "Post-reveal construction from m0..m6; m7/m8 are unused same-source readouts.", ""]
    for area, row in sources.items():
        lines += [f"## N{area}", "", row["status"], ""]
        if "estimate" not in row:
            lines += [row["reason"], ""]
            continue
        lines += ["| Parameter/readout | Estimate | Batch SE |", "|---|---:|---:|"]
        for name, value, se in zip(LABELS, row["estimate"], row.get("se", [None]*len(LABELS))):
            lines.append(f"| {name} | {value:.10g} | {se:.7g} |" if se is not None else f"| {name} | {value:.10g} | not resolved |")
        if "unused_moment_readout" in row:
            readout = row["unused_moment_readout"]
            lines += ["", f"Unused-moment residual: chi2={readout['chi2']:.9g}/2; nominal p={readout['nominal_p']:.8g}.", ""]
    lines += [result["boundary"], ""]
    (args.output/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
