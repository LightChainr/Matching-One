#!/usr/bin/env python3
"""Use the new N100 stream to distinguish low-order thermal transport.

This is post-reveal mechanism analysis, with no new sampling. Moment anchors
fix a polynomial transport velocity; higher moments are overidentifications,
not statistically independent holdouts. All source uncertainty is jackknifed.
"""
from __future__ import annotations

import csv
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path

import numpy as np
from scipy.stats import binom, chi2

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/etop-n100-three-modulus"
OUT = ROOT / "results/etop-thermal-transport"


def load_fields():
    source = json.loads((SOURCE / "score.json").read_text())
    contract = source["contract"]
    n, batches = contract["area"], contract["batches"]
    count = contract["samples_per_shape_pair"] // batches
    fields = np.zeros((batches, 3, 2, n+1))
    for s, shape in enumerate(contract["shapes"]):
        delta = float(Fraction(shape["delta_cos4"]))
        with (SOURCE / "raw" / (shape["name"] + ".hist.csv")).open() as stream:
            for row in csv.DictReader(stream):
                b, k = int(row["batch"]), int(row["k"])
                sign = 1 if row["orientation"] == "first" else -1
                value = sign*int(row["count"])/(count*delta)
                fields[b, s, 0, k] += value
                fields[b, s, 1, k] += value*(1 if row["kind"] == "plus" else -1)
    return source, fields


def moment_kernel(n, pc, scale, degree=9):
    k = np.arange(n+1, dtype=float)
    rising = np.ones(n+1)
    beta_moments = []
    for j in range(1, degree+2):
        rising = rising*(k+j-1)/(n+j)
        beta_moments.append(rising.copy())
    return np.array([sum(comb(j, m)*(-pc)**(j-m)*(-beta_moments[m]/(m+1))
                         for m in range(j+1))*scale**j for j in range(degree+1)]).T


def jack_cov(leave_one):
    b = len(leave_one)
    return np.atleast_2d(np.cov(leave_one, rowvar=False, ddof=1))*(b-1)**2/b


def describe(value, leave_one, with_score=False):
    covariance = jack_cov(leave_one)
    se = np.sqrt(np.maximum(np.diag(covariance), 0))
    result = {"mean": np.asarray(value).tolist(), "covariance": covariance.tolist(),
              "se": se.tolist()}
    if with_score:
        statistic = float(value@np.linalg.solve(covariance, value))
        result.update({"chi_square": statistic, "df": len(value),
                       "nominal_p_value": float(chi2.sf(statistic, len(value)))})
    return result


def fit_moments(mean_moments, degree):
    d, u = mean_moments[1]-mean_moments[0], mean_moments[2]-mean_moments[0]
    r = u[:, 0]/d[:, 0]
    # Each parity has its own area-normalizing amplitude, no shared-area claim.
    residual = u-r[:, None]*d
    anchors = np.arange(1, degree+2)
    matrix = np.array([[-j*d[0, j-1+l] for l in range(degree+1)] for j in anchors])
    theta = np.linalg.solve(matrix, residual[0, anchors])
    # R_A=(1/s) d_p[v(z) D_A], v=sum theta_l z^l.
    odd_j = np.arange(degree+2, 7)
    odd_pred = np.array([-j*sum(theta[l]*d[0, j-1+l] for l in range(degree+1)) for j in odd_j])
    odd_remainder = residual[0, odd_j]-odd_pred
    # A common geometric velocity has coefficient theta/r_A. The even
    # amplitude r_E is fixed by its own area, and supplies no new warp fit.
    even_j = np.arange(1, 7)
    even_pred = np.array([-j*(r[1]/r[0])*sum(theta[l]*d[1, j-1+l] for l in range(degree+1)) for j in even_j])
    even_remainder = residual[1, even_j]-even_pred
    return {"theta": theta, "area_amplitudes": r, "odd_orders": odd_j,
            "odd_observed": residual[0, odd_j], "odd_prediction": odd_pred,
            "odd_remainder": odd_remainder, "even_orders": even_j,
            "even_observed": residual[1, even_j], "even_prediction": even_pred,
            "even_remainder": even_remainder, "anchor_condition_number": float(np.linalg.cond(matrix))}


def main():
    source, fields = load_fields()
    contract = source["contract"]
    n, b, pc = contract["area"], contract["batches"], contract["fixed_p"]
    scale = n**(3/8)
    kernel = moment_kernel(n, pc, scale)
    batch_moments = fields @ kernel
    mean_moments = batch_moments.mean(axis=0)
    loo = (b*mean_moments-batch_moments)/(b-1)
    field_mean = fields.mean(axis=0)
    d, u = field_mean[1]-field_mean[0], field_mean[2]-field_mean[0]
    grid = np.linspace(0., 1., 401)
    z = scale*(grid-pc)
    k = np.arange(n+1)
    tails = binom.sf(k[:, None]-1, n, grid[None, :])
    derivatives = n*binom.pmf(k[:, None]-1, n-1, grid[None, :])
    dcurve, ucurve, dprime = d@tails, u@tails, d@derivatives
    models = {}
    for degree, name in enumerate(("translation_tangent", "affine_velocity_tangent", "quadratic_velocity_tangent")):
        point = fit_moments(mean_moments, degree)
        replicates = [fit_moments(m, degree) for m in loo]
        r, theta = point["area_amplitudes"], point["theta"]
        velocity = sum(theta[l]*z**l for l in range(degree+1))
        dv = sum(l*theta[l]*z**(l-1) for l in range(1, degree+1))
        predicted = (dv*dcurve + (velocity/scale)*dprime) * np.array([1., r[1]/r[0]])[:, None]
        observed = ucurve-r[:, None]*dcurve
        row = {"degree": degree, "anchor_orders": list(range(1, degree+2)),
               "coefficient_meaning": "R_A=(1/s)*d_p[(sum theta_l z^l)D_A]; normalized common velocity is theta/r_A",
               "theta": describe(theta, np.array([a["theta"] for a in replicates])),
               "area_amplitudes_A_E": describe(r, np.array([a["area_amplitudes"] for a in replicates])),
               "anchor_condition_number": point["anchor_condition_number"],
               "normalized_velocity_coefficients": (theta/r[0]).tolist(),
               "p_grid": grid.tolist(), "odd_predicted_profile": predicted[0].tolist(),
               "even_predicted_profile": predicted[1].tolist(),
               "odd_observed_profile": observed[0].tolist(), "even_observed_profile": observed[1].tolist()}
        for kind in ("odd", "even"):
            row[kind+"_remaining_moments"] = {
                "orders": point[kind+"_orders"].tolist(),
                "observed": point[kind+"_observed"].tolist(),
                "prediction": point[kind+"_prediction"].tolist(),
                **describe(point[kind+"_remainder"],
                           np.array([a[kind+"_remainder"] for a in replicates]), with_score=True)}
        models[name] = row
    result = {"source_commit": "7b30648", "source_path": str((SOURCE/"score.json").relative_to(ROOT)),
              "source_sha256": hashlib.sha256((SOURCE/"score.json").read_bytes()).hexdigest(),
              "status": "post-reveal exploratory tangent-transport analysis; zero new samples",
              "dependency": source["contract"]["sampling"], "thermal_coordinate": "z=N^(3/8)(p-p_ref), a finite-N coordinate only",
              "moment_identity": "R_j=-j sum_l theta_l D_(j-1+l)",
              "models": models,
              "boundaries": ["Moment anchors and remaining moments share data; remaining means overidentifying restrictions, not independent holdouts.",
                             "Delete-one paired-batch covariance propagates clock normalization, moment anchors, source and target jointly.",
                             "These are tangent deformations. Rejecting finite-degree tangent models does not reject all finite monotone coordinate maps.",
                             "Jacobian transport preserves the signed profile integral; an ordinary observable coordinate change has no Jacobian. D is not assumed to be a physical probability density.",
                             "Even cross-channel prediction uses the same normalized velocity after its separately measured area amplitude; it is an extra mechanism hypothesis, not parity algebra."]}
    OUT.mkdir(exist_ok=True)
    (OUT/"transport.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# New N100 thermal deformation: translation, dilation, or another shape?", "",
        "The lowest odd signed moments determine a tangent-transport velocity. The remaining odd moments and the entire even profile are then genuine additional model constraints, with all same-stream source uncertainty retained. This is post-reveal mechanism analysis, not fresh evidence or a test-suite rerun.", "",
        "## Odd profile: low moments fix the transport, higher moments constrain it", "",
        "| candidate | anchored moments | theta | remaining chi-square / df | nominal p |",
        "|---|---|---|---:|---:|"]
    for name, row in models.items():
        test = row["odd_remaining_moments"]
        lines.append(f"| {name} | {row['anchor_orders']} | {row['theta']['mean']} | {test['chi_square']:.7g} / {test['df']} | {test['nominal_p_value']:.6g} |")
    lines += ["", "## Even profile: transfer the same normalized velocity, not a new fit", "",
              "Each parity has a separate area amplitude r=U_0/D_0. Odd anchors determine theta/r_A; the even moments use r_E times this same velocity.", "",
              "| candidate | even remaining chi-square / df | nominal p |", "|---|---:|---:|"]
    for name, row in models.items():
        test = row["even_remaining_moments"]
        lines.append(f"| {name} | {test['chi_square']:.7g} / {test['df']} | {test['nominal_p_value']:.6g} |")
    lines += ["", "All uncertainty is estimated by deleting each of the 200 aligned batches and refitting its area amplitudes and low-moment velocity. Gaussian-reference scores are exploratory; the candidates are not independent evidence blocks.", "", "## What the model means", "",
        "Write D=Y(4i)-Y(2i), U=Y(1/2+i)-Y(2i). After fixing the area ratio r, the candidate is U-rD=(1/s) d_p[v(z)D], with z=s(p-p_ref). Integration by parts gives R_j=-j sum_l theta_l D_(j-1+l). Constant v is a translation tangent; linear v adds dilation with the area-preserving amplitude correction; quadratic v introduces a non-affine tangent.", "",
        "This is a signed-profile transport model, not an assumption that A_top is a probability density. Ordinary observable reparameterization omits the Jacobian. The present scores do not settle finite nonlinear monotone transport, and cannot identify a continuous field.", ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
