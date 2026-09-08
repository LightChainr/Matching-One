#!/usr/bin/env python3
"""Post-reveal two equal-width Gaussian lobes: fit moments 3/4, predict 5/6."""
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.optimize import brentq
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/etop-two-lobe-moment-closure"
SOURCE_REV = "6d8a3ed9"
SOURCE_PATH = "results/p267-standardized-rank-shape/score.json"
LABELS = ["right_lobe_weight", "between_lobes_variance_fraction",
          "within_lobe_sigma", "left_lobe_center", "right_lobe_center",
          "predicted_mu5", "predicted_mu6", "mu5_residual", "mu6_residual"]


def closure(moments):
    m3, m4, m5, m6 = moments[:4]
    c4 = m4 - 3
    # r is the fraction of unit variance in the two-point lobe centers.
    # c3=t*r^(3/2), c4=(t^2-2)*r^2, hence the unique positive root below.
    lower = np.sqrt(max(0., -c4/2))
    equation = lambda r: 2*r**3+c4*r-m3*m3
    if lower >= 1 or equation(1) < 0:
        raise ValueError("Moment pair does not admit positive equal-width Gaussian lobes")
    r = brentq(equation, max(lower, 1e-14), 1, xtol=1e-14)
    t = m3/r**1.5
    w = (1-t/np.sqrt(t*t+4))/2
    c5 = t*(t*t-8)*r**2.5
    c6 = (t**4-22*t*t+16)*r**3
    predicted5 = c5 + 10*m3
    predicted6 = c6 + 15*c4 + 10*m3*m3 + 15
    centers = [-np.sqrt(r*w/(1-w)), np.sqrt(r*(1-w)/w)]
    return np.array([w, r, np.sqrt(1-r), *centers,
                     predicted5, predicted6, m5-predicted5, m6-predicted6])


def covariance(loo):
    centered = loo-loo.mean(axis=0)
    return (len(loo)-1)/len(loo) * centered.T@centered


SYMMETRIC_LABELS = ["right_lobe_weight", "between_lobes_variance_fraction",
                    "kernel_variance", "kernel_fourth_moment",
                    "kernel_sixth_moment", "kernel_moment_Hankel_determinant"]


def symmetric_kernel(moments):
    """Any common symmetric positive kernel, not necessarily Gaussian."""
    m3, m4, m5, m6 = moments[:4]
    c4, c5 = m4-3, m5-10*m3
    c6 = m6-15*m4-10*m3*m3+30
    if abs(m3) < 1e-14:
        raise ValueError("This odd-moment identification needs nonzero skewness")
    # c5/c3 = t^2*r-8*r = c3^2/r^2-8*r; monotone in r.
    equation = lambda r: 8*r+c5/m3-m3*m3/r**2
    if equation(1) < 0:
        raise ValueError("Odd moments require between-center variance above total variance")
    r = brentq(equation, 1e-10, 1, xtol=1e-14)
    t, v = m3/r**1.5, 1-r
    w = (1-t/np.sqrt(t*t+4))/2
    kernel_c4 = c4-(t*t-2)*r*r
    kernel_c6 = c6-(t**4-22*t*t+16)*r**3
    kernel_m4 = kernel_c4+3*v*v
    kernel_m6 = kernel_c6+15*kernel_c4*v+15*v**3
    return np.array([w, r, v, kernel_m4, kernel_m6, v*kernel_m6-kernel_m4**2])


def score(delta, cov):
    q = float(delta@np.linalg.solve(cov, delta))
    return {"chi2": q, "df": len(delta), "nominal_p": float(chi2.sf(q, len(delta)))}


def main():
    revision = subprocess.check_output(["git", "rev-parse", SOURCE_REV], cwd=ROOT, text=True).strip()
    blob = subprocess.check_output(["git", "show", revision+":"+SOURCE_PATH], cwd=ROOT)
    source = json.loads(blob)
    result = {"source_commit": revision, "source_path": SOURCE_PATH,
              "source_sha256": hashlib.sha256(blob).hexdigest(), "labels": LABELS,
              "interpretation": "Post-reveal deterministic mechanism model chosen after observing two peaks and standardized moments; no independent new evidence, no Gaussian assumption imposed on signed data. Moments 3/4 identify two lobe parameters and moments 5/6 are the overidentifying readout. Source fitting repeated inside every aligned delete-one replicate.",
              "sources": {}, "common_symmetric_kernel": {"labels": SYMMETRIC_LABELS, "sources": {},
                  "interpretation": "Structural extension after the Gaussian closure result: moments 3/5 identify the two-center component for any common symmetric lobe kernel. Its reconstructed sixth even moment and 2x2 moment Hankel determinant must be nonnegative. No extra fitted kernel family."}}
    for n, row in source["sources"].items():
        point = closure(row["estimate"])
        loo = np.array([closure(v) for v in row["loo_vectors"]])
        cov = covariance(loo)
        result["sources"][n] = {"estimate": point.tolist(), "se": np.sqrt(cov.diagonal()).tolist(),
                                 "covariance": cov.tolist(), "loo_vectors": loo.tolist(),
                                 "residual_score": score(point[-2:], cov[-2:, -2:])}
        sym_point = symmetric_kernel(row["estimate"])
        sym_loo = np.array([symmetric_kernel(v) for v in row["loo_vectors"]])
        sym_cov = covariance(sym_loo)
        result["common_symmetric_kernel"]["sources"][n] = {
            "estimate": sym_point.tolist(), "se": np.sqrt(sym_cov.diagonal()).tolist(),
            "covariance": sym_cov.tolist(), "loo_vectors": sym_loo.tolist()}
    a, b = (result["sources"][n] for n in ("100", "400"))
    diff = np.array(b["estimate"])-np.array(a["estimate"])
    cov = np.array(a["covariance"])+np.array(b["covariance"])
    result["N400_minus_N100"] = {"estimate": diff.tolist(), "se": np.sqrt(cov.diagonal()).tolist(),
                                  "covariance": cov.tolist(),
                                  "residual_change_score": score(diff[-2:], cov[-2:, -2:])}
    OUT.mkdir(exist_ok=True)
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Two-lobe Gaussian moment closure", "", result["interpretation"], "",
             "| Quantity | N100 | N400 | change | change SE |", "|---|---:|---:|---:|---:|"]
    for k, label in enumerate(LABELS):
        lines.append(f"| {label} | {a['estimate'][k]:.9g} | {b['estimate'][k]:.9g} | {diff[k]:.8g} | {np.sqrt(cov[k,k]):.6g} |")
    for n, row in result["sources"].items():
        lines += ["", f"N{n} residual moments 5/6: {row['residual_score']}"]
    lines += ["", f"Cross-scale residual change: {result['N400_minus_N100']['residual_change_score']}", ""]
    lines += ["## Any common symmetric lobe kernel", "",
              result["common_symmetric_kernel"]["interpretation"], "",
              "| Derived necessary quantity | N100 | SE | N400 | SE |", "|---|---:|---:|---:|---:|"]
    sa, sb = (result["common_symmetric_kernel"]["sources"][n] for n in ("100", "400"))
    for k, label in enumerate(SYMMETRIC_LABELS):
        lines.append(f"| {label} | {sa['estimate'][k]:.9g} | {sa['se'][k]:.6g} | {sb['estimate'][k]:.9g} | {sb['se'][k]:.6g} |")
    lines.append("")
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
