#!/usr/bin/env python3
"""Compare the independent N100/N400 shape blocks without imposing an area law."""
import json

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2

from etop_thermal_transport import ROOT, load_fields, moment_kernel, describe

OUT = ROOT/"results/etop-n100-n400-scale-transport"


def profile_ray(x, sx, y, sy):
    def objective(theta):
        c, s = np.cos(theta), np.sin(theta)
        r = c*y-s*x
        cov = c*c*sy+s*s*sx
        return float(r@np.linalg.solve(cov, r))
    grid = np.linspace(-np.pi/2, np.pi/2, 2049)
    values = np.array([objective(t) for t in grid])
    fits = [(float(values[0]), float(grid[0]))]
    for i in range(1, len(grid)-1):
        if values[i] <= values[i-1] and values[i] <= values[i+1]:
            fit = minimize_scalar(objective, bounds=(grid[i-1], grid[i+1]), method="bounded")
            fits.append((float(fit.fun), float(fit.x)))
    q, theta = min(fits)
    return {"amplitude_N400_over_N100": float(np.tan(theta)), "chi_square": q,
            "df": len(x)-1, "nominal_p_value": float(chi2.sf(q, len(x)-1)),
            "scope": "One freely profiled signed amplitude; source and target covariances both retained, with independent scale blocks."}


def clock_moments(directory):
    source, fields = load_fields(directory)
    n, pc = source["contract"]["area"], source["contract"]["fixed_p"]
    b = len(fields)
    moments = fields @ moment_kernel(n, pc, n**(3/8), 4)
    mean = moments.mean(axis=0)
    loo = (b*mean-moments)/(b-1)
    def value(q):
        d, u = q[1, 0]-q[0, 0], q[2, 0]-q[0, 0]
        md, mu = d[1]/d[0], u[1]/u[0]
        vd, vu = d[2]/d[0]-md**2, u[2]/u[0]-mu**2
        return np.array([d[0], u[0], md, mu, mu-md, vd, vu])
    labels = ["D_A_area", "U_A_area", "D_mean_z", "U_mean_z", "U_minus_D_mean_z",
              "D_centered_second_z", "U_centered_second_z"]
    return {"area": n, "labels": labels, **describe(value(mean), np.array([value(q) for q in loo])),
            "scope": "Area-normalized signed-profile moments, with z=N^(3/8)(p-p_ref); no exponent was fitted or proved by choosing this coordinate."}


def main():
    sources = {}
    for n, directory in ((100, "etop-finite-transport-invariants"),
                          (400, "etop-n400-finite-transport-invariants")):
        sources[n] = json.loads((ROOT/"results"/directory/"invariants.json").read_text())
    old, new = (sources[n]["finite_transport_remainders"] for n in (100, 400))
    x, y = np.array(old["mean"]), np.array(new["mean"])
    sx, sy = np.array(old["covariance"]), np.array(new["covariance"])
    ray = profile_ray(x, sx, y, sy)
    clocks = {str(n): clock_moments(ROOT/"results"/f"etop-n{n}-three-modulus") for n in (100, 400)}
    delta = np.array(clocks["400"]["mean"])[2:] - np.array(clocks["100"]["mean"])[2:]
    cov = (np.array(clocks["400"]["covariance"])+np.array(clocks["100"]["covariance"]))[2:, 2:]
    result = {"status": "Independent new N400 scale compared with the exploratory N100 source",
        "source_freezes": {str(n): sources[n].get("input_prediction_freeze_commit", "4c1ec50") for n in (100, 400)},
        "source_score_sha256": {str(n): sources[n]["source_sha256"] for n in (100, 400)},
        "profile_remainder_direction_transfer": ray, "clock_moments": clocks,
        "clock_scale_difference": {"labels": clocks["100"]["labels"][2:], "mean": delta.tolist(),
                                   "covariance": cov.tolist(), "se": np.sqrt(np.diag(cov)).tolist()},
        "boundary": "Two areas do not establish an asymptotic exponent. A surviving direction-transfer comparison is not proof of one field; a weak target is not recovery of the common-coordinate model. No p-values from source and target are added."}
    OUT.mkdir(exist_ok=True)
    (OUT/"scale_transport.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# N100 to N400: an actual new scale, not a rescore", "",
        "All periods are doubled. The three moduli and rational rotation are fixed, while the Smith pair changes uniformly from(1,100)/(5,20) to(2,200)/(10,40). The N400 seed/block is independent of N100.", "",
        "## Does the coordinate-free E deformation keep its direction?", "",
        f"A freely signed amplitude gives N400/N100={ray['amplitude_N400_over_N100']:.9g}, chi-square={ray['chi_square']:.8g}/5, nominal p={ray['nominal_p_value']:.6g}. Both source and target covariance matrices are used. No area exponent is imposed.", "",
        "## Full-p odd clock profiles in the declared thermal coordinate", "",
        "| readout | N100 | N400 |", "|---|---:|---:|"]
    for j, label in enumerate(clocks["100"]["labels"]):
        row = [f"{clocks[str(n)]['mean'][j]:.8g} +/- {clocks[str(n)]['se'][j]:.3g}" for n in (100, 400)]
        lines.append("| "+label+" | "+" | ".join(row)+" |")
    lines += ["", "These are signed-profile moment ratios. Their finite-scale change is measured directly, without claiming a continuum collapse from two sizes.", "", result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
