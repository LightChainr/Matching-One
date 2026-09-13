#!/usr/bin/env python3
"""Compare the independent N100/N400 shape blocks without imposing an area law."""
import json

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import roots_legendre
from scipy.stats import binom
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
    cutoff = q + chi2.ppf(.95, 1)
    edges = [-np.pi/2]
    for a, z, fa, fz in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if (fa-cutoff)*(fz-cutoff) < 0:
            edges.append(float(brentq(lambda t: objective(t)-cutoff, a, z)))
    edges.append(np.pi/2)
    intervals = [[None if a == -np.pi/2 else float(np.tan(a)),
                  None if z == np.pi/2 else float(np.tan(z))]
                 for a, z in zip(edges[:-1], edges[1:]) if objective((a+z)/2) <= cutoff]
    return {"amplitude_N400_over_N100": float(np.tan(theta)), "chi_square": q,
            "df": len(x)-1, "nominal_p_value": float(chi2.sf(q, len(x)-1)),
            "approximate_95_profile_intervals": intervals,
            "zero_amplitude_delta_chi_square": float(objective(0)-q),
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


def window_moments(directory):
    """Full-p versus shrinking central windows, on the same exact polynomials."""
    source, fields = load_fields(directory)
    n, pc = source["contract"]["area"], source["contract"]["fixed_p"]
    scale, b = n**(3/8), len(fields)
    windows = [("full_p", 0., 1.)] + [
        (f"abs_z_le_{h:g}", max(0., pc-h/scale), min(1., pc+h/scale))
        for h in (.5, 1., 2.)]
    x, w = roots_legendre((n+4)//2)
    kernels = []
    for _, low, high in windows:
        p = low+(x+1)*(high-low)/2
        tails = binom.sf(np.arange(n+1)[:, None]-1, n, p[None, :])
        for j in range(3):
            kernels.append(tails @ (w*(high-low)/2*(scale*(p-pc))**j))
    moments = fields @ np.array(kernels).T
    mean = moments.mean(axis=0)
    loo = (b*mean-moments)/(b-1)
    def value(q):
        d, u = (q[1]-q[0]).reshape(2, len(windows), 3), (q[2]-q[0]).reshape(2, len(windows), 3)
        r = u[0, 0, 0]/d[0, 0, 0]
        residual = u-r*d
        return np.concatenate([d[0, :, 0]/d[0, 0, 0], u[0, :, 0]/u[0, 0, 0],
                               residual[0].reshape(-1), residual[1].reshape(-1)])
    labels = (["D_A_fraction:"+v[0] for v in windows] + ["U_A_fraction:"+v[0] for v in windows]
              + [f"R_{s}:{v[0]}:z{j}" for s in ("A", "E") for v in windows for j in range(3)])
    return {"area": n, "windows": [{"name": v[0], "p_low": v[1], "p_high": v[2]} for v in windows],
            "labels": labels, **describe(value(mean), np.array([value(q) for q in loo])),
            "scope": "Signed area fractions and clock-calibrated residual moments; window widths 0.5,1,2 are post-reveal descriptive choices. All windows share their source block. Full odd residual area is zero by construction."}


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
    windows = {str(n): window_moments(ROOT/"results"/f"etop-n{n}-three-modulus") for n in (100, 400)}
    delta = np.array(clocks["400"]["mean"])[2:] - np.array(clocks["100"]["mean"])[2:]
    cov = (np.array(clocks["400"]["covariance"])+np.array(clocks["100"]["covariance"]))[2:, 2:]
    result = {"status": "Independent new N400 scale compared with the exploratory N100 source",
        "source_freezes": {str(n): sources[n].get("input_prediction_freeze_commit", "4c1ec50") for n in (100, 400)},
        "source_score_sha256": {str(n): sources[n]["source_sha256"] for n in (100, 400)},
        "profile_remainder_direction_transfer": ray, "clock_moments": clocks, "window_moments": windows,
        "clock_scale_difference": {"labels": clocks["100"]["labels"][2:], "mean": delta.tolist(),
                                   "covariance": cov.tolist(), "se": np.sqrt(np.diag(cov)).tolist()},
        "boundary": "Two areas do not establish an asymptotic exponent. A surviving direction-transfer comparison is not proof of one field; a weak target is not recovery of the common-coordinate model. No p-values from source and target are added."}
    OUT.mkdir(exist_ok=True)
    (OUT/"scale_transport.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# N100 to N400: an actual new scale, not a rescore", "",
        "All periods are doubled. The three moduli and rational rotation are fixed, while the Smith pair changes uniformly from(1,100)/(5,20) to(2,200)/(10,40). The N400 seed/block is independent of N100.", "",
        "## Does the coordinate-free E deformation keep its direction?", "",
        f"A freely signed amplitude gives N400/N100={ray['amplitude_N400_over_N100']:.9g}, chi-square={ray['chi_square']:.8g}/5, nominal p={ray['nominal_p_value']:.6g}. Both source and target covariance matrices are used. No area exponent is imposed.", "",
        f"Approximate 95% profile intervals for this amplitude: {ray['approximate_95_profile_intervals']}. Zero amplitude costs only delta-chi-square={ray['zero_amplitude_delta_chi_square']:.6g}; the best-fit sign must not be described as a measured reversal.", "",
        "## Full-p odd clock profiles in the declared thermal coordinate", "",
        "| readout | N100 | N400 |", "|---|---:|---:|"]
    for j, label in enumerate(clocks["100"]["labels"]):
        row = [f"{clocks[str(n)]['mean'][j]:.8g} +/- {clocks[str(n)]['se'][j]:.3g}" for n in (100, 400)]
        lines.append("| "+label+" | "+" | ".join(row)+" |")
    lines += ["", "These are signed-profile moment ratios. Their finite-scale change is measured directly, without claiming a continuum collapse from two sizes.", "",
              "## Does the odd signed area concentrate inside a critical-width window?", "",
              "The declared coordinate is z=N^(3/8)(p-p_ref). Window half-widths below are exploratory, not extra independent observations.", "",
              "| profile and window | N100 signed area fraction | N400 signed area fraction |", "|---|---:|---:|"]
    for j, label in enumerate(windows["100"]["labels"][:8]):
        row = [f"{windows[str(n)]['mean'][j]:.8g} +/- {windows[str(n)]['se'][j]:.3g}" for n in (100, 400)]
        lines.append("| "+label+" | "+" | ".join(row)+" |")
    lines += ["", "The JSON also retains every full/window A/E residual moment through z^2 and their joint covariance. These are signed profiles, so an area fraction is not automatically a probability.", "", result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
