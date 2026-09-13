#!/usr/bin/env python3
"""New fixed-N modulus transport inference from the completed P267 archive.

No raw replay or new simulation. Two independent, jointly observed A/E/C/W
vectors are profiled as noisy vectors, not as fixed regression predictors.
The new-shape predictions are conditional hypotheses, not field identifications.
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "087b07cb69f2481cbcd55fe2194150d7620835e5"
INPUT = "results/p267-etop-tau-topology-factorial/score.json"
FIELDS = ["A_top", "E_top", "C", "W"]
OUT = ROOT / "results/etop-modulus-survivors"


def e4(tau: complex) -> float:
    """Real E4 on the declared real-part 0 or 1/2 target shapes."""
    mp.mp.dps = 50
    q = mp.exp(2j * mp.pi * mp.mpc(tau))
    value = 1 + 240 * mp.fsum(n**3 * q**n / (1 - q**n) for n in range(1, 80))
    if abs(mp.im(value)) > mp.mpf("1e-40"):
        raise ValueError("This real covector is only declared on reflection-real shapes")
    return float(mp.re(value))


def shape_g(tau: complex) -> float:
    return tau.imag**2 * e4(tau) / e4(1j)


def profile_angle(theta, x, y, sx, sy):
    # Homogeneous parameterization includes lambda=+/-infinity. The covariance
    # determinant is NOT added: this profiles a fixed-covariance joint Gaussian
    # experiment over its latent source vector, not a new residual likelihood.
    c, s = np.cos(theta), np.sin(theta)
    residual = c * y - s * x
    variance = c*c * sy + s*s * sx
    return float(residual @ np.linalg.solve(variance, residual))


def fit_ray(x, y, sx, sy):
    grid = np.linspace(-np.pi / 2, np.pi / 2, 4097)
    values = np.array([profile_angle(t, x, y, sx, sy) for t in grid])
    candidates = [(float(values[0]), float(grid[0]))]
    for i in range(1, len(grid)-1):
        if values[i] <= values[i-1] and values[i] <= values[i+1]:
            r = minimize_scalar(lambda t: profile_angle(t, x, y, sx, sy),
                                bounds=(grid[i-1], grid[i+1]), method="bounded",
                                options={"xatol": 1e-14})
            candidates.append((float(r.fun), float(r.x)))
    value, theta = min(candidates)
    lam = float(np.tan(theta))
    ix, iy = np.linalg.inv(sx), np.linalg.inv(sy)
    latent = np.linalg.solve(ix + lam*lam*iy, ix @ x + lam*iy @ y)
    joint_cov = np.zeros((8, 8))
    joint_cov[:4, :4], joint_cov[4:, 4:] = sx, sy
    jac = np.zeros((8, 5))
    jac[:4, :4], jac[4:, :4], jac[4:, 4] = np.eye(4), lam*np.eye(4), latent
    parameter_cov = np.linalg.inv(jac.T @ np.linalg.solve(joint_cov, jac))
    profiles = {}
    for coverage in (.95, .99):
        cutoff = value + chi2.ppf(coverage, 1)
        crossings = []
        for i in range(len(grid)-1):
            if (values[i]-cutoff) * (values[i+1]-cutoff) < 0:
                crossings.append(brentq(
                    lambda t: profile_angle(t, x, y, sx, sy)-cutoff,
                    grid[i], grid[i+1], xtol=1e-14))
        bounds = [float(grid[0]), *crossings, float(grid[-1])]
        intervals = []
        for lo, hi in zip(bounds[:-1], bounds[1:]):
            if profile_angle((lo+hi)/2, x, y, sx, sy) <= cutoff:
                intervals.append([
                    None if lo == grid[0] else float(np.tan(lo)),
                    None if hi == grid[-1] else float(np.tan(hi)),
                ])
        profiles[str(coverage)] = intervals
    residual = y - lam*x
    v = sy + lam*lam*sx
    profiled_joint_distance = float(
        (x-latent) @ ix @ (x-latent)
        + (y-lam*latent) @ iy @ (y-lam*latent))
    return {
        "lambda": lam, "chi_square": value, "df": len(x)-1,
        "p_value": float(chi2.sf(value, len(x)-1)),
        "lambda_profile_intervals": profiles,
        "latent_tau_i": latent.tolist(),
        "latent_tau_2i": (lam*latent).tolist(),
        "parameter_order": FIELDS + ["lambda"],
        "local_parameter_covariance": parameter_cov.tolist(),
        "remaining_direction_residual": residual.tolist(),
        "remaining_direction_covariance": v.tolist(),
        "profile_vs_joint_distance_error": abs(value-profiled_joint_distance),
        "compact_boundary_chi_square": float(values[0]),
        "stationary_basin_candidates": len(candidates),
    }


def fit_parity(x, y, sx, sy, lam):
    def objective(z):
        scale = np.array([z[0], z[1], z[0], z[1]])
        residual = y-scale*x
        variance = sy+scale[:, None]*sx*scale[None, :]
        return float(residual @ np.linalg.solve(variance, residual))
    fits = [minimize(objective, start, method="Nelder-Mead",
                     options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 3000})
            for start in [(lam, lam), (3., 5.), (-3., 5.), (3., -5.)]]
    result = min(fits, key=lambda r: r.fun)
    return {"lambda_odd_A_C": float(result.x[0]),
            "lambda_even_E_W": float(result.x[1]),
            "chi_square": float(result.fun), "df": 2,
            "p_value": float(chi2.sf(result.fun, 2)),
            "optimizer_success": bool(result.success)}


def interpolation_weight(model, tau):
    """Independent offset/slope vectors: no surviving-common-ray assumption."""
    if model == "affine_log_height":
        return float(np.log(tau.imag)/np.log(2))
    if model == "affine_height_squared":
        return (tau.imag**2-1)/3
    if model == "affine_height_E4":
        return (shape_g(1j*tau.imag)-1)/(shape_g(2j)-1)
    if model == "affine_E4":
        return (shape_g(tau)-1)/(shape_g(2j)-1)
    raise ValueError(model)


def fieller_ratio(x, y, vx, vy, coverage=.95):
    """Confidence set for x/y with independent Gaussian estimates.

    None denotes an infinite endpoint, not a missing/error observation.
    Unlike delta intervals, these sets retain weak-denominator nonidentification.
    """
    cutoff = float(chi2.ppf(coverage, 1))
    a, b, c = y*y-cutoff*vy, -2*x*y, x*x-cutoff*vx
    scale = max(abs(a), abs(b), abs(c))
    if abs(a) <= 1e-13*scale:
        if abs(b) <= 1e-13*scale:
            intervals = [[None, None]] if c <= 0 else []
        else:
            root = -c/b
            intervals = [[None, root]] if b > 0 else [[root, None]]
    else:
        disc = b*b-4*a*c
        if disc < 0:
            intervals = [[None, None]] if a < 0 else []
        else:
            lo, hi = sorted(((-b-np.sqrt(disc))/(2*a),
                             (-b+np.sqrt(disc))/(2*a)))
            intervals = [[lo, hi]] if a > 0 else [[None, lo], [hi, None]]
    return {"coverage": coverage, "intervals": intervals,
            "bounded": bool(intervals) and all(v is not None for row in intervals for v in row),
            "quadratic_coefficients": [a, b, c]}


def affine_intervals(intervals, offset, slope):
    if slope == 0:
        return [[offset, offset]] if intervals else []
    answer = []
    for lo, hi in intervals:
        left = None if lo is None else float(offset+slope*lo)
        right = None if hi is None else float(offset+slope*hi)
        answer.append([left, right] if slope > 0 else [right, left])
    return sorted(answer, key=lambda row: -np.inf if row[0] is None else row[0])


def vector_predictions(x, y, sx, sy):
    ratios = x/y
    jx, jy = np.diag(1/y), np.diag(-x/y**2)
    ratio_cov = jx@sx@jx.T+jy@sy@jy.T
    fieller = {field: {str(c): fieller_ratio(x[i], y[i], sx[i, i], sy[i, i], c)
                        for c in (.95, .99)} for i, field in enumerate(FIELDS)}
    predictions = {}
    for model in ("affine_log_height", "affine_height_squared", "affine_height_E4", "affine_E4"):
        predictions[model] = {}
        for name, tau in {"tau_4i": 4j, "tau_half_plus_i": .5+1j}.items():
            t = interpolation_weight(model, tau)
            point = t+(1-t)*ratios
            predictions[model][name] = {
                "t": t,
                "old_N50_prototype_not_absolute_N100_forecast": ((1-t)*x+t*y).tolist(),
                "prototype_covariance": ((1-t)**2*sx+t*t*sy).tolist(),
                "same_N100_ratio_to_2i": dict(zip(FIELDS, point.tolist())),
                "ratio_delta_covariance_diagnostic_only": ((1-t)**2*ratio_cov).tolist(),
                "source_fieller_sets": {field: {c: affine_intervals(
                    fieller[field][c]["intervals"], t, 1-t) for c in ("0.95", "0.99")}
                    for field in FIELDS},
            }
    return {"source_ratio_i_over_2i": dict(zip(FIELDS, ratios.tolist())),
            "source_ratio_fieller": fieller,
            "source_ratio_delta_covariance_diagnostic_only": ratio_cov.tolist(),
            "models": predictions}


def calculate():
    raw = subprocess.check_output(["git", "show", f"{SOURCE}:{INPUT}"], cwd=ROOT)
    data = json.loads(raw)
    first, second = [data["projected_rows"][k] for k in ("tau_i", "tau_2i")]
    x, y = [np.array(row["estimate"], dtype=float) for row in (first, second)]
    sx, sy = [np.array(row["covariance"], dtype=float) for row in (first, second)]
    ray = fit_ray(x, y, sx, sy)
    lam = ray["lambda"]
    fixed = {}
    for name, value in [("no_modulus_response", 1.), ("height", 2.),
                         ("height_squared", 4.), ("pure_area_E4", shape_g(2j))]:
        q = profile_angle(np.arctan(value), x, y, sx, sy)
        fixed[name] = {"lambda": value, "chi_square": q, "df": 4,
                       "p_value": float(chi2.sf(q, 4)),
                       "penalty_vs_free_ray": q-ray["chi_square"]}
    parity = fit_parity(x, y, sx, sy, lam)
    parity["improvement_over_common_ray"] = ray["chi_square"]-parity["chi_square"]
    parity["improvement_p_value_1df"] = float(chi2.sf(
        max(0., parity["improvement_over_common_ray"]), 1))
    predictions = vector_predictions(x, y, sx, sy)
    return {
        "schema": "matching-one/etop-modulus-survivors/v2",
        "status": "retrospective_new_transport_analysis_no_new_samples",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input": {"commit": SOURCE, "path": INPUT,
                  "sha256": hashlib.sha256(raw).hexdigest(),
                  "field_order": FIELDS, "fixed_p": .59274605079,
                  "area": 50, "rows": {"tau_i": first, "tau_2i": second},
                  "independence": "old P205 N50 and new factorial block have disjoint streams; preserve full within-row covariance",
                  "dependency_warning": "P205 N50 is reused by other Etop summaries; their p-values are not pooled here"},
        "free_common_amplitude": ray, "fixed_amplitudes": fixed,
        "two_parity_amplitudes": parity,
        "shape_values": {str(t): shape_g(t) for t in (1j, 2j, 4j, .5+1j)},
        "conditional_N100_predictions": predictions,
        "prediction_contract": {
            "same_area": 100, "reference_modulus": "2i",
            "targets": ["4i", "1/2+i"],
            "source_shape_assumption": "Each coordinate has independent offset/slope in its shape coordinate. Both coefficients receive the same unknown coordinate-specific gain c_j(N) between N50 and N100. This shape/area separability is new and testable, not measured scaling.",
            "amplitude": "Arbitrary separate gain c_j(N100) for each A/E/C/W cancels in the SAME-N100 target/reference ratio; no N^-13/8 or other area law is imposed",
            "models": {
                "affine_log_height": "v(tau)=(1-t)x+t*y, t=log(Im tau)/log 2",
                "affine_height_squared": "v(tau)=(1-t)x+t*y, t=((Im tau)^2-1)/3",
                "affine_height_E4": "v(tau)=(1-t)x+t*y, t=(g(i Im tau)-1)/(g(2i)-1); height-only E4 adversary",
                "affine_E4": "v(tau)=(1-t)x+t*y, t=(g(tau)-1)/(g(2i)-1), g=Im(tau)^2 Re E4(tau)/E4(i), on reflection-real shapes",
            },
            "uncertainty": "Coordinatewise 95%/99% Fieller sets propagated from x_j/y_j; these are not simultaneous sets. Full delta covariance is diagnostic only for unresolved denominators. Future target noise and shape/area discrepancy are additional.",
            "chronology": "Models formulated after source reveal; do not label this a prospective source fit. No target data acquired.",
        },
        "geometry_source": "b9e4ea19bc585cbed18ec6ba1d13e85f2b5accc7:results/p267-third-geometry-feasibility/certificate.json",
        "boundary": "The common four-vector amplitude is rejected. Its lambda interval describes that rejected family and is NOT used for forecasts. Two source shapes saturate each vector-affine family, so none is source-selected. Pure E4 common-ray rejection does not exclude all thermal Q4 couplings. Common/parity scores use approximate Gaussian covariance inference.",
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "mpmath": mp.__version__},
    }


def report(data):
    ray, parity = data["free_common_amplitude"], data["two_parity_amplitudes"]
    def intervals_text(intervals):
        return " U ".join("["+("-inf" if lo is None else f"{lo:.6g}")+", "+
                          ("+inf" if hi is None else f"{hi:.6g}")+"]"
                          for lo, hi in intervals) or "empty"
    lines = ["# E_top: modulus interaction versus direction transport", "",
             "## Geometry changes the response direction, not only its amplitude", "",
             "The completed N50 archive rejects one common amplitude for the four coordinates A_top,E_top,C,W. This is a new retrospective model comparison using zero new samples; it is not an E_top-only detection or field identification.", "",
             f"Free common amplitude: lambda={ray['lambda']:.9g}; chi2={ray['chi_square']:.6g}/3, p={ray['p_value']:.6g}.",
             "The lambda profile interval belongs to the rejected model and is not used in the predictions below.", "",
             "A and E are fixed-p rank probabilities; C and W are integrated clock/lifetime coordinates. Every input is a first-minus-second direction contrast divided by its exact cos(4 theta) contrast, at N50 and p_ref=.59274605079. Two random streams are independent, while each four-vector retains its full covariance.", "",
             "| Declared common-vector shape | lambda | chi2/4 | p |", "|---|---:|---:|---:|"]
    for name, row in data["fixed_amplitudes"].items():
        lines.append(f"| {name} | {row['lambda']:.7g} | {row['chi_square']:.6g} | {row['p_value']:.6g} |")
    lines += ["", "## Is a second transport amplitude needed?", "",
              f"Odd A/C lambda={parity['lambda_odd_A_C']:.7g}; even E/W lambda={parity['lambda_even_E_W']:.7g}.",
              f"Two-amplitude chi2={parity['chi_square']:.6g}/2, p={parity['p_value']:.6g}; improvement={parity['improvement_over_common_ray']:.6g}/1, p={parity['improvement_p_value_1df']:.6g}.", "",
              "The parity split is a narrow survivor at alpha=.01 and fails at .05. It is a useful direction, not a selected mechanism. Earlier A+C-plane compatibility is not contradicted: a four-vector can change direction while still satisfying E=beta*A+gamma*C. That earlier test used different geometries/centers and is not pooled here.", "",
              "## A sheared geometry separates the conditional shape hypotheses", "",
              "A vector-affine family uses v(tau)=(1-t)x_i+t*x_2i with independent offset and slope vectors, so the source vectors need not share a ray. The four choices of t are log(height), height squared, height-only E4, and the actual complex-modulus E4 covector. Each exactly interpolates the two source means; none is selected by those two points.", "",
              "For N100, assume both affine coefficients for coordinate j receive the same unknown gain c_j(N100). This is an explicit new shape/area-separability hypothesis. In the SAME-N100 ratio to a new 2i bridge, that gain cancels:", "",
              "`R_j(tau)/R_j(2i) = t+(1-t)*x_j(50,i)/x_j(50,2i)`.", "",
              "No N^-13/8 or other area exponent is imposed. The old-N50 affine vectors in JSON are prototypes, not absolute N100 amplitude predictions. Source-only 95% Fieller sets below retain weak-denominator nonidentification; they are coordinatewise, not simultaneous.", "",
              "| Shape family | target | coordinate | same-N100 ratio | source 95% Fieller set |", "|---|---|---|---:|---|"]
    for model, p in data["conditional_N100_predictions"]["models"].items():
        for target, row in p.items():
            for field in FIELDS:
                lines.append(f"| {model} | {target} | {field} | {row['same_N100_ratio_to_2i'][field]:.7g} | {intervals_text(row['source_fieller_sets'][field]['0.95'])} |")
    e4_shear = data["conditional_N100_predictions"]["models"]["affine_E4"]["tau_half_plus_i"]
    height_shear = data["conditional_N100_predictions"]["models"]["affine_height_squared"]["tau_half_plus_i"]
    lines += ["", f"The cleanest source-resolved contrast is C: actual-modulus affine-E4 predicts shear/2i={e4_shear['same_N100_ratio_to_2i']['C']:.6g}, whereas every listed height-only family predicts {height_shear['same_N100_ratio_to_2i']['C']:.6g}. Their source uncertainty sets are distinct; future target noise and model discrepancy remain to be measured. E_top's denominator is weak, so its unbounded Fieller set cannot be presented as a precise forecast. A possible W sign reversal is a hypothesis, not a resolved sign prediction.", "",
              "N50 with Smith(1,50)/(5,10) admits no third modulus. N100 with Smith(1,100)/(5,20) admits 2i,4i,1/2+i. The exact proof/matrices are at commit b9e4ea1. The bridge-and-shear pair is a smaller source-informed test; all three N100 shapes permit a new-area offset/amplitude-free test that drops the cross-area separability assumption.", "",
              "## Dependence and scope", "",
              data["boundary"], "",
              "The likelihood profiles the latent source vector with both 4x4 covariance blocks. It does not treat the source point as exact, add a residual log-determinant, combine prior p-values, or rerun the factorial score. The newly selected model comparisons and shape forecasts are retrospective, exploratory outputs, not preregistered source evidence.", "",
              "Source: `"+SOURCE+":"+INPUT+"`, SHA256 `"+data['input']['sha256']+"`.", "",
              "Reproduce with `/Users/lc/python-envs/research-py311/bin/python scripts/etop_modulus_survivors.py`. The pinned source Git object must be present (fetch the PR267 source branch if needed).", ""]
    return "\n".join(lines)


if __name__ == "__main__":
    result = calculate()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "latest.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    (OUT / "REPORT.md").write_text(report(result))
    print(report(result))
