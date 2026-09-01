#!/usr/bin/env python3
"""Compare amplitude-only and direction-changing P267 responses on saved data."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import time

import numpy as np
import scipy
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def profile(x, y, cs, cr, points, confidence):
    """Profile latent true x from the joint likelihood; no log(det S) term."""
    def objective(phi):
        c, s = math.cos(phi), math.sin(phi)
        residual = c * y - s * x
        variance = c * c * cr + s * s * cs
        return float(residual @ np.linalg.solve(variance, residual))

    grid = np.linspace(-math.pi / 2, math.pi / 2, points)
    values = np.array([objective(float(phi)) for phi in grid])
    candidates = [(float(grid[0]), float(values[0])),
                  (float(grid[-1]), float(values[-1]))]
    for i in range(1, points - 1):
        if values[i] <= values[i - 1] and values[i] <= values[i + 1]:
            fit = minimize_scalar(objective, method="bounded",
                                  bounds=(grid[i - 1], grid[i + 1]),
                                  options={"xatol": 1e-14})
            candidates.append((float(fit.x), float(fit.fun)))
    phi, distance = min(candidates, key=lambda pair: pair[1])
    threshold = distance + float(chi2.ppf(confidence, 1))
    crossings = []
    for i in range(points - 1):
        if (values[i] - threshold) * (values[i + 1] - threshold) < 0:
            crossings.append(float(brentq(lambda z: objective(z) - threshold,
                                          grid[i], grid[i + 1], xtol=1e-14)))
    edges = [-math.pi / 2] + crossings + [math.pi / 2]
    intervals = []
    for left, right in zip(edges[:-1], edges[1:]):
        if objective((left + right) / 2) <= threshold:
            intervals.append({
                "phi": [left, right],
                "lambda": [None if left == -math.pi / 2 else math.tan(left),
                           None if right == math.pi / 2 else math.tan(right)],
                "null_endpoint_means": "negative_or_positive_infinity",
            })
    finite = abs(math.cos(phi)) > 1e-10
    result = {
        "lambda": math.tan(phi) if finite else None,
        "phi": phi, "chi_square": distance, "nominal_df": len(x) - 1,
        "nominal_p": float(chi2.sf(distance, len(x) - 1)),
        "lambda_one_chi_square": objective(math.pi / 4),
        "projective_infinity_chi_square": objective(math.pi / 2),
        "profile_confidence": confidence, "profile_intervals": intervals,
        "profile_interval_boundary": "conditional_parameter_range_inside_a_model_that_may_be_rejected",
        "refined_candidates_phi_distance": candidates,
        "search": "dense_compact_projective_grid_then_local_refinement_not_an_exact_global_certificate",
    }
    if finite:
        lam = result["lambda"]
        residual = y - lam * x
        variance = cr + lam * lam * cs
        latent = x + lam * cs @ np.linalg.solve(variance, residual)
        result.update({"residual": residual.tolist(), "residual_covariance_at_fit": variance.tolist(),
                       "fitted_latent_x": latent.tolist(), "fitted_y": (lam * latent).tolist()})
    return result


def schur_allocation(residual, variance, first):
    """Exact quadratic partition at the fitted lambda; order-dependent diagnostic."""
    second = [i for i in range(len(residual)) if i not in first]
    va = variance[np.ix_(first, first)]
    cross = variance[np.ix_(second, first)]
    ra, rb = residual[first], residual[second]
    conditional_r = rb - cross @ np.linalg.solve(va, ra)
    conditional_v = variance[np.ix_(second, second)] - cross @ np.linalg.solve(va, cross.T)
    return {"first_indices": first, "conditional_indices": second,
            "marginal_distance": float(ra @ np.linalg.solve(va, ra)),
            "conditional_distance": float(conditional_r @ np.linalg.solve(conditional_v, conditional_r)),
            "conditional_residual": conditional_r.tolist(),
            "conditional_covariance": conditional_v.tolist()}


def future_geometry():
    """Concrete integer period matrices, not a production request or scale extrapolation."""
    matrices = [
        (1, [[13, -9], [9, 13]], [[5, -15], [15, 5]]),
        (2, [[11, 4], [-2, 22]], [[10, -10], [5, 20]]),
        (5, [[7, -5], [1, 35]], [[5, -25], [5, 25]]),
    ]
    def det(m):
        return m[0][0] * m[1][1] - m[0][1] * m[1][0]
    def cos4(m):
        a, b = m[0][0], m[1][0]
        return Fraction(a**4 - 6*a*a*b*b + b**4, (a*a+b*b)**2)
    rows = []
    for rho, left, right in matrices:
        rotated = [[Fraction(4*left[0][j]-3*left[1][j], 5) for j in range(2)],
                   [Fraction(3*left[0][j]+4*left[1][j], 5) for j in range(2)]]
        assert rotated == right
        smith = []
        for m in [left, right]:
            d1 = math.gcd(*(v for row in m for v in row))
            assert det(m) == 250
            assert sum(m[i][0]*m[i][1] for i in range(2)) == 0
            assert sum(m[i][1]**2 for i in range(2)) == rho*rho*sum(m[i][0]**2 for i in range(2))
            smith.append([d1, det(m)//d1])
        rows.append({"rho": rho, "cyclic": left, "noncyclic": right, "smith": smith,
                     "delta_cos4_cyclic_minus_noncyclic": str(cos4(left)-cos4(right))})
    return {"N": 250, "common_map": "(1/5)[[4,-3],[3,4]]", "rows": rows,
            "status": "exact_feasible_future_design_not_acquired_or_powered",
            "N50_obstruction": "Smith(5,10) implies L=5B with detB=2; the three index-two sublattices have only square and rectangle shapes up to square-lattice symmetry",
            "conditional_shape_candidates": {
                "log_even": {"definition": "f(rho)=(log(rho)/log(2))^2", "f5": (math.log(5)/math.log(2))**2},
                "geometric_even": {"definition": "f(rho)=2*(rho+1/rho-2)", "f5": 6.4},
                "prediction": "u5=(1-f5)*u1+f5*u2, using anchors AT N250, not current N50 vectors",
                "training_covariance_if_independent": "(1-f5)^2*C1+f5^2*C2; shared streams require cross terms",
                "boundary": "two explicitly post-reveal phenomenological hypotheses, not identified laws, a model ranking or a sampling authorization"}}


def analyze(manifest):
    start = time.perf_counter()
    source = manifest["source"]
    raw = subprocess.check_output(["git", "show", source["commit"]+":"+source["path"]], cwd=ROOT)
    if sha(raw) != source["sha256"]:
        raise ValueError("Pinned input hash differs")
    parent = json.loads(raw)
    fields = parent["field_order"]
    if fields != manifest["field_order"]:
        raise ValueError("Coordinate order differs")
    a, b = parent["projected_rows"]["tau_i"], parent["projected_rows"]["tau_2i"]
    x, y = np.asarray(a["estimate"]), np.asarray(b["estimate"])
    cs, cr = np.asarray(a["covariance"]), np.asarray(b["covariance"])
    for covariance in (cs, cr):
        np.linalg.cholesky(covariance)
        if not np.allclose(covariance, covariance.T, atol=1e-20, rtol=1e-12):
            raise ValueError("Covariance is not symmetric")
    kwargs = {"points": manifest["projective_grid_points"], "confidence": manifest["profile_confidence"]}
    full = profile(x, y, cs, cr, **kwargs)
    diagnostics = []
    for names in manifest["coordinate_diagnostics"]:
        ix = [fields.index(name) for name in names]
        diagnostics.append({"fields": names, **profile(x[ix], y[ix], cs[np.ix_(ix,ix)], cr[np.ix_(ix,ix)], **kwargs)})
    residual, variance = np.asarray(full["residual"]), np.asarray(full["residual_covariance_at_fit"])
    partitions = [schur_allocation(residual, variance, [0,1]), schur_allocation(residual, variance, [2,3])]
    checks = {"fixed_lambda_parent_absolute_difference": abs(full["lambda_one_chi_square"]-parent["primary_character_normalized_interaction"]["chi_square"]),
              "largest_Schur_sum_difference": max(abs(p["marginal_distance"]+p["conditional_distance"]-full["chi_square"]) for p in partitions)}
    if max(checks.values()) > 1e-7:
        raise ArithmeticError("Quadratic distance identities failed")
    joint = np.zeros((8,8)); joint[:4,:4]=cs; joint[4:,4:]=cr
    return {"schema": "matching-one.p267-response-ray.v1", "source": source,
            "analysis_role": manifest["analysis_role"], "dependency": manifest["dependency"],
            "field_order": fields, "x_tau_i": x.tolist(), "y_tau_2i": y.tolist(),
            "joint_covariance_order": ["tau_i:"+f for f in fields]+["tau_2i:"+f for f in fields],
            "joint_covariance": joint.tolist(), "full_ray": full,
            "coordinate_diagnostics": diagnostics, "Schur_diagnostics_at_full_lambda": partitions,
            "amplitude_improvement": {"chi_square_difference": full["lambda_one_chi_square"]-full["chi_square"], "nominal_df": 1},
            "checks": checks, "future_geometry": future_geometry(),
            "scope": manifest["scope"], "inference": manifest["inference"],
            "environment": {"python": platform.python_version(), "machine": platform.machine(), "numpy": np.__version__, "scipy": scipy.__version__},
            "script_sha256": sha(Path(__file__).read_bytes()), "elapsed_seconds": time.perf_counter()-start,
            "generated_at_utc": datetime.now(timezone.utc).isoformat()}


def report(r):
    f = r["full_ray"]
    lines = ["# P267: changing amplitude does not close the four-coordinate response", "",
             f"The best common amplitude is `{f['lambda']:.7g}`, but the two fixed-p response vectors still disagree at nominal `chi2={f['chi_square']:.7g}/3, p={f['nominal_p']:.7g}`.",
             "This excludes the declared finite four-coordinate amplitude-only model. It does not count microscopic sources or identify a continuum operator.", "",
             "## Definitions and evidence", "",
             "x and y are the exact cos4-normalized `(A_top,E_top,C,W)` contrasts at tau=i and 2i, fixed N50 and p=.59274605079. C and W are normalized integrated clocks. The old 12M square block and new 100k rectangular block are independent; each retains its internally paired direction covariance.", "",
             "| coordinate subset | fitted amplitude | residual chi2 / nominal df | nominal p |", "|---|---:|---:|---:|",
             f"| A/E/C/W | {f['lambda']:.7g} | {f['chi_square']:.7g} / 3 | {f['nominal_p']:.7g} |"]
    for row in r["coordinate_diagnostics"]:
        lines.append(f"| {'/'.join(row['fields'])} | {row['lambda']:.7g} | {row['chi_square']:.7g} / {row['nominal_df']} | {row['nominal_p']:.7g} |")
    lines += ["", "These subsets localize one shared result, not independent discoveries. Their nominal p-values are retrospective and not multiplicity-adjusted. The source vectors and complete joint covariance are saved in score.json.", "",
              "## Model and uncertainty", "",
              "The joint Gaussian model is `x~N(mu,Cs), y~N(lambda*mu,Cr)`. Profiling mu gives `D(lambda)=(y-lambda*x)^T(Cr+lambda^2*Cs)^-1(y-lambda*x)`. The original joint covariance determinant is constant, so no contrast log-determinant is added. A compact projective-angle search includes lambda=infinity and refines the grid minima.", "",
              f"Fixed lambda=1 gives `{f['lambda_one_chi_square']:.8g}` on nominal df4. Allowing one amplitude improves distance by `{r['amplitude_improvement']['chi_square_difference']:.8g}` on nominal df1, but leaves `{f['chi_square']:.8g}` on df3. The infinity boundary has distance `{f['projective_infinity_chi_square']:.8g}`.", "",
              f"The 95% conditional amplitude profile interval is `{f['profile_intervals'][0]['lambda']}`. This is a parameter range inside a rejected model, not an identified physical transport coefficient. All inference conditions on the saved estimated covariance; no exact finite-sample coverage is claimed.", "",
              "## Where the four-coordinate mismatch remains", ""]
    for p in r["Schur_diagnostics_at_full_lambda"]:
        first='/'.join(r['field_order'][i] for i in p['first_indices'])
        second='/'.join(r['field_order'][i] for i in p['conditional_indices'])
        lines.append(f"At the full fitted amplitude, {first} marginal distance is `{p['marginal_distance']:.7g}`; {second} conditional on it contributes `{p['conditional_distance']:.7g}`. Their sum is the full distance.")
    lines += ["", "These order-dependent Schur partitions describe correlated residuals; they are not causal shares or separate tests. Two observed vectors always fit some rank-two plane, so a saturated rank-two fit would add no mechanism information.", "",
              "## A third modulus needs a realizable geometry", "",
              "At N50, Smith(5,10) forces L=5B with det B=2. Up to square-lattice symmetry, the index-two sublattices provide only the square and 1x2 rectangle. A third same-N, same-Smith modulus is therefore not available. A new map must change a concrete endpoint/readout, not merely rename a Smith class.", "",
              "The saved exact future design uses N250 and Smith(1,250)/(5,50), with the same O=(1/5)[[4,-3],[3,4]] at rho=1,2,5. Its integer matrices, determinants, orthogonal columns and exact direction factors are recorded in score.json. It is not acquired, powered or a request to start production.", "",
              "Two explicit exploratory shape candidates are `f_log(rho)=(log(rho)/log(2))^2` and `f_geom(rho)=2(rho+1/rho-2)`. Both interpolate the two anchors and are even under axis exchange. At rho=5 they predict coefficients `f_log="+f"{(math.log(5)/math.log(2))**2:.8g}"+"` versus `f_geom=6.4` in `u5=(1-f5)u1+f5*u2`. These are distinguishable phenomenological hypotheses, not laws selected by the old two points.", "",
              "Any N250 prediction needs N250 anchors or an explicitly justified coordinate-wise scale law. Current N50 vectors are not numerical predictions at N250. Likewise, fixed-p rows must not be concatenated with intrinsic-center norm-4/norm-5 production rows without recomputing the same observable definition and covariance.", "",
              "## Next scientific output", "",
              "Preserve the amplitude-only exclusion. Use the coordinate localization to define the next physical response model; compare its explicit predictions on a realizable design. Do not treat an arbitrary rank-two interpolation, a free surface, another compiler or more N50 factorial replicas as identification.", "",
              "## Reproduce and provenance", "",
              f"Input: `{r['source']['commit']}:{r['source']['path']}`; SHA256 `{r['source']['sha256']}`. This is a dependent, post-reveal reanalysis of the existing factorial. No Monte Carlo or test suite was run.", "",
              "```sh", "python3 scripts/analyze_p267_response_ray.py", "```", "",
              "One in-run arithmetic check recovers the parent's fixed-amplitude distance and both Schur sums. No chart is needed: the small exact-value table exposes the coordinate comparison directly. All raw covariance, fitted residuals, profile intervals and prospective integer designs are in the companion JSON."]
    return "\n".join(lines)+"\n"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, default=ROOT/"analysis/p267_response_ray_manifest.json")
    p.add_argument("--output-dir", type=Path, default=ROOT/"results/p267-response-ray")
    a = p.parse_args()
    manifest = json.loads(a.manifest.read_text())
    result = analyze(manifest)
    a.output_dir.mkdir(parents=True, exist_ok=True)
    (a.output_dir/"score.json").write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
    (a.output_dir/"REPORT.md").write_text(report(result))
    print(json.dumps({"lambda": result['full_ray']['lambda'], "chi_square": result['full_ray']['chi_square'], "nominal_p": result['full_ray']['nominal_p'], "checks": result['checks'], "elapsed_seconds": result['elapsed_seconds']}))


if __name__ == "__main__":
    main()
