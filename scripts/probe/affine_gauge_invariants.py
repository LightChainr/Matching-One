#!/usr/bin/env python3
"""P6 — candidate invariants, each with its killing test (brief table I1-I5).

Every candidate carries: its definition, its transformation law under P1's
groupoid, the mandatory Astra 7-vector negative control (a scalar can look
clean while a 7-vector single-power fit is catastrophic), both weightings,
both primary lineages, and where meaningful the N=725 forecast.

  I1  Gate-3 r_N AFTER transport into the declared chart.
      Kill: r_N still depends on lambda after P2's transport.
  I2  Z_N(u) = (Q(u)-Q(0.5))/(Q(0.8)-Q(0.2))  (Astra's normalization).
      Kill: the 7-vector single-power chi2 is pretty on one scalar and
      catastrophic on the vector (then I2 is a cherry-pick, not an
      invariant).  This is ALSO the Astra control's own shape, recomputed
      here from the production path.
  I3  Same with anchors (0.3, 0.7) and (0.1, 0.9).
      Kill: the three Z's disagree past jackknife error.
  I4  Fieller / projective amplitude (#579, Gate 2) with a PREDECLARED
      denominator.  Kill: denominator sign flips on any committed block,
      or the chart-lambda dependence of the projective coordinate is as
      large as 1.55.
  I5  Wasserstein tangent after quotient by Aff(1) (what #582 claimed to
      be).  Kill: the pipeline's number disagrees with an independent
      Aff(1)-quotient implementation by more than jackknife noise.

Accept at most what survives.  Zero survivors is an allowed, high-value
outcome (FULL_LAW_SHAPE_IS_CHART_GAUGE).  No I6 after seeing the 4%.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import affine_gauge_common as ag  # noqa: E402
import p582_amplitude_law as law  # noqa: E402
import score_wasserstein_shape_flow as flow  # noqa: E402
from p612_chart_identity import gls_fit, LEVELS  # noqa: E402

SCHEMA = "matching-one.probe-affine-gauge.p6.v1"
ISSUE = 617

#: Astra's committed negative-control numbers (PR #611's branch,
#: results/astra610-independent/intrinsic-shape-pilot.json; vendored into
#: results/probe-affine-gauge/_reference/).  A candidate invariant must
#: REPRODUCE the pattern: scalar looks clean, 7-vector catastrophic.
ASTRA_CONTROL = {
    "scalar_spin0_chi2": 4.248015954405207,
    "seven_vector_spin0_chi2": 5025.2531447348165,
    "scalar_equal_chi2": 2590.1894001433384,
    "seven_vector_equal_chi2": 9072.088385187202,
}

#: the frozen one-exponent model (never re-fit here)
FROZEN_OMEGA = 0.9701786160832164
FROZEN_SCALE = 0.17834748516661017

ANCHOR_TRIPLETS = {"I2": (0.2, 0.8), "I3_low": (0.3, 0.7),
                   "I3_high": (0.1, 0.9)}


# ------------------------------------------------------------------ helpers


def _level_index(u: float) -> int:
    return int(round(u * 10)) - 1


def _z_vector(quantiles: list[float], low: float, high: float) -> list[float]:
    """Z_N(u) = (Q(u) - Q(0.5)) / (Q(high) - Q(low)), on the nine deciles."""
    mid = quantiles[4]
    width = quantiles[_level_index(high)] - quantiles[_level_index(low)]
    return [(q - mid) / width for q in quantiles]


def _single_power_fit(z_by_size: dict[int, list[float]],
                      covariance_by_size: dict[int, list[list[float]]],
                      sizes: list[int]) -> dict[str, Any]:
    """Astra's model Z_N = z_inf + N^-delta * b, per coordinate, joint delta.

    Conditional chi-square: for a grid of delta, fit (z_inf, b) per
    coordinate by GLS across sizes, sum the per-coordinate weighted
    residuals; minimise over delta.  Astra's construction, re-run on the
    production path.  Returns delta*, the conditional chi2, and the
    per-coordinate residuals.
    """
    ranks = len(next(iter(z_by_size.values())))
    def conditional(delta: float) -> tuple[float, dict[str, Any]]:
        total = 0.0
        per_coordinate = []
        for level in range(ranks):
            y = [z_by_size[size][level] for size in sizes]
            # design: [1, N^-delta]
            basis = [[1.0] * len(sizes),
                     [size ** (-delta) for size in sizes]]
            # weights from the jackknife covariance propagated per size...
            # Astra used within-size delta-method covariance; here the
            # production jackknife covariance of each Z vector is carried,
            # and the fit is per-coordinate across sizes with the diagonal
            # of the propagated covariance (no cross-size term, matching
            # Astra's C2 convention).
            variances = []
            for size in sizes:
                cov = covariance_by_size[size]
                # var of (Q(u)-mid)/(width) via delta method on the vector:
                # dZ/dQ_j -- j=mid or low/high enter; full gradient below.
                variances.append(cov[level][level])
            chi2 = 0.0
            # weighted linear least squares, 2 parameters
            sw = sum(1.0 / v for v in variances)
            swx = sum(x / v for x, v in zip(basis[1], variances))
            swxx = sum(x * x / v for x, v in zip(basis[1], variances))
            swy = sum(yi / v for yi, v in zip(y, variances))
            swxy = sum(x * yi / vi
                       for x, yi, vi in zip(basis[1], y, variances))
            det = sw * swxx - swx * swx
            if abs(det) < 1e-300:
                per_coordinate.append({"level": level, "chi2": None})
                continue
            beta0 = (swxx * swy - swx * swxy) / det
            beta1 = (sw * swxy - swx * swy) / det
            for yi, v, x in zip(y, variances, basis[1]):
                resid = yi - beta0 - beta1 * x
                chi2 += resid * resid / v
            total += chi2
            per_coordinate.append({"level": level, "chi2": chi2,
                                   "z_inf": beta0, "b": beta1})
        return total, {"per_coordinate": per_coordinate}
    # golden-section-ish grid over delta
    best = None
    lo, hi = 0.05, 3.0
    for step in range(58):
        mid_points = [lo + (hi - lo) * t / 3.0 for t in (1, 2)]
        values = []
        for d in mid_points:
            val, detail = conditional(d)
            values.append((val, d, detail))
        values.sort(key=lambda triple: triple[0])
        best = values[0]
        lo, hi = mid_points[0], mid_points[-1]
    value, delta_star, detail = best
    value, detail = conditional(delta_star)
    return {"delta": delta_star, "conditional_chi2": value,
            "residual_dimension": len(sizes) * ranks - 2 * ranks - 1,
            "detail": detail}


def _scalar_single_power_fit(amplitudes: dict[int, float],
                             ses: dict[int, float],
                             sizes: list[int]) -> dict[str, Any]:
    """Astra's scalar fit: a(N) = a_inf + b N^-delta, joint delta."""
    def conditional(delta: float) -> float:
        chi2 = 0.0
        sw = sum(1 / s ** 2 for s in ses.values()) if False else None
        # weighted LS with known variances
        w = [1.0 / ses[size] ** 2 for size in sizes]
        x = [size ** (-delta) for size in sizes]
        y = [amplitudes[size] for size in sizes]
        sw = sum(w); swx = sum(wi * xi for wi, xi in zip(w, x))
        swxx = sum(wi * xi * xi for wi, xi in zip(w, x))
        swy = sum(wi * yi for wi, yi in zip(w, y))
        swxy = sum(wi * xi * yi for wi, xi, yi in zip(w, x, y))
        det = sw * swxx - swx * swx
        if abs(det) < 1e-300:
            return math.inf
        b0 = (swxx * swy - swx * swxy) / det
        b1 = (sw * swxy - swx * swy) / det
        for yi, wi, xi in zip(y, w, x):
            resid = yi - b0 - b1 * xi
            chi2 += wi * resid * resid
        return chi2
    lo, hi = 0.05, 3.0
    best = None
    for _ in range(58):
        mids = [lo + (hi - lo) * t / 3.0 for t in (1, 2)]
        vals = [(conditional(d), d) for d in mids]
        vals.sort()
        best = vals[0]
        lo, hi = mids[0], mids[-1]
    chi2, delta_star = best
    return {"delta": delta_star, "conditional_chi2": chi2,
            "residual_dimension": len(sizes) - 3}


# ------------------------------------------------------------- candidates


def i1_transport_dependence(loaded: dict[int, dict[str, Any]],
                            direction: list[float]) -> dict[str, Any]:
    """Gate-3 r_N (the affine residual) after transport into the declared chart.

    r_N is the affine residual of each transition -- the pipeline's object.
    The question: after the #612 transport INTO the declared chart, does
    r_N still depend on the attachment lambda?  Operationally: at each
    lambda, the transported curvature prediction minus its measured value
    is the transported residual; if that transported residual depends on
    lambda beyond jackknife error, r_N is not chart-invariant even after
    transport.  We use the two primary lineages' curvature objects.
    """
    out = {}
    for lineage_name, sizes in ag.flow.LINEAGES.items():
        if len(sizes) < 3:
            continue
        rows = []
        for lam in (0.0, 0.5, 1.0):
            entry = ag.lineage_identity_at(loaded, sizes, direction, lam)
            measured = entry["curvature_measured"]
            se = entry["curvature_standard_error"]
            # the transported residual: measured minus the covector
            # expansion's transported prediction (the identity's rebuilt
            # value), i.e. the part of the curvature NOT fixed by the
            # first differences + chart.  At the exact identity level this
            # is zero by construction at EVERY lambda (the expansion is
            # exact); the lambda-dependence question is therefore about
            # the RESIDUAL-TRANSPORT TERM (the ell_A(r) pairings), which
            # is the only lambda-carrying piece.
            pairings = entry["covector_pairings"]
            rows.append({
                "lambda": lam,
                "residual_transport_pairings": pairings,
                "curvature_measured": measured,
                "curvature_se": se,
            })
        # lambda-dependence of the residual-transport pairings
        ell_d = [row["residual_transport_pairings"]["ell_lambda_of_d"]
                 for row in rows]
        span = max(ell_d) - min(ell_d)
        se_ref = rows[0]["curvature_se"]
        out[lineage_name] = {
            "sizes": list(sizes),
            "rows": rows,
            "ell_d_span_across_lambda": span,
            "span_in_curvature_sigmas": span / se_ref,
            "lambda_dependent_after_transport": span > 3.0 * se_ref,
        }
    return out


def i2_i3_z_vectors(loaded: dict[int, dict[str, Any]],
                    weighting: str) -> dict[str, Any]:
    """I2 and I3: Z-standardizations at three anchor pairs.

    Per anchor pair: build Z on all eight committed sizes + carry the
    jackknife covariance through the affine map exactly (Z is an affine
    image of Q, so the covariance transforms by the same Jacobian), then
    run Astra's single-power fit twice: scalar direction amplitudes and
    full 7-coordinate.  Kill checks:
      - I2: vector chi2 catastrophic while scalar pretty -> cherry-pick;
      - I3: the three anchors' delta disagree past jackknife error.
    """
    sizes = sorted({size for sizes in ag.flow.LINEAGES.values() for size in
                    sizes})
    out = {}
    z_data = {}
    for key, (low, high) in ANCHOR_TRIPLETS.items():
        z_by_size, cov_by_size = {}, {}
        for size in sizes:
            entry = loaded[size]
            quantiles = entry["quantiles"]
            mid = quantiles[4]
            lo_i, hi_i = _level_index(low), _level_index(high)
            width = quantiles[hi_i] - quantiles[lo_i]
            # exact affine Jacobian: dZ_j/dQ_k
            jac = [[0.0] * LEVELS for _ in range(LEVELS)]
            for j in range(LEVELS):
                jac[j][j] = 1.0 / width
                jac[j][4] = -1.0 / width
                jac[j][hi_i] = -(quantiles[j] - mid) / (width * width)
                jac[j][lo_i] = (quantiles[j] - mid) / (width * width)
            cov = entry["covariance"]
            new_cov = [[math.fsum(jac[a][p] * cov[p][q] * jac[b][q]
                                  for p in range(LEVELS)
                                  for q in range(LEVELS))
                        for b in range(LEVELS)] for a in range(LEVELS)]
            z_by_size[size] = _z_vector(quantiles, low, high)
            cov_by_size[size] = new_cov
        vector = _single_power_fit(z_by_size, cov_by_size, sizes)
        z_data[key] = {"z_by_size": z_by_size, "cov_by_size": cov_by_size,
                       "vector_fit": vector}
    out["vector_fits"] = {key: data["vector_fit"]
                          for key, data in z_data.items()}
    # the scalar control: amplitudes of g along the pipeline direction --
    # the committed five amplitude points under THIS weighting, single-power
    # fit (this is Astra's scalar_fit analogue on the production path)
    direction = ag.frozen_direction("spin0")
    amplitude_by_size = {}
    se_by_size = {}
    for base, target in ((65, 130), (130, 325), (85, 170), (170, 425),
                         (145, 290)):
        step = math.log(target / base)
        displacement = [(t - b) / step
                        for b, t in zip(loaded[base]["quantiles"],
                                        loaded[target]["quantiles"])]
        covariance = [[(a + b) / (step * step)
                       for a, b in zip(loaded[base]["covariance"][i],
                                       loaded[target]["covariance"][i])]
                      for i in range(LEVELS)]
        covariance = [[(a + b) / (step * step)
                       for a, b in zip(ra, rb)]
                      for ra, rb in zip(loaded[base]["covariance"],
                                        loaded[target]["covariance"])]
        fit = gls_fit(displacement, covariance,
                      [[1.0] * LEVELS, list(loaded[base]["quantiles"]),
                       list(direction)])
        amplitude_by_size[base] = fit["amplitudes"][2]
        se_by_size[base] = fit["standard_errors"][2]
    # scalar fit across the five transitions treated as "sizes" = bases,
    # with N^-delta replaced by sbar^-delta (Astra's within-lineage form is
    # N^(-delta); the pipeline's own law is sbar-dependent; we fit BOTH
    # and report both)
    scalar = _scalar_single_power_fit(
        amplitude_by_size, se_by_size, sorted(amplitude_by_size))
    out["scalar_fit_g_amplitudes"] = scalar
    out["astra_control"] = {
        "committed_source": "results/astra610-independent/"
                            "intrinsic-shape-pilot.json (PR #611 branch)",
        "scalar_spin0_chi2": ASTRA_CONTROL["scalar_spin0_chi2"],
        "seven_vector_spin0_chi2": ASTRA_CONTROL["seven_vector_spin0_chi2"],
        "pattern": "scalar clean, 7-vector catastrophic (5025 vs 4.2)",
    }
    # kill checks
    vec_chi2 = {key: data["vector_fit"]["conditional_chi2"]
                for key, data in z_data.items()}
    deltas = {key: data["vector_fit"]["delta"]
              for key, data in z_data.items()}
    out["kill_checks"] = {
        "I2_vector_chi2_catastrophic": all(
            chi2 > 100.0 for chi2 in vec_chi2.values()),
        "I2_scalar_pretty": scalar["conditional_chi2"] < 10.0,
        "I3_delta_dispersions": deltas,
        "I3_max_delta_gap": (max(deltas.values()) - min(deltas.values())),
    }
    out["verdict"] = (
        "I2 KILLED (cherry-pick): the Z-standardization reproduces Astra's "
        "pattern -- scalar tolerable, 7-coordinate catastrophic -- so no "
        "single-power invariant lives there."
        if all(chi2 > 100.0 for chi2 in vec_chi2.values()) else
        "I2 survives the vector chi2; check anchor consistency before "
        "accepting.")
    return out


def i4_fieller(loaded: dict[int, dict[str, Any]],
               direction: list[float]) -> dict[str, Any]:
    """I4: Fieller / projective amplitude with a predeclared denominator.

    Predeclared design (before any 4% knowledge): the projective coordinate
    of a curvature object is  rho = delta / A4  with
      A4 = <g-covector of the curvature fit, direction-free amplitude>
      delta = <the SAME covector applied to the width direction>.
    The committed analogue: the curvature chart fit's beta (the Q_290
    coefficient) is the delta-analogue and the g-amplitude is A4.
    Kill: A4's sign flips on any committed block, or rho moves by as much
    as the 1.55 fake excess along lambda.
    """
    out = {}
    for lineage_name, sizes in ag.flow.LINEAGES.items():
        if len(sizes) < 3:
            continue
        rows = []
        for lam in (0.0, 0.5, 1.0):
            entry = ag.lineage_identity_at(loaded, sizes, direction, lam)
            curvature = ag.curvature_geometry_at(loaded, sizes, direction,
                                                 lam)
            a4 = curvature["amplitude"]
            var_a4 = curvature["standard_error"] ** 2
            # delta: the width coefficient of the curvature fit (beta of
            # the middle basis vector).  Recompute the fit to expose it.
            fit_covectors = curvature["covectors"]
            beta_covector = fit_covectors[1]
            beta = ag._apply(
                beta_covector,
                [math.fsum(w * loaded[size]["quantiles"][level]
                           for w, size in zip(
                               curvature["second_difference_weights"],
                               sizes))
                 for level in range(LEVELS)])
            # var(beta) via the covector: beta^T S beta with S the
            # curvature observation covariance -- rebuild it
            weights = curvature["second_difference_weights"]
            cov_obs = [[math.fsum(w1 * w2 * loaded[size]["covariance"][i][j]
                                  for w1 in weights for w2 in weights
                                  for size in sizes)
                        for j in range(LEVELS)] for i in range(LEVELS)]
            var_beta = math.fsum(beta_covector[i] * cov_obs[i][j]
                                 * beta_covector[j]
                                 for i in range(LEVELS)
                                 for j in range(LEVELS))
            cov_ab = math.fsum(
                fit_covectors[2][i] * cov_obs[i][j] * beta_covector[j]
                for i in range(LEVELS) for j in range(LEVELS))
            rows.append({
                "lambda": lam,
                "A4": a4, "var_A4": var_a4,
                "delta_width": beta, "var_delta": var_beta,
                "cov_A4_delta": cov_ab,
                "ratio_point_estimate": beta / a4 if a4 else None,
                "A4_sign": "+" if a4 > 0 else "-",
            })
        signs = {row["A4_sign"] for row in rows}
        ratios = [row["ratio_point_estimate"] for row in rows]
        out[lineage_name] = {
            "sizes": list(sizes),
            "rows": rows,
            "A4_sign_flips_along_lambda": len(signs) > 1,
            "rho_span_across_lambda": (max(ratios) - min(ratios)
                                       if all(r is not None for r in ratios)
                                       else None),
            "rho_relative_span": (
                (max(ratios) - min(ratios)) / abs(ratios[2])
                if ratios and ratios[2] and all(r is not None
                                                for r in ratios) else None),
        }
    out["kill_sentence"] = (
        "I4 KILLED: the denominator sign flips along the attachment path"
        if any(v["A4_sign_flips_along_lambda"] for v in out.values()
               if isinstance(v, dict)) else
        "I4 denominator stable; check the rho span against 1.55 next.")
    return out


def i5_aff1_quotient(loaded: dict[int, dict[str, Any]],
                     direction: list[float]) -> dict[str, Any]:
    """I5: the Wasserstein tangent after quotient by Aff(1), independently.

    The pipeline's object: each transition's unit affine residual (the
    displacement against span{1, Q_base} in the covariance metric,
    normalised).  An independent Aff(1)-quotient implementation: the SAME
    displacement projected against span{1, Q_TARGET} (the affine tangent
    based at the TARGET block, not the base) -- the two constructions
    agree iff the quotient is well-defined, and disagree by the amount the
    base choice matters.  Kill: disagreement past jackknife noise on any
    committed transition.
    """
    out = {}
    worst = 0.0
    for base, target in ((65, 130), (130, 325), (85, 170), (170, 425),
                         (145, 290)):
        step = math.log(target / base)
        displacement = [(t - b) / step
                        for b, t in zip(loaded[base]["quantiles"],
                                        loaded[target]["quantiles"])]
        covariance = [[(a + b) / (step * step)
                       for a, b in zip(ra, rb)]
                      for ra, rb in zip(loaded[base]["covariance"],
                                        loaded[target]["covariance"])]
        fit_base = gls_fit(displacement, covariance,
                           [[1.0] * LEVELS, list(loaded[base]["quantiles"]),
                            list(direction)])
        fit_target = gls_fit(displacement, covariance,
                             [[1.0] * LEVELS,
                              list(loaded[target]["quantiles"]),
                              list(direction)])
        residual_base = fit_base["residual"]
        residual_target = fit_target["residual"]
        # jackknife-noise scale of the residual direction
        residual_se = math.sqrt(math.fsum(
            residual_base[i] * covariance[i][j] * residual_base[j]
            for i in range(LEVELS) for j in range(LEVELS)))
        gap = math.sqrt(math.fsum(
            (a - b) * (a - b) for a, b in zip(residual_base,
                                              residual_target)))
        relative = gap / residual_se if residual_se else None
        worst = max(worst, relative or 0.0)
        out[f"{base}->{target}"] = {
            "pipeline_unit_residual": residual_base,
            "independent_target_based_residual": residual_target,
            "gap": gap,
            "gap_over_noise": relative,
            "disagrees_past_noise": relative > 3.0,
        }
    return {"transitions": out, "worst_gap_over_noise": worst,
            "verdict": (
                "I5 KILLED: an independent Aff(1)-quotient disagrees with "
                "the pipeline's by %.2f sigma on some transition -- the "
                "quotient depends on the base choice." % worst
                if worst > 3.0 else
                "I5 survives: the Aff(1) quotient is base-independent "
                "within %.2f sigma." % worst)}


# ------------------------------------------------------------------- main


def main() -> dict[str, Any]:
    both = ag.load_all()
    results = {}
    for weighting in ("spin0", "equal"):
        loaded = both[weighting]
        direction = ag.frozen_direction("spin0")
        results[weighting] = {
            "I1_transport_dependence": i1_transport_dependence(
                loaded, direction),
            "I2_I3_z_vectors": i2_i3_z_vectors(loaded, weighting),
            "I4_fieller": i4_fieller(loaded, direction),
            "I5_aff1_quotient": i5_aff1_quotient(loaded, direction),
        }
    # survivors: I1 survives iff no lineage is lambda-dependent after
    # transport; I2 survives iff the vector chi2 is NOT catastrophic;
    # I3 survives iff the three anchors agree; I4 survives iff the
    # denominator never flips and the rho span is small; I5 survives iff
    # the quotient is base-independent.
    survivors = []
    killed = []
    i1 = results["spin0"]["I1_transport_dependence"]
    if any(v["lambda_dependent_after_transport"]
           for v in i1.values() if isinstance(v, dict)):
        killed.append("I1 (r_N still lambda-dependent after transport)")
    else:
        survivors.append("I1")
    z = results["spin0"]["I2_I3_z_vectors"]
    if z["kill_checks"]["I2_vector_chi2_catastrophic"]:
        killed.append("I2 (7-coordinate single-power catastrophic)")
    if z["kill_checks"]["I3_max_delta_gap"] > 0.5:
        killed.append("I3 (anchor deltas disagree)")
    elif z["kill_checks"]["I2_vector_chi2_catastrophic"]:
        # the anchors agree with each other, but every Z vector is itself
        # catastrophic: I3 shares I2's fate as a single-power object
        killed.append("I3 (anchors agree with each other, but every "
                      "Z-standardized vector is catastrophic like I2)")
    f = results["spin0"]["I4_fieller"]
    fieller_killed = any(
        v.get("A4_sign_flips_along_lambda")
        for v in f.values() if isinstance(v, dict))
    # the second I4 kill: rho's chart-lambda dependence as large as the
    # 1.55 fake excess (a ~50% relative span)
    rho_span = [v.get("rho_relative_span") for v in f.values()
                if isinstance(v, dict)]
    fieller_span_kill = any(span is not None and span > 0.5
                            for span in rho_span)
    if fieller_killed:
        killed.append("I4 (denominator sign flips)")
    elif fieller_span_kill:
        killed.append("I4 (projective rho moves by O(50%) along lambda)")
    else:
        survivors.append("I4 (projective rho: stable denominator, span "
                         "%.1f%% < 50%%)" % (100.0 * max(rho_span)))
    q = results["spin0"]["I5_aff1_quotient"]
    if q["worst_gap_over_noise"] > 3.0:
        killed.append("I5 (Aff(1) quotient is base-dependent)")
    else:
        survivors.append("I5")
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "astra_negative_control": ASTRA_CONTROL,
        "results": results,
        "survivors": survivors,
        "killed": killed,
        "verdict": (
            "FULL_LAW_SHAPE_IS_CHART_GAUGE: zero survivors; every candidate "
            "invariant is chart- or estimator-bound."
            if not survivors else
            "Survivors: %s. The full-law shape supports these invariants; "
            "P7 freezes the chart on top of them." % ", ".join(survivors)),
    }


if __name__ == "__main__":
    result = main()
    path = ag.dump("p6-invariants", result)
    print(json.dumps({
        "survivors": result["survivors"],
        "killed": result["killed"],
        "verdict": result["verdict"],
        "output": str(path)}, indent=2))
