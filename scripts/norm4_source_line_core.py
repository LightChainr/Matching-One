"""Conditional primitive-line O4 responses on already saved matching roots.

Pure numerical core: no file I/O, root finder, simulation or uncertainty driver.
The caller keeps the original paired/source-group delete-one alignment and
updates p0, source root motion and the E loading in the same omitted sample.

O4 is the sign/scale-invariant physical-lift quantity
    (dx + i*dy)^4 / (dx^2 + dy^2)^2
inside ambient rank one.  Both complex coordinates are retained.  Source s is
the BULK black-plus-white matching cluster count, without division by N.
"""
from __future__ import annotations

import itertools
import math
from fractions import Fraction

import numpy as np

from analyze_norm4_source_thermal import binomial_moments

INPUT_COLUMNS = (
    "sum_rank1", "sum_rank1_s", "sum_line4_re", "sum_line4_im",
    "sum_line4_s_re", "sum_line4_s_im",
)
COMPLEX_COMPONENTS = ("first_re", "first_im", "second_re", "second_im")
DIRECTIONS = ("first", "second")
QUANTITIES = ("mu", "mu_p", "conditional_cov_s", "nu_s", "pred_E", "resid_E")


def _real_vector(values):
    values = np.asarray(values, dtype=complex)
    if values.shape != (2,):
        raise ValueError("a complex conditional vector must contain both geometries")
    return np.asarray([values[0].real, values[0].imag, values[1].real, values[1].imag])


def common_clock_wedges(conditional_cov_s, mu_p):
    """All six unnormalized minors of (conditional source covariance, mu_p).

    For X=a_g+b*K+c*q+d*E with b,c,d real and geometry-common, q=E=0
    within rank one, hence
      Cov(O4,X | rank1,g) = b*p*(1-p)*d_p E[O4 | rank1,g].
    The four-real covariance vector and mu_p vector must be collinear with
    the same real multiplier.  Every minor is necessary; none divides by a
    small slope.  This is a conditional-response null, not a claim that the
    microscopic source functions are pointwise identical.  If mu_p vanishes,
    these minors alone lose power, so the raw covariance vector is retained.
    """
    covariance = _real_vector(conditional_cov_s)
    slope = _real_vector(mu_p)
    if not np.isfinite(covariance).all() or not np.isfinite(slope).all():
        raise ValueError("conditional source covariance and slopes must be finite")
    return {
        f"clock_wedge.{COMPLEX_COMPONENTS[i]}__{COMPLEX_COMPONENTS[j]}":
        float(covariance[i] * slope[j] - covariance[j] * slope[i])
        for i, j in itertools.combinations(range(4), 2)
    }


def evaluate(sums, samples, n, p0, raw_scalar, baseline_direction):
    """Return (flat_values, diagnostics) at one supplied central/delete-one root.

    sums has shape (2,N+1,6), holding aggregate sums in INPUT_COLUMNS order.
    samples is the retained number of permutations PER geometry, at every K.
    Each original permutation supplies all K; K rows are not independent.

    raw_scalar uses unprefixed saved scalar keys:
      rootdot_fugacity, root_comoving_rank1_fugacity, rank1_common_E.
    The first two are derivatives under exp(t*s); rank1_common_E uses exp(eta*E)
    and has NO factor N.  baseline_direction is old direction_values' mapping
    first/second -> q,E,q_p,..., evaluated on the same saved root/sample state.

    The E+clock coefficient is alpha_E = r_s/r_E, common to both geometries.
    Since E=0 in rank one, its conditional covariance is zero; its conditional
    moving-root response is p0dot_E*mu_p.  The prediction and residual therefore
    require neither a new fit nor source/E mixed moments.
    """
    profiles = np.asarray(sums, dtype=float)
    samples, p0 = float(samples), float(p0)
    if profiles.shape != (2, n + 1, 6) or not np.isfinite(profiles).all():
        raise ValueError("expected finite (2,N+1,6) bulk-source conditional sum profiles")
    if samples <= 0 or not math.isfinite(samples) or not 0 < p0 < 1:
        raise ValueError("positive retained sample count and a supplied root in (0,1) required")
    if n in (65, 130, 260):
        delta_exact = Fraction("1152/845")
    elif n in (85, 170, 340):
        delta_exact = Fraction("2304/1445")
    else:
        raise ValueError("only the six frozen norm-4 lineage sizes are supported")
    delta = float(delta_exact)

    rootdot_s = float(raw_scalar["rootdot_fugacity"])
    rank1dot_s = float(raw_scalar["root_comoving_rank1_fugacity"])
    rank1dot_e = float(raw_scalar["rank1_common_E"])
    rows = [baseline_direction[direction] for direction in DIRECTIONS]
    m = np.asarray([row["q"] for row in rows], dtype=float)
    e = np.asarray([row["E"] for row in rows], dtype=float)
    matching_slope = float(np.mean([row["q_p"] for row in rows]))
    if (not all(map(math.isfinite, (rootdot_s, rank1dot_s, rank1dot_e, matching_slope)))
            or matching_slope <= 0 or not np.isfinite(m).all() or not np.isfinite(e).all()):
        raise ValueError("saved scalar responses and matching jets must be finite with positive slope")
    if rank1dot_e == 0:
        raise ValueError("E-source rank1 loading is undefined because its composition response is zero")
    alpha_e = rank1dot_s / rank1dot_e
    rootdot_e = float(-np.mean(m * (1 - e)) / matching_slope)

    packets = [binomial_moments(profiles[g], samples, p0, n) for g in range(2)]
    quantities = {name: np.zeros(2, dtype=complex) for name in QUANTITIES}
    values = {
        "p0": p0, "alpha_E": float(alpha_e), "rootdot_s": rootdot_s,
        "rootdot_E": rootdot_e, "rank1rootdot_s": rank1dot_s,
        "rank1rootdot_E": rank1dot_e,
    }
    direction_diagnostics = {}
    for g, direction in enumerate(DIRECTIONS):
        means, first, second, integration = packets[g]
        a = float(means[0])
        if not math.isfinite(a) or a <= 0:
            raise ValueError(f"N{n} {direction}: conditional rank1 probability is not positive")
        b = complex(means[2], means[3])
        c = float(means[1])
        ds = complex(means[4], means[5])
        ap = float(first[0])
        bp = complex(first[2], first[3])
        mu = b / a
        mu_p = (bp - mu * ap) / a
        conditional_source_mean = c / a
        conditional_cov_s = ds / a - mu * conditional_source_mean
        nu_s = conditional_cov_s + rootdot_s * mu_p
        pred_e = alpha_e * rootdot_e * mu_p
        resid_e = nu_s - pred_e
        for name, value in zip(QUANTITIES, (mu, mu_p, conditional_cov_s, nu_s, pred_e, resid_e)):
            quantities[name][g] = value
            values[f"{direction}.{name}_re"] = float(value.real)
            values[f"{direction}.{name}_im"] = float(value.imag)
        values[f"{direction}.rank1_probability"] = a
        values[f"{direction}.rank1_probability_p"] = ap
        values[f"{direction}.conditional_source_mean_bulk"] = conditional_source_mean
        direction_diagnostics[direction] = {
            "A": a, "A_p": ap,
            "B_re": float(b.real), "B_im": float(b.imag),
            "B_p_re": float(bp.real), "B_p_im": float(bp.imag),
            "C_bulk": c, "D_bulk_re": float(ds.real), "D_bulk_im": float(ds.imag),
            "raw_means": means.tolist(), "raw_p_derivatives": first.tolist(),
            "raw_p_second_derivatives": second.tolist(), "integration": integration,
        }

    # Condition in each geometry BEFORE either direction subtraction or P4.
    # Averaging B or A across geometries first would define a different readout.
    for name in QUANTITIES:
        contrast = quantities[name][0] - quantities[name][1]
        for part, value in (("re", contrast.real), ("im", contrast.imag)):
            values[f"pair.{name}_{part}"] = float(value)
            values[f"P4.{name}_{part}"] = float(value / delta)
    values.update(common_clock_wedges(quantities["conditional_cov_s"], quantities["mu_p"]))
    if not all(map(math.isfinite, values.values())):
        raise ValueError("conditional primitive-line response contains nonfinite values")
    diagnostics = {
        "N": int(n), "p0": p0, "samples_per_geometry": samples,
        "input_columns": list(INPUT_COLUMNS), "delta_cos4_exact": str(delta_exact),
        "conditioning": "I1 within each geometry, followed by first-minus-second; imaginary parts retained",
        "source_units": "bulk s=CB+CW; common E is dimensionless and has no factor N",
        "matching_slope": matching_slope,
        "baseline_q": m.tolist(), "baseline_E": e.tolist(),
        "directions": direction_diagnostics,
        "conditional_covariance_vector_order": list(COMPLEX_COMPONENTS),
        "conditional_covariance_bulk_vector": _real_vector(quantities["conditional_cov_s"]).tolist(),
        "mu_p_vector": _real_vector(quantities["mu_p"]).tolist(),
        "moving_source_vector": _real_vector(quantities["nu_s"]).tolist(),
        "common_E_prediction_vector": _real_vector(quantities["pred_E"]).tolist(),
        "common_E_residual_vector": _real_vector(quantities["resid_E"]).tolist(),
        "clock_null": "a_g+b*K+c*q+d*E, with common real b,c,d: conditional_cov_s=b*p0*(1-p0)*mu_p in every geometry",
        "clock_null_thermal_factor": p0 * (1 - p0),
        "wedge_definition": "covariance[i]*mu_p[j]-covariance[j]*mu_p[i], all six i<j, no slope division or selected component",
        "wedge_dependence": "the six minors are correlated; with nonzero mu_p the collinearity null has at most three independent first-order constraints, not six automatic chi-square degrees of freedom",
        "small_slope_boundary": "if mu_p is zero the wedge alone is uninformative although the null predicts zero conditional covariance; retain both four-real vectors",
        "wedge_root_invariance": "adding the same real multiple of mu_p, including root motion or the E-loading prediction, leaves all six wedges unchanged",
        "loading_definition": "alpha_E=rank1rootdot_s/rank1rootdot_E; shared across directions, recomputed by caller's aligned omissions",
        "interpretation": "These are one conditional-observer experiment and correlated component readouts; null compatibility is not a pointwise microscopic-source identity or energy-field identification",
    }
    return values, diagnostics
