"""Two-phase norm-4 source response, with a disjoint q/E production anchor.

Numerical core only: no input/output, sampling, covariance orchestration or
analysis entry point.  Input source profiles use read_raw's density units:
(q, E, S=s/N, qS, ES).  Returned bulk responses differentiate exp(t*s).
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

from analyze_norm4_source_thermal import binomial_moments

P_REF = 0.59274605079
BATCHES = 100
SAMPLES_PER_BATCH = 1000
FOLDS = 5
COEFFICIENT_NAMES = (
    "a_first", "a_second", "b_common_K_over_N",
    "c_first_q", "c_second_q", "d_first_E", "d_second_E",
)
TERM_NAMES = ("direct", "rootmotion", "slope_source", "slope_root")


def baseline(qe_sums, samples, n, delta, bracket):
    """Build the unmarked complement anchor from its two-direction q/E sums.

    qe_sums: (2,N+1,2), with integer-derived sums rather than means.
    samples: number of complement permutations in EACH direction.
    Each call refinds the local pooled root; baseline delete-groups therefore
    propagate changes of p0, D, B, H and T rather than holding them fixed.
    """
    profiles = np.asarray(qe_sums, dtype=float)
    samples, delta = float(samples), float(delta)
    if profiles.shape != (2, n + 1, 2) or not np.isfinite(profiles).all():
        raise ValueError("baseline requires finite (2,N+1,2) q/E sum profiles")
    if samples <= 0 or not math.isfinite(samples) or delta == 0 or not math.isfinite(delta):
        raise ValueError("positive samples and finite nonzero exact-direction delta required")
    low, high = map(float, bracket)
    if not 0 < low < high < 1:
        raise ValueError("baseline root bracket must be inside (0,1)")

    def packets(p):
        return [binomial_moments(profiles[g], samples, p, n) for g in range(2)]

    def pooled_q(p):
        return sum(packet[0][0] for packet in packets(p)) / 2

    endpoint_values = [float(pooled_q(low)), float(pooled_q(high))]
    if endpoint_values[0] * endpoint_values[1] > 0:
        raise ValueError(f"N{n}: complement pooled root is not bracketed")
    root = float(brentq(pooled_q, low, high, xtol=5e-14, rtol=5e-14))
    packet = packets(root)
    mean = np.asarray([row[0] for row in packet])
    first = np.asarray([row[1] for row in packet])
    second = np.asarray([row[2] for row in packet])
    d = float(first[:, 0].mean())
    b = float((first[0, 1] - first[1, 1]) / delta)
    h = float((second[0, 1] - second[1, 1]) / delta)
    t = float(second[:, 0].mean())
    if not all(map(math.isfinite, (d, b, h, t))) or d <= 0:
        raise ValueError(f"N{n}: complement has nonpositive or nonfinite local slope")
    return {
        "root": root, "p0": root,
        "m": mean[:, 0].tolist(), "e": mean[:, 1].tolist(),
        "m_p": first[:, 0].tolist(), "e_p": first[:, 1].tolist(),
        "m_pp": second[:, 0].tolist(), "e_pp": second[:, 1].tolist(),
        "D": d, "B": b, "H": h, "T": t,
        "U": float(n ** (13 / 8) * b / (2 * d)),
        "samples": samples, "n": int(n), "delta": delta,
        "bracket": [low, high], "bracket_pooled_q": endpoint_values,
        "pooled_q": float(mean[:, 0].mean()), "integration": packet[0][3],
    }


def _response(anchor, n, delta, jq, je, jq_p, je_p):
    """The fixed-anchor linear functional L on one density-normalized source."""
    jq, je, jq_p, je_p = [np.asarray(value, dtype=float) for value in (jq, je, jq_p, je_p)]
    d, b, h, t = (float(anchor[key]) for key in ("D", "B", "H", "T"))
    jq_bar, jqp_bar = float(jq.mean()), float(jq_p.mean())
    je_p4 = float((je_p[0] - je_p[1]) / delta)
    prefactor = n ** (13 / 8) / 2
    pieces = [prefactor * je_p4 / d,
              -prefactor * h * jq_bar / d**2,
              -prefactor * b * jqp_bar / d**2,
              prefactor * b * t * jq_bar / d**3]
    value = math.fsum(pieces)
    if not math.isfinite(value):
        raise ValueError(f"N{n}: nonfinite two-phase source response")
    return {
        "density": float(value), "bulk": float(n * value),
        "terms_density": dict(zip(TERM_NAMES, map(float, pieces))),
        "terms_bulk": dict(zip(TERM_NAMES, (float(n * piece) for piece in pieces))),
        "rootdot_density": float(-jq_bar / d), "rootdot_bulk": float(-n * jq_bar / d),
        "Jq_density": jq.tolist(), "JE_density": je.tolist(),
        "Jq_p_density": jq_p.tolist(), "JE_p_density": je_p.tolist(),
    }


def _source_response(source_sums, samples, anchor, n, delta):
    """Center held-out S,qS,ES means using the independent anchor's q/E jets."""
    packet = [binomial_moments(source_sums[g], samples, anchor["root"], n) for g in range(2)]
    mean = np.asarray([row[0] for row in packet])
    first = np.asarray([row[1] for row in packet])
    m, e, mp, ep = [np.asarray(anchor[key], dtype=float) for key in ("m", "e", "m_p", "e_p")]
    s, qs, es = mean.T
    sp, qsp, esp = first.T
    jq = qs - m * s
    je = es - e * s
    jq_p = qsp - mp * s - m * sp
    je_p = esp - ep * s - e * sp
    result = _response(anchor, n, delta, jq, je, jq_p, je_p)
    result["source_means_density"] = mean.tolist()
    result["source_first_derivatives_density"] = first.tolist()
    return result


def _topology_response(coefficients, anchor, n, delta):
    """Use qE=q and E^2=E to evaluate c_g*q+d_g*E on the large anchor."""
    c = coefficients[[3, 4]]
    d = coefficients[[5, 6]]
    m, e, mp, ep = [np.asarray(anchor[key], dtype=float) for key in ("m", "e", "m_p", "e_p")]
    odd_even = m * (1 - e)
    odd_even_p = mp * (1 - e) - m * ep
    jq = c * (e - m * m) + d * odd_even
    je = c * odd_even + d * e * (1 - e)
    jq_p = c * (ep - 2 * m * mp) + d * odd_even_p
    je_p = c * odd_even_p + d * ep * (1 - 2 * e)
    return _response(anchor, n, delta, jq, je, jq_p, je_p)


def _fit(training_sums, samples, n):
    """Fit all seven frozen columns at p_ref, without requiring S^2.

    The common thermal column is K/N, keeping coefficient magnitudes in the
    same density-source units as the other columns.  Normal equations use the
    exact three-state algebra and are solved after diagonal column scaling;
    no ridge penalty, column dropping or outcome-selected basis is used.
    """
    k_over_n = np.arange(n + 1, dtype=float) / n
    gram = np.zeros((7, 7), dtype=float)
    rhs = np.zeros(7, dtype=float)
    for g in range(2):
        q, e, s, qs, es = training_sums[g].T
        packed = np.column_stack((
            np.full(n + 1, samples, dtype=float),
            samples * k_over_n, samples * k_over_n**2,
            q, e, k_over_n * q, k_over_n * e,
            s, k_over_n * s, qs, es,
        ))
        values = binomial_moments(packed, samples, P_REF, n)[0]
        one, k, kk, m, ev, kq, ke, sv, ks, qsv, esv = values
        # Local column order a_g, shared K/N, c_g*q, d_g*E.
        indices = [g, 2, 3 + g, 5 + g]
        block = np.asarray([
            [one, k, m, ev],
            [k, kk, kq, ke],
            [m, kq, ev, m],
            [ev, ke, m, ev],
        ])
        gram[np.ix_(indices, indices)] += block / 2
        rhs[indices] += np.asarray([sv, ks, qsv, esv]) / 2
    scaling = np.sqrt(np.diag(gram))
    if np.any(scaling <= 0) or not np.isfinite(gram).all() or not np.isfinite(rhs).all():
        raise ValueError(f"N{n}: seven-column source Gram is not finite with positive column norms")
    scaled_gram = gram / np.outer(scaling, scaling)
    scaled_rhs = rhs / scaling
    coefficients = np.linalg.solve(scaled_gram, scaled_rhs) / scaling
    eigenvalues = np.linalg.eigvalsh(scaled_gram)
    if eigenvalues[0] <= 0:
        raise np.linalg.LinAlgError(f"N{n}: frozen seven-column source Gram is not positive definite")
    norm_rhs = float(np.linalg.norm(rhs))
    return coefficients, {
        "p_ref": P_REF,
        "coefficient_names": list(COEFFICIENT_NAMES),
        "coefficients_density": coefficients.tolist(),
        "coefficients_by_name": dict(zip(COEFFICIENT_NAMES, map(float, coefficients))),
        "scaled_gram_condition": float(eigenvalues[-1] / eigenvalues[0]),
        "scaled_gram_eigenvalues": eigenvalues.tolist(),
        "scaled_gram_rank": int(np.linalg.matrix_rank(scaled_gram)),
        "column_scales": scaling.tolist(),
        "normal_equation_relative_residual": float(np.linalg.norm(gram @ coefficients - rhs) / max(norm_rhs, np.finfo(float).tiny)),
        "source_variance_or_R_squared": "not estimated; S^2 is neither present nor needed",
    }


def _residual_sums(heldout_sums, samples, coefficients, n):
    k_over_n = np.arange(n + 1, dtype=float) / n
    residual = np.empty((2, n + 1, 3), dtype=float)
    for g in range(2):
        a, b, c, d = coefficients[g], coefficients[2], coefficients[3 + g], coefficients[5 + g]
        q, e, s, qs, es = heldout_sums[g].T
        clock = b * k_over_n
        residual[g, :, 0] = s - (a + clock) * samples - c * q - d * e
        residual[g, :, 1] = qs - (a + clock) * q - c * e - d * q
        residual[g, :, 2] = es - (a + clock) * e - c * q - d * e
    return residual


def estimate(marked_profiles, anchor, n, delta, omitted_batch=None):
    """Anchor-only and cross-fitted two-phase estimates of the same bulk source.

    marked_profiles is (100,2,N+1,5), each batch holding 1000 old permutations.
    The anchor belongs to the disjoint unmarked complement.  Deleting a marked
    batch refits each fold from its remaining training batches; fold identities
    remain original batch_id % 5, and held-out weights use actual sample counts.
    """
    marked = np.asarray(marked_profiles, dtype=float)
    delta = float(delta)
    if marked.shape != (BATCHES, 2, n + 1, 5) or not np.isfinite(marked).all():
        raise ValueError("marked profiles must be finite (100,2,N+1,5) read_raw density sums")
    if anchor["n"] != n or not math.isclose(float(anchor["delta"]), delta, rel_tol=1e-14, abs_tol=0):
        raise ValueError("anchor size or frozen direction projection differs")
    if omitted_batch is not None and (
            not isinstance(omitted_batch, (int, np.integer)) or not 0 <= omitted_batch < BATCHES):
        raise ValueError("omitted_batch must be an original batch index 0..99")
    retained = np.arange(BATCHES)
    if omitted_batch is not None:
        retained = retained[retained != omitted_batch]
    total_samples = len(retained) * SAMPLES_PER_BATCH
    total = marked[retained].sum(axis=0)
    anchor_only = _source_response(total[:, :, 2:], total_samples, anchor, n, delta)
    folds = []
    topology_parts, residual_parts = [], []
    for fold in range(FOLDS):
        heldout_ids = retained[retained % FOLDS == fold]
        training_ids = retained[retained % FOLDS != fold]
        heldout_samples = len(heldout_ids) * SAMPLES_PER_BATCH
        training_samples = len(training_ids) * SAMPLES_PER_BATCH
        heldout = marked[heldout_ids].sum(axis=0)
        training = total - heldout
        coefficients, diagnostic = _fit(training, training_samples, n)
        topology = _topology_response(coefficients, anchor, n, delta)
        residual_profiles = _residual_sums(heldout, heldout_samples, coefficients, n)
        residual = _source_response(residual_profiles, heldout_samples, anchor, n, delta)
        weight = heldout_samples / total_samples
        topology_parts.append(weight * topology["density"])
        residual_parts.append(weight * residual["density"])
        folds.append({
            "fold": fold, "training_batch_ids": training_ids.tolist(),
            "heldout_batch_ids": heldout_ids.tolist(), "training_samples": training_samples,
            "heldout_samples": heldout_samples, "weight": float(weight),
            "fit": diagnostic, "topology": topology, "residual": residual,
            "clock_bulk": 0.0, "constant_bulk": 0.0,
            "two_phase_bulk": float(n * math.fsum((topology["density"], residual["density"]))),
        })
    topology_density = math.fsum(topology_parts)
    residual_density = math.fsum(residual_parts)
    two_phase_density = math.fsum((topology_density, residual_density))
    return {
        "anchor_only_bulk": anchor_only["bulk"],
        "two_phase_bulk": float(n * two_phase_density),
        "topology_bulk": float(n * topology_density),
        "residual_bulk": float(n * residual_density),
        "anchor_only_density": anchor_only["density"],
        "two_phase_density": float(two_phase_density),
        "topology_density": float(topology_density),
        "residual_density": float(residual_density),
        "clock_bulk": 0.0, "constant_bulk": 0.0,
        "anchor_only": anchor_only,
        "crossfit": {
            "fold_rule": "original batch_id modulo 5", "fit_p_ref": P_REF,
            "retained_marked_samples": total_samples,
            "omitted_batch": None if omitted_batch is None else int(omitted_batch),
            "coefficient_names": list(COEFFICIENT_NAMES), "folds": folds,
            "max_scaled_gram_condition": max(fold["fit"]["scaled_gram_condition"] for fold in folds),
            "min_scaled_gram_rank": min(fold["fit"]["scaled_gram_rank"] for fold in folds),
            "thermal_derivative_rule": "training coefficients fixed in p; no coefficient p-derivatives",
            "clock_rule": "the shared b*(K/N) term has exactly zero L_anchor response and is omitted analytically",
        },
        "source_coordinate": "input density S=(CB+CW)/N; all bulk values multiply the density derivative by N",
        "interpretation": "topology and residual allocate the same source response; cross-fitted components and alternative estimators are not independent evidence",
    }
