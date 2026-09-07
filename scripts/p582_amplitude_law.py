#!/usr/bin/env python3
"""The amplitude law of #582's frozen shape direction -- issue #584, after Gate 3.

Gate 3 (``scripts/score_type582_residual.py``) asked whether a pre-existing
discrete label indexes the small structured remainder that survives #582's
dominant shape direction.  No label beat the exact permutation null, and the
verdict was to keep the dominant direction as the robust finite object.

This module asks the question Gate 3 left on the table: **the dominant direction
was kept -- what law does its amplitude obey?**

That is a different object from the one #584 rules out.  #584 says

    do not fit a free exponent to five residual *directions*.

Nothing here fits a direction.  The direction ``g`` is frozen first, by the same
power iteration Gate 3 uses, and one exponent is then fitted to the five
*amplitudes* of that already-frozen direction -- five numbers spanning a factor
of 2.8, each measured to better than 0.3%.  The fit is then asked to predict a
transition it never saw, and finally confronted with a functional of the data it
was not fitted to at all.

**The model.**  If the shape part of the threshold law carries a single
correction-to-scaling term,

    Q(s, u) = A(u) + lambda * exp(-omega * s) * g(u) + ...,     s = log N,

then the finite difference #582 actually measures is not the derivative but its
exact image,

    v = [Q(s_t) - Q(s_b)] / h
      = -lambda * omega * exp(-omega * sbar) * sinh(omega*h/2) / (omega*h/2) * g,

with ``sbar`` the log-midpoint and ``h`` the log-step.  The ``sinh(x)/x`` factor
is 1.9% at ``h = log 2`` and 3.3% at ``h = log 2.5`` -- step-size dependent, and
the same size as the misfit being measured -- so it is carried exactly rather
than dropped.  Dropping it would manufacture a difference between the m=2 and
m=2.5 transitions, which is precisely the ``multiplier`` label Gate 3 screened.

**The independent check that decides the result.**  A second divided difference
across a three-size lineage is a different functional of the same productions.
Applied to the model above it is an exact linear combination, ``sum_k c_k
exp(-omega s_k)``, with weights that annihilate constants.  The parameters
fitted on the five *first* differences are frozen and asked to predict it.  They
do not: the measured second difference is about 55% larger, in **both**
three-size lineages, and the two ratios agree with each other to about 1%.  A
single exponent therefore describes the tangent well and the curvature badly,
and it fails reproducibly rather than noisily.  GOVERNANCE section 2 minimum A
asks for one check by independent means; this is that check, and it fired.

Nothing is sampled here.  Every histogram is a committed production artifact.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import mpmath as mp

try:  # pragma: no cover - import shape depends on how the script is invoked
    from scripts.projective_inference import (spectral_pseudo_inverse,
                                              subspace_residual)
    from scripts import score_wasserstein_shape_flow as flow
    from scripts import threshold_quantile_lineage as lineage
except ModuleNotFoundError:  # pragma: no cover
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from projective_inference import spectral_pseudo_inverse, subspace_residual
    import score_wasserstein_shape_flow as flow
    import threshold_quantile_lineage as lineage

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "p582-amplitude-law" / "latest.json"
COMMITTED_FLOW = ROOT / "results" / "wasserstein-shape-flow" / "latest.json"
SCHEMA = "matching-one.p582-amplitude-law.v1"
ISSUE = 584

LEVELS = 9

#: The orientation weighting #582 declared primary, and the naive one it warns
#: about.  Both are run: the difference between the exponents they give is the
#: orientation systematic, and it is a real part of the error budget.
PRIMARY_WEIGHTING = "spin0"
CONTROL_WEIGHTING = "equal"

#: Search window for the exponent.  Wide enough to contain any plausible answer
#: and narrow enough that the profile is unimodal; the fitted value is reported
#: with the window so a boundary hit would be visible rather than silent.
OMEGA_WINDOW = (0.25, 2.50)
GOLDEN_STEPS = 200

#: Sizes proposed for the crossed/lever-arm acquisition.  Declared here, before
#: any scoring, so the candidate list cannot be chosen after seeing which one
#: the analysis would prefer.  ``650`` is the size #583/#589 already design for.
CANDIDATE_SIZES = (338, 442, 650, 725, 845, 1105, 1690)

#: A production needs two orientations, so a size with a single primitive
#: Gaussian representative cannot enter this pipeline at all.
MINIMUM_PRIMITIVE_REPRESENTATIVES = 2


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


# --------------------------------------------------------------------------
# transitions, covariances, and the control that we are looking at #582's object
# --------------------------------------------------------------------------

def transition_rows(loaded: Mapping[int, Mapping[str, Any]]) -> list[dict[str, Any]]:
    """#582's five transitions, with the log geometry this module needs.

    ``sbar`` and ``h`` are the log-midpoint and the log-step.  Neither appears in
    #582's artifact, and both are what turn five amplitudes into a law.
    """
    rows = flow.transitions()
    for row in rows:
        base, target = row["base"], row["target"]
        scale = math.log(row["multiplier"])
        row["sbar"] = (math.log(base) + math.log(target)) / 2.0
        row["h"] = math.log(target) - math.log(base)
        row["covariance"] = [
            [(a + b) / (scale * scale) for a, b in zip(left, right)]
            for left, right in zip(loaded[base]["covariance"],
                                   loaded[target]["covariance"])
        ]
        row["affine_only"] = flow.shape_flow(loaded[base], loaded[target],
                                             row["multiplier"])
    return rows


def reproduction_control(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Re-derived affine statistics against the committed #582 artifact.

    The wrong number this stops us believing is an amplitude law fitted to a
    re-read of the histograms that is not the object #582 published -- a changed
    quantile grid, a changed jackknife, a changed orientation weighting.  The
    agreement here is exact, not approximate.
    """
    committed = json.loads(COMMITTED_FLOW.read_text())
    published = {row["label"]: row["affine_only"]
                 for row in committed["scored"][PRIMARY_WEIGHTING]["transitions"]}
    out = []
    for row in rows:
        mine = row["affine_only"]
        theirs = published[row["label"]]
        out.append({
            "label": row["label"],
            "statistic": mine["statistic"],
            "published_statistic": theirs["statistic"],
            "largest_residual_difference":
                max(abs(a - b) for a, b in zip(mine["residual"], theirs["residual"])),
            "statistic_matches_exactly": mine["statistic"] == theirs["statistic"],
        })
    return {
        "per_transition": out,
        "every_statistic_matches_exactly": all(item["statistic_matches_exactly"]
                                               for item in out),
    }


def consensus_direction(rows: Sequence[Mapping[str, Any]]) -> list[float]:
    """The frozen dominant shape direction: Gate 3's ``g_N``, rebuilt identically.

    Power iteration on the sum of outer products of the five unit affine
    residuals.  The direction is fixed here and never touched again; every
    amplitude below is read along it.
    """
    residuals = [flow._normalise(row["affine_only"]["residual"]) for row in rows]
    return flow._leading_direction(residuals)


# --------------------------------------------------------------------------
# amplitudes along the frozen direction, with their exact standard errors
# --------------------------------------------------------------------------

def amplitude(row: Mapping[str, Any], base_quantiles: Sequence[float],
              observation: Sequence[float], covariance: Sequence[Sequence[float]],
              direction: Sequence[float]) -> dict[str, Any]:
    """Amplitude of ``direction`` in ``observation``, beside the affine nuisance.

    The standard error is the ``(2, 2)`` entry of the inverse weighted Gram of
    the *full* basis, not ``1 / (g^T S^+ g)``.  Those two differ by the affine
    nuisance's leverage on ``g``, and using the second would understate the
    amplitude's uncertainty -- which is precisely the number the exponent fit is
    weighted by.
    """
    basis = [[1.0] * LEVELS, list(base_quantiles), list(direction)]
    fit = subspace_residual(observation, covariance, basis)
    pinv, _, _, _ = spectral_pseudo_inverse(covariance)
    columns = [mp.matrix(vector) for vector in basis]
    gram = mp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            gram[i, j] = (columns[i].T * (pinv * columns[j]))[0]
    inverse = gram ** -1
    return {
        "amplitude": fit["amplitudes"][2],
        "standard_error": float(mp.sqrt(inverse[2, 2])),
        "residual_statistic": fit["statistic"],
        "residual_degrees_of_freedom": fit["degrees_of_freedom"],
    }


def amplitude_points(rows: Sequence[Mapping[str, Any]],
                     loaded: Mapping[int, Mapping[str, Any]],
                     direction: Sequence[float]) -> list[dict[str, Any]]:
    """One ``(sbar, h, amplitude, standard error)`` point per transition."""
    points = []
    for row in rows:
        measured = amplitude(row, loaded[row["base"]]["quantiles"],
                             row["affine_only"]["displacement"],
                             row["covariance"], direction)
        points.append({
            "label": row["label"],
            "lineage": row["lineage"],
            "base": row["base"],
            "target": row["target"],
            "sbar": row["sbar"],
            "h": row["h"],
            "multiplier": row["multiplier"],
            "amplitude": measured["amplitude"],
            "standard_error": measured["standard_error"],
            "significance": measured["amplitude"] / measured["standard_error"],
            "statistic_after_removing_the_direction": measured["residual_statistic"],
            "degrees_of_freedom_after_removing_the_direction":
                measured["residual_degrees_of_freedom"],
        })
    return points


# --------------------------------------------------------------------------
# the one-exponent model and its exact finite-difference image
# --------------------------------------------------------------------------

def first_difference_image(scale_amplitude: float, omega: float,
                           sbar: float, step: float) -> float:
    """Exact finite-difference image of ``lambda * exp(-omega * s)``.

    Not ``-lambda * omega * exp(-omega * sbar)``: that is the derivative, and it
    is wrong here by ``sinh(x)/x`` with ``x = omega*h/2``, which at the fitted
    exponent is 1.9% at ``h = log 2`` and 3.3% at ``h = log 2.5``.  Those are
    step-size dependent and the same size as the
    misfit being measured, so dropping the factor would manufacture a spurious
    step-size dependence and, with it, a spurious label effect on ``multiplier``.
    """
    _require(step > 0.0, "a transition needs a positive log-step")
    _require(omega != 0.0, "the exponential image is singular at omega = 0")
    half = omega * step / 2.0
    return -scale_amplitude * omega * math.exp(-omega * sbar) * math.sinh(half) / half


def second_difference_weights(nodes: Sequence[float]) -> list[float]:
    """Weights of ``2 * f[s0, s1, s2]`` on an unequal grid.

    They sum to zero, which is what makes the second difference blind to the
    additive limit shape ``A(u)``; a weight set that did not would let the whole
    threshold law leak into a curvature amplitude.
    """
    _require(len(nodes) == 3, "a second divided difference needs exactly three nodes")
    first, middle, last = nodes
    left, right = middle - first, last - middle
    _require(left > 0.0 and right > 0.0, "nodes must increase")
    span = left + right
    return [2.0 / (left * span), -2.0 / (left * right), 2.0 / (right * span)]


def second_difference_image(scale_amplitude: float, omega: float,
                            nodes: Sequence[float]) -> float:
    """Exact second-divided-difference image of ``lambda * exp(-omega * s)``."""
    weights = second_difference_weights(nodes)
    return scale_amplitude * sum(weight * math.exp(-omega * node)
                                 for weight, node in zip(weights, nodes))


def _profiled_scale(omega: float, points: Sequence[Mapping[str, Any]]) -> float:
    """``lambda`` is linear given ``omega``, so it is profiled out exactly."""
    numerator = sum(point["amplitude"] * first_difference_image(1.0, omega,
                                                                point["sbar"],
                                                                point["h"])
                    / point["standard_error"] ** 2 for point in points)
    denominator = sum(first_difference_image(1.0, omega, point["sbar"], point["h"]) ** 2
                      / point["standard_error"] ** 2 for point in points)
    _require(denominator > 0.0, "the model direction has no weight")
    return numerator / denominator


def _chi_square(scale_amplitude: float, omega: float,
                points: Sequence[Mapping[str, Any]]) -> float:
    return sum(((point["amplitude"]
                 - first_difference_image(scale_amplitude, omega,
                                          point["sbar"], point["h"]))
                / point["standard_error"]) ** 2 for point in points)


def fit_exponent(points: Sequence[Mapping[str, Any]],
                 window: tuple[float, float] = OMEGA_WINDOW) -> dict[str, Any]:
    """One exponent, ``lambda`` profiled out, by golden-section on the window."""
    low, high = window
    for _ in range(GOLDEN_STEPS):
        left = low + (high - low) / 3.0
        right = high - (high - low) / 3.0
        if _chi_square(_profiled_scale(left, points), left, points) \
                < _chi_square(_profiled_scale(right, points), right, points):
            high = right
        else:
            low = left
    omega = (low + high) / 2.0
    scale_amplitude = _profiled_scale(omega, points)
    margin = min(omega - window[0], window[1] - omega)
    return {
        "omega": omega,
        "scale_amplitude": scale_amplitude,
        "statistic": _chi_square(scale_amplitude, omega, points),
        "degrees_of_freedom": len(points) - 2,
        "search_window": list(window),
        "distance_to_the_nearest_window_edge": margin,
        "window_binds": margin < 1e-3,
    }


def fit_report(points: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """The fit, its per-point pulls, and the unit-exponent comparison."""
    fitted = fit_exponent(points)
    per_point = []
    for point in points:
        model = first_difference_image(fitted["scale_amplitude"], fitted["omega"],
                                       point["sbar"], point["h"])
        per_point.append({
            "label": point["label"],
            "measured": point["amplitude"],
            "model": model,
            "pull": (point["amplitude"] - model) / point["standard_error"],
            "relative_error": point["amplitude"] / model - 1.0,
        })
    unit_scale = _profiled_scale(1.0, points)
    unit = _chi_square(unit_scale, 1.0, points)
    return {
        **fitted,
        "per_transition": per_point,
        "largest_relative_error": max(abs(item["relative_error"])
                                      for item in per_point),
        "unit_exponent_statistic": unit,
        "unit_exponent_excess_statistic": unit - fitted["statistic"],
    }


def leave_one_transition_out(points: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Fit on four transitions, predict the fifth's amplitude.

    The exponent is never allowed to see the transition it predicts.  The wrong
    number this stops us believing is an in-sample exponent quoted as if it
    forecast a new size: the spread of the four-point exponents below is the
    honest statistical part of the error budget on ``omega``.
    """
    folds = []
    for index, held in enumerate(points):
        training = [point for offset, point in enumerate(points) if offset != index]
        fitted = fit_exponent(training)
        model = first_difference_image(fitted["scale_amplitude"], fitted["omega"],
                                       held["sbar"], held["h"])
        folds.append({
            "held_out": held["label"],
            "omega_trained_on_the_other_four": fitted["omega"],
            "predicted": model,
            "measured": held["amplitude"],
            "relative_error": held["amplitude"] / model - 1.0,
            "pull": (held["amplitude"] - model) / held["standard_error"],
        })
    exponents = [fold["omega_trained_on_the_other_four"] for fold in folds]
    return {
        "folds": folds,
        "omega_range": [min(exponents), max(exponents)],
        "omega_half_spread": (max(exponents) - min(exponents)) / 2.0,
        "largest_relative_prediction_error":
            max(abs(fold["relative_error"]) for fold in folds),
    }


# --------------------------------------------------------------------------
# the independent check: a functional the fit never saw
# --------------------------------------------------------------------------

def second_difference_check(loaded: Mapping[int, Mapping[str, Any]],
                            direction: Sequence[float],
                            scale_amplitude: float, omega: float) -> dict[str, Any]:
    """Predict each three-size lineage's curvature from the first-difference fit.

    ``scale_amplitude`` and ``omega`` are frozen: they were fitted on the five
    first differences and are not refitted here, so this is a prediction and not
    a description.  A one-exponent law would return a ratio of one.  The wrong
    number this stops us believing is a single correction-to-scaling term read
    off the tangent alone -- which fits the tangent to 2% and misses the
    curvature by half.
    """
    out = []
    for name, sizes in flow.LINEAGES.items():
        if len(sizes) < 3:
            out.append({"lineage": name, "sizes": list(sizes),
                        "usable": False,
                        "reason": "a second difference needs three sizes"})
            continue
        nodes = [math.log(size) for size in sizes]
        weights = second_difference_weights(nodes)
        observation = [sum(weight * loaded[size]["quantiles"][level]
                           for weight, size in zip(weights, sizes))
                       for level in range(LEVELS)]
        covariance = [[sum(weight * weight * loaded[size]["covariance"][i][j]
                           for weight, size in zip(weights, sizes))
                       for j in range(LEVELS)] for i in range(LEVELS)]
        middle = sizes[1]
        measured = amplitude({}, loaded[middle]["quantiles"], observation,
                             covariance, direction)
        predicted = second_difference_image(scale_amplitude, omega, nodes)
        out.append({
            "lineage": name,
            "sizes": list(sizes),
            "usable": True,
            "weights": weights,
            "weights_annihilate_constants": abs(sum(weights)) < 1e-9,
            "measured": measured["amplitude"],
            "standard_error": measured["standard_error"],
            "predicted_from_the_first_differences": predicted,
            "ratio": measured["amplitude"] / predicted,
            "pull": (measured["amplitude"] - predicted) / measured["standard_error"],
        })
    usable = [item for item in out if item["usable"]]
    ratios = [item["ratio"] for item in usable]
    spread = (max(ratios) - min(ratios)) / (sum(ratios) / len(ratios)) if ratios else None
    return {
        "per_lineage": out,
        "lineages_checked": len(usable),
        "ratio_range": [min(ratios), max(ratios)] if ratios else None,
        "relative_spread_between_lineages": spread,
        "the_discrepancy_reproduces_across_lineages":
            bool(ratios) and spread is not None and spread < 0.05
            and min(abs(value - 1.0) for value in ratios) > 0.10,
    }


# --------------------------------------------------------------------------
# the orientation-weighting systematic
# --------------------------------------------------------------------------

def weighting_systematic(both: Mapping[str, Mapping[int, Mapping[str, Any]]]
                         ) -> dict[str, Any]:
    """Refit everything under the naive equal-orientation weighting.

    #582 warns that the equal-weight residue alternates in sign along a lineage.
    The exponent moves when the weighting does, and that shift is a systematic,
    not a discovery.  The wrong number this stops us believing is a
    statistics-only error bar on ``omega`` -- which would exclude ``omega = 1``
    on a systematic that is larger than it.

    This control also carries the second-difference check, for a reason that is
    not about ``omega`` at all.  Under the primary weighting the spin-0
    combination **extrapolates** at 325 and 425 (weights ``1.278 / -0.278`` and
    ``-0.026 / 1.026``) and interpolates everywhere else -- so rung 3 of *both*
    three-size lineages is an extrapolation, and the two lineages share that
    structure exactly, along with a 5x larger per-batch sample count at the same
    two sizes.  The close agreement of their two ratios is therefore **not**
    independent evidence on its own.  The equal weighting is ``0.5 / 0.5`` at
    every size and never extrapolates, so running the check under both is what
    separates a shared systematic from a shared law.
    """
    report = {}
    directions = {}
    for name in (PRIMARY_WEIGHTING, CONTROL_WEIGHTING):
        loaded = both[name]
        rows = transition_rows(loaded)
        direction = consensus_direction(rows)
        directions[name] = direction
        points = amplitude_points(rows, loaded, direction)
        fitted = fit_exponent(points)
        check = second_difference_check(loaded, direction,
                                        fitted["scale_amplitude"], fitted["omega"])
        report[name] = {
            "omega": fitted["omega"],
            "scale_amplitude": fitted["scale_amplitude"],
            "statistic": fitted["statistic"],
            "degrees_of_freedom": fitted["degrees_of_freedom"],
            "amplitudes": [point["amplitude"] for point in points],
            "second_difference_ratios":
                {item["lineage"]: item["ratio"]
                 for item in check["per_lineage"] if item["usable"]},
        }
    cosine = sum(a * b for a, b in zip(directions[PRIMARY_WEIGHTING],
                                       directions[CONTROL_WEIGHTING]))
    cosine = max(-1.0, min(1.0, cosine))
    shift = abs(report[PRIMARY_WEIGHTING]["omega"] - report[CONTROL_WEIGHTING]["omega"])
    cells = {f"{name}/{lineage}": ratio
             for name, block in report.items()
             for lineage, ratio in block["second_difference_ratios"].items()}
    values = list(cells.values())
    return {
        "per_weighting": report,
        "angle_between_the_two_frozen_directions_degrees":
            math.degrees(math.acos(abs(cosine))),
        "omega_shift": shift,
        "second_difference_ratio_cells": cells,
        "second_difference_ratio_range": [min(values), max(values)],
        "second_difference_ratio_spread_over_all_cells":
            (max(values) - min(values)) / (sum(values) / len(values)),
        "the_discrepancy_survives_the_weighting_change":
            min(values) > 1.2 and max(values) < 2.0,
    }


def error_budget(fit: Mapping[str, Any], folds: Mapping[str, Any],
                 weighting: Mapping[str, Any]) -> dict[str, Any]:
    """``omega`` with both parts of its uncertainty, and what that excludes.

    The wrong number this stops us believing is ``omega != 1`` claimed from the
    fit's own curvature: the fit is rejected, so its internal error bar is not
    the uncertainty on ``omega``.  The honest bar is the leave-one-out spread
    added to the weighting shift, and it reaches one.
    """
    central = fit["omega"]
    statistical = folds["omega_half_spread"]
    systematic = weighting["omega_shift"]
    total = math.sqrt(statistical ** 2 + systematic ** 2)
    return {
        "omega": central,
        "leave_one_out_half_spread": statistical,
        "orientation_weighting_shift": systematic,
        "combined_uncertainty": total,
        "distance_from_unity_in_combined_uncertainties":
            abs(central - 1.0) / total if total > 0.0 else None,
        "unit_exponent_is_excluded": abs(central - 1.0) > 2.0 * total,
        "the_fit_is_itself_rejected": fit["statistic"] > 20.0,
    }


# --------------------------------------------------------------------------
# what to measure next
# --------------------------------------------------------------------------

def primitive_representatives(size: int) -> list[tuple[int, int]]:
    """Primitive Gaussian representatives ``a > b > 0``, ``gcd(a, b) = 1``."""
    out = []
    for first in range(math.isqrt(size), 0, -1):
        remainder = size - first * first
        if remainder < 1:
            continue
        second = math.isqrt(remainder)
        if second * second == remainder and second <= first and math.gcd(first, second) == 1:
            out.append((first, second))
    return out


def five_adic_valuation(size: int) -> int:
    value, count = size, 0
    while value % 5 == 0:
        value //= 5
        count += 1
    return count


def candidate_sizes(sizes: Sequence[int] = CANDIDATE_SIZES) -> list[dict[str, Any]]:
    """Exact arithmetic on the declared candidates; nothing is fitted here.

    Gate 3 found that ``multiplier``, ``5-adic valuation of the target`` and
    ``target interpolation flag`` are the *same* partition of the existing five
    transitions, so no screen on this design can tell them apart.  A candidate
    breaks that degeneracy when its valuation and its interpolation flag
    disagree with the pattern the existing sizes set.  A candidate is unusable
    outright when it has fewer than two primitive representatives, because the
    pipeline needs two orientations.
    """
    out = []
    for size in sizes:
        representatives = primitive_representatives(size)
        cos4 = [lineage.cos_four_theta(a, b) for a, b in representatives]
        usable = len(representatives) >= MINIMUM_PRIMITIVE_REPRESENTATIVES
        straddles = usable and min(cos4) < 0.0 < max(cos4)
        valuation = five_adic_valuation(size)
        out.append({
            "size": size,
            "primitive_representatives": [list(item) for item in representatives],
            "cos_four_theta": cos4,
            "usable_in_this_pipeline": usable,
            "reason_if_unusable":
                None if usable else "fewer than two primitive representatives",
            "five_adic_valuation": valuation,
            "an_interpolating_spin_zero_combination_exists": straddles,
            "breaks_the_valuation_interpolation_degeneracy":
                bool(usable and ((valuation >= 2) == straddles)),
        })
    return out


def acquisition(candidates: Sequence[Mapping[str, Any]],
                second_difference: Mapping[str, Any]) -> dict[str, Any]:
    """The single next production, and the three things it settles at once.

    Chosen by the declared criteria, not by preference: usable in the pipeline,
    breaks the degeneracy Gate 3 exposed, and extends a two-size lineage to
    three so that a *third* independent curvature exists.
    """
    two_size = [name for name, sizes in flow.LINEAGES.items() if len(sizes) < 3]
    extends = {}
    for name in two_size:
        last = flow.LINEAGES[name][-1]
        for candidate in candidates:
            if candidate["usable_in_this_pipeline"] and candidate["size"] % last == 0:
                extends.setdefault(name, []).append(candidate["size"])
            elif candidate["usable_in_this_pipeline"] and \
                    abs(candidate["size"] / last - 2.5) < 1e-9:
                extends.setdefault(name, []).append(candidate["size"])
    preferred = [candidate for candidate in candidates
                 if candidate["usable_in_this_pipeline"]
                 and candidate["breaks_the_valuation_interpolation_degeneracy"]
                 and any(candidate["size"] in sizes for sizes in extends.values())]
    return {
        "two_size_lineages": two_size,
        "candidates_that_extend_a_two_size_lineage": extends,
        "preferred": [candidate["size"] for candidate in preferred],
        "what_it_settles": [
            "a third independent curvature, so the second-difference discrepancy "
            f"(ratio {second_difference['ratio_range']}) is tested outside the two "
            "lineages that produced it",
            "the multiplier / 5-adic-valuation / interpolation-flag degeneracy that "
            "Gate 3 showed makes every label screen on the current five transitions "
            "observationally identical to the production block",
            "a longer lever arm in N for the exponent itself",
        ],
    }


# --------------------------------------------------------------------------
# verdict
# --------------------------------------------------------------------------

def decide(reproduction: Mapping[str, Any], fit: Mapping[str, Any],
           folds: Mapping[str, Any], budget: Mapping[str, Any],
           second: Mapping[str, Any], weighting: Mapping[str, Any]) -> dict[str, Any]:
    """Verdict letters, each tied to a number above rather than to a reading."""
    if not reproduction["every_statistic_matches_exactly"]:
        return {"verdict": "REPRODUCTION_CONTROL_FAILED",
                "reason": "the re-derived affine statistics are not #582's"}
    findings = []
    if folds["largest_relative_prediction_error"] < 0.05:
        findings.append("AMPLITUDE_LAW_PREDICTS_A_HELD_OUT_TRANSITION")
    if not budget["unit_exponent_is_excluded"]:
        findings.append("EXPONENT_CONSISTENT_WITH_UNITY")
    if second["the_discrepancy_reproduces_across_lineages"] and \
            weighting["the_discrepancy_survives_the_weighting_change"]:
        findings.append("ONE_EXPONENT_FALSIFIED_BY_THE_SECOND_DIFFERENCE")
    verdict = ("ONE_EXPONENT_DESCRIBES_THE_TANGENT_AND_FAILS_THE_CURVATURE"
               if "ONE_EXPONENT_FALSIFIED_BY_THE_SECOND_DIFFERENCE" in findings
               and "AMPLITUDE_LAW_PREDICTS_A_HELD_OUT_TRANSITION" in findings
               else "AMPLITUDE_LAW_NOT_ESTABLISHED")
    return {"verdict": verdict, "findings": findings}


# --------------------------------------------------------------------------

def assemble() -> dict[str, Any]:
    both = flow.load_sizes(sorted(flow.SOURCES))
    loaded = both[PRIMARY_WEIGHTING]
    rows = transition_rows(loaded)
    reproduction = reproduction_control(rows)
    direction = consensus_direction(rows)
    points = amplitude_points(rows, loaded, direction)
    fit = fit_report(points)
    folds = leave_one_transition_out(points)
    second = second_difference_check(loaded, direction,
                                     fit["scale_amplitude"], fit["omega"])
    weighting = weighting_systematic(both)
    budget = error_budget(fit, folds, weighting)
    candidates = candidate_sizes()
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "primary_weighting": PRIMARY_WEIGHTING,
        "reproduction_control": reproduction,
        "frozen_direction": direction,
        "amplitudes": points,
        "one_exponent_fit": fit,
        "leave_one_transition_out": folds,
        "second_difference_check": second,
        "weighting_systematic": weighting,
        "error_budget": budget,
        "candidate_sizes": candidates,
        "acquisition": acquisition(candidates, second),
        "decision": decide(reproduction, fit, folds, budget, second, weighting),
        "not_established": [
            "The exponent is not identified with any named percolation "
            "correction-to-scaling exponent. omega is measured in the site count "
            "N, so N^-1 is L^-2 on a square torus; that is a coincidence of "
            "numbers until something independent fixes the operator.",
            "The second-difference discrepancy is measured, not explained. Two "
            "exponents would produce it; so would a step-size-dependent bias in "
            "the quantile reconstruction that the first differences happen to "
            "cancel. Only a third three-size lineage separates those.",
            "The two three-size lineages are NOT structurally independent at "
            "their third rung. Under the primary weighting both extrapolate "
            "there (325 and 425 are the only extrapolating spin-0 combinations "
            "in the tree) and both carry a 5x larger per-batch sample count "
            "there. So the closeness of their two ratios under that weighting "
            "is not evidence on its own; what carries weight is that the "
            "discrepancy survives the equal weighting, which never "
            "extrapolates, with a wider spread (1.44 to 1.58 over all four "
            "cells rather than 1.54 to 1.56 over two).",
            "The five transitions are one correlated evidence block with #582 "
            "and with Gate 3. They are the same histograms read three ways.",
        ],
    }


def render(report: Mapping[str, Any]) -> str:
    lines = [f"#{report['issue']} -- the amplitude law of #582's frozen shape direction",
             ""]
    control = report["reproduction_control"]
    lines.append("reproduction control: "
                 + ("every affine statistic matches #582 exactly"
                    if control["every_statistic_matches_exactly"] else "FAILED"))
    lines.append("")
    lines.append(f"{'transition':>10} {'sbar':>8} {'h':>7} {'amplitude':>13} "
                 f"{'se':>10} {'z':>8}")
    for point in report["amplitudes"]:
        lines.append(f"{point['label']:>10} {point['sbar']:8.5f} {point['h']:7.5f} "
                     f"{point['amplitude']:+13.6e} {point['standard_error']:10.3e} "
                     f"{point['significance']:8.1f}")
    fit = report["one_exponent_fit"]
    lines += ["", f"one exponent: omega = {fit['omega']:.4f}, "
                  f"lambda = {fit['scale_amplitude']:.5e}, "
                  f"chi2 = {fit['statistic']:.1f} on {fit['degrees_of_freedom']} df",
              f"  largest relative misfit {100 * fit['largest_relative_error']:.2f}%"]
    folds = report["leave_one_transition_out"]
    lines.append(f"  held-out amplitude prediction worst case "
                 f"{100 * folds['largest_relative_prediction_error']:.2f}%, "
                 f"omega across folds {folds['omega_range'][0]:.4f}"
                 f"-{folds['omega_range'][1]:.4f}")
    budget = report["error_budget"]
    lines.append(f"  omega = {budget['omega']:.3f} +- {budget['combined_uncertainty']:.3f} "
                 f"(loo {budget['leave_one_out_half_spread']:.3f}, "
                 f"weighting {budget['orientation_weighting_shift']:.3f}); "
                 f"unity excluded: {budget['unit_exponent_is_excluded']}")
    lines += ["", "independent check -- second divided difference, parameters frozen:"]
    for item in report["second_difference_check"]["per_lineage"]:
        if not item["usable"]:
            lines.append(f"  {item['lineage']:<12} {item['reason']}")
            continue
        lines.append(f"  {item['lineage']:<12} measured {item['measured']:+.5e} "
                     f"predicted {item['predicted_from_the_first_differences']:+.5e} "
                     f"ratio {item['ratio']:.4f}")
    check = report["second_difference_check"]
    if check["relative_spread_between_lineages"] is not None:
        lines.append(f"  the two ratios agree to "
                     f"{100 * check['relative_spread_between_lineages']:.1f}% "
                     f"under the primary weighting alone")
    weighting = report["weighting_systematic"]
    lines.append("  all four (weighting x lineage) cells, the honest spread:")
    for cell, ratio in sorted(weighting["second_difference_ratio_cells"].items()):
        lines.append(f"    {cell:<24} {ratio:.4f}")
    lines.append(f"    range {weighting['second_difference_ratio_range'][0]:.4f}"
                 f"-{weighting['second_difference_ratio_range'][1]:.4f}, spread "
                 f"{100 * weighting['second_difference_ratio_spread_over_all_cells']:.1f}%"
                 f", survives the weighting change: "
                 f"{weighting['the_discrepancy_survives_the_weighting_change']}")
    lines += ["", f"next production: {report['acquisition']['preferred']}",
              "", f"verdict: {report['decision']['verdict']}"]
    for finding in report["decision"]["findings"]:
        lines.append(f"  - {finding}")
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    report = assemble()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(render(report))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
