#!/usr/bin/env python3
"""#589 option 2: calibrate the Smith nuisance in the derivative channel.

#589 shows that the N=650 three-angle design changes two coordinates at once --
the angle and the quotient/Smith class -- so the third fitted coefficient is not
a clean ``A8`` unless the noncyclic loading is separately controlled.  #591
quantifies the leakage exactly: a unit offset ``delta`` on the one noncyclic row
maps into the saturated ``(C, A4, A8)`` as

    square:     (+0.45933, +0.77170, +0.32494) delta
    rectangle:  (+0.45933, -0.77170, +0.32494) delta

so a quotient-only response shows up as a same-sign apparent ``H8`` in both
modulus families.  The nuisance therefore has to be measured before the third
coefficient is read.

#205 already ran the closest control the repository has: at N=325 and N=425 it
put one *noncyclic* same-N ``C`` node against two cyclic ``A``/``B`` nodes under
a common priority field, and scored a frozen one-harmonic affine null.  It found
no gross breakdown -- but it scored the **fixed-p scalar** ``M``, and #583 is
motivated by a derivative/slope channel.  A quotient response can be small in
``M`` and large in ``M'``; nothing in the fixed-p result rules that out.

This script rescores the identical frozen residuals in the derivative channels
``Sp`` and ``Dp`` at the same frozen ``p_ref``, with the same aligned 3x3
delete-one covariance and the same integer nulls.  Nothing is refitted, no new
sample is drawn, and no harmonic is selected.

**This is the same #205 evidence block viewed in another channel, not a new
independent vote** (GOVERNANCE 2E).  It cannot move the H4 verdict in either
direction.  Its only job is to say whether the noncyclic row loads the
derivative channels enough to contaminate an ``A8`` reading at N=650.

The output that matters is dimensionless:

    rho = delta / A4      (noncyclic offset per unit angular amplitude,
                           both measured in the same channel at the same size)

because ``delta`` in a derivative channel at N=325 and an ``A8`` at N=650 have
neither the same units nor the same size.  Under the *assumption* that ``rho``
transfers -- stated as an assumption, not shown -- the spurious ``A8/A4`` that a
quotient response would inject at N=650 is ``0.32494 * rho``.
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import mpmath as mp

if __package__ in (None, ""):  # ``python3 scripts/score_p205_derivative_smith_loading.py``
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from analyze_matching_parity_derivatives_fast import H, combine, obs, remove  # noqa: E402
from score_p205_norm5_conjugate_coalescence import (  # noqa: E402
    EXPECTED,
    GEOMETRY_ORDER,
    MODEL_ORDER,
    SIZES,
    aligned_rows,
    jackknife_covariance,
    load_contract,
    load_pair,
    pair_spec,
    residual_score,
    score_size,
    sha256,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "p205-derivative-smith-loading" / "latest.json"
SCHEMA = "matching-one.p205-derivative-smith-loading.v1"
ISSUE = 589

#: ``M`` is included to reproduce the published #205 score exactly, as the
#: control that this pipeline is the same pipeline.  ``Sp`` and ``Dp`` are the
#: derivative channels #583 actually cares about.  ``S`` is the fixed-p even
#: combination.  ``D`` is carried for a different reason: ``obs`` defines
#: ``D = (rg - rh)/2`` and ``M = rg - rh``, so ``D`` is *exactly* ``M/2`` and
#: must therefore return a bit-for-bit identical ``rho`` -- a free scale
#: invariance check on the ratio machinery, asserted below rather than assumed.
CHANNELS: Tuple[str, ...] = ("M", "S", "D", "Sp", "Dp")
CONTROL_CHANNEL = "M"

#: #591's N=650 leakage of a unit noncyclic offset into the saturated three
#: coefficient fit, as exact rationals.  Recomputed here from scratch rather
#: than quoted, because the whole bridge from this measurement to #583 rests on
#: these three numbers.
EXPECTED_LEAKAGE = {
    "C": Fraction(477911411, 1040449536),
    "A4": Fraction(946294375, 1226244096),
    "A8": Fraction(11156640625, 34334834688),
}

N650_SQUARE = ((25, 5), (23, 11), (19, 17))
N650_RECTANGLE = ((18, 1), (17, 6), (15, 10))
#: The noncyclic row in each family, by index into the tuples above.
N650_NONCYCLIC_INDEX = {"square": 0, "rectangle": 2}

#: The N=325 / N=425 orientations as Gaussian generators, in ``GEOMETRY_ORDER``.
#: Taken from the frozen #205 contract (``EXPECTED``) rather than restated, so a
#: change of geometry there breaks this script instead of silently rescaling it.
GEOMETRY_GENERATORS = {n: tuple(EXPECTED[n][name] for name in GEOMETRY_ORDER) for n in SIZES}

#: Pre-registered reading of the *induced* spurious ``|A8/A4|`` at N=650.  These
#: bands were written before any derivative-channel number was computed; the
#: only #205 numbers already seen at that moment are the published ``M``-channel
#: scores, which is why ``M`` is carried through only as a reproduction control
#: and is excluded from the verdict (GOVERNANCE 2C, 2E).
CLEAN_BAND = mp.mpf("0.05")
BANDED_BAND = mp.mpf("0.20")
#: Channels the verdict is actually read from.
VERDICT_CHANNELS: Tuple[str, ...] = ("Sp", "Dp")
#: Two-sided coverage multiple used for the reported bound.
SIGMA_MULTIPLE = mp.mpf(2)

#: Channels that are exact positive rescalings of one another.  ``rho`` is a
#: ratio of two linear functionals of the same channel, so it must be invariant;
#: if these ever disagreed, the covariance and the functionals would have fallen
#: out of alignment with each other.
RESCALING_PAIRS = (("M", "D"),)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


# --------------------------------------------------------------------------
# Exact recomputation of the N=650 leakage
# --------------------------------------------------------------------------


def cos_harmonic(a: int, b: int, spin: int) -> Fraction:
    """``cos(spin * theta)`` as an exact rational, from ``(a+ib)^spin``."""

    _require(spin > 0 and spin % 4 == 0, "spin must be a positive multiple of four")
    real, imaginary = 1, 0
    for _ in range(spin):
        real, imaginary = real * a - imaginary * b, real * b + imaginary * a
    return Fraction(real, (a * a + b * b) ** (spin // 2))


def invert3(matrix: Sequence[Sequence[Fraction]]) -> List[List[Fraction]]:
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    _require(determinant != 0, "harmonic design matrix is singular")
    adjugate = [
        [e * i - f * h, c * h - b * i, b * f - c * e],
        [f * g - d * i, a * i - c * g, c * d - a * f],
        [d * h - e * g, b * g - a * h, a * e - b * d],
    ]
    return [[value / determinant for value in row] for row in adjugate]


def leakage(family: str) -> Dict[str, Fraction]:
    """Where a unit offset on the one noncyclic row lands in ``(C, A4, A8)``.

    The design matrix has rows ``[1, cos4, cos8]`` per orientation; its inverse
    column for the noncyclic row *is* the leakage.  Recomputed in exact rational
    arithmetic rather than read from #591, because everything downstream of this
    measurement is multiplied by these three numbers.
    """

    rows = N650_SQUARE if family == "square" else N650_RECTANGLE
    index = N650_NONCYCLIC_INDEX[family]
    design = [
        [Fraction(1), cos_harmonic(a, b, 4), cos_harmonic(a, b, 8)] for a, b in rows
    ]
    inverse = invert3(design)
    return {
        "C": inverse[0][index],
        "A4": inverse[1][index],
        "A8": inverse[2][index],
    }


def verify_leakage() -> Dict[str, Any]:
    square = leakage("square")
    rectangle = leakage("rectangle")
    agree = (
        square["C"] == rectangle["C"] == EXPECTED_LEAKAGE["C"]
        and abs(square["A4"]) == abs(rectangle["A4"]) == EXPECTED_LEAKAGE["A4"]
        and square["A8"] == rectangle["A8"] == EXPECTED_LEAKAGE["A8"]
        and square["A4"] == -rectangle["A4"]
    )
    return {
        "square": {key: str(value) for key, value in square.items()},
        "rectangle": {key: str(value) for key, value in rectangle.items()},
        "square_float": {key: float(value) for key, value in square.items()},
        "matches_issue_591": agree,
        "note": (
            "The A4 leakage flips sign between the two modulus families and the "
            "A8 leakage does not, which is why a quotient-only response reads as "
            "a same-sign apparent H8 in both."
        ),
    }


# --------------------------------------------------------------------------
# Channel-parameterized rescore of the frozen #205 block
# --------------------------------------------------------------------------


def score_size_channel(ca: Any, cb: Any, p_ref: mp.mpf, channel: str) -> Dict[str, Any]:
    """``score_size`` with the observable channel left open.

    Every guard in the published scorer still runs -- the common-field
    signature, the byte-identity of the duplicated ``C`` stream, the batch grid.
    Only the scalar read off each combined histogram changes.  The batch
    alignment of the delete-one jackknife is preserved, so the 3x3 covariance is
    the same block structure the published score used, just in another channel.
    """

    _require(channel in CHANNELS, f"unknown channel {channel!r}")
    if ca.n != cb.n or ca.partner != "A" or cb.partner != "B":
        raise ValueError("score_size_channel requires aligned C-A and C-B runs")
    base = score_size(ca, cb, p_ref)  # re-runs every provenance guard
    geometry_rows = {
        "C": aligned_rows(ca, "first"),
        "A": aligned_rows(ca, "second"),
        "B": aligned_rows(cb, "second"),
    }
    combined = {name: combine(rows) for name, rows in geometry_rows.items()}
    point = {name: obs(combined[name], p_ref)[channel] for name in GEOMETRY_ORDER}
    deleted = []
    for batch in range(len(geometry_rows["C"])):
        deleted.append([
            obs(remove(combined[name], geometry_rows[name][batch]), p_ref)[channel]
            for name in GEOMETRY_ORDER
        ])
    covariance = jackknife_covariance(deleted)
    return {
        "N": ca.n,
        "channel": channel,
        "point": point,
        "covariance": covariance,
        "common_C_histogram_rows_sha256": base["common_C_histogram_rows_sha256"],
        "common_C_moments_rows_sha256": base["common_C_moments_rows_sha256"],
    }


# --------------------------------------------------------------------------
# The dimensionless loading
# --------------------------------------------------------------------------


def angular_functional(n: int) -> List[mp.mpf]:
    """The ``(C, A, B)`` functional whose value is the two-cyclic-row ``A4``.

    With two cyclic rows a one-harmonic model ``x = C + A4 cos(4 theta)`` is
    saturated, so ``A4 = (x_A - x_B) / (cos4_A - cos4_B)``.  The noncyclic ``C``
    row gets weight zero on purpose: the denominator of ``rho`` must not itself
    be contaminated by the offset the numerator is trying to measure.
    """

    (_, a_gen, b_gen) = GEOMETRY_GENERATORS[n]
    gap = cos_harmonic(*a_gen, 4) - cos_harmonic(*b_gen, 4)
    _require(gap != 0, f"N={n}: the two cyclic rows are H4-degenerate")
    scale = mp.mpf(gap.numerator) / gap.denominator
    return [mp.mpf(0), 1 / scale, -1 / scale]


def h4_offset_functional(n: int, integer_weights: Sequence[int]) -> List[mp.mpf]:
    """The frozen H4 null, normalized so the ``C`` coefficient is exactly one.

    ``residual_score`` already normalizes this way, which is why its ``residual``
    field *is* ``delta``: the residual reads ``x_C`` minus the value the two
    cyclic rows predict for the noncyclic row under the same one-harmonic model.
    Recomputed here so the identity can be checked rather than assumed.
    """

    scale = Fraction(1, int(integer_weights[0]))
    return [
        mp.mpf(weight * scale.numerator) / scale.denominator
        for weight in integer_weights
    ]


def check_h4_null_is_the_cos4_interpolant(n: int, integer_weights: Sequence[int]) -> mp.mpf:
    """How far the frozen integer null is from the exact ``cos4`` interpolation.

    If this is not zero the whole reading of ``residual`` as a *noncyclic offset*
    collapses: the residual would be a mixture of an offset and a mis-specified
    angular model, and ``rho`` would not be dimensionless in the intended sense.
    """

    generators = GEOMETRY_GENERATORS[n]
    c4 = [cos_harmonic(a, b, 4) for a, b in generators]
    lambda_a = (c4[0] - c4[2]) / (c4[1] - c4[2])
    lambda_b = (c4[1] - c4[0]) / (c4[1] - c4[2])
    exact = [Fraction(1), -lambda_a, -lambda_b]
    scale = Fraction(1, int(integer_weights[0]))
    frozen = [Fraction(weight) * scale for weight in integer_weights]
    return max(abs(mp.mpf(str(float(x - y)))) for x, y in zip(exact, frozen))


def quadratic_form(
    left: Sequence[mp.mpf], covariance: Sequence[Sequence[mp.mpf]],
    right: Sequence[mp.mpf],
) -> mp.mpf:
    return mp.fsum(
        left[i] * covariance[i][j] * right[j] for i in range(3) for j in range(3)
    )


def smith_loading(
    scored: Mapping[str, Any], integer_weights: Sequence[int]
) -> Dict[str, Any]:
    """``rho = delta / A4`` with a delta-method standard error.

    ``delta`` and ``A4`` are two linear functionals of the *same* three
    correlated numbers, so their covariance is not optional: the two functionals
    share the ``A`` and ``B`` rows and the common priority field correlates all
    three.  Ignoring the cross term would misstate the error in either direction
    depending on the sign of the correlation.
    """

    n = scored["N"]
    point = [scored["point"][name] for name in GEOMETRY_ORDER]
    covariance = scored["covariance"]
    offset_weights = h4_offset_functional(n, integer_weights)
    angular_weights = angular_functional(n)

    delta = mp.fsum(offset_weights[i] * point[i] for i in range(3))
    amplitude = mp.fsum(angular_weights[i] * point[i] for i in range(3))
    var_delta = quadratic_form(offset_weights, covariance, offset_weights)
    var_amplitude = quadratic_form(angular_weights, covariance, angular_weights)
    cross = quadratic_form(offset_weights, covariance, angular_weights)
    _require(var_delta > 0 and var_amplitude > 0, f"N={n}: degenerate covariance")

    rho = delta / amplitude
    variance = (
        var_delta / amplitude**2
        - 2 * delta * cross / amplitude**3
        + delta**2 * var_amplitude / amplitude**4
    )
    # The delta method can return a negative variance when the ratio is badly
    # determined; that is a signal, not a number to take a square root of.
    well_posed = variance > 0
    error = mp.sqrt(variance) if well_posed else mp.mpf("nan")
    bound = abs(rho) + SIGMA_MULTIPLE * error if well_posed else mp.mpf("inf")
    return {
        "N": n,
        "channel": scored["channel"],
        "delta": delta,
        "delta_standard_error": mp.sqrt(var_delta),
        "delta_z": delta / mp.sqrt(var_delta),
        "A4": amplitude,
        "A4_standard_error": mp.sqrt(var_amplitude),
        "A4_z": amplitude / mp.sqrt(var_amplitude),
        "delta_A4_correlation": cross / mp.sqrt(var_delta * var_amplitude),
        "rho": rho,
        "rho_standard_error": error,
        "rho_bound_2sigma": bound,
        "rho_well_posed": bool(well_posed),
    }


def induced_a8_over_a4(rho: mp.mpf, family: str) -> Dict[str, mp.mpf]:
    """Spurious ``A8/A4`` a quotient response of size ``rho * A4`` injects at N=650.

    Under the *assumption* that the dimensionless ``rho`` transfers from N=325 /
    N=425 to N=650 -- an assumption this script states and does not establish --
    the offset on the noncyclic N=650 row is ``delta' = rho * A4_true``.  #591's
    leakage then gives

        A8_fitted = leak_A8 * delta',
        A4_fitted = A4_true + leak_A4 * delta',

    and ``A4_true`` cancels out of the ratio, which is the point of using a
    dimensionless ``rho``.
    """

    leak = leakage(family)
    leak_a8 = mp.mpf(leak["A8"].numerator) / leak["A8"].denominator
    leak_a4 = mp.mpf(leak["A4"].numerator) / leak["A4"].denominator
    denominator = 1 + leak_a4 * rho
    return {
        "leak_A8": leak_a8,
        "leak_A4": leak_a4,
        # ``rho`` at which the *fitted* A4 is annihilated by the quotient
        # response.  A confidence interval that reaches this point does not
        # merely give a large ratio -- it contains the case where N=650's second
        # coefficient is entirely nuisance, so no ratio is defined at all.
        "pole": -1 / leak_a4,
        "induced_A8_over_A4": (
            leak_a8 * rho / denominator if denominator != 0 else mp.mpf("inf")
        ),
    }


def rho_window(band: mp.mpf, family: str) -> Tuple[mp.mpf, mp.mpf]:
    """The ``rho`` interval on which the induced ``|A8/A4|`` stays under ``band``.

    ``induced(rho) = leak_A8 rho / (1 + leak_A4 rho)`` is monotone on the branch
    containing zero, so the constraint is an interval, and it is *asymmetric*:
    an offset that shrinks the fitted A4 buys a larger ratio than one that
    inflates it.  The tighter side is the one a sample-size plan has to clear.
    """

    leak = leakage(family)
    leak_a8 = mp.mpf(leak["A8"].numerator) / leak["A8"].denominator
    leak_a4 = mp.mpf(leak["A4"].numerator) / leak["A4"].denominator
    upper = band / (leak_a8 - band * leak_a4)
    lower = -band / (leak_a8 + band * leak_a4)
    return (min(upper, lower), max(upper, lower))


def sample_multiple(error: mp.mpf, band: mp.mpf) -> mp.mpf:
    """Sample-count multiple that would shrink a 2-sigma window into ``band``.

    Stated assumption, and it is load-bearing: this asks what happens *if the
    true offset is zero*, so that a larger sample shrinks the point estimate
    toward zero along with its error.  If the offset is real and merely
    unresolved, no sample size makes the induced ratio small -- it makes it
    measurable instead.  The multiple is therefore a floor on the cost of a
    clean answer, not a promise of a small one.
    """

    if not mp.isfinite(error) or error <= 0:
        return mp.mpf("inf")
    # Take the tighter (smaller-magnitude) side of the window, over both families.
    reach = min(
        min(abs(edge) for edge in rho_window(band, family))
        for family in ("square", "rectangle")
    )
    return (SIGMA_MULTIPLE * error / reach) ** 2


# --------------------------------------------------------------------------
# Control: reproduce the published #205 M-channel score
# --------------------------------------------------------------------------


def control_scores(
    scored: Mapping[int, Mapping[str, Any]], contract: Mapping[str, Any]
) -> Dict[str, Any]:
    """Re-derive #205's published ``M`` scores through *this* code path.

    The point of the control is narrow and worth stating: it does not confirm
    #205, it confirms that ``score_size_channel`` and ``h4_offset_functional``
    are the published pipeline with one substitution.  If these ``signed_z``
    values did not match the archived ``analysis/score.json``, every derivative
    number below would be measuring a different experiment.
    """

    models = []
    for name in MODEL_ORDER:
        by_size = {}
        chi_square = mp.mpf(0)
        for n in SIZES:
            row = residual_score(
                scored[n][CONTROL_CHANNEL]["point"],
                scored[n][CONTROL_CHANNEL]["covariance"],
                contract["weights"][name][n],
            )
            by_size[str(n)] = {
                "residual": mp.nstr(mp.mpf(row["residual"]), 20),
                "standard_error": row["standard_error"],
                "signed_z": row["signed_z"],
            }
            chi_square += mp.mpf(row["signed_z"]) ** 2
        models.append({
            "name": name,
            "by_size": by_size,
            "joint_chi_square": mp.nstr(chi_square, 20),
            "chi_square_survival_df2": mp.nstr(mp.exp(-chi_square / 2), 18),
        })
    return {"model_order": list(MODEL_ORDER), "models": models}


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------


def _num(value: mp.mpf, digits: int = 12) -> str:
    return mp.nstr(value, digits)


def channel_block(
    scored: Mapping[int, Mapping[str, Any]], contract: Mapping[str, Any], channel: str
) -> Dict[str, Any]:
    rows = {}
    for n in SIZES:
        loading = smith_loading(scored[n][channel], contract["weights"]["H4"][n])
        induced = {
            family: induced_a8_over_a4(loading["rho"], family)
            for family in ("square", "rectangle")
        }
        induced_bound = {
            family: induced_a8_over_a4(
                mp.sign(loading["rho"]) * loading["rho_bound_2sigma"], family
            )
            for family in ("square", "rectangle")
            if loading["rho_well_posed"]
        }
        interval = (
            (loading["rho"] - SIGMA_MULTIPLE * loading["rho_standard_error"],
             loading["rho"] + SIGMA_MULTIPLE * loading["rho_standard_error"])
            if loading["rho_well_posed"] else None
        )
        crosses_pole = {
            family: bool(
                interval is not None
                and interval[0] <= value["pole"] <= interval[1]
            )
            for family, value in induced.items()
        }
        rows[str(n)] = {
            "raw": {
                name: _num(scored[n][channel]["point"][name], 16)
                for name in GEOMETRY_ORDER
            },
            "delta": _num(loading["delta"], 16),
            "delta_standard_error": _num(loading["delta_standard_error"], 12),
            "delta_z": _num(loading["delta_z"], 8),
            "A4": _num(loading["A4"], 16),
            "A4_standard_error": _num(loading["A4_standard_error"], 12),
            "A4_z": _num(loading["A4_z"], 8),
            "delta_A4_correlation": _num(loading["delta_A4_correlation"], 8),
            "rho": _num(loading["rho"], 12),
            "rho_standard_error": _num(loading["rho_standard_error"], 12),
            "rho_bound_2sigma": _num(loading["rho_bound_2sigma"], 12),
            "rho_well_posed": loading["rho_well_posed"],
            "induced_A8_over_A4_at_point": {
                family: _num(value["induced_A8_over_A4"], 12)
                for family, value in induced.items()
            },
            "induced_A8_over_A4_at_2sigma": {
                family: _num(abs(value["induced_A8_over_A4"]), 12)
                for family, value in induced_bound.items()
            },
            "two_sigma_interval_contains_fitted_A4_pole": crosses_pole,
            "samples_multiple_for_banded": _num(
                sample_multiple(loading["rho_standard_error"], BANDED_BAND), 6
            ),
            "samples_multiple_for_clean": _num(
                sample_multiple(loading["rho_standard_error"], CLEAN_BAND), 6
            ),
        }
    return rows


def rescaling_control(blocks: Mapping[str, Any]) -> Dict[str, Any]:
    """Channels that differ by a constant factor must give the same ``rho``."""

    checks = []
    for left, right in RESCALING_PAIRS:
        for n in SIZES:
            a = mp.mpf(blocks[left][str(n)]["rho"])
            b = mp.mpf(blocks[right][str(n)]["rho"])
            checks.append({
                "pair": f"{left}/{right}",
                "N": n,
                "absolute_difference": _num(abs(a - b), 6),
                "agrees": bool(abs(a - b) <= mp.mpf("1e-20") * (1 + abs(a))),
            })
    return {
        "checks": checks,
        "all_agree": all(row["agrees"] for row in checks),
        "why": (
            "obs() defines D = M/2 exactly, so rho -- a ratio of two linear "
            "functionals of the same channel -- must be identical.  A "
            "disagreement would mean the functionals and the jackknife "
            "covariance were built from differently scaled quantities."
        ),
    }


def worst_induced(block: Mapping[str, Any]) -> mp.mpf:
    worst = mp.mpf(0)
    for row in block.values():
        if not row["rho_well_posed"]:
            return mp.mpf("inf")
        for value in row["induced_A8_over_A4_at_2sigma"].values():
            worst = max(worst, abs(mp.mpf(value)))
    return worst


def decide(blocks: Mapping[str, Any]) -> Dict[str, Any]:
    """Read the pre-registered bands.

    The verdict deliberately uses only ``Sp``/``Dp``.  ``M`` was published in
    #205 before this script existed, so folding it into the verdict would score
    a number that was already seen; ``S``/``D`` are reported because they
    separate "the parity split matters" from "differentiation matters", but they
    are not the channel #583 is motivated by.
    """

    per_channel = {name: worst_induced(blocks[name]) for name in VERDICT_CHANNELS}
    worst = max(per_channel.values())
    pole_crossings = [
        f"{name} N{n} {family}"
        for name in VERDICT_CHANNELS
        for n in SIZES
        for family, hit in blocks[name][str(n)][
            "two_sigma_interval_contains_fitted_A4_pole"
        ].items()
        if hit
    ]
    cheapest_banded = min(
        mp.mpf(blocks[name][str(n)]["samples_multiple_for_banded"])
        for name in VERDICT_CHANNELS for n in SIZES
    )
    if worst < CLEAN_BAND:
        verdict = "DERIVATIVE_SMITH_LOADING_SMALL__N650_THIRD_COEFFICIENT_READABLE_AS_A8"
    elif worst < BANDED_BAND:
        verdict = "DERIVATIVE_SMITH_LOADING_BOUNDED__N650_MUST_REPORT_A_CONTAMINATION_BAND"
    else:
        verdict = "DERIVATIVE_SMITH_LOADING_UNCONSTRAINED__N650_NEEDS_AN_INDEPENDENT_QUOTIENT_CONTROL"
    return {
        "verdict_channels": list(VERDICT_CHANNELS),
        "worst_induced_A8_over_A4_at_2sigma": {
            name: _num(value, 12) for name, value in per_channel.items()
        },
        "worst_overall": _num(worst, 12),
        "bands": {
            "clean_below": _num(CLEAN_BAND, 4),
            "banded_below": _num(BANDED_BAND, 4),
        },
        "verdict": verdict,
        "two_sigma_intervals_containing_the_fitted_A4_pole": pole_crossings,
        "cheapest_banded_sample_multiple": _num(cheapest_banded, 6),
        "what_this_cannot_do": (
            "This is the #205 evidence block re-read in another channel, not a "
            "new independent vote (GOVERNANCE 2E).  It cannot move #205's H4 "
            "verdict, and the transfer of rho from N=325/425 to N=650 is an "
            "assumption stated here, not a result shown here."
        ),
    }


def assemble(
    pairs: Mapping[int, Any], contract: Mapping[str, Any],
    prediction_path: Path, experiment_path: Path,
) -> Dict[str, Any]:
    scored: Dict[int, Dict[str, Any]] = {}
    for n in SIZES:
        scored[n] = {
            channel: score_size_channel(pairs[n][0], pairs[n][1], contract["p_ref"], channel)
            for channel in CHANNELS
        }
    blocks = {channel: channel_block(scored, contract, channel) for channel in CHANNELS}
    null_defect = {
        str(n): _num(check_h4_null_is_the_cos4_interpolant(n, contract["weights"]["H4"][n]), 6)
        for n in SIZES
    }
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "status": "re-view of the frozen #205 block in additional channels; no new sampling, no refit",
        "fixed_probability": mp.nstr(contract["p_ref"], 30),
        "geometry_order": list(GEOMETRY_ORDER),
        "channels": list(CHANNELS),
        "control_channel": CONTROL_CHANNEL,
        "generators": {
            str(n): {name: list(EXPECTED[n][name]) for name in GEOMETRY_ORDER}
            for n in SIZES
        },
        "frozen_H4_null_is_the_cos4_interpolant_defect": null_defect,
        "n650_leakage": verify_leakage(),
        "published_M_control": control_scores(scored, contract),
        "exact_rescaling_control": rescaling_control(blocks),
        "by_channel": blocks,
        "decision": decide(blocks),
        "provenance": {
            "prediction": str(prediction_path.relative_to(ROOT)),
            "prediction_sha256": sha256(prediction_path),
            "experiment": str(experiment_path.relative_to(ROOT)),
            "experiment_sha256": sha256(experiment_path),
            "inputs": [
                {
                    "N": run.n,
                    "pair": f"C-{run.partner}",
                    "histogram": str(run.histogram_path),
                    "histogram_sha256": sha256(run.histogram_path),
                    "metadata_sha256": sha256(run.metadata_path),
                    "git_commit": run.metadata["git_commit"],
                }
                for n in SIZES for run in pairs[n]
            ],
        },
    }


def render(payload: Mapping[str, Any]) -> str:
    lines = [
        f"#{ISSUE} P205 derivative-channel Smith loading (p_ref = {payload['fixed_probability'][:14]})",
        "",
        "control: published #205 M-channel scores through this code path",
    ]
    for model in payload["published_M_control"]["models"]:
        z = "  ".join(
            f"N{n}: {model['by_size'][str(n)]['signed_z'][:9]:>9}" for n in SIZES
        )
        lines.append(
            f"  {model['name']:<4} chi2_2 = {model['joint_chi_square'][:8]:>8}"
            f"  p = {model['chi_square_survival_df2'][:7]:>7}   {z}"
        )
    lines += ["", "loading by channel (delta = noncyclic offset, A4 = cyclic-pair amplitude)", ""]
    lines.append(
        f"  {'ch':<3} {'N':>4} {'delta_z':>9} {'A4_z':>9} {'rho':>11} {'+-':>11}"
        f" {'|A8/A4|2s':>10} {'pole?':>6} {'x samples':>10}"
    )
    for channel in payload["channels"]:
        for n in SIZES:
            row = payload["by_channel"][channel][str(n)]
            worst_induced_value = max(
                (abs(mp.mpf(value)) for value in row["induced_A8_over_A4_at_2sigma"].values()),
                default=mp.mpf("inf"),
            )
            pole = "YES" if any(row["two_sigma_interval_contains_fitted_A4_pole"].values()) else "-"
            lines.append(
                f"  {channel:<3} {n:>4} {row['delta_z'][:9]:>9} {row['A4_z'][:9]:>9}"
                f" {row['rho'][:11]:>11} {row['rho_standard_error'][:11]:>11}"
                f" {_num(worst_induced_value, 5)[:10]:>10} {pole:>6}"
                f" {row['samples_multiple_for_banded'][:10]:>10}"
            )
    decision = payload["decision"]
    lines += [
        "",
        "verdict channels: " + ", ".join(decision["verdict_channels"]),
        "worst induced |A8/A4| at 2 sigma: " + ", ".join(
            f"{name} {value}" for name, value in
            decision["worst_induced_A8_over_A4_at_2sigma"].items()
        ),
        f"bands: clean < {decision['bands']['clean_below']}, banded < {decision['bands']['banded_below']}",
        "2-sigma intervals reaching the fitted-A4 pole: "
        + (", ".join(decision["two_sigma_intervals_containing_the_fitted_A4_pole"]) or "none"),
        "cheapest x-samples to reach the banded window (assuming a zero true offset): "
        + decision["cheapest_banded_sample_multiple"],
        "exact-rescaling control (rho(D) == rho(M)): "
        + ("PASS" if payload["exact_rescaling_control"]["all_agree"] else "FAIL"),
        f"VERDICT: {decision['verdict']}",
    ]
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair", action="append", type=pair_spec, required=True)
    parser.add_argument(
        "--prediction", type=Path,
        default=ROOT / "predictions/norm5_conjugate_coalescence_20260829.yaml",
    )
    parser.add_argument(
        "--experiment", type=Path,
        default=ROOT / "experiments/p205_norm5_conjugate_coalescence_20260829.yaml",
    )
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    mp.mp.dps = args.dps

    contract = load_contract(args.prediction, args.experiment)
    loaded = {}
    for n, partner, histogram, moments, metadata in args.pair:
        key = (n, partner)
        if key in loaded:
            raise SystemExit(f"duplicate pair {key}")
        loaded[key] = load_pair(n, partner, histogram, moments, metadata, contract)
    expected_keys = {(n, partner) for n in SIZES for partner in ("A", "B")}
    if set(loaded) != expected_keys:
        raise SystemExit(f"pairs must be exactly {sorted(expected_keys)}")
    pairs = {n: (loaded[(n, "A")], loaded[(n, "B")]) for n in SIZES}

    payload = assemble(pairs, contract, args.prediction, args.experiment)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(render(payload))
    print()
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
