#!/usr/bin/env python3
"""#595 v2 / #596 Gate 2: projective channel design for the Smith nuisance.

#595 v1 maximized ``|A4| / se(A4)``.  That is the wrong objective and this
script replaces it.  The decision #583 waits on is not "is the angular amplitude
large" but "can the angular direction be separated from the Smith/noncyclic
nuisance well enough to bound the induced contamination", and a channel can
improve the first while damaging the second.

Forming ``rho = delta / A4`` during design is also invalid on this block:
``|A4_z|`` runs 0.06 to 2.18, so a delta-method interval on the ratio silently
replaces a confidence set that is genuinely unbounded.  #589 hit that from the
other side -- all ten channel-size cells there have a 2-sigma rho interval
reaching the pole where the quotient response annihilates the fitted ``A4``.

So the design object here is the joint estimator

    y = (A4, delta)

with its full 2x2 covariance, never reduced to a ratio, and the score is the
smallest future sample multiplier ``n`` for which the alpha-level **Fieller**
set for ``rho`` lies wholly inside the declared contamination window

    |A8_fit / A4_fit| = |0.324936 rho / (1 +- 0.771701 rho)| < epsilon,

which also excludes both nuisance poles ``rho = -/+1.29584`` by construction.

WHAT THE PRIMARY SCORE MAY USE (GOVERNANCE 2D, 2E)
--------------------------------------------------
The primary score is computed **under the design null ``delta = 0``**.  That is
not a claim the offset is zero; it is what makes the number a *design* quantity:
it uses the covariance and the angular amplitude and never the measured offset.
A channel therefore cannot be selected because its offset happened to look small
in the archive.

Reported separately, and explicitly not a design score, is the same calculation
at the measured point.  If the measured ``rho`` is real, no sample size reaches
the band -- more sampling would make the offset *measurable*, not small -- and
saying so is the honest complement to the planning number.

A leading-order identity falls out of the Fieller algebra and is worth stating
because it is the whole correction to v1: for large ``n`` the interval
half-width is ``z sqrt(sigma_delta_delta / n) / |A4|``, so the figure of merit is

    A4^2 / var(delta)        not     A4^2 / var(A4).

The optimal linear combination is then a matched filter against the *offset*
covariance, ``v* ∝ S^-1 a``, with value ``a^T S^-1 a`` -- a different filter from
v1's, which whitened by the ``A4`` covariance instead.

Selection is nested three-fold cross-fitting over the archived batches: choose on
two folds, score on the held-out third, three times.  A channel discovered here
is a design recommendation until it is bought prospectively; no Smith verdict may
be read from this block in any channel.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import mpmath as mp

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from analyze_matching_parity_derivatives_fast import combine, obs, remove  # noqa: E402
from p205_channel_leverage import (  # noqa: E402
    BASELINE_CHANNEL,
    BASIS_CHANNELS,
    COMBINATION_FLOOR,
    P_GRID_POINTS,
    P_HALF_WIDTH,
    P_REF,
    baseline_key,
    probability_grid,
    symmetric_eigen,
)
from score_p205_derivative_smith_loading import (  # noqa: E402
    GEOMETRY_ORDER,
    angular_functional,
    h4_offset_functional,
    leakage,
    quadratic_form,
    rho_window,
)
from score_p205_norm5_conjugate_coalescence import (  # noqa: E402
    SIZES,
    aligned_rows,
    jackknife_covariance,
    load_contract,
    load_pair,
    pair_spec,
    sha256,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "p205-projective-channel-design" / "latest.json"
SCHEMA = "matching-one.p205-projective-channel-design.v1"
ISSUE = 595

#: Declared contamination bands, unchanged from #589 so the two are comparable.
BANDS: Tuple[str, ...] = ("0.20", "0.05")

#: Two-sided normal quantile for the Fieller set.  Frozen at 95% before any
#: cost was computed; the cost scales as z^2, so changing it later would rescale
#: every reported number by a known factor and must be declared, not tuned.
Z_LEVEL = mp.mpf("1.959963984540054")

#: Sample multiplier search bounds.  A design needing more than 10^6 times the
#: archived block is not a design, and is reported as unreachable rather than as
#: a large number that invites rounding.
MIN_MULTIPLIER = mp.mpf("1e-3")
MAX_MULTIPLIER = mp.mpf("1e6")
BISECTION_STEPS = 200

#: Three-fold nested cross-fitting over the archived production batches.
FOLDS = 3

#: A cost within this factor of the fold minimum counts as "near-optimal" when
#: measuring how broad the optimum is in p.  Breadth matters because a razor-thin
#: peak is unlikely to survive the move to N=650 or N=2210.
BREADTH_FACTOR = mp.mpf(2)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


# --------------------------------------------------------------------------
# The joint (A4, delta) estimator
# --------------------------------------------------------------------------


def joint_estimator(
    point: Mapping[str, Sequence[mp.mpf]],
    covariance: Mapping[str, Sequence[Sequence[mp.mpf]]],
    n: int,
    integer_weights: Sequence[int],
    channels: Sequence[str],
    weights: Optional[Sequence[mp.mpf]] = None,
) -> Dict[str, Any]:
    """``(A4, delta)`` and their 2x2 covariance for one readout or combination.

    ``weights`` combines the basis channels linearly; ``None`` means a single
    channel.  Both functionals are applied to the *same* three correlated
    geometry values, so the cross term is not optional -- it is what distinguishes
    a channel that resolves the angle from one that merely resolves something.
    """

    angular = angular_functional(n)
    offset = h4_offset_functional(n, integer_weights)
    if weights is None:
        _require(len(channels) == 1, "a single readout takes one channel")
        weights = [mp.mpf(1)]
    _require(len(weights) == len(channels), "weight/channel length mismatch")

    combined_point = [
        mp.fsum(weights[j] * point[channels[j]][i] for j in range(len(channels)))
        for i in range(3)
    ]
    combined_cov = [
        [
            mp.fsum(
                weights[a] * weights[b] * covariance[channels[a]][channels[b]][i][j]
                for a in range(len(channels)) for b in range(len(channels))
            )
            for j in range(3)
        ]
        for i in range(3)
    ]
    amplitude = mp.fsum(angular[i] * combined_point[i] for i in range(3))
    delta = mp.fsum(offset[i] * combined_point[i] for i in range(3))
    return {
        "A4": amplitude,
        "delta": delta,
        "var_A4": quadratic_form(angular, combined_cov, angular),
        "var_delta": quadratic_form(offset, combined_cov, offset),
        "cov_A4_delta": quadratic_form(angular, combined_cov, offset),
    }


def fieller_set(
    estimator: Mapping[str, mp.mpf], multiplier: mp.mpf, design_null: bool
) -> Optional[Tuple[mp.mpf, mp.mpf]]:
    """The Fieller confidence set for ``rho = delta / A4`` at ``n x`` the sample.

    Returns the interval, or ``None`` when the set is unbounded -- which is the
    honest answer whenever ``A4`` is not resolved at the declared level, and is
    exactly the case a delta-method interval would have hidden.
    """

    amplitude = estimator["A4"]
    delta = mp.mpf(0) if design_null else estimator["delta"]
    scale = 1 / multiplier
    var_a = estimator["var_A4"] * scale
    var_d = estimator["var_delta"] * scale
    cov = estimator["cov_A4_delta"] * scale
    z2 = Z_LEVEL**2

    a = amplitude**2 - z2 * var_a
    b = -2 * (amplitude * delta - z2 * cov)
    c = delta**2 - z2 * var_d
    if a <= 0:
        return None  # A4 unresolved: the set is unbounded or all of R.
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        return None  # Empty set; treated as uninformative rather than as a pass.
    root = mp.sqrt(discriminant)
    return ((-b - root) / (2 * a), (-b + root) / (2 * a))


def window(band: str) -> Tuple[mp.mpf, mp.mpf]:
    """The rho interval on which BOTH modulus families stay inside the band."""

    lows, highs = [], []
    for family in ("square", "rectangle"):
        low, high = rho_window(mp.mpf(band), family)
        lows.append(low)
        highs.append(high)
    return (max(lows), min(highs))


def poles() -> Tuple[mp.mpf, mp.mpf]:
    values = []
    for family in ("square", "rectangle"):
        leak = leakage(family)
        leak_a4 = mp.mpf(leak["A4"].numerator) / leak["A4"].denominator
        values.append(-1 / leak_a4)
    return (min(values), max(values))


def sample_multiplier(
    estimator: Mapping[str, mp.mpf], band: str, design_null: bool = True
) -> mp.mpf:
    """Smallest ``n`` for which the Fieller set fits inside the band's window.

    Monotone in ``n`` -- more sample shrinks the set and never grows it -- so a
    bisection is exact up to its own tolerance rather than a search over a
    landscape.  ``inf`` means no multiplier inside the declared bounds works,
    which for the measured-point variant is the expected and informative answer.
    """

    low, high = window(band)

    def fits(multiplier: mp.mpf) -> bool:
        interval = fieller_set(estimator, multiplier, design_null)
        if interval is None:
            return False
        return interval[0] >= low and interval[1] <= high

    if fits(MIN_MULTIPLIER):
        return MIN_MULTIPLIER
    if not fits(MAX_MULTIPLIER):
        return mp.mpf("inf")
    lower, upper = MIN_MULTIPLIER, MAX_MULTIPLIER
    for _ in range(BISECTION_STEPS):
        middle = mp.sqrt(lower * upper)
        if fits(middle):
            upper = middle
        else:
            lower = middle
        if upper / lower < mp.mpf("1.000001"):
            break
    return upper


# --------------------------------------------------------------------------
# The channel table on one batch subset
# --------------------------------------------------------------------------


def subset_table(
    rows: Mapping[str, Sequence[Any]], batches: Sequence[int], p: mp.mpf
) -> Dict[str, Any]:
    """Point values and the full cross-channel jackknife covariance at one ``p``.

    The covariance is kept across *channels as well as geometries* because the
    optimal combination needs it; reducing to per-channel blocks first would
    throw away exactly the correlations that make one combination better than
    another.
    """

    combined = {name: combine([rows[name][b] for b in batches]) for name in GEOMETRY_ORDER}
    point_obs = {name: obs(combined[name], p) for name in GEOMETRY_ORDER}
    point = {
        channel: [point_obs[name][channel] for name in GEOMETRY_ORDER]
        for channel in BASIS_CHANNELS
    }
    replicates = []
    for b in batches:
        deleted = {
            name: obs(remove(combined[name], rows[name][b]), p) for name in GEOMETRY_ORDER
        }
        replicates.append([
            deleted[name][channel]
            for channel in BASIS_CHANNELS for name in GEOMETRY_ORDER
        ])
    flat = jackknife_covariance(replicates)
    width = len(GEOMETRY_ORDER)
    covariance = {
        left: {
            right: [
                [flat[a * width + i][b * width + j] for j in range(width)]
                for i in range(width)
            ]
            for b, right in enumerate(BASIS_CHANNELS)
        }
        for a, left in enumerate(BASIS_CHANNELS)
    }
    return {"point": point, "covariance": covariance}


def optimal_weights(
    table: Mapping[str, Any], n: int, integer_weights: Sequence[int]
) -> Optional[List[mp.mpf]]:
    """``v* ∝ S^-1 a``: matched filter against the OFFSET covariance.

    This is the whole correction to #595 v1.  The Fieller half-width at large
    ``n`` is ``z sqrt(var(delta)/n) / |A4|``, so the quantity to maximize is
    ``A4^2 / var(delta)`` -- the angular amplitude against the variance of the
    *offset* estimator.  v1 whitened by ``var(A4)`` instead, which optimizes a
    precision nobody needs.
    """

    channels = list(BASIS_CHANNELS)
    singles = [
        joint_estimator(table["point"], table["covariance"], n, integer_weights, [channel])
        for channel in channels
    ]
    amplitudes = [row["A4"] for row in singles]

    angular = angular_functional(n)
    offset = h4_offset_functional(n, integer_weights)
    size = len(channels)
    offset_covariance = [
        [
            quadratic_form(offset, table["covariance"][channels[a]][channels[b]], offset)
            for b in range(size)
        ]
        for a in range(size)
    ]
    _ = angular
    values, vectors = symmetric_eigen(offset_covariance)
    largest = max(values)
    if largest <= 0:
        return None
    weights = [mp.mpf(0)] * size
    for value, vector in zip(values, vectors):
        if value <= COMBINATION_FLOOR * largest:
            continue
        overlap = mp.fsum(vector[i] * amplitudes[i] for i in range(size))
        for i in range(size):
            weights[i] += vector[i] * overlap / value
    norm = mp.sqrt(mp.fsum(w**2 for w in weights))
    if norm == 0:
        return None
    return [w / norm for w in weights]


def readout_costs(
    table: Mapping[str, Any], n: int, integer_weights: Sequence[int],
    weights: Optional[Sequence[mp.mpf]] = None,
    amplitude_override: Optional[Mapping[str, mp.mpf]] = None,
) -> Dict[str, Any]:
    """Design cost of every basis channel plus one combination, at one ``p``.

    ``amplitude_override`` replaces the fold-local ``A4`` with a full-sample one.
    It exists only for the stabilized cross-fit variant, where the point is to
    cross-fit the covariance structure without letting a 33-batch amplitude
    estimate dominate the comparison.
    """

    out: Dict[str, Any] = {}
    for channel in BASIS_CHANNELS:
        estimator = dict(joint_estimator(
            table["point"], table["covariance"], n, integer_weights, [channel]
        ))
        if amplitude_override is not None and channel in amplitude_override:
            estimator["A4"] = amplitude_override[channel]
        out[channel] = _score(estimator)
    if weights is not None:
        estimator = dict(joint_estimator(
            table["point"], table["covariance"], n, integer_weights,
            list(BASIS_CHANNELS), weights,
        ))
        if amplitude_override is not None and "combination" in amplitude_override:
            estimator["A4"] = amplitude_override["combination"]
        out["combination"] = _score(estimator)
        out["combination"]["weights"] = [mp.nstr(w, 6) for w in weights]
    return out


def full_sample_amplitudes(
    rows: Mapping[str, Sequence[Any]], batches: Sequence[int],
    grid: Sequence[mp.mpf], n: int, integer_weights: Sequence[int],
) -> Dict[str, Dict[str, mp.mpf]]:
    """``A4`` per readout per grid point, from the whole archived block."""

    out: Dict[str, Dict[str, mp.mpf]] = {}
    for p in grid:
        key = mp.nstr(p, 12)
        table = subset_table(rows, batches, p)
        weights = optimal_weights(table, n, integer_weights)
        entry = {
            channel: joint_estimator(
                table["point"], table["covariance"], n, integer_weights, [channel]
            )["A4"]
            for channel in BASIS_CHANNELS
        }
        if weights is not None:
            entry["combination"] = joint_estimator(
                table["point"], table["covariance"], n, integer_weights,
                list(BASIS_CHANNELS), weights,
            )["A4"]
        out[key] = entry
    return out


def _score(estimator: Mapping[str, mp.mpf]) -> Dict[str, Any]:
    """Design quantities only.  ``delta`` never enters the primary score."""

    projective_information = (
        estimator["A4"] ** 2 / estimator["var_delta"]
        if estimator["var_delta"] > 0 else mp.mpf(0)
    )
    row: Dict[str, Any] = {
        "projective_information": projective_information,
        "A4_snr": (
            abs(estimator["A4"]) / mp.sqrt(estimator["var_A4"])
            if estimator["var_A4"] > 0 else mp.mpf(0)
        ),
        "correlation": (
            estimator["cov_A4_delta"]
            / mp.sqrt(estimator["var_A4"] * estimator["var_delta"])
            if estimator["var_A4"] > 0 and estimator["var_delta"] > 0 else mp.mpf(0)
        ),
    }
    for band in BANDS:
        row[f"multiplier_{band}"] = sample_multiplier(estimator, band, design_null=True)
        row[f"multiplier_at_measured_point_{band}"] = sample_multiplier(
            estimator, band, design_null=False
        )
    return row


# --------------------------------------------------------------------------
# Nested three-fold discovery / validation
# --------------------------------------------------------------------------


def folds(count: int) -> List[List[int]]:
    return [[b for b in range(count) if b % FOLDS == r] for r in range(FOLDS)]


def scan_subset(
    rows: Mapping[str, Sequence[Any]], batches: Sequence[int],
    grid: Sequence[mp.mpf], n: int, integer_weights: Sequence[int],
    weights_by_p: Optional[Mapping[str, Sequence[mp.mpf]]] = None,
    amplitudes: Optional[Mapping[str, Mapping[str, mp.mpf]]] = None,
) -> Dict[str, Any]:
    """Design costs across the grid on one batch subset.

    ``weights_by_p`` lets a validation fold be scored with the *discovery* fold's
    combination weights.  Re-optimizing the weights on the validation fold would
    reintroduce exactly the selection this design is cross-fitting away.
    """

    out = {}
    for p in grid:
        key = mp.nstr(p, 12)
        table = subset_table(rows, batches, p)
        weights = (
            list(weights_by_p[key]) if weights_by_p is not None and key in weights_by_p
            else optimal_weights(table, n, integer_weights)
        )
        override = amplitudes.get(key) if amplitudes is not None else None
        out[key] = readout_costs(table, n, integer_weights, weights, override)
    return out


def cheapest(scan: Mapping[str, Any], band: str) -> Tuple[str, str, mp.mpf]:
    best: Tuple[str, str, mp.mpf] = ("", "", mp.mpf("inf"))
    for key, entry in scan.items():
        for name, row in entry.items():
            value = row[f"multiplier_{band}"]
            if value < best[2]:
                best = (name, key, value)
    return best


def cross_fit(
    rows: Mapping[str, Sequence[Any]], count: int, grid: Sequence[mp.mpf],
    n: int, integer_weights: Sequence[int], bands: Sequence[str],
    full_amplitudes: Optional[Mapping[str, Mapping[str, mp.mpf]]] = None,
) -> Dict[str, Any]:
    """Choose on two folds, score on the held-out third; all three rotations.

    Two variants are computed and both are reported, because the raw one has a
    defect worth naming.  The design cost is proportional to
    ``var(delta) / A4^2``, and on a 33-batch fold the ``A4^2`` in the denominator
    is barely resolved -- so the *validation statistic itself* becomes a ratio
    with a weak denominator, which is the exact disease this whole exercise
    exists to avoid.  The stabilized variant takes every ``A4`` from the full
    sample and cross-fits only the covariance.

    That is defensible and its limitation is stated: ``A4`` is a property of the
    experiment, identical for the discovery and validation folds by construction,
    so its fold-level noise adds variance to the comparison without checking
    anything about selection; but sharing it also removes the one route by which
    a channel could be selected for an upward ``A4`` fluctuation.  The stabilized
    variant was introduced **after** seeing the raw variant's instability, so it
    is labelled as such, and the verdict is read conservatively: if the two
    disagree, the disagreement is the finding.
    """

    partition = folds(count)
    anchor = baseline_key(grid)
    rounds: Dict[str, List[Dict[str, Any]]] = {band: [] for band in bands}
    stabilized: Dict[str, List[Dict[str, Any]]] = {band: [] for band in bands}

    for held_out in range(FOLDS):
        discovery = [b for r, fold in enumerate(partition) if r != held_out for b in fold]
        validation = partition[held_out]
        discovery_scan = scan_subset(rows, discovery, grid, n, integer_weights)
        discovery_stable = (
            scan_subset(rows, discovery, grid, n, integer_weights,
                        amplitudes=full_amplitudes)
            if full_amplitudes is not None else discovery_scan
        )
        for band in bands:
            for label, scan, store in (
                ("raw", discovery_scan, rounds),
                ("stabilized", discovery_stable, stabilized),
            ):
                if label == "stabilized" and full_amplitudes is None:
                    continue
                name, key, discovery_cost = cheapest(scan, band)
                weights = None
                if name == "combination":
                    table = subset_table(rows, discovery, mp.mpf(key))
                    weights = optimal_weights(table, n, integer_weights)
                amplitudes = full_amplitudes if label == "stabilized" else None
                validation_scan = scan_subset(
                    rows, validation, [mp.mpf(key)], n, integer_weights,
                    {key: weights} if weights is not None else None,
                    amplitudes=amplitudes,
                )
                baseline_scan = scan_subset(
                    rows, validation, [mp.mpf(anchor)], n, integer_weights,
                    amplitudes=amplitudes,
                )
                held_cost = validation_scan[key][name][f"multiplier_{band}"]
                baseline_cost = baseline_scan[anchor][BASELINE_CHANNEL][f"multiplier_{band}"]
                store[band].append({
                    "held_out_fold": held_out,
                    "chosen_readout": name,
                    "chosen_p": key,
                    "discovery_multiplier": discovery_cost,
                    "validation_multiplier": held_cost,
                    "validation_baseline_multiplier": baseline_cost,
                    "validation_saving_vs_baseline": (
                        baseline_cost / held_cost
                        if held_cost > 0 and mp.isfinite(held_cost) else mp.mpf(0)
                    ),
                    "beats_baseline_held_out": bool(
                        mp.isfinite(held_cost) and mp.isfinite(baseline_cost)
                        and held_cost < baseline_cost
                    ),
                })

    def summarize(store: Mapping[str, List[Dict[str, Any]]], band: str) -> Dict[str, Any]:
        entries = store[band]
        return {
            "band": band,
            "rounds": entries,
            "all_rounds_beat_baseline": bool(entries) and all(
                row["beats_baseline_held_out"] for row in entries
            ),
            "selection_is_stable": len(
                {(row["chosen_readout"], row["chosen_p"]) for row in entries}
            ) == 1,
            "held_out_multiplier_spread": (
                max(float(row["validation_multiplier"]) for row in entries)
                / min(float(row["validation_multiplier"]) for row in entries)
                if entries and all(
                    mp.isfinite(row["validation_multiplier"]) and row["validation_multiplier"] > 0
                    for row in entries
                ) else float("inf")
            ),
        }

    return {
        band: {
            "raw": summarize(rounds, band),
            "stabilized": summarize(stabilized, band) if full_amplitudes is not None else None,
        }
        for band in bands
    }


def breadth(scan: Mapping[str, Any], readout: str, band: str) -> Dict[str, Any]:
    """How wide the near-optimal region in ``p`` is for one readout.

    A razor-thin optimum is a worse design than a slightly costlier broad one:
    the grid moves when the size does, and #583 is at N=650, not N=325.
    """

    costs = {key: entry[readout][f"multiplier_{band}"] for key, entry in scan.items()
             if readout in entry}
    finite = {key: value for key, value in costs.items() if mp.isfinite(value)}
    if not finite:
        return {"near_optimal_points": 0, "grid_points": len(costs), "best_p": None}
    best_key = min(finite, key=lambda key: finite[key])
    floor = finite[best_key] * BREADTH_FACTOR
    near = [key for key, value in finite.items() if value <= floor]
    return {
        "best_p": best_key,
        "best_multiplier": mp.nstr(finite[best_key], 6),
        "near_optimal_points": len(near),
        "grid_points": len(costs),
        "near_optimal_fraction": float(len(near)) / len(costs),
        "near_optimal_p_range": [min(near), max(near)] if near else None,
    }


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------


def _num(value: mp.mpf, digits: int = 8) -> str:
    if not mp.isfinite(value):
        return "inf"
    return mp.nstr(value, digits)


def size_report(
    ca: Any, cb: Any, contract: Mapping[str, Any], grid: Sequence[mp.mpf]
) -> Dict[str, Any]:
    n = ca.n
    integer_weights = contract["weights"]["H4"][n]
    rows = {
        "C": aligned_rows(ca, "first"),
        "A": aligned_rows(ca, "second"),
        "B": aligned_rows(cb, "second"),
    }
    count = len(rows["C"])
    full = list(range(count))
    full_scan = scan_subset(rows, full, grid, n, integer_weights)
    anchor = baseline_key(grid)

    curves = {
        readout: {
            key: {
                "multiplier_0.20": _num(entry[readout]["multiplier_0.20"], 6),
                "multiplier_0.05": _num(entry[readout]["multiplier_0.05"], 6),
                "projective_information": _num(entry[readout]["projective_information"], 6),
                "A4_snr": _num(entry[readout]["A4_snr"], 6),
                "correlation": _num(entry[readout]["correlation"], 5),
            }
            for key, entry in full_scan.items() if readout in entry
        }
        for readout in list(BASIS_CHANNELS) + ["combination"]
    }
    return {
        "N": n,
        "batches": count,
        "baseline": f"{BASELINE_CHANNEL} @ {anchor}",
        "baseline_multiplier": {
            band: _num(full_scan[anchor][BASELINE_CHANNEL][f"multiplier_{band}"], 6)
            for band in BANDS
        },
        "baseline_multiplier_at_measured_point": {
            band: _num(
                full_scan[anchor][BASELINE_CHANNEL][f"multiplier_at_measured_point_{band}"], 6
            )
            for band in BANDS
        },
        "full_sample_curves": curves,
        "breadth": {
            readout: breadth(full_scan, readout, "0.20")
            for readout in list(BASIS_CHANNELS) + ["combination"]
        },
        "cross_fit": cross_fit(
            rows, count, grid, n, integer_weights, BANDS,
            full_sample_amplitudes(rows, full, grid, n, integer_weights),
        ),
    }


def decide(reports: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    """Read #596's Gate-2 branches A / B / C, conservatively.

    Two cross-fit variants are available and they answer slightly different
    questions.  If they disagree, that disagreement is the result -- the honest
    branch is then "not established", not whichever variant is friendlier.
    """

    def variant(label: str) -> Dict[str, Any]:
        transports = all(
            report["cross_fit"]["0.20"][label]["all_rounds_beat_baseline"]
            for report in reports
        )
        savings = [
            mp.mpf(row["validation_saving_vs_baseline"])
            for report in reports
            for row in report["cross_fit"]["0.20"][label]["rounds"]
        ]
        costs = [
            mp.mpf(row["validation_multiplier"])
            for report in reports
            for row in report["cross_fit"]["0.20"][label]["rounds"]
            if mp.isfinite(row["validation_multiplier"])
        ]
        clean = [
            mp.mpf(row["validation_multiplier"])
            for report in reports
            for row in report["cross_fit"]["0.05"][label]["rounds"]
            if mp.isfinite(row["validation_multiplier"])
        ]
        spreads = [
            report["cross_fit"]["0.20"][label]["held_out_multiplier_spread"]
            for report in reports
        ]
        return {
            "all_rounds_beat_baseline": transports,
            "worst_held_out_saving_vs_baseline": _num(
                min(savings) if savings else mp.mpf(0), 6
            ),
            "median_held_out_saving_vs_baseline": _num(
                sorted(savings)[len(savings) // 2] if savings else mp.mpf(0), 6
            ),
            "cheapest_held_out_multiplier_band_0.20": _num(
                min(costs) if costs else mp.mpf("inf"), 6
            ),
            "cheapest_held_out_multiplier_band_0.05": _num(
                min(clean) if clean else mp.mpf("inf"), 6
            ),
            "worst_fold_to_fold_spread": max(spreads) if spreads else float("inf"),
            "selection_stable_everywhere": all(
                report["cross_fit"][band][label]["selection_is_stable"]
                for report in reports for band in BANDS
            ),
        }

    raw = variant("raw")
    stable = variant("stabilized")

    def branch(summary: Mapping[str, Any]) -> str:
        cheapest_value = summary["cheapest_held_out_multiplier_band_0.20"]
        cheapest_number = mp.mpf(cheapest_value) if cheapest_value != "inf" else mp.mpf("inf")
        if not summary["all_rounds_beat_baseline"]:
            return "B"
        if not mp.isfinite(cheapest_number) or cheapest_number > mp.mpf(50):
            return "C"
        return "A"

    raw_branch, stable_branch = branch(raw), branch(stable)
    if raw_branch == stable_branch:
        verdict = {
            "A": "GATE2_A__USE_THE_SELECTED_CHANNEL_PROSPECTIVELY",
            "B": "GATE2_B__IMPROVEMENT_DOES_NOT_TRANSPORT__STOP_OPTIMIZING_ON_THIS_BLOCK",
            "C": "GATE2_C__NUISANCE_UNBOUNDED_AT_REALISTIC_COST__DO_NOT_BUY_N650_AS_A_PURE_A8_EXPERIMENT",
        }[raw_branch]
    else:
        verdict = (
            f"GATE2_NOT_ESTABLISHED__RAW_SAYS_{raw_branch}_STABILIZED_SAYS_{stable_branch}"
            "__THE_VALIDATION_STATISTIC_IS_ITSELF_UNDERPOWERED"
        )

    return {
        "raw_cross_fit": raw,
        "stabilized_cross_fit": stable,
        "raw_branch": raw_branch,
        "stabilized_branch": stable_branch,
        "verdict": verdict,
        "primary_score_is_delta_free": True,
        "why_two_variants": (
            "The design cost is proportional to var(delta)/A4^2, and on a "
            "33-batch fold the A4^2 denominator is barely resolved, so the "
            "validation statistic is itself a weak-denominator ratio.  The "
            "stabilized variant shares the full-sample A4 across discovery and "
            "validation and cross-fits only the covariance.  It was introduced "
            "after seeing the raw variant's instability and is labelled as such."
        ),
        "what_this_cannot_do": (
            "Every primary number here is computed under the design null "
            "delta = 0, so no channel can be selected because its measured "
            "offset looked small.  The measured-point variant is reported "
            "separately and is not a design score.  A channel selected here is a "
            "design recommendation until it is bought prospectively; no Smith "
            "verdict may be read from this archived block in any channel "
            "(GOVERNANCE 2E)."
        ),
        "note_on_poles": (
            "Both declared windows lie strictly inside the poles at "
            "rho = -/+1.29584, so a Fieller set that fits the band excludes the "
            "poles automatically.  Pole exclusion is therefore implied by, not "
            "additional to, the band requirement."
        ),
    }


def assemble(
    pairs: Mapping[int, Any], contract: Mapping[str, Any],
    half_width: str, points: int,
) -> Dict[str, Any]:
    predeclared = (half_width == P_HALF_WIDTH and points == P_GRID_POINTS)
    grid = probability_grid(half_width, points)
    reports = [size_report(pairs[n][0], pairs[n][1], contract, grid) for n in SIZES]
    low, high = window("0.20")
    tight_low, tight_high = window("0.05")
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "status": "projective design scan; design quantities only, no Smith verdict",
        "supersedes": "the |A4|/se(A4) screen in scripts/p205_channel_leverage.py",
        "basis_channels": list(BASIS_CHANNELS),
        "bands": list(BANDS),
        "confidence_level_z": mp.nstr(Z_LEVEL, 10),
        "rho_windows": {
            "0.20": [mp.nstr(low, 10), mp.nstr(high, 10)],
            "0.05": [mp.nstr(tight_low, 10), mp.nstr(tight_high, 10)],
        },
        "nuisance_poles": [mp.nstr(value, 10) for value in poles()],
        "grid": {
            "half_width": half_width,
            "points": points,
            "predeclared": predeclared,
        },
        "figure_of_merit": (
            "A4^2 / var(delta).  The Fieller half-width at large n is "
            "z sqrt(var(delta)/n) / |A4|, so the offset variance is the "
            "denominator that matters -- not var(A4), which #595 v1 used."
        ),
        "by_size": {str(report["N"]): report for report in reports},
        "decision": decide(reports),
        "provenance": {
            "inputs": [
                {
                    "N": run.n,
                    "pair": f"C-{run.partner}",
                    "histogram_sha256": sha256(run.histogram_path),
                    "git_commit": run.metadata["git_commit"],
                }
                for n in SIZES for run in pairs[n]
            ],
        },
    }


def render(payload: Mapping[str, Any]) -> str:
    lines = [
        f"#{ISSUE} v2 projective channel design "
        f"(bands {', '.join(payload['bands'])}; rho window 0.20 = "
        f"{payload['rho_windows']['0.20'][0][:8]} .. {payload['rho_windows']['0.20'][1][:8]})",
        "",
    ]
    for key, report in payload["by_size"].items():
        lines.append(
            f"N = {key}   baseline {report['baseline']}  "
            f"x{report['baseline_multiplier']['0.20']} (band 0.20), "
            f"x{report['baseline_multiplier']['0.05']} (band 0.05)"
        )
        lines.append(
            f"  {'readout':<12} {'best p':>13} {'x0.20':>10} {'x0.05':>10}"
            f" {'A4^2/var(d)':>12} {'A4 snr':>8} {'broad':>7}"
        )
        for readout, curve in report["full_sample_curves"].items():
            spread = report["breadth"][readout]
            if spread["best_p"] is None:
                lines.append(f"  {readout:<12} {'--':>13} {'inf':>10}")
                continue
            best = curve[spread["best_p"]]
            lines.append(
                f"  {readout:<12} {spread['best_p']:>13} {best['multiplier_0.20'][:10]:>10}"
                f" {best['multiplier_0.05'][:10]:>10} {best['projective_information'][:12]:>12}"
                f" {best['A4_snr'][:8]:>8}"
                f" {spread['near_optimal_points']}/{spread['grid_points']:<4}"
            )
        for band in payload["bands"]:
            for label in ("raw", "stabilized"):
                fit = report["cross_fit"][band][label]
                if fit is None:
                    continue
                lines.append(
                    f"  cross-fit [{label}], band {band}  (choose on 2 folds, score on the third)"
                )
                for row in fit["rounds"]:
                    lines.append(
                        f"    hold fold {row['held_out_fold']}: {row['chosen_readout']} @ {row['chosen_p']}"
                        f"  discovery x{_num(row['discovery_multiplier'], 5)}"
                        f"  held-out x{_num(row['validation_multiplier'], 5)}"
                        f"  vs baseline x{_num(row['validation_baseline_multiplier'], 5)}"
                        f"  {'BEATS' if row['beats_baseline_held_out'] else 'does not beat'}"
                    )
                lines.append(
                    f"    selection stable: {fit['selection_is_stable']}"
                    f"   fold-to-fold cost spread: x{fit['held_out_multiplier_spread']:.4g}"
                )
        lines.append("")
    decision = payload["decision"]
    for label in ("raw_cross_fit", "stabilized_cross_fit"):
        summary = decision[label]
        lines.append(
            f"{label:<22} all rounds beat baseline: {str(summary['all_rounds_beat_baseline']):<5}"
            f"  worst saving x{summary['worst_held_out_saving_vs_baseline']}"
            f"  median saving x{summary['median_held_out_saving_vs_baseline']}"
            f"  cheapest x{summary['cheapest_held_out_multiplier_band_0.20']} (0.20)"
            f" / x{summary['cheapest_held_out_multiplier_band_0.05']} (0.05)"
            f"  spread x{summary['worst_fold_to_fold_spread']:.4g}"
        )
    lines.append(f"VERDICT: {decision['verdict']}")
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
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
    parser.add_argument("--dps", type=int, default=40)
    parser.add_argument("--half-width", default=P_HALF_WIDTH)
    parser.add_argument("--points", type=int, default=P_GRID_POINTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    mp.mp.dps = args.dps

    contract = load_contract(args.prediction, args.experiment)
    loaded = {}
    for n, partner, histogram, moments, metadata in args.pair:
        loaded[(n, partner)] = load_pair(n, partner, histogram, moments, metadata, contract)
    expected = {(n, partner) for n in SIZES for partner in ("A", "B")}
    if set(loaded) != expected:
        raise SystemExit(f"pairs must be exactly {sorted(expected)}")
    pairs = {n: (loaded[(n, "A")], loaded[(n, "B")]) for n in SIZES}

    payload = assemble(pairs, contract, args.half_width, args.points)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(render(payload))
    print()
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
