#!/usr/bin/env python3
"""Per-batch quantile vectors of the finite threshold law, for issue #582.

#582 needs the quantile function of the reconstructed threshold law, one vector
per batch, so that an aligned delete-one jackknife can give the covariance of the
whole grid rather than of one level at a time.  This module supplies that, and
nothing else: no shape flow is computed here.

**The law.**  ``analysis/threshold_histogram_profile_contract.json`` fixes the
channel: a transition at the k-th occupied site is Beta(k, N+1-k), the mixture
weights are the equal-weight average of the K_minus and K_plus rank laws, and

    F(p) = sum_k w_k I_p(k, N+1-k).

That is the ``F_N(p) = [1 + M_N(p)] / 2`` of #582, and this module computes the
same object as ``threshold_histogram_profile`` rather than a second one.

**Why it does not use the exact path in production.**  The exact route builds a
degree-N rational polynomial whose coefficients carry ``N * C(N-1, k-1)``, and
then isolates a root by Sturm sequences.  At N=325, across 100 batches and nine
levels and six sizes, that is not affordable.  The identity

    F(p) = sum_j B(N, j, p) G(j),      G(j) = P(K <= j),

collapses the whole mixture into one binomial pass, because
``I_p(k, N+1-k) = P(Bin(N,p) >= k)`` and the order of summation exchanges.  The
binomial weights are computed in log space and anchored at the mode, which is
the underflow the P49 pilot hit at N=1300 and is fixed here by construction
rather than by a tolerance.

``exactness_control`` checks the fast path against the exact rational one, which
is what #582's second gate asks for.

**Orientation, and why equal weight is wrong.**  Each production carries two
lattice orientations.  Averaging them with equal weight does *not* give the
spin-0 law, because two Gaussian integers of the same norm do not have opposite
``cos 4theta``: at N=65 the pair is ``(+0.8788, -0.4845)`` and the equal-weight
average retains ``+0.1972`` of a spin-4 component.  Worse, that residue
alternates in sign along a lineage and has opposite sign between the two Gaussian
lineages, which is a pattern easily mistaken for a shape flow that differs
between lineages.

So the law is reconstructed twice.  ``Q_bar`` is the equal-weight average, kept
because it is what a reader would build by default and because the difference
between the two is the size of the systematic.  ``Q_spin0`` removes the residue
exactly to first order: writing ``Q_i(u) = Q_0(u) + c_i D(u)`` for the two
orientations, the orientation difference measures ``D = (Q_1 - Q_2)/(c_1 - c_2)``
and

    Q_0 = Q_bar - cbar * (Q_1 - Q_2) / (c_1 - c_2),    cbar = (c_1 + c_2)/2,

which is the linear combination ``w_1 Q_1 + w_2 Q_2`` with
``w_1 = -c_2/(c_1-c_2)``.  Where the two orientations carry the same sign of
``cos 4theta`` -- N=325 and N=425 do -- those weights leave the unit interval and
the correction is an extrapolation rather than a mixture.  That is recorded per
size rather than hidden, because an extrapolated correction is a weaker object
than an interpolated one.

This is a first-order removal of spin 4.  Spin 8 and the second order in ``D``
survive it, and no claim here says otherwise.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

ROOT = Path(__file__).resolve().parents[1]

#: Frozen interior levels.  Nine symmetric deciles, chosen before any lineage was
#: loaded, so that #582's warning against picking a window by whichever value
#: minimises the rank cannot apply to them.
FROZEN_LEVELS: tuple[float, ...] = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

#: Binomial terms are dropped once their log-weight falls this far below the
#: mode.  exp(-745) is the smallest normal double, so 700 discards only terms
#: that cannot change a sum of order one.
LOG_WEIGHT_FLOOR = 700.0

QUANTILE_TOLERANCE = 1e-13


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


ORIENTATIONS = ("first", "second")


def cos_four_theta(a: int, b: int) -> float:
    """``cos 4theta`` for the Gaussian integer ``a + bi``, exactly.

    ``Re(w^4)/|w|^4 = (a^4 - 6 a^2 b^2 + b^4) / (a^2 + b^2)^2``.  This is the
    coefficient that decides how much spin 4 an orientation contributes, and it
    is arithmetic about the period lattice rather than anything measured.
    """
    norm = a * a + b * b
    _require(norm > 0, "the Gaussian integer must be nonzero")
    return (a ** 4 - 6 * a * a * b * b + b ** 4) / (norm * norm)


def load_batch_histograms(path: Path) -> dict[str, Any]:
    """Per-batch, per-orientation K_minus and K_plus rank counts.

    Returns ``{"n": N, "orientation_cos4theta": {...}, "batches": {batch:
    {orientation: {"minus": {...}, "plus": {...}}}}}``.  Orientations are kept
    apart because averaging them with equal weight leaves a spin-4 residue; see
    the module docstring.
    """
    counts: dict[int, dict[str, dict[str, dict[int, int]]]] = defaultdict(
        lambda: defaultdict(lambda: {"minus": defaultdict(int), "plus": defaultdict(int)}))
    sizes: set[int] = set()
    reps: dict[str, tuple[int, int]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            sizes.add(int(row["n"]))
            orientation = row["orientation"]
            _require(orientation in ORIENTATIONS,
                     f"{path} carries an unexpected orientation {orientation!r}")
            rep = (int(row["a"]), int(row["b"]))
            if orientation in reps:
                _require(reps[orientation] == rep,
                         f"{path} gives orientation {orientation!r} two representatives")
            reps[orientation] = rep
            counts[int(row["batch"])][orientation][row["kind"]][int(row["k"])] += int(row["count"])
    _require(len(sizes) == 1, f"{path} mixes site counts: {sorted(sizes)}")
    size = sizes.pop()
    _require(set(reps) == set(ORIENTATIONS), f"{path} is missing an orientation")
    for batch, payload in counts.items():
        _require(set(payload) == set(ORIENTATIONS),
                 f"{path} batch {batch} is missing an orientation")
    return {
        "n": size,
        "orientation_cos4theta": {name: cos_four_theta(*reps[name]) for name in ORIENTATIONS},
        "orientation_representative": {name: list(reps[name]) for name in ORIENTATIONS},
        "batches": {
            batch: {name: {"minus": dict(payload[name]["minus"]),
                           "plus": dict(payload[name]["plus"])}
                    for name in ORIENTATIONS}
            for batch, payload in sorted(counts.items())
        },
    }


def spin_zero_weights(cos4: Mapping[str, float]) -> dict[str, float]:
    """Weights on the two orientations that cancel spin 4 exactly, to first order.

    With ``Q_i = Q_0 + c_i D``, the combination ``w_1 Q_1 + w_2 Q_2`` removes
    ``D`` when ``w_1 c_1 + w_2 c_2 = 0`` and ``w_1 + w_2 = 1``.  When the two
    ``c_i`` share a sign the weights leave ``[0, 1]`` and the correction is an
    extrapolation; ``is_interpolation`` below records that rather than hiding it.
    """
    first, second = cos4["first"], cos4["second"]
    gap = first - second
    _require(abs(gap) > 1e-12,
             "the two orientations carry the same cos4theta; spin 4 is not separable")
    weights = {"first": -second / gap, "second": first / gap}
    return weights


def is_interpolation(weights: Mapping[str, float]) -> bool:
    return all(0.0 <= value <= 1.0 for value in weights.values())


def rank_cdf(minus: Mapping[int, int], plus: Mapping[int, int], n: int) -> list[float]:
    """``G(j) = P(K <= j)`` for the equal-weight K_minus / K_plus mixture, j = 0..n."""
    minus_total = sum(minus.values())
    plus_total = sum(plus.values())
    _require(minus_total > 0 and plus_total > 0, "both births must carry samples")
    _require(minus_total == plus_total,
             f"K_minus and K_plus sample counts differ: {minus_total} vs {plus_total}")
    denominator = 2.0 * minus_total
    cumulative = 0.0
    out = [0.0]
    for rank in range(1, n + 1):
        cumulative += (minus.get(rank, 0) + plus.get(rank, 0)) / denominator
        out.append(cumulative)
    _require(abs(out[-1] - 1.0) < 1e-9, f"mixture weights sum to {out[-1]}, not one")
    return out


def _binomial_weights(n: int, p: float) -> tuple[int, list[float]]:
    """Binomial pmf over a window around the mode, in log space.

    Returns the window's first index and the weights on it.  Anchoring at the
    mode rather than at j = 0 is what stops ``(1-p)**n`` underflowing to exactly
    zero and taking the whole recurrence with it.
    """
    _require(0.0 < p < 1.0, "binomial weight needs an interior probability")
    log_p, log_q = math.log(p), math.log1p(-p)
    log_fact_n = math.lgamma(n + 1)

    def log_pmf(j: int) -> float:
        return (log_fact_n - math.lgamma(j + 1) - math.lgamma(n - j + 1)
                + j * log_p + (n - j) * log_q)

    mode = min(n, max(0, int((n + 1) * p)))
    peak = log_pmf(mode)
    low = mode
    while low > 0 and peak - log_pmf(low - 1) < LOG_WEIGHT_FLOOR:
        low -= 1
    high = mode
    while high < n and peak - log_pmf(high + 1) < LOG_WEIGHT_FLOOR:
        high += 1
    return low, [math.exp(log_pmf(j) - peak) for j in range(low, high + 1)]


def profile_cdf(rank_cumulative: Sequence[float], n: int, p: float) -> float:
    """``F(p) = sum_j B(N, j, p) G(j)``, the threshold law's CDF at ``p``."""
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    low, weights = _binomial_weights(n, p)
    total = 0.0
    mass = 0.0
    for offset, weight in enumerate(weights):
        mass += weight
        total += weight * rank_cumulative[low + offset]
    _require(mass > 0.0, "binomial window carried no mass")
    return total / mass


def quantile(rank_cumulative: Sequence[float], n: int, level: float,
             tolerance: float = QUANTILE_TOLERANCE) -> float:
    """``Q(u)`` by bisection on a CDF that is monotone by construction."""
    _require(0.0 < level < 1.0, "quantile level must be interior")
    low, high = 0.0, 1.0
    while high - low > tolerance:
        middle = 0.5 * (low + high)
        if profile_cdf(rank_cumulative, n, middle) < level:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def quantile_vector(minus: Mapping[int, int], plus: Mapping[int, int], n: int,
                    levels: Sequence[float] = FROZEN_LEVELS) -> list[float]:
    cumulative = rank_cdf(minus, plus, n)
    return [quantile(cumulative, n, level) for level in levels]


def pooled_histograms(batches: Mapping[int, Mapping[str, Mapping[str, Mapping[int, int]]]]
                      ) -> dict[str, dict[str, dict[int, int]]]:
    """The full-sample histogram per orientation, as the sum over every batch."""
    pooled: dict[str, dict[str, dict[int, int]]] = {
        name: {"minus": defaultdict(int), "plus": defaultdict(int)}
        for name in ORIENTATIONS}
    for payload in batches.values():
        for name in ORIENTATIONS:
            for kind in ("minus", "plus"):
                for rank, count in payload[name][kind].items():
                    pooled[name][kind][rank] += count
    return {name: {kind: dict(pooled[name][kind]) for kind in ("minus", "plus")}
            for name in ORIENTATIONS}


def _subtract(total: Mapping[int, int], part: Mapping[int, int]) -> dict[int, int]:
    out = dict(total)
    for rank, count in part.items():
        remaining = out.get(rank, 0) - count
        _require(remaining >= 0, "delete-one removed more samples than were present")
        if remaining:
            out[rank] = remaining
        else:
            out.pop(rank, None)
    return out


def combined_quantiles(per_orientation: Mapping[str, Mapping[str, Mapping[int, int]]],
                       n: int, weights: Mapping[str, float],
                       levels: Sequence[float] = FROZEN_LEVELS) -> list[float]:
    """``sum_i w_i Q_i(u)``: each orientation's quantile vector, then combined.

    The combination is applied to the quantile functions and not to the mixture
    weights, so a weight outside ``[0, 1]`` -- an extrapolation, which N=325 and
    N=425 need -- stays meaningful as a linear correction even though it is not a
    probability mixture.
    """
    vectors = {
        name: quantile_vector(per_orientation[name]["minus"],
                              per_orientation[name]["plus"], n, levels)
        for name in ORIENTATIONS}
    return [sum(weights[name] * vectors[name][index] for name in ORIENTATIONS)
            for index in range(len(levels))]


EQUAL_WEIGHTS = {"first": 0.5, "second": 0.5}


def jackknife_quantiles(loaded: Mapping[str, Any],
                        weights: Mapping[str, float] | None = None,
                        levels: Sequence[float] = FROZEN_LEVELS) -> dict[str, Any]:
    """Full-sample quantile vector, and one delete-one vector per batch.

    The deletion removes batch ``b`` from **every** level and **both**
    orientations at once, which is what makes the resulting covariance a
    covariance of one random object across the grid rather than a collection of
    per-level variances.
    """
    size = loaded["n"]
    batches = loaded["batches"]
    if weights is None:
        weights = EQUAL_WEIGHTS
    pooled = pooled_histograms(batches)
    full = combined_quantiles(pooled, size, weights, levels)
    deleted = {}
    for batch, payload in batches.items():
        reduced = {
            name: {kind: _subtract(pooled[name][kind], payload[name][kind])
                   for kind in ("minus", "plus")}
            for name in ORIENTATIONS}
        deleted[batch] = combined_quantiles(reduced, size, weights, levels)
    return {"n": size, "levels": list(levels), "full": full, "deleted": deleted,
            "batches": len(batches), "weights": dict(weights)}


def jackknife_covariance(full: Sequence[float],
                         deleted: Mapping[int, Sequence[float]]) -> list[list[float]]:
    """Delete-one covariance of the quantile vector, batches aligned."""
    count = len(deleted)
    _require(count > 1, "need at least two batches")
    size = len(full)
    pseudo = [[count * full[j] - (count - 1) * vector[j] for j in range(size)]
              for vector in deleted.values()]
    means = [sum(row[j] for row in pseudo) / count for j in range(size)]
    scale = 1.0 / (count * (count - 1))
    return [[scale * sum((row[i] - means[i]) * (row[j] - means[j]) for row in pseudo)
             for j in range(size)] for i in range(size)]


def exactness_control(minus: Mapping[int, int], plus: Mapping[int, int], n: int,
                      levels: Sequence[float] = FROZEN_LEVELS) -> dict[str, Any]:
    """The fast path against the exact rational one, on the same histogram.

    #582's second gate.  The exact route is the repository's committed
    ``threshold_histogram_profile`` mixture, integrated in ``Fraction``
    arithmetic, and its quantiles are Sturm brackets from
    ``exact_threshold_quantile_certificate``.  The wrong number this control
    stops us believing is a quantile grid whose whole shape flow is an artifact
    of the binomial window or of the bisection tolerance.
    """
    try:  # pragma: no cover - import shape depends on how the caller runs
        from scripts.threshold_histogram_profile import (
            density_coefficients, evaluate_polynomial, integrate_density, mixture_weights)
        from scripts.exact_threshold_quantile_certificate import quantile_bracket
    except ModuleNotFoundError:  # pragma: no cover
        from threshold_histogram_profile import (
            density_coefficients, evaluate_polynomial, integrate_density, mixture_weights)
        from exact_threshold_quantile_certificate import quantile_bracket

    weights = mixture_weights(dict(minus), dict(plus), n)
    exact_cdf = integrate_density(density_coefficients(weights))
    cumulative = rank_cdf(minus, plus, n)

    cdf_error = 0.0
    for index in range(1, 20):
        point = Fraction(index, 20)
        exact_value = float(evaluate_polynomial(exact_cdf, point))
        cdf_error = max(cdf_error, abs(profile_cdf(cumulative, n, float(point)) - exact_value))

    quantile_error = 0.0
    for level in levels:
        left, right = quantile_bracket(exact_cdf, Fraction(level).limit_denominator(10 ** 6))
        fast = quantile(cumulative, n, level)
        inside = float(left) - 1e-9 <= fast <= float(right) + 1e-9
        quantile_error = max(quantile_error,
                             0.0 if inside else min(abs(fast - float(left)),
                                                    abs(fast - float(right))))
    return {
        "n": n,
        "largest_cdf_difference": cdf_error,
        "largest_quantile_excursion_outside_the_exact_bracket": quantile_error,
        "levels": list(levels),
    }
