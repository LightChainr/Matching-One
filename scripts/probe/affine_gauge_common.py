"""Shared machinery for probe HIGH (#617): the affine gauge of the full law.

Everything here rebuilds on the committed pipeline, never replacing it:

- ``load_all()`` loads the eight committed histogram blocks (plus N=725 from
  PR #614's branch, when present on disk) under both orientation weightings,
  through the production path ``score_wasserstein_shape_flow.load_sizes``.
- ``frozen_direction()`` returns the published consensus ``g`` from the
  committed #582 artifact -- the published vector is the probe's fixed
  reference; nothing here re-fits it.
- ``lineage_identity_at`` is the #612 identity generalised from the historical
  charts to an arbitrary curvature-chart attachment ``Q_lambda``.

Chart bookkeeping (matching the committed pipeline exactly):

- the LOWER first-difference fit lives in ``span{1, Q_base, g}``;
- the UPPER first-difference fit lives in ``span{1, Q_middle, g}``;
- the SECOND-DIFFERENCE (curvature) fit lives in ``span{1, Q_middle, g}`` --
  the "historical middle chart".

The probe's attachment family moves ONLY the curvature chart's middle basis
vector along the lower transition's segment,

    Q_lambda = (1-lambda) Q_base + lambda Q_middle,   lambda in [0, 1],

with the first-difference charts held at their historical anchors.  lambda = 1
is the historical middle chart (the attachment that produced the 1.55 fake
excess); lambda = 0 attaches the curvature to the LOWER chart, where the
curvature and ``a0`` share one basis and the naive difference is consistent.
The two endpoints are P2's expected 1.55 and 1.00.

The general-lambda identity.  GLS amplitudes are linear functionals of the
observation at fixed basis and covariance, so for ANY attachment A:

    a_curv(A) = 2/(h0+h1) * [ a1 - a0
                              + ell_A(r1) - ell_A(r0)
                              - beta0 * e_star(A) ],
    e_star(A) = ell_A(Q_base)  -- the g-coefficient of Q_base in the
                                 curvature chart at attachment A,

with ``a0, a1, r0, r1, beta0`` all at their historical anchors and ``ell_A``
the g-covector of the curvature fit in ``span{1, A, g}``.  This is exact
(zero approximation beyond mpmath): at A = Q_middle it must collapse onto the
committed #612 formula, which P0 asserts to 1e-12; the collapse is itself the
proof that the committed identity is the A = Q_middle member of this family.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

import score_wasserstein_shape_flow as flow  # noqa: E402
import threshold_quantile_lineage as lineage  # noqa: E402
import p582_amplitude_law as law  # noqa: E402

# p612_chart_identity lives on PR #614 (claude/p612-n725-chart), not yet on the
# frontier.  The probe vendors it verbatim at commit f0981a98 into
# scripts/probe/_vendored/ and imports from there; provenance is recorded in
# the atlas.  If the frontier ever merges #614, the vendored copy is a pinned
# reference, not a fork.
_PROBE_DIR = Path(__file__).resolve().parent
if str(_PROBE_DIR / "_vendored") not in sys.path:
    sys.path.insert(0, str(_PROBE_DIR / "_vendored"))
from p612_chart_identity import (  # noqa: E402
    DPS,
    LEVELS,
    gls_fit,
    transition_geometry,
)

N725 = 725
N725_SOURCE = "results/server-20260907/P612-n725-fullcurve/raw/n725_100m.hist.csv"

COMMITTED_FLOW = ROOT / "results" / "wasserstein-shape-flow" / "latest.json"
COMMITTED_LAW = ROOT / "results" / "p582-amplitude-law" / "latest.json"
COMMITTED_N725 = ROOT / "results" / "p612-n725-score" / "latest.json"
COMMITTED_CHART = ROOT / "results" / "p612-chart-identity" / "latest.json"

OUT_DIR = ROOT / "results" / "probe-affine-gauge"


# ---------------------------------------------------------------- loading


def load_all(with_725: bool = False) -> dict[str, dict[int, dict[str, Any]]]:
    """Both weightings for every committed size; 725 added only when asked.

    The production path is deterministic on committed inputs, so the result
    is cached on disk (first run ~100 s, later runs < 1 s).  The cache key is
    the set of sources; if SOURCES ever changes upstream, delete the cache
    file and it rebuilds.  A stale cache is worse than no cache, so the file
    records every source path and the production validates them.
    """
    import pickle
    tag = "all725" if with_725 else "all8"
    cache = OUT_DIR / f"_cache-{tag}.pkl"
    if cache.exists():
        with cache.open("rb") as handle:
            return pickle.load(handle)
    sizes = sorted({size for sizes in flow.LINEAGES.values() for size in sizes})
    if with_725:
        patched = dict(flow.SOURCES)
        patched[N725] = N725_SOURCE
        original = flow.SOURCES
        flow.SOURCES = patched
        try:
            result = flow.load_sizes(sizes + [N725])
        finally:
            flow.SOURCES = original
    else:
        result = flow.load_sizes(sizes)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with cache.open("wb") as handle:
        pickle.dump(result, handle)
    return result


def frozen_direction(weighting: str = "spin0") -> list[float]:
    """The published #582/#584 consensus g, straight from the artifact."""
    committed = json.loads(COMMITTED_LAW.read_text())
    direction = list(committed["frozen_direction"])
    if committed.get("primary_weighting") != weighting:
        raise ValueError(
            f"the committed frozen direction was frozen under "
            f"{committed.get('primary_weighting')!r}, not {weighting!r}; "
            "re-freeze explicitly or pass the matching weighting")
    return direction


def published_amplitudes(weighting: str = "spin0") -> dict[str, float]:
    """The five published first-difference amplitudes of g, by label."""
    committed = json.loads(COMMITTED_LAW.read_text())
    if "weightings" in committed and weighting in committed["weightings"]:
        points = committed["weightings"][weighting]["amplitude_points"]
    else:
        points = committed["amplitudes"]
    return {point["label"]: point["amplitude"] for point in points}


# ------------------------------------------------------------ chart algebra


def attach_point(base_q: Sequence[float], target_q: Sequence[float],
                 lam: float) -> list[float]:
    """Q_lambda on the LOWER transition's segment, lambda in [0, 1].

    ``lambda = 1`` is ``Q_middle`` (the historical curvature chart);
    ``lambda = 0`` is ``Q_base`` (the lower chart).
    """
    return [(1.0 - lam) * b + lam * t for b, t in zip(base_q, target_q)]


def curvature_geometry_at(loaded: Mapping[int, Mapping[str, Any]],
                          sizes: Sequence[int],
                          direction: Sequence[float],
                          lam: float) -> dict[str, Any]:
    """The curvature object scored against ``span{1, Q_lambda, g}``.

    ``sizes`` is a three-size lineage.  The basis middle vector is the
    attachment point ``Q_lambda`` on the LOWER transition's segment;
    everything else (weights on log N, the observation, the covariance) is
    the committed construction, unchanged.
    """
    first, middle, last = sizes
    weights = law.second_difference_weights([math.log(s) for s in sizes])
    observation = [math.fsum(weight * loaded[size]["quantiles"][level]
                             for weight, size in zip(weights, sizes))
                   for level in range(LEVELS)]
    covariance = [[math.fsum(weight * weight * loaded[size]["covariance"][i][j]
                             for weight, size in zip(weights, sizes))
                   for j in range(LEVELS)] for i in range(LEVELS)]
    attach = attach_point(loaded[first]["quantiles"],
                          loaded[middle]["quantiles"], lam)
    fit = gls_fit(observation, covariance,
                  [[1.0] * LEVELS, attach, list(direction)])
    return {
        "sizes": list(sizes), "lambda": lam,
        "amplitude": fit["amplitudes"][2],
        "standard_error": fit["standard_errors"][2],
        "residual": fit["residual"],
        "covectors": fit["covectors"],
        "statistic": fit["statistic"],
        "covariance_rank": fit["covariance_rank"],
        "second_difference_weights": weights,
    }


def lineage_identity_at(loaded: Mapping[int, Mapping[str, Any]],
                        sizes: Sequence[int],
                        direction: Sequence[float],
                        lam: float) -> dict[str, Any]:
    """The #612 identity generalised to curvature-chart attachment lambda.

    Two exact statements live here; both are asserted numerically every call.

    (1) THE COVECTOR EXPANSION (exact at every lambda, 1e-13 verified).  The
    second-difference weights sum to zero, so the curvature observation
    ``y = w0 Q1 + w1 Q2 + w2 Q3`` satisfies, at attachment
    ``A(λ) = (1-λ) Q1 + λ Q2`` with segment ``d = Q2 - Q1``:

        a_curv(λ) = (−w0 λ + w1 (1−λ)) ⟨ell_λ, d⟩
                    + w2 ⟨ell_λ, Q3 − A(λ)⟩,

    because ⟨ell_λ, A(λ)⟩ = 0 by the GLS covector's defining biorthogonality.
    Every attachment dependence of the measured curvature runs through two
    scalar covector pairings — this is the explicit "parallel transport of the
    g-covector" that P1 asks to exhibit, and it is a matrix identity, not a
    metaphor.

    (2) THE COMMITTED /k FORM AS THE λ = 1 MEMBER.  At λ = 1
    (attachment = Q_middle) the committed #612 identity holds:

        a_curv(1) = 2/(h0+h1) [a1 − a0/k + ell(r1) − ell(r0)/k],
        k = 1 + h0 β0,

    and is here recomputed through the same lower/upper fits (never re-fit).
    The probe's λ-sweep therefore interpolates between the historical chart
    (λ=1, where the naive ratio is 1.55) and the consistent chart (λ=0,
    where the curvature shares the LOWER chart) with an exact formula at
    every point between.
    """
    first, middle, last = sizes
    h0 = math.log(middle) - math.log(first)
    h1 = math.log(last) - math.log(middle)
    factor = 2.0 / (h0 + h1)

    lower = transition_geometry(loaded, first, middle, direction)
    upper = transition_geometry(loaded, middle, last, direction)
    curvature = curvature_geometry_at(loaded, sizes, direction, lam)

    a0, a1 = lower["amplitude"], upper["amplitude"]
    beta0 = lower["beta"]
    k = 1.0 + h0 * beta0

    ell_c = curvature["covectors"][2]
    weights = curvature["second_difference_weights"]
    segment = [(t - b) for b, t in zip(loaded[first]["quantiles"],
                                       loaded[middle]["quantiles"])]
    attach = attach_point(loaded[first]["quantiles"],
                          loaded[middle]["quantiles"], lam)

    # (1) covector expansion
    ell_d = _apply(ell_c, segment)
    ell_last_minus_attach = _apply(
        ell_c, [q - a for q, a in zip(loaded[last]["quantiles"], attach)])
    covector_expansion = ((-weights[0] * lam + weights[1] * (1.0 - lam)) * ell_d
                          + weights[2] * ell_last_minus_attach)

    # (2) committed /k form at the historical anchors
    transported = a1 - a0 / k
    residual_transport = (_apply(ell_c, upper["residual"])
                          - _apply(ell_c, lower["residual"]) / k)
    committed_rebuilt = factor * (transported + residual_transport)

    measured = curvature["amplitude"]
    collapse_error = abs(committed_rebuilt - measured)

    return {
        "lineage_sizes": list(sizes),
        "lambda": lam,
        "h0": h0, "h1": h1, "chart_factor": factor,
        "lower_amplitude": a0, "upper_amplitude": a1, "beta0": beta0,
        "k": k,
        "e_star": _apply(ell_c, list(loaded[first]["quantiles"])),
        "second_difference_weights": list(weights),
        "covector_pairings": {
            "ell_lambda_of_d": ell_d,
            "ell_lambda_of_Q3_minus_attach": ell_last_minus_attach,
            "ell_lambda_of_Qbase": _apply(ell_c, list(loaded[first]["quantiles"])),
        },
        "curvature_measured": measured,
        "curvature_standard_error": curvature["standard_error"],
        "identity": {
            "covector_expansion": covector_expansion,
            "expansion_absolute_error": abs(covector_expansion - measured),
            "committed_k_form_at_lambda_1": committed_rebuilt,
            "committed_collapse_absolute_error": collapse_error,
            "naive_no_transport": factor * (a1 - a0),
            "ratio_measured_over_naive": (
                measured / (factor * (a1 - a0))
                if abs(factor * (a1 - a0)) > 0 else None),
        },
    }


def _apply(covector: Sequence[float], vector: Sequence[float]) -> float:
    return math.fsum(left * right for left, right in zip(covector, vector))


# ------------------------------------------------------------------ output


def dump(program: str, payload: Mapping[str, Any]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / f"{program}.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n",
                    encoding="utf-8")
    return path
