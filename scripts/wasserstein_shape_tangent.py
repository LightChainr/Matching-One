#!/usr/bin/env python3
"""Affine-null controls for the Wasserstein shape tangent of issue #582.

Issue #582 asks whether a finite-size change of a whole threshold law is only a
center/width drift or a genuine change of shape.  In one dimension the quantile
function is an isometric coordinate for W2, so the finite-size displacement is
the function ``v(u) = Q_M(u) - Q_N(u)``, and an infinitesimal affine map
``p -> a p + b`` moves the quantile function inside

    span{ 1, Q_N }  subset  L2(0,1).

The shape flow is what is left after projecting that out.

**The projection is the statistic of #579, not a new one.**  Writing the affine
tangent as a two-column basis ``V = [1, Q_N]`` and asking for the
covariance-weighted distance from ``v`` to ``span(V)`` is exactly
``projective_inference.subspace_residual`` with ``dim(V) = 2`` instead of the
``dim(V) = 1`` of a model ray.  So the machinery already exists, and this module
supplies the basis, the controls that show it behaves, and one thing #582's own
formulation does not have.

**What the weighting buys.**  #582 proposes ``||v_shape||_L2``.  That norm treats
every quantile coordinate alike, and the same issue's tail section notes that a
quantile's uncertainty carries a factor ``1 / f(Q(u))``, so coordinates in a thin
tail are both the noisiest and, in an unweighted norm, potentially the loudest.
Weighting by ``S^+`` fixes that automatically: a coordinate with large variance
contributes little, which is the whole reason the frozen central window ``u_min``
exists in the unweighted formulation.  ``tail_dominance_control`` below exhibits
a case where the unweighted norm is dominated by one tail coordinate while the
weighted statistic is not, so the choice of window stops being load-bearing.

Nothing here touches percolation data.  Every control is synthetic and exact.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable, Sequence

try:  # pragma: no cover - import shape depends on how the script is invoked
    from scripts.projective_inference import subspace_residual
except ModuleNotFoundError:  # pragma: no cover
    from projective_inference import subspace_residual

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "wasserstein-shape-tangent-controls" / "latest.json"
SCHEMA = "matching-one.wasserstein-shape-tangent-controls.v1"
ISSUE = 582
DEFAULT_LEVELS = 9


def quantile_levels(count: int = DEFAULT_LEVELS) -> list[float]:
    """The frozen interior grid ``u = 1/(n+1), ..., n/(n+1)``.

    Interior by construction, so no control here can accidentally rest on a
    quantile at 0 or 1 where the inverse CDF is not finite.
    """
    if count < 3:
        raise ValueError("need at least three levels to have a shape direction at all")
    return [(index + 1) / (count + 1) for index in range(count)]


def affine_tangent(base_quantile: Sequence[float]) -> list[list[float]]:
    """The two directions an affine map on the variable can move a quantile along.

    ``p -> a p + b`` sends ``Q`` to ``a Q + b``, so the tangent at ``Q`` is
    ``span{1, Q}``.  Returning it as a basis, rather than as a formula, is what
    lets the same residual routine score it.
    """
    ones = [1.0] * len(base_quantile)
    return [ones, list(base_quantile)]


def shape_residual(displacement: Sequence[float],
                   covariance: Sequence[Sequence[float]],
                   base_quantile: Sequence[float],
                   extra_directions: Sequence[Sequence[float]] = ()) -> dict[str, Any]:
    """Covariance-weighted shape flow: the part of ``displacement`` outside the affine tangent.

    ``extra_directions`` lets a candidate shape generator be added to the basis,
    which is how ``rank_recovery_control`` asks whether one direction closes the
    residual.
    """
    basis = affine_tangent(base_quantile) + [list(row) for row in extra_directions]
    return subspace_residual(list(displacement), covariance, basis)


def unweighted_shape_norm(displacement: Sequence[float],
                          base_quantile: Sequence[float]) -> float:
    """#582's ``||v_shape||_L2``: the same projection with the identity metric.

    Kept so the two can be compared on the same displacement, which is the point
    of ``tail_dominance_control``.
    """
    size = len(displacement)
    identity = [[1.0 if i == j else 0.0 for j in range(size)] for i in range(size)]
    fit = shape_residual(displacement, identity, base_quantile)
    return math.sqrt(max(fit["statistic"], 0.0))


def _base_quantile(levels: Sequence[float]) -> list[float]:
    """A strictly increasing base quantile with no affine relation to ``u``.

    A base proportional to ``u`` would make ``span{1, Q}`` equal ``span{1, u}``
    and every control below would pass for the wrong reason.
    """
    return [math.tan(math.pi * (level - 0.5) / 2.0) for level in levels]


def _identity(size: int) -> list[list[float]]:
    return [[1.0 if i == j else 0.0 for j in range(size)] for i in range(size)]


def affine_null_controls(levels: Sequence[float]) -> list[dict[str, Any]]:
    """Pure translation, pure scale and both: the shape flow must be exactly zero.

    These are the controls that would catch a projection that is not actually a
    projection.  If any of them returns a nonzero statistic, every shape norm
    the method later reports is contaminated by a location or width change.
    """
    base = _base_quantile(levels)
    covariance = _identity(len(base))
    cases = {
        "pure_translation": [0.37 for _ in base],
        "pure_scale": [0.21 * value for value in base],
        "translation_and_scale": [0.21 * value + 0.37 for value in base],
        "negative_scale": [-1.4 * value for value in base],
    }
    rows = []
    for name, displacement in cases.items():
        fit = shape_residual(displacement, covariance, base)
        rows.append({
            "case": name,
            "statistic": fit["statistic"],
            "degrees_of_freedom": fit["degrees_of_freedom"],
            "unweighted_shape_norm": unweighted_shape_norm(displacement, base),
        })
    return rows


def rank_recovery_control(levels: Sequence[float]) -> list[dict[str, Any]]:
    """Declared shape deformations, and whether the residual finds their rank.

    A method that returns zero on affine displacements is only half tested; it
    must also return something nonzero on a deformation that is genuinely
    outside the tangent, and must fall back to zero once that deformation is
    added to the basis.  The wrong result this catches is a projector that
    annihilates everything.
    """
    base = _base_quantile(levels)
    size = len(base)
    covariance = _identity(size)
    first = [math.sin(2.0 * math.pi * level) for level in levels]
    second = [math.sin(4.0 * math.pi * level) for level in levels]
    rows = []
    for name, generators in (("one_shape_direction", [first]),
                             ("two_shape_directions", [first, second])):
        displacement = [0.21 * base[i] + 0.37
                        + sum(0.5 * generator[i] for generator in generators)
                        for i in range(size)]
        alone = shape_residual(displacement, covariance, base)
        with_one = shape_residual(displacement, covariance, base, [generators[0]])
        with_all = shape_residual(displacement, covariance, base, generators)
        rows.append({
            "case": name,
            "generators": len(generators),
            "statistic_affine_basis_only": alone["statistic"],
            "statistic_after_one_generator": with_one["statistic"],
            "statistic_after_all_generators": with_all["statistic"],
            "degrees_of_freedom_affine_basis_only": alone["degrees_of_freedom"],
        })
    return rows


def tail_dominance_control(levels: Sequence[float]) -> dict[str, Any]:
    """One tail coordinate, noisy because its density is small, and two verdicts.

    The displacement is affine plus a fluctuation confined to the extreme
    quantile, whose variance is inflated by the ``1 / f(Q(u))`` factor #582
    names.  The unweighted L2 shape norm sees a large shape flow.  The
    covariance-weighted statistic sees a residual consistent with its degrees of
    freedom.  The wrong number this control stops us believing is a shape rank
    read off an unweighted norm that a thin tail supplied.
    """
    base = _base_quantile(levels)
    size = len(base)
    inflation = 400.0
    covariance = _identity(size)
    covariance[size - 1][size - 1] = inflation ** 2
    excursion = 1.0 * inflation
    displacement = [0.21 * value + 0.37 for value in base]
    displacement[size - 1] += excursion
    weighted = shape_residual(displacement, covariance, base)
    return {
        "tail_standard_deviation": inflation,
        "tail_excursion_in_units_of_its_own_sigma": excursion / inflation,
        "unweighted_shape_norm": unweighted_shape_norm(displacement, base),
        "weighted_statistic": weighted["statistic"],
        "degrees_of_freedom": weighted["degrees_of_freedom"],
        "weighted_equivalent_sigma": weighted["equivalent_sigma"],
        "what_it_shows": (
            "an affine displacement carrying a one-sigma excursion in the "
            "noisiest quantile produces a large unweighted shape norm and a "
            "weighted statistic consistent with its degrees of freedom, so the "
            "frozen u_min window that the unweighted formulation needs is not "
            "load-bearing for the weighted one"
        ),
    }


def assemble(count: int = DEFAULT_LEVELS) -> dict[str, Any]:
    levels = quantile_levels(count)
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "quantile_levels": levels,
        "base_quantile": _base_quantile(levels),
        "affine_tangent_is_the_579_statistic": (
            "the affine tangent span{1, Q} is a two-column basis for "
            "projective_inference.subspace_residual, the same routine that "
            "scores a model ray with a one-column basis"
        ),
        "affine_null": affine_null_controls(levels),
        "rank_recovery": rank_recovery_control(levels),
        "tail_dominance": tail_dominance_control(levels),
        "not_established": [
            "anything about a percolation threshold law: every control here is synthetic",
            "that a covariance-weighted shape norm removes the need to declare a "
            "quantile window; it removes the need for the window to carry the "
            "conditioning, which is a smaller claim",
            "that the shape subspace of any real lineage is low-dimensional",
        ],
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--levels", type=int, default=DEFAULT_LEVELS)
    args = parser.parse_args(argv)
    payload = assemble(args.levels)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
