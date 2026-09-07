#!/usr/bin/env python3
"""Shape flow of the finite threshold law across existing lineages -- issue #582.

The question, in #582's own words: after projecting out the two-dimensional
affine tangent that a center shift and a width change span, is the finite-size
displacement of the whole threshold law zero, one-dimensional, two-dimensional,
or broad and rotating?

Nothing is sampled here.  Every histogram is a committed production artifact and
every number below is a reanalysis of it.

**The statistic is #579's.**  The affine tangent at ``Q_N`` is ``span{1, Q_N}``,
a two-column basis, so the shape flow is ``projective_inference.subspace_residual``
with ``dim(V) = 2``.  Nine frozen levels leave 7 degrees of freedom.  Weighting by
``S^+`` rather than using an unweighted L2 norm is what keeps a thin tail, whose
uncertainty carries ``1 / f(Q(u))``, from supplying the answer; see
``notes/wasserstein-shape-tangent-controls-20260906.md``.

**The covariance across a transition.**  ``v = Q_M - Q_N`` mixes two productions.
Each size's own ``S`` is an aligned delete-one across its batches, deleting batch
``b`` from every level at once.  The cross term is not assumed away silently: for
every transition this module measures the batch-by-batch correlation between the
two sizes' quantile vectors and reports it, so a coupling would be visible rather
than assumed absent.  Two of the lineages reuse a seed across sizes, which is
exactly why that check exists.

**Held-out prediction.**  A shape generator fitted on four transitions is frozen
and then asked to close the residual of the fifth, scored in that transition's
own covariance and degrees of freedom.  The falsification target is the full
quantile vector, not a selected moment.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path
from typing import Any, Mapping, Sequence

try:  # pragma: no cover - import shape depends on how the script is invoked
    from scripts.projective_inference import subspace_residual
    from scripts import threshold_quantile_lineage as lineage
except ModuleNotFoundError:  # pragma: no cover
    from projective_inference import subspace_residual
    import threshold_quantile_lineage as lineage

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "wasserstein-shape-flow" / "latest.json"
SCHEMA = "matching-one.wasserstein-shape-flow.v1"
ISSUE = 582

#: Committed productions, one per size.  Every one is a square torus at two
#: orientations, reconstructed through the same frozen channel, which is what
#: makes them one flow rather than several.
SOURCES: dict[int, str] = {
    65: "results/server-20260828/P45-root-amplitude/n65.hist.csv",
    85: "results/server-20260828/P45-root-amplitude/n85.hist.csv",
    130: "results/server-20260828/P49-fullcurve-doubling-100m/raw/n130.hist.csv",
    145: "results/server-20260829/P50-n145-n290-fullcurve/raw/n145_100m.hist.csv",
    170: "results/server-20260828/P49-fullcurve-doubling-100m/raw/n170.hist.csv",
    290: "results/server-20260829/P50-n145-n290-fullcurve/raw/n290_100m.hist.csv",
    325: "results/server-20260829/P57-norm5-500m/raw/n325_500m.hist.csv",
    425: "results/server-20260829/P57-norm5-500m/raw/n425_500m.hist.csv",
}

#: The three lineages named in the directive that opened this work.
LINEAGES: dict[str, tuple[int, ...]] = {
    "gaussian_13": (65, 130, 325),
    "gaussian_17": (85, 170, 425),
    "p50": (145, 290),
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def transitions() -> list[dict[str, Any]]:
    """Adjacent pairs within each lineage, with their scale multiplier."""
    out = []
    for name, sizes in LINEAGES.items():
        for base, target in zip(sizes, sizes[1:]):
            out.append({"lineage": name, "base": base, "target": target,
                        "multiplier": target / base,
                        "label": f"{base}->{target}"})
    return out


WEIGHTINGS = ("equal", "spin0")


def load_sizes(sizes: Sequence[int]) -> dict[str, dict[int, dict[str, Any]]]:
    """Quantile vectors and covariance per size, under both orientation weightings.

    ``equal`` is the naive average of the two orientations; ``spin0`` is the
    combination that cancels ``cos 4theta`` exactly to first order.  Both are
    computed because the difference between them is the size of the systematic,
    and because the naive one is what a reader would build by default.
    """
    out: dict[str, dict[int, dict[str, Any]]] = {name: {} for name in WEIGHTINGS}
    for size in sizes:
        loaded = lineage.load_batch_histograms(ROOT / SOURCES[size])
        _require(loaded["n"] == size, f"artifact for {size} carries n={loaded['n']}")
        spin0 = lineage.spin_zero_weights(loaded["orientation_cos4theta"])
        cos4 = loaded["orientation_cos4theta"]
        for name, weights in (("equal", lineage.EQUAL_WEIGHTS), ("spin0", spin0)):
            jack = lineage.jackknife_quantiles(loaded, weights)
            covariance = lineage.jackknife_covariance(jack["full"], jack["deleted"])
            out[name][size] = {
                "source": SOURCES[size],
                "batches": jack["batches"],
                "orientation_cos4theta": cos4,
                "orientation_representative": loaded["orientation_representative"],
                "weights": jack["weights"],
                "net_cos4theta": sum(jack["weights"][k] * cos4[k] for k in cos4),
                "spin4_correction_is_an_interpolation": lineage.is_interpolation(spin0),
                "quantiles": jack["full"],
                "deleted": jack["deleted"],
                "covariance": covariance,
                "standard_errors": [math.sqrt(covariance[i][i])
                                    for i in range(len(jack["full"]))],
            }
    return out


def _pearson(left: Sequence[float], right: Sequence[float]) -> float:
    count = len(left)
    mean_left = sum(left) / count
    mean_right = sum(right) / count
    cov = sum((a - mean_left) * (b - mean_right) for a, b in zip(left, right))
    var_left = sum((a - mean_left) ** 2 for a in left)
    var_right = sum((b - mean_right) ** 2 for b in right)
    if var_left <= 0.0 or var_right <= 0.0:
        return 0.0
    return cov / math.sqrt(var_left * var_right)


def cross_size_coupling(base: Mapping[str, Any], target: Mapping[str, Any]) -> dict[str, Any]:
    """Batch-by-batch correlation between the two sizes' quantile vectors.

    Two of the lineages reuse a seed across sizes.  A shared seed driving
    different-sized systems is not a coupling -- the streams diverge as soon as
    the two runs consume different numbers of draws -- but that is an argument,
    and this measures it instead.  Independent runs give a correlation of order
    ``1/sqrt(batches)``; anything much larger means the cross term of
    ``Var(Q_M - Q_N)`` cannot be dropped.
    """
    shared = sorted(set(base["deleted"]) & set(target["deleted"]))
    _require(len(shared) > 1, "sizes share fewer than two batch labels")
    levels = len(base["quantiles"])
    correlations = []
    for level in range(levels):
        left = [base["deleted"][b][level] for b in shared]
        right = [target["deleted"][b][level] for b in shared]
        correlations.append(_pearson(left, right))
    noise = 1.0 / math.sqrt(len(shared))
    return {
        "batches_compared": len(shared),
        "per_level_correlation": correlations,
        "largest_absolute_correlation": max(abs(value) for value in correlations),
        "correlation_expected_from_noise_alone": noise,
        "cross_term_droppable": max(abs(value) for value in correlations) < 3.0 * noise,
    }


def _add(left: Sequence[Sequence[float]],
         right: Sequence[Sequence[float]]) -> list[list[float]]:
    return [[a + b for a, b in zip(row_a, row_b)] for row_a, row_b in zip(left, right)]


def shape_flow(base: Mapping[str, Any], target: Mapping[str, Any],
               multiplier: float,
               extra_directions: Sequence[Sequence[float]] = ()) -> dict[str, Any]:
    """Residual of the transition displacement against the affine tangent.

    ``extra_directions`` adds candidate shape generators to the basis, which is
    how the held-out prediction asks whether a frozen generator closes what the
    affine tangent alone leaves.
    """
    scale = math.log(multiplier)
    displacement = [(t - b) / scale
                    for b, t in zip(base["quantiles"], target["quantiles"])]
    covariance = [[value / (scale * scale) for value in row]
                  for row in _add(base["covariance"], target["covariance"])]
    basis = [[1.0] * len(displacement), list(base["quantiles"])]
    basis += [list(row) for row in extra_directions]
    fit = subspace_residual(displacement, covariance, basis)
    fitted = [sum(fit["amplitudes"][i] * basis[i][j] for i in range(len(basis)))
              for j in range(len(displacement))]
    return {
        "displacement": displacement,
        "residual": [displacement[j] - fitted[j] for j in range(len(displacement))],
        "statistic": fit["statistic"],
        "degrees_of_freedom": fit["degrees_of_freedom"],
        "equivalent_sigma": fit["equivalent_sigma"],
        "p_value": fit["p_value"],
        "amplitudes": fit["amplitudes"],
        "covariance_rank": fit["covariance_rank"],
        "covariance_condition_number": fit["covariance_condition_number"],
    }


def _normalise(vector: Sequence[float]) -> list[float]:
    norm = math.sqrt(sum(value * value for value in vector))
    _require(norm > 0.0, "cannot normalise a zero direction")
    return [value / norm for value in vector]


def _leading_direction(vectors: Sequence[Sequence[float]],
                       iterations: int = 500) -> list[float]:
    """Leading eigenvector of the sum of outer products of unit residuals.

    Power iteration; the matrix is small and symmetric positive semi-definite.
    """
    size = len(vectors[0])
    gram = [[sum(v[i] * v[j] for v in vectors) for j in range(size)] for i in range(size)]
    current = [1.0 / math.sqrt(size)] * size
    for _ in range(iterations):
        nxt = [sum(gram[i][j] * current[j] for j in range(size)) for i in range(size)]
        norm = math.sqrt(sum(value * value for value in nxt))
        if norm == 0.0:
            break
        current = [value / norm for value in nxt]
    return current


def _deflate(vectors: Sequence[Sequence[float]],
             direction: Sequence[float]) -> list[list[float]]:
    out = []
    for vector in vectors:
        overlap = sum(a * b for a, b in zip(vector, direction))
        out.append([a - overlap * b for a, b in zip(vector, direction)])
    return out


def held_out_prediction(loaded: Mapping[int, Mapping[str, Any]],
                        rows: Sequence[Mapping[str, Any]],
                        held_out: str, rank: int) -> dict[str, Any]:
    """Freeze a shape subspace on every other transition, then score this one.

    The generator is extracted from the *training* residuals only and is not
    refitted on the held-out transition; only its amplitude is free there.  The
    wrong result this guards against is a shape direction that was chosen after
    seeing the transition it then explains.
    """
    training = [row for row in rows if row["label"] != held_out]
    _require(len(training) >= rank, "not enough training transitions for that rank")
    residuals = [_normalise(row["affine_only"]["residual"]) for row in training]
    directions = []
    pool = [list(vector) for vector in residuals]
    for _ in range(rank):
        direction = _leading_direction(pool)
        directions.append(direction)
        pool = _deflate(pool, direction)
    target_row = next(row for row in rows if row["label"] == held_out)
    base = loaded[target_row["base"]]
    target = loaded[target_row["target"]]
    with_generators = shape_flow(base, target, target_row["multiplier"], directions)
    alone = target_row["affine_only"]
    return {
        "held_out": held_out,
        "rank": rank,
        "trained_on": [row["label"] for row in training],
        "generators": [list(direction) for direction in directions],
        "statistic_affine_only": alone["statistic"],
        "degrees_of_freedom_affine_only": alone["degrees_of_freedom"],
        "statistic_with_frozen_generators": with_generators["statistic"],
        "degrees_of_freedom_with_frozen_generators": with_generators["degrees_of_freedom"],
        "equivalent_sigma_with_frozen_generators": with_generators["equivalent_sigma"],
        "generator_amplitudes_on_the_held_out_transition":
            with_generators["amplitudes"][2:],
    }


def pairwise_angles(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Angles between the transitions' unit shape residuals, in degrees.

    A stable one-dimensional shape generator would put every pair near 0 or 180.
    A rotating flow would not.  Sign is not meaningful -- a generator and its
    negative are the same direction -- so the reported angle is folded into
    ``[0, 90]``.
    """
    out = []
    for index, first in enumerate(rows):
        for second in rows[index + 1:]:
            left = _normalise(first["affine_only"]["residual"])
            right = _normalise(second["affine_only"]["residual"])
            cosine = max(-1.0, min(1.0, sum(a * b for a, b in zip(left, right))))
            angle = math.degrees(math.acos(abs(cosine)))
            out.append({"pair": [first["label"], second["label"]],
                        "cosine": cosine,
                        "acute_angle_degrees": angle})
    return out


def same_production_null_control(sizes: Sequence[int] = (65, 145, 290, 325)) -> list[dict[str, Any]]:
    """Two disjoint halves of one production, scored as if they were a transition.

    Nothing changes between them -- same site count, same geometry, same engine --
    so the shape flow must be consistent with zero.  This is the control that
    makes a chi-square of 10^5 on a real transition mean something: without it,
    the number could be the pipeline rather than the physics.  The wrong result
    it stops us believing is a shape flow manufactured by the quantile
    reconstruction, the binomial window, the bisection tolerance or the
    jackknife.
    """
    out = []
    for size in sizes:
        loaded = lineage.load_batch_histograms(ROOT / SOURCES[size])
        weights = lineage.spin_zero_weights(loaded["orientation_cos4theta"])
        labels = sorted(loaded["batches"])
        middle = len(labels) // 2
        halves = []
        for keep in (labels[:middle], labels[middle:]):
            subset = {"n": loaded["n"],
                      "batches": {b: loaded["batches"][b] for b in keep}}
            jack = lineage.jackknife_quantiles(subset, weights)
            halves.append((jack["full"],
                           lineage.jackknife_covariance(jack["full"], jack["deleted"])))
        (first, first_cov), (second, second_cov) = halves
        displacement = [b - a for a, b in zip(first, second)]
        covariance = _add(first_cov, second_cov)
        fit = subspace_residual(displacement, covariance,
                               [[1.0] * len(first), list(first)])
        out.append({
            "n": size,
            "batches_per_half": [middle, len(labels) - middle],
            "displacement_norm": math.sqrt(sum(x * x for x in displacement)),
            "statistic": fit["statistic"],
            "degrees_of_freedom": fit["degrees_of_freedom"],
            "p_value": fit["p_value"],
            "equivalent_sigma": fit["equivalent_sigma"],
        })
    return out


ALTERNATIVE_GRIDS: dict[str, tuple[float, ...]] = {
    "narrow_0.25_to_0.75": tuple(0.25 + 0.5 * i / 8 for i in range(9)),
    "wide_0.05_to_0.95": tuple(0.05 + 0.90 * i / 8 for i in range(9)),
    "fine_15_levels": tuple((i + 1) / 16 for i in range(15)),
}


def level_grid_robustness(sizes: Sequence[int]) -> list[dict[str, Any]]:
    """The same analysis on grids that were not the frozen one.

    #582 warns against choosing a quantile window by whichever value makes the
    rank smallest.  The frozen grid was fixed before any lineage was loaded, and
    this reruns the whole thing on a narrow window, a wide one and a finer one so
    that the reader can see the answer does not depend on it.  These are reported
    alongside the frozen grid, never in place of it.
    """
    out = []
    for name, levels in ({"frozen_9_deciles": lineage.FROZEN_LEVELS}
                         | ALTERNATIVE_GRIDS).items():
        loaded = {}
        for size in sizes:
            raw = lineage.load_batch_histograms(ROOT / SOURCES[size])
            weights = lineage.spin_zero_weights(raw["orientation_cos4theta"])
            jack = lineage.jackknife_quantiles(raw, weights, levels)
            loaded[size] = {
                "quantiles": jack["full"], "deleted": jack["deleted"],
                "covariance": lineage.jackknife_covariance(jack["full"], jack["deleted"])}
        rows = [{**meta,
                 "affine_only": shape_flow(loaded[meta["base"]], loaded[meta["target"]],
                                           meta["multiplier"])}
                for meta in transitions()]
        removed = []
        for row in rows:
            fit = held_out_prediction(loaded, rows, row["label"], 1)
            removed.append(1.0 - fit["statistic_with_frozen_generators"]
                           / fit["statistic_affine_only"])
        statistics = [row["affine_only"]["statistic"] for row in rows]
        out.append({
            "grid": name,
            "levels": list(levels),
            "degrees_of_freedom": rows[0]["affine_only"]["degrees_of_freedom"],
            "affine_statistic_range": [min(statistics), max(statistics)],
            "rank1_held_out_fraction_removed_range": [min(removed), max(removed)],
        })
    return out


def score_one_weighting(loaded: Mapping[int, Mapping[str, Any]]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for meta in transitions():
        base, target = loaded[meta["base"]], loaded[meta["target"]]
        rows.append({
            **meta,
            "coupling": cross_size_coupling(base, target),
            "affine_only": shape_flow(base, target, meta["multiplier"]),
        })
    held_out = [held_out_prediction(loaded, rows, row["label"], rank)
                for row in rows for rank in (1, 2)]
    return {
        "transitions": rows,
        "pairwise_angles": pairwise_angles(rows),
        "held_out": held_out,
        "random_direction_control": random_direction_control(loaded, rows),
    }


def random_direction_control(loaded: Mapping[int, Mapping[str, Any]],
                             rows: Sequence[Mapping[str, Any]],
                             draws: int = 200, seed: int = 20260906) -> list[dict[str, Any]]:
    """What a direction that was never fitted removes from the same residual.

    The held-out numbers only mean something against this.  A frozen generator
    that removes 99 percent of a chi-square is a finding if a random direction
    removes 2 percent, and is nothing if a random direction removes 90.
    """
    generator = random.Random(seed)
    out = []
    for row in rows:
        base, target = loaded[row["base"]], loaded[row["target"]]
        statistics = []
        for _ in range(draws):
            raw = [generator.gauss(0.0, 1.0) for _ in range(len(row["affine_only"]["residual"]))]
            statistics.append(
                shape_flow(base, target, row["multiplier"], [_normalise(raw)])["statistic"])
        statistics.sort()
        affine = row["affine_only"]["statistic"]
        out.append({
            "label": row["label"],
            "draws": draws,
            "statistic_affine_only": affine,
            "median_random": statistics[draws // 2],
            "best_random": statistics[0],
            "median_random_fraction_removed": 1.0 - statistics[draws // 2] / affine,
            "best_random_fraction_removed": 1.0 - statistics[0] / affine,
        })
    return out


def assemble() -> dict[str, Any]:
    sizes = sorted({size for pair in LINEAGES.values() for size in pair})
    loaded = load_sizes(sizes)
    scored = {name: score_one_weighting(loaded[name]) for name in WEIGHTINGS}

    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "levels": list(lineage.FROZEN_LEVELS),
        "lineages": {name: list(sizes) for name, sizes in LINEAGES.items()},
        "orientation_weightings": {
            "equal": "the naive average of the two orientations",
            "spin0": ("the combination cancelling cos4theta exactly to first order; "
                      "primary, because the equal-weight residue alternates in sign "
                      "along a lineage and has opposite sign between the two Gaussian "
                      "lineages, which is a pattern easily mistaken for a shape flow "
                      "that differs between lineages"),
        },
        "primary_weighting": "spin0",
        "sizes": {
            str(size): {
                "source": payload["source"],
                "batches": payload["batches"],
                "orientation_representative": payload["orientation_representative"],
                "orientation_cos4theta": payload["orientation_cos4theta"],
                "equal_weight_net_cos4theta":
                    loaded["equal"][size]["net_cos4theta"],
                "spin0_weights": loaded["spin0"][size]["weights"],
                "spin4_correction_is_an_interpolation":
                    payload["spin4_correction_is_an_interpolation"],
                "quantiles_equal": loaded["equal"][size]["quantiles"],
                "quantiles_spin0": loaded["spin0"][size]["quantiles"],
                "standard_errors_spin0": loaded["spin0"][size]["standard_errors"],
            } for size, payload in loaded["spin0"].items()
        },
        "same_production_null_control": same_production_null_control(),
        "level_grid_robustness": level_grid_robustness(sizes),
        "scored": {
            name: {
                "transitions": [{k: v for k, v in row.items() if k != "deleted"}
                                for row in scored[name]["transitions"]],
                "pairwise_angles": scored[name]["pairwise_angles"],
                "held_out": scored[name]["held_out"],
                "random_direction_control": scored[name]["random_direction_control"],
            } for name in WEIGHTINGS
        },
        "not_established": [
            "any mechanism. A shape direction is a direction; naming it needs #581",
            "that a transition's cross-size covariance is exactly zero. It is measured "
            "batch by batch and reported, not assumed",
            "freedom from spin 8. The correction here removes spin 4 to first order "
            "and nothing else; N=325 and N=425 need an extrapolation rather than a "
            "mixture to get even that, which is recorded per size",
            "an exponent for the shape amplitude. Three base sizes do not determine "
            "one and this analysis deliberately does not fit it",
        ],
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    payload = assemble()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
