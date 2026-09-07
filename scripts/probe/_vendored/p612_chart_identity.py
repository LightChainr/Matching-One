#!/usr/bin/env python3
"""#612 step 1: is the 55% curvature excess an affine change of chart?

The curvature check in ``scripts/p582_amplitude_law.py`` fits each three-size
lineage's second divided difference against the basis ``[1, Q_middle, g]``.
The five first differences that the exponent was fitted to are measured against
``[1, Q_base, g]``.  Those are two different affine charts.  ``g`` being frozen
does not make the *extracted amplitude* a fixed covector applied to ``Q``:
the covector is the row of a weighted Gram inverse, and the Gram depends on
which ``Q`` sits in the basis.

Write the first step as

    v0 = alpha0 1 + beta0 Q0 + a0 g + r0,        k = 1 + h0 beta0,

so that, exactly,

    Q1 = k Q0 + h0 alpha0 1 + h0 a0 g + h0 r0,

and therefore, still exactly,

    v0 = (alpha0/k) 1 + (beta0/k) Q1 + (a0/k) g + r0/k.

Let ``ell_C`` be the curvature estimator's own covector: the row that extracts
the ``g`` coefficient from a GLS fit against ``[1, Q1, g]`` with the curvature
covariance.  It satisfies ``ell_C(1) = ell_C(Q1) = 0`` and ``ell_C(g) = 1``.
The curvature observation is ``C * (v12 - v01)`` with ``C = 2/(h0+h1)``, so

    a_curvature = C [a1 - a0/k] + C [ell_C(r1) - ell_C(r0)/k].

The first bracket is the chart-transported comparison; setting ``k = 1``
recovers the one the exponent was checked against.  The second bracket is the
residual transport, which is small but not zero and is reported rather than
assumed away.

This module recomputes every term from the eight committed productions, using
the production mpmath path and the published ``g``, and closes the identity.
It deliberately does **not** refit anything: if the identity fails, the
discrepant term is returned and scoring stops.

Nothing is sampled.  Every histogram is a committed production artifact.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from mpmath import mp

try:  # pragma: no cover - import shape depends on how the script is invoked
    from scripts import p582_amplitude_law as law
    from scripts import score_wasserstein_shape_flow as flow
    from scripts.projective_inference import (spectral_pseudo_inverse,
                                              subspace_residual)
except ModuleNotFoundError:  # pragma: no cover
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import p582_amplitude_law as law
    import score_wasserstein_shape_flow as flow
    from projective_inference import spectral_pseudo_inverse, subspace_residual

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "p612-chart-identity" / "latest.json"
SCHEMA = "matching-one.p612-chart-identity.v1"
ISSUE = 612

LEVELS = 9
PRIMARY_WEIGHTING = "spin0"
CONTROL_WEIGHTING = "equal"
DPS = 60

#: #612 step 1 asks for the identity to close at this level.  Anything larger is
#: a discrepancy to be reported, not a rounding artefact to be averaged away.
IDENTITY_TOLERANCE = 1e-9


def gls_fit(observation: Sequence[float],
            covariance: Sequence[Sequence[float]],
            basis: Sequence[Sequence[float]],
            dps: int = DPS) -> dict[str, Any]:
    """GLS fit, its residual, and the covector of every coefficient.

    ``covectors[i]`` is the row vector ``ell_i`` with ``ell_i(b_j) = delta_ij``:
    ``ell_i(y)`` is the coefficient of ``basis[i]`` in the fit of ``y``.  The
    curvature covector used below is ``covectors[2]``.
    """
    with mp.workdps(dps):
        y = mp.matrix([mp.mpf(str(v)) for v in observation])
        size = y.rows
        columns = [mp.matrix([mp.mpf(str(v)) for v in direction]) for direction in basis]
        for column in columns:
            if column.rows != size:
                raise ValueError("a basis direction has the wrong length")
        pinv, rank, condition, _ = spectral_pseudo_inverse(covariance)
        width = len(columns)
        weighted = [pinv * column for column in columns]
        gram = mp.zeros(width, width)
        for i in range(width):
            for j in range(width):
                gram[i, j] = (columns[i].T * weighted[j])[0]
        inverse = gram ** -1
        rhs = mp.matrix([(columns[i].T * (pinv * y))[0] for i in range(width)])
        amplitudes = inverse * rhs
        fitted = sum((amplitudes[i] * columns[i] for i in range(width)),
                     mp.zeros(size, 1))
        residual = y - fitted
        # ell_i = S^+ V (V^T S^+ V)^{-1} e_i, as a column of that matrix product.
        covectors = []
        for i in range(width):
            column = sum((inverse[j, i] * weighted[j] for j in range(width)),
                         mp.zeros(size, 1))
            covectors.append([float(v) for v in column])
        statistic = (residual.T * (pinv * residual))[0]
        return {
            "amplitudes": [float(amplitudes[i]) for i in range(width)],
            "residual": [float(residual[i]) for i in range(size)],
            "covectors": covectors,
            "standard_errors": [float(mp.sqrt(inverse[i, i])) for i in range(width)],
            "statistic": float(statistic),
            "covariance_rank": rank,
            "covariance_condition_number": float(condition),
        }


def _apply(covector: Sequence[float], vector: Sequence[float]) -> float:
    return math.fsum(left * right for left, right in zip(covector, vector))


def transition_geometry(loaded: Mapping[int, Mapping[str, Any]],
                        base: int, target: int,
                        direction: Sequence[float]) -> dict[str, Any]:
    """Everything the identity needs about one first difference."""
    step = math.log(target) - math.log(base)
    displacement = [ (t - b) / step
                     for b, t in zip(loaded[base]["quantiles"],
                                     loaded[target]["quantiles"]) ]
    covariance = [[(a + b) / (step * step)
                   for a, b in zip(row_a, row_b)]
                  for row_a, row_b in zip(loaded[base]["covariance"],
                                          loaded[target]["covariance"])]
    fit = gls_fit(displacement, covariance,
                  [[1.0] * LEVELS, list(loaded[base]["quantiles"]),
                   list(direction)])
    return {
        "base": base,
        "target": target,
        "step": step,
        "alpha": fit["amplitudes"][0],
        "beta": fit["amplitudes"][1],
        "amplitude": fit["amplitudes"][2],
        "standard_error": fit["standard_errors"][2],
        "residual": fit["residual"],
        "residual_statistic": fit["statistic"],
        "residual_degrees_of_freedom": fit["covariance_rank"] - 3,
    }


def lineage_identity(loaded: Mapping[int, Mapping[str, Any]],
                     sizes: Sequence[int],
                     direction: Sequence[float],
                     scale_amplitude: float, omega: float) -> dict[str, Any]:
    """Close the identity on one three-size lineage, and decompose the excess."""
    first, middle, last = sizes
    h0 = math.log(middle) - math.log(first)
    h1 = math.log(last) - math.log(middle)
    factor = 2.0 / (h0 + h1)

    lower = transition_geometry(loaded, first, middle, direction)
    upper = transition_geometry(loaded, middle, last, direction)

    # The curvature observation, its covariance, and its own covector.
    weights = law.second_difference_weights([math.log(s) for s in sizes])
    observation = [math.fsum(weight * loaded[size]["quantiles"][level]
                             for weight, size in zip(weights, sizes))
                   for level in range(LEVELS)]
    covariance = [[math.fsum(weight * weight * loaded[size]["covariance"][i][j]
                             for weight, size in zip(weights, sizes))
                   for j in range(LEVELS)] for i in range(LEVELS)]
    curvature_fit = gls_fit(observation, covariance,
                            [[1.0] * LEVELS, list(loaded[middle]["quantiles"]),
                             list(direction)])
    ell_c = curvature_fit["covectors"][2]

    a0, a1 = lower["amplitude"], upper["amplitude"]
    k = 1.0 + h0 * lower["beta"]

    transported = a1 - a0 / k
    residual_transport = _apply(ell_c, upper["residual"]) \
        - _apply(ell_c, lower["residual"]) / k
    rebuilt = factor * (transported + residual_transport)
    measured = curvature_fit["amplitudes"][2]
    absolute_error = abs(rebuilt - measured)
    relative_error = absolute_error / abs(measured) if measured else float("nan")

    # The one-exponent model's own first amplitudes, and the two predictions.
    nodes = [math.log(size) for size in sizes]
    model_lower = law.first_difference_image(scale_amplitude, omega,
                                             (nodes[0] + nodes[1]) / 2.0, h0)
    model_upper = law.first_difference_image(scale_amplitude, omega,
                                             (nodes[1] + nodes[2]) / 2.0, h1)
    predicted = law.second_difference_image(scale_amplitude, omega, nodes)
    predicted_chart_corrected = factor * (model_upper - model_lower / k)
    chart_term = predicted_chart_corrected - predicted
    excess = measured - predicted

    return {
        "lineage_sizes": list(sizes),
        "h0": h0,
        "h1": h1,
        "chart_factor": factor,
        "second_difference_weights": weights,
        "lower_transition": {key: lower[key] for key in
                             ("base", "target", "step", "alpha", "beta",
                              "amplitude", "standard_error",
                              "residual_statistic",
                              "residual_degrees_of_freedom")},
        "upper_transition": {key: upper[key] for key in
                             ("base", "target", "step", "alpha", "beta",
                              "amplitude", "standard_error",
                              "residual_statistic",
                              "residual_degrees_of_freedom")},
        "k": k,
        "k_is_the_width_contraction": True,
        "curvature_measured": measured,
        "curvature_standard_error": curvature_fit["standard_errors"][2],
        "curvature_residual_statistic": curvature_fit["statistic"],
        "curvature_residual_degrees_of_freedom": curvature_fit["covariance_rank"] - 3,
        "identity": {
            "a1_minus_a0_over_k": transported,
            "ell_C_of_r1": _apply(ell_c, upper["residual"]),
            "ell_C_of_r0_over_k": _apply(ell_c, lower["residual"]) / k,
            "residual_transport": residual_transport,
            "rebuilt_curvature": rebuilt,
            "measured_curvature": measured,
            "absolute_error": absolute_error,
            "relative_error": relative_error,
            "closes_within_tolerance": absolute_error < IDENTITY_TOLERANCE,
            "tolerance": IDENTITY_TOLERANCE,
        },
        "decomposition": {
            "naive_C_times_a1_minus_a0": factor * (a1 - a0),
            "chart_transported_C_times_a1_minus_a0_over_k": factor * transported,
            "residual_transport_term": factor * residual_transport,
            "model_first_amplitude_lower": model_lower,
            "model_first_amplitude_upper": model_upper,
            "one_exponent_prediction": predicted,
            "chart_corrected_prediction": predicted_chart_corrected,
            "chart_term": chart_term,
            "excess_over_the_one_exponent_prediction": excess,
            "fraction_of_excess_from_the_chart_term":
                chart_term / excess if excess else None,
            "ratio_measured_over_one_exponent": measured / predicted,
            "ratio_measured_over_chart_corrected":
                measured / predicted_chart_corrected,
        },
    }


def assemble() -> dict[str, Any]:
    """The identity on both weightings, on every three-size lineage."""
    sizes = sorted({size for sizes in flow.LINEAGES.values() for size in sizes})
    both = flow.load_sizes(sizes)
    out: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": ISSUE,
        "weightings": {},
        "identity_tolerance": IDENTITY_TOLERANCE,
    }
    for name in (PRIMARY_WEIGHTING, CONTROL_WEIGHTING):
        loaded = both[name]
        rows = law.transition_rows(loaded)
        reproduction = law.reproduction_control(rows)
        direction = law.consensus_direction(rows)
        points = law.amplitude_points(rows, loaded, direction)
        fit = law.fit_exponent(points)
        per_lineage = []
        for lineage_name, lineage_sizes in flow.LINEAGES.items():
            if len(lineage_sizes) < 3:
                per_lineage.append({"lineage": lineage_name,
                                    "sizes": list(lineage_sizes),
                                    "usable": False,
                                    "reason": "a second difference needs three sizes"})
                continue
            entry = lineage_identity(loaded, lineage_sizes, direction,
                                     fit["scale_amplitude"], fit["omega"])
            entry["lineage"] = lineage_name
            entry["usable"] = True
            per_lineage.append(entry)
        usable = [entry for entry in per_lineage if entry.get("usable")]
        out["weightings"][name] = {
            "reproduction_control": reproduction,
            "consensus_direction": direction,
            "exponent_fit": {key: fit[key] for key in
                             ("omega", "scale_amplitude", "statistic",
                              "degrees_of_freedom", "search_window")},
            "amplitude_points": [{key: point[key] for key in
                                  ("label", "sbar", "h", "amplitude",
                                   "standard_error", "significance")}
                                 for point in points],
            "per_lineage": per_lineage,
            "every_identity_closes": all(entry["identity"]["closes_within_tolerance"]
                                         for entry in usable),
            "largest_absolute_identity_error": max(
                entry["identity"]["absolute_error"] for entry in usable),
            "chart_fraction_of_excess": [entry["decomposition"][
                "fraction_of_excess_from_the_chart_term"] for entry in usable],
        }
    return out


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    report = assemble()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"wrote {args.output}")
    for name, block in report["weightings"].items():
        print(f"\n[{name}]  omega={block['exponent_fit']['omega']:.6f}  "
              f"lambda={block['exponent_fit']['scale_amplitude']:.6e}")
        print(f"  identities close: {block['every_identity_closes']}  "
              f"largest |error| = {block['largest_absolute_identity_error']:.3e}")
        for entry in block["per_lineage"]:
            if not entry.get("usable"):
                print(f"  {entry['lineage']}: {entry['reason']}")
                continue
            identity = entry["identity"]
            decomp = entry["decomposition"]
            print(f"  {entry['lineage']}: k={entry['k']:.7f}  "
                  f"M={identity['measured_curvature']:+.6e}  "
                  f"rebuilt={identity['rebuilt_curvature']:+.6e}  "
                  f"|err|={identity['absolute_error']:.3e}")
            print(f"      P ={decomp['one_exponent_prediction']:+.6e}  "
                  f"P'={decomp['chart_corrected_prediction']:+.6e}  "
                  f"chart/E = {decomp['fraction_of_excess_from_the_chart_term']:.4f}  "
                  f"M/P={decomp['ratio_measured_over_one_exponent']:.5f}  "
                  f"M/P'={decomp['ratio_measured_over_chart_corrected']:.5f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
