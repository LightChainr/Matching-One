#!/usr/bin/env python3
"""#612 steps 3-5: score the one N725 block against the declared forecast.

This is the *scorer* that `notes/p612-preregistration-20260907.md` froze before
the N725 histograms were read.  Everything it fixes is fixed here and not
refit:

* the direction `g` — rebuilt from the five **committed** transitions only;
  725 is never allowed into the consensus;
* the nine deciles and the `spin0` weighting (`equal` is a declared
  sensitivity from the same block, not a second attempt);
* the estimator for the 290 -> 725 first difference — basis `[1, Q_290, g]`,
  observation `(Q_725 - Q_290)/log 2.5`, covariance
  `(S_290 + S_725)/log(2.5)^2`, and the reported number is the **third
  coefficient** of that fit, not an amplitude read off a differently
  normalized curve;
* cross-size covariance: none.  The three productions use distinct seeds, so
  `S_290 + S_725` is the declared convention, and this module *measures* the
  batch-level correlation rather than assuming it away.

The targets are #612's own and are already public, so this is a **reanalysis**,
not a blind forecast test; the module says so in its output rather than
implying otherwise.

The secondary p50 curvature is compared against the **conditional chart
calculation**, not against `6.70567e-4`, because the chart result says a raw
curvature near `1.03e-3` is already what one transported exponent predicts.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

try:  # pragma: no cover - import shape depends on how the script is invoked
    from scripts import p582_amplitude_law as law
    from scripts import score_wasserstein_shape_flow as flow
    from scripts.p612_chart_identity import gls_fit, lineage_identity, transition_geometry
except ModuleNotFoundError:  # pragma: no cover
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import p582_amplitude_law as law
    import score_wasserstein_shape_flow as flow
    from p612_chart_identity import gls_fit, lineage_identity, transition_geometry

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "p612-n725-score" / "latest.json"
SCHEMA = "matching-one.p612-n725-score.v1"
ISSUE = 612

N725 = 725
N725_SOURCE = "results/server-20260907/P612-n725-fullcurve/raw/n725_100m.hist.csv"

PRIMARY_WEIGHTING = "spin0"
CONTROL_WEIGHTING = "equal"

#: #612's declared forecast and its 5% tolerance band.  Public before this run,
#: which is why the run is a reanalysis (see the preregistration note).
FORECAST_AMPLITUDE = -0.0004680842487166532
TOLERANCE_BAND = (-0.0004914884611524859, -0.0004446800362808205)
SIGMA_MULTIPLE = 3.0

#: The p50 second-difference weights #609 pre-registered.
PRE_REGISTERED_P50_WEIGHTS = [1.79279, -3.14899, 1.35620]
#: The two curvature numbers #609 said the run would separate.  After the chart
#: result they are not separated, so the second one is reported, not tested.
ONE_EXPONENT_CURVATURE_TARGET = 6.70567e-04
FACTOR_1547_CURVATURE_TARGET = 1.03745e-03


def _load_with_725():
    """Load every committed size plus 725, leaving ``flow.LINEAGES`` untouched.

    ``flow.LINEAGES`` must stay as committed so that ``flow.transitions()``
    returns the five transitions that froze ``g``; only ``SOURCES`` grows.
    """
    patched = dict(flow.SOURCES)
    patched[N725] = N725_SOURCE
    original = flow.SOURCES
    flow.SOURCES = patched
    try:
        sizes = sorted({size for sizes in flow.LINEAGES.values() for size in sizes}
                       | {N725})
        return flow.load_sizes(sizes)
    finally:
        flow.SOURCES = original


def decision(measured: float, standard_error: float) -> dict[str, Any]:
    """#612's declared three-standard-error rule, applied literally."""
    low = measured - SIGMA_MULTIPLE * standard_error
    high = measured + SIGMA_MULTIPLE * standard_error
    band_low, band_high = TOLERANCE_BAND
    if high < band_low or low > band_high:
        outcome = "stops_this_forecast"
    elif low >= band_low and high <= band_high:
        outcome = "supports_this_finite_forecast"
    else:
        outcome = "unresolved"
    return {
        "measured": measured,
        "standard_error": standard_error,
        "sigma_multiple": SIGMA_MULTIPLE,
        "measurement_interval": [low, high],
        "tolerance_band": [band_low, band_high],
        "displacement_from_the_forecast": measured - FORECAST_AMPLITUDE,
        "relative_displacement": (measured - FORECAST_AMPLITUDE)
                                 / abs(FORECAST_AMPLITUDE),
        "displacement_in_standard_errors": (measured - FORECAST_AMPLITUDE)
                                           / standard_error,
        "outcome": outcome,
    }


def score_weighting(both: Mapping[str, Mapping[int, Mapping[str, Any]]],
                    name: str) -> dict[str, Any]:
    loaded = both[name]
    rows = law.transition_rows(loaded)
    reproduction = law.reproduction_control(rows)

    # g is frozen here and never sees 725.
    direction = law.consensus_direction(rows)
    points = law.amplitude_points(rows, loaded, direction)
    fit = law.fit_exponent(points)

    # ---- primary: the 290 -> 725 first difference, in the declared chart ----
    step = math.log(N725) - math.log(290)
    sbar = (math.log(290) + math.log(N725)) / 2.0
    new_point = transition_geometry(loaded, 290, N725, direction)
    predicted_amplitude = law.first_difference_image(fit["scale_amplitude"],
                                                     fit["omega"], sbar, step)
    verdict = decision(new_point["amplitude"], new_point["standard_error"])

    # ---- the label crossing #609 said 725 would break ----
    orientation = loaded[N725]
    five_adic = law.five_adic_valuation(N725)
    interpolation = orientation["spin4_correction_is_an_interpolation"]
    crossing = {
        "size": N725,
        "five_adic_valuation": five_adic,
        "spin0_is_an_interpolation": interpolation,
        "orientation_cos4theta": orientation["orientation_cos4theta"],
        "spin0_weights": orientation["weights"],
        "net_cos4theta": orientation["net_cos4theta"],
        "valuation_two_and_interpolating": five_adic == 2 and interpolation,
        "breaks_the_committed_identification": five_adic == 2 and interpolation,
        "note": ("on the committed sizes v_5=2, 'no interpolating spin-0 "
                 "combination' and 'multiplier 2.5' were the same partition; "
                 "725 has v_5=2 and an interpolating combination, so they differ"),
    }

    # ---- secondary: the p50 curvature now that there are three rungs ----
    nodes = [math.log(145), math.log(290), math.log(N725)]
    weights = law.second_difference_weights(nodes)
    h0 = nodes[1] - nodes[0]
    h1 = nodes[2] - nodes[1]
    factor = 2.0 / (h0 + h1)
    observation = [math.fsum(w * loaded[size]["quantiles"][level]
                             for w, size in zip(weights, (145, 290, N725)))
                   for level in range(law.LEVELS)]
    covariance = [[math.fsum(w * w * loaded[size]["covariance"][i][j]
                             for w, size in zip(weights, (145, 290, N725)))
                   for j in range(law.LEVELS)] for i in range(law.LEVELS)]
    curvature_fit = gls_fit(observation, covariance,
                            [[1.0] * law.LEVELS,
                             list(loaded[290]["quantiles"]), list(direction)])
    ell_c = curvature_fit["covectors"][2]

    lower = transition_geometry(loaded, 145, 290, direction)
    upper = new_point
    k = 1.0 + h0 * lower["beta"]
    a0, a1 = lower["amplitude"], upper["amplitude"]
    transported = a1 - a0 / k
    residual_transport = (math.fsum(c * v for c, v in zip(ell_c, upper["residual"]))
                          - math.fsum(c * v for c, v in zip(ell_c, lower["residual"])) / k)
    rebuilt = factor * (transported + residual_transport)
    measured_curvature = curvature_fit["amplitudes"][2]

    model_lower = law.first_difference_image(fit["scale_amplitude"], fit["omega"],
                                             (nodes[0] + nodes[1]) / 2.0, h0)
    model_upper = law.first_difference_image(fit["scale_amplitude"], fit["omega"],
                                             (nodes[1] + nodes[2]) / 2.0, h1)
    naive_prediction = factor * (model_upper - model_lower)
    chart_corrected = factor * (model_upper - model_lower / k)

    curvature = {
        "sizes": [145, 290, N725],
        "second_difference_weights": weights,
        "pre_registered_p50_weights": PRE_REGISTERED_P50_WEIGHTS,
        "weights_agree_with_pre_registration":
            max(abs(a - b) for a, b in zip(weights, PRE_REGISTERED_P50_WEIGHTS)) < 5e-5,
        "chart_factor": factor,
        "measured": measured_curvature,
        "standard_error": curvature_fit["standard_errors"][2],
        "residual_statistic": curvature_fit["statistic"],
        "residual_degrees_of_freedom": curvature_fit["covariance_rank"] - 3,
        "k_145_to_290": k,
        "identity": {
            "a1_minus_a0_over_k": transported,
            "ell_C_of_r1": math.fsum(c * v for c, v in zip(ell_c, upper["residual"])),
            "ell_C_of_r0_over_k": math.fsum(c * v for c, v in zip(ell_c, lower["residual"])) / k,
            "residual_transport": residual_transport,
            "rebuilt": rebuilt,
            "absolute_error": abs(rebuilt - measured_curvature),
        },
        "predictions": {
            "naive_one_exponent": naive_prediction,
            "chart_corrected": chart_corrected,
            "one_exponent_target_from_609": ONE_EXPONENT_CURVATURE_TARGET,
            "factor_1547_target_from_609": FACTOR_1547_CURVATURE_TARGET,
            "chart_corrected_minus_factor_1547_target":
                chart_corrected - FACTOR_1547_CURVATURE_TARGET,
            "relative_separation_of_the_two_609_targets":
                (FACTOR_1547_CURVATURE_TARGET - ONE_EXPONENT_CURVATURE_TARGET)
                / ONE_EXPONENT_CURVATURE_TARGET,
            "the_two_targets_are_separated_by_the_chart":
                abs(FACTOR_1547_CURVATURE_TARGET - chart_corrected)
                < 0.1 * abs(FACTOR_1547_CURVATURE_TARGET
                            - ONE_EXPONENT_CURVATURE_TARGET),
        },
        "measured_over_chart_corrected": measured_curvature / chart_corrected,
        "measured_over_naive": measured_curvature / naive_prediction,
    }

    coupling = flow.cross_size_coupling(loaded[290], loaded[N725])

    return {
        "weighting": name,
        "reproduction_control": reproduction,
        "consensus_direction_frozen_from_five_committed_transitions": direction,
        "transitions_used_to_freeze_g": [row["label"] for row in rows],
        "exponent_fit": {key: fit[key] for key in
                         ("omega", "scale_amplitude", "statistic",
                          "degrees_of_freedom", "search_window")},
        "amplitude_points": [{key: point[key] for key in
                              ("label", "sbar", "h", "amplitude", "standard_error")}
                             for point in points],
        "new_transition": {
            "label": "290->725",
            "basis": ["1", "Q_290", "g_frozen"],
            "observation": "(Q_725 - Q_290)/log(2.5)",
            "covariance_convention": "(S_290 + S_725)/log(2.5)^2, "
                                     "delete-one jackknife, no cross-size term",
            "amplitude": new_point["amplitude"],
            "standard_error": new_point["standard_error"],
            "significance": new_point["amplitude"] / new_point["standard_error"],
            "beta_Q290": new_point["beta"],
            "residual_statistic": new_point["residual_statistic"],
            "residual_degrees_of_freedom": new_point["residual_degrees_of_freedom"],
            "predicted_from_the_five_amplitude_fit": predicted_amplitude,
            "forecast_target": FORECAST_AMPLITUDE,
            "forecast_reproduced_by_the_frozen_fit":
                abs(predicted_amplitude - FORECAST_AMPLITUDE) < 1e-12,
        },
        "forecast_verdict": verdict,
        "label_crossing": crossing,
        "p50_curvature": curvature,
        "cross_size_coupling_290_725": coupling,
    }


def assemble() -> dict[str, Any]:
    both = _load_with_725()
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "design": "reanalysis: every target in #612 was public before this run",
        "preregistration": "notes/p612-preregistration-20260907.md",
        "production": {
            "size": N725,
            "orientations": [[26, 7], [23, 14]],
            "batches": 100,
            "samples_per_batch": 1_000_000,
            "paired_configurations": 100_000_000,
            "samples_per_orientation": 100_000_000,
            "sample_count_correction": ("#609 wrote '100 batches x 100M', which "
                                        "is 100x this; n290 metadata says "
                                        "samples_per_pair=100000000, and we match it"),
        },
        "weightings": {name: score_weighting(both, name)
                       for name in (PRIMARY_WEIGHTING, CONTROL_WEIGHTING)},
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    report = assemble()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"wrote {args.output}\n")
    for name, block in report["weightings"].items():
        new = block["new_transition"]
        verdict = block["forecast_verdict"]
        curv = block["p50_curvature"]
        print(f"[{name}]")
        print(f"  a_hat(290->725) = {new['amplitude']:+.10e}  "
              f"se = {new['standard_error']:.4e}  "
              f"({new['significance']:.1f} sigma)")
        print(f"  a0 target       = {FORECAST_AMPLITUDE:+.10e}   "
              f"frozen fit predicts {new['predicted_from_the_five_amplitude_fit']:+.10e}")
        print(f"  a_hat - a0      = {verdict['displacement_from_the_forecast']:+.4e}  "
              f"= {verdict['displacement_in_standard_errors']:+.2f} se")
        print(f"  3se interval    = [{verdict['measurement_interval'][0]:+.6e}, "
              f"{verdict['measurement_interval'][1]:+.6e}]")
        print(f"  tolerance band  = [{TOLERANCE_BAND[0]:+.6e}, {TOLERANCE_BAND[1]:+.6e}]")
        print(f"  OUTCOME         = {verdict['outcome']}")
        print(f"  curvature: measured {curv['measured']:+.6e}  "
              f"chart-corrected {curv['predictions']['chart_corrected']:+.6e}  "
              f"ratio {curv['measured_over_chart_corrected']:.5f}")
        print(f"             identity |err| = {curv['identity']['absolute_error']:.3e}"
              f"   coupling droppable = "
              f"{block['cross_size_coupling_290_725']['cross_term_droppable']}")
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
