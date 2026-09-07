#!/usr/bin/env python3
"""P8 — N=725 is a forecast, not a sixth transition.

In the declared chart, for every surviving P6 invariant (I4, the projective
rho, is the only survivor), report:

    a_hat(290->725), se, a0, (a_hat - a0)/se,
    the 3-se interval vs the pre-registered +/-5% band,
    spin0 AND equal.

The committed forecast target (PR #614, never re-fit here):
    FORECAST_AMPLITUDE = -4.680842487165320e-04  (from #609's five-amplitude
    fit, preregistered); the scorer's band is +/-5% around it widened by
    3 sigma of the measurement.  #612's verdict stands: the spin0
    measurement supports the finite forecast (displacement -5.56 sigma of
    the FORECAST, inside the pre-registered band, which is ~10.2 sigma
    wide) -- supporting a finite forecast is NOT the model being exact.

Do not fold 725 into a three-size curvature: that would use 725 twice (as
freeze and as test), which is forbidden (GOVERNANCE section 2C).
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import affine_gauge_common as ag  # noqa: E402
import threshold_quantile_lineage as lineage  # noqa: E402
import score_wasserstein_shape_flow as flow  # noqa: E402
from p612_chart_identity import gls_fit, LEVELS  # noqa: E402

SCHEMA = "matching-one.probe-affine-gauge.p8.v1"
ISSUE = 617

#: committed, preregistered (do not re-fit)
FORECAST_AMPLITUDE = -4.680842487165320e-04
BAND_RELATIVE = 0.05          # the pre-registered +/-5% band
SIGMA_MULTIPLE = 3.0

P50_NODES = (145, 290, 725)


def _raw_blocks() -> dict[int, dict[str, Any]]:
    patched = dict(flow.SOURCES)
    patched[725] = ag.N725_SOURCE
    return {size: lineage.load_batch_histograms(ag.ROOT / patched[size])
            for size in P50_NODES}


def _transition_fit(q_a: list[float], q_b: list[float],
                    s_a: list[list[float]], s_b: list[list[float]],
                    multiplier: float, direction: list[float]
                    ) -> dict[str, Any]:
    step = math.log(multiplier)
    displacement = [(t - b) / step for b, t in zip(q_a, q_b)]
    covariance = [[(x + y) / (step * step)
                   for x, y in zip(ra, rb)]
                  for ra, rb in zip(s_a, s_b)]
    fit = gls_fit(displacement, covariance,
                  [[1.0] * LEVELS, list(q_a), list(direction)])
    return {"amplitude": fit["amplitudes"][2],
            "standard_error": fit["standard_errors"][2],
            "beta": fit["amplitudes"][1],
            "residual": fit["residual"],
            "statistic": fit["statistic"],
            "degrees_of_freedom": fit["covariance_rank"] - 3}


def main() -> dict[str, Any]:
    raw = _raw_blocks()
    committed = json.loads(
        (ag.OUT_DIR / "_reference" / "p612-n725-score-latest.json")
        .read_text())
    frozen_g = list(
        committed["weightings"]["spin0"]
        ["consensus_direction_frozen_from_five_committed_transitions"])

    weightings = {}
    for weighting in ("spin0", "equal"):
        readouts = {}
        for size, block in raw.items():
            if weighting == "spin0":
                w = lineage.spin_zero_weights(block["orientation_cos4theta"])
            else:
                w = {"first": 0.5, "second": 0.5}
            jack = lineage.jackknife_quantiles(
                block, w, levels=[0.1 * i for i in range(1, 10)])
            readouts[size] = {
                "quantiles": jack["full"],
                "covariance": lineage.jackknife_covariance(jack["full"],
                                                           jack["deleted"]),
            }
        # a_hat: the 290->725 first difference in the declared chart
        a_hat = _transition_fit(readouts[290]["quantiles"],
                                readouts[725]["quantiles"],
                                readouts[290]["covariance"],
                                readouts[725]["covariance"],
                                2.5, frozen_g)
        # a0: the 145->290 amplitude (the last committed first difference
        # NOT involving 725)
        a0 = _transition_fit(readouts[145]["quantiles"],
                             readouts[290]["quantiles"],
                             readouts[145]["covariance"],
                             readouts[290]["covariance"],
                             2.0, frozen_g)
        displacement = (a_hat["amplitude"]
                        - FORECAST_AMPLITUDE)
        low = a_hat["amplitude"] - SIGMA_MULTIPLE * a_hat["standard_error"]
        high = a_hat["amplitude"] + SIGMA_MULTIPLE * a_hat["standard_error"]
        band_low = FORECAST_AMPLITUDE * (1.0 - BAND_RELATIVE)
        band_high = FORECAST_AMPLITUDE * (1.0 + BAND_RELATIVE)
        # the forecast is NEGATIVE: (1-5%) is the smaller-magnitude end,
        # so the interval convention is [min, max] of the two
        band_low, band_high = min(band_low, band_high), max(band_low,
                                                            band_high)
        # the #612 verdict convention: the band is the pre-registered +/-5%
        # AROUND THE FORECAST, and the 3-se measurement interval must lie
        # inside it for "supports this finite forecast"
        inside = low >= band_low and high <= band_high
        # I4 (the only P6 survivor): the projective rho of the new
        # transition, denominator predeclared as the width coefficient
        rho = a_hat["beta"] / a_hat["amplitude"] if a_hat["amplitude"] else None
        weightings[weighting] = {
            "a_hat_290_725": a_hat["amplitude"],
            "a_hat_standard_error": a_hat["standard_error"],
            "a0_145_290": a0["amplitude"],
            "a_hat_minus_forecast": displacement,
            "displacement_in_forecast_sigmas": (
                displacement / abs(FORECAST_AMPLITUDE)),
            "displacement_in_measurement_sigmas": (
                displacement / a_hat["standard_error"]),
            "measurement_interval_3se": [low, high],
            "preregistered_band_5pct": [band_low, band_high],
            "inside_band": inside,
            "i4_projective_rho": rho,
            "residual_statistic": a_hat["statistic"],
            "residual_degrees_of_freedom": a_hat["degrees_of_freedom"],
        }
    spin0 = weightings["spin0"]
    equal = weightings["equal"]
    both_inside = spin0["inside_band"] and equal["inside_band"]
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "forecast_target_preregistered": FORECAST_AMPLITUDE,
        "weightings": weightings,
        "both_weightings_inside_band": both_inside,
        "verdict": (
            "BOTH weightings sit inside the pre-registered +/-5%% band "
            "around the frozen forecast at 3 se: the finite forecast is "
            "supported.  This does NOT make the model exact (GOVERNANCE "
            "2C); the displacement is %.1f forecast-sigmas (spin0) and the "
            "residual statistic is %.0f on %d df, not zero."
            % (spin0["displacement_in_forecast_sigmas"],
               spin0["residual_statistic"],
               spin0["residual_degrees_of_freedom"])
            if both_inside else
            "at least one weighting misses the pre-registered band: the "
            "forecast is NOT supported under the declared chart."),
    }


if __name__ == "__main__":
    result = main()
    path = ag.dump("p8-n725", result)
    print(json.dumps({
        "both_inside": result["both_weightings_inside_band"],
        "verdict": result["verdict"],
        "output": str(path)}, indent=2))
