#!/usr/bin/env python3
"""P4 — weighting as a covector on the two-orientation fibre.

spin0 and equal are two points in the affine line of orientation weightings
w(t) = (1-t) w_spin0 + t w_equal, t in [0,1].  For each t this script
rebuilds the jackknife quantiles under w(t) (straight from the per-batch,
per-orientation histograms -- the raw orientations ARE stored for every
committed block, and for N=725), and measures:

- the first-difference amplitude of g on every transition;
- the N=725 forecast residual (290->725 amplitude minus the frozen forecast);
- the transported curvature residual.

Deliverable: the 12% / 25-sigma spin0->equal shift as a DIRECTIONAL
DERIVATIVE along t (and one more direction: the segment is the only
physically available one, since the fibre is one-dimensional in the two
stored orientations -- an "orthogonal complement" direction does not exist
as an independent weighting; the honest orthogonal object is recorded as
absent with a reason).
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
from p612_chart_identity import LEVELS  # noqa: E402

SCHEMA = "matching-one.probe-affine-gauge.p4.v1"
ISSUE = 617

TS = [round(i / 10, 1) for i in range(11)]     # 0.0 ... 1.0
COMMITTED_FORECAST = -4.680842487165320e-04


def load_raw_sizes(sources: dict[int, str]) -> dict[int, dict[str, Any]]:
    """Per-batch histograms for each size, from the raw committed blocks."""
    out = {}
    for size, rel in sources.items():
        out[size] = lineage.load_batch_histograms(ag.ROOT / rel)
    return out


def jackknife_under(weights_fn, raw: dict[str, Any]) -> dict[str, Any]:
    weights = weights_fn(raw["orientation_cos4theta"])
    jack = lineage.jackknife_quantiles(raw, weights)
    covariance = lineage.jackknife_covariance(jack["full"], jack["deleted"])
    return {
        "weights": weights,
        "quantiles": jack["full"],
        "deleted": jack["deleted"],
        "covariance": covariance,
        "batches": jack["batches"],
    }


def main() -> dict[str, Any]:
    raw = load_raw_sizes({**flow.SOURCES, ag.N725: ag.N725_SOURCE})
    g = ag.frozen_direction("spin0")
    import p582_amplitude_law as law

    def w_at(t: float, cos4: dict[str, float]) -> dict[str, float]:
        spin0 = lineage.spin_zero_weights(cos4)
        return {"first": (1 - t) * spin0["first"] + t * 0.5,
                "second": (1 - t) * spin0["second"] + t * 0.5}

    rows = []
    cache: dict[tuple[float, int], dict[str, Any]] = {}
    for t in TS:
        loaded: dict[int, dict[str, Any]] = {}
        for size, payload in raw.items():
            loaded[size] = jackknife_under(lambda cos4, t=t: w_at(t, cos4),
                                           payload)
        # first-difference amplitudes on every committed transition
        amps = {}
        for name, sizes in flow.LINEAGES.items():
            for base, target in zip(sizes, sizes[1:]):
                geo = ag.transition_geometry(loaded, base, target, g)
                amps[f"{base}->{target}"] = geo["amplitude"]
        # N=725 forecast residual: the 290->725 amplitude in the same chart
        geo725 = ag.transition_geometry(loaded, 290, ag.N725, g)
        forecast_resid = geo725["amplitude"] - COMMITTED_FORECAST
        # transported curvature residual at lambda=1, per lineage
        curv = {}
        for name, sizes in flow.LINEAGES.items():
            if len(sizes) < 3:
                continue
            entry = ag.lineage_identity_at(loaded, sizes, g, 1.0)
            curv[name] = {
                "a_curv": entry["curvature_measured"],
                "ratio_naive": entry["identity"]["ratio_measured_over_naive"],
            }
        rows.append({
            "t": t,
            "weights_example": loaded[65]["weights"],
            "first_amplitudes": amps,
            "n725_amplitude": geo725["amplitude"],
            "n725_forecast_residual": forecast_resid,
            "curvature": curv,
        })

    # directional derivative at t=0 (spin0): d(a_290->725)/dt
    derivatives = {}
    for key in rows[0]["first_amplitudes"]:
        derivatives[key] = (rows[-1]["first_amplitudes"][key]
                            - rows[0]["first_amplitudes"][key])
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "t_grid": TS,
        "rows": rows,
        "directional_derivatives_spin0_to_equal": derivatives,
        "orthogonal_complement_note": (
            "the weighting fibre is ONE-dimensional: the two stored "
            "orientations span it, w0+w1=1 fixes the affine gauge, and there "
            "is no independent 'orthogonal weighting' to test.  The honest "
            "record is that the orthogonal direction does not exist in the "
            "stored data."),
        "summary": {
            "amplitude_shift_290_725_percent": 100.0 * abs(
                rows[-1]["n725_amplitude"] - rows[0]["n725_amplitude"])
                / abs(rows[0]["n725_amplitude"]),
            "n725_residual_at_spin0": rows[0]["n725_forecast_residual"],
            "n725_residual_at_equal": rows[-1]["n725_forecast_residual"],
        },
    }


if __name__ == "__main__":
    result = main()
    path = ag.dump("p4-weighting", result)
    print(json.dumps({"summary": result["summary"],
                      "derivatives": result["directional_derivatives_spin0_to_equal"],
                      "output": str(path)}, indent=2))
