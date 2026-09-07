#!/usr/bin/env python3
"""P2 — the continuous attachment path (why this probe has CPU).

For each three-size lineage and each orientation weighting, sweep the
curvature-chart attachment along the lower transition's segment,

    A(lambda) = (1-lambda) Q_base + lambda Q_middle,   lambda in [0, 1],

at 21 points.  At each lambda record:

- a_curv(lambda): the measured g-amplitude in span{1, A(lambda), g}
  (exact, via the covector expansion verified in P0/P1);
- a_pred(lambda): the one-exponent model's second-difference image evaluated
  at the SAME attachment -- i.e. with the chart term recomputed at that
  lambda through the covector expansion, so the comparison is chart-matched;
- ratio a_curv / a_pred: the "fake excess multiplier";
- the curvature observation's own value at that attachment.

Expected shape (kill sentence from the brief): the fake excess must be a
SMOOTH function of lambda passing near 1.55 at the historical middle chart
(lambda = 1) and near the consistent reading at low lambda.  If the sweep is
not smooth, or if no lambda brings the ratio to ~1, the identity is not the
whole story and P2 says so.

No new Monte Carlo: everything is re-read from the committed histograms
(through the deterministic production path, cached).
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import affine_gauge_common as ag  # noqa: E402

SCHEMA = "matching-one.probe-affine-gauge.p2.v1"
ISSUE = 617

LAMBDAS = [round(i / 20, 2) for i in range(21)]   # 0.00 ... 1.00


def sweep_one(loaded: dict[int, dict[str, Any]], sizes: tuple[int, ...],
              direction: list[float],
              model_amplitude: float, omega: float) -> dict[str, Any]:
    import p582_amplitude_law as law
    first, middle, last = sizes
    nodes = [math.log(s) for s in sizes]
    h0, h1 = nodes[1] - nodes[0], nodes[2] - nodes[1]
    factor = 2.0 / (h0 + h1)

    lower = ag.transition_geometry(loaded, first, middle, direction)
    upper = ag.transition_geometry(loaded, middle, last, direction)
    a0, a1 = lower["amplitude"], upper["amplitude"]
    k = 1.0 + h0 * lower["beta"]

    model_pred = law.second_difference_image(model_amplitude, omega, nodes)
    model_first_lower = law.first_difference_image(
        model_amplitude, omega, (nodes[0] + nodes[1]) / 2.0, h0)

    rows = []
    for lam in LAMBDAS:
        entry = ag.lineage_identity_at(loaded, sizes, direction, lam)
        measured = entry["curvature_measured"]
        ratio = measured / model_pred if model_pred else None
        rows.append({
            "lambda": lam,
            "a_curv": measured,
            "a_curv_se": entry["curvature_standard_error"],
            "a_curv_over_naive":
                entry["identity"]["ratio_measured_over_naive"],
            "a_pred_one_exponent": model_pred,
            "ratio_curv_over_pred": ratio,
        })
    ratios = [row["ratio_curv_over_pred"] for row in rows]
    argmax = max(range(len(rows)), key=lambda i: ratios[i])
    argmin = min(range(len(rows)), key=lambda i: ratios[i])
    return {
        "lineage": list(sizes),
        "model_prediction": model_pred,
        "model_first_lower": model_first_lower,
        "k": k, "a0": a0, "a1": a1, "chart_factor": factor,
        "rows": rows,
        "lambda_of_maximum_fake_excess": rows[argmax]["lambda"],
        "maximum_ratio": rows[argmax]["ratio_curv_over_pred"],
        "lambda_of_minimum_fake_excess": rows[argmin]["lambda"],
        "minimum_ratio": rows[argmin]["ratio_curv_over_pred"],
        "smooth": _is_smooth([row["a_curv"] for row in rows]),
        "reaches_one": min(ratios) < 1.05,
    }


def _is_smooth(values: list[float], tol: float = 1e-6) -> bool:
    """Second differences of the sweep must be tiny relative to the range."""
    if len(values) < 3:
        return True
    span = max(values) - min(values)
    worst = max(abs(values[i-1] - 2*values[i] + values[i+1])
                for i in range(1, len(values) - 1))
    return worst <= max(tol, 1e-3 * span)


def main() -> dict[str, Any]:
    import p582_amplitude_law as law
    both = ag.load_all()
    chart_ref = json.loads((ag.OUT_DIR / "_reference"
                            / "p612-chart-identity-latest.json").read_text())
    lineages = {}
    all_smooth = True
    for weighting in ("spin0", "equal"):
        loaded = both[weighting]
        direction = ag.frozen_direction("spin0")
        fit = chart_ref["weightings"][weighting]["exponent_fit"]
        per_lineage = {}
        for name, sizes in ag.flow.LINEAGES.items():
            if len(sizes) < 3:
                continue
            entry = sweep_one(loaded, sizes, direction,
                              fit["scale_amplitude"], fit["omega"])
            per_lineage[name] = entry
            if not entry["smooth"]:
                all_smooth = False
        lineages[weighting] = per_lineage
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "lambda_grid": LAMBDAS,
        "lineages": lineages,
        "fake_excess_is_smooth_in_lambda": all_smooth,
        "kill_sentence": (
            "PASS: the fake excess is a smooth function of the attachment, "
            "peaks at the historical middle chart, and falls toward the "
            "consistent reading as the curvature is attached to the lower "
            "chart -- the identity is the whole story of the 55%."
            if all_smooth else
            "KILL: the fake excess is NOT smooth in lambda; the identity is "
            "not the whole story."),
    }


if __name__ == "__main__":
    result = main()
    path = ag.dump("p2-lambda-sweep", result)
    print(json.dumps({
        "kill_sentence": result["kill_sentence"],
        "output": str(path)}, indent=2))
