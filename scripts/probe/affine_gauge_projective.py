#!/usr/bin/env python3
"""P10 — Gate 2 projective design vs this groupoid.

#595/#596 Gate 2 ranked channels by a projective figure of merit
(Fieller / A4^2/var(delta)); #602 forbade retargeting Gate 3's coordinate
with that ranking.  On the same eight committed blocks:

  - the Gate-2 projective amplitude of g as a function of lambda (P2's
    path): the projective coordinate rho(lambda) = delta(lambda)/A4(lambda)
    where A4 is the curvature's g-amplitude and delta its width
    coefficient (the predeclared denominator pair of P6-I4);
  - whether projective coordinates KILL the 1.55 fake excess (they should,
    if they are attachment-invariant) or merely hide it in a denominator.

Kill: if the projective rho still moves by O(50%) along lambda, Gate 2 did
not buy chart invariance and #602's warning is stronger than written.
Report that.  Do not retarget Gate 3.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import affine_gauge_common as ag  # noqa: E402
from p612_chart_identity import LEVELS  # noqa: E402

SCHEMA = "matching-one.probe-affine-gauge.p10.v1"
ISSUE = 617

LAMBDAS = [round(i / 10, 2) for i in range(11)]   # 0.0 ... 1.0
O50_THRESHOLD = 0.5


def _rho_at(loaded: dict[int, dict[str, Any]], sizes: tuple[int, ...],
            direction: list[float], lam: float) -> dict[str, Any]:
    curvature = ag.curvature_geometry_at(loaded, sizes, direction, lam)
    a4 = curvature["amplitude"]
    var_a4 = curvature["standard_error"] ** 2
    beta_covector = curvature["covectors"][1]
    weights = curvature["second_difference_weights"]
    observation = [math.fsum(w * loaded[size]["quantiles"][level]
                             for w, size in zip(weights, sizes))
                   for level in range(LEVELS)]
    delta = ag._apply(beta_covector, observation)
    cov_obs = [[math.fsum(w1 * w2 * loaded[size]["covariance"][i][j]
                          for w1 in weights for w2 in weights
                          for size in sizes)
                for j in range(LEVELS)] for i in range(LEVELS)]
    var_delta = math.fsum(beta_covector[i] * cov_obs[i][j] * beta_covector[j]
                          for i in range(LEVELS) for j in range(LEVELS))
    cov_ad = math.fsum(
        curvature["covectors"][2][i] * cov_obs[i][j] * beta_covector[j]
        for i in range(LEVELS) for j in range(LEVELS))
    # Fieller set for rho = delta / A4 at 95% (the Gate 2 figure of merit)
    z2 = 1.959963984540054 ** 2
    a = a4 * a4 - z2 * var_a4
    b = -2.0 * (a4 * delta - z2 * cov_ad)
    c = delta * delta - z2 * var_delta
    fieller = None
    if a > 0:
        disc = b * b - 4 * a * c
        if disc >= 0:
            root = math.sqrt(disc)
            fieller = [(-b - root) / (2 * a), (-b + root) / (2 * a)]
    return {
        "lambda": lam,
        "A4": a4,
        "delta_width": delta,
        "rho_point": delta / a4 if a4 else None,
        "fieller_set_95": fieller,
        "A4_sign": "+" if a4 > 0 else "-",
    }


def main() -> dict[str, Any]:
    both = ag.load_all()
    out = {}
    for weighting in ("spin0", "equal"):
        loaded = both[weighting]
        direction = ag.frozen_direction("spin0")
        per_lineage = {}
        for lineage_name, sizes in ag.flow.LINEAGES.items():
            if len(sizes) < 3:
                continue
            rows = [_rho_at(loaded, sizes, direction, lam)
                    for lam in LAMBDAS]
            rhos = [row["rho_point"] for row in rows]
            span = (max(rhos) - min(rhos)) / abs(rhos[-1]) if rhos[-1] else None
            signs = {row["A4_sign"] for row in rows}
            fieller_all_bounded = all(row["fieller_set_95"] is not None
                                      for row in rows)
            per_lineage[lineage_name] = {
                "sizes": list(sizes),
                "rows": rows,
                "rho_relative_span_across_lambda": span,
                "rho_moves_O50": span is not None and span > O50_THRESHOLD,
                "A4_sign_flips": len(signs) > 1,
                "fieller_bounded_everywhere": fieller_all_bounded,
            }
        out[weighting] = per_lineage
    spans = [v["rho_relative_span_across_lambda"]
             for w in out.values() for v in w.values()
             if v["rho_relative_span_across_lambda"] is not None]
    any_O50 = any(v["rho_moves_O50"] for w in out.values()
                  for v in w.values())
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "lambda_grid": LAMBDAS,
        "lineages": out,
        "max_relative_span": max(spans) if spans else None,
        "projective_moves_O50_along_lambda": any_O50,
        "verdict": (
            "KILL: the projective rho still moves by O(50%%) along the "
            "attachment path (max relative span %.1f%%) -- Gate 2 did NOT "
            "buy chart invariance, and #602's warning is stronger than "
            "written.  Do not retarget Gate 3."
            % (100.0 * max(spans) if spans else 0.0)
            if any_O50 else
            "PASS: the projective coordinate is attachment-stable (max "
            "relative span %.1f%% < 50%%): the projective design kills "
            "the fake excess rather than hiding it, and I4's survival in "
            "P6 is confirmed along the full path."
            % (100.0 * max(spans) if spans else 0.0)),
    }


if __name__ == "__main__":
    result = main()
    path = ag.dump("p10-projective", result)
    print(json.dumps({
        "max_relative_span": result["max_relative_span"],
        "verdict": result["verdict"],
        "output": str(path)}, indent=2))
