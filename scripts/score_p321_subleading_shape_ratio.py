#!/usr/bin/env python3
"""Score the empirical P321 subleading width ratio with Fieller sets.

For equal-area rectangles the fitted coefficients obey

    C_width = C_N / rho**2,
    D_width = D_N / rho**3,

so the requested ratio is D_N / (rho*C_N).  Several fitted C_N values are
near zero, making symmetric delta-method intervals unreliable.  This scorer
therefore reports both the local delta scale and the full Fieller confidence
set implied by the joint GLS covariance of (C_N,D_N).
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from statistics import NormalDist
from typing import Mapping, Sequence


RHO_ORDER = ("1", "16/9", "9/4", "4", "9")
RHO_VALUES = (1.0, 16.0 / 9.0, 9.0 / 4.0, 4.0, 9.0)


def fieller_set(
    x: float,
    y: float,
    var_x: float,
    var_y: float,
    cov_xy: float,
    *,
    alpha: float,
) -> dict:
    """Return the Gaussian Fieller set for the ratio y/x."""

    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0,1)")
    z = NormalDist().inv_cdf(1.0 - alpha / 2.0)
    a = x * x - z * z * var_x
    b = -2.0 * x * y + 2.0 * z * z * cov_xy
    c = y * y - z * z * var_y
    discriminant = b * b - 4.0 * a * c
    scale = max(1.0, abs(x * x), abs(z * z * var_x))
    tolerance = 1e-14 * scale
    result = {
        "alpha": alpha,
        "z": z,
        "quadratic": {"a": a, "b": b, "c": c, "discriminant": discriminant},
        "denominator_excludes_zero": a > 0.0,
    }
    if abs(a) <= tolerance:
        if abs(b) <= tolerance:
            result.update({"kind": "all_real" if c <= 0.0 else "empty"})
        else:
            bound = -c / b
            result.update(
                {
                    "kind": "half_line",
                    "bound": bound,
                    "direction": "at_most" if b > 0.0 else "at_least",
                }
            )
        return result
    if discriminant < 0.0:
        result.update({"kind": "all_real" if a < 0.0 else "empty"})
        return result
    root = math.sqrt(max(0.0, discriminant))
    lower, upper = sorted(((-b - root) / (2.0 * a), (-b + root) / (2.0 * a)))
    if a > 0.0:
        result.update({"kind": "bounded", "lower": lower, "upper": upper})
    else:
        result.update(
            {
                "kind": "disjoint_unbounded",
                "excluded_open_interval": [lower, upper],
            }
        )
    return result


def analyze(payload: Mapping) -> dict:
    order = list(payload["parameter_order"])
    expected = ["pc"] + [f"C_N[{rho}]" for rho in RHO_ORDER] + [
        f"D_N[{rho}]" for rho in RHO_ORDER
    ]
    if order != expected:
        raise ValueError("unexpected GLS parameter order")
    estimate = [float(value) for value in payload["estimate"]]
    covariance = [[float(value) for value in row] for row in payload["covariance"]]
    if len(covariance) != 11 or any(len(row) != 11 for row in covariance):
        raise ValueError("expected an 11 by 11 GLS covariance")

    rows = []
    for index, (rho_label, rho) in enumerate(zip(RHO_ORDER, RHO_VALUES)):
        c_index = 1 + index
        d_index = 6 + index
        x = estimate[c_index]
        y = estimate[d_index] / rho
        var_x = covariance[c_index][c_index]
        var_y = covariance[d_index][d_index] / (rho * rho)
        cov_xy = covariance[c_index][d_index] / rho
        ratio = y / x
        grad_x = -y / (x * x)
        grad_y = 1.0 / x
        delta_variance = (
            grad_x * grad_x * var_x
            + grad_y * grad_y * var_y
            + 2.0 * grad_x * grad_y * cov_xy
        )
        rows.append(
            {
                "rho": rho_label,
                "C_N": x,
                "D_N": estimate[d_index],
                "C_denominator_z": x / math.sqrt(var_x),
                "D_width_over_C_width": ratio,
                "delta_standard_error_diagnostic": math.sqrt(max(0.0, delta_variance)),
                "fieller_95": fieller_set(
                    x, y, var_x, var_y, cov_xy, alpha=0.05
                ),
                "fieller_99": fieller_set(
                    x, y, var_x, var_y, cov_xy, alpha=0.01
                ),
            }
        )

    bounded_99 = [row["rho"] for row in rows if row["fieller_99"]["kind"] == "bounded"]
    nonzero_99 = [
        row["rho"]
        for row in rows
        if row["fieller_99"]["kind"] == "bounded"
        and not row["fieller_99"]["lower"] <= 0.0 <= row["fieller_99"]["upper"]
    ]
    return {
        "schema": "matching-one/p321-subleading-shape-ratio/v1",
        "ratio_definition": "D_width/C_width = D_N/(rho*C_N)",
        "rows": rows,
        "decision": {
            "bounded_at_99_percent": bounded_99,
            "ratio_excludes_zero_at_99_percent": nonzero_99,
            "status": (
                "all_geometries_identified"
                if len(nonzero_99) == len(rows)
                else "one_bounded_set_but_no_nonzero_ratio"
                if bounded_99 and not nonzero_99
                else "partial_ratio_identification"
                if nonzero_99
                else "ratio_curve_not_identified"
            ),
            "interpretation": (
                "Do not compare a five-point D/C curve to a continuum model unless "
                "the corresponding Fieller set is bounded. Near-zero C_N is a "
                "denominator-identification failure, not evidence for a large ratio."
            ),
        },
        "claim_boundary": (
            "This is an empirical coefficient-ratio audit under the fixed N^-2+N^-3 "
            "fit. It neither supplies the missing thermal one-point F_t nor selects "
            "identity dressing or a Jordan extension."
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gls", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = json.loads(args.gls.read_text(encoding="utf-8"))
    result = analyze(payload)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
