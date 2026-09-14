#!/usr/bin/env python3
"""Evaluate the fixed-width dual-odd charge scaling function.

Uses the transparent safe transfer from fixed_width_charge_transfer.py and
reports

    F_w(X) = w * Theta_w(h_w + X*w^(-3/4)),

where h is Bernoulli logit and h_w is the charge-coexistence root.  It also
fits the antisymmetric part to a local odd polynomial and estimates the
quadratic analytic thermal-metric coefficient from the symmetric residual.

This is a small-width deterministic transfer control, not a continuum proof.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from fixed_width_charge_transfer import SafeTransfer


def logistic(h: float) -> float:
    return 1.0 / (1.0 + math.exp(-h))


def logit(p: float) -> float:
    return math.log(p / (1.0 - p))


def run_width(width: int, x_grid: list[float]) -> dict[str, object]:
    g4 = SafeTransfer(width, matching=False)
    g8 = SafeTransfer(width, matching=True)

    def theta_p(p: float) -> float:
        return g4.perron(p)["I0"] - g8.perron(1.0 - p)["I0"]

    root = brentq(theta_p, 0.5000000001, 0.999999999, xtol=2e-14)
    h_root = logit(root)

    def scaled(x: float) -> float:
        h = h_root + x * width ** (-0.75)
        return width * theta_p(logistic(h))

    curve = {format(x, ".12g"): scaled(x) for x in x_grid}

    positive = sorted({abs(x) for x in x_grid if x != 0.0})
    fit_x = []
    fit_y = []
    for x in positive:
        if x <= 1.5:
            odd = 0.5 * (scaled(x) - scaled(-x))
            fit_x.extend([x, -x])
            fit_y.extend([odd, -odd])
    fit_x = np.asarray(fit_x)
    fit_y = np.asarray(fit_y)
    design = np.column_stack([fit_x, fit_x**3, fit_x**5])
    a1, a3, a5 = np.linalg.lstsq(design, fit_y, rcond=None)[0]

    def odd_fit_prime(x: float) -> float:
        return a1 + 3.0 * a3 * x * x + 5.0 * a5 * x**4

    metric_points = [x for x in (0.5, 0.75, 1.0, 1.5) if x in positive]
    numerator = 0.0
    denominator = 0.0
    point_estimates = {}
    for x in metric_points:
        symmetric_sum = scaled(x) + scaled(-x)
        target = 2.0 * x * x * odd_fit_prime(x) * width ** (-0.75)
        estimate = symmetric_sum / target
        point_estimates[format(x, ".12g")] = estimate
        numerator += symmetric_sum * target
        denominator += target * target
    c2 = numerator / denominator

    return {
        "width": width,
        "charge_root_p": root,
        "charge_root_h": h_root,
        "scaled_curve": curve,
        "odd_polynomial": {
            "a1_X": float(a1),
            "a3_X3": float(a3),
            "a5_X5": float(a5),
            "fifth_derivative_120a5": float(120.0 * a5),
        },
        "quadratic_metric_estimate": float(c2),
        "pointwise_metric_estimates": point_estimates,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-width", type=int, default=4)
    parser.add_argument("--max-width", type=int, default=9)
    parser.add_argument(
        "--x-grid",
        default="-1.5,-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1,1.5",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    x_grid = [float(item) for item in args.x_grid.split(",")]
    result = {
        "schema": "dual-odd-charge-scaling-collapse-v1",
        "definition": "F_w(X)=w*Theta_w(h_w+X*w^(-3/4))",
        "claim_boundary": [
            "deterministic small-width transfer control",
            "odd polynomial is only a local descriptive fit",
            "quadratic metric interpretation is a scaling hypothesis",
        ],
        "records": [
            run_width(width, x_grid)
            for width in range(args.min_width, args.max_width + 1)
        ],
    }
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
