#!/usr/bin/env python3
"""Parameter-free fixed-power diagnostics for a finite-width threshold sequence.

The primary statistic is the dyadic difference ratio

    (p(4n)-p(2n)) / (p(2n)-p(n)),

which tends to 2**(-delta) when the leading shift is A*n**(-delta).
No infinite-volume threshold is inserted or fitted for this statistic.
"""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal, getcontext
import json
import math
from pathlib import Path
from statistics import median
from typing import Iterable


getcontext().prec = 80


def read_sequence(path: Path, width_column: str, value_column: str) -> dict[int, Decimal]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle)
        values = {
            int(row[width_column]): Decimal(row[value_column])
            for row in rows
            if row.get(value_column, "").strip()
        }
    if len(values) < 3:
        raise ValueError("at least three non-empty sequence values are required")
    return values


def dyadic_diagnostics(values: dict[int, Decimal], delta: int = 4) -> list[dict[str, str | int]]:
    target = Decimal(2) ** Decimal(-delta)
    output: list[dict[str, str | int]] = []
    for n in sorted(values):
        if 2 * n not in values or 4 * n not in values:
            continue
        d_n = values[2 * n] - values[n]
        d_2n = values[4 * n] - values[2 * n]
        ratio = d_2n / d_n
        effective_delta = -math.log(abs(float(ratio)), 2)
        output.append(
            {
                "n": n,
                "d_n": str(d_n),
                "d_2n": str(d_2n),
                "ratio": str(ratio),
                "fixed_delta_target": str(target),
                "ratio_minus_target": str(ratio - target),
                "effective_delta": format(effective_delta, ".15g"),
            }
        )
    return output


def two_point_prediction_residuals(
    values: dict[int, Decimal], delta: int
) -> list[dict[str, str | int]]:
    """Fit p_inf+A*n^-delta to two consecutive widths and predict the next."""
    widths = sorted(values)
    output: list[dict[str, str | int]] = []
    power = Decimal(delta)
    for n0, n1, n2 in zip(widths, widths[1:], widths[2:]):
        if n1 != n0 + 1 or n2 != n1 + 1:
            continue
        x0 = Decimal(n0) ** (-power)
        x1 = Decimal(n1) ** (-power)
        x2 = Decimal(n2) ** (-power)
        p0, p1, p2 = values[n0], values[n1], values[n2]
        amplitude = (p1 - p0) / (x1 - x0)
        intercept = p0 - amplitude * x0
        prediction = intercept + amplitude * x2
        output.append(
            {
                "fit_width_0": n0,
                "fit_width_1": n1,
                "test_width": n2,
                "intercept": str(intercept),
                "prediction": str(prediction),
                "observed": str(p2),
                "residual": str(p2 - prediction),
            }
        )
    return output


def residual_summary(rows: Iterable[dict[str, str | int]], tail: int = 8) -> dict[str, str | int]:
    selected = list(rows)[-tail:]
    if not selected:
        return {
            "tail_count": 0,
            "tail_first_test_width": 0,
            "tail_last_test_width": 0,
            "median_absolute_residual": "NaN",
            "maximum_absolute_residual": "NaN",
        }
    absolute = [abs(Decimal(str(row["residual"]))) for row in selected]
    return {
        "tail_count": len(selected),
        "tail_first_test_width": int(selected[0]["test_width"]),
        "tail_last_test_width": int(selected[-1]["test_width"]),
        "median_absolute_residual": str(median(absolute)),
        "maximum_absolute_residual": str(max(absolute)),
    }


def build_score(values: dict[int, Decimal], source: Path) -> dict[str, object]:
    powers: dict[str, object] = {}
    for delta in (3, 4, 5, 6):
        rows = two_point_prediction_residuals(values, delta)
        powers[str(delta)] = {
            "tail_summary": residual_summary(rows),
            "rolling_predictions": rows,
        }
    return {
        "schema": "matching-one.p321-fixed-power-sequence-score.v1",
        "source": str(source),
        "width_min": min(values),
        "width_max": max(values),
        "primary_hypothesis": {
            "delta": 4,
            "dyadic_target": "0.0625",
            "uses_fitted_p_infinity": False,
        },
        "dyadic_diagnostics": dyadic_diagnostics(values, delta=4),
        "fixed_power_rolling_controls": powers,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--width-column", default="n")
    parser.add_argument("--value-column", default="value")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    values = read_sequence(args.input, args.width_column, args.value_column)
    score = build_score(values, args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(score, indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
