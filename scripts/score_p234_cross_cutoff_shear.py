#!/usr/bin/env python3
"""Joint cross-cutoff score for the P234 logarithmic top-field shear."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

try:
    from scripts.score_p234_fixed_delta_continuum import _gls, read_size
except ModuleNotFoundError:  # Direct `python scripts/...` execution.
    from score_p234_fixed_delta_continuum import _gls, read_size


def realized_cutoff_row(path: Path) -> dict:
    row = read_size(path)
    ratio = (row["realized_delta"] / row["declared_delta"]) ** (-25.0 / 24.0)
    scaling = [1.0, ratio, ratio * ratio]
    point = [value * scale for value, scale in zip(row["normalized_point"], scaling)]
    covariance = row["normalized_covariance_delta_method"]
    covariance = [
        [covariance[i][j] * scaling[i] * scaling[j] for j in range(3)]
        for i in range(3)
    ]
    return {
        **row,
        "natural_realized_cutoff_point": point,
        "natural_realized_cutoff_covariance": covariance,
    }


def cross_cutoff_score(rows: Sequence[dict]) -> dict:
    points = [row["natural_realized_cutoff_point"] for row in rows]
    covariances = [row["natural_realized_cutoff_covariance"] for row in rows]
    width = 3 * len(rows)
    columns = [[0.0] * width for _ in range(6)]
    for block, row in enumerate(rows):
        inverse_size = 1.0 / float(row["L"])
        inverse_radius = 1.0 / (float(row["L"]) * float(row["realized_delta"]))
        log_cutoff = math.log(2.0 * float(row["realized_delta"]))
        columns[0][3 * block] = inverse_size
        columns[1][3 * block + 1] = 1.0
        columns[2][3 * block + 1] = inverse_radius
        columns[3][3 * block + 2] = 1.0
        columns[4][3 * block + 2] = log_cutoff
        columns[5][3 * block + 2] = inverse_radius
    score = _gls(
        points,
        covariances,
        columns,
        [
            "LL_1_over_L",
            "LD_continuum",
            "LD_1_over_lattice_radius",
            "DD_intercept",
            "DD_log_2delta_slope",
            "DD_1_over_lattice_radius",
        ],
        "LL=c/L; LD=B+c/(L delta); DD=D+s log(2 delta)+c/(L delta)",
    )
    mixed = score["coefficients"][1]
    slope = score["coefficients"][4]
    kappa = -slope / (2.0 * mixed)
    gradient_mixed = slope / (2.0 * mixed**2)
    gradient_slope = -1.0 / (2.0 * mixed)
    covariance = score["coefficient_covariance"]
    variance = (
        gradient_mixed**2 * covariance[1][1]
        + gradient_slope**2 * covariance[4][4]
        + 2.0 * gradient_mixed * gradient_slope * covariance[1][4]
    )
    return {
        "gls": score,
        "kappa_proxy": {
            "definition": "-d(DD)/dlog(2delta)/(2*LD)",
            "estimate": kappa,
            "standard_error_delta_method": math.sqrt(max(0.0, variance)),
        },
    }


def render(paths: Sequence[Path]) -> dict:
    rows = sorted(
        (realized_cutoff_row(path) for path in paths),
        key=lambda row: (-row["declared_delta"], row["L"]),
    )
    cutoffs = sorted({row["declared_delta"] for row in rows}, reverse=True)
    if len(cutoffs) < 2:
        raise ValueError("at least two declared cutoffs are required")
    return {
        "schema": "matching-one.p234-cross-cutoff-shear.v1",
        "issue": 234,
        "status": "scorer_frozen_after_L64_preview_before_L96_L128_L192_phaseB_reveal",
        "rows": rows,
        "declared_cutoffs": cutoffs,
        "joint_score": cross_cutoff_score(rows),
        "scope": [
            "The natural realized cutoff removes the deterministic nearest-vertex prefactor before fitting.",
            "The leading finite-radius correction is 1/(L*delta_realized), not a free correction per cutoff.",
            "The LCFT relation B_delta=hat_phi-kappa*log(2delta)*phi predicts a common LD and a DD slope -2*kappa*LD.",
            "kappa_proxy remains in the connection-probability field gauge; cross-lattice universality needs an additional amplitude invariant.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = render(args.inputs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
