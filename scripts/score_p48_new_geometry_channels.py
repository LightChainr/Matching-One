#!/usr/bin/env python3
"""Score the frozen P48 pure laws on prospective N=185/265 geometries.

The source amplitudes are read from the retrospective P48 training fit
(N=65,85,130). The target values come from the independent N=185/265
full-curve production. No target parameter is fitted.

For each scaled channel y_N the two target sizes are independent, while the
single frozen source-amplitude uncertainty is common to both predictions:

    Cov(y_target - A) = diag(target_se**2) + Var(A) * 11^T.

This script scores the intrinsic-center P48 projectors. It must not be confused
with the separate Issue #43 fixed-coordinate wrapping-channel DeltaS contract.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List


CHANNELS = {
    "P4_S": {
        "power": 1.0,
        "scaled_key": "A_S_N1",
        "law": "P4[S] ~ N^-1",
    },
    "P4_D": {
        "power": 13.0 / 8.0,
        "scaled_key": "A_D_N13_8",
        "law": "P4[D] ~ N^-13/8",
    },
    "P4_S_prime": {
        "power": 5.0 / 4.0,
        "scaled_key": "A_Sprime_N5_4",
        "law": "P4[S'] ~ N^-5/4",
    },
    "P4_D_prime": {
        "power": 5.0 / 8.0,
        "scaled_key": "A_Dprime_N5_8",
        "law": "P4[D'] ~ N^-5/8",
    },
}


def _load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _two_by_two_chi_square(residual: List[float], covariance: List[List[float]]) -> float:
    c00, c01 = covariance[0]
    c10, c11 = covariance[1]
    if abs(c01 - c10) > 1e-15 * max(1.0, abs(c01), abs(c10)):
        raise ValueError("covariance is not symmetric")
    determinant = c00 * c11 - c01 * c10
    if determinant <= 0.0:
        raise ValueError("residual covariance is not positive definite")
    r0, r1 = residual
    return (r0 * r0 * c11 - 2.0 * r0 * r1 * c01 + r1 * r1 * c00) / determinant


def score(source_summary: dict, target_payloads: Dict[int, dict]) -> dict:
    sizes = sorted(target_payloads)
    if sizes != [185, 265]:
        raise ValueError("this frozen score requires exactly N=185 and N=265")

    output = {
        "schema": "P48 prospective new-geometry pure-law score v1",
        "source_training_sizes": [65, 85, 130],
        "target_sizes": sizes,
        "target_refit_parameters": 0,
        "target_independence": "N=185 and N=265 use disjoint counter domains",
        "important_observable_boundary": (
            "intrinsic-center P48 projectors; distinct from the Issue #43 "
            "fixed-coordinate cross/either DeltaS contract"
        ),
        "channels": {},
    }

    source_scores = source_summary["channel_scores"]
    for channel, specification in CHANNELS.items():
        frozen = source_scores[channel]
        amplitude = float(frozen["amplitude"])
        amplitude_se = float(frozen["amplitude_se"])
        observed = []
        observed_se = []

        for n in sizes:
            row = target_payloads[n]["by_N"][str(n)]
            scaled_value = float(row["scaled"][specification["scaled_key"]])
            point_se = float(row["se"][channel])
            scaled_se = point_se * n ** specification["power"]
            observed.append(scaled_value)
            observed_se.append(scaled_se)

        residual = [value - amplitude for value in observed]
        common_variance = amplitude_se * amplitude_se
        covariance = [
            [observed_se[0] ** 2 + common_variance, common_variance],
            [common_variance, observed_se[1] ** 2 + common_variance],
        ]
        chi_square = _two_by_two_chi_square(residual, covariance)
        zero_chi_square = sum((value / se) ** 2 for value, se in zip(observed, observed_se))
        marginal_z = [
            difference / math.sqrt(se * se + common_variance)
            for difference, se in zip(residual, observed_se)
        ]

        output["channels"][channel] = {
            "law": specification["law"],
            "frozen_source_amplitude": amplitude,
            "frozen_source_amplitude_se": amplitude_se,
            "observed_scaled": observed,
            "observed_scaled_se": observed_se,
            "residual": residual,
            "residual_covariance": covariance,
            "marginal_signed_z": marginal_z,
            "chi_square": chi_square,
            "df": 2,
            "zero_chi_square": zero_chi_square,
        }

    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("results/server-20260828/P48-retrospective/summary.json"),
    )
    parser.add_argument(
        "--n185",
        type=Path,
        default=Path("results/server-20260828/P43-heldout-fullcurve-500m/analysis/n185.p48.json"),
    )
    parser.add_argument(
        "--n265",
        type=Path,
        default=Path("results/server-20260828/P43-heldout-fullcurve-500m/analysis/n265.p48.json"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    payload = score(
        _load_json(args.source),
        {185: _load_json(args.n185), 265: _load_json(args.n265)},
    )
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
