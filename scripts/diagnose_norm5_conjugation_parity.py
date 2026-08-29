#!/usr/bin/env python3
"""Project the norm-5 Jordan residual onto Gaussian-conjugation parity.

This is a post-reveal mechanism diagnostic.  It consumes the already-scored
ten-component residual and its full covariance; it does not create an
independent evidence row.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Sequence

import mpmath as mp

from score_norm5_thermal_jet import generalized_chi_square


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def project_parity(
    residual: Sequence[mp.mpf], covariance: Sequence[Sequence[mp.mpf]]
) -> tuple[
    list[mp.mpf],
    list[mp.mpf],
    list[list[mp.mpf]],
    list[list[mp.mpf]],
    list[list[mp.mpf]],
]:
    if len(residual) % 2:
        raise ValueError("parity projection requires two equal lineage blocks")
    width = len(residual) // 2
    if len(covariance) != len(residual) or any(
        len(row) != len(residual) for row in covariance
    ):
        raise ValueError("residual covariance dimension mismatch")
    half = mp.mpf("0.5")
    even = [half * (residual[i] + residual[i + width]) for i in range(width)]
    odd = [half * (residual[i + width] - residual[i]) for i in range(width)]

    quarter = mp.mpf("0.25")
    even_cov = [[quarter * (
        covariance[i][j] + covariance[i][j + width]
        + covariance[i + width][j] + covariance[i + width][j + width]
    ) for j in range(width)] for i in range(width)]
    odd_cov = [[quarter * (
        covariance[i][j] - covariance[i][j + width]
        - covariance[i + width][j] + covariance[i + width][j + width]
    ) for j in range(width)] for i in range(width)]
    even_odd_cov = [
        [quarter * (
            -covariance[i][j] + covariance[i][j + width]
            - covariance[i + width][j] + covariance[i + width][j + width]
        ) for j in range(width)] for i in range(width)
    ]
    return even, odd, even_cov, odd_cov, even_odd_cov


def strings(values: Sequence[mp.mpf]) -> list[str]:
    return [mp.nstr(value, 25) for value in values]


def matrix_strings(values: Sequence[Sequence[mp.mpf]]) -> list[list[str]]:
    return [strings(row) for row in values]


def render(source: dict, input_path: Path) -> dict:
    model = source["secondary_multiplier_cocycles"]["rank2_Jordan"]
    labels = list(model["labels"])
    midpoint = len(labels) // 2
    expected_left = [f"N65_to_N325_r{rank}" for rank in range(2, 7)]
    expected_right = [f"N85_to_N425_r{rank}" for rank in range(2, 7)]
    if labels[:midpoint] != expected_left or labels[midpoint:] != expected_right:
        raise ValueError("unexpected norm-5 lineage or rank ordering")
    residual = [mp.mpf(value) for value in model["residual"]]
    covariance = [
        [mp.mpf(value) for value in row] for row in model["covariance"]
    ]
    even, odd, even_cov, odd_cov, cross_cov = project_parity(residual, covariance)
    correlations = [
        cross_cov[i][i] / mp.sqrt(even_cov[i][i] * odd_cov[i][i])
        for i in range(midpoint)
    ]
    sine_scale = mp.mpf(25) / 24
    sine_quadrature = [sine_scale * value for value in odd]
    return {
        "schema": "matching-one/norm5-conjugation-parity/v1",
        "status": "post-reveal correlated mechanism diagnostic; not additive evidence",
        "input": str(input_path),
        "input_sha256": sha256(input_path),
        "model": "rank2_Jordan",
        "ranks": list(range(2, 7)),
        "lineages": {
            "minus": {"sizes": [65, 325], "multiplier": "2-i"},
            "plus": {"sizes": [85, 425], "multiplier": "2+i"},
        },
        "spin4_phase": {
            "cosine_both": "-7/25",
            "sine_minus": "-24/25",
            "sine_plus": "24/25",
        },
        "definition": {
            "even": "(residual_(2-i) + residual_(2+i))/2",
            "odd": "(residual_(2+i) - residual_(2-i))/2",
        },
        "even": {
            "residual": strings(even),
            "covariance": matrix_strings(even_cov),
            "score": generalized_chi_square(even, even_cov),
        },
        "odd": {
            "residual": strings(odd),
            "covariance": matrix_strings(odd_cov),
            "score": generalized_chi_square(odd, odd_cov),
            "sine_quadrature_estimate": strings(sine_quadrature),
        },
        "even_odd_cross_covariance": matrix_strings(cross_cov),
        "same_rank_even_odd_correlations": strings(correlations),
        "source_joint_score": model["score"],
        "locked_phase_node_prediction": {
            "multipliers": ["1-i", "1+i"],
            "cosine_both": "-1",
            "sine_both": "0",
            "target": "N145_to_N290_fullcurve",
            "prediction": (
                "a simple spin-4 sine-quadrature remainder has zero "
                "conjugation-odd amplitude at the norm-2 phase node"
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    mp.mp.dps = 80
    source = json.loads(args.input.read_text(encoding="utf-8"))
    payload = render(source, args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
