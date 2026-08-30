#!/usr/bin/env python3
"""Pool independent P321 campaigns by precision before the fixed scale fit."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import mpmath as mp

from score_p321_equal_area_rectangles import (
    DEFAULT_DESIGN,
    DEFAULT_ORACLE,
    RHO_ORDER,
    fit_fixed_model,
    score_campaign,
)


def _contrast_covariance(covariance: Sequence[Sequence[float]]) -> list[list[float]]:
    rows = [
        [-1.0 if column == 0 else (1.0 if column == index else 0.0) for column in range(5)]
        for index in range(1, 5)
    ]
    return [
        [
            math.fsum(
                rows[i][a] * covariance[a][b] * rows[j][b]
                for a in range(5)
                for b in range(5)
            )
            for j in range(4)
        ]
        for i in range(4)
    ]


def combine_scored_campaigns(
    campaigns: Sequence[Mapping[str, Any]], sources: Sequence[str]
) -> dict[str, Any]:
    """Combine independent root vectors with their full within-run covariance."""

    if not campaigns or len(campaigns) != len(sources):
        raise ValueError("nonempty campaigns and one source per campaign required")
    n_values = {int(row["N"]) for row in campaigns}
    if len(n_values) != 1:
        raise ValueError("only campaigns at one common N can be combined")
    precision = mp.matrix(5, 5)
    rhs = mp.matrix(5, 1)
    for row in campaigns:
        covariance = mp.matrix(row["root_covariance"])
        inverse = mp.inverse(covariance)
        precision += inverse
        rhs += inverse * mp.matrix(row["roots"])
    covariance_mp = mp.inverse(precision)
    roots_mp = covariance_mp * rhs
    covariance = [[float(covariance_mp[i, j]) for j in range(5)] for i in range(5)]
    roots = [float(roots_mp[i]) for i in range(5)]
    contrast_covariance = _contrast_covariance(covariance)
    return {
        "N": n_values.pop(),
        "root_order": list(RHO_ORDER),
        "roots": roots,
        "root_covariance": covariance,
        "root_standard_errors": [math.sqrt(max(0.0, covariance[i][i])) for i in range(5)],
        "contrast_order": list(RHO_ORDER[1:]),
        "contrasts_to_square": [roots[i] - roots[0] for i in range(1, 5)],
        "contrast_covariance": contrast_covariance,
        "contrast_standard_errors": [
            math.sqrt(max(0.0, contrast_covariance[i][i])) for i in range(4)
        ],
        "batches": sum(int(row["batches"]) for row in campaigns),
        "samples_per_shape": sum(int(row["samples_per_shape"]) for row in campaigns),
        "square_histograms_byte_identical_within_each_component": all(
            bool(row["square_histograms_byte_identical"]) for row in campaigns
        ),
        "square_moments_byte_identical_within_each_component": all(
            bool(row["square_moments_byte_identical"]) for row in campaigns
        ),
        "elapsed_seconds_all_four_pairs": math.fsum(
            float(row["elapsed_seconds_all_four_pairs"]) for row in campaigns
        ),
        "independent_precision_pooling": True,
        "components": [
            {
                "source": source,
                "samples_per_shape": int(row["samples_per_shape"]),
                "batches": int(row["batches"]),
            }
            for source, row in zip(sources, campaigns)
        ],
    }


def score_independent(
    campaign_dirs: Sequence[Path],
    design_path: Path = DEFAULT_DESIGN,
    oracle_path: Path = DEFAULT_ORACLE,
) -> dict[str, Any]:
    design = json.loads(design_path.read_text(encoding="utf-8"))
    oracle = json.loads(oracle_path.read_text(encoding="utf-8"))
    scored = [(path, score_campaign(path, design)) for path in campaign_dirs]
    grouped: dict[int, list[tuple[Path, Mapping[str, Any]]]] = {}
    for path, row in scored:
        grouped.setdefault(int(row["N"]), []).append((path, row))
    combined = [
        combine_scored_campaigns(
            [row for _, row in rows], [str(path) for path, _ in rows]
        )
        for _, rows in sorted(grouped.items())
    ]
    if sorted(grouped) != [144, 576, 1296]:
        raise ValueError("independent P321 score requires N=144,576,1296")
    return {
        "schema": "matching-one/p321-independent-campaign-pool/v1",
        "data_class": "independent RNG domains pooled by full root precision",
        "campaigns": combined,
        "scale_fit": fit_fixed_model(combined, oracle),
        "pooling_boundary": (
            "Independence is asserted by disjoint declared RNG domains; within-run "
            "cross-rho covariance is retained and no cross-run covariance is fitted."
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-dir", action="append", required=True, type=Path)
    parser.add_argument("--design", type=Path, default=DEFAULT_DESIGN)
    parser.add_argument("--oracle", type=Path, default=DEFAULT_ORACLE)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    mp.mp.dps = 60
    result = score_independent(args.campaign_dir, args.design, args.oracle)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, allow_nan=False, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
