#!/usr/bin/env python3
"""Numerically scaled entry point for the adjacent-axis annihilator scorer.

The physical models in :mod:`score_axis_pair_annihilator` use design columns
whose magnitudes can differ by many orders (`L^-w` versus `1`, or
`Delta L^-q` versus `Delta L^4`).  This module keeps the exact same model and
output contract but rescales each weighted-least-squares design column before
forming normal equations, then maps coefficients/covariances back to the
original basis.

No scientific parameter, exponent, train/heldout split or target changes.
"""
from __future__ import annotations

from typing import Sequence

import score_axis_pair_annihilator as base


_original_weighted_linear_fit = base.weighted_linear_fit


def scaled_weighted_linear_fit(
    features: Sequence[Sequence[float]],
    values: Sequence[float],
    variances: Sequence[float],
):
    if not features or not features[0]:
        raise ValueError("empty design matrix")
    width = len(features[0])
    if any(len(row) != width for row in features):
        raise ValueError("ragged design matrix")

    scales = []
    for j in range(width):
        scale = max(abs(row[j]) for row in features)
        if not scale > 0.0:
            raise ValueError("zero design column")
        scales.append(scale)

    normalized = [
        [row[j] / scales[j] for j in range(width)]
        for row in features
    ]
    beta_scaled, cov_scaled, residual, chi2 = _original_weighted_linear_fit(
        normalized, values, variances
    )

    beta = [beta_scaled[j] / scales[j] for j in range(width)]
    covariance = [
        [cov_scaled[i][j] / (scales[i] * scales[j]) for j in range(width)]
        for i in range(width)
    ]
    return beta, covariance, residual, chi2


# fit_f_shape / fit_root_power resolve this module global from base at call time.
base.weighted_linear_fit = scaled_weighted_linear_fit

fit_f_shape = base.fit_f_shape
fit_root_power = base.fit_root_power
calculate = base.calculate
report = base.report


def main() -> int:
    return base.main()


if __name__ == "__main__":
    raise SystemExit(main())
