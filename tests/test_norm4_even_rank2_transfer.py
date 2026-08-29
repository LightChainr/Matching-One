#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_norm4_even_rank2_transfer import (  # noqa: E402
    next_point,
    rank_one_gls,
    render,
)


class Norm4EvenRank2TransferTests(unittest.TestCase):
    def test_exact_rank_one_residual(self) -> None:
        first = [1.0, -2.0, 3.0, -4.0, 5.0]
        residual = first + [0.25 * value for value in first]
        covariance = [[float(i == j) for j in range(10)] for i in range(10)]
        chi_square, direction, error = rank_one_gls(
            residual, covariance, mp.mpf("0.25")
        )
        self.assertAlmostEqual(float(chi_square), 0.0)
        self.assertEqual([float(value) for value in direction], first)
        self.assertTrue(all(abs(float(value)) < 1e-12 for value in error))

    def test_half_eigenvalue_recurrence(self) -> None:
        self.assertEqual(next_point([1.0], [3.0], [2.0], 0.5), [6.0])

    def test_committed_fit_regression(self) -> None:
        result = render(
            ROOT / "results/server-20260829/P154-norm4-production/analysis/scalar_score.json",
            ROOT / "results/server-20260829/P154-norm4-production/analysis/thermal_jet_score.json",
            0.5,
        )
        fit = result["covariance_aware_rank_one_fit"]
        self.assertAlmostEqual(fit["lineage_85_over_65_amplitude_ratio"], 0.17775649, places=6)
        self.assertAlmostEqual(fit["chi_square"], 4.75632317, places=6)
        self.assertEqual(fit["degrees_of_freedom"], 4)


if __name__ == "__main__":
    unittest.main()
