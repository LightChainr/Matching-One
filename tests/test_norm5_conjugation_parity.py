#!/usr/bin/env python3


from __future__ import annotations
from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from diagnose_norm5_conjugation_parity import project_parity  # noqa: E402


class Norm5ConjugationParityTests(unittest.TestCase):
    def test_projection_separates_even_and_odd_lineage_vectors(self) -> None:
        residual = [mp.mpf(1), mp.mpf(2), mp.mpf(3), mp.mpf(0)]
        covariance = [
            [mp.mpf(int(i == j)) for j in range(4)] for i in range(4)
        ]
        even, odd, even_cov, odd_cov, cross_cov = project_parity(
            residual, covariance
        )
        self.assertEqual(even, [mp.mpf(2), mp.mpf(1)])
        self.assertEqual(odd, [mp.mpf(1), mp.mpf(-1)])
        self.assertEqual(even_cov, [[mp.mpf("0.5"), 0], [0, mp.mpf("0.5")]])
        self.assertEqual(odd_cov, [[mp.mpf("0.5"), 0], [0, mp.mpf("0.5")]])
        self.assertEqual(cross_cov, [[0, 0], [0, 0]])

    def test_projection_rejects_ragged_covariance(self) -> None:
        with self.assertRaises(ValueError):
            project_parity([mp.mpf(1), mp.mpf(2)], [[mp.mpf(1)]])


if __name__ == "__main__":
    unittest.main()
