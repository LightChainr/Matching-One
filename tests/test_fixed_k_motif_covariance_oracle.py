from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import fixed_k_motif_covariance_oracle as oracle  # noqa: E402


class FixedKMotifCovarianceOracleTests(unittest.TestCase):
    def test_inclusion_probability_is_exact(self) -> None:
        self.assertEqual(oracle.inclusion_probability(6, 3, 2), Fraction(1, 5))
        self.assertEqual(oracle.inclusion_probability(6, 1, 2), 0)
        with self.assertRaises(ValueError):
            oracle.inclusion_probability(6, 7, 2)

    def test_overlap_oracle_matches_every_fixed_k_enumeration(self) -> None:
        n, families, contrasts = oracle.fixture()
        for k in range(n + 1):
            self.assertEqual(
                oracle.signed_moments(n, k, families, contrasts),
                oracle.brute_force_signed_moments(n, k, families, contrasts),
            )

    def test_equal_multiplicity_contrasts_have_zero_mean_for_every_k(self) -> None:
        n, families, contrasts = oracle.fixture()
        for k in range(n + 1):
            means, _ = oracle.signed_moments(n, k, families, contrasts)
            self.assertEqual(means, {"edge_difference": 0, "triangle_difference": 0})

    def test_covariance_matrix_is_symmetric_positive_semidefinite(self) -> None:
        n, families, contrasts = oracle.fixture()
        _, covariance = oracle.signed_moments(n, 3, families, contrasts)
        a = covariance["edge_difference"]["edge_difference"]
        b = covariance["edge_difference"]["triangle_difference"]
        c = covariance["triangle_difference"]["triangle_difference"]
        self.assertEqual(b, covariance["triangle_difference"]["edge_difference"])
        self.assertGreater(a, 0)
        self.assertGreater(c, 0)
        self.assertGreaterEqual(a * c - b * b, 0)

    def test_identical_signed_counts_have_zero_variance(self) -> None:
        family = oracle.normalize_family(4, ((0, 1), (1, 2)))
        contrasts = {"zero": {"left": Fraction(1), "right": Fraction(-1)}}
        means, covariance = oracle.signed_moments(
            4, 2, {"left": family, "right": family}, contrasts
        )
        self.assertEqual(means["zero"], 0)
        self.assertEqual(covariance["zero"]["zero"], 0)

    def test_invalid_embeddings_and_unknown_families_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "distinct"):
            oracle.normalize_family(4, ((0, 0),))
        with self.assertRaisesRegex(ValueError, "outside"):
            oracle.normalize_family(4, ((0, 4),))
        with self.assertRaisesRegex(ValueError, "unknown"):
            oracle.signed_moments(4, 2, {}, {"bad": {"missing": Fraction(1)}})

    def test_checked_artifact_reproduces(self) -> None:
        checked = json.loads(
            (ROOT / "results/fixed-k-motif-covariance/latest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(checked, oracle.build_artifact())
        self.assertEqual(checked["oracle_vs_enumeration_failures"], 0)


if __name__ == "__main__":
    unittest.main()
