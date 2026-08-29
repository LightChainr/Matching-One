from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p123_correlated_hyperedge_closure as oracle  # noqa: E402


class CorrelatedHyperedgeClosureTests(unittest.TestCase):
    def test_duality_is_an_involution_and_normalization_is_preserved(self) -> None:
        tensor = oracle.independent_two_block_family()
        self.assertEqual(oracle.dual(oracle.dual(tensor)), tensor)
        self.assertEqual(oracle.normalization(tensor), oracle.poly(1))
        self.assertEqual(oracle.normalization(oracle.dual(tensor)), oracle.poly(1))

    def test_all_none_slice_has_only_half_selfdual(self) -> None:
        all_none = (oracle.poly(1, -1), oracle.poly(0), oracle.poly(0), oracle.poly(0, 1), oracle.poly(0), oracle.poly(0))
        residual = oracle.duality_residual(all_none)
        self.assertEqual(residual[0], oracle.poly(-1, 2))
        self.assertEqual(oracle.evaluate(residual[0], Fraction(1, 2)), 0)

    def test_composition_generates_unpaired_partial_state(self) -> None:
        tensor = oracle.independent_two_block_family()
        self.assertEqual(tensor[1], oracle.poly(0, Fraction(1, 2), Fraction(-1, 2)))
        residual = oracle.duality_residual(tensor)
        self.assertEqual(residual[0], oracle.poly(-1, 2))
        self.assertEqual(residual[1], oracle.poly(0, Fraction(-1, 2), Fraction(1, 2)))
        for point in (Fraction(0), Fraction(1)):
            self.assertNotEqual(oracle.evaluate(residual[0], point), 0)
        self.assertNotEqual(oracle.evaluate(residual[1], Fraction(1, 2)), 0)

    def test_positive_selfdual_simplex(self) -> None:
        for u, v in ((Fraction(0), Fraction(0)), (Fraction(1, 32), Fraction(1, 32))):
            point = oracle.selfdual_simplex_point(u, v)
            self.assertEqual(point, (point[3], point[2], point[1], point[0], point[5], point[4]))
            self.assertEqual(sum(m * x for m, x in zip(oracle.MULTIPLICITIES, point)), 1)
            self.assertLessEqual(point[3], Fraction(1, 2))

    def test_committed_artifact_is_reproducible(self) -> None:
        committed = json.loads(
            (ROOT / "results" / "exact-correlated-hyperedge-closure" / "latest.json").read_text()
        )
        self.assertEqual(committed, oracle.build_oracle())


if __name__ == "__main__":
    unittest.main()
