from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from production_pair_motif_covariance_certificate import (  # noqa: E402
    CONTROL_ORDER,
    DECLARED_PAIRS,
    build_artifact,
    determinant,
    matrix_rank,
)


class ProductionPairMotifCovarianceCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = build_artifact()

    def test_declared_pair_inventory(self) -> None:
        self.assertEqual(len(self.artifact["declared_pairs"]), len(DECLARED_PAIRS))
        self.assertEqual(self.artifact["control_order"], list(CONTROL_ORDER))
        self.assertEqual([item["N"] for item in self.artifact["declared_pairs"]], [65, 85, 130, 145, 170])

    def test_all_exact_gates_pass(self) -> None:
        totals = self.artifact["totals"]
        self.assertEqual(totals["checked_K_values"], 600)
        for key, value in totals.items():
            if key not in {"pairs", "checked_K_values"}:
                self.assertEqual(value, 0, key)

    def test_equal_multiplicities(self) -> None:
        for item in self.artifact["declared_pairs"]:
            self.assertEqual(item["embedding_multiplicities"]["nn_edge"], [2 * item["N"], 2 * item["N"]])
            self.assertEqual(item["embedding_multiplicities"]["diagonal_pair"], [2 * item["N"], 2 * item["N"]])
            self.assertEqual(item["embedding_multiplicities"]["face"], [item["N"], item["N"]])
            self.assertEqual(item["embedding_multiplicities"]["right_angle"], [item["N"], item["N"]])

    def test_representative_covariances_are_nontrivial(self) -> None:
        for item in self.artifact["declared_pairs"]:
            middle = item["representative_rows"][1]
            self.assertGreaterEqual(middle["rank"], 2)
            self.assertTrue(any(value != "0" for row in middle["covariance"] for value in row))

    def test_exact_linear_algebra_helpers(self) -> None:
        matrix = [[Fraction(2), Fraction(1)], [Fraction(1), Fraction(2)]]
        self.assertEqual(determinant(matrix), 3)
        self.assertEqual(matrix_rank(matrix), 2)
        self.assertEqual(determinant([[Fraction(1), Fraction(1)], [Fraction(1), Fraction(1)]]), 0)
        self.assertEqual(matrix_rank([[Fraction(1), Fraction(1)], [Fraction(1), Fraction(1)]]), 1)

    def test_checked_artifact_reproduces(self) -> None:
        path = ROOT / "results" / "production-pair-motif-covariance-certificate.json"
        self.assertEqual(json.loads(path.read_text(encoding="utf-8")), self.artifact)


if __name__ == "__main__":
    unittest.main()

