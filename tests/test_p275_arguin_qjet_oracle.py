from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p275_arguin_qjet_oracle import build_oracle  # noqa: E402


class ArguinQJetOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_oracle(10)

    def test_all_stirling_conversions_give_pascal_rows(self) -> None:
        self.assertTrue(self.result["all_stirling_checks_pass"])
        for row in self.result["orders"]:
            self.assertEqual(
                row["coefficients_from_ordinary_identity_via_stirling"],
                row["log_Q_identity"]["coefficients"],
            )

    def test_all_polynomial_basis_checks_pass(self) -> None:
        self.assertTrue(self.result["all_polynomial_basis_checks_pass"])

    def test_first_four_logarithmic_rows(self) -> None:
        rows = self.result["orders"]
        self.assertEqual(rows[0]["log_Q_identity"]["coefficients"], [1])
        self.assertEqual(rows[1]["log_Q_identity"]["coefficients"], [1, 2])
        self.assertEqual(rows[2]["log_Q_identity"]["coefficients"], [1, 3, 3])
        self.assertEqual(rows[3]["log_Q_identity"]["coefficients"], [1, 4, 6, 4])

    def test_third_order_basis_warning_is_explicit(self) -> None:
        warning = self.result["basis_warning"]
        self.assertIn("6 partial_Q P", warning["ordinary_third_order"])
        self.assertIn("3 D P", warning["logarithmic_third_order"])
        self.assertIn("not the D-derivative basis", warning["do_not_mix"])

    def test_residual_is_the_only_new_typed_information(self) -> None:
        consequence = self.result["research_consequence"]
        self.assertIn("R_n", consequence["residual_definition"])
        self.assertIn("Only a typed nonzero residual", consequence["interpretation"])

    def test_committed_artifact_is_reproducible(self) -> None:
        committed = json.loads(
            (ROOT / "results" / "arguin-qjet" / "latest.json").read_text()
        )
        self.assertEqual(committed, self.result)


if __name__ == "__main__":
    unittest.main()
