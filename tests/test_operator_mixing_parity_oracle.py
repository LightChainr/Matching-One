from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import operator_mixing_parity_oracle as oracle  # noqa: E402


class OperatorMixingParityOracleTests(unittest.TestCase):
    def test_matching_parity_alternates_exactly(self) -> None:
        even, odd = oracle.FIELDS
        for order in range(6):
            self.assertEqual(bool(oracle.selection_coefficient(even, "S", order)), order % 2 == 0)
            self.assertEqual(bool(oracle.selection_coefficient(even, "D", order)), order % 2 == 1)
            self.assertEqual(bool(oracle.selection_coefficient(odd, "S", order)), order % 2 == 1)
            self.assertEqual(bool(oracle.selection_coefficient(odd, "D", order)), order % 2 == 0)

    def test_primary_exponents_are_exact(self) -> None:
        artifact = oracle.build_artifact()
        found = {row["channel"]: row["contributions"][0]["N_exponent"]
                 for row in artifact["primary_rows"]}
        self.assertEqual(found, {"P4_S": "-1", "P4_D": "-13/8",
                                 "P4_S_prime": "-5/4", "P4_D_prime": "-5/8"})

    def test_primary_structural_map_is_full_rank(self) -> None:
        artifact = oracle.build_artifact()
        self.assertEqual(artifact["primary_map"]["rank"], 4)
        self.assertEqual(artifact["primary_map"]["matrix"], [
            ["1", "0", "0", "0"], ["0", "0", "1", "0"],
            ["0", "0", "0", "1"], ["0", "1", "0", "0"],
        ])

    def test_fractional_rank_handles_dependent_rows(self) -> None:
        self.assertEqual(oracle.matrix_rank([[Fraction(1), Fraction(2)],
                                             [Fraction(2), Fraction(4)]]), 1)

    def test_checked_in_artifacts_reproduce(self) -> None:
        artifact = oracle.build_artifact()
        checked_json = json.loads((ROOT / "results/operator-mixing-parity/latest.json").read_text())
        checked_md = (ROOT / "results/operator-mixing-parity/latest.md").read_text()
        self.assertEqual(checked_json, artifact)
        self.assertEqual(checked_md, oracle.render_markdown(artifact))


if __name__ == "__main__":
    unittest.main()
