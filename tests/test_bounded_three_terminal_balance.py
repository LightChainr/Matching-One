from __future__ import annotations

import json
from pathlib import Path
import unittest

from scripts.score_bounded_three_terminal_balance import (
    DEFAULT_OUTPUT,
    build_artifact,
    connected_orbit_catalog,
    evaluate_power,
    validate_artifact,
)


class BoundedThreeTerminalBalanceTests(unittest.TestCase):
    def test_complete_connected_candidate_space(self) -> None:
        catalog = connected_orbit_catalog()
        self.assertEqual(len(catalog), 11)
        self.assertEqual(sum(1 for _, graph in catalog.values() if sum(3 in edge for edge in graph) >= 3), 4)

    def test_star_polynomial_and_known_root(self) -> None:
        artifact = build_artifact()
        star = next(
            row
            for row in artifact["ranking"]
            if sorted(row["edges"]) == [[0, 3], [1, 3], [2, 3]]
        )
        self.assertEqual(star["primitive_balance_power_coefficients_low_to_high"], [-1, 0, 3, -1])
        self.assertAlmostEqual(float(star["root_decimal"]), 0.6527036446661393, places=15)

    def test_every_isolation_brackets_a_sign_change(self) -> None:
        artifact = build_artifact()
        for row in artifact["ranking"]:
            numerator, denominator = map(int, row["root_isolation"]["lower"].split("/"))
            from fractions import Fraction

            lower = Fraction(numerator, denominator)
            numerator, denominator = map(int, row["root_isolation"]["upper"].split("/"))
            upper = Fraction(numerator, denominator)
            polynomial = row["primitive_balance_power_coefficients_low_to_high"]
            self.assertLess(evaluate_power(polynomial, lower), 0)
            self.assertGreaterEqual(evaluate_power(polynomial, upper), 0)

    def test_checked_artifact_reproduces(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_artifact())
        self.assertEqual(
            validate_artifact(checked)["status"],
            "valid_exact_bounded_three_terminal_balance_screen",
        )

    def test_claim_boundary_stays_open(self) -> None:
        artifact = build_artifact()
        self.assertEqual(artifact["decision"]["exact_square_site_claims_certified"], 0)
        self.assertEqual(artifact["claim_boundary"]["parent_issue"], "remain open")


if __name__ == "__main__":
    unittest.main()
