
from __future__ import annotations
import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import composite_spin4_algebra as algebra  # noqa: E402


class CompositeSpin4AlgebraTest(unittest.TestCase):
    def test_exact_cosine_powers(self) -> None:
        self.assertEqual(algebra.cosine_harmonics(0), {0: Fraction(1)})
        self.assertEqual(algebra.cosine_harmonics(1), {4: Fraction(1)})
        self.assertEqual(
            algebra.cosine_harmonics(2),
            {0: Fraction(1, 2), 8: Fraction(1, 2)},
        )
        self.assertEqual(
            algebra.cosine_harmonics(3),
            {4: Fraction(3, 4), 12: Fraction(1, 4)},
        )

    def test_cosine_expansions_sum_to_one_at_zero_angle(self) -> None:
        for power in range(11):
            self.assertEqual(sum(algebra.cosine_harmonics(power).values()), 1)

    def test_named_selection_rules(self) -> None:
        scalar = algebra.classify_monomial({"T4": 1, "S0": 1})
        one_i4 = algebra.classify_monomial({"T4": 1, "I4": 1})
        cubic = algebra.classify_monomial({"T4": 1, "I4": 2})
        self.assertEqual((scalar["relative_q"], scalar["accelerated_w"]), ("2", "6"))
        self.assertEqual(one_i4["harmonics"], {"H0": "1/2", "H8": "1/2"})
        self.assertFalse(one_i4["has_H4"])
        self.assertEqual(cubic["harmonics"], {"H4": "3/4", "H12": "1/4"})
        self.assertEqual(cubic["H12_over_H4"], "1/3")

    def test_conditional_v4_rows(self) -> None:
        mixed = algebra.classify_monomial({"T4": 1, "I4": 1, "V4": 1})
        doubled = algebra.classify_monomial({"T4": 1, "V4": 2})
        self.assertEqual((mixed["relative_q"], mixed["accelerated_w"]), ("14/3", "26/3"))
        self.assertEqual((doubled["relative_q"], doubled["accelerated_w"]), ("16/3", "28/3"))
        self.assertEqual(mixed["harmonics"], {"H4": "3/4", "H12": "1/4"})

    def test_matching_odd_h4_selection_is_counting_rule(self) -> None:
        for values in algebra.exponent_tuples(6):
            counts = dict(zip(algebra.GENERATOR_ORDER, values))
            if counts["T4"] == 0:
                continue
            row = algebra.classify_monomial(counts)
            expected_parity = -1 if counts["T4"] % 2 else 1
            expected_h4 = (counts["T4"] + counts["I4"] + counts["V4"]) % 2 == 1
            self.assertEqual(row["matching_parity"], expected_parity)
            self.assertEqual(row["has_H4"], expected_h4)

    def test_q3_is_absent_beyond_frozen_enumeration(self) -> None:
        for row in algebra.enumerate_matching_odd_h4(10):
            self.assertNotEqual(row["relative_q"], "3")
        certificate = algebra.q3_diophantine_certificate()
        self.assertIn("no monomial at any degree", certificate["conclusion"])

    def test_q6_is_not_exponent_unique(self) -> None:
        rows = [
            row["name"]
            for row in algebra.enumerate_matching_odd_h4(5)
            if row["relative_q"] == "6"
        ]
        self.assertIn("T4*S0^3", rows)
        self.assertIn("T4*I4^2*S0", rows)

    def test_checked_in_results_are_reproducible(self) -> None:
        artifact = algebra.build_artifact()
        checked_json = json.loads(
            (ROOT / "results/composite-spin4-algebra/latest.json").read_text(encoding="utf-8")
        )
        checked_markdown = (
            ROOT / "results/composite-spin4-algebra/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(checked_json, artifact)
        self.assertEqual(checked_markdown, algebra.render_markdown(artifact))


if __name__ == "__main__":
    unittest.main()
