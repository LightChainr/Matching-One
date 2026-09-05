#!/usr/bin/env python3
"""Regression tests for the A-lattice complexity certification.

The wrong number these tests exist to stop us believing is the *height* of the
historical complexity range.  The manuscript's §6.6 censuses `C(<=6, <=4)`
because this polynomial has height 4; if the height here drifts to 3, the paper
silently reverts to a class that does not cover the record, and says it does.
"""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from ziff_a_lattice_complexity import (  # noqa: E402
    DEFAULT_OUTPUT,
    POLYNOMIAL,
    build_result,
    exact_divide,
    monic_integer_factorizations,
    validate_result,
)


def evaluate(coefficients, point: Fraction) -> Fraction:
    total = Fraction(0)
    for coefficient in reversed(coefficients):
        total = total * point + coefficient
    return total


class ZiffALatticeComplexityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_result()
        cls.committed = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))

    def test_committed_artifact_reproduces_exactly(self) -> None:
        self.assertEqual(self.committed, self.result)
        validate_result(self.committed)

    def test_the_height_is_four_and_the_degree_is_five(self) -> None:
        self.assertEqual(self.result["height"], 4)
        self.assertEqual(self.result["degree"], 5)

    def test_the_bracket_really_brackets_a_root(self) -> None:
        """Sign change at the certified endpoints, in exact arithmetic.

        Isolation code that returned a plausible but wrong interval would still
        pass the decimal cross-check to six places.  This does not depend on
        the isolator being right about anything except the endpoints.
        """
        low, high = (Fraction(text) for text in self.result["root_in_unit_interval"]["isolating_bracket"])
        self.assertLess(low, high)
        below, above = evaluate(POLYNOMIAL, low), evaluate(POLYNOMIAL, high)
        self.assertNotEqual(below, 0)
        self.assertNotEqual(above, 0)
        self.assertLess(below * above, 0)

    def test_the_root_agrees_with_the_independently_quoted_decimal(self) -> None:
        """The two halves of the sourcing must corroborate each other.

        The polynomial and the decimal `0.625457` came from separate summaries
        of the same paper, neither of which we could read.  Their agreement is
        the whole of the evidence that either was transcribed correctly, so it
        is checked here rather than left in a prose note.
        """
        check = self.result["cross_check_against_the_quoted_decimal"]
        self.assertTrue(check["agrees_with_the_isolated_root"])
        low, high = (Fraction(text) for text in self.result["root_in_unit_interval"]["isolating_bracket"])
        window_low, window_high = (Fraction(text) for text in check["quoted_window"])
        self.assertLessEqual(window_low, high)
        self.assertLessEqual(low, window_high)

    def test_the_factorization_search_finds_a_factorization_when_one_exists(self) -> None:
        """A positive control on the irreducibility claim.

        `irreducible_over_Q` is reported as "the search found nothing", which is
        also what a broken search reports.  Feed it a monic quintic that does
        factor, and require it to say so.
        """
        reducible = [1, 0, 0, 1, -1, 1]  # (x^2 + 1)(x^3 - x^2 + 1), ascending
        self.assertEqual(evaluate(reducible, Fraction(2)), Fraction(1 + 4) * Fraction(8 - 4 + 1))
        found = monic_integer_factorizations(reducible)
        self.assertTrue(found)
        for entry in found:
            product = [0] * (len(entry["factor"]) + len(entry["cofactor"]) - 1)
            for i, a in enumerate(entry["factor"]):
                for j, b in enumerate(entry["cofactor"]):
                    product[i + j] += a * b
            self.assertEqual(product, reducible)

    def test_exact_divide_rejects_a_division_with_a_remainder(self) -> None:
        self.assertIsNone(exact_divide(POLYNOMIAL, [1, 1]))
        self.assertIsNone(exact_divide(POLYNOMIAL, [-1, 1]))

    def test_the_sourcing_status_is_not_quietly_upgraded(self) -> None:
        """This row sets the height bound of a published claim.

        It is corroborated, not read.  If someone changes the status to a
        primary reading, the note saying what is still owed has to go with it,
        and this test is where that trade is enforced.
        """
        self.assertEqual(self.result["verification_status"], "CORROBORATED_NOT_PRIMARY")
        self.assertIn("primary text", self.result["what_is_still_owed"])


if __name__ == "__main__":
    unittest.main()
