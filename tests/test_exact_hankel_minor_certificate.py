from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_hankel_minor_certificate import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    compile_minors,
    determinant,
    hankel_matrix,
    validate_result,
)


class ExactHankelMinorCertificateTests(unittest.TestCase):
    def test_fraction_determinant_and_hankel_shape(self) -> None:
        matrix = hankel_matrix([Fraction(1, 2), Fraction(2, 3), Fraction(3, 4)], 2, 2)
        self.assertEqual(matrix, [[Fraction(1, 2), Fraction(2, 3)], [Fraction(2, 3), Fraction(3, 4)]])
        self.assertEqual(determinant(matrix), Fraction(-5, 72))

    def test_compiler_enumerates_every_declared_minor(self) -> None:
        matrix = hankel_matrix([Fraction(n) for n in range(1, 6)], 2, 4)
        records = compile_minors(matrix, 2)
        self.assertEqual(len(records), 6)
        self.assertEqual([row["determinant"] for row in records], ["-1", "-2", "-3", "-1", "-2", "-1"])

    def test_checked_in_certificate_reproduces_exactly(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_result())
        summary = validate_result(checked)
        self.assertEqual(summary["minor_count"], 6)
        self.assertEqual(summary["rank_lower_bound"], 2)

    def test_shape_and_order_errors_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "insufficient moments"):
            hankel_matrix([Fraction(1)], 2, 2)
        with self.assertRaisesRegex(ValueError, "outside matrix shape"):
            compile_minors([[Fraction(1), Fraction(2)]], 2)


if __name__ == "__main__":
    unittest.main()
