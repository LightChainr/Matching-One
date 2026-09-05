import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_linear_ideal_certificate import (  # noqa: E402
    DEFAULT_OUTPUT,
    add,
    build_result,
    multiply,
    validate_result,
    verify_ideal_witness,
)


class ExactLinearIdealCertificateTests(unittest.TestCase):
    def test_polynomial_arithmetic_is_low_to_high_and_exact(self) -> None:
        self.assertEqual(add([1, 2], [-1, 0, 3]), [0, 2, 3])
        self.assertEqual(multiply([1, 1], [1, -1]), [1, 0, -1])

    def test_generic_polynomial_multiplier_witness(self) -> None:
        verified = verify_ideal_witness([["0", "1"], ["1", "-1"]], [["1"], ["1"]])
        self.assertEqual(verified["result_coefficients_low_to_high"], ["1"])
        self.assertEqual(verified["status"], "exact_ideal_contains_one")
        self.assertTrue(verified["primitive_after_common_denominator"])

    def test_checked_in_certificate_reproduces_exactly(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_result())
        summary = validate_result(checked)
        self.assertEqual(summary["constraint_count"], 2)
        self.assertEqual(summary["result"], ["1"])
        self.assertIs(summary["solver_invoked"], False)


if __name__ == "__main__":
    unittest.main()
