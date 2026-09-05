from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_recurrence_certificate import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    compile_first_recurrence,
    recurrence_candidate,
    solve_linear_system,
    validate_result,
)


class ExactRecurrenceCertificateTests(unittest.TestCase):
    def test_overdetermined_exact_solver_distinguishes_statuses(self) -> None:
        status, solution, rank = solve_linear_system(
            [[Fraction(1), Fraction(2)], [Fraction(2), Fraction(4)]], 1
        )
        self.assertEqual((status, solution, rank), ("supported_unique", [Fraction(2)], 1))
        status, solution, rank = solve_linear_system(
            [[Fraction(1), Fraction(2)], [Fraction(2), Fraction(5)]], 1
        )
        self.assertEqual((status, solution, rank), ("inconsistent", None, 1))

    def test_first_supported_recurrence_rejects_order_one(self) -> None:
        attempts, selected = compile_first_recurrence([Fraction(n) for n in range(1, 6)], 2)
        self.assertEqual([row["status"] for row in attempts], ["inconsistent", "supported_unique"])
        self.assertEqual(selected["coefficients_low_to_high"], ["-1", "2"])
        self.assertEqual(selected["residuals"], ["0", "0", "0"])

    def test_short_sequences_fail_before_overclaiming_uniqueness(self) -> None:
        with self.assertRaisesRegex(ValueError, "insufficient sequence length"):
            recurrence_candidate([Fraction(1), Fraction(2), Fraction(3)], 2)

    def test_checked_in_certificate_reproduces_exactly(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_result())
        summary = validate_result(checked)
        self.assertEqual(summary["selected_order"], 2)
        self.assertEqual(summary["lower_orders_rejected"], [1])
        self.assertEqual(summary["discriminant"], "0")


if __name__ == "__main__":
    unittest.main()
