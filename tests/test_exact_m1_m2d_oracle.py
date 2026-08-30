import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_m1_m2d_oracle import DEFAULT_OUTPUT, build_result, validate_result  # noqa: E402


class ExactM1M2dOracleTests(unittest.TestCase):
    def test_scalar_character_is_excluded_by_nonzero_hankel_minor(self) -> None:
        result = build_result()
        self.assertEqual(result["m1_infeasibility_certificate"]["evaluated_minor"], "1")
        self.assertEqual(result["m1_infeasibility_certificate"]["status"], "exactly_infeasible")

    def test_diagonal_rank_two_realization_reproduces_all_rows(self) -> None:
        realization = build_result()["m2d_extracted_realization"]
        self.assertEqual(realization["eigenvalues"], ["1", "2"])
        self.assertEqual(realization["reproduced_moments"], ["2", "3", "5"])
        self.assertEqual(realization["status"], "exactly_feasible")

    def test_checked_in_certificate_reproduces_exactly(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_result())
        self.assertEqual(validate_result(checked)["m1_minor"], "1")

    def test_tampering_fails_closed(self) -> None:
        tampered = copy.deepcopy(build_result())
        tampered["m1_infeasibility_certificate"]["evaluated_minor"] = "0"
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(tampered)


if __name__ == "__main__":
    unittest.main()
