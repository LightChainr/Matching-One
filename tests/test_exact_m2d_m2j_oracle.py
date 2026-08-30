import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_m2d_m2j_oracle import DEFAULT_OUTPUT, build_result, validate_result  # noqa: E402


class ExactM2dM2jOracleTests(unittest.TestCase):
    def test_minimal_recurrence_has_a_repeated_root_and_rank_two(self) -> None:
        certificate = build_result()["m2d_infeasibility_certificate"]
        self.assertEqual(certificate["nonzero_rank_two_hankel_minor"], "-1")
        self.assertEqual(certificate["recurrence_residuals"], ["0", "0", "0"])
        self.assertEqual(certificate["discriminant"], "0")
        self.assertEqual(certificate["status"], "exactly_infeasible")

    def test_common_nilpotent_jordan_realization_is_exact(self) -> None:
        realization = build_result()["m2j_extracted_realization"]
        self.assertEqual(realization["nilpotent_rank"], 1)
        self.assertEqual(realization["nilpotent_square"], [["0", "0"], ["0", "0"]])
        self.assertEqual(realization["reproduced_moments"], ["1", "2", "3", "4", "5"])

    def test_checked_in_certificate_reproduces_exactly(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_result())
        self.assertEqual(validate_result(checked)["discriminant"], "0")

    def test_tampering_fails_closed(self) -> None:
        tampered = copy.deepcopy(build_result())
        tampered["m2j_extracted_realization"]["nilpotent_rank"] = 0
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(tampered)


if __name__ == "__main__":
    unittest.main()
