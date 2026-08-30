import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_morphism_rank3_oracle import DEFAULT_OUTPUT, build_result, validate_result  # noqa: E402


class ExactMorphismRank3OracleTests(unittest.TestCase):
    def test_endpoint_rows_are_compatible_with_jordan_rank_two(self) -> None:
        endpoint = build_result()["endpoint_m2j_positive_control"]
        self.assertEqual(endpoint["reproduced_moments"], ["1", "2", "3", "4"])
        self.assertEqual(endpoint["status"], "exactly_feasible")

    def test_morphism_row_forces_common_rank_three(self) -> None:
        certificate = build_result()["common_rank_two_infeasibility_certificate"]
        self.assertEqual(certificate["determinant"], "-1")
        self.assertEqual(certificate["certified_minimum_predictive_rank"], 3)
        self.assertEqual(certificate["status"], "exactly_infeasible")

    def test_extracted_rank_three_realization_reproduces_both_row_types(self) -> None:
        realization = build_result()["rank3_extracted_realization"]
        self.assertEqual(realization["reproduced_endpoint_moments"], ["1", "2", "3", "4"])
        self.assertEqual(realization["reproduced_morphism_moments"], ["1", "2", "4"])

    def test_checked_in_certificate_reproduces_and_tampering_fails_closed(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_result())
        self.assertEqual(validate_result(checked)["minimum_rank"], 3)
        tampered = copy.deepcopy(checked)
        tampered["common_rank_two_infeasibility_certificate"]["determinant"] = "0"
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(tampered)


if __name__ == "__main__":
    unittest.main()
