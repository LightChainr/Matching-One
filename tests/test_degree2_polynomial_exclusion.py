import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from degree2_polynomial_exclusion import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    primitive_quadratic_count,
    validate_result,
)


class DegreeTwoPolynomialExclusionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_small_primitive_count(self):
        self.assertEqual(primitive_quadratic_count(1), 9)

    def test_checked_in_result_reproduces_exactly(self):
        checked_in = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked_in, self.result)
        summary = validate_result(checked_in)
        self.assertTrue(summary["all_method_intervals_excluded"])
        self.assertEqual(summary["interval_count"], 4)

    def test_every_interval_is_excluded_and_sturm_checked(self):
        for row in self.result["interval_results"]:
            self.assertTrue(row["excluded"])
            self.assertEqual(row["root_containing_polynomials"], 0)
            certificate = row["closest_polynomial"]["independent_sturm_certificate"]
            self.assertEqual(certificate["sturm_open_root_count_in_method_interval"], 0)
            self.assertEqual(certificate["isolation_bits"], 60)

    def test_result_and_boundary_tampering_fail_closed(self):
        altered = copy.deepcopy(self.result)
        altered["interval_results"][0]["excluded"] = False
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(altered)
        boundary = copy.deepcopy(self.result)
        boundary["claim_boundary"]["parent_issue"] = "closed"
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(boundary)


if __name__ == "__main__":
    unittest.main()
