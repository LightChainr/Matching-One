import copy
import json
from pathlib import Path
import unittest

from scripts.degree1_rational_exclusion import (
    DEFAULT_CONTRACT,
    ROOT,
    build_result,
    primitive_degree_one_coefficients,
    validate_result,
)


RESULT = ROOT / "results" / "pslq-degree1-rational-exclusion" / "latest.json"


class DegreeOneRationalExclusionTests(unittest.TestCase):
    def test_primitive_enumeration_is_sign_normalized(self):
        values = primitive_degree_one_coefficients(2)
        self.assertEqual(len(values), len(set(values)))
        self.assertTrue(all(a1 > 0 for _, a1 in values))
        self.assertIn((-1, 2), values)
        self.assertNotIn((-2, 2), values)

    def test_checked_in_result_reproduces_exactly(self):
        checked_in = json.loads(RESULT.read_text(encoding="utf-8"))
        self.assertEqual(checked_in, build_result(DEFAULT_CONTRACT))
        summary = validate_result(checked_in)
        self.assertTrue(summary["all_method_intervals_excluded"])
        self.assertEqual(summary["interval_count"], 4)

    def test_every_interval_has_no_degree_one_root(self):
        result = build_result()
        for interval in result["interval_results"]:
            self.assertTrue(interval["excluded"])
            self.assertEqual(interval["zero_containing_residuals"], 0)
            self.assertGreater(interval["closest_polynomial"]["slope_abs"], 0)

    def test_result_tampering_fails_closed(self):
        result = build_result()
        tampered = copy.deepcopy(result)
        tampered["interval_results"][0]["excluded"] = False
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(tampered)

    def test_contract_digest_tampering_fails_closed(self):
        result = build_result()
        tampered = copy.deepcopy(result)
        tampered["contract"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(tampered)


if __name__ == "__main__":
    unittest.main()
