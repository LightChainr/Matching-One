import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from degree4_interval_exclusion import build_result, output_path, primitive_quartic_count, validate_result  # noqa: E402


INTERVAL_ID = "jacobsen-2015-eigenvalue"


class DegreeFourJacobsenExclusionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result(INTERVAL_ID)

    def test_primitive_count(self):
        self.assertEqual(primitive_quartic_count(), 157_309_446_881)

    def test_committed_result_reproduces(self):
        checked = json.loads(output_path(INTERVAL_ID).read_text())
        self.assertEqual(checked, self.result)
        self.assertTrue(validate_result(checked, INTERVAL_ID)["excluded"])

    def test_certified_screen_and_closest_witness(self):
        row = self.result["interval_result"]
        self.assertEqual(row["near_candidates_exactly_checked"], 1543)
        self.assertEqual(row["root_filter_candidates"], 0)
        self.assertEqual(row["root_containing_polynomials"], 0)
        self.assertEqual(row["near_candidates_with_stationary_point"], 0)
        self.assertEqual(row["closest_polynomial"]["coefficients_ascending"], [-84, 99, -7, 99, 58])

    def test_tampering_fails_closed(self):
        changed = copy.deepcopy(self.result)
        changed["interval_result"]["excluded"] = False
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(changed, INTERVAL_ID)


if __name__ == "__main__":
    unittest.main()
