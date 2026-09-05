import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from degree4_interval_exclusion import build_result, output_path, validate_result  # noqa: E402


INTERVAL_ID = "mertens-2022-p-med"


class DegreeFourMertensMedCensusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result(INTERVAL_ID)

    def test_committed_result_reproduces(self):
        checked = json.loads(output_path(INTERVAL_ID).read_text())
        self.assertEqual(checked, self.result)
        self.assertFalse(validate_result(checked, INTERVAL_ID)["excluded"])

    def test_one_quartic_survives_the_method_interval(self):
        row = self.result["interval_result"]
        self.assertEqual(row["near_candidates_exactly_checked"], 1548)
        self.assertEqual(row["root_filter_candidates"], 3)
        self.assertEqual(row["root_containing_polynomials"], 1)
        self.assertEqual(row["distinct_roots_in_interval"], 1)
        self.assertEqual(row["near_candidates_with_stationary_point"], 0)
        self.assertEqual(row["root_witnesses"][0]["coefficients_ascending"], [-84, 99, -7, 99, 58])


if __name__ == "__main__":
    unittest.main()
