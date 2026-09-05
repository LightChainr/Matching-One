import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from degree3_interval_exclusion import build_result, output_path, validate_result  # noqa: E402


INTERVAL_ID = "jacobsen-2015-eigenvalue"


class DegreeThreeJacobsenExclusionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result(INTERVAL_ID)

    def test_committed_result_reproduces(self):
        checked = json.loads(output_path(INTERVAL_ID).read_text())
        self.assertEqual(checked, self.result)
        self.assertTrue(validate_result(checked, INTERVAL_ID)["excluded"])

    def test_exact_exclusion_and_closest_sturm_control(self):
        row = self.result["interval_result"]
        self.assertEqual(row["primitive_cubics_checked"], 749_507_743)
        self.assertEqual(row["derivative_stationary_fibers"], 0)
        self.assertEqual(row["root_containing_polynomials"], 0)
        self.assertEqual(row["closest_polynomial"]["coefficients_ascending"], [-63, 40, 65, 79])
        self.assertEqual(
            row["closest_polynomial"]["independent_sturm_certificate"]["sturm_open_root_count_in_method_interval"],
            0,
        )


if __name__ == "__main__":
    unittest.main()
