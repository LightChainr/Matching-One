from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location(
    "p205_etop", ROOT / "scripts" / "p205_quotient_etop_ray_test.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class P205QuotientETopRayTests(unittest.TestCase):
    def test_common_ray_synthetic_rows_have_zero_penalty(self):
        rows = []
        for index, amplitude in enumerate((1, 2, -3)):
            rows.append({"id": str(index), "N": index + 1,
                         "dependency_group": "x",
                         "estimate": [mp.mpf(amplitude), mp.mpf(amplitude) / 2],
                         "covariance": [[mp.mpf(".01"), 0],
                                        [0, mp.mpf(".02")]]})
        fit = MODULE.fit_ray(rows)
        self.assertLess(fit["min_chi2"], 1e-20)

    def test_prism_sources_are_exactly_the_completed_three_sizes(self):
        self.assertEqual(MODULE.SIZES, (25, 50, 125))
        self.assertTrue(MODULE.SOURCE_COMMIT.startswith("fc14817"))

    def test_rotation_pair_is_rejected_in_synthetic_control(self):
        rows = [
            {"id": "a", "N": 1, "dependency_group": "x",
             "estimate": [mp.mpf(1), mp.mpf(0)],
             "covariance": [[mp.mpf("1e-4"), 0],
                            [0, mp.mpf("1e-4")]]},
            {"id": "b", "N": 2, "dependency_group": "x",
             "estimate": [mp.mpf(0), mp.mpf(1)],
             "covariance": [[mp.mpf("1e-4"), 0],
                            [0, mp.mpf("1e-4")]]},
        ]
        self.assertGreater(MODULE.fit_ray(rows)["min_chi2"], 1000)


if __name__ == "__main__":
    unittest.main()
