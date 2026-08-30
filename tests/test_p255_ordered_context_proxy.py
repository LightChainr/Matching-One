from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "score_p255_ordered_context_proxy",
    ROOT / "scripts" / "score_p255_ordered_context_proxy.py",
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class OrderedContextProxyTest(unittest.TestCase):
    def test_psd_score_and_even_survival(self) -> None:
        covariance = np.diag([1.0, 2.0, 3.0, 4.0])
        inverse, rank = MODULE.psd_inverse(covariance)
        self.assertEqual(rank, 4)
        self.assertTrue(np.allclose(inverse, np.diag([1, 1 / 2, 1 / 3, 1 / 4])))
        self.assertAlmostEqual(MODULE.chi2_survival(0.0, 4), 1.0)
        self.assertAlmostEqual(MODULE.chi2_survival(0.0, 3), 1.0)

    def test_freeze_names(self) -> None:
        self.assertEqual(MODULE.CHANNELS, ("ES", "ED", "OS", "OD"))
        self.assertEqual(MODULE.COUNTER, (18000000000, 18000020000))
        self.assertEqual(MODULE.SAMPLES, 20000)


if __name__ == "__main__":
    unittest.main()
