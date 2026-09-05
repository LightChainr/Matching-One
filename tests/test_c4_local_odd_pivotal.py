#!/usr/bin/env python3


from __future__ import annotations
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from c4_local_odd_pivotal import exact_n10_report  # noqa: E402


class C4LocalOddPivotalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = exact_n10_report()

    def test_exact_response_matrix(self) -> None:
        self.assertEqual(
            self.report["response_matrix"],
            [["15/8", "5/4"], ["-3/64", "11/64"]],
        )
        self.assertEqual(self.report["determinant"], "195/512")

    def test_observable_is_nontrivial_and_symmetric(self) -> None:
        self.assertEqual(
            self.report["local_twice_observable_counts"],
            {"-1": 88, "0": 848, "1": 88},
        )
        self.assertFalse(any(self.report["symmetry_violations"].values()))
        self.assertTrue(self.report["passed"])


if __name__ == "__main__":
    unittest.main()
