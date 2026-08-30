from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p321_subleading_shape_ratio import fieller_set  # noqa: E402


class P321SubleadingShapeRatioTests(unittest.TestCase):
    def test_precise_denominator_gives_bounded_interval(self) -> None:
        result = fieller_set(2.0, 4.0, 0.01, 0.04, 0.0, alpha=0.05)
        self.assertEqual(result["kind"], "bounded")
        self.assertLess(result["lower"], 2.0)
        self.assertGreater(result["upper"], 2.0)
        self.assertTrue(result["denominator_excludes_zero"])

    def test_weak_denominator_is_not_rendered_as_symmetric_error(self) -> None:
        result = fieller_set(0.1, 1.0, 1.0, 0.1, 0.0, alpha=0.05)
        self.assertIn(result["kind"], {"disjoint_unbounded", "all_real"})
        self.assertFalse(result["denominator_excludes_zero"])

    def test_invalid_alpha_rejected(self) -> None:
        with self.assertRaises(ValueError):
            fieller_set(1.0, 1.0, 1.0, 1.0, 0.0, alpha=0.0)


if __name__ == "__main__":
    unittest.main()
