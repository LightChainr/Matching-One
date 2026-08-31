from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from p334_twelve_canonical_crossings import root_certificate


class CanonicalCrossingTests(unittest.TestCase):
    def test_one_coefficient_sign_change_forces_one_root(self):
        intervals, _ = root_certificate([3, 1, -2, -4])
        self.assertEqual(len(intervals), 1)

    def test_two_sign_changes_can_survive_or_disappear(self):
        retained, _ = root_certificate([-1, 3, -1])
        smoothed, _ = root_certificate([-10, 1, -10])
        self.assertEqual(len(retained), 2)
        self.assertEqual(len(smoothed), 0)


if __name__ == "__main__":
    unittest.main()
