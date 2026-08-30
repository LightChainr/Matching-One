from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from diagnose_p28_complement_clock_modes import complement_coordinates  # noqa: E402


class P28ComplementClockModeTests(unittest.TestCase):
    def test_exact_even_odd_reconstruction(self):
        row = {"log_density": {"left": [1, 2, 4, 8, 16], "right": [2, 3, 5, 9, 17]}}
        profiles = {
            "K1": {"first": row, "second": row},
            "K2": {"first": row, "second": row},
        }
        coordinates, labels, pairs = complement_coordinates(
            profiles, [2.5, 2.75, 3.0, 3.25, 3.5]
        )
        self.assertEqual(len(coordinates["even"]), 12)
        self.assertEqual(len(coordinates["odd"]), 12)
        self.assertEqual(len(labels), 12)
        self.assertEqual(len(pairs), 4)
        scale = 2**0.5
        for index, pair in enumerate(pairs):
            for mode in range(3):
                even = coordinates["even"][3 * index + mode]
                odd = coordinates["odd"][3 * index + mode]
                self.assertAlmostEqual((even + odd) / scale, pair["matching_residual"][mode])
                self.assertAlmostEqual((even - odd) / scale, pair["primal_residual"][mode])


if __name__ == "__main__":
    unittest.main()
