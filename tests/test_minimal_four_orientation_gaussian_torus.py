
from __future__ import annotations
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from minimal_four_orientation_gaussian_torus import (  # noqa: E402
    minimal_n_for_at_least_four_orientations,
    primitive_first_octant_representations,
    primitive_orientation_count_formula,
)


class MinimalFourOrientationGaussianTorusTests(unittest.TestCase):
    def test_formula_matches_bruteforce_through_1105(self) -> None:
        for n in range(1, 1106):
            self.assertEqual(
                primitive_orientation_count_formula(n),
                len(primitive_first_octant_representations(n)),
                n,
            )

    def test_1105_is_first_four_orientation_size(self) -> None:
        self.assertEqual(minimal_n_for_at_least_four_orientations(), 1105)
        self.assertEqual(
            primitive_first_octant_representations(1105),
            [(33, 4), (32, 9), (31, 12), (24, 23)],
        )


if __name__ == "__main__":
    unittest.main()
