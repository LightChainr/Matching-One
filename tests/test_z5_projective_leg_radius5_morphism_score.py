from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_z5_projective_leg_radius5_morphism import ROWS_3, transformed_rows  # noqa: E402


class Z5ProjectiveLegRadius5MorphismScoreTests(unittest.TestCase):
    def test_transformed_degree_three_rows_stay_on_degree_three(self) -> None:
        for alexander in (False, True):
            for power in range(4):
                rows = transformed_rows(alexander, power)
                self.assertEqual(len(set(rows)), len(ROWS_3))
                self.assertTrue(all(abs(a) + abs(b) == 3 for a, b in rows))


if __name__ == "__main__":
    unittest.main()
