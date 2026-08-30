from __future__ import annotations

from pathlib import Path
import sys
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_z5_projective_leg_annihilator_bridge import transformed_basis  # noqa: E402
from score_z5_projective_leg_radius5_morphism import (  # noqa: E402
    ROWS_3,
    extension_matrix,
    transformed_rows,
)


class Z5ProjectiveLegRadius5MorphismScoreTests(unittest.TestCase):
    def test_old_only_extension_matrix_has_two_dimensional_empty_append(self) -> None:
        basis = transformed_basis(alexander_reflection=False, rotation_power=0)
        with patch("score_z5_projective_leg_radius5_morphism.old_pair", return_value=1.0 + 0.0j):
            matrix = extension_matrix({}, {}, "plus", basis, ())
        self.assertEqual(matrix.shape, (2 * len(basis), len(basis)))

    def test_transformed_degree_three_rows_stay_on_degree_three(self) -> None:
        for alexander in (False, True):
            for power in range(4):
                rows = transformed_rows(alexander, power)
                self.assertEqual(len(set(rows)), len(ROWS_3))
                self.assertTrue(all(abs(a) + abs(b) == 3 for a, b in rows))


if __name__ == "__main__":
    unittest.main()
