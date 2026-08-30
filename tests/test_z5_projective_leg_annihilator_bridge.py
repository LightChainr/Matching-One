from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_z5_projective_leg_annihilator_bridge import (  # noqa: E402
    annihilator_line,
    canonical,
    exact_map_gate,
    projective_residual,
    transformed_basis,
)


class Z5ProjectiveLegAnnihilatorBridgeTests(unittest.TestCase):
    def test_projective_residual_ignores_complex_scale(self) -> None:
        line = np.asarray([1 + 2j, -0.3j, 0.7, 1.2j, -0.9, 0.2 + 0.1j])
        mapped = (2.7 - 1.4j) * line
        self.assertLess(max(abs(value) for value in projective_residual(line, mapped, 0)), 1e-12)

    def test_minimal_covector_is_invariant_under_charge_row_mixing(self) -> None:
        rng = np.random.default_rng(250)
        matrix = rng.normal(size=(12, 6)) + 1j * rng.normal(size=(12, 6))
        matrix[:, -1] = matrix[:, 0] - 2j * matrix[:, 2]
        line, _ = annihilator_line(matrix)
        mixing = rng.normal(size=(12, 12)) + 1j * rng.normal(size=(12, 12))
        while abs(np.linalg.det(mixing)) < 1e-8:
            mixing += np.eye(12)
        mixed_line, _ = annihilator_line(mixing @ matrix)
        pivot = int(np.argmax(np.abs(line)))
        self.assertLess(np.max(np.abs(canonical(line, pivot) - canonical(mixed_line, pivot))), 1e-10)

    def test_basis_family_is_C4_closed(self) -> None:
        bases = [transformed_basis(alexander_reflection=False, rotation_power=value) for value in range(4)]
        self.assertEqual(len({tuple(row) for row in bases}), 4)
        self.assertTrue(all(max(abs(a) + abs(b) for a, b in row) == 2 for row in bases))

    def test_same_parent_hands_are_not_exact_D4_reflections(self) -> None:
        gate = exact_map_gate()
        self.assertTrue(gate["fiber_multipliers_are_inverse"])
        self.assertFalse(gate["same_parent_children_are_exact_D4_reflections"])


if __name__ == "__main__":
    unittest.main()
