
from __future__ import annotations
import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import discrete_holomorphic_spin4_alias_gate as gate  # noqa: E402


class DiscreteHolomorphicSpin4AliasGateTests(unittest.TestCase):
    def test_axial_spin4_phase_is_scalar_character(self) -> None:
        self.assertEqual(
            tuple(gate.spin4_phase(item) for item in gate.c4_orbit(gate.AXIS)),
            ((Fraction(1), Fraction(0)),) * 4,
        )
        self.assertEqual(gate.real_character_rank(gate.AXIS), 1)

    def test_diagonal_spin4_phase_is_negative_scalar_character(self) -> None:
        self.assertEqual(
            tuple(gate.spin4_phase(item) for item in gate.c4_orbit(gate.DIAGONAL)),
            ((Fraction(-1), Fraction(0)),) * 4,
        )
        self.assertEqual(gate.real_character_rank(gate.DIAGONAL), 1)

    def test_naive_single_orbit_moment_cannot_remove_constant_defect(self) -> None:
        axis = gate.orbit_averages(gate.AXIS, (Fraction(7),) * 4)
        diagonal = gate.orbit_averages(gate.DIAGONAL, (Fraction(7),) * 4)
        self.assertEqual(axis["spin4"], (axis["scalar"], Fraction(0)))
        self.assertEqual(diagonal["spin4"], (-diagonal["scalar"], Fraction(0)))

    def test_axis_diagonal_character_matrix_is_invertible(self) -> None:
        artifact = gate.build_artifact()
        separation = artifact["two_orbit_separation"]
        self.assertEqual(separation["response_matrix"], [["1", "1"], ["1", "-1"]])
        self.assertEqual(separation["determinant"], "-2")
        self.assertEqual(separation["rank"], 2)

    def test_invalid_direction_and_orbit_length_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "nonzero"):
            gate.spin4_phase((0, 0))
        with self.assertRaisesRegex(ValueError, "exactly four"):
            gate.orbit_averages(gate.AXIS, (Fraction(1),) * 3)

    def test_checked_in_artifacts_reproduce(self) -> None:
        artifact = gate.build_artifact()
        checked_json = json.loads((ROOT / "results/discrete-holomorphic-spin4-alias/latest.json").read_text())
        checked_md = (ROOT / "results/discrete-holomorphic-spin4-alias/latest.md").read_text()
        self.assertEqual(checked_json, artifact)
        self.assertEqual(checked_md, gate.render_markdown(artifact))


if __name__ == "__main__":
    unittest.main()
