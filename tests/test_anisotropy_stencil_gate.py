from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import anisotropy_stencil_gate as gate  # noqa: E402


class AnisotropyStencilGateTest(unittest.TestCase):
    def test_spin4_integer_orbits_are_exact(self) -> None:
        self.assertEqual(gate.spin4(1, 0), (1, 0))
        self.assertEqual(gate.spin4(0, 1), (1, 0))
        self.assertEqual(gate.spin4(1, 1), (-4, 0))
        self.assertEqual(gate.spin4(1, -1), (-4, 0))
        self.assertEqual(gate.spin4(2, 1), (-7, 24))

    def test_probability_weighted_critical_axis_family_is_constant(self) -> None:
        for denominator in range(1, 17):
            for numerator in range(denominator + 1):
                t = Fraction(numerator, denominator)
                self.assertEqual(gate.axis_probability_moment(t), (Fraction(2), Fraction(0)))

    def test_variance_weighted_axis_family_has_only_degenerate_zeros(self) -> None:
        self.assertEqual(gate.axis_variance_moment(Fraction(0)), (Fraction(0), Fraction(0)))
        self.assertEqual(gate.axis_variance_moment(Fraction(1)), (Fraction(0), Fraction(0)))
        for denominator in range(2, 17):
            for numerator in range(1, denominator):
                real, imag = gate.axis_variance_moment(Fraction(numerator, denominator))
                self.assertGreater(real, 0)
                self.assertEqual(imag, 0)

    def test_nonnegative_axis_weights_cannot_cancel(self) -> None:
        for horizontal in range(5):
            for vertical in range(5):
                moment = gate.weighted_moment(
                    [
                        gate.EdgeOrbit(1, 0, Fraction(horizontal), multiplicity=2),
                        gate.EdgeOrbit(0, 1, Fraction(vertical), multiplicity=2),
                    ]
                )
                self.assertEqual(moment, (Fraction(2 * (horizontal + vertical)), Fraction(0)))
                if horizontal or vertical:
                    self.assertGreater(moment[0], 0)

    def test_axis_diagonal_zero_ratio_and_signs(self) -> None:
        self.assertEqual(
            gate.axis_diagonal_moment(Fraction(4), Fraction(1)),
            (Fraction(0), Fraction(0)),
        )
        self.assertGreater(gate.axis_diagonal_moment(Fraction(5), Fraction(1))[0], 0)
        self.assertLess(gate.axis_diagonal_moment(Fraction(3), Fraction(1))[0], 0)

    def test_reflection_orbit_cancels_imaginary_part(self) -> None:
        moment = gate.reflected_orbit_moment(2, 1, Fraction(3, 7))
        self.assertEqual(moment, (Fraction(-12), Fraction(0)))

    def test_checked_in_results_are_reproducible(self) -> None:
        artifact = gate.build_artifact()
        checked_json = json.loads(
            (ROOT / "results/anisotropy-stencil-gate/latest.json").read_text(encoding="utf-8")
        )
        checked_markdown = (
            ROOT / "results/anisotropy-stencil-gate/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(checked_json, artifact)
        self.assertEqual(checked_markdown, gate.render_markdown(artifact))


if __name__ == "__main__":
    unittest.main()
