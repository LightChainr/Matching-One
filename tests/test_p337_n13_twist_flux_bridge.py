from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p337_n13_twist_flux_bridge import build_certificate, matmul2  # noqa: E402


def complex_fraction(payload: dict[str, str]) -> tuple[Fraction, Fraction]:
    return Fraction(payload["real"]), Fraction(payload["imag"])


class N13TwistFluxBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_certificate()

    def test_fixed_character_bridge_at_every_state_coefficient(self) -> None:
        payload = self.payload
        self.assertTrue(payload["gates"]["all_pass"])
        self.assertTrue(payload["gates"]["quarter_turn_commutes_with_period_matrix"])
        self.assertEqual(complex_fraction(payload["z_axis"]),
                         (Fraction(-119, 169), Fraction(120, 169)))
        z = complex_fraction(payload["z_axis"])
        for row in payload["state_coefficient_rows"]:
            h = Fraction(row["H_F3_orbit_coefficient"])
            a4 = complex_fraction(row["primitive_A4_coefficient"])
            self.assertEqual(a4, (z[0] * h, z[1] * h))
            self.assertEqual(Fraction(row["H_F3_unit_coefficient"]), h / 2)
            self.assertEqual(row["A_axis_odd_numerator"], 0)
            self.assertEqual(row["D_diagonal_odd_numerator"], 0)

    def test_quarter_turn_commutes_with_generic_gaussian_periods(self) -> None:
        rotation = [[0, -1], [1, 0]]
        for a, b in ((1, 0), (3, 2), (8, 1), (7, 4), (5, -3)):
            periods = [[a, -b], [b, a]]
            self.assertEqual(matmul2(rotation, periods), matmul2(periods, rotation))

    def test_source_sink_derivative_bridge_is_coefficientwise(self) -> None:
        z = complex_fraction(self.payload["z_axis"])
        for row in self.payload["flux_coefficient_rows"]:
            self.assertTrue(row["all_flux_bridge_identities_pass"])
            source = Fraction(row["H_twist_source_coefficient"])
            sink = Fraction(row["H_twist_sink_coefficient"])
            derivative = Fraction(row["dH_dp_coefficient"])
            self.assertEqual(derivative, source - sink)
            d_a4 = complex_fraction(row["primitive_dA4_dp_coefficient"])
            self.assertEqual(d_a4, (z[0] * derivative, z[1] * derivative))
            for name in (
                "A_birth_odd_numerator", "D_birth_odd_numerator",
                "A_exit_odd_numerator", "D_exit_odd_numerator",
            ):
                self.assertEqual(row[name], 0)

    def test_reference_evaluation_preserves_exact_bridge(self) -> None:
        reference = self.payload["reference_evaluation"]
        z = complex_fraction(self.payload["z_axis"])
        h = Fraction(reference["H_F3"])
        d_h = Fraction(reference["dH_dp"])
        self.assertEqual(
            complex_fraction(reference["primitive_A4"]),
            (z[0] * h, z[1] * h),
        )
        self.assertEqual(
            complex_fraction(reference["primitive_dA4_dp"]),
            (z[0] * d_h, z[1] * d_h),
        )
        self.assertAlmostEqual(reference["decimal"]["H_F3"], 0.2997725442393157)


if __name__ == "__main__":
    unittest.main()
