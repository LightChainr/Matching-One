
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from observer_bandwidth_high_pass import (  # noqa: E402
    apply_high_pass,
    build_report,
    filter_coefficients,
    raw_monomial_values,
    spectral_multiplier,
)
from observer_bandwidth_product_walsh import centered_basis  # noqa: E402


class ObserverBandwidthHighPassTests(unittest.TestCase):
    def test_six_level_polynomial_has_exact_low_degree_roots(self) -> None:
        rho = Fraction(1, 2)
        coefficients = filter_coefficients(rho, 4)
        self.assertEqual(len(coefficients), 6)
        for degree in range(5):
            self.assertEqual(spectral_multiplier(rho, 4, degree), 0)

    def test_raw_degree_four_monomial_is_annihilated(self) -> None:
        n = 5
        filtered = apply_high_pass(
            raw_monomial_values(n, 0b01111),
            n,
            Fraction(2, 5),
            Fraction(1, 2),
            4,
        )
        self.assertTrue(all(value == 0 for value in filtered.values()))

    def test_centered_degree_five_control_survives(self) -> None:
        n = 5
        p = Fraction(2, 5)
        rho = Fraction(1, 2)
        values = {mask: centered_basis(mask, 0b11111, p) for mask in range(1 << n)}
        filtered = apply_high_pass(values, n, p, rho, 4)
        multiplier = spectral_multiplier(rho, 4, 5)
        self.assertGreater(multiplier, 0)
        self.assertTrue(any(filtered.values()))
        self.assertEqual(
            filtered,
            {mask: multiplier * value for mask, value in values.items()},
        )

    def test_checked_report_is_exactly_reproducible(self) -> None:
        manifest = json.loads(
            (ROOT / "analysis/observer_bandwidth_high_pass_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        checked = json.loads(
            (ROOT / "results/observer-bandwidth-high-pass/latest.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(checked, build_report(manifest))
        self.assertEqual(checked["euler_control"]["filtered_l2_energy"], "0")
        self.assertEqual(checked["degree_five_control"]["nonzero_output_points"], 32)


if __name__ == "__main__":
    unittest.main()
