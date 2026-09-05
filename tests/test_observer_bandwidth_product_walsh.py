
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from observer_bandwidth_product_walsh import (  # noqa: E402
    apply_noise,
    build_report,
    centered_basis,
    direct_covariance,
    evaluate_multilinear,
    popcount,
    source_fixture,
    walsh_degree_coefficients,
)


class ProductWalshBandwidthTests(unittest.TestCase):
    def test_centered_basis_has_exact_noise_eigenvalue(self) -> None:
        n = 4
        p = Fraction(2, 5)
        rho = Fraction(3, 7)
        for subset in range(1 << n):
            values = {
                mask: centered_basis(mask, subset, p) for mask in range(1 << n)
            }
            transformed = apply_noise(values, n, p, rho)
            for mask in range(1 << n):
                self.assertEqual(
                    transformed[mask],
                    rho ** popcount(subset) * values[mask],
                )

    def test_degree_two_observer_has_no_higher_covariance_modes(self) -> None:
        n = 4
        p = Fraction(1, 3)
        observer = evaluate_multilinear(
            {0: Fraction(5), 1 << 0: Fraction(2), (1 << 1) | (1 << 3): Fraction(-3)},
            n,
        )
        source = source_fixture(n)
        coefficients = walsh_degree_coefficients(observer, source, n, p)
        self.assertEqual(coefficients[3], 0)
        self.assertEqual(coefficients[4], 0)

    def test_direct_covariance_equals_degree_polynomial(self) -> None:
        n = 4
        p = Fraction(2, 5)
        rho = Fraction(5, 8)
        observer = evaluate_multilinear(
            {1 << 0: Fraction(1), (1 << 1) | (1 << 2): Fraction(4)}, n
        )
        source = source_fixture(n)
        coefficients = walsh_degree_coefficients(observer, source, n, p)
        self.assertEqual(
            direct_covariance(observer, source, n, p, rho),
            sum(rho**degree * value for degree, value in coefficients.items()),
        )

    def test_checked_report_is_exactly_reproducible(self) -> None:
        manifest = json.loads(
            (ROOT / "analysis/observer_bandwidth_product_walsh_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        checked = json.loads(
            (ROOT / "results/observer-bandwidth-product-walsh/latest.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(checked, build_report(manifest))
        self.assertEqual([row["residual"] for row in checked["comparisons"]], ["0"] * 3)
        self.assertLessEqual(max(checked["active_nonconstant_degrees"]), 4)


if __name__ == "__main__":
    unittest.main()
