from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gaussian_semigroup_design import (  # noqa: E402
    Gaussian,
    default_catalog,
    lineage_payload,
)


class GaussianSemigroupDesignTests(unittest.TestCase):
    def test_exact_harmonics(self) -> None:
        value = Gaussian(8, 1)
        self.assertEqual(value.cos4m(1), Fraction(3713, 4225))
        self.assertEqual(value.cos4m(2), 2 * value.cos4() ** 2 - 1)
        self.assertEqual(
            value.cos4m(3), 4 * value.cos4() ** 3 - 3 * value.cos4()
        )

    def test_doubling_reverses_H4_delta(self) -> None:
        payload = lineage_payload(
            Gaussian(8, 1), Gaussian(7, 4), Gaussian(1, 1)
        )
        prediction = payload["harmonic_predictions"]["H4"]
        self.assertEqual(prediction["parent_delta"]["numerator"], 1152)
        self.assertEqual(prediction["parent_delta"]["denominator"], 845)
        self.assertEqual(prediction["child_delta"]["numerator"], -1152)
        self.assertEqual(prediction["child_delta"]["denominator"], 845)
        self.assertEqual(prediction["angular_ratio"]["numerator"], -1)
        self.assertEqual(prediction["angular_ratio"]["denominator"], 1)
        self.assertEqual(
            payload["child"]["first_canonical"]["pair"], [9, 7]
        )
        self.assertEqual(
            payload["child"]["second_canonical"]["pair"], [11, 3]
        )

    def test_third_doubling_lineage(self) -> None:
        payload = lineage_payload(
            Gaussian(12, 1), Gaussian(9, 8), Gaussian(1, 1)
        )
        self.assertEqual(payload["child"]["first_canonical"]["pair"], [13, 11])
        self.assertEqual(payload["child"]["second_canonical"]["pair"], [17, 1])
        prediction = payload["harmonic_predictions"]["H4"]
        self.assertEqual(prediction["angular_ratio"]["numerator"], -1)
        self.assertEqual(prediction["angular_ratio"]["denominator"], 1)

    def test_norm5_H4_H8_H12_ratios(self) -> None:
        for first, second, multiplier in (
            (Gaussian(8, 1), Gaussian(7, 4), Gaussian(2, -1)),
            (Gaussian(9, 2), Gaussian(7, 6), Gaussian(2, 1)),
        ):
            payload = lineage_payload(first, second, multiplier)
            predictions = payload["harmonic_predictions"]
            self.assertEqual(
                (
                    predictions["H4"]["angular_ratio"]["numerator"],
                    predictions["H4"]["angular_ratio"]["denominator"],
                ),
                (-14, 25),
            )
            self.assertEqual(
                (
                    predictions["H8"]["angular_ratio"]["numerator"],
                    predictions["H8"]["angular_ratio"]["denominator"],
                ),
                (-1054, 625),
            )
            self.assertEqual(
                (
                    predictions["H12"]["angular_ratio"]["numerator"],
                    predictions["H12"]["angular_ratio"]["denominator"],
                ),
                (23506, 15625),
            )

    def test_N1105_catalog_has_four_orientations(self) -> None:
        catalog = default_catalog()
        orientations = set()
        for edge in catalog["N1105_edges"].values():
            orientations.add(tuple(edge["child"]["first_canonical"]["pair"]))
            orientations.add(tuple(edge["child"]["second_canonical"]["pair"]))
            self.assertEqual(edge["child"]["first_canonical"]["N"], 1105)
            self.assertEqual(edge["child"]["second_canonical"]["N"], 1105)
        self.assertEqual(
            orientations,
            {(33, 4), (32, 9), (31, 12), (24, 23)},
        )


if __name__ == "__main__":
    unittest.main()
