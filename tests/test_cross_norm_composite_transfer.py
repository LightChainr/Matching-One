from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import unittest

import mpmath as mp
import yaml


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "predictions"
    / "gaussian_norm2_norm5_residual_transfer_20260828.yaml"
)


def cos4(pair: tuple[int, int]) -> Fraction:
    a, b = pair
    n = a * a + b * b
    return Fraction(a**4 - 6 * a * a * b * b + b**4, n * n)


class CrossNormCompositeTransferTests(unittest.TestCase):
    def test_cubic_sideband_uses_lineage_specific_angular_ratio(self) -> None:
        mp.mp.dps = 80
        with ARTIFACT.open(encoding="utf-8") as handle:
            frozen = yaml.safe_load(handle)

        alpha = mp.mpf(13) / 8
        beta = mp.mpf(2)
        r4_norm5 = -mp.mpf(14) / 25
        r4_norm2 = -mp.mpf(1)
        rf_norm2 = -mp.mpf(1)

        cases = {
            "N65_to_N325": {
                "parent": ((8, 1), (7, 4)),
                "child": ((17, 6), (18, 1)),
                "source_key": "N65",
            },
            "N85_to_N425": {
                "parent": ((9, 2), (7, 6)),
                "child": ((16, 13), (19, 8)),
                "source_key": "N85",
            },
        }

        model = frozen["transfer_laws"]["nonlinear_T4_I4_squared_aligned"]
        for lineage, case in cases.items():
            parent_first, parent_second = case["parent"]
            child_first, child_second = case["child"]
            parent_cube = cos4(parent_first) ** 3 - cos4(parent_second) ** 3
            child_cube = cos4(child_first) ** 3 - cos4(child_second) ** 3
            rf_norm5_fraction = child_cube / parent_cube

            reported_fraction = Fraction(
                model["lineages"][lineage]["norm5_cubic_angular_ratio_exact"]
            )
            self.assertEqual(reported_fraction, rf_norm5_fraction)

            rf_norm5 = (
                mp.mpf(rf_norm5_fraction.numerator)
                / rf_norm5_fraction.denominator
            )
            expected_transfer = (
                mp.power(5, -alpha)
                * (rf_norm5 * mp.power(5, -beta) - r4_norm5)
                / (
                    mp.power(2, -alpha)
                    * (rf_norm2 * mp.power(2, -beta) - r4_norm2)
                )
            )
            reported_transfer = mp.mpf(
                model["lineages"][lineage]["E5_over_E2"]
            )
            self.assertLess(
                abs(reported_transfer - expected_transfer), mp.mpf("1e-48")
            )

            source = frozen["source_norm2_residuals"][case["source_key"]]
            target = model["lineages"][lineage]["prospective_E5_from_source"]
            self.assertLess(
                abs(
                    mp.mpf(target["mean"])
                    - reported_transfer * mp.mpf(source["E2"])
                ),
                mp.mpf("1e-34"),
            )
            self.assertLess(
                abs(
                    mp.mpf(target["source_only_standard_error"])
                    - reported_transfer * mp.mpf(source["E2_standard_error"])
                ),
                mp.mpf("1e-34"),
            )


if __name__ == "__main__":
    unittest.main()
