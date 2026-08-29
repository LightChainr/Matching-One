from __future__ import annotations

from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_p263_rectangle_conformal_secondary import (  # noqa: E402
    LAMBDA_ORDER,
    mapped_geometry,
    rectangle_modulus,
)


class P263RectangleConformalSecondaryTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.mp.dps = 60
        self.modulus = rectangle_modulus()

    def test_modulus_has_frozen_rectangle_aspect_ratio(self) -> None:
        self.assertLess(
            abs(self.modulus["K_prime"] / self.modulus["K"] - mp.mpf(4) / 3),
            mp.mpf("1e-50"),
        )
        self.assertLess(
            abs(
                self.modulus["parameter_m"]
                - mp.mpf(
                    "0.21549970429193269707024331706636106943156632191124"
                )
            ),
            mp.mpf("1e-49"),
        )

    def test_effective_cross_ratios_are_high_precision_regressions(self) -> None:
        expected = (
            "0.23161561099460535818383083600611699355844419221062",
            "0.31280861631947778502067149135200723065244349073853",
            "0.64998175749754660990705108769548529186112608100857",
            "0.73654763309312352113775685892170357505078814310356",
        )
        actual = [
            mapped_geometry(lam, self.modulus)["effective_lambda"]
            for lam in LAMBDA_ORDER
        ]
        for value, target in zip(actual, expected):
            self.assertLess(abs(value - mp.mpf(target)), mp.mpf("1e-49"))

    def test_map_preserves_boundary_order(self) -> None:
        for lam in LAMBDA_ORDER:
            geometry = mapped_geometry(lam, self.modulus)
            images = geometry["uhp_images"]
            self.assertEqual(images, sorted(images))
            self.assertGreater(geometry["mapped_K"], 0)
            self.assertGreater(geometry["effective_lambda"], 0)
            self.assertLess(geometry["effective_lambda"], 1)


if __name__ == "__main__":
    unittest.main()
