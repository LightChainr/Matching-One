import unittest

from scripts.score_threshold_shape_collapse_production import (
    effective_power_from_cover_ratio,
    projection_scale,
)


class ProductionShapeHelpersTest(unittest.TestCase):
    def test_projection_scale(self):
        self.assertAlmostEqual(projection_scale([1.0, 2.0], [3.0, 6.0]), 3.0)

    def test_cover_power_round_trip(self):
        power = 0.625
        ratio = (1.0 - 5.0 ** (-power)) / (1.0 - 2.0 ** (-power))
        self.assertAlmostEqual(effective_power_from_cover_ratio(ratio), power, places=12)


if __name__ == "__main__":
    unittest.main()
