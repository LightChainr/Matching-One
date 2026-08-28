from __future__ import annotations

from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_pivotal_russo import (  # noqa: E402
    russo_audit,
    total_pivotal_mass,
    wrapping_event,
)
from integer_period_torus import (  # noqa: E402
    axis_integer_torus,
    diamond_integer_torus,
    gaussian_integer_torus,
)


class ExactPivotalRussoTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.mp.dps = 70

    def test_cross_derivative_matches_pivotal_and_threshold_rank_oracles(self) -> None:
        cases = (
            (axis_integer_torus(2), mp.mpf("0.317")),
            (diamond_integer_torus(2), mp.mpf("0.431")),
            (gaussian_integer_torus(2, 1), mp.mpf("0.619")),
        )
        for geometry, probability in cases:
            with self.subTest(geometry=geometry.name):
                result = russo_audit(geometry, probability)
                self.assertLess(
                    abs(result["analytic_minus_pivotal"]), mp.mpf("1e-60")
                )
                self.assertLess(
                    abs(result["threshold_rank_minus_pivotal"]), mp.mpf("1e-60")
                )
                self.assertGreater(result["primal_total_pivotal_mass"], 0)
                self.assertGreater(
                    result["matching_total_pivotal_mass_at_complement"], 0
                )

    def test_either_channel_russo_identity_is_independent_of_rank_engine(self) -> None:
        result = russo_audit(
            axis_integer_torus(2),
            mp.mpf("0.61"),
            channel="either",
            include_threshold_rank=False,
        )
        self.assertLess(abs(result["analytic_minus_pivotal"]), mp.mpf("1e-60"))
        self.assertNotIn("threshold_rank_derivative", result)

    def test_pivotal_mass_is_nonnegative_for_every_supported_channel(self) -> None:
        geometry = axis_integer_torus(2)
        for channel in ("cross", "either", "both", "direction_0", "direction_1"):
            with self.subTest(channel=channel):
                mass = total_pivotal_mass(
                    geometry,
                    mp.mpf("0.5"),
                    matching=False,
                    channel=channel,
                )
                self.assertGreaterEqual(mass, 0)

    def test_unknown_channel_fails_closed(self) -> None:
        geometry = axis_integer_torus(2)
        with self.assertRaises(ValueError):
            wrapping_event(
                geometry,
                [False] * geometry.n,
                matching=False,
                channel="some_wrap",
            )


if __name__ == "__main__":
    unittest.main()
