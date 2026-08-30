from __future__ import annotations

import cmath
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_z5_projective_leg_pair_transfer import channel_metrics  # noqa: E402
from z5_projective_leg_pair_transfer_mc import FIELD_ORDER, exact_gate, run  # noqa: E402


def synthetic_values(function) -> dict[str, float]:
    output = {}
    for hand in ("plus", "minus"):
        for charge in (1, 2):
            for separation in (1, 2, 3):
                value = function(separation)
                output[f"d{separation}_T{charge}_{hand}_re"] = value.real
                output[f"d{separation}_T{charge}_{hand}_im"] = value.imag
    return output


class Z5ProjectiveLegPairTransferTests(unittest.TestCase):
    def test_exact_gate_covers_six_separations(self) -> None:
        gate = exact_gate()
        self.assertTrue(gate["passed"])
        self.assertEqual(gate["separations"], [1, 2, 3, 4, 5, 6])
        self.assertEqual(gate["coordinates"], 96)
        self.assertEqual(gate["anchor_collapse_failures"], 0)

    def test_tiny_stream_has_full_96_coordinate_covariance(self) -> None:
        rows, analysis = run(4, 2, 1, 0.59274605079, 25025033720260831, 0)
        self.assertEqual(len(rows), 2)
        self.assertEqual(len(FIELD_ORDER), 96)
        self.assertEqual(len(analysis["covariance_of_mean"]), 96)

    def test_exponential_invariant_is_exact(self) -> None:
        lam = 0.7 * cmath.exp(0.2j)
        values = synthetic_values(lambda separation: (1.3 - 0.4j) * lam**separation)
        row = channel_metrics(values, "plus", 1)
        self.assertAlmostEqual(row["exponential_mass_delta"], 0.0, places=13)
        self.assertAlmostEqual(row["phase_step_delta"], 0.0, places=13)
        self.assertGreater(abs(row["power_eta_delta"]), 0.1)

    def test_power_invariant_is_exact(self) -> None:
        values = synthetic_values(
            lambda separation: (1.3 - 0.4j) * separation**-2.0 * cmath.exp(0.2j * separation)
        )
        row = channel_metrics(values, "minus", 2)
        self.assertAlmostEqual(row["power_eta_delta"], 0.0, places=13)
        self.assertAlmostEqual(row["phase_step_delta"], 0.0, places=13)
        self.assertGreater(abs(row["exponential_mass_delta"]), 0.1)


if __name__ == "__main__":
    unittest.main()
