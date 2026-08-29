from __future__ import annotations

import sys
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_two_clock_generator import exact_oracles, score  # noqa: E402


class TwoClockGeneratorTests(unittest.TestCase):
    def test_three_exact_common_generator_classes(self) -> None:
        cases = exact_oracles()["cases"]
        for name in ("ordinary", "Jordan", "complex_pair"):
            self.assertTrue(cases[name]["A_squared_equals_U"])
            self.assertTrue(cases[name]["commutator_zero"])

    def test_same_spectrum_does_not_align_jordan_direction(self) -> None:
        bad = exact_oracles()["cases"]["same_Jordan_spectrum_wrong_alignment"]
        self.assertFalse(bad["A_squared_equals_U"])
        self.assertFalse(bad["commutator_zero"])

    def test_spin_removal_and_commensurate_power(self) -> None:
        payload = {
            "gaussian_transfer": [[0, -3], [2, 0]],
            "spin_rotation": [[0, -1], [1, 0]],
            "annulus_transfer": [[4, 0], [0, 9]],
            "commensurate_relation": {"gaussian_power": 2, "annulus_power": 1},
        }
        rendered = score(payload)
        self.assertEqual(rendered["verdict"], "compatible_with_one_generator_at_declared_times")
        self.assertEqual(rendered["similarity_invariants"]["gaussian_class"], "distinct_real_or_complex_pair")

    def test_power_gate_rejects_commuting_wrong_clock(self) -> None:
        rendered = score({
            "gaussian_transfer": [[2, 0], [0, 3]],
            "annulus_transfer": [[8, 0], [0, 27]],
            "commensurate_relation": {"gaussian_power": 2, "annulus_power": 1},
        })
        self.assertTrue(rendered["mixed_context_gate"]["passes"])
        self.assertFalse(rendered["commensurate_time_gate"]["passes"])


if __name__ == "__main__":
    unittest.main()
