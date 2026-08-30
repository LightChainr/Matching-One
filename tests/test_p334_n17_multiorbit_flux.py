from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p334_n17_multiorbit_flux as flux  # noqa: E402


class P334N17MultiorbitFluxTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = flux.build_certificate()

    def test_n17_gate_has_two_opposite_character_orbits(self) -> None:
        census = self.payload["n17_census"]
        self.assertEqual(census["geometry"]["N"], 17)
        self.assertEqual(census["geometry"]["subset_states"], 131072)
        self.assertEqual(census["geometry"]["directed_addition_edges"], 1114112)
        self.assertEqual(census["gates"]["orbit_count"], 2)
        self.assertTrue(census["gates"]["characters_are_exact_opposites"])
        self.assertTrue(self.payload["orbit_gate"]["passed"])
        self.assertFalse(self.payload["orbit_gate"]["fallback_scan_used"])

    def test_n17_exact_flow_counts(self) -> None:
        census = self.payload["n17_census"]
        axis = census["orbits"]["axis_orbit"]
        diagonal = census["orbits"]["diagonal_orbit"]
        self.assertEqual(axis["chi4"], {"real": "161/289", "imag": "240/289"})
        self.assertEqual(
            diagonal["chi4"], {"real": "-161/289", "imag": "-240/289"}
        )
        self.assertEqual(
            (axis["birth_edge_count"], axis["exit_edge_count"]),
            (150824, 81600),
        )
        self.assertEqual(
            (diagonal["birth_edge_count"], diagonal["exit_edge_count"]),
            (16218, 9418),
        )
        self.assertEqual(census["direct_rank2_edge_count"], 8823)

    def test_continuity_is_coefficientwise_exact(self) -> None:
        census = self.payload["n17_census"]
        self.assertTrue(census["gates"]["coefficientwise_dA4_equals_birth_minus_exit"])
        self.assertEqual(census["gates"]["coefficient_failures"], 0)
        for row in census["coefficient_rows"]:
            for label in flux.LABELS:
                self.assertTrue(row[label]["coefficient_identity_pass"])

    def test_cross_quotient_shares_reinforce_and_are_close(self) -> None:
        comparison = self.payload["cross_geometry_signed_share"]
        axis13 = float(comparison["n13"]["axis_orbit"])
        axis17 = float(comparison["n17"]["axis_orbit"])
        self.assertTrue(self.payload["gates"]["both_quotients_reinforce_at_p_ref"])
        self.assertGreater(axis17, axis13)
        self.assertLess(abs(axis17 - axis13), 0.01)
        self.assertLess(float(comparison["l1_shift"]), 0.02)

    def test_checked_in_artifacts_reproduce(self) -> None:
        expected_json = json.loads(
            (ROOT / "results/p334-n17-multiorbit-flux/latest.json").read_text(
                encoding="utf-8"
            )
        )
        expected_markdown = (
            ROOT / "results/p334-n17-multiorbit-flux/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(expected_json, self.payload)
        self.assertEqual(expected_markdown, flux.render_markdown(self.payload) + "\n")


if __name__ == "__main__":
    unittest.main()
