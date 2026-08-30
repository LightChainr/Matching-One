from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p334_two_orbit_exact_atlas as atlas  # noqa: E402


class P334TwoOrbitExactAtlasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = atlas.build_certificate(ROOT)

    def test_geometry_gate_is_exhaustive_and_flux_blind(self) -> None:
        gate = self.payload["geometry_gate"]
        self.assertEqual(gate["maximum_order"], 12)
        self.assertEqual(gate["matrices_scanned"], 119)
        self.assertIn("Only geometry, line support and chi4 enter", gate["gate"])
        matrices = [row["geometry"]["period_matrix"] for row in self.payload["atlas"]]
        self.assertEqual(
            matrices,
            [
                [[7, 2], [0, 1]],
                [[7, 3], [0, 1]],
                [[7, 4], [0, 1]],
                [[7, 5], [0, 1]],
                [[8, 3], [0, 1]],
                [[8, 5], [0, 1]],
            ],
        )

    def test_root_pair_is_universal_but_frozen_closeness_is_not(self) -> None:
        summary = self.payload["summary"]
        self.assertEqual(summary["included_hnf_count"], 6)
        self.assertEqual(summary["one_simple_root_per_orbit_count"], 6)
        self.assertIsNone(summary["minimal_root_count_counterexample"])
        self.assertEqual(summary["close_pair_count"], 4)
        self.assertEqual(summary["minimal_close_pair_counterexample"], "hnf-8-3-0-1")
        self.assertEqual(summary["exact_mechanism_signature_count"], 2)
        self.assertEqual(summary["total_subset_states"], 1024)
        self.assertEqual(summary["total_directed_boundary_edges"], 3840)

    def test_gram_sign_exactly_stratifies_cooperation_topology(self) -> None:
        for row in self.payload["atlas"]:
            topology = [interval["topology"] for interval in row["phase_intervals"]]
            if row["character_gram_sign"] == "positive":
                self.assertEqual(row["character_gram"], "527/625")
                self.assertEqual(topology, ["reinforce", "cancel", "reinforce"])
            else:
                self.assertEqual(row["character_gram"], "-7/25")
                self.assertEqual(topology, ["cancel", "reinforce", "cancel"])
            self.assertTrue(row["cross_term_polynomial_identity"]["pass"])
            self.assertEqual(
                row["cross_term_polynomial_identity"]["residual_coefficients"], ["0"]
            )

    def test_multiline_orbits_compress_coefficientwise(self) -> None:
        multiline = 0
        for row in self.payload["atlas"]:
            self.assertTrue(row["coefficientwise_continuity"])
            self.assertTrue(row["orbit_compression_gate"]["all_pass"])
            for orbit in row["orbit_compression_gate"]["orbits"]:
                self.assertTrue(orbit["line_boundary_coefficients_equal"])
                self.assertTrue(orbit["effective_character_average_exact"])
                multiline += orbit["line_count"] > 1
        self.assertEqual(multiline, 2)

    def test_checked_in_artifacts_reproduce(self) -> None:
        expected_json = json.loads(
            (ROOT / "results/p334-two-orbit-exact-atlas/latest.json").read_text(
                encoding="utf-8"
            )
        )
        expected_markdown = (
            ROOT / "results/p334-two-orbit-exact-atlas/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(expected_json, self.payload)
        self.assertEqual(expected_markdown, atlas.render_markdown(self.payload) + "\n")


if __name__ == "__main__":
    unittest.main()
