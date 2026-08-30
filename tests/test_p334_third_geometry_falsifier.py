from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p334_third_geometry_falsifier as falsifier  # noqa: E402


class P334ThirdGeometryFalsifierTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = falsifier.build_certificate(ROOT)

    def test_flux_blind_scan_selects_first_asymmetric_two_orbit_hnf(self) -> None:
        selection = self.payload["selection"]
        self.assertEqual(selection["selected_matrix"], [[7, 2], [0, 1]])
        self.assertEqual(selection["selected_N"], 7)
        self.assertEqual(selection["selected_orbits"], [[[0, 1]], [[1, -3]]])
        self.assertTrue(selection["selection_was_flux_blind"])
        self.assertFalse(selection["gaussian_similarity"])
        selected_row = selection["rows_examined"][-1]
        self.assertFalse(selected_row["quarter_turn_lattice_symmetry"])
        self.assertEqual(selected_row["decision"], "select_first_lexicographic_candidate")

    def test_exact_boundary_census_and_continuity(self) -> None:
        census = self.payload["census"]
        self.assertEqual(census["geometry"]["subset_states"], 128)
        self.assertEqual(census["geometry"]["directed_addition_edges"], 448)
        self.assertEqual(census["direct_rank2_edge_count"], 14)
        self.assertTrue(census["gates"]["coefficientwise_continuity"])
        self.assertEqual(census["gates"]["coefficient_failures"], 0)
        self.assertEqual(
            (census["orbits"]["orbit_0"]["birth_edge_count"],
             census["orbits"]["orbit_0"]["exit_edge_count"]),
            (84, 49),
        )
        self.assertEqual(
            (census["orbits"]["orbit_1"]["birth_edge_count"],
             census["orbits"]["orbit_1"]["exit_edge_count"]),
            (28, 21),
        )

    def test_root_pair_is_close_and_one_root_is_exact(self) -> None:
        score = self.payload["score"]
        self.assertTrue(score["paired_zero_test"]["close_pair_pass"])
        exact = score["root_sets"]["orbit_1"]["roots"][0]
        self.assertEqual(exact["kind"], "exact_rational")
        self.assertEqual(exact["root"], "4/7")
        other = score["root_sets"]["orbit_0"]["roots"][0]
        self.assertLess(Fraction(other["lower"]), Fraction(5928, 10000))
        self.assertGreater(Fraction(other["upper"]), Fraction(5927, 10000))

    def test_frozen_reinforcement_topology_is_falsified(self) -> None:
        score = self.payload["score"]
        self.assertEqual(score["character_geometry"]["chi4_gram_real"], "527/625")
        self.assertTrue(
            score["character_geometry"]["characters_linearly_independent_over_R"]
        )
        self.assertEqual(
            [row["character_contribution_alignment"] for row in score["phase_intervals"]],
            ["reinforce", "cancel", "reinforce"],
        )
        frozen = score["frozen_prediction_score"]
        self.assertTrue(frozen["paired_net_zeros_close"])
        self.assertFalse(frozen["reinforcement_only_between_zeros"])
        self.assertEqual(
            frozen["verdict"],
            "paired_timing_survives_but_reinforcement_topology_is_falsified",
        )
        theorem = self.payload["mechanism_update"]["exact_two_orbit_alignment_theorem"]
        self.assertIn("Gram(chi1,chi2) J1 J2", theorem["formula"])

    def test_checked_in_artifacts_reproduce(self) -> None:
        expected_json = json.loads(
            (ROOT / "results/p334-third-geometry-falsifier/latest.json").read_text(
                encoding="utf-8"
            )
        )
        expected_markdown = (
            ROOT / "results/p334-third-geometry-falsifier/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(expected_json, self.payload)
        self.assertEqual(expected_markdown, falsifier.render_markdown(self.payload) + "\n")


if __name__ == "__main__":
    unittest.main()
