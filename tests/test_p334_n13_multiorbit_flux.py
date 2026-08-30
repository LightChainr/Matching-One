from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p334_n13_multiorbit_flux as flux  # noqa: E402


class P334N13MultiorbitFluxTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = flux.build_certificate()

    def test_n13_frontier_has_two_opposite_character_orbits(self) -> None:
        census = self.payload["census"]
        self.assertEqual(census["geometry"]["N"], 13)
        self.assertEqual(census["geometry"]["subset_states"], 8192)
        self.assertEqual(census["geometry"]["directed_addition_edges"], 53248)
        self.assertEqual(census["gates"]["orbit_count"], 2)
        self.assertTrue(census["gates"]["characters_are_exact_opposites"])

    def test_birth_exit_identity_is_coefficientwise_exact(self) -> None:
        census = self.payload["census"]
        self.assertTrue(census["gates"]["coefficientwise_dA4_equals_birth_minus_exit"])
        self.assertEqual(census["gates"]["coefficient_failures"], 0)
        for row in census["coefficient_rows"]:
            for label in census["orbits"]:
                self.assertTrue(row[label]["coefficient_identity_pass"])

    def test_both_orbits_have_birth_and_exit_flux(self) -> None:
        for row in self.payload["census"]["orbits"].values():
            self.assertGreater(row["birth_edge_count"], 0)
            self.assertGreater(row["exit_edge_count"], 0)

    def test_reference_contributions_reinforce(self) -> None:
        reference = self.payload["reference_evaluation"]
        self.assertTrue(reference["shares_sum_to_one"])
        self.assertTrue(reference["both_orbits_reinforce_total"])
        shares = [float(value) for value in reference["signed_collinear_share"].values()]
        self.assertAlmostEqual(sum(shares), 1.0, places=15)
        self.assertTrue(all(0.0 < value < 1.0 for value in shares))

    def test_checked_in_artifacts_reproduce(self) -> None:
        expected_json = json.loads(
            (ROOT / "results/p334-n13-multiorbit-flux/latest.json").read_text(
                encoding="utf-8"
            )
        )
        expected_markdown = (
            ROOT / "results/p334-n13-multiorbit-flux/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(expected_json, self.payload)
        self.assertEqual(expected_markdown, flux.render_markdown(self.payload) + "\n")


if __name__ == "__main__":
    unittest.main()
