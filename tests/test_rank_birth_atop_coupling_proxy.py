from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import rank_birth_atop_coupling_proxy as oracle  # noqa: E402


class RankBirthATopCouplingProxyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = oracle.build_artifact()
        cls.rows = {row["id"]: row for row in cls.artifact["tiny_exact"]}

    def test_connected_complement_sign_controls_close(self) -> None:
        for row in self.rows.values():
            self.assertTrue(
                all(value == "0" for value in row["complement_covariance_residuals"].values())
            )

    def test_line_resolved_D_coupling_is_finitely_nonzero(self) -> None:
        expected = {
            "axis-L2-degenerate": ("1", "0"),
            "gaussian-2-1": ("-49/128", "21/16"),
            "axis-L4-fixed-root": ("5013127/4194304", "0"),
        }
        for name, pair in expected.items():
            covariance = self.rows[name]["primal"]["connected_covariance_with_A_top"]
            self.assertEqual((covariance["line_cos4_D"], covariance["line_sin4_D"]), pair)
            self.assertTrue(self.rows[name]["finite_lattice_nonzero"]["line_D_complex"])

    def test_axis_landing_proxy_and_birth_normalization(self) -> None:
        axis = self.rows["axis-L4-fixed-root"]["primal"]
        covariance = axis["connected_covariance_with_A_top"]
        normalized = axis["normalized_D_H4_proxy_over_birth_mass"]
        self.assertEqual(covariance["landing_h4_D"], "6977235/33554432")
        self.assertEqual(normalized["line_cos4_D"], "5013127/17240064")

    def test_scaling_exponents_keep_the_two_normalizations_separate(self) -> None:
        scaling = self.artifact["scaling_if_Q4_epsilon"]
        self.assertEqual(scaling["candidate_dimension_x"], "21/4")
        self.assertEqual(
            scaling["direct_CFT_density_normalization"],
            {
                "one_local_insertion_N_power": "-21/8",
                "sum_over_N_sites_N_power": "-13/8",
                "derivation": "N^(-x/2) locally and N^(1-x/2) after summing N sites",
            },
        )
        rank_birth = scaling["rank_birth_measure_normalization"]
        self.assertEqual(rank_birth["normalized_proxy_C_over_B_N_power"], "-13/8")
        self.assertEqual(rank_birth["raw_connected_sum_C_N_power"], "-5/4")
        self.assertEqual(rank_birth["raw_per_site_contribution_N_power"], "-9/4")

    def test_archive_snapshot_stops_before_the_connected_proxy(self) -> None:
        snapshot = self.artifact["archive_reconstructibility"]["bedc94b_snapshot"]
        self.assertEqual(snapshot["source_N"], 64)
        self.assertEqual(snapshot["source_sample_count"], 100000)
        self.assertIsNotNone(snapshot["available_unmarked_birth_mass_B_equals_M_prime"])
        self.assertIsNotNone(snapshot["available_unmarked_D_equals_f12_minus_f01"])
        self.assertIsNone(snapshot["connected_D_H4_proxy"])

    def test_checked_in_artifacts_reproduce(self) -> None:
        directory = ROOT / "results/rank-birth-atop-coupling"
        checked_json = json.loads((directory / "latest.json").read_text(encoding="utf-8"))
        checked_markdown = (directory / "latest.md").read_text(encoding="utf-8")
        self.assertEqual(checked_json, self.artifact)
        self.assertEqual(checked_markdown, oracle.render_markdown(self.artifact))


if __name__ == "__main__":
    unittest.main()

