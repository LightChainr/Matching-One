from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import projective_essential_birth_oracle as oracle  # noqa: E402


class ProjectiveEssentialBirthOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.certificate = oracle.build_certificate()

    def test_three_minimal_controls_and_exact_path_counts(self) -> None:
        rows = {row["id"]: row for row in self.certificate["geometries"]}
        self.assertEqual(
            set(rows),
            {"axis-L2", "gaussian-2-1", "c4-self-matching-3-1"},
        )
        self.assertEqual(rows["axis-L2"]["all_paths"], 24)
        self.assertEqual(rows["gaussian-2-1"]["all_paths"], 120)
        self.assertEqual(rows["c4-self-matching-3-1"]["all_paths"], 3_628_800)

    def test_projective_mark_and_integral_contract(self) -> None:
        for row in self.certificate["geometries"]:
            self.assertTrue(row["canonicality"]["all_pass"])
            self.assertEqual(row["subgroup_contract"]["iota"], 1)
            self.assertEqual(row["subgroup_contract"]["after_line_birth"], "Z ell with ell primitive")
            self.assertTrue(row["sum_over_ell_plus_direct_atom_recovers_K1"])
            self.assertGreater(row["terminal_marked_state_count"], 0)
            self.assertTrue(row["birth_site_pair_counts"])

    def test_complement_basis_and_D4_covariance(self) -> None:
        for row in self.certificate["geometries"]:
            self.assertTrue(row["complement_Alexander"]["all_pass"])
        for audit in self.certificate["SL2Z_basis_covariance"].values():
            self.assertTrue(audit["all_pass"])
        self.assertTrue(self.certificate["D4_chi4_action"]["all_pass"])
        self.assertEqual(
            {row["chi4_action"] for row in self.certificate["D4_chi4_action"]["elements"]},
            {"identity", "conjugation"},
        )

    def test_direction_is_extra_but_tiny_chi4_is_not(self) -> None:
        for row in self.certificate["geometries"]:
            self.assertTrue(row["ell_not_determined_by_K1_K2"])
            self.assertTrue(row["ell_independent_of_K1_K2_conditional_on_line_birth"])
            self.assertFalse(row["chi4_adds_information_beyond_K1_K2_on_this_quotient"])
        rows = {row["id"]: row for row in self.certificate["geometries"]}
        self.assertGreater(rows["axis-L2"]["direct_rank2_paths_without_projective_line"], 0)
        self.assertEqual(rows["gaussian-2-1"]["direct_rank2_paths_without_projective_line"], 0)
        self.assertGreater(
            rows["c4-self-matching-3-1"]["direct_rank2_paths_without_projective_line"],
            0,
        )

    def test_A4_is_issue156_character_and_flux_split_is_new(self) -> None:
        for row in self.certificate["geometries"]:
            crosswalk = row["primitive_sector_crosswalk"]
            self.assertTrue(crosswalk["all_pass"])
            self.assertEqual(crosswalk["state_crosswalk_path_failures"], 0)
            self.assertEqual(crosswalk["derivative_coefficient_failures"], 0)
        self.assertIn("Issue 156", self.certificate["scientific_crosswalk"]["A4"])

    def test_checked_in_artifacts_reproduce(self) -> None:
        expected_json = json.loads(
            (ROOT / "results/projective-essential-birth/latest.json").read_text(
                encoding="utf-8"
            )
        )
        expected_markdown = (
            ROOT / "results/projective-essential-birth/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(expected_json, self.certificate)
        self.assertEqual(expected_markdown, oracle.render_markdown(self.certificate) + "\n")


if __name__ == "__main__":
    unittest.main()
