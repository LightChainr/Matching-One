from __future__ import annotations

import json
import sys
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import paired_motif_microcanonical_oracle as oracle  # noqa: E402
from integer_period_torus import gaussian_integer_torus  # noqa: E402


class PairedMotifMicrocanonicalOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = oracle.build_artifact()

    def test_frozen_motif_sizes_and_multiplicities(self) -> None:
        self.assertEqual(oracle.MOTIF_SIZES, {"nn_edge": 2, "diagonal_pair": 2, "face": 4, "right_angle": 3})
        for a, b, n in ((2, 1, 5), (3, 2, 13)):
            motifs = oracle.motif_embeddings(gaussian_integer_torus(a, b))
            self.assertEqual({name: len(items) for name, items in motifs.items()},
                             {"nn_edge": 2 * n, "diagonal_pair": 2 * n, "face": n, "right_angle": n})

    def test_hypergeometric_expectation_is_exact(self) -> None:
        self.assertEqual(oracle.expected_count(13, 5, 26, 2), Fraction(10, 3))
        self.assertEqual(oracle.expected_count(13, 2, 13, 4), Fraction(0))
        with self.assertRaises(ValueError):
            oracle.expected_count(5, 6, 10, 2)

    def test_both_gaussian_pairs_have_zero_fixed_K_difference(self) -> None:
        pairs = self.artifact["paired_gaussian_certificates"]
        self.assertEqual([item["N"] for item in pairs], [5, 13])
        self.assertTrue(all(item["all_fixed_K_difference_sums_zero"] for item in pairs))
        self.assertTrue(all(max(item["max_absolute_fixed_K_difference_sum"].values()) == 0 for item in pairs))
        self.assertTrue(all(item["configurationwise_nonzero_masks"] == 0 for item in pairs))

    def test_declared_production_pairs_are_nontrivial_but_zero_mean(self) -> None:
        pairs = self.artifact["declared_production_pair_gates"]
        self.assertEqual([item["N"] for item in pairs], [65, 85, 130, 145, 170])
        self.assertTrue(all(item["multiplicities_equal"] for item in pairs))
        self.assertTrue(all(item["all_witness_differences_nonzero"] for item in pairs))
        self.assertTrue(all(item["incremental_witness_agrees"] for item in pairs))
        self.assertTrue(all(set(item["fixed_K_mean_difference"].values()) == {"0 for every K=0,...,N"} for item in pairs))

    def test_direct_and_incremental_counters_agree_everywhere(self) -> None:
        all_summaries = self.artifact["paired_gaussian_certificates"] + self.artifact["independent_geometry_controls"]
        self.assertTrue(all(item["incremental_failures"] == 0 for item in all_summaries))

    def test_axis_and_diamond_controls_match_formula(self) -> None:
        controls = self.artifact["independent_geometry_controls"]
        self.assertEqual([(item["name"], item["N"]) for item in controls], [("axis", 9), ("diamond", 8)])
        self.assertTrue(all(item["formula_failures"] == 0 for item in controls))

    def test_unimodular_basis_histograms_are_invariant(self) -> None:
        checks = self.artifact["unimodular_basis_checks"]
        self.assertEqual(len(checks), 6)
        self.assertTrue(all(item["determinant"] == 1 and item["joint_histograms_equal"] for item in checks))

    def test_checked_in_artifacts_reproduce(self) -> None:
        checked_json = json.loads((ROOT / "results/paired-motif-microcanonical/latest.json").read_text())
        checked_md = (ROOT / "results/paired-motif-microcanonical/latest.md").read_text()
        self.assertEqual(checked_json, self.artifact)
        self.assertEqual(checked_md, oracle.render_markdown(self.artifact))


if __name__ == "__main__":
    unittest.main()
