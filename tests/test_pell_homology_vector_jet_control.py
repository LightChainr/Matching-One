from __future__ import annotations

import sys
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from pell_homology_vector_jet_control import (  # noqa: E402
    PILOT_DESIGNS,
    analyze_design,
    category_key,
    exact_n4_jet,
    run_batches,
)


class ExactVectorJetTests(unittest.TestCase):
    def test_n4_exact_jet_obeys_probability_sum_rules(self) -> None:
        payload = exact_n4_jet()
        self.assertEqual(payload["configurations"], 256)
        self.assertAlmostEqual(payload["sum_rules"]["value"], 1.0)
        self.assertAlmostEqual(payload["sum_rules"]["d_dp"], 0.0)
        self.assertAlmostEqual(payload["sum_rules"]["d2_dp2"], 0.0)
        values = {
            row["category"]: row["jet"]["value"]
            for row in payload["coordinates"]
        }
        self.assertAlmostEqual(values["rank0"], 75 / 256)
        self.assertAlmostEqual(values["rank2"], 75 / 256)
        self.assertAlmostEqual(values["rank1:1,0"], 57 / 256)
        self.assertAlmostEqual(values["rank1:1,-2"], 1 / 256)

    def test_rank_one_categories_retain_full_direction(self) -> None:
        self.assertEqual(category_key("l0", (1, 0)), "rank1:1,0")
        self.assertEqual(category_key("rank1_other", (2, -1)), "rank1:2,-1")


class SmallVectorJetPilotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.rows = run_batches(
            samples=400,
            batches=4,
            seed=19,
            workers=1,
            p0=0.5,
            h=0.01,
        )

    def test_full_vectors_conserve_counts_at_all_three_points(self) -> None:
        self.assertEqual(len(self.rows), 8)
        for identifier, matrix in PILOT_DESIGNS:
            result = analyze_design(identifier, matrix, self.rows, h=0.01, dps=50)
            self.assertEqual(result["count_conservation"], {
                "minus": True,
                "center": True,
                "plus": True,
            })
            self.assertEqual(result["invariant_failures"], {
                "minus": 0,
                "center": 0,
                "plus": 0,
            })
            self.assertIn("rank0", result["support"])
            self.assertIn("rank2", result["support"])
            self.assertIn("rank1:1,0", result["support"])
            self.assertIn("d_dp", result["jet_covariance_of_mean"])
            self.assertIn("d2_dp2", result["three_shortest_line_contrasts"])


if __name__ == "__main__":
    unittest.main()
