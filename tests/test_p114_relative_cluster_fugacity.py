from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p114_relative_cluster_fugacity import build_oracle  # noqa: E402


class RelativeClusterFugacityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_oracle()

    def test_tiny_source_factorizes_configurationwise(self) -> None:
        tiny = self.result["tiny_exact_oracle"]
        self.assertEqual(tiny["N"], 5)
        self.assertEqual(tiny["configuration_count"], 32)
        self.assertEqual(tiny["factorization_failures"], 0)
        self.assertEqual(tiny["topology_channel_failures"], 0)

    def test_sector_derivative_is_matching_charge(self) -> None:
        tiny = self.result["tiny_exact_oracle"]
        sectors = tiny["sector_Bernstein_coefficients_by_q"]
        expected = [plus - minus for plus, minus in zip(sectors["1"], sectors["-1"])]
        self.assertEqual(tiny["matching_Bernstein_coefficients"], expected)
        self.assertEqual(tiny["p_half"]["first_logQ_derivative"]["text"], "-5/16")

    def test_one_colour_plaquette_local_class_has_exact_counterexample(self) -> None:
        witness = self.result["ordinary_FK_local_obstruction"]
        self.assertTrue(witness["same_black_cluster_count"])
        self.assertTrue(witness["same_complete_plaquette_histogram"])
        self.assertTrue(witness["different_matching_charge"])
        rows = witness["configurations"]
        self.assertEqual([row["q_matching"] for row in rows], [0, -1])
        self.assertEqual([row["white_matching_clusters"] for row in rows], [1, 2])

    def test_three_source_terms_are_not_collapsed(self) -> None:
        terms = self.result["derivative_ledger"]["fixed_p_relative_Q_direction"]
        self.assertEqual(len(terms), 3)
        self.assertTrue(any("white" in term for term in terms))
        self.assertTrue(any("Euler" in term for term in terms))

    def test_committed_artifact_is_reproducible(self) -> None:
        committed = json.loads(
            (ROOT / "results" / "exact-relative-cluster-fugacity" / "latest.json").read_text()
        )
        self.assertEqual(committed, self.result)


if __name__ == "__main__":
    unittest.main()
