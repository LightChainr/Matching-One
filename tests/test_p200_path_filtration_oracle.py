from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p200_path_filtration_oracle import render  # noqa: E402


class P200PathFiltrationOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = render()

    def test_real_gaussian_hnf_lift(self) -> None:
        geometry = self.payload["geometry"]
        self.assertEqual(geometry["final_column_HNF"], [[10, 3], [0, 1]])
        self.assertEqual(geometry["N2_column_HNF"], [[2, 1], [0, 1]])
        self.assertEqual(geometry["N5_column_HNF"], [[5, 3], [0, 1]])

    def test_first_homology_activation_is_path_ordered(self) -> None:
        activation = self.payload["exhaustive_checks"]["activation"]
        self.assertEqual(activation["black_NN"]["nonzero_configurations"], 75)
        self.assertEqual(activation["white_matching"]["nonzero_configurations"], 50)
        self.assertEqual(
            activation["typed_complement_even"]["mean_at_p_half"]["exact"],
            "-117/2048",
        )
        self.assertEqual(
            activation["typed_complement_odd"]["mean_at_p_half"]["exact"],
            "-33/2048",
        )
        witness = self.payload["witnesses"]["balanced_first_H1_activation"]
        self.assertEqual(witness["mask"], 62)
        self.assertEqual(witness["white_path"]["tau_R2_then_R5"], 2)
        self.assertEqual(witness["white_path"]["tau_R5_then_R2"], 1)
        self.assertEqual(witness["C_typed_even"], "-1/2")
        self.assertEqual(witness["C_typed_odd"], "1/2")

    def test_Doob_contrast_is_nonzero_but_has_exact_zero_mean(self) -> None:
        doob = self.payload["exhaustive_checks"]["Doob_quadratic_variation"]
        self.assertEqual(doob["black_NN"]["nonzero_configurations"], 669)
        self.assertEqual(doob["white_matching"]["nonzero_configurations"], 382)
        self.assertEqual(doob["black_NN"]["mean_at_p_half"]["exact"], "0")
        self.assertEqual(doob["white_matching"]["mean_at_p_half"]["exact"], "0")
        collision = self.payload["witnesses"][
            "same_endpoint_different_Doob_contrast"
        ]
        self.assertEqual(collision["shared_endpoint_ranks_h0_h25"], [0, 2])
        self.assertNotEqual(
            collision["first_C_Doob_Q"], collision["second_C_Doob_Q"]
        )

    def test_scalar_partition_Rc_filtration_is_an_exact_no_go(self) -> None:
        no_go = self.payload["exhaustive_checks"][
            "partition_Rc_scalar_rank_no_go"
        ]
        self.assertTrue(no_go["C_Doob_Q_zero_all_black"])
        self.assertTrue(no_go["C_Doob_Q_zero_all_white"])
        self.assertIn(
            "marked cluster lineage",
            self.payload["minimal_extra_data"]["for_partition_Rc_chronology"],
        )

    def test_frozen_artifact_matches(self) -> None:
        artifact = json.loads(
            (
                ROOT
                / "results"
                / "exact-cover-character-oracles"
                / "p200_path_filtration.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(artifact, self.payload)


if __name__ == "__main__":
    unittest.main()
