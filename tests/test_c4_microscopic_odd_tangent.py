from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from c4_microscopic_odd_tangent import build_oracle  # noqa: E402


class C4MicroscopicOddTangentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_oracle()

    def test_single_theory_intertwiner_is_exact(self) -> None:
        row = self.result["single_theory_intertwiner"]
        self.assertTrue(row["C_squared_is_identity"])
        self.assertEqual(row["configurationwise_weight_signature_violations"], 0)

    def test_microscopic_scores_are_odd_and_fisher_orthogonal(self) -> None:
        row = self.result["microscopic_scores"]
        self.assertEqual(row["complement_odd_violations"], {"score_t": 0, "score_lambda": 0})
        self.assertEqual(row["Fisher_Gram"], [["40", "0"], ["0", "40"]])
        self.assertTrue(row["exact_thermal_orthogonality"])

    def test_unique_local_thermal_null_has_nonzero_staggered_matrix_element(self) -> None:
        row = self.result["local_thermal_orthogonalization"]
        self.assertEqual(row["alpha_star"], "3/64")
        self.assertTrue(row["unique_t_response_zero"])
        self.assertEqual(row["surviving_lambda_matrix_element"], "11/64")
        self.assertEqual(row["complement_odd_violations"], 0)
        self.assertEqual(row["response_rows_columns_t_lambda"]["local_thermal_null"], ["0", "11/64"])

    def test_coupling_null_is_locally_visible(self) -> None:
        row = self.result["coupling_space_zero"]
        self.assertEqual(row["direction_delta_t_delta_lambda"], ["2", "-3"])
        self.assertEqual(row["global_cross_matrix_element"], "0")
        self.assertEqual(row["local_h4_matrix_element"], "-39/64")

    def test_committed_artifact_is_reproducible(self) -> None:
        committed = json.loads(
            (ROOT / "results" / "exact-c4-microscopic-odd-tangent" / "latest.json").read_text()
        )
        self.assertEqual(committed, self.result)


if __name__ == "__main__":
    unittest.main()
