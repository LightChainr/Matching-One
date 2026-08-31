import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from p398_rooted_second_jet import build_result


class TestRootedSecondJet(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_second_jet_is_uniquely_zero(self):
        row = self.result["second_jet_affine_certificate"]
        self.assertEqual(row["unique_solution"], "X2=0")
        self.assertEqual(row["total_rank_in_remaining_variables"], 322)
        self.assertEqual(row["full_matrix_variable_count"], 529)

    def test_involution_and_no_next_layer(self):
        self.assertEqual(self.result["fixed_radical_translation"]["involution_residual_rank"], 0)
        self.assertEqual(self.result["Jantzen_certificate"]["valuation_multiplicities"], {"0": 19, "1": 4})
        self.assertEqual(self.result["Jantzen_certificate"]["J2_dimension"], 0)
        for row in self.result["Gram_Taylor_coefficients"]:
            self.assertEqual(row["C4_covariance_residual_rank"], 0)
            self.assertEqual(row["projected_Gram_skew_rank"], 0)

    def test_artifact_recomputes(self):
        self.assertEqual(self.result, json.loads((ROOT / "results/p398-rooted-second-jet/latest.json").read_text()))


if __name__ == "__main__":
    unittest.main()
