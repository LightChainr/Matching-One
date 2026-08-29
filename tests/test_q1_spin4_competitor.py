import json
from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import q1_spin4_competitor as oracle  # noqa: E402


class Q1Spin4CompetitorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads((ROOT / "analysis" / "q1_spin4_competitor_manifest.yaml").read_text())
        cls.result = oracle.analyze(cls.config)

    def test_exact_q1_weights_dimension_and_spin(self):
        field = self.result["competitor"]
        self.assertEqual(field["left_weight"]["text"], "1/8")
        self.assertEqual(field["right_weight"]["text"], "33/8")
        self.assertEqual(field["scaling_dimension"]["text"], "17/4")
        self.assertEqual(field["conformal_spin"]["text"], "-4")
        self.assertEqual(field["leg_count"], 4)

    def test_central_charge_and_dimension_gap(self):
        self.assertEqual(self.result["central_charge"]["text"], "0")
        self.assertEqual(self.result["dimension_gap_thermal_minus_competitor"]["text"], "1")

    def test_chirality_conjugation_only_flips_spin(self):
        first = self.result["competitor"]
        second = self.result["chirality_conjugate"]
        self.assertEqual(first["scaling_dimension"], second["scaling_dimension"])
        self.assertEqual(Fraction(first["conformal_spin"]["text"]), -Fraction(second["conformal_spin"]["text"]))

    def test_relative_two_field_drift_is_q_minus_one_half(self):
        for row in self.result["continuum_dilation_oracle"]:
            expected = row["area_multiplier"] ** -0.5
            self.assertAlmostEqual(row["thermal_to_competitor_relative_factor"], expected, places=15)

    def test_physical_sector_claims_remain_unresolved(self):
        statuses = {row["id"]: row["status"] for row in self.result["unresolved_gates"]}
        self.assertTrue(all(value.startswith("unresolved") for value in statuses.values()))
        self.assertIn("does not establish", self.result["claim_boundary"])


if __name__ == "__main__":
    unittest.main()
