import json
import math
from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import q1_spin4_velocity_oracle as oracle  # noqa: E402


class Q1Spin4VelocityOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads((ROOT / "analysis" / "q_velocity_spin4_manifest.yaml").read_text())
        cls.result = oracle.analyze(cls.config)
        cls.fields = {field["id"]: field for field in cls.result["fields"]}

    def test_primary_exact_velocities(self):
        four_leg = self.fields["loop_V_2_2"]
        thermal = self.fields["thermal_Q4_epsilon"]
        self.assertEqual(four_leg["scaling_dimension"]["text"], "17/4")
        self.assertEqual(thermal["scaling_dimension"]["text"], "21/4")
        self.assertEqual(four_leg["dx_dQ_at_q1"]["coefficient_of_sqrt3_over_pi"]["text"], "-5/16")
        self.assertEqual(thermal["dx_dQ_at_q1"]["coefficient_of_sqrt3_over_pi"]["text"], "-9/16")
        self.assertEqual(self.result["primary_pair"]["velocity_four_leg_minus_thermal"]["coefficient_of_sqrt3_over_pi"]["text"], "1/4")

    def test_chirality_flips_spin_not_dimension_or_velocity(self):
        left = self.fields["loop_V_2_2"]
        right = self.fields["loop_V_2_minus2"]
        self.assertEqual(left["scaling_dimension"], right["scaling_dimension"])
        self.assertEqual(left["dx_dQ_at_q1"], right["dx_dQ_at_q1"])
        self.assertEqual(Fraction(left["conformal_spin"]["text"]), -Fraction(right["conformal_spin"]["text"]))

    def test_same_r_loop_spin8_spin12_controls(self):
        spin8 = self.fields["loop_V_2_4_spin8_control"]
        spin12 = self.fields["loop_V_2_6_spin12_control"]
        self.assertEqual(spin8["scaling_dimension"]["text"], "53/4")
        self.assertEqual(spin8["dx_dQ_at_q1"]["coefficient_of_sqrt3_over_pi"]["text"], "-41/16")
        self.assertEqual(spin12["scaling_dimension"]["text"], "113/4")
        self.assertEqual(spin12["dx_dQ_at_q1"]["coefficient_of_sqrt3_over_pi"]["text"], "-101/16")

    def test_chain_rule_matches_finite_difference(self):
        u0 = 2 / 3
        step = 1e-7
        def q(u):
            return 4 * math.cos(math.pi * u) ** 2
        def x22(u):
            return 1 + 1.5 * (u + 1 / u)
        numeric = (x22(u0 + step) - x22(u0 - step)) / (q(u0 + step) - q(u0 - step))
        exact = -5 * math.sqrt(3) / (16 * math.pi)
        self.assertAlmostEqual(numeric, exact, places=8)

    def test_angular_aliases_remain_unassigned(self):
        aliases = self.result["angular_aliases_without_q_family"]
        self.assertEqual(len(aliases), 2)
        self.assertTrue(all(row["status"].startswith("unresolved") for row in aliases))
        self.assertIn("forbidden", self.result["selection_boundary"])


if __name__ == "__main__":
    unittest.main()
