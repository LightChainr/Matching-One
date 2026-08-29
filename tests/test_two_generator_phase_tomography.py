import unittest
from fractions import Fraction
from pathlib import Path

import yaml

from scripts.verify_two_generator_phase_tomography import chi4, mul


ROOT = Path(__file__).resolve().parents[1]


class TwoGeneratorPhaseTomographyTests(unittest.TestCase):
    def test_phase_nodes_and_norm5_norm10_rows(self):
        self.assertEqual(chi4((1, 1)), (Fraction(-1), Fraction(0)))
        self.assertEqual(chi4((1, -1)), (Fraction(-1), Fraction(0)))
        self.assertEqual(chi4((0, 2)), (Fraction(1), Fraction(0)))
        self.assertEqual(chi4((2, 1)), (Fraction(-7, 25), Fraction(24, 25)))
        self.assertEqual(chi4((2, -1)), (Fraction(-7, 25), Fraction(-24, 25)))
        self.assertEqual(chi4((3, 1)), (Fraction(7, 25), Fraction(24, 25)))
        self.assertEqual(chi4((3, -1)), (Fraction(7, 25), Fraction(-24, 25)))

    def test_character_composition(self):
        product = mul((1, 1), (2, -1))
        self.assertEqual(product, (3, 1))
        self.assertEqual(chi4(product), chi4((3, 1)))

        left = chi4((1, 1))
        right = chi4((2, -1))
        composed = (
            left[0] * right[0] - left[1] * right[1],
            left[0] * right[1] + left[1] * right[0],
        )
        self.assertEqual(composed, chi4((3, 1)))

    def test_quadrature_reconstruction(self):
        cosine = Fraction(7, 13)
        sine = Fraction(-5, 11)
        y5 = Fraction(-7, 25) * cosine + Fraction(24, 25) * sine
        y10 = Fraction(7, 25) * cosine + Fraction(24, 25) * sine

        self.assertEqual(Fraction(25, 14) * (y10 - y5), cosine)
        self.assertEqual(Fraction(25, 48) * (y10 + y5), sine)

    def test_phase_matrix_condition_number(self):
        # The columns are orthogonal. Squared singular values are 98/625 and
        # 1152/625, so their ratio of singular values is exactly 24/7.
        cos_sq = Fraction(98, 625)
        sin_sq = Fraction(1152, 625)
        self.assertEqual(sin_sq / cos_sq, Fraction(576, 49))
        self.assertEqual(Fraction(24, 7) ** 2, sin_sq / cos_sq)

    def test_norm4_guardrail_is_frozen_before_target(self):
        prediction = yaml.safe_load(
            (ROOT / "predictions" / "norm4_two_generator_transfer_20260829.yaml").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(prediction["status"], "frozen_before_N260_N340_target_production")
        guardrail = prediction["quotient_arithmetic_guardrail"]
        channels = {row["observable"]: row for row in guardrail["joint_channels"]}
        self.assertEqual(channels["P4[D]"]["target_ratio_target_over_parent"], "4^(-13/8)")
        self.assertEqual(channels["P4[D_prime]"]["target_ratio_target_over_parent"], "4^(-5/8)")
        self.assertEqual(channels["DeltaRoot"]["target_ratio_target_over_parent"], "1/16")
        self.assertIn("5%", guardrail["decision_rule"])
        self.assertEqual(
            prediction["execution_gate"][0],
            "score_the_pre_reveal_N290_conjugation_phase_node_first",
        )
        self.assertIn(
            "freeze_the_source_trained_full_jet_even_generator_if_claimed_prospectively",
            prediction["execution_gate"],
        )


if __name__ == "__main__":
    unittest.main()
