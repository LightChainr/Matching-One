import sys
import unittest
from fractions import Fraction
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify_norm5_conjugate_coalescence import (  # noqa: E402
    DESIGNS,
    affine_weights,
    d4_canonical,
    mul,
    norm,
    smith_invariants,
    verify_design,
)


class Norm5ConjugateCoalescenceTests(unittest.TestCase):
    def test_exact_products_coalesce_under_d4(self) -> None:
        for n, design in DESIGNS.items():
            parents = design["parents"]
            observed = design["observed_multiplier"]
            conjugate = design["conjugate_multiplier"]
            a, b, c = design["A"], design["B"], design["C"]
            self.assertEqual(
                tuple(d4_canonical(mul(parent, observed)) for parent in parents),
                (a, b),
            )
            self.assertEqual(
                tuple(d4_canonical(mul(parent, conjugate)) for parent in parents),
                (c, c),
            )
            self.assertEqual(norm(a), n)
            self.assertEqual(norm(b), n)
            self.assertEqual(norm(c), n)

    def test_noncyclic_child_is_the_only_smith_class_change(self) -> None:
        for n, design in DESIGNS.items():
            self.assertEqual(smith_invariants(design["A"]), (1, n))
            self.assertEqual(smith_invariants(design["B"]), (1, n))
            self.assertEqual(smith_invariants(design["C"]), (5, n // 5))

    def test_frozen_harmonic_weights(self) -> None:
        expected = {
            325: {
                4: (Fraction(11, 5), Fraction(-6, 5)),
                8: (Fraction(22517, 44795), Fraction(22278, 44795)),
                12: (Fraction(363263, 7144145), Fraction(6780882, 7144145)),
            },
            425: {
                4: (Fraction(-13, 20), Fraction(33, 20)),
                8: (Fraction(89531, 242420), Fraction(152889, 242420)),
                12: (Fraction(181189, 68620), Fraction(-112569, 68620)),
            },
        }
        for n, by_spin in expected.items():
            design = DESIGNS[n]
            for spin, weights in by_spin.items():
                self.assertEqual(
                    affine_weights(design["A"], design["B"], design["C"], spin),
                    weights,
                )
                self.assertEqual(sum(weights), 1)

    def test_h4_residuals_cancel_scalar_background(self) -> None:
        self.assertEqual(DESIGNS[325]["h4_integer_residual_C_A_B"], (5, -11, 6))
        self.assertEqual(DESIGNS[425]["h4_integer_residual_C_A_B"], (20, 13, -33))
        self.assertEqual(sum(DESIGNS[325]["h4_integer_residual_C_A_B"]), 0)
        self.assertEqual(sum(DESIGNS[425]["h4_integer_residual_C_A_B"]), 0)
        self.assertEqual(DESIGNS[325]["h4_conjugate_ratio"], Fraction(6, 11))
        self.assertEqual(DESIGNS[425]["h4_conjugate_ratio"], Fraction(33, 13))

    def test_prediction_and_execution_contract(self) -> None:
        prediction = yaml.safe_load(
            (ROOT / "predictions/norm5_conjugate_coalescence_20260829.yaml").read_text()
        )
        experiment = yaml.safe_load(
            (ROOT / "experiments/p205_norm5_conjugate_coalescence_20260829.yaml").read_text()
        )
        self.assertEqual(prediction["status"], "prospective_before_noncyclic_C_target_reveal")
        self.assertEqual(prediction["observable"]["p_ref"], 0.59274605079)
        self.assertEqual(prediction["primary_H4"]["N325"]["integer_residual_C_A_B"], [5, -11, 6])
        self.assertEqual(prediction["primary_H4"]["N425"]["integer_residual_C_A_B"], [20, 13, -33])
        self.assertEqual(prediction["fixed_angular_adversaries"]["frozen_order"], ["H12", "H8"])
        self.assertEqual(
            experiment["score_order"][:3],
            [
                "joint_H4_fixed_p_ref_affine_residual",
                "fixed_H12_affine_residual",
                "fixed_H8_affine_residual",
            ],
        )
        self.assertEqual(experiment["pilot"]["samples_per_pair"], 10000000)
        self.assertEqual(experiment["pilot"]["batches"], 100)
        for key in ("N325", "N425"):
            pair_a, pair_b = experiment[key]["pairs"]
            self.assertEqual(pair_a["first_representation"], experiment[key]["C_representation"])
            self.assertEqual(pair_b["first_representation"], experiment[key]["C_representation"])
            self.assertEqual(pair_a["seed"], pair_b["seed"])
            self.assertEqual(pair_a["replica_counter"], pair_b["replica_counter"])
        self.assertNotEqual(experiment["N325"]["pairs"][0]["seed"], experiment["N425"]["pairs"][0]["seed"])

    def test_verifier_runs_all_contracts(self) -> None:
        for n, design in DESIGNS.items():
            self.assertEqual(verify_design(n, design)["N"], n)


if __name__ == "__main__":
    unittest.main()
