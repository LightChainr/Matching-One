from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from flat_twist_representation import (  # noqa: E402
    R,
    T_SHEAR,
    action_permutation,
    build_certificate,
    compose,
    observed_f3,
    representation_certificate,
)


def power(permutation: tuple[int, ...], exponent: int) -> tuple[int, ...]:
    result = tuple(range(len(permutation)))
    for _ in range(exponent):
        result = compose(permutation, result)
    return result


class FlatTwistRepresentationTests(unittest.TestCase):
    def test_exact_finite_group_images_and_projectors(self) -> None:
        f2 = representation_certificate(2)
        f3 = representation_certificate(3)
        self.assertEqual(f2["D4_projective_image_order"], 2)
        self.assertEqual(f2["modular_projective_image_order"], 6)
        self.assertEqual(f2["D4_decomposition"], "2 A1 + B1")
        self.assertEqual(f3["D4_projective_image_order"], 4)
        self.assertEqual(f3["modular_projective_image_order"], 12)
        self.assertEqual(f3["D4_decomposition"], "2 A1 + B1 + B2")
        for certificate in (f2, f3):
            self.assertTrue(all(certificate["projector_gates"].values()))
            self.assertEqual(certificate["standard_character_inner_product"], "1")
            self.assertTrue(certificate["standard_is_irreducible"])

    def test_projective_modular_relations(self) -> None:
        for q in (2, 3):
            s = action_permutation(R, q)
            t = action_permutation(T_SHEAR, q)
            identity = tuple(range(q + 1))
            self.assertEqual(power(s, 2), identity)
            self.assertEqual(power(compose(s, t), 3), identity)
            self.assertEqual(power(t, q), identity)

    def test_n65_observation_is_a_charged_doublet_coordinate(self) -> None:
        score = ROOT / (
            "results/local-20260830/P334-projective-birth-N65-smoke/"
            "flat_twist_score.json"
        )
        observed = observed_f3(score)
        self.assertAlmostEqual(
            observed["orientation_contrast"][2], 0.0020775612521830144
        )
        self.assertAlmostEqual(
            observed["T_charged_doublet_quadratic"], 5.905028130859152
        )
        self.assertEqual(observed["T_charged_doublet_df"], 2)
        self.assertAlmostEqual(
            observed["no_fit_T_shear_prediction_in_D4_basis"][2],
            -0.0028238455708840906,
        )
        self.assertIn("charged doublet", observed["interpretation"])

    def test_machine_certificate_closes(self) -> None:
        score = ROOT / (
            "results/local-20260830/P334-projective-birth-N65-smoke/"
            "flat_twist_score.json"
        )
        certificate = build_certificate(score)
        self.assertTrue(certificate["all_exact_gates_pass"])
        self.assertIn("roughly 2-sigma", certificate["claim_boundary"])


if __name__ == "__main__":
    unittest.main()
