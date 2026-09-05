
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_profile_distance_certificate import (  # noqa: E402
    build_artifact,
    exact_profile_distance,
    polynomial_inner_product,
)


class ExactProfileDistanceCertificateTests(unittest.TestCase):
    def test_checked_artifact_reproduces(self) -> None:
        artifact = build_artifact()
        checked = json.loads(
            (ROOT / "analysis" / "exact_profile_distance_certificate.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(artifact, checked)

    def test_frozen_uniform_distances_are_exact(self) -> None:
        artifact = build_artifact()
        self.assertEqual(
            artifact["distance"],
            {
                "density_L2_squared": "1/20",
                "cdf_Cramer_von_Mises_squared": "1/840",
            },
        )

    def test_symmetry_and_identity(self) -> None:
        frozen = [Fraction(1, 2), 3, -3]
        uniform = [1]
        self.assertEqual(
            exact_profile_distance(frozen, uniform),
            exact_profile_distance(uniform, frozen),
        )
        self.assertEqual(
            exact_profile_distance(frozen, frozen),
            {"density_L2_squared": 0, "cdf_Cramer_von_Mises_squared": 0},
        )

    def test_inner_product_is_exact_and_bilinear(self) -> None:
        first = [Fraction(1), Fraction(2)]
        second = [Fraction(-1), Fraction(3)]
        self.assertEqual(
            polynomial_inner_product([2 * value for value in first], second),
            2 * polynomial_inner_product(first, second),
        )
        self.assertEqual(
            polynomial_inner_product(first, second),
            polynomial_inner_product(second, first),
        )

    def test_unnormalized_inputs_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "integrate exactly"):
            exact_profile_distance([2], [1])
        with self.assertRaisesRegex(ValueError, "must not be empty"):
            exact_profile_distance([], [1])

    def test_control_gram_is_positive_semidefinite(self) -> None:
        artifact = build_artifact()
        gram = artifact["control_gram"]
        self.assertEqual(
            artifact["exact_affine_mixture_identity"],
            "frozen = (uniform + Beta(2,2))/2",
        )
        self.assertTrue(gram["positive_semidefinite_certified"])
        self.assertEqual(Fraction(gram["determinant"]), 0)


if __name__ == "__main__":
    unittest.main()
