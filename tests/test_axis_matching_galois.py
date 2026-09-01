from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from certify_axis_matching_galois import (  # noqa: E402
    AXIS_POWER_POLYNOMIALS,
    FROBENIUS_CERTIFICATES,
    factor_degree_partition,
    l2_splitting_field_degree,
    pairwise_gcds,
    run_suite,
)
from exact_matching_polynomial import bernstein_counts, bernstein_to_power  # noqa: E402
from matched_torus_reference import axis_geometry  # noqa: E402


class AxisMatchingGaloisTests(unittest.TestCase):
    def test_locked_coefficients_match_committed_enumeration_for_l2_l3(self) -> None:
        for L in (2, 3):
            enumerated = bernstein_to_power(bernstein_counts(axis_geometry(L)))
            self.assertEqual(enumerated, AXIS_POWER_POLYNOMIALS[L], L)

    def test_pairwise_gcds_are_one(self) -> None:
        gcds = pairwise_gcds()
        self.assertEqual(gcds, {"L2_L3": "1", "L2_L4": "1", "L3_L4": "1"})

    def test_l2_is_c4_not_s4(self) -> None:
        tower = l2_splitting_field_degree()
        self.assertTrue(tower["even_in_p2"])
        self.assertEqual(tower["quadratic_in_u_discriminant"], 8)
        self.assertFalse(tower["quadratic_in_u_disc_is_square"])
        self.assertEqual(tower["splitting_field_degree"], 4)
        self.assertEqual(tower["group"], "C4")
        self.assertEqual(factor_degree_partition(AXIS_POWER_POLYNOMIALS[2], 3), (4,))
        self.assertEqual(factor_degree_partition(AXIS_POWER_POLYNOMIALS[2], 7), (2, 2))

    def test_l3_frobenius_witnesses_force_s9(self) -> None:
        poly = AXIS_POWER_POLYNOMIALS[3]
        spec = FROBENIUS_CERTIFICATES[3]
        self.assertEqual(
            factor_degree_partition(poly, spec["irreducible_prime"]),
            spec["irreducible_partition"],
        )
        self.assertEqual(
            factor_degree_partition(poly, spec["primitive_prime"]),
            spec["primitive_partition"],
        )
        self.assertEqual(
            factor_degree_partition(poly, spec["transposition_prime"]),
            spec["transposition_partition"],
        )

    def test_l4_frobenius_witnesses_force_s16(self) -> None:
        poly = AXIS_POWER_POLYNOMIALS[4]
        spec = FROBENIUS_CERTIFICATES[4]
        self.assertEqual(
            factor_degree_partition(poly, spec["irreducible_prime"]),
            spec["irreducible_partition"],
        )
        self.assertEqual(
            factor_degree_partition(poly, spec["primitive_prime"]),
            spec["primitive_partition"],
        )
        self.assertEqual(
            factor_degree_partition(poly, spec["transposition_prime"]),
            spec["transposition_partition"],
        )

    def test_suite_does_not_touch_l5(self) -> None:
        payload = run_suite()
        self.assertTrue(payload["excludes_L5"])
        self.assertTrue(payload["excludes_PR84"])
        groups = {row["L"]: row["group"] for row in payload["certificates"]}
        self.assertEqual(groups, {2: "C4", 3: "S9", 4: "S16"})
        self.assertNotIn(5, groups)
        for row in payload["certificates"]:
            self.assertTrue(row["not_a_statement_about_infinite_pc"])


if __name__ == "__main__":
    unittest.main()
