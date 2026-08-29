from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from null_field_residue_ward import (  # noqa: E402
    C,
    H,
    Dual,
    build_artifact,
    covectors,
    fixed_c_jordan_null_jet,
    kac_tangent_null_jet,
    level2_kac_central_charge,
    ordinary_null_jet,
    q4_scalar_ward_vector,
    singular_levels,
)


Q = Fraction


class NullFieldResidueWardTests(unittest.TestCase):
    def test_level2_kac_curve_and_tangent(self) -> None:
        value = level2_kac_central_charge(Dual(H, Q(1)))
        self.assertEqual(value.value, 0)
        self.assertEqual(value.tangent, Q(-40, 9))
        self.assertEqual(singular_levels()[:3], [2, 10, 16])

    def test_inhomogeneous_null_jets_separate_three_mechanisms(self) -> None:
        self.assertEqual(ordinary_null_jet(H), (Q(0), Q(0), Q(0)))
        self.assertEqual(fixed_c_jordan_null_jet(), (Q(-8, 3), Q(-10, 3), Q(0)))
        self.assertEqual(kac_tangent_null_jet(), (Q(-8, 3), Q(-10, 3), Q(-20, 9)))
        self.assertEqual(kac_tangent_null_jet()[2] / kac_tangent_null_jet()[1], Q(2, 3))

    def test_four_leg_and_negative_dimension_controls_have_exact_signatures(self) -> None:
        self.assertEqual(ordinary_null_jet(Q(33, 8)), (Q(-28, 3), Q(-77), Q(0)))
        self.assertEqual(ordinary_null_jet(Q(1, 8)), (Q(4, 3), Q(1, 3), Q(0)))
        negative = ordinary_null_jet(Q(-1, 24))
        self.assertGreater(negative[0], 0)
        self.assertLess(negative[1], 0)

    def test_q4_positive_mode_vector_and_jordan_residue(self) -> None:
        base = tuple(value.value for value in q4_scalar_ward_vector(Dual(H), Dual(C)))
        self.assertEqual(base, (Q(40), Q(-60), Q(30)))
        fixed = tuple(
            value.tangent
            for value in q4_scalar_ward_vector(Dual(H, 1), Dual(C, 0))
        )
        self.assertEqual(fixed, (Q(864), Q(-546), Q(48)))
        self.assertEqual(covectors(base), (Q(0), Q(0)))
        self.assertEqual(covectors(fixed), (Q(1500), Q(-450)))
        self.assertEqual(covectors(fixed)[0] / covectors(fixed)[1], Q(-10, 3))

    def test_kac_curve_q4_residue_is_distinct(self) -> None:
        c_path = level2_kac_central_charge(Dual(H, 1))
        tangent = tuple(
            value.tangent for value in q4_scalar_ward_vector(Dual(H, 1), c_path)
        )
        self.assertEqual(tangent, (Q(-616), Q(362, 3), Q(248)))
        self.assertEqual(covectors(tangent), (Q(-4820, 3), Q(1850, 3)))
        self.assertEqual(covectors(tangent)[0] / covectors(tangent)[1], Q(-482, 185))

    def test_artifact_freezes_primary_before_ratio_scoring(self) -> None:
        artifact = build_artifact()
        prediction = artifact["frozen_high_risk_prediction"]
        self.assertIn("4:-6:3", prediction["stage_1_operator_gate"])
        self.assertIn("-10/3", prediction["stage_2_module_gate"])
        self.assertIn("linear constraints before ratios", artifact["minimal_execution"]["covariance"])


if __name__ == "__main__":
    unittest.main()
