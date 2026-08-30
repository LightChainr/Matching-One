from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p321_graded_closure_extension import (  # noqa: E402
    DEFAULT_CERTIFICATE,
    block_diagonal,
    build_certificate,
    exact_nullity,
    identity,
    intertwiner_constraints,
)


class P321GradedClosureExtensionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.certificate = build_certificate()
        cls.rows = {row["width"]: row for row in cls.certificate["widths"]}

    def test_frozen_certificate_matches(self) -> None:
        frozen = json.loads(DEFAULT_CERTIFICATE.read_text(encoding="utf-8"))
        self.assertEqual(self.certificate, frozen)

    def test_block_helper_and_constraint_nullity(self) -> None:
        one = identity(1)
        doubled = block_diagonal(one, one)
        self.assertEqual(doubled, ((1, 0), (0, 1)))
        self.assertEqual(exact_nullity(intertwiner_constraints((one,), (one,))), 1)

    def test_graded_block_pull_through_is_exact(self) -> None:
        keys = (
            "J_constant_commutator",
            "J_Q_coefficient_commutator",
            "closure_seam_pull_through",
            "graded_defect_constant_pull_through",
            "graded_defect_Q_coefficient_pull_through",
        )
        for row in self.rows.values():
            residuals = row["block_pull_through"]
            self.assertTrue(residuals["all_sites_equal"])
            for key in keys:
                self.assertEqual(residuals[key]["rank_over_Q"], 0)

    def test_regular_Q_derivative_is_not_nilpotent(self) -> None:
        for width, dimension in ((2, 2), (3, 5), (4, 14)):
            derivative = self.rows[width]["regular_Q_derivative"]
            self.assertEqual(derivative["rank_J_prime"], dimension)
            self.assertEqual(derivative["rank_J_prime_squared"], dimension)
            self.assertEqual(derivative["W_prime_idempotency_residual"]["rank_over_Q"], 0)

    def test_Jordan_intertwiner_exists_but_is_not_unique(self) -> None:
        expected_dimensions = {2: (2, 2), 3: (5, 3), 4: (20, 9)}
        for width, (join_dimension, affine_dimension) in expected_dimensions.items():
            row = self.rows[width]
            dimensions = row["intertwiner_dimensions_over_Q"]
            self.assertEqual(
                dimensions["join_sigma_intertwiners_Xe_i=e_(i+1)X"], join_dimension
            )
            self.assertEqual(
                dimensions["affine_sigma_intertwiners_also_XT=TX"], affine_dimension
            )
            self.assertGreater(affine_dimension, 1)
            self.assertEqual(row["canonical_off_diagonal"]["N_X_squared"]["rank_over_Q"], 0)
            self.assertEqual(row["canonical_off_diagonal"]["pull_through"]["rank_over_Q"], 0)
            self.assertEqual(row["singular_confluence"]["residual_rank"], 0)

    def test_scientific_decision_prefers_scalar_module_amplitudes(self) -> None:
        decision = self.certificate["decision"]
        self.assertIn("independent scalar amplitudes", decision["F_t_classification"])
        self.assertIn("not forced", decision["F_t_classification"])
        self.assertIn("-pi_0D", self.certificate["q_lift_crosswalk"])


if __name__ == "__main__":
    unittest.main()
