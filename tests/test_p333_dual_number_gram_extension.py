#!/usr/bin/env python3

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p333_dual_number_gram_extension import (  # noqa: E402
    analyze,
    expected_inertia,
    mobius_pivot_jet,
    restricted_first_form,
    symmetric_inertia,
    sharp_jordan_gate_oracle,
    join_semilattice_oracle,
    detach_join_semigroup_oracle,
    weighted_detach_join_jordan_oracle,
)


class DualNumberGramExtensionTests(unittest.TestCase):
    def test_pivot_jet_closed_form(self):
        self.assertEqual([mobius_pivot_jet(k) for k in range(2, 7)],
                         [1, -1, 2, -6, 24])

    def test_low_leg_inertia_is_exact(self):
        for n in range(2, 6):
            self.assertEqual(symmetric_inertia(restricted_first_form(n)),
                             expected_inertia(n))

    def test_no_first_order_null_direction_remains(self):
        payload = analyze(5)
        self.assertTrue(all(row["prediction_exact"] for row in payload["checks"]))
        self.assertTrue(all(row["first_radical_form_inertia"]["zero"] == 0
                            for row in payload["checks"]))

    def test_isotropic_jordan_gate_is_sharp(self):
        oracle = sharp_jordan_gate_oracle()
        self.assertTrue(oracle["gram_self_adjoint"])
        self.assertEqual(oracle["bottom_norm"], "0")

    def test_indefinite_radical_starts_at_three_marks(self):
        self.assertEqual(expected_inertia(2), (1, 0, 0))
        self.assertEqual(expected_inertia(3), (3, 1, 0))

    def test_join_only_algebra_is_semisimple(self):
        oracle = join_semilattice_oracle(4)
        self.assertTrue(oracle["all_idempotent"])
        self.assertTrue(oracle["all_commuting"])
        self.assertTrue(oracle["all_first_jet_gram_self_adjoint"])

    def test_detach_is_minimal_noncoarsening_positive_control(self):
        self.assertEqual(detach_join_semigroup_oracle(2)["defective_deterministic_elements"], 0)
        n3 = detach_join_semigroup_oracle(3)
        self.assertEqual(n3["first_defective_word"], ["D0", "J01", "D1"])
        self.assertEqual(n3["first_defective_height"], 2)
        self.assertEqual(n3["gram_self_adjoint_defective_elements"], 0)
        self.assertEqual(detach_join_semigroup_oracle(4)["gram_self_adjoint_defective_elements"], 0)

    def test_signed_history_jordan_control(self):
        oracle = weighted_detach_join_jordan_oracle()
        self.assertTrue(oracle["formula_exact"])
        self.assertEqual(oracle["rank"], 1)
        self.assertTrue(oracle["nonzero"])
        self.assertTrue(oracle["square_zero"])
        self.assertTrue(oracle["first_jet_gram_self_adjoint"])
        self.assertTrue(oracle["chain"]["K_partner_equals_bottom"])
        self.assertTrue(oracle["chain"]["K_bottom_zero"])
        self.assertEqual(oracle["chain"]["bottom_H_norm"], "0")


if __name__ == "__main__":
    unittest.main()
