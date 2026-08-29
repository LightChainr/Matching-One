from __future__ import annotations

import sys
import unittest
from fractions import Fraction
from pathlib import Path

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from design_p275_atop_q4_field_identity import (  # noqa: E402
    REPRESENTATIVES,
    SIZES,
    build_payload,
    determinant,
    exact_modulus_oracle,
    period_matrix,
    phase4,
    smith_invariants,
)


class P275AtopQ4FieldIdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_payload(70)

    def test_selected_periods_have_exact_N_and_cyclic_smith_type(self) -> None:
        modulus = {"i": (1, 1), "2i": (2, 1), "5i_over_2": (5, 2)}
        for n in SIZES:
            for name, (numerator, denominator) in modulus.items():
                matrix = period_matrix(REPRESENTATIVES[n][name], numerator, denominator)
                self.assertEqual(determinant(matrix), n)
                self.assertEqual(smith_invariants(matrix), (1, n))

    def test_complex_frame_characters_are_exact_unit_phases(self) -> None:
        for rows in REPRESENTATIVES.values():
            for representation in rows.values():
                real, imag = phase4(representation)
                self.assertEqual(real * real + imag * imag, Fraction(1))
        self.assertEqual(phase4((2, 1)), (Fraction(-7, 25), Fraction(24, 25)))

    def test_modulus_oracle_keeps_exact_hecke_ratio_and_eta_cocycle(self) -> None:
        oracle = exact_modulus_oracle(70)
        self.assertEqual(oracle["exact_E4hat_2i_over_i"], "11/4")
        ratio = mp.mpf(oracle["eta_cocycle_ratio_2i_over_5i_over_2"])
        self.assertGreater(ratio, mp.mpf("0.66"))
        self.assertLess(ratio, mp.mpf("0.67"))
        self.assertGreater(
            mp.mpf(oracle["values"]["5i_over_2"]["E4hat_over_i"]),
            mp.mpf(oracle["values"]["2i"]["E4hat_over_i"]),
        )

    def test_model_dimensions_and_archive_boundary_are_frozen(self) -> None:
        models = self.payload["frozen_model_subspaces"]
        self.assertEqual(models["Q4_epsilon_ordinary"]["real_observations_parameters_df"], [18, 2, 16])
        self.assertEqual(models["Q4_energy_Jordan"]["real_observations_parameters_df"], [18, 6, 12])
        self.assertEqual(
            self.payload["archive_audit"]["conclusion"],
            "no committed archive can reconstruct the field-identity score",
        )
        self.assertFalse(self.payload["acquisition"]["production_authorized"])

    def test_scientific_card_is_exactly_five_lines(self) -> None:
        self.assertEqual(len(self.payload["scientific_card"]), 5)
        self.assertTrue(all(line.strip() for line in self.payload["scientific_card"]))


if __name__ == "__main__":
    unittest.main()
