
from __future__ import annotations
import unittest

from scripts.landing_registry_spin_alias import (
    axis_minus_diagonal,
    build_contract,
    cosine_response,
    matrix_rank,
    sine_response,
    validate_contract,
)


class LandingRegistrySpinAliasTests(unittest.TestCase):
    def test_checked_in_contract_closes_exactly(self) -> None:
        self.assertEqual(validate_contract(), build_contract())

    def test_axis_minus_diagonal_selects_four_modulo_eight(self) -> None:
        for spin in range(0, 101, 4):
            expected = 2 if spin % 8 == 4 else 0
            self.assertEqual(axis_minus_diagonal(spin), expected)

    def test_sine_quadrature_is_zero_for_every_c4_allowed_spin(self) -> None:
        for spin in range(0, 101, 4):
            self.assertEqual(sine_response(spin), [0, 0])
        with self.assertRaisesRegex(ValueError, "requires spin 4k"):
            cosine_response(6)


if __name__ == "__main__":
    unittest.main()
