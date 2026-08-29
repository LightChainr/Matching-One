from __future__ import annotations

from fractions import Fraction
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from square_bond_duality_exact import (  # noqa: E402
    CHANNELS,
    enumerate_exact,
    geometric_dual_mask,
    primal_dual_wrapping,
    square_bond_pairs,
)

L2_EVEN_EXPECTATIONS = {
    "cross": Fraction(69, 256),
    "both": Fraction(73, 256),
    "either": Fraction(187, 256),
    "direction_0": Fraction(65, 128),
    "direction_1": Fraction(65, 128),
}
L3_EVEN_EXPECTATIONS = {
    "cross": Fraction(18865, 65536),
    "both": Fraction(20785, 65536),
    "either": Fraction(46671, 65536),
    "direction_0": Fraction(527, 1024),
    "direction_1": Fraction(527, 1024),
}


class SquareBondDualityExactTests(unittest.TestCase):
    def test_l2_geometric_dual_transport_swaps_primal_and_dual(self) -> None:
        result = enumerate_exact(2)
        self.assertTrue(result["passed"])
        self.assertEqual(result["N_bonds"], 8)
        self.assertEqual(result["configurations"], 256)
        self.assertEqual(result["geometric_dual_transport_failures"], 0)
        self.assertEqual(result["even_odd_involution_failures"], 0)
        self.assertTrue(result["geometric_dual_map_bijective"])

    def test_l2_naive_complement_is_not_the_duality_map(self) -> None:
        result = enumerate_exact(2)
        self.assertEqual(result["naive_complement_swap_failures"], 138)

    def test_l2_odd_sector_vanishes_exactly_at_half(self) -> None:
        result = enumerate_exact(2)
        for name in CHANNELS:
            row = result["channels"][name]
            self.assertTrue(row["D_vanishes_at_half"], name)
            self.assertEqual(row["E_D"]["fraction"], "0")
            self.assertTrue(row["primal_dual_equidistributed"], name)

    def test_l2_even_sector_locked_rationals(self) -> None:
        result = enumerate_exact(2)
        for name, expected in L2_EVEN_EXPECTATIONS.items():
            row = result["channels"][name]
            self.assertEqual(Fraction(row["E_S"]["fraction"]), expected, name)
            self.assertFalse(row["S_identically_zero"], name)
        self.assertEqual(
            L2_EVEN_EXPECTATIONS["either"],
            L2_EVEN_EXPECTATIONS["direction_0"]
            + L2_EVEN_EXPECTATIONS["direction_1"]
            - L2_EVEN_EXPECTATIONS["both"],
        )

    def test_single_mask_geometric_transport_swaps_wrapping(self) -> None:
        pairs = square_bond_pairs(2)
        mask = 0b00101101
        primal, dual = primal_dual_wrapping(2, mask, pairs)
        transported = geometric_dual_mask(mask, pairs)
        self.assertNotEqual(transported, mask ^ 255)
        swapped_primal, swapped_dual = primal_dual_wrapping(2, transported, pairs)
        self.assertEqual(primal.as_dict(), swapped_dual.as_dict())
        self.assertEqual(dual.as_dict(), swapped_primal.as_dict())

    def test_l3_identities(self) -> None:
        result = enumerate_exact(3)
        self.assertTrue(result["passed"])
        self.assertEqual(result["N_bonds"], 18)
        self.assertEqual(result["configurations"], 1 << 18)
        self.assertEqual(result["geometric_dual_transport_failures"], 0)
        self.assertEqual(result["naive_complement_swap_failures"], 147560)
        for name, expected in L3_EVEN_EXPECTATIONS.items():
            row = result["channels"][name]
            self.assertEqual(row["E_D"]["fraction"], "0", name)
            self.assertTrue(row["primal_dual_equidistributed"], name)
            self.assertEqual(Fraction(row["E_S"]["fraction"]), expected, name)
        self.assertEqual(
            L3_EVEN_EXPECTATIONS["either"],
            L3_EVEN_EXPECTATIONS["direction_0"]
            + L3_EVEN_EXPECTATIONS["direction_1"]
            - L3_EVEN_EXPECTATIONS["both"],
        )


if __name__ == "__main__":
    unittest.main()
