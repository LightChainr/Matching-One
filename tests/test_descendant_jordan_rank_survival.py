
from __future__ import annotations
from fractions import Fraction
import unittest

from scripts.descendant_jordan_rank_survival import (
    apply,
    basis_vector,
    bottom_survives,
    build_contract,
    descendant_intertwines,
    dilatation,
    identity,
    jordan_chain_length,
    jordan_nilpotent,
    validate_contract,
)


class DescendantJordanRankSurvivalTests(unittest.TestCase):
    def test_checked_in_contract_closes_exactly(self) -> None:
        self.assertEqual(validate_contract(), build_contract())

    def test_nonzero_bottom_image_preserves_the_full_chain(self) -> None:
        for rank in range(1, 7):
            nilpotent = jordan_nilpotent(rank)
            descendant = identity(rank)
            top_image = apply(descendant, basis_vector(rank, rank - 1))
            self.assertTrue(bottom_survives(descendant))
            self.assertEqual(jordan_chain_length(nilpotent, top_image), rank)

    def test_q4_identity_has_the_required_level_four_intertwining(self) -> None:
        source = dilatation(2, Fraction(5, 4))
        target = dilatation(2, Fraction(21, 4))
        self.assertTrue(descendant_intertwines(source, target, identity(2), 4))

    def test_commuting_map_can_collapse_rank_when_it_kills_the_bottom(self) -> None:
        nilpotent = jordan_nilpotent(3)
        source = dilatation(3, Fraction(7, 3))
        target = dilatation(3, Fraction(13, 3))
        self.assertTrue(descendant_intertwines(source, target, nilpotent, 2))
        self.assertFalse(bottom_survives(nilpotent))
        top_image = apply(nilpotent, basis_vector(3, 2))
        self.assertEqual(jordan_chain_length(nilpotent, top_image), 2)


if __name__ == "__main__":
    unittest.main()
