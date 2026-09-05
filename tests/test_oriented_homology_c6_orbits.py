
from __future__ import annotations
import unittest

from scripts.oriented_homology_c6_orbits import (
    build_contract,
    c6_character_inner,
    c6_orbits,
    descends_to_unoriented_lines,
    negate,
    oriented_primitive_vectors,
    oriented_spin_charge,
    quotient_c3_charge,
    rotate,
    validate_contract,
)


class OrientedHomologyC6OrbitTests(unittest.TestCase):
    def test_checked_in_contract_closes_exactly(self) -> None:
        self.assertEqual(validate_contract(), build_contract())

    def test_rotation_has_order_six_and_cube_is_sign(self) -> None:
        for vector in oriented_primitive_vectors(31):
            third = vector
            sixth = vector
            for _ in range(3):
                third = rotate(third)
            for _ in range(6):
                sixth = rotate(sixth)
            self.assertEqual(third, negate(vector))
            self.assertEqual(sixth, vector)

    def test_c6_orbits_partition_every_oriented_vector(self) -> None:
        vectors = oriented_primitive_vectors(31)
        flattened = [vector for orbit in c6_orbits(31) for vector in orbit]
        self.assertEqual(sorted(flattened), vectors)
        self.assertEqual(len(flattened), len(set(flattened)))

    def test_character_quotient_is_exact(self) -> None:
        gram = [[c6_character_inner(a, b) for b in range(6)] for a in range(6)]
        self.assertEqual(gram, [[6 if a == b else 0 for b in range(6)] for a in range(6)])
        self.assertEqual(
            [quotient_c3_charge(oriented_spin_charge(spin)) for spin in (4, 8, 12)],
            [2, 1, 0],
        )
        self.assertFalse(descends_to_unoriented_lines(oriented_spin_charge(3)))
        with self.assertRaisesRegex(ValueError, "odd C6 charge"):
            quotient_c3_charge(3)


if __name__ == "__main__":
    unittest.main()
