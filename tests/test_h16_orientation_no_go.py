
from __future__ import annotations
import unittest

from scripts.h16_orientation_no_go import (
    FIRST_FIVE_ORIENTATION_N,
    N1105,
    build_contract,
    first_full_rank_subset,
    harmonic_rank,
    minimum_primitive_layer,
    validate_contract,
)
from scripts.minimal_four_orientation_gaussian_torus import (
    primitive_first_octant_representations,
)


class H16OrientationNoGoTests(unittest.TestCase):
    def test_checked_in_contract_closes_exactly(self) -> None:
        self.assertEqual(validate_contract(), build_contract())

    def test_n1105_has_only_four_rows_for_five_harmonics(self) -> None:
        orientations = primitive_first_octant_representations(N1105)
        self.assertEqual(len(orientations), 4)
        self.assertEqual(harmonic_rank(N1105, orientations), 4)

    def test_first_layer_with_five_orbits_is_32045_with_eight(self) -> None:
        self.assertEqual(minimum_primitive_layer(4), (1105, 4, [5, 13, 17]))
        self.assertEqual(
            minimum_primitive_layer(5),
            (FIRST_FIVE_ORIENTATION_N, 8, [5, 13, 17, 29]),
        )
        self.assertEqual(
            len(primitive_first_octant_representations(FIRST_FIVE_ORIENTATION_N)),
            8,
        )

    def test_first_five_orientation_subset_has_exact_rank_five(self) -> None:
        subset = first_full_rank_subset(FIRST_FIVE_ORIENTATION_N)
        self.assertEqual(len(subset), 5)
        self.assertEqual(harmonic_rank(FIRST_FIVE_ORIENTATION_N, subset), 5)


if __name__ == "__main__":
    unittest.main()
