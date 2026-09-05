
from __future__ import annotations
import json
from pathlib import Path
import tempfile
import unittest

from scripts.cyclic_deck_invariant_monomials import (
    build_contract,
    invariant_census,
    is_primitive_neutral,
    neutral_monomials,
    primitive_neutral_generators,
    total_charge,
    validate_contract,
)


class CyclicDeckInvariantMonomialTests(unittest.TestCase):
    def test_checked_in_contract_closes_exactly(self) -> None:
        self.assertEqual(validate_contract(), build_contract())

    def test_hilbert_counts_through_degree_five_are_exact(self) -> None:
        self.assertEqual(invariant_census(2, 5)["hilbert_counts"], [1, 0, 1, 0, 1, 0])
        self.assertEqual(invariant_census(5, 5)["hilbert_counts"], [1, 0, 2, 4, 7, 12])

    def test_c5_quadratics_are_exactly_the_conjugate_pairs(self) -> None:
        self.assertEqual(neutral_monomials(5, 2), [(0, 1, 1, 0), (1, 0, 0, 1)])
        for exponent in neutral_monomials(5, 5):
            self.assertEqual(total_charge(5, exponent), 0)

    def test_primitive_generators_are_indecomposable(self) -> None:
        generators = primitive_neutral_generators(5, 5)
        self.assertEqual(len(generators), 14)
        self.assertTrue(all(is_primitive_neutral(5, exponent) for exponent in generators))
        self.assertEqual(
            [exponent for exponent in generators if sum(exponent) == 5],
            [(0, 0, 0, 5), (0, 0, 5, 0), (0, 5, 0, 0), (5, 0, 0, 0)],
        )

    def test_contract_drift_and_invalid_shape_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "every nontrivial"):
            total_charge(5, (1, 4))
        frozen = build_contract()
        frozen["c5"]["hilbert_counts"][5] = 11
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "drift.json"
            path.write_text(json.dumps(frozen), encoding="utf-8")
            with self.assertRaisesRegex(AssertionError, "contract drifted"):
                validate_contract(path)


if __name__ == "__main__":
    unittest.main()
