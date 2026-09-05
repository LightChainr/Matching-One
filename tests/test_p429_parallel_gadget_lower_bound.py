import unittest
from fractions import Fraction

from scripts.p429_parallel_gadget_lower_bound import (
    EXIT_A,
    EXIT_B,
    fork_probability,
    fork_probability_by_successor_sum,
    safe_counts,
    survival_law,
    verify,
)


class ParallelGadgetLowerBoundTests(unittest.TestCase):
    def test_base_successor_moments(self):
        self.assertEqual(sum(EXIT_A), 13)
        self.assertEqual(sum(EXIT_B), 13)
        self.assertEqual(sum(x * x for x in EXIT_A), 29)
        self.assertEqual(sum(x * x for x in EXIT_B), 25)

    def test_base_witness(self):
        self.assertEqual(fork_probability(1, 0), Fraction(93, 196))
        self.assertEqual(fork_probability(1, 1), Fraction(95, 196))
        self.assertEqual(
            fork_probability(1, 1) - fork_probability(1, 0), Fraction(1, 98)
        )

    def test_parallel_closed_form(self):
        for k in range(1, 9):
            for a in range(k + 1):
                self.assertEqual(
                    fork_probability(k, a),
                    fork_probability_by_successor_sum(k, a),
                )

    def test_strict_class_splitting(self):
        for k in range(1, 9):
            expected = Fraction(1, 2 * k * (8 * k - 1) ** 2)
            values = [fork_probability(k, a) for a in range(k + 1)]
            self.assertEqual(len(set(values)), k + 1)
            for left, right in zip(values, values[1:]):
                self.assertEqual(right - left, expected)

    def test_survival_polynomial_examples(self):
        self.assertEqual(
            safe_counts(2)[:9],
            (1, 14, 85, 292, 620, 832, 688, 320, 64),
        )
        self.assertEqual(survival_law(1)[:5], (
            Fraction(1, 1),
            Fraction(7, 8),
            Fraction(9, 14),
            Fraction(5, 14),
            Fraction(4, 35),
        ))

    def test_full_verifier(self):
        verify(12)


if __name__ == "__main__":
    unittest.main()
