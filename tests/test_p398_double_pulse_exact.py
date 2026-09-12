"""Small checks that prevent a false zero, false rank, or a pulse normalizer error."""
from fractions import Fraction
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import p398_double_pulse_exact as p


def determinant_fraction(a):
    b = [[Fraction(x) for x in row] for row in a]
    result = Fraction(1)
    for k in range(len(b)):
        pivot = next((i for i in range(k, len(b)) if b[i][k]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != k:
            b[k], b[pivot] = b[pivot], b[k]
            result = -result
        v = b[k][k]
        result *= v
        for i in range(k + 1, len(b)):
            factor = b[i][k] / v
            for j in range(k + 1, len(b)):
                b[i][j] -= factor * b[k][j]
    return result


class DoublePulseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = {w: p.analyse(w) for w in (4, 5)}

    def test_existing_model_counts_and_parity(self):
        for w, n in ((4, 14), (5, 42)):
            row = self.rows[w]
            self.assertEqual(len(row['states']), n)
            self.assertEqual(row['double_pulse']['odd_dimension'], n - row['baseline']['minimal_order'])
            # H is NOT a stochastic generator on its own; only G+epsilon H is.
            self.assertTrue(any(x < 0 for i, r in enumerate(row['H']) for j, x in enumerate(r) if i != j))

    def test_exact_first_order_zero_and_second_order_nonzero(self):
        for row in self.rows.values():
            src, g, h, f = (row[k] for k in ('sources_scaled', 'G', 'H', 'readouts'))
            first = p.multiply(p.multiply(src, h), f)
            self.assertTrue(all(x == 0 for r in first for x in r))
            second = p.multiply(p.multiply(p.multiply(src, h), h), f)
            self.assertEqual(second[0][0], 2 * row['source_denominator'])
            self.assertEqual(second[1], [0, 0, 0])

    def test_width4_scalar_resolvent_and_independent_determinant(self):
        scalar = self.rows[4]['double_pulse']['scalar_witness']
        self.assertEqual(scalar['moments'][:8], [0, 2, -20, 156, -1122, 7822, -53932, 371172])
        self.assertEqual(scalar['laplace_numerator_high_first'], [2, 22, 54])
        self.assertEqual(scalar['laplace_denominator_high_first'], [1, 21, 159, 513, 598])
        self.assertEqual(determinant_fraction(scalar['hankel_minor']), -16)

    def test_all_odd_modes_accessible(self):
        for w, d in ((4, 4), (5, 16)):
            data = self.rows[w]['double_pulse']
            cert = data['rank_certificate']
            self.assertEqual(cert['minimal_order'], d)
            self.assertNotEqual(p.determinant_mod(cert['integer_hankel_minor']), 0)
            self.assertNotEqual(p.determinant_mod(data['scalar_witness']['hankel_minor']), 0)

    def test_controlled_language_requires_more_than_baseline(self):
        for w, base, total in ((4, 10, 14), (5, 26, 42)):
            row = self.rows[w]
            self.assertEqual(row['baseline']['minimal_order'], base)
            self.assertEqual(row['arbitrary_controlled_word']['minimal_order'], total)
            self.assertEqual(row['controlled_source_contrasts']['minimal_order'], total - 1)
            for name in ('baseline', 'arbitrary_controlled_word', 'controlled_source_contrasts'):
                cert = row[name]
                self.assertEqual(p.determinant_mod(cert['integer_hankel_minor']),
                                 cert['hankel_minor_determinant_mod_prime'])

    def test_sign_not_fixed_by_parity(self):
        data = self.rows[4]['double_pulse']
        n = data['source_denominator']
        self.assertEqual(data['markov_moments_scaled'][0][0][0], 2 * n)
        self.assertEqual(data['markov_moments_scaled'][1][0][2], -2 * n)
        self.assertEqual(data['markov_moments_scaled'][1][2][2], 2 * n)

    def test_archived_half_strength_gives_quarter_kernel(self):
        row = self.rows[4]
        h = [[Fraction(x, 2) for x in r] for r in row['H']]
        for lag_power in range(3):
            right = p.multiply(h, row['readouts'])
            for _ in range(lag_power):
                right = p.multiply(row['G'], right)
            got = p.multiply(p.multiply(row['sources_scaled'], h), right)
            expected = row['double_pulse']['markov_moments_scaled'][lag_power]
            self.assertEqual(got, [[Fraction(x, 4) for x in r] for r in expected])


if __name__ == '__main__':
    unittest.main()
