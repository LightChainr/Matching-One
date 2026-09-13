"""Focused exact regressions; no statistical or continuum claim is tested here."""
import hashlib
import itertools
import json
from collections import Counter
from fractions import Fraction
from math import comb
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import p401_cooperative_continuation as cc
from p334_birth_age_collision_review_20260830 import enumerate_states, memory_witness


class CooperativeContinuationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificate = cc.build_certificate()
        _, cls.s10 = enumerate_states(3, 1)
        _, cls.s13 = enumerate_states(3, 2)

    def test_inherited_topology_blob_is_unchanged(self):
        data = (ROOT/'scripts/p334_birth_age_collision_review_20260830.py').read_bytes()
        blob = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        self.assertEqual(blob, '62e06795fdfa91a956aedd62b7344e84aa5efc5c')

    def test_all_bounded_gates(self):
        gates = self.certificate['gates']
        self.assertEqual([g['rank_one_states'] for g in gates], [10, 162, 310, 2340])
        self.assertEqual(sum(g['two_step_identity_cases'] for g in gates), 2822)
        self.assertEqual(sum(g['three_step_identity_cases'] for g in gates), 2812)

    def test_new_history_gaps(self):
        witnesses = self.certificate['history_witnesses']
        self.assertEqual([w['latest_minus_earliest'] for w in witnesses],
                         ['1/44', '1/154', '2/315'])
        self.assertEqual([r['prefixes'] for r in witnesses[0]['rows']], [960, 2640])
        self.assertEqual([r['probability'] for r in witnesses[0]['rows']], ['13/20', '37/55'])
        self.assertEqual(self.certificate['independent_prefix_check']['ordered_prefixes_checked'], 30240)

    def test_original_one_step_witness_is_preserved(self):
        self.assertEqual(memory_witness(10, self.s10)['hazard_difference'], '1/57')

    def test_configuration_and_overlap_witnesses(self):
        a = cc.continuation(10, self.s10, 155)
        b = cc.continuation(10, self.s10, 157)
        self.assertEqual((a['x'], b['x'], a['c2'], b['c2']), (1, 1, 2, 3))
        a = cc.continuation(13, self.s13, 655)
        b = cc.continuation(13, self.s13, 693)
        self.assertEqual((a['x'], b['x'], a['c2'], b['c2']), (0, 0, 5, 5))
        self.assertEqual((a['pair_wedges'], b['pair_wedges']), (4, 6))
        self.assertEqual((a['minimal_triples'], b['minimal_triples']), (0, 3))
        self.assertEqual((a['exit_probability']['3'], b['exit_probability']['3']), ('3/5', '22/35'))

    def test_two_future_variance(self):
        row = self.certificate['two_future_control']
        self.assertEqual(Fraction(row['mean_q2_squared'])-Fraction(row['mean_q2'])**2,
                         Fraction(1, 450))
        self.assertEqual(row['variance_after_also_controlling_c2'], '0')

    def test_full_survival_signature_by_independent_future_orders(self):
        for mask in (155, 157):
            vacant = [v for v in range(10) if not mask >> v & 1]
            counts = Counter()
            orders = list(itertools.permutations(vacant))
            for order in orders:
                current = mask
                counts[0] += 1
                for m, v in enumerate(order, 1):
                    current |= 1 << v
                    counts[m] += self.s10[current][0] == 1
            signature = cc.survival_signature(10, self.s10, mask)
            for m, value in enumerate(signature):
                self.assertEqual(Fraction(counts[m], len(orders)), Fraction(value, comb(len(vacant), m)))

    def test_second_derivative_by_full_bernstein_expansion(self):
        for n, states, mask in ((10, self.s10, 155), (10, self.s10, 157),
                                 (13, self.s13, 655), (13, self.s13, 693)):
            b = cc.survival_signature(n, states, mask)
            d = len(b)-1
            power = [0]*(d+1)
            for m in range(d+1):
                successful = comb(d, m)-b[m]
                for j in range(d-m+1):
                    power[m+j] += successful*comb(d-m, j)*(-1)**j
            row = cc.continuation(n, states, mask)
            self.assertEqual(power[0], 0)
            self.assertEqual(power[1], row['x'])
            self.assertEqual(2*power[2], row['bernoulli_second_derivative_at_zero'])

    def test_n13_three_step_history_by_conditioned_prefix_enumeration(self):
        totals = Counter(); weighted = Counter(); checked = 0
        for mask, state in enumerate(self.s13):
            if mask.bit_count() != 6 or state != (1, (0, 1)):
                continue
            row = cc.continuation(13, self.s13, mask)
            if (row['x'], row['c2']) != (0, 5):
                continue
            occupied = [v for v in range(13) if mask >> v & 1]
            vacant = [v for v in range(13) if not mask >> v & 1]
            # Direct terminal ranks, not the three-step overlap formula.
            successes = sum(self.s13[mask | sum(1 << v for v in triple)][0] == 2
                            for triple in itertools.combinations(vacant, 3))
            for prefix in itertools.permutations(occupied):
                checked += 1
                current = 0; birth = None
                for j, v in enumerate(prefix, 1):
                    current |= 1 << v
                    if birth is None and self.s13[current][0] >= 1:
                        birth = j
                totals[birth] += 1; weighted[birth] += successes
        self.assertEqual(checked, 37440)
        self.assertEqual(totals, {5: 9360, 6: 28080})
        self.assertEqual(weighted, {5: 199680, 6: 605280})
        self.assertEqual(Fraction(weighted[6], 35*totals[6])-Fraction(weighted[5], 35*totals[5]),
                         Fraction(2, 315))

    def test_reject_invalid_shape_and_start(self):
        with self.assertRaises(ValueError): cc.validate_states(10, self.s10[:-1])
        with self.assertRaises(ValueError): cc.continuation(10, self.s10, 0)
        with self.assertRaises(ValueError): cc.survival_signature(10, self.s10, 1 << 10)

    def test_reject_nonmonotone_cache(self):
        bad = list(self.s10); bad[-2] = (0, None)
        with self.assertRaisesRegex(ValueError, 'Nonmonotone'): cc.validate_states(10, bad)

    def test_reject_line_change(self):
        bad = list(self.s10); bad[155] = (1, (1, 1))
        with self.assertRaisesRegex(ValueError, 'changed its projective line'):
            cc.validate_states(10, bad)

    def test_committed_result_reproduces(self):
        expected = json.loads((ROOT/'results/p401-cooperative-continuation/exact.json').read_text())
        actual = json.loads(json.dumps(self.certificate))
        self.assertEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
