import os
from pathlib import Path
import sys
import unittest
from fractions import Fraction
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import width4_site_sources as ws

CERTIFICATE = Path(os.environ.get('MATCHING_ONE_RANK_CERTIFICATE', str(ws.DEFAULT_CERTIFICATE)))

class SiteSourcesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = ws.load_certificate(CERTIFICATE)

    def test_all_profile_quotients_and_group_actions(self):
        r = ws.profile_report(self.source)
        self.assertEqual([x['classes'] for x in r['profiles'].values()], list(ws.EXPECTED_COUNTS))
        self.assertTrue(all(x['equals_orbit_partition_blockwise'] for x in r['profiles'].values()))

    def test_seven_profiles_exhaust_column_partitions(self):
        words = [(0,)]
        for _ in range(3):
            words = [word+(k,) for word in words for k in range(max(word)+2)]
        def groups(word):
            return tuple(tuple(j for j,v in enumerate(word) if v==k) for k in range(max(word)+1))
        def canon(gs):
            return tuple(sorted(tuple(sorted(g)) for g in gs))
        representatives = {min(canon([[pi[j] for j in g] for g in groups(word)]) for pi in ws.D4)
                           for word in words}
        expected = {min(canon([[pi[j] for j in g] for g in gs]) for pi in ws.D4)
                    for gs in ws.PROFILES.values()}
        self.assertEqual(len(words), 15)
        self.assertEqual(len(representatives), 7)
        self.assertEqual(representatives, expected)

    def test_first_odd_response_vanishes(self):
        for p in (Fraction(1, 3), Fraction(1, 2), Fraction(3, 5)):
            self.assertEqual(ws.dipole_response(self.source, 5, (2,), p), 0)

    def test_finite_amplitude_is_exact_not_an_extrapolation(self):
        for p, a, b in ((Fraction(1, 2), Fraction(1, 3), Fraction(1, 4)),
                        (Fraction(2, 5), Fraction(1, 10), Fraction(1, 8))):
            value, _ = ws.four_sign(self.source, 4, 1, 2, a, b, p)
            self.assertEqual(value, ws.dipole_response(self.source, 4, (1, 2), p))
        self.assertEqual(ws.dipole_response_half(self.source, 4, (1, 2)), Fraction(327, 1024))

    def test_homogeneous_block_does_not_carry_local_input(self):
        s = self.source
        labels = ws.refine(s, ws.PROFILES['homogeneous'])['labels']
        i, j = ws.follow(s, [1, 1]), ws.follow(s, [2, 2])
        self.assertEqual(labels[i], labels[j])
        p,q = Fraction(1, 3), Fraction(3, 5)
        weights, den = ws.row_integer_weights([q,p,p,p])
        out=s['quotient_rank_output'];t=s['quotient_transitions']
        difference=Fraction(sum(w*(out[t[i][b]]-out[t[j][b]]) for b,w in enumerate(weights)),den)
        self.assertEqual(difference, q-p)

    def test_invalid_inputs_rejected(self):
        with self.assertRaises(ValueError):ws.validate_groups(((0,1),(1,2,3)))
        with self.assertRaises(ValueError):ws.row_integer_weights([Fraction(5,4)]*4)
        with self.assertRaises(ValueError):ws.dipole_response_half(self.source,4,(1,1))

if __name__ == '__main__':unittest.main()
