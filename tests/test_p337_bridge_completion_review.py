"""Focused exact checks; no Monte Carlo or repository production imports."""
import json
import sys
import unittest
from fractions import Fraction as F
from itertools import product
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import p337_bridge_completion_review as m


class BridgeCompletionTests(unittest.TestCase):
    def test_bell_counts(self):
        self.assertEqual(len(m.partitions(4)),15)
        self.assertEqual(len(m.partitions(8)),4140)

    def test_invalid_partition_size(self):
        for n in (0,-1,True,1.5):
            with self.assertRaises(ValueError):
                m.partitions(n)

    def test_invalid_patterns(self):
        for pi in ((0,1,1),(-1,0,1,1),(0,2,2,0)):
            with self.assertRaises(ValueError):
                m.endpoint_signature(pi)
        with self.assertRaises(ValueError):
            m.bridge_endpoints((0,1,2,3,0,4,5,6))

    def test_endpoint_count_and_sharp_envelope(self):
        ps=m.distinguished_patterns()
        self.assertEqual(len(ps),62)
        self.assertEqual(set(m.endpoint_signature(p) for p in ps),
                         {F(-1,2),F(-1,4),F(0),F(1,2),F(3,4),F(1),F(3,2)})
        self.assertEqual(m.endpoint_signature((0,2,1,2)),F(3,2))

    def test_existing_two_component_witness(self):
        self.assertEqual(m.endpoint_signature((0,0,1,1)),F(-1,4))
        self.assertEqual(m.pair_activation((0,0,1,1,0,0,1,1)),F(1,16))

    def test_all_1922_two_bridge_factorizations(self):
        n=0
        for pi in m.partitions(8):
            if len(set(pi[:4]) & set(pi[4:]))==2:
                x,y=m.bridge_endpoints(pi)
                self.assertEqual(m.pair_activation(pi),m.endpoint_signature(x)*m.endpoint_signature(y),pi)
                n+=1
        self.assertEqual(n,1922)

    def test_bridge_label_exchange(self):
        for p in m.distinguished_patterns():
            swap=tuple(1-x if x in (0,1) else x for x in p)
            self.assertEqual(m.endpoint_signature(p),m.endpoint_signature(swap))

    def test_c4_rotation(self):
        for p in m.distinguished_patterns():
            self.assertEqual(m.endpoint_signature(p),m.endpoint_signature(p[1:]+p[:1]))
            self.assertEqual(m.endpoint_signature(p),m.endpoint_signature((p[0],p[3],p[2],p[1])))

    def test_one_colour_zero_for_every_completion(self):
        for a in (F(-10),F(0),F(1,2),F(5,3)):
            self.assertEqual(m.kernel_at_one((0,0,0,0),a),0)

    def test_literal_integer_two_colour_factorization(self):
        for q in (4,5,6):
            for a in (F(-1),F(0),F(1,2)):
                row=m.integer_two_bridge_check(q,(0,1,2,3),(0,2,1,2),a)
                self.assertEqual(F(row['residual']),0)

    def test_literal_four_wire_matches_independent_rational_norm(self):
        for q in (4,5,6):
            qf=F(q)
            n2=qf*(q-3)*(3*q*q-9*q+8)/(8*(q-1)*(q-2))
            n0=F(2*q*q-4*q+3,2*q*(q-1))
            cross=F(q-3,4*(q-1))
            for a in (F(-1),F(0),F(1,2),F(1)):
                c=1+a*(q-1)
                expected=n2+2*c*cross+c*c*n0
                actual=sum((m.kernel_at_integer(q,tuple(v),a)**2 for v in product(range(q),repeat=4)),F(0))
                self.assertEqual(actual,expected,(q,a))

    def test_bivariate_activation(self):
        pi=(0,1,2,3,0,1,2,3)
        for a in (F(-3),F(0),F(1,2),F(1),F(5,3)):
            for b in (F(-3),F(0),F(1,2),F(1),F(5,3)):
                self.assertEqual(m.pair_activation(pi,a,b),F(3,2)+(a-F(1,2))*(b-F(1,2))/2)
                self.assertEqual(m.pair_activation((0,1,2,3,0,3,2,1),a,b),m.pair_activation(pi,a,b))

    def test_same_completion_sharp_minimum(self):
        pi=(0,1,2,3,0,1,2,3)
        self.assertEqual(m.pair_activation(pi,F(1,2),F(1,2)),F(3,2))
        self.assertEqual(m.pair_activation(pi),F(13,8))
        self.assertGreater(m.pair_activation(pi,F(10),F(10)),F(3,2))

    def test_different_completion_is_not_positive_theorem(self):
        self.assertLess(m.pair_activation((0,1,2,3,0,1,2,3),F(-3),F(3)),0)

    def test_gram_schur_complement(self):
        g00,g01,g11=F(13,8),F(-1,4),F(1,2)
        self.assertEqual(g00*g11-g01*g01,F(3,4))
        self.assertEqual(g00-g01*g01/g11,F(3,2))

    def test_signed_two_bridge_is_not_event_probability(self):
        # The exact two-bridge kernel can be negative despite genuine bridges.
        pi=(0,0,1,1,0,1,2,3)
        self.assertEqual(m.pair_activation(pi),F(-1,8))

    def test_certificate_replays(self):
        report=m.build_report()
        stored=json.loads((ROOT/'results/p337-p0-bridge-review/exact.json').read_text())
        self.assertEqual(json.loads(json.dumps(report)),stored)

if __name__=='__main__':
    unittest.main()
