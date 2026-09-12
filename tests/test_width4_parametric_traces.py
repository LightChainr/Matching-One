"""Mathematical controls: all-p identity, specialization and physical counts."""
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import width4_parametric_traces as s
from rank_two_onset import rank

class ParametricTraceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source=s.SOURCE
        if not source.exists():source=ROOT.parent/'source_inputs/width4-rank-closure-certificate.json'
        cls.data=s.verify(source)
        cls.definition=json.loads(s.DEFINITION.read_text())

    def test_every_parameter_polynomial_is_certified(self):
        c=self.data['certificate']
        self.assertEqual(c['zero_polynomials_checked_each'],126)
        self.assertTrue(c['resultant_sylvester_identity'])
        self.assertEqual(self.data['weighted_quotient']['counts'],[3,35,94,94])

    def test_independent_lift_census_at_rectangular_controls(self):
        for m in (2,3):
            adj=[]
            for y in range(m):
                for x in range(4):
                    adj.append([(((y+dy)%m)*4+(x+dx)%4,dx,dy)
                                for dx,dy in ((1,0),(-1,0),(0,1),(0,-1))])
            counts=[[0]*(4*m+1) for _ in range(3)]
            for mask in range(1<<(4*m)):counts[rank(mask,adj)][mask.bit_count()]+=1
            row=self.data['small_integer_sector_polynomials'][m-1]
            self.assertEqual(counts[0],row['rank0_coefficients'])
            self.assertEqual(counts[2],row['rank2_coefficients'])

    def test_exceptional_parameter_is_not_generic_order(self):
        orders={(r['probability'],r['channel']):r['minimum_scalar_order'] for r in self.data['minimality']}
        self.assertEqual(orders['1/2','M'],15)
        self.assertEqual(orders['2/3','M'],14)
        self.assertEqual(orders['3/4','M'],16)
        f=self.definition['factors']
        def at(name,t):return [s.polynomial_value(c,t) for c in f[name]['coefficients_desc_x_ascending_t']]
        self.assertEqual(at('closed_b',2),at('negative_square',2))
        self.assertEqual(at('closed_b',1),[1,0])

    def test_shared_modes_cancel_only_in_the_difference(self):
        weights=self.data['trace_factor_weights']
        for name in ('shared_3','shared_4a','shared_4b'):
            self.assertEqual(weights['P0'][name],weights['P2'][name])
            self.assertNotIn(name,weights['M'])
        self.assertEqual(self.data['generic_minimum_orders'],{'P0':17,'P2':23,'M':16})

    def test_diagonal_thermal_derivative_has_a_length_factor(self):
        # Exact coefficient of e in two diagonal eigenvalue powers.
        from math import comb
        lam=Fraction(3,5)
        for m in range(1,21):
            coefficient=comb(m,1)*lam**(m-1)-comb(m,1)*lam**(m-1)*(-1)
            self.assertEqual(coefficient,2*m*lam**(m-1))

if __name__=='__main__':unittest.main()
