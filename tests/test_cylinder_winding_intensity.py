import sys
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import cylinder_winding_intensity as c

class CylinderIntensityTests(unittest.TestCase):
    def test_complete_closure_counts(self):
        for w,n,k in ((2,6,3),(3,14,4),(4,38,7)):
            for matching in (False,True):
                states,table=c.build_transfer(w,matching)
                red,blocks=c.reward_lump(table)
                self.assertEqual((len(states),len(red)),(n,k))
                self.assertTrue(all(len(row)==1<<w for row in table))
                self.assertTrue(all(row[0][0]==0 for row in table))

    def test_path_information_not_dropped(self):
        _,table=c.build_transfer(4)
        def bits(rows): return sum(m<<(4*j) for j,m in enumerate(rows))
        self.assertEqual(c.reward_count(bits([13,5,13,0]),4,4,table),0)
        self.assertEqual(c.reward_count(bits([7,5,13,0]),4,4,table),1)

    def test_count_components_not_rows(self):
        for w in (2,3,4):
            _,table=c.build_transfer(w)
            full=(1<<w)-1
            mask=full+(full<<w)+(full<<(3*w))
            self.assertEqual(c.reward_count(mask,w,4,table),2)

    def test_independent_small_open_graph(self):
        for w,h in ((2,3),(3,2)):
            for matching in (False,True):
                _,table=c.build_transfer(w,matching)
                for mask in range(1<<(w*h)):
                    expected=c.graph_components(mask,w,h,matching)[1]
                    self.assertEqual(c.reward_count(mask,w,h,table),expected)

    def test_stochastic_not_configurationwise_lumping(self):
        for matching in (False,True):
            _,table=c.build_transfer(4,matching)
            red,blocks=c.reward_lump(table)
            self.assertEqual(c.reward_histogram(4,5,table),c.reward_histogram(4,5,red))
            signatures={}
            for i,row in enumerate(table):
                sig=Counter((mask.bit_count(),r,blocks[j]) for mask,(j,r) in enumerate(row))
                if blocks[i] in signatures:self.assertEqual(signatures[blocks[i]],sig)
                signatures[blocks[i]]=sig

    def test_exact_density_and_variance(self):
        expected={2:(F(7,48),F(343,6912)),3:(F(169,1984),F(4769011,122023936)),
                  4:(F(323849,5576960),F(186754153229427053,6098108298338304000))}
        for w,pair in expected.items():
            _,table=c.build_transfer(w)
            red,_=c.reward_lump(table)
            out=c.stationary_reward(red,F(1,2))
            self.assertEqual((out['mean'],out['variance_rate']),pair)

    def test_unreduced_stationary_check(self):
        _,full=c.build_transfer(3)
        red,_=c.reward_lump(full)
        a=c.stationary_reward(full,F(1,3));b=c.stationary_reward(red,F(1,3))
        self.assertEqual(a['mean'],b['mean'])
        self.assertEqual(a['variance_rate'],b['variance_rate'])

    def test_matching_density_and_variance_duality(self):
        for w in (2,3,4):
            _,a=c.build_transfer(w);_,b=c.build_transfer(w,True)
            a,_=c.reward_lump(a);b,_=c.reward_lump(b)
            x=c.stationary_reward(a,F(1,4));y=c.stationary_reward(b,F(3,4))
            self.assertEqual(x['mean'],y['mean'])
            self.assertEqual(x['variance_rate'],y['variance_rate'])

    def test_topological_count_identity(self):
        full=(1<<9)-1
        for mask in range(1<<9):
            r,a=c.graph_components(mask,3,3,False,True)
            s,b=c.graph_components(full^mask,3,3,True,True)
            self.assertEqual(r+s,2)
            self.assertEqual(a-b,r-1)
            self.assertNotEqual((a,b),(0,0))

    def test_residual_certificate(self):
        _,table=c.build_transfer(4);red,_=c.reward_lump(table)
        exact=c.stationary_reward(red,F(1,4))
        approx=[F(int(x*1024),1024) for x in exact['stationary']]
        approx[0]+=1-sum(approx)
        bound=c.stationary_certificate(red,F(1,4),approx)
        self.assertLessEqual(bound['lower'],exact['mean'])
        self.assertGreaterEqual(bound['upper'],exact['mean'])
        self.assertGreater(bound['stationarity_l1_residual'],0)

    def test_invalid_input_is_rejected(self):
        with self.assertRaises(ValueError):c.empty_state(1)
        with self.assertRaises(RuntimeError):c.build_transfer(4,state_cap=3)
        with self.assertRaises(ValueError):c.solve_fraction([[F(0)]],[F(1)])
        _,table=c.build_transfer(2)
        with self.assertRaises(ValueError):c.stationary_reward(table,F(1))
        with self.assertRaises(ValueError):c.stationary_certificate(table,F(1,2),[F(0)]*len(table))

    def test_renewal_coefficients_and_prefactor(self):
        self.assertEqual([c.central_trinomial(n) for n in range(8)],[1,1,3,7,19,51,141,393])
        out=c.renewal_controls()['controls']
        self.assertLess(abs(out[-1]['effective_beta']-.5),.002)
        self.assertLess(abs(out[-1]['normalized_amplitude']-1),.002)
        self.assertTrue(all(x['effective_beta']<y['effective_beta'] for x,y in zip(out,out[1:])))

if __name__=='__main__':unittest.main()
