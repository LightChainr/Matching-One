"""Small mathematical controls, not document-wording or publication gates."""
import importlib.util
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal as D, localcontext
import unittest

FILE=Path(__file__).resolve().parents[1]/'scripts'/'cluster_sewing_identity.py'
spec=importlib.util.spec_from_file_location('sewing',FILE)
s=importlib.util.module_from_spec(spec); spec.loader.exec_module(s)

class SewingTests(unittest.TestCase):
    def test_column_boundary_all_small_masks(self):
        for matching in (False,True):
            for mask in range(64):
                c=s.vertices(mask,3,2)
                _,b=s.boundary_columns(c,3,matching)
                self.assertEqual(s.boundary_direct(c,3,matching),
                                 frozenset((i,y) for i,t in enumerate(b) for y in t))
    def test_halo_is_not_incidence(self):
        e=s.exact_marking_counterexample()
        self.assertEqual((e['external_void_sites'],e['open_to_void_incidences']),(8,12))
    def test_two_independent_winding_detectors(self):
        for matching in (False,True):
            for mask in range(256):
                c=s.vertices(mask,4,2)
                self.assertEqual(any(w for _,w in s.components(c,4,matching)),s.dsu_winding(c,4,matching))
    def test_straight_component_one_mark(self):
        c=frozenset((x,0) for x in range(5))
        for matching in (False,True):
            self.assertEqual(len(s.seam_marks(c,5,matching)),1)
    def test_random_guards_match_activity(self):
        for matching in (False,True):
            t,_=s.shape_table(3,1,matching)
            a=F(s.evaluate_table(t,F(1,3))['nu_truncated'])
            b,_=s.guard_enumeration(3,1,matching,F(1,3))
            self.assertEqual(a,b)
    def test_palm_reciprocal(self):
        t,_=s.shape_table(4,3,False)
        a=s.evaluate_table(t,F(1,2))
        self.assertEqual(F(a['nu_truncated']),F(9087,1048576))
        self.assertEqual(F(a['mark_mean_inverse_c']),1/F(a['component_mean_c']))
        self.assertGreater(F(a['wrong_component_reciprocal_estimate']),F(a['nu_truncated']))
    def test_two_row_trace(self):
        for matching in (False,True):
            t,_=s.shape_table(4,2,matching)
            p=F(1,4)
            a=F(s.evaluate_table(t,p)['nu_truncated'])
            b=s.two_row_transfer_trace(4,matching,p)-(p*(1-p)**2)**4
            self.assertEqual(a,b)
    def test_height_difference_monotonicity(self):
        a=[]
        for h in (1,2,3):
            t,_=s.shape_table(3,h,False)
            a.append(F(s.evaluate_table(t,F(1,2))['nu_truncated']))
        self.assertLess(a[0],a[1]); self.assertLess(a[1],a[2])
    def test_cyclic_gauge_telescopes(self):
        word=(1,3,2,3,1)
        gauge=F(1)
        for i,b in enumerate(word):
            a=word[i-1]; c=word[(i+1)%len(word)]
            gauge*=F(1+a+3*b,1+b+3*c)
        self.assertEqual(gauge,1)
    def test_all_width_contrast_annihilates_amplitude_mass(self):
        with localcontext() as ctx:
            ctx.prec=65
            for widths in ((2,4,8),(4,8,12),(3,11,17)):
                c=s.contrast_weights(widths)
                self.assertEqual(sum(c),0)
                self.assertEqual(sum(a*w for a,w in zip(c,widths)),0)
                y={w:D('1.7')-D('2.3')*w-D('.5')*D(w).ln() for w in widths}
                self.assertLess(abs(s.contrast(y,widths)-D('.5')),D('1e-55'))
    def test_returned_unequal_window_is_not_negative_seven(self):
        r=s.contrast_report()['models']
        self.assertAlmostEqual(float(r['NN']['correct_beta_2_4_8']),.71824578691523,12)
        self.assertAlmostEqual(float(r['matching']['correct_beta_2_4_8']),.53337743380313,12)
    def test_brownian_target_is_not_site_data(self):
        vals=[s.range_cdf(D(str(h))) for h in (.5,1,1.5,2,3)]
        self.assertTrue(all(0<a<1 for a in vals))
        self.assertEqual(vals,sorted(vals))
        self.assertAlmostEqual(float(vals[1]),.177923355643,11)
    def test_unequal_window_validation(self):
        with self.assertRaises(ValueError): s.contrast_weights((2,2,8))

if __name__=='__main__': unittest.main()
