from fractions import Fraction as F
import json
from pathlib import Path
import sys,unittest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import conditional_odds_integration as ti
import invisible_spatial_marks as im
import verify_source_hessian as verify
DATA=ROOT/'results/research-control-20260912'

class MarksAndIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.h=json.loads((DATA/'full-site-hessian.json').read_text())
        cls.cert=verify.read_cert(DATA/'width4-rank-closure-certificate.json')
        cls.ranks=verify.all_ranks(cls.cert)
    def test_orbit_sizes(self):
        self.assertEqual([len(im.orbit(x)) for x in (51,275,39)],[16,128,64])
    def test_all_orbit_support_is_rank0_K4(self):
        for shape in (51,275,39):
            for mask in im.orbit(shape):
                self.assertEqual(self.ranks[mask],0);self.assertEqual(mask.bit_count(),4)
    def test_mark_cellwise_invisibility(self):
        for shape,fa,den in ((275,8,16),(39,4,8)):
            mark={x:fa for x in im.orbit(51)};mark.update({x:-1 for x in im.orbit(shape)})
            self.assertEqual(sum(mark.values()),0)
            self.assertGreaterEqual(min(1-F(v,den) for v in mark.values()),F(1,2))
            self.assertGreaterEqual(min(1+F(v,den) for v in mark.values()),F(1,2))
    def test_two_mark_rank_certificate(self):
        matrix=[]
        for k in ((2,0),(2,2)):
            field=verify.pattern(*k);row=[]
            for shape,fa in ((275,8),(39,4)):
                mark={x:fa for x in im.orbit(51)};mark.update({x:-1 for x in im.orbit(shape)})
                row.append(sum(v*sum(field[i] for i in range(16) if (x>>i)&1)**2 for x,v in mark.items()))
            matrix.append(row)
        self.assertEqual(matrix,[[-512,-128],[0,-256]])
        self.assertEqual(matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0],131072)
    def test_L_checkerboard_null_finite_amplitude(self):
        mark={x:8 for x in im.orbit(51)};mark.update({x:-1 for x in im.orbit(275)})
        for p in (F(1,3),F(3,5)):
            self.assertEqual(im.weighted_field_sum(mark,verify.pattern(2,2),p,F(1,20)),0)
    def test_conditional_integration_onset(self):
        s=self.h['rank_bernstein_counts']
        self.assertEqual(ti.normalize_sector(s[0])[:2],(0,1))
        self.assertEqual(ti.normalize_sector(s[2])[:2],(7,16))
    def test_conditional_mean_identity_exact(self):
        a,b=self.h['rank_bernstein_counts'][0],self.h['rank_bernstein_counts'][2]
        ka,aa,ra=ti.normalize_sector(a);kb,ab,rb=ti.normalize_sector(b)
        for t in (F(1,2),F(1),F(2)):
            mu0,_=ti.occupation_moments(a,t);mu2,_=ti.occupation_moments(b,t)
            derivative=ti.value(ti.deriv(rb),t)/ti.value(rb,t)-ti.value(ti.deriv(ra),t)/ti.value(ra,t)
            self.assertEqual(mu2-mu0,(kb-ka)+t*derivative)
    def test_anchor_bound_majorizes_exact_polynomials(self):
        a,b=self.h['rank_bernstein_counts'][0],self.h['rank_bernstein_counts'][2]
        _,_,ra=ti.normalize_sector(a);_,_,rb=ti.normalize_sector(b)
        for t in (F(1,10**4),F(1,10**8)):
            bound=ti.anchor_bound(16,7,16,t)
            self.assertLessEqual(ti.value(ra,t)-1,bound)
            self.assertLessEqual(ti.value(rb,t)-1,bound)
    def test_source_fields_unit_RMS(self):
        for k in ((1,0),(2,0),(1,1),(2,1),(2,2)):
            field=verify.pattern(*k)
            self.assertEqual(sum(field),0);self.assertEqual(sum(x*x for x in field),16)

if __name__=='__main__':unittest.main()
