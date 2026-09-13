from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import tagged_winding_span as T
import tagged_span_controls as C

class TaggedSpanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cache={}
        for g in (False,True):
            for w in (2,3,4):
                states,tr,src=T.build(w,g); rr,src,blocks=T.lump(tr,src)
                cls.cache[(w,g)]=(states,tr,rr,src)

    def calc(self,w,g,p=F(1,2),bins=8):
        _,_,tr,src=self.cache[(w,g)]
        return T.moments(tr,src,p,bins=bins)

    def test_known_density_and_full_moments_nn2(self):
        r=self.calc(2,False)
        self.assertEqual(r['nu'],F(7,48));self.assertEqual(r['mean'],F(76,21))
        self.assertEqual(r['variance'],F(2068,441));self.assertEqual(r['cv2'],F(517,1444))

    def test_equal_density_does_not_equal_shape(self):
        a=self.calc(2,False);b=self.calc(2,True)
        self.assertEqual(a['nu'],b['nu']);self.assertEqual(b['mean'],F(100,21))
        self.assertEqual(b['variance'],F(5596,441));self.assertNotEqual(a['d_h'],b['d_h'])

    def test_isolated_row_all_controls(self):
        for g in (False,True):
            for w in (2,3,4):
                for p in (F(1,4),F(1,2),F(3,4)):
                    self.assertEqual(self.calc(w,g,p,1)['d_h'][0],p**w*(1-p)**(2*w))

    def test_physical_height_three(self):
        r=self.calc(4,False,bins=3)
        self.assertEqual(sum(r['d_h']),F(9087,1048576))
        self.assertEqual(r['nu'],F(323849,5576960))

    def test_all_source_words(self):
        for g in (False,True):
            n,mass=C.enumerate_tag_paths(2,2,g)
            self.assertEqual(n,256);self.assertEqual(mass,self.calc(2,g,bins=2)['d_h'][1])

    def test_late_merge_tie_break(self):
        # Two candidates at time zero eventually join one winding component.
        paths=[x for x in T.source_paths(4,False) if x[:2]==(0,5)]
        accepted=[]
        for prev,mask,anchor,s in paths:
            outcome=0
            for mask in (7,15,0):
                s,outcome=T.step(s,mask)
                if outcome:break
            if outcome==1:accepted.append(anchor)
        self.assertEqual(accepted,[0])

    def test_past_contact_rejects_new_tag(self):
        paths=[x for x in T.source_paths(4,False) if x[:2]==(4,5)]
        self.assertEqual(len(paths),1)
        _,_,_,s=paths[0];s,outcome=T.step(s,7)
        self.assertEqual(outcome,-1)

    def test_empty_mask_terminates(self):
        for states,tr,_,_ in self.cache.values():
            self.assertTrue(all(row[0][0]<0 for row in tr))
            self.assertTrue(all(s.colours.count(T.TAG)==1 for s in states))

    def test_exact_tail_moments(self):
        r=self.calc(4,False,bins=3)
        self.assertEqual(sum(r['d_h'])+r['tail'],r['nu'])
        self.assertEqual(sum((i+1)*x for i,x in enumerate(r['d_h']))+r['tail_first'],r['nu']*r['mean'])
        self.assertEqual(sum((i+1)**2*x for i,x in enumerate(r['d_h']))+r['tail_second'],r['nu']*r['second'])

    def test_residual_intervals_cover_exact(self):
        _,_,tr,src=self.cache[(4,False)];r=T.certified_moments(tr,src,F(1,4));ex=self.calc(4,False,F(1,4))
        for interval,key in [('density_interval','nu'),('mean_interval','mean'),('variance_interval','variance'),('cv2_interval','cv2')]:
            lo,hi=r[interval];self.assertLessEqual(lo,ex[key]);self.assertGreaterEqual(hi,ex[key])

    def test_two_parameter_activity_polynomials(self):
        for g in (False,True):
            st,tr,src=T.activity_transfer(3,g);cs=C.activity_coefficients(tr,src,3)
            for h,c in enumerate(cs,1):
                physical,_=C.shape_activity_counts(3,h,g);self.assertEqual(c,physical)

    def test_activity_and_tag_resolvents_agree(self):
        for g in (False,True):
            _,tr,src=T.activity_transfer(4,g);r=T.activity_joint_moments(tr,src,F(1,4));e=self.calc(4,g,F(1,4))
            self.assertEqual(r['nu'],e['nu']);self.assertEqual(r['mean_span'],e['mean']);self.assertEqual(r['variance_span'],e['variance'])

    def test_joint_moment_psd(self):
        for g in (False,True):
            _,tr,src=T.activity_transfer(3,g);r=T.activity_joint_moments(tr,src,F(1,4))
            self.assertLessEqual(r['covariance_span_occupation']**2,r['variance_span']*r['variance_occupation'])
            self.assertGreaterEqual(r['mean_occupation'],3)

    def test_exact_palm_samples(self):
        for g,p in ((False,F(1,4)),(True,F(1,8))):
            for seed in (0,1):
                r=T.sample_component(4,g,p,seed);C.verify_sample(4,g,r)
                self.assertEqual(r['path_probability'],r['unconditioned_row_weight']/r['nu'])

    def test_range_not_total_variation(self):
        y=[0,1]*50+[0];range_y=max(y)-min(y);variation=sum(abs(a-b) for a,b in zip(y,y[1:]))
        self.assertEqual(range_y,1);self.assertEqual(variation,100)
        # Decorations within B of a skeleton change its range by at most 2B.
        for B in (0,1,2):
            full=[v+d for v in y for d in range(-B,B+1)]
            self.assertLessEqual((max(full)-min(full))-range_y,2*B)

    def test_parameter_boundary_is_rejected(self):
        _,_,tr,src=self.cache[(2,False)]
        for p in (F(0),F(1),F(-1),F(3,2)):
            with self.assertRaises(ValueError):T.numeric_system(tr,src,p)

    def test_symbolic_pgf(self):
        r=C.symbolic_width_two();self.assertEqual(len(r),2)
        self.assertIn('p**2',r[0]['density_pgf'])

if __name__=='__main__':unittest.main()
