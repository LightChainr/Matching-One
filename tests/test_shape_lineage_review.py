import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import shape_lineage_review as r

class ShapeReview(unittest.TestCase):
    def test_width_uses_difference_variance(self):
        c=np.eye(9)*2;c[1,7]=c[7,1]=.5
        self.assertAlmostEqual(r.width_se(c),np.sqrt(3))
    def test_delete_one_is_not_independent_samples(self):
        deleted=np.arange(100)/10000
        correct=np.sqrt(99/100*np.sum((deleted-deleted.mean())**2))
        wrong=np.std(deleted,ddof=1)/np.sqrt(100)
        self.assertAlmostEqual(correct/wrong,99)
    def test_shared_middle_scalar(self):
        v=[np.array([0.]),np.array([1.]),np.array([3.])]
        c=[np.eye(1)]*3;got=r.interval_norm_change(v,c)
        self.assertAlmostEqual(got['value'],1)
        self.assertAlmostEqual(got['se'],np.sqrt(6))
        self.assertAlmostEqual(got['old_independent_se'],2)
    def test_complex_step_width_and_covariance(self):
        q=np.linspace(.4,.7,9);c=np.eye(9)*1e-8
        _,cov,j=r.propagated(lambda x:np.array([x[7]-x[1]]),q,c)
        self.assertAlmostEqual(cov[0,0],2e-8)
        self.assertEqual(j[0,7],1);self.assertEqual(j[0,1],-1)
    def test_quadratic_symmetric_law_has_exact_normal_form(self):
        z=np.array([-.8,-.5,-.3,-.15,0,.15,.3,.5,.8])
        for scale in (.15,.07):
            t=.6+scale*z;q=t+.2*t*t
            self.assertLess(np.linalg.norm(r.normal_form_residual(q)),1e-12)
    def test_smooth_nonquadratic_normal_form_remainder_is_order_width_squared(self):
        z=np.array([-.8,-.5,-.3,-.15,0,.15,.3,.5,.8]);errors=[]
        for scale in (.05,.025):
            t=.6+scale*z;q=t+.2*t*t+.1*t**3
            errors.append(np.linalg.norm(r.normal_form_residual(q)))
        self.assertLess(abs(errors[0]/errors[1]-4),.02)
    def test_source_beta_uncertainty_enters_target(self):
        q=np.concatenate([np.linspace(.4,.8,9)]*3)
        j=r.derivative(r.common_quadratic_residual,q)
        self.assertGreater(np.linalg.norm(j[3:,0:9]),0)

if __name__=='__main__':unittest.main()
