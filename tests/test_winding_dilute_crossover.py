import importlib.util
from pathlib import Path
from fractions import Fraction as F
from math import comb, factorial
import unittest

P=Path(__file__).resolve().parents[1]/'scripts/winding_dilute_crossover.py'
spec=importlib.util.spec_from_file_location('dilute',P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class DiluteTests(unittest.TestCase):
    def test_physical_winding_and_no_winding(self):
        self.assertEqual(m.winding_components({(x,0) for x in range(6)},6),1)
        self.assertEqual(m.winding_components({(x,0) for x in range(5)},6),0)
        self.assertEqual(m.winding_components({(x,y) for x in range(6) for y in (0,2)},6),2)
    def test_cyclic_gap_counts(self):
        for w in range(6,13):
            for k in range(0,min(4,w)+1):
                self.assertEqual(len(list(m.gap_subsets(w,k))),m.gap_count(w,k))
    def test_collision_union_bound(self):
        for w in range(6,65):
            for k in range(0,8):
                self.assertGreaterEqual(F(m.gap_count(w,k)*factorial(k),w**k),
                                        1-F(5*k*(k-1),2*w))
    def test_cycle_geometry(self):
        x=m.cycle_geometry_checks(12)
        self.assertGreater(x['cycles'],100)
        self.assertLessEqual(x['max_contact_over_r'],16)
    def test_first_nonrow_cycle_count(self):
        for row in m.minimal_nonrow_census(6):
            w=row['width'];self.assertEqual(row['nonrow_winding_counts'][w+2],w*(w-3))
    def test_retained_exact_small_width_series(self):
        for w in (3,4):
            a=m.rational_series(w,10)
            self.assertEqual(a[w],1);self.assertEqual(a[w+1],0)
            self.assertEqual(a[w+2],w*(w-3))
        self.assertEqual(m.evaluate_small_nu(2,F(1,2)),F(7,48))
        self.assertEqual(m.evaluate_small_nu(4,F(1,2)),F(323849,5576960))
    def test_finite_lower_probability_coefficients(self):
        for w in (6,8,12):
            p=F(1,128)
            lower=m.finite_cycle_lower(w,p)
            self.assertGreater(lower,0)
            self.assertLess(lower,2)
            self.assertGreater(m.universal_lower_factor(w,p),0)
    def test_matrix_trace_is_independent_of_logdet_recurrence(self):
        self.assertEqual(m.matrix_loop_coefficients(10),m.direct_matrix_loops(10))
    def test_matrix_diffusion_poisson_equation(self):
        P=((F(3,4),F(1,4)),(F(1,4),F(3,4)))
        Q1=tuple(tuple(P[i][j]*F(1 if j==0 else -1,2) for j in range(2)) for i in range(2))
        h=(F(1,2),F(-1,2));mu=F(3,2)
        for i in range(2):
            self.assertEqual(h[i]-sum(P[i][j]*h[j] for j in range(2)),sum(Q1[i]))
        sigma=F(1,2)+sum(Q1[i][j]*h[j] for i,j in ((0,0),(0,1),(1,0),(1,1)))
        self.assertEqual(sigma,1);self.assertEqual(sigma/mu,F(2,3))
    def test_mass_and_amplitude_cancellation(self):
        # exact nu_w=A*t^w/w for beta=1 gives R=4/3.
        for w in (4,7,12):
            f=lambda n:F(7,5)*F(2,5)**n/n
            self.assertEqual(f(w)*f(3*w)/f(2*w)**2,F(4,3))
            # An exact log-density c/w correction has coefficient c/(3w).
            self.assertEqual(F(1,w)+F(1,3*w)-2*F(1,2*w),F(1,3*w))
    def test_bessel_control(self):
        import mpmath as mp
        with mp.workdps(50):
            self.assertLess(m.bessel_contrast('0.1',mp),mp.mpf('.07'))
            self.assertGreater(m.bessel_contrast('0.5',mp),mp.mpf('.60'))
            x=mp.mpf('.7')
            series=mp.fsum(x**(2*r)/factorial(r)**2 for r in range(50))
            self.assertLess(abs(series-mp.besseli(0,2*x)),mp.mpf('1e-45'))
            moment=mp.fsum(r*r*x**(2*r)/factorial(r)**2 for r in range(50))
            self.assertLess(abs(moment-x*x*series),mp.mpf('1e-45'))
    def test_matching_minimal_pattern_coefficient(self):
        from itertools import product
        for w in range(2,8):
            n=sum(sum(steps)==0 for steps in product((-1,0,1),repeat=w))
            self.assertEqual(n,m.central_trinomial(w))
        self.assertEqual([m.central_trinomial(w) for w in (4,8,12)],[19,1107,73789])
    def test_uniform_walk_bound_numerically(self):
        import mpmath as mp
        with mp.workdps(45):
            for w,p in ((6,F(1,128)),(64,F(1,256))):
                a=m.walk_upper_normalized(w,p,mp)
                lam=w*m.mp_number(p,mp);b=mp.besseli(0,2*lam)
                self.assertLessEqual(a/b,mp.exp(6*w*m.mp_number(p,mp)**2))
                self.assertLessEqual(m.mp_number(m.finite_cycle_lower(w,p),mp),a)

if __name__=='__main__': unittest.main()
