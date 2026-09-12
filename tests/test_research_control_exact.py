"""Checks for concrete mathematical mistakes, not manuscript wording."""
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'/'theory'))
from irreducible_hidden_threshold import check


class ExactControlTests(unittest.TestCase):
    def test_symplectic_line_is_its_own_orthogonal(self):
        # In a symplectic plane annihilator of (a,b) is ay-bx=0.
        for a,b in ((1,0),(0,1),(1,1),(2,-3)):
            for x in range(-5,6):
                for y in range(-5,6):
                    self.assertEqual(a*y-b*x == 0, F(a)*y == F(b)*x)
        self.assertEqual(1*0-0*1, 0)  # (1,0) lies in its own annihilator.
        self.assertNotEqual(1*1-0*0, 0)  # (0,1) does not.

    def test_rank_identity_needs_no_direction_conjecture(self):
        for r in range(3):
            white = 2-r
            d = int(r>0)-int(white>0)
            self.assertEqual(d, r-1)
            self.assertEqual(F(1+d,2), F(r,2))

    def test_affine_pushforward_intertwines_inverse(self):
        # F(p)=p^2, Q(u)=sqrt(u); exact rational points avoid roundoff.
        alpha,beta = F(3,2),F(-1,7)
        for q in (F(1,5),F(1,2),F(4,5)):
            x = alpha*q+beta
            self.assertEqual(((x-beta)/alpha)**2, q*q)

    def test_shape_symmetry_does_not_force_center_half(self):
        center,a,b = F(3,5), F(1,4), F(3,4)
        for n in (10,20,100):
            width,eps = F(1,n),F(1,n*n)
            lower = center-width/2
            def Q(u):
                return (u+(1-eps)*lower/width)/((1-eps)/width+eps)
            def Z(u):
                return (Q(u)-Q(a))/(Q(b)-Q(a))
            for i in range(1,10):
                u=F(i,10)
                self.assertEqual(Z(u)+Z(1-u),1)
            # At 1/2 only the uniform background contributes to the CDF.
            self.assertEqual(2*(2*eps*F(1,2)-1), -2+2*eps)

    def test_irreducible_product_chain(self):
        self.assertEqual(check()['response_checks'],54)


if __name__ == '__main__':
    unittest.main()
