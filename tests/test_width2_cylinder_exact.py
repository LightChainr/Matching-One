import sys
from pathlib import Path
from fractions import Fraction
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import width2_cylinder_exact as w

class WidthTwo(unittest.TestCase):
    def test_all_coefficients_from_independent_winding(self):
        for m in (2,3,4,5):
            self.assertEqual(w.matching_polynomial(m),w.enumerated_polynomial(m)[0])
    def test_known_l2_control(self):
        self.assertEqual(w.matching_polynomial(2),[-1,0,4,0,-2])
    def test_probability_trace_at_half(self):
        # At p=1/2, T=K/4. Check trace directly with integer powers.
        k=[[1,0,1],[0,1,1],[1,1,1]]
        a=[[int(i==j) for j in range(3)] for i in range(3)]
        for m in range(1,8):
            a=[[sum(a[i][l]*k[l][j] for l in range(3)) for j in range(3)] for i in range(3)]
            if m>=2:
                self.assertEqual(w.evaluate(w.matching_polynomial(m),Fraction(1,2)),
                                 Fraction(sum(a[i][i] for i in range(3))-3**m,4**m))

if __name__=='__main__':unittest.main()
