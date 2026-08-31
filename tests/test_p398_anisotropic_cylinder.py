from collections import defaultdict
from fractions import Fraction as F
from math import comb
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from p398_anisotropic_cylinder import derive_sector, evaluate, factored_certificate, signed_jordan_witness


class AnisotropicCylinderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.poly=derive_sector()

    def test_polynomial_and_isotropic_specialization(self):
        factored_certificate(self.poly)
        self.assertEqual([evaluate(self.poly[k],F(1,2),F(1,2)) for k in ("a","b","c","d")],
                         [F(1,16),-F(1,32),-F(1,64),F(1,32)])

    def test_critical_discriminant_and_signed_exceptional_point(self):
        actual=defaultdict(F)
        for (ph,pv),coefficient in self.poly["discriminant"].items():
            for j in range(pv+1):
                actual[ph+j]+=coefficient*(-1)**j*comb(pv,j)
        expected=defaultdict(F)
        # 4 h^2 (1-h)^8 (2-2h+h^2), strictly positive for 0<h<1.
        for j in range(9):
            for k,c in enumerate((2,-2,1)):
                expected[2+j+k]+=4*(-1)**j*comb(8,j)*c
        self.assertEqual({k:v for k,v in actual.items() if v},
                         {k:v for k,v in expected.items() if v})
        point=signed_jordan_witness(self.poly)
        self.assertEqual(point["eigenvalue"],"243/1024")
        self.assertTrue(point["nilpotent_square_zero"] and point["nonzero_nilpotent"])
        self.assertFalse(point["physical"])


if __name__=="__main__":
    unittest.main()
