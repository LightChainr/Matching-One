import sys
from pathlib import Path
import unittest
from math import log,exp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from shape_transport_order import composition_difference

class TransportOrder(unittest.TestCase):
    def test_affine_maps_need_not_commute(self):
        a=lambda x:2*x+1;b=lambda x:3*x-2
        for p in (-1,0,1):self.assertEqual(composition_difference(a,b,p),-4)
    def test_common_nonlinear_conjugacy_preserves_sign_not_value(self):
        a=lambda x:exp(2*log(x)+1);b=lambda x:exp(3*log(x)-2)
        vals=[composition_difference(a,b,p) for p in (.5,1.,2.)]
        self.assertTrue(all(v<0 for v in vals));self.assertNotEqual(vals[0],vals[-1])
    def test_crossing_orientation_preserving_maps_cannot_be_common_affine(self):
        a=lambda x:2*x;b=lambda x:x+x**3
        self.assertGreater(composition_difference(a,b,-.5),0)
        self.assertLess(composition_difference(a,b,.5),0)

if __name__=='__main__':unittest.main()
