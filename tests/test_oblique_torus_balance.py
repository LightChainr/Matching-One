import sys, unittest
from pathlib import Path
from fractions import Fraction as F
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import oblique_torus_balance as o

class TestObliqueGeometry(unittest.TestCase):
    def test_gaussian_basis(self):
        for h,n in [((10,3,1),10),((13,5,1),13),((4225,268,1),4225)]:
            t=o.Torus(*h);u,v=t.basis()
            self.assertEqual(o.norm2(u),n);self.assertEqual(o.det(u,v),t.n)
    def test_hnf_identifications(self):
        t=o.Torus(5,2,3)
        for x in range(-10,11):
            for y in range(-7,8):
                self.assertEqual(t.reduce(x,y),t.reduce(x+5,y))
                self.assertEqual(t.reduce(x,y),t.reduce(x+2,y+3))
    def test_nonprimitive_ambient_shortest(self):
        u,v=o.Torus(64,0,80).basis();self.assertEqual(o.norm2(u),64**2)
    def test_slabs_disjoint(self):
        for h in ((65,17,130),(4225,268,1)):
            t=o.Torus(*h);g=o.slab_geometry(t,o.MATCHING,coarse=True)
            nodes=[x for b in g['slabs'] for x in b['nodes']]
            self.assertEqual(len(nodes),len(set(nodes)))
            self.assertGreaterEqual(g['count']*g['systole_squared'],4*t.n)
    def test_all_and_empty_crossings(self):
        t=o.Torus(13,5,1);g=o.slab_geometry(t,o.MATCHING)
        self.assertTrue(all(o.has_crossing((1<<t.n)-1,b) for b in g['slabs']))
        self.assertTrue(all(not o.has_crossing(0,b) for b in g['slabs']))
    def test_exact_sign_certificate(self):
        t=o.Torus(4096,113,4096)
        self.assertTrue(o.sign_certificate(t,4,F(1,10))['strict_sign_certified'])
        self.assertFalse(o.sign_certificate(t,4,F(1,2))['strict_sign_certified'])
    def test_exhaustive_tilted_small(self):
        r=o.check_small(o.Torus(10,3,1));self.assertTrue(r['complementary_rank_identity'])
    def test_invalid_geometry(self):
        with self.assertRaises(ValueError):o.reduced_basis((1,0),(2,0))
        with self.assertRaises(ValueError):o.Torus(2,2,1)

if __name__=='__main__':unittest.main()
