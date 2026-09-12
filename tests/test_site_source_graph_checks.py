import os
from pathlib import Path
import sys
import unittest
from fractions import Fraction
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import width4_site_sources as ws
import site_source_graph_checks as gc
CERTIFICATE=Path(os.environ.get('MATCHING_ONE_RANK_CERTIFICATE',str(ws.DEFAULT_CERTIFICATE)))

class PhysicalSourceTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.source=ws.load_certificate(CERTIFICATE)

    def test_spiral_is_1D_not_rank2(self):
        rank, images=gc.lifted_graph([3,6,12,9])
        self.assertEqual(rank,1)
        self.assertTrue(any(x and y for im in images for x,y in im['cycle_gains']))

    def test_existing_topological_memory_witnesses(self):
        for word,r in (([13,5,13,0],0),([7,5,13,0],1),([0,7,5,7],0),([0,13,5,7],1),([11,14],2)):
            rank,images=gc.lifted_graph(word)
            self.assertEqual(rank,r)
            self.assertTrue(all(im['saturated'] for im in images))

    def test_full_physical_polynomial_and_root_certificate(self):
        result=gc.physical_controls(self.source)
        self.assertEqual(result['physical_configurations_checked'],69888)
        self.assertEqual(result['tables']['4']['mixed_at_half'][1],'327/1024')
        self.assertEqual(result['tables']['4']['mixed_at_half'][2],'633/2048')
        for row in result['root_source_responses']:
            lo,hi=map(Fraction,row['mixed_root_interval'])
            self.assertLess(lo,hi)
            self.assertLess(hi,0)

    def test_bernstein_derivative_and_interval(self):
        # p^2+(1-p)^2 represented without binomial factors.
        c=[1,0,1]
        self.assertEqual(gc.bernstein_derivative(c),[-2,2])
        lo,hi=gc.bernstein_interval(c,Fraction(1,3),Fraction(1,2))
        for p in (Fraction(1,3),Fraction(2,5),Fraction(1,2)):
            v=gc.unnormalized_bernstein(c,p)
            self.assertLessEqual(lo,v);self.assertLessEqual(v,hi)

if __name__ == '__main__':unittest.main()
