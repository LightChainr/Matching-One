import sys
import unittest
from fractions import Fraction as F
from pathlib import Path
from itertools import product

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import dilute_fragmentation_clock as d

class FragmentationTests(unittest.TestCase):
    def test_minimal_classification(self):
        r=d.physical_minimal_control()
        self.assertEqual(r['configurations'],6196)
        self.assertEqual(r['winding_sets'],73)
        self.assertEqual(r['hole_centres'],12)

    def test_matching_minimal_rings(self):
        ms=d.motifs(4,True)
        self.assertEqual(sum(m.kind=='barrier' for m in ms),19)
        for m in ms:
            self.assertEqual(d.dsu_winds(m.sites,4,True),m.kind=='barrier')

    def test_nn_minimal_rings(self):
        ms=d.motifs(8,False)
        self.assertEqual(sum(m.kind=='barrier' for m in ms),1)
        for m in ms:
            self.assertEqual(d.dsu_winds(m.sites,8,False),m.kind=='barrier')

    def test_no_same_order_overlap(self):
        for w,m in ((4,True),(8,False)):
            r=d.dependency_polynomials(w,m)
            self.assertGreater(r['smallest_union'],w)

    def test_dependency_polynomials(self):
        self.assertEqual(d.dependency_polynomials(4,True)['b2'],{'5':152,'6':428,'7':712})
        self.assertEqual(d.dependency_polynomials(8,False)['b2'],{'12':32,'13':64,'14':112,'15':32})

    def test_column_control_four(self):
        r=d.intensities(4,True,(F(2),F(1),F(1,2),F(1)))
        self.assertEqual(r['barrier'],19)
        self.assertEqual(r['ratio'],F(25,76))

    def test_column_control_eight(self):
        r=d.intensities(8,False,(F(2),F(1),F(1),F(1),F(1,2),F(1),F(1),F(1)))
        self.assertEqual(r['barrier'],1)
        self.assertEqual(r['ratio'],F(45,2))

    def test_opposite_pairs_can_be_invisible_at_four(self):
        r=d.intensities(4,True,(F(2),F(1,2),F(2),F(1,2)))
        self.assertEqual(r['barrier'],19)
        self.assertEqual(r['ratio'],F(4,19))

    def test_two_time_integral(self):
        for k,r,z,v in product(map(F,[1,2,5]),(F(4,19),F(8)),(F(0),F(1,2),F(1)),(F(0),F(1,3),F(1))):
            self.assertEqual(d.two_time(k,r,z,v,F(1,2),F(1,4)),
                             d.two_time_by_integrals(k,r,z,v,F(1,2),F(1,4)))

    def test_two_time_marginals(self):
        for k,r,z in product(map(F,[1,2,5]),(F(4,19),F(8)),(F(0),F(1,2),F(1))):
            self.assertEqual(d.two_time(k,r,z,F(1)),1/(1+r*(1-z))**2)
            self.assertEqual(d.two_time(k,r,F(1),z),1/(1+r*(1-z))**2)
            self.assertEqual(d.two_time(F(1),r,z,z),1/(1+r*(1-z*z))**2)

    def test_mark_covariance(self):
        for r,k in product((F(4,19),F(8)),map(F,[1,2,4])):
            self.assertEqual(d.mark_moments(k,r)['correlation'],1/k)

    def test_generator_stationarity(self):
        self.assertTrue(all(d.generator_moments(n)==0 for n in range(1,10)))

    def test_time_reverse(self):
        for m,n in product(range(7),repeat=2):
            self.assertEqual(*d.adjoint_moment(m,n))

    def test_common_labels(self):
        self.assertEqual(d.motif_multitime_check(),610)

    def test_bad_parameters(self):
        with self.assertRaises(ValueError):d.motifs(3,True)
        with self.assertRaises(ValueError):d.two_time(F(1,2),F(1),F(1),F(1))
        with self.assertRaises(ValueError):d.intensities(4,True,(F(0),)*4)

if __name__=='__main__':unittest.main()
