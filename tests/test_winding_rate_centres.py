"""Mathematical controls, not tests of document wording."""
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import winding_rate_centres as w

class WindingRateControls(unittest.TestCase):
    def test_seed_enumeration_against_independent_formulas(self):
        for p in (F(1,100),F(1,5),F(1,4),F(1,2),F(4,5),F(1)):
            self.assertEqual(w.bernstein_count_value(w.seed_counts('NN'),p),
                             p*(p+(1-p)*(2*p**3-p**6)))
            self.assertEqual(w.bernstein_count_value(w.seed_counts('matching'),p),
                             p*(1-(1-p)**3))

    def test_seed_dyadic_signs(self):
        for kind in ('NN','matching'):
            counts = w.seed_counts(kind)
            lo,hi = w.seed_root_interval(counts)
            self.assertEqual(hi-lo,F(1,2**64))
            self.assertLess(w.bernstein_count_value(counts,lo),F(1,16))
            self.assertGreater(w.bernstein_count_value(counts,hi),F(1,16))
            self.assertLess(hi,F(1,4))

    def test_all_tiny_cut_witnesses(self):
        for kind in ('NN','matching'):
            data=w.census_control(3,3,kind,ring_check=True)
            self.assertEqual(data['configurations'],512)
            self.assertEqual(data['first_span_cut_witnesses_checked'],
                             data['nonzero_winding_configurations'])

    def test_seam_ring_is_sufficient_not_identical(self):
        for kind in ('NN','matching'):
            # The top full row winds, but the seed requires the middle source.
            top=7
            self.assertIsNotNone(w.winding_walk(top,3,3,kind))
            self.assertFalse(w.seed_ring(top,3,3,kind))
            middle=7<<3
            self.assertTrue(w.seed_ring(middle,3,3,kind))
            self.assertIsNotNone(w.winding_walk(middle,3,3,kind))

    def test_parallel_edges_keep_short_winding(self):
        # Circumference two has distinct lifted edges between the same vertices.
        for kind in ('NN','matching'):
            walk=w.winding_walk(3,2,3,kind)
            self.assertIsNotNone(walk)
            w.verify_cut(3,2,3,kind,walk)

    def test_normalized_seed_counts_not_double_binomial(self):
        for kind in ('NN','matching'):
            counts=w.seed_counts(kind)
            self.assertEqual(w.bernstein_count_value(counts,F(1,2)),F(sum(counts),256))
            self.assertEqual(w.bernstein_count_value(counts,F(0)),0)
            self.assertEqual(w.bernstein_count_value(counts,F(1)),1)

    def test_finite_cluster_likelihood_bound(self):
        data=w.seed_and_cluster_controls()
        for kind in ('NN','matching'):
            row=data[kind]['cluster_comparison']
            lhs=F(int(row['lhs']['numerator']),int(row['lhs']['denominator']))
            rhs=F(int(row['rhs']['numerator']),int(row['rhs']['denominator']))
            self.assertLessEqual(lhs,rhs)

    def test_backtracking_first_span(self):
        path=[(0,0),(1,0),(0,0),(0,1),(1,1),(2,1),(3,1),(4,1),(4,0)]
        prefix,origin,axis=w.cut_witness(path,4)
        self.assertEqual(axis,'x')
        self.assertEqual(prefix[-1],(3,1))
        self.assertEqual(origin,(0,0))

    def test_invalid_arguments(self):
        with self.assertRaises(ValueError): w.steps('bond')
        with self.assertRaises(ValueError): w.bernstein_count_value([1,1],F(3,2))
        with self.assertRaises(ValueError): w.census_control(4,3,'NN')
        with self.assertRaises(ValueError): w.seed_ring(0,2,3,'NN')

if __name__=='__main__': unittest.main()
