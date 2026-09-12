"""Small checks preventing wrong finite polynomials and missing frontier memory."""
from fractions import Fraction
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import width3_cylinder_exact as w


class WidthThreeExactTests(unittest.TestCase):
    def test_existing_three_by_three_bernstein_rung(self):
        self.assertEqual(w.bernstein_counts(3),[-1,-9,-36,-78,-90,-36,36,36,9,1])

    def test_small_configurationwise_rank_dictionary(self):
        for m in (2,3):
            self.assertEqual(w.check_census(m)['classification_failures'],0)

    def test_full_seven_state_trace_matches_block_recurrence(self):
        for m in range(2,9):
            self.assertEqual(w.overlap_trace_half(3,m),w.evaluate(w.power_coefficients(m),Fraction(1,2)))

    def test_all_length_positive_defect_certificate(self):
        c=w.interval_certificate()
        self.assertTrue(c['all_m_at_least_2_defect_positive'])
        self.assertTrue(c['hprime_positive'])
        self.assertEqual(len(c['eigenvalue_intervals']),5)

    def test_width_four_overlap_can_both_overcount_and_undercount_rank(self):
        a,b=w.width4_witnesses()
        self.assertEqual((a['actual_ambient_rank'],a['naive_rank']),(0,1))
        self.assertEqual((b['actual_ambient_rank'],b['naive_rank']),(2,1))

    def test_equal_frontier_masks_do_not_preserve_continuation(self):
        witness=w.continuation_memory_witness()
        self.assertEqual(witness['completed_torus_ranks'],[0,1])
        self.assertEqual(witness['completed_matching_D'],[-1,0])

    def test_ordinary_partition_and_current_rank_still_lose_topological_memory(self):
        witness=w.topology_memory_witness()
        self.assertEqual(witness['completed_torus_ranks'],[0,1])
        signatures=witness['prefix_signatures']
        self.assertEqual(signatures[0]['ordinary_frontier_partition'],signatures[1]['ordinary_frontier_partition'])
        self.assertEqual(signatures[0]['occupied_count'],5)
        self.assertEqual(signatures[1]['occupied_count'],5)


if __name__=='__main__':unittest.main()
