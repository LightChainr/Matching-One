"""Concrete topology/count errors these controls prevent; no stochastic tests."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from lifted_boundary_rank import (State, initial_state, advance, close_rank, torus_rank,
    state_for_rows, reachable_states, rank_count_polynomials, minimize_rank_machine)
from verify_width4_rank_closure import graph_rank, minimality_signature
from width4_half_probability_spectrum import report as spectral_report


class BoundaryRankTests(unittest.TestCase):
    def test_parallel_periodic_edges_are_not_collapsed(self):
        self.assertEqual(torus_rank(2,[3,0]),1)
        self.assertEqual(torus_rank(2,[1,1]),1)
        self.assertEqual(torus_rank(2,[3,3]),2)

    def test_old_partition_only_memory_collision(self):
        a,b=state_for_rows(4,[13,5]),state_for_rows(4,[7,5])
        self.assertNotEqual(a,b)
        self.assertEqual((torus_rank(4,[13,5,13,0]),torus_rank(4,[7,5,13,0])),(0,1))
        for rows in ([13,5,13,0],[7,5,13,0],[1,5,4,5],[11,14]):
            self.assertEqual(torus_rank(4,rows),graph_rank(4,rows))

    def test_forgotten_horizontal_cycle_remains_in_global_span(self):
        state=state_for_rows(4,[0,15,0])
        self.assertTrue(state.horizontal)
        self.assertTrue(all(r<0 for r,_ in state.entries))
        self.assertEqual(close_rank(state),1)

    def test_common_shear_preserves_all_short_future_ranks(self):
        # Switch one representative by a deliberately large shear. The four
        # occupancy labels are fixed; the lift integer is not a new site.
        state=state_for_rows(4,[3,6,12],shear=False)
        entries=tuple((r,g+17*(int(j>=4)-int(r>=4))) if r>=0 else (r,g)
                      for j,(r,g) in enumerate(state.entries))
        shifted=State(state.horizontal,entries)
        self.assertNotEqual(state,shifted)
        for a in range(16):
            for b in range(16):
                old=advance(advance(state,a,False),b,False)
                new=advance(advance(shifted,a,False),b,False)
                self.assertEqual(close_rank(old),close_rank(new))

    def test_complete_row_closure_and_continuation_separation(self):
        states,transition,initial=reachable_states(4)
        labels,reps,q,rank,steps=minimize_rank_machine(states,transition)
        self.assertEqual((len(states),sum(map(len,transition))),(1448,23168))
        self.assertEqual(steps,[3,164,509,509])
        self.assertEqual(len(set(minimality_signature(q,rank,i) for i in range(len(q)))),509)

    def test_exact_four_square_polynomial(self):
        counts,_=rank_count_polynomials(4,4)
        self.assertEqual(list(map(sum,counts)),[36559,19932,9045])
        self.assertEqual([b-a for a,b in zip(counts[0],counts[2])],
            [-1,-16,-120,-560,-1812,-4272,-7448,-9424,-7874,-2896,1720,2832,1660,560,120,16,1])

    def test_spectral_identity_has_a_physical_tail_certificate(self):
        data=spectral_report()
        self.assertEqual(data['minimal_scalar_recurrence_order'],15)
        self.assertNotEqual(data['hankel_determinant_15_start_m2'],0)
        self.assertEqual(data['annihilator_zero_scalar_moments_verified'],509)
        # The scalar recurrence is NOT a claim p(A)c=0 as a whole vector.
        self.assertEqual(data['pA_c_nonzero_entries'],82)
        counts,_=rank_count_polynomials(4,4)
        self.assertEqual(data['first_31_s_m'][3],sum(counts[2])-sum(counts[0]))


if __name__=='__main__':unittest.main()
