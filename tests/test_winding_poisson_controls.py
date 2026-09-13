"""Small exact mathematical controls; not tests of document wording."""
from fractions import Fraction
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import winding_poisson_controls as wp


class WindingAnchors(unittest.TestCase):
    def graph(self, mask, w=3, m=7, matching=False):
        return wp.components(mask, w, m, wp.adjacency(w,m,matching))

    def test_two_components_not_rank_two(self):
        w,m,h=3,7,1
        mask=7 | (7 << (2*w))
        cc=self.graph(mask,w,m)
        self.assertEqual(wp.ambient_rank(cc),1)
        self.assertEqual(wp.global_anchors(cc,w,m,h),(0,2*w))

    def test_full_guard_closure_prevents_duplicate_anchor(self):
        # Adding one occupied neighbour kills the one-row anchor. An anchor
        # event is not increasing and must not be used directly in BK.
        self.assertEqual(wp.global_anchors(self.graph(7),3,7,1),(0,))
        self.assertEqual(wp.global_anchors(self.graph(7 | (1<<3)),3,7,1),())
        self.assertEqual(wp.global_anchors(self.graph(7 | (1<<3)),3,7,2),(0,))

    def test_vertical_seam_is_not_a_new_component(self):
        w,m,h=3,7,2
        mask=(7 << (6*w)) | 1
        cc=self.graph(mask,w,m)
        self.assertEqual(wp.global_anchors(cc,w,m,h),(6*w,))
        la=wp.adjacency(w,h+2,False,False)
        self.assertEqual(wp.local_all_anchors(mask,w,m,h,la),(6*w,))

    def test_matching_diagonal_winding(self):
        mask=(1<<0)|(1<<2)|(1<<4)
        self.assertEqual(wp.ambient_rank(self.graph(mask)),0)
        cc=self.graph(mask,matching=True)
        self.assertEqual(wp.ambient_rank(cc),1)
        self.assertEqual(wp.global_anchors(cc,3,7,2),(0,))

    def test_long_component_excluded(self):
        mask=7 | (1<<3) | (1<<6)
        cc=self.graph(mask)
        self.assertEqual(wp.global_anchors(cc,3,7,2),())
        self.assertEqual(wp.global_anchors(cc,3,7,3),(0,))

    def test_one_row_exact_intensity(self):
        for w in (2,3):
            for matching in (False,True):
                cs=wp.local_intensity_coefficients(w,1,matching)
                self.assertEqual(cs[w],1)
                self.assertEqual(sum(cs),1)
                p=Fraction(2,7)
                self.assertEqual(wp.bernstein_value(cs,p),p**w*(1-p)**(2*w))

    def test_two_row_independent_exact_coefficients(self):
        self.assertEqual(wp.local_intensity_coefficients(2,2,False),[0,0,1,6,6,0,0,0,0])
        self.assertEqual(wp.local_intensity_coefficients(2,2,True),[0,0,3,6,2,0,0,0,0])

    def test_horizontal_tiebreak_not_uniform(self):
        # The leftmost label is a convention, not a translation invariant site
        # mark. For a full row only column zero is the chosen representative.
        self.assertEqual(wp.global_anchors(self.graph(7),3,7,1),(0,))

    def test_closed_guard_positive_contact(self):
        d=wp.exact_row_contact_control()
        self.assertGreater(Fraction(d['two_anchors_two_rows_apart']),
                           Fraction(d['product_of_marginals']))
        self.assertEqual(Fraction(d['ratio']),Fraction(64,27))

    def test_exact_finite_law_and_mean(self):
        t=wp.enumerate_one(2,5,1,False,(Fraction(1,4),))
        c=t['cases'][0]
        self.assertEqual(Fraction(c['lambda']),Fraction(405,4096))
        self.assertEqual(t['local_global_anchor_failures'],0)
        self.assertGreater(Fraction(c['void_mismatch_probability']),0)
        self.assertLessEqual(c['poisson_tv_float_diagnostic'],float(Fraction(c['agg_process_tv_bound'])))

    def test_common_labels_are_not_finitely_independent(self):
        d=wp.two_colour_control()
        self.assertEqual(Fraction(d['count_covariance']),Fraction(4831,1048576))
        self.assertLessEqual(Fraction(d['opposite_monotone_joint']),
                             Fraction(d['opposite_monotone_product']))
        self.assertGreater(Fraction(d['count_covariance']),0)

    def test_neighbourhood_uses_whole_windows(self):
        w,m,h=3,12,2
        dd=wp.dependency(w,m,h)
        self.assertIn(0,dd[0])
        self.assertIn(3*w,dd[0])
        self.assertNotIn(5*w,dd[0])
        self.assertIn(11*w,dd[0])
        for i in range(w*m):
            self.assertLessEqual(len(dd[i]),w*(2*h+3))

    def test_component_activity_uses_distinct_boundary_sites(self):
        for matching in (False, True):
            self.assertEqual(wp.component_activity_counts(2,1,matching),{(2,4):1})
            self.assertEqual(wp.component_activity_counts(3,1,matching),{(3,6):1})

    def test_component_and_full_window_derivatives_agree(self):
        for matching in (False,True):
            d=wp.component_activity_control(3,2,matching)
            self.assertTrue(all(c['two_differentiation_routes_agree'] for c in d['cases']))
            for c in d['cases']:
                self.assertGreaterEqual(Fraction(c['log_nu_second_derivative']),
                                        Fraction(c['curvature_lower_bound']))

    def test_one_row_log_curvature_is_contact_term(self):
        d=wp.component_activity_control(2,1,False)
        for c in d['cases']:
            self.assertEqual(c['log_nu_second_derivative'],c['curvature_lower_bound'])

    def test_invalid_geometry_rejected(self):
        with self.assertRaises(ValueError):
            wp.adjacency(1,8)
        with self.assertRaises(ValueError):
            wp.local_all_anchors(0,3,4,2,wp.adjacency(3,4,False,False))


if __name__ == '__main__':
    unittest.main()
