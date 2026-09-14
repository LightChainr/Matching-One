from fractions import Fraction as F
from pathlib import Path
import importlib.util
import unittest

P=Path(__file__).resolve().parents[1]/'scripts'/'query_cavity_controls.py'
spec=importlib.util.spec_from_file_location('cav',P)
cav=importlib.util.module_from_spec(spec);spec.loader.exec_module(cav)

class QueryCavityTests(unittest.TestCase):
    def test_wired_all_black_has_hidden_centre(self):
        self.assertEqual(cav.wired_statistics(0),(0,8,1))
    def test_hidden_centre_colour_is_not_queried(self):
        i=list(cav.product(range(-1,2),repeat=2)).index((0,0))
        self.assertEqual(cav.wired_statistics(1<<i),(0,8,1))
    def test_all_white(self): self.assertEqual(cav.wired_statistics(511),(9,0,0))
    def test_all_wired_fraction_moments(self):
        for q in (F(1,3),F(1,2),F(3,4),F(9,10)):
            self.assertGreater(F(cav.shield_control(q)['conditional_residual']['fraction']),0)
    def test_positive_polynomial_factorization(self):
        for p in (F(1,100),F(1,4),F(1,2),F(9,10),F(99,100)):
            poly=sum(F(c)*p**i for i,c in enumerate(cav.P_COEFF))
            self.assertEqual((1-p**8)*(9-p**8)-64*(1-p)*p**7,(1-p)**2*poly)
    def test_independent_nonparallel_increment_vectors(self):
        # Ring all black versus one open gate, with the centre respectively black/white.
        sites=list(cav.product(range(-1,2),repeat=2))
        gate=1<<sites.index((1,0)); centre=1<<sites.index((0,0))
        base=cav.wired_statistics(0); one=cav.wired_statistics(gate)
        two=cav.wired_statistics(gate|centre)
        a=(one[0]-base[0],one[1]-base[1])
        b=(two[0]-base[0],two[1]-base[1])
        self.assertEqual(a[0]*b[1]-a[1]*b[0],-1)
    def test_small_cylinder_cages(self):
        r=cav.animal_control();self.assertEqual(r['checked_classes'],3792)
    def test_winding_not_generic_graph_cycle(self):
        self.assertTrue(cav.winds({(x,0) for x in range(8)},8))
        self.assertFalse(cav.winds(cav.KING,8))
    def test_double_ring_defect_count(self):
        black=set(cav.KING)|{((x+4)%8,y) for x,y in cav.KING}
        self.assertEqual(len(cav.hidden_sites(black,8)),2)
    def test_distinct_rare_patterns_do_not_duplicate(self):
        self.assertEqual(cav.pattern_control()['least_distinct_union'],12)
    def test_singleton_and_dimer_boundaries(self):
        self.assertEqual(len(cav.boundary_of({(0,0)})),8)
        self.assertEqual(len(cav.boundary_of({(0,0),(1,0)})),10)
        self.assertEqual(len(cav.boundary_of({(0,0),(1,1)})),12)
    def test_boolean_low_density_coefficients(self):
        self.assertEqual(cav.sparse_defect_series()['J_coefficients'],{8:1,9:-2,10:9,11:-20})
    def test_geometric_limit_normalization(self):
        r=cav.limiting_laws();self.assertEqual(F(r['variance_H']['fraction']),72)
        self.assertEqual(F(r['correlation_E_H_squared']['fraction']),F(8,9))
    def test_rational_clock_limit(self):
        r=cav.independent_block_control()
        for j in range(3):
            errors=[F(r[k+j]['absolute_difference']['fraction']) for k in (0,3,6)]
            self.assertLessEqual(errors[2],errors[1]); self.assertLessEqual(errors[1],errors[0])
    def test_invalid_probability(self):
        with self.assertRaises(ValueError):cav.shield_control(F(0))

if __name__=='__main__':unittest.main()
