import importlib.util
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

P=Path(__file__).resolve().parents[1]/'scripts/giant_cluster_random_information.py'
spec=importlib.util.spec_from_file_location('score_field_tests',P)
m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m)

class ScoreFieldTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.report=m.run()
    def test_wired_multilinear(self):
        self.assertEqual(sum(r['zero_multilinear_moments'] for r in self.report['wired_controls']),314)
    def test_wired_normalization(self):
        self.assertTrue(all(r['finite_tilt_normalizer']['fraction']=='1' for r in self.report['wired_controls']))
    def test_source_engine_pinned(self):
        self.assertEqual(self.report['engine_blob'],m.ENGINE_BLOB)
    def test_physical_coefficients(self):
        self.assertEqual(sum(r['complete_shapes'] for r in self.report['physical_controls']),2956)
    def test_spatial_quotient_not_uniform(self):
        self.assertEqual([r['colored_states'] for r in self.report['physical_controls'][-2:]],[23,23])
    def test_off_uniform_duality(self):
        self.assertEqual(len(self.report['inhomogeneous_duality']),6)
    def test_exact_cross_covariance(self):
        self.assertTrue(all(r['cross_covariance']['fraction']=='0' for r in self.report['all_height_scores']))
    def test_shared_exposure_not_independence(self):
        self.assertEqual(self.report['limit_model']['orthogonal_square_correlation']['fraction'],'1/5')
    def test_student_pooling(self):
        self.assertEqual(self.report['limit_model']['pooled_estimator_limit_df'],{'1':2,'2':4,'3':6,'4':8})
    def test_variance_gamma_inverse(self):
        self.assertEqual(self.report['limit_model']['pooled_limit_variance_times_c']['3']['fraction'],'1/2')
    def test_reject_bad_width(self):
        with self.assertRaises(ValueError):m.colored_activity(m.load_engine(),9)
    def test_inverse(self):
        G=m.inverse_fraction([[2,1],[1,2]])
        self.assertEqual(G,[[F(2,3),F(-1,3)],[F(-1,3),F(2,3)]])

if __name__=='__main__': unittest.main()
