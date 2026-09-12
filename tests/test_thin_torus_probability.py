"""Mathematical checks: wrong tail constant, collapsed births, or lost modes."""
from __future__ import annotations
from fractions import Fraction
from itertools import product
import json
from math import comb
import os
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import thin_torus_extremes as geo
import thin_torus_spectral_controls as spec
import width4_topological_source_spectrum as source_spec
from mpmath import mp

DEFINITION=Path(os.environ.get('MATCHING_PARAMETRIC_DEFINITION',str(spec.DEFAULT_DEFINITION)))
CERTIFICATE=Path(os.environ.get('MATCHING_WIDTH4_CERTIFICATE',str(ROOT/'results/research-control-20260912/width4-rank-closure-certificate.json')))


def rational_trace(co,m):
    """Newton recurrence independent of the numerical companion-matrix path."""
    n=len(co)-1
    result=[Fraction(n)]
    for k in range(1,m+1):
        if k<=n: result.append(-sum(co[i]*result[k-i] for i in range(1,k))-k*co[k])
        else: result.append(-sum(co[i]*result[k-i] for i in range(1,n+1)))
    return result[m]


class ThinTorusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # These are explicit dependencies, not a reason to silently skip validation.
        cls.definition=json.loads(DEFINITION.read_text())

    def test_central_trinomial_is_a_shape_count(self):
        for w in range(2,8):
            counted=sum(sum(word)==0 for word in product((-1,0,1),repeat=w))
            self.assertEqual(geo.central_trinomial(w),counted)
        self.assertEqual([geo.central_trinomial(w) for w in (2,3,4)],[3,7,19])

    def test_minimal_matching_cycles_are_not_only_empty_rows(self):
        self.assertEqual(geo.motif_control(3)['minimal_cycle_sets'],49)

    def test_full_and_empty_rows_force_rank_one(self):
        report=geo.barrier_control(2,3)
        self.assertTrue(report['rank_one_in_every_case'])
        self.assertGreater(report['full_and_empty_row_cases'],0)

    def test_birth_subset_DP_against_explicit_permutations(self):
        self.assertTrue(geo.brute_permutation_control()['DP_matches_explicit_permutations'])

    def test_finite_births_are_not_claimed_independent(self):
        report=geo.first_two_birth_counts(2,3)
        self.assertEqual(Fraction(report['continuous_birth_covariance']),Fraction(1123,58800))
        self.assertEqual(sum(c for _,_,c in report['joint_K1_K2']),720)

    def test_endpoint_coefficients_are_exact(self):
        for w in (2,3,4):
            rows=spec.WidthSpectrum(w,self.definition).endpoint_certificates()
            self.assertEqual([r['first_nonzero_degree'] for r in rows],[w,w])
            self.assertEqual([r['first_coefficient'] for r in rows],[1,geo.central_trinomial(w)])

    def test_all_finite_sector_weights_against_physical_graphs(self):
        for w,m in ((2,3),(3,2),(4,2)):
            model=spec.WidthSpectrum(w,self.definition)
            n=w*m
            ranks=[geo.ambient_rank(mask,w,m) for mask in range(1 << n)]
            for p in (Fraction(1,3),Fraction(2,3)):
                measured={r:sum(p**mask.bit_count()*(1-p)**(n-mask.bit_count())
                                for mask,rank in enumerate(ranks) if rank==r) for r in (0,1,2)}
                for channel,expected in (('P0',measured[0]),('P2',measured[2]),
                                         ('M',measured[2]-measured[0])):
                    got=sum(c*rational_trace(model.coefficients(name,p),m)
                            for name,c in model.weights[channel].items())
                    self.assertEqual(got,expected)

    def test_bulk_rank_CDF_flattens_not_sharpens(self):
        with mp.workdps(80):
            model=spec.WidthSpectrum(2)
            for p in (mp.mpf('.25'),mp.mpf('.75')):
                self.assertLess(abs(model.cdf(512,p)-mp.mpf('.5')),mp.mpf('1e-12'))

    def test_conditional_logistic_is_not_unconditional_CDF(self):
        with mp.workdps(80):
            model=spec.WidthSpectrum(2)
            q,l,a0,a2=model.branch_data()
            p=q+mp.mpf('.5')/256
            p0=model.value('P0',256,p);p2=model.value('P2',256,p)
            conditional=p2/(p0+p2)
            target=1/(1+mp.exp(-(a2-a0)/2))
            self.assertLess(abs(conditional-target),mp.mpf('.002'))
            self.assertGreater(conditional,mp.mpf('.8'))
            self.assertLess(abs(model.cdf(256,p)-mp.mpf('.5')),mp.mpf('1e-35'))

    def test_rank_source_formula_has_correct_normalizer(self):
        p0,p2,z=Fraction(1,8),Fraction(1,4),Fraction(2)
        p1=1-p0-p2;M=p2-p0;E=p0+p2
        cosh,sinh=(z+1/z)/2,(z-1/z)/2
        Z=p0/z+p1+p2*z
        self.assertEqual(Z,1+M*sinh+E*(cosh-1))
        self.assertEqual((p2*z-p0/z)/Z,(M*cosh+E*sinh)/(1+M*sinh+E*(cosh-1)))

    def test_source_exposes_28_and_29_modes(self):
        report=source_spec.verify(CERTIFICATE,DEFINITION)
        orders={r['channel']:r['generic_minimum_scalar_order'] for r in report['channels']}
        self.assertEqual(orders,{'M':16,'E':28,'twice_tilted_mean_numerator':28,'twice_Z_at_exp_s_2':29})
        for row in report['channels']:
            self.assertNotEqual(int(row['exact_Hankel_determinant']),0)
            self.assertTrue(any(row['independent_modular_determinant_checks'].values()))

    def test_one_event_probability_is_not_normal_root_precision(self):
        # At activity 1/4, 11 snapshots suffice for one event with probability .95,
        # while ten do not. This is not a Gaussian root-precision calculation.
        self.assertGreater(Fraction(3,4)**10,Fraction(1,20))
        self.assertLess(Fraction(3,4)**11,Fraction(1,20))

if __name__=='__main__':unittest.main()
