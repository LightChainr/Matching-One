"""Independent exact algebra, input rejection, and normalized-response checks."""
import copy
import csv
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('p337_anti_review',ROOT/'scripts/p337_antisymmetric_trace_review.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def projector(Q):
    pairs=[(i,j) for i in range(Q) for j in range(Q)]
    J=lambda i,k:F(int(i==k))-F(1,Q)
    return pairs,[[F(1,2)*(J(i,k)*J(j,l)-J(i,l)*J(j,k))
                   for k,l in pairs] for i,j in pairs]


def matvec(A,x):
    return [sum(a*b for a,b in zip(row,x)) for row in A]


def decimal_from_fraction(x):
    return D(x.numerator)/D(x.denominator)

class AntisymmetricReviewTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=m.build_result()
        cls.ordinary,cls.trace,cls.p,cls.sources=m.load_inputs()

    def test_saved_exact_certificate_reproduces(self):
        self.assertEqual(self.result,json.loads((m.DATA/'result.json').read_text()))

    def test_antisymmetric_activation_certified_negative(self):
        v=self.result['responses']['antisymmetric_activation']['V_over_A']
        self.assertLess(F(v['hi']),0)
        self.assertLess(F(v['hi'])-F(v['lo']),F(1,10**35))

    def test_symmetric_control_matches_published_value(self):
        self.assertAlmostEqual(self.result['responses']['symmetric_control']['V_display_only'],
                               -0.001904836180602413,places=16)

    def test_linearity_and_two_connection_signs(self):
        r=self.result['responses']
        va=r['type_A']['V_display_only'];vb=r['type_B']['V_display_only']
        self.assertGreater(va,0);self.assertLess(vb,0)
        self.assertAlmostEqual(r['antisymmetric_activation']['V_display_only'],(vb-va)/2,places=16)
        self.assertAlmostEqual(r['symmetric_control']['V_display_only'],-va-vb,places=16)

    def test_all_colour_character_checks(self):
        for Q in range(4,10):
            self.assertEqual(m.character_coefficients(Q),
                {'constant':0,'one_port':0,'A':1,'B':-1,'norm':1})

    def test_projector_idempotent_and_closure_traces(self):
        for Q in (4,5):
            pairs,P=projector(Q);n=len(pairs)
            self.assertEqual([[sum(P[i][k]*P[k][j] for k in range(n))
                               for j in range(n)] for i in range(n)],P)
            d=F((Q-1)*(Q-2),2)
            self.assertEqual(sum(P[i][i] for i in range(n)),d)
            self.assertEqual(sum(P[i][pairs.index((b,a))] for i,(a,b) in enumerate(pairs)),-d)

    def test_even_endpoint_vector_annihilated(self):
        pairs,P=projector(4)
        for i in range(4):
            for j in range(4):
                v=[int((a,b)==(i,j))+int((a,b)==(j,i)) for a,b in pairs]
                self.assertEqual(matvec(P,v),[0]*len(pairs))

    def test_nonzero_odd_standard_endpoint_survives(self):
        pairs,P=projector(4)
        a=[F(1),F(0),F(-1),F(0)];b=[F(0),F(1),F(-1),F(0)]
        v=[a[i]*b[j]-b[i]*a[j] for i,j in pairs]
        self.assertTrue(any(v));self.assertEqual(matvec(P,v),v)

    def test_even_operator_is_not_even_endpoint_vector(self):
        pairs,P=projector(4);swap=[pairs.index((j,i)) for i,j in pairs]
        self.assertEqual([[P[swap[i]][swap[j]] for j in range(16)] for i in range(16)],P)
        self.assertEqual(sum(P[i][i] for i in range(16)),3)

    def test_colour_endpoint_zero_and_fixed_derivative(self):
        dim=lambda Q:F((Q-1)*(Q-2),2)
        self.assertEqual(dim(1),0)
        dprime=F(2*1-3,2)
        self.assertEqual(dprime/F(1)**2,F(-1,2))
        self.assertEqual(-dprime/F(1),F(1,2))
        self.assertEqual(dim(4)/16,F(3,16))
        self.assertEqual(-dim(4)/4,F(-3,4))

    def test_prefactor_invariance_of_vanishing_packet(self):
        # For F(1)=0, (cF)'/(cZ)=F'/Z at Q1 for any nonzero c(1).
        fprime,z=F(7,11),F(13,17)
        for c0,cprime in ((F(2),F(9)),(F(1,3),F(-25,2))):
            self.assertEqual((cprime*0+c0*fprime)/(c0*z),fprime/z)

    def test_bad_digest_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            target=Path(folder)/'data';shutil.copytree(m.DATA,target)
            p=target/'inputs/trace_types.csv';p.write_text(p.read_text().replace('axis,A,10,17,5','axis,A,10,17,6'))
            with self.assertRaisesRegex(ValueError,'digest mismatch'):
                m.load_inputs(target)

    def test_bad_root_bracket_rejected_after_digest_update(self):
        with tempfile.TemporaryDirectory() as folder:
            target=Path(folder)/'data';shutil.copytree(m.DATA,target)
            p=target/'inputs/root_bracket.json';r=json.loads(p.read_text());r['root_interval']={'lo':'1/10','hi':'1/5'}
            p.write_text(json.dumps(r))
            s=target/'SOURCES.json';a=json.loads(s.read_text());a['input_sha256'][p.name]=hashlib.sha256(p.read_bytes()).hexdigest();s.write_text(json.dumps(a))
            with self.assertRaisesRegex(ValueError,'root bracket'):
                m.build_result(target)

    def test_interval_rejects_zero_divisor(self):
        with self.assertRaises(ZeroDivisionError):
            1/m.Interval(-1,1)

    def test_bernstein_normalization(self):
        for rows in self.ordinary.values():
            self.assertEqual(m.polynomial([r['count'] for r in rows]),[1]+[0]*25)

    def test_independent_bernstein_finite_difference_at_comoving_root(self):
        # Direct Bernstein sums and nonlinear normalized ratios, not power-basis jets.
        # A linearized root path suffices for the first derivative; no new root search.
        with localcontext() as ctx:
            ctx.prec=80
            p=decimal_from_fraction((self.p.lo+self.p.hi)/2)
            eps=D('1e-17')
            def data_at(p):
                ans=[]
                for g in ('axis','tilted'):
                    rows=self.ordinary[g];c=[F(0)]*26
                    for row in self.trace:
                        if row['geometry']==g:
                            c[row['k']]+=row['count']*(F(-1,2) if row['type']=='A' else F(1,2))
                    w=[p**k*(1-p)**(25-k) for k in range(26)]
                    wp=[w[k]*(D(k)/p-D(25-k)/(1-p)) for k in range(26)]
                    q=sum(D(r['sum_q'])*w[k] for k,r in enumerate(rows))
                    qp=sum(D(r['sum_q'])*wp[k] for k,r in enumerate(rows))
                    e=sum(D(r['sum_e'])*w[k] for k,r in enumerate(rows))
                    ep=sum(D(r['sum_e'])*wp[k] for k,r in enumerate(rows))
                    b=sum(decimal_from_fraction(c[k])*w[k] for k in range(26))
                    bp=sum(decimal_from_fraction(c[k])*wp[k] for k in range(26))
                    ans.append((q,qp,e,ep,b,bp))
                return ans
            base=data_at(p);Dp=sum(a[1] for a in base)/2
            pdot=sum(a[0]*a[4] for a in base)/(2*Dp)
            def u_over_A(s):
                ans=data_at(p+s*pdot);qp=[];ep=[]
                for q,q0p,e,e0p,b,bp in ans:
                    Z=1+s*b
                    qp.append(q0p/Z-q*s*bp/Z**2)
                    ep.append(e0p/Z-e*s*bp/Z**2)
                return ((ep[0]-ep[1])/decimal_from_fraction(m.DELTA))/(sum(qp)/2)
            observed=(u_over_A(eps)-u_over_A(-eps))/(2*eps)
            cert=self.result['responses']['antisymmetric_activation']['V_over_A']
            expected=decimal_from_fraction((F(cert['lo'])+F(cert['hi']))/2)
            self.assertLess(abs(observed-expected),D('1e-30'))

if __name__=='__main__':
    unittest.main()
