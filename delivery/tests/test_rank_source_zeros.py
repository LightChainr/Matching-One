from pathlib import Path
from fractions import Fraction as F
import json,sys,unittest
from mpmath import mp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import topological_source_zeros as z
DATA=ROOT/'results/research-control-20260912'

class RankSourceZerosTests(unittest.TestCase):
    def test_double_zero_without_operator(self):
        d=z.source_jet_identity_controls()[0]
        self.assertEqual(d['discriminant'],'0');self.assertEqual(d['u_sum'],'-2');self.assertEqual(d['u_product'],'1')
    def test_complex_and_real_classification(self):
        with mp.workdps(90):
            a=z.source_zeros([mp.mpf(1)/3]*3)
            b=z.source_zeros([mp.mpf('.01'),mp.mpf('.98'),mp.mpf('.01')])
            self.assertEqual(a['type'],'complex_conjugate_pair');self.assertEqual(b['type'],'negative_real_pair')
            for row in a['zeros']+b['zeros']:
                self.assertGreaterEqual(abs(mp.mpf(row['s_imag'])),mp.pi/2)
    def test_cosine_zero_free_identity(self):
        with mp.workdps(90):
            a,b,c=mp.mpf('.2'),mp.mpf('.3'),mp.mpf('.5')
            for x in (-5,0,5):
                for y in (-mp.pi/3,mp.mpf(0),mp.pi/3):
                    actual=mp.re(a*mp.exp(-x-1j*y)+b+c*mp.exp(x+1j*y))
                    exact=b+mp.cos(y)*(a*mp.exp(-x)+c*mp.exp(x))
                    self.assertLess(abs(actual-exact),mp.mpf('1e-80'));self.assertGreater(actual,0)
    def test_newton_trace_vs_matrix_power(self):
        with mp.workdps(90):
            c=[mp.mpf(1),mp.mpf(-3),mp.mpf(1)]
            A=mp.matrix([[3,-1],[1,0]]);seq=z.traces(c,10)
            for m in range(1,11):
                v=A**m;self.assertEqual(seq[m],v[0,0]+v[1,1])
    def test_full_spectrum_vs_independent_enumeration(self):
        h=json.loads((DATA/'full-site-hessian.json').read_text())
        d=z.read_definition(DATA/'width4-parametric-definition.json')
        with mp.workdps(100):
            model=z.Spectra(d)
            for p in (mp.mpf('.5'),mp.mpf('.7')):
                probs=model.probabilities(p,4)[4]
                for i in range(3):
                    direct=sum(a*p**k*(1-p)**(16-k) for k,a in enumerate(h['rank_bernstein_counts'][i]))
                    self.assertLess(abs(direct-probs[i]),mp.mpf('1e-80'))
    def test_extensive_source_three_regions(self):
        d=json.loads((DATA/'topological-source-zeros.json').read_text())
        row=next(x for x in d['parameters'] if x['parameter']=='q4')
        vals=[x for x in row['extensive_source_controls'] if x['length']==256]
        self.assertEqual(len(vals),3)
        for i,x in enumerate(vals):
            self.assertGreater(float(x['tilted_rank_probabilities'][i]),1-1e-9)
    def test_fixed_width_limit_not_uniform_in_aspect(self):
        # Exact analytic control; no percolation exponent is being inferred.
        with mp.workdps(60):
            gamma,a,omega=mp.mpf(2),mp.mpf(1),2
            for w in (4,8,16):
                m=w**(1+omega)
                log_l=-gamma/w+a/w**(1+omega)
                ratio=mp.exp(m*log_l+gamma*m/w)
                self.assertLess(abs(ratio-mp.e),mp.mpf('1e-50'))
    def test_source_probability_validation(self):
        with self.assertRaises(ValueError):z.source_zeros([0,1,0])

if __name__=='__main__':unittest.main()
