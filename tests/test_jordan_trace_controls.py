import sys,unittest,tempfile,json
from pathlib import Path
from fractions import Fraction as F
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import sympy as s
import jordan_trace_controls as j
import verify_jordan_trace_controls as v

class TestJordanControls(unittest.TestCase):
    def test_nilpotent_and_commutation(self):
        P,Q,D,N=j.matrices()
        self.assertEqual(N*N,s.zeros(3));self.assertNotEqual(N,s.zeros(3))
        self.assertEqual(D*N,N);self.assertEqual(N*D,-N)
    def test_positive_markov_and_shared_stationary(self):
        for t in [-s.Rational(1,4),0,s.Rational(1,4)]:
            for G in j.markov_pair(t):
                self.assertEqual(G*s.ones(3,1),s.zeros(3,1))
                self.assertEqual(s.ones(1,3)*G,s.zeros(1,3))
                self.assertGreater(min(G[i,k] for i in range(3) for k in range(3) if i!=k),0)
    def test_jordan_nullities(self):
        A,B=j.markov_pair(0)
        self.assertEqual(len((A+s.eye(3)).nullspace()),2)
        self.assertEqual(len((B+s.eye(3)).nullspace()),1)
    def test_probability_readout_not_trace(self):
        A,B=j.markov_pair(0);A=s.eye(3)+A/2;B=s.eye(3)+B/2
        self.assertEqual(s.trace(A*A),s.trace(B*B))
        self.assertNotEqual((A*A)[0,2],(B*B)[0,2])
    def test_affine_nonanalytic_pencil(self):
        t,x=s.symbols('t x');C=s.Matrix([[0,t,0],[0,0,t],[t,0,1]])
        self.assertEqual(s.expand(C.charpoly(x).as_expr()-(x*x*(x-1)-t**3)),0)
        self.assertEqual(len(C.subs(t,0).nullspace()),2)
    def test_independent_fraction_validator(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'report.json';p.write_text(json.dumps(j.run_checks()))
            r=v.verify(p);self.assertTrue(r['success'])
            self.assertEqual(r['hankel_determinants'],['-1/36','-1/6912'])

if __name__=='__main__':unittest.main()
