"""Checks of conditional algebra and old-data arithmetic, not SITE scaling."""
import sys
import unittest
from pathlib import Path
from fractions import Fraction as F
from math import comb
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import thermal_ward_root_control as w

class WardControls(unittest.TestCase):
    def test_virasoro_commutator(self):
        for word in [(),(1,),(2,1),(3,)]:
            v={word:F(1)}
            for a,b in [(1,-2),(2,-1),(2,-2),(3,-2)]:
                lhs=w.add(w.act(a,w.act(b,v)),w.scale(-F(1),w.act(b,w.act(a,v))))
                rhs=w.scale(F(a-b),w.act(a+b,v))
                self.assertEqual(lhs,rhs)

    def test_null_positive_modes(self):
        chi,_,_=w.descendants()
        for n in range(1,6): self.assertEqual(w.act(n,chi),{})

    def test_quasiprimary(self):
        _,q,U=w.descendants()
        self.assertEqual(w.act(1,U),{})
        self.assertEqual(w.act(0,U),w.scale(w.H+4,U))

    def test_normalized_quotient_class(self):
        _,_,U=w.descendants(); basis=list(w.partitions(4)); sub=w.null_derivative_subspace()
        self.assertEqual(w.rational_rank(sub,basis),4)
        self.assertEqual(w.rational_rank(sub+[w.add(U,{(4,):-F(1)})],basis),4)
        self.assertEqual(w.rational_rank(sub+[{(4,):F(1)}],basis),5)

    def test_retained_null_anomaly(self):
        chi,_,U=w.descendants(); basis=list(w.partitions(4))
        ds=[w.act(-1,{v:F(1)}) for v in w.partitions(3)]
        rem=w.add(U,{(4,):-F(1)},w.scale(-F(80,87),w.act(-2,chi)))
        self.assertEqual(w.rational_rank(ds,basis),w.rational_rank(ds+[rem],basis))

    def test_magnetic_fusion_and_zero_mode(self):
        m=F(2,3)*w.H*(w.H+1)-w.H
        self.assertEqual(m,F(5,96))
        self.assertEqual(m-w.H/12,0)
        self.assertEqual(w.H/240,F(1,384))
        self.assertEqual(3*w.H/45,F(1,24))

    def test_eta_two_routes(self):
        self.assertEqual(w.eta_product_coeffs(2*w.H,14),w.eta_ward_coeffs(w.H,14))
        self.assertEqual(w.eta_product_coeffs(2*w.H,3),[F(1),-F(5,4),-F(35,32),F(45,128)])

    def test_modular_covariance(self):
        import mpmath as mp
        with mp.workdps(65):
            tau=mp.mpc('.21','1.17'); E=w.eisenstein(tau)
            self.assertLess(abs(w.eisenstein(tau+1)-E),mp.mpf('1e-55'))
            self.assertLess(abs(w.eisenstein(-1/tau)-tau**4*E),mp.mpf('1e-55'))
            # Physical basis invariance: u'=u*tau, tau'=-1/tau.
            u=mp.mpc('2.3','.4')
            self.assertLess(abs((u*tau)**-4*w.eisenstein(-1/tau)-u**-4*E),mp.mpf('1e-55'))

    def test_hex_simple_zero(self):
        import mpmath as mp
        with mp.workdps(65):
            tau=mp.mpc(mp.mpf('.5'),mp.sqrt(3)/2)
            self.assertLess(abs(w.eisenstein(tau)),mp.mpf('1e-55'))
            expected=-(2j*mp.pi/3)*w.eisenstein(tau,6)
            self.assertGreater(abs(expected),1)
            self.assertLess(abs(mp.diff(lambda z:w.eisenstein(z),tau)-expected),mp.mpf('1e-50'))

    def test_supplied_counts(self):
        for L,d in w.RANK_INPUT.items():
            for k in range(L*L+1):
                self.assertEqual(sum(d[f'rank{j}'][k] for j in range(3)),comb(L*L,k))
            self.assertEqual(sum(map(sum,[d[f'rank{j}'] for j in range(3)])),2**(L*L))

    def test_root_residuals(self):
        import mpmath as mp
        with mp.workdps(60):
            for record in w.root_diagnostics():
                data=w.RANK_INPUT[record['L']]; p=mp.mpf(record['torus_root'])
                self.assertLess(abs(w.rank_probability(data['rank2'],p)-w.rank_probability(data['rank0'],p)),mp.mpf('1e-33'))

    def test_report_precision(self):
        a=w.report(60); b=w.report(90)
        self.assertEqual(a['exact'],b['exact'])
        self.assertEqual(a['existing_data_diagnostics'],b['existing_data_diagnostics'])
        for field in ['E4_i','E4_2i']:
            self.assertEqual(a['modular_numerical'][field],b['modular_numerical'][field])

if __name__=='__main__': unittest.main()
