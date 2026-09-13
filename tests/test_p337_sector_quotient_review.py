"""Focused independent checks for the P337 review, no Monte Carlo."""
import csv
import importlib.util
import json
import math
import sys
import tempfile
import unittest
from decimal import Decimal, localcontext
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "p337_sector_quotient_review", ROOT / "scripts/p337_sector_quotient_review.py")
M = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = M
SPEC.loader.exec_module(M)


def midpoint(entry):
    return (F(entry["lower"]) + F(entry["upper"])) / 2


def dec(frac):
    return Decimal(frac.numerator) / Decimal(frac.denominator)


def direct_row_moments(rows, h, m, score_shift=0):
    """Independent decimal row sums, NOT the integer polynomial evaluator."""
    weights = [Decimal(n) * h**k / Decimal(m)**g for k,g,q,n in rows]
    total = sum(weights)
    ps = [w/total for w in weights]
    K = [Decimal(k) for k,g,q,n in rows]
    q = [Decimal(q) for k,g,q,n in rows]
    E = [v*v for v in q]
    S = [Decimal(-g+score_shift*k) for k,g,q,n in rows]
    def mean(v):
        return sum(p*x for p,x in zip(ps,v))
    def covariance(v,w):
        a,b=mean(v),mean(w)
        return sum(p*(x-a)*(y-b) for p,x,y in zip(ps,v,w))
    def third(v,w,u):
        a,b,c=mean(v),mean(w),mean(u)
        return sum(p*(x-a)*(y-b)*(z-c) for p,x,y,z in zip(ps,v,w,u))
    out={"K":mean(K)}
    for name,o in (("q",q),("E",E)):
        out[name]=mean(o);out[name+"x"]=covariance(o,K)
        out[name+"xx"]=third(o,K,K)
        out["J"+name]=covariance(o,S);out["J"+name+"x"]=third(o,S,K)
    for suffix,condition in (("full",lambda g:True),
            ("remainder",lambda g:g>min(r[1] for r in rows if r[2]==0))):
        active=[(p,k) for p,k,row in zip(ps,K,rows) if row[2]==0 and condition(row[1])]
        mass=sum(p for p,k in active)
        C=sum(p*(k-out["K"]) for p,k in active)
        absolute=sum(p*abs(k-out["K"]) for p,k in active)
        second=sum(p*(k-out["K"])**2 for p,k in active)
        km=sum(p*k for p,k in active)/mass
        conditional=mass*sum(p*(k-km)**2 for p,k in active)
        out[suffix] = {"mass":mass,"mean":C,"absolute":absolute,
            "iid_variance":second-C*C,"one_proposal_floor":absolute*absolute-C*C,
            "conditional_oracle_variance":conditional}
    return out


class SectorReviewTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows=[M.load_rows(M.PACKAGE/"inputs"/(n+".csv"),M.HASHES[n])
                  for n in ("axis","tilted")]
        # A stored report is optional: the connector did not publish the
        # generated JSON. Clean checkouts reconstruct it from the pinned CSVs.
        stored = M.PACKAGE / "result.json"
        cls.report = (json.loads(stored.read_text()) if stored.exists()
                      else M.analyze(M.PACKAGE / "inputs"))

    def test_source_multiplicities(self):
        for rows in self.rows:
            self.assertEqual(sum(row[3] for row in rows),2**25)
            self.assertEqual(sum(row[3] for row in rows if row[0]==12),math.comb(25,12))

    def test_corrupted_source_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/"bad.csv"
            p.write_bytes((M.PACKAGE/"inputs/axis.csv").read_bytes()+b"\n")
            with self.assertRaises(ValueError):M.load_rows(p,M.HASHES["axis"])

    def test_exact_horner(self):
        c=[i*i-30 for i in range(26)]
        for h in (F(1,3),F(1),F(7,5)):
            direct=sum(F(a)*h**k for k,a in enumerate(c))
            self.assertEqual(F(M.homogeneous_eval(c,h),h.denominator**25),direct)

    def test_interval_arithmetic(self):
        a,b=M.Interval(F(-2,3),F(7,5)),M.Interval(F(1,3),F(3,2))
        for x in (F(-2,3),F(0),F(7,5)):
            for y in (F(1,3),F(1),F(3,2)):
                self.assertTrue((a+b).contains(x+y))
                self.assertTrue((a*b).contains(x*y))
                self.assertTrue((a/b).contains(x/y))
        with self.assertRaises(ZeroDivisionError):a/M.Interval(-1,1)
        with self.assertRaises(ValueError):M.Interval(2,1)

    def test_signed_polynomial_interval(self):
        c=[(-1)**i*(i+1) for i in range(26)]
        h=M.Interval(F(1,2),F(7,10));bound=M.polynomial_bound(c,h)
        for x in (F(1,2),F(3,5),F(7,10)):
            self.assertTrue(bound.contains(sum(F(a)*x**k for k,a in enumerate(c))))

    def test_rank_tilt_preserves_projective_invariant(self):
        z=[F(2),F(7),F(3)]
        for v in (F(1,4),F(2),F(64)):
            zp=[z[r]*v**r for r in range(3)]
            self.assertEqual(z[1]**2/(z[0]*z[2]),zp[1]**2/(zp[0]*zp[2]))
            self.assertEqual(zp[2]/zp[0],v*v*z[2]/z[0])

    def test_rank_one_insertion_changes_only_middle_coordinate(self):
        z=[F(2),F(7),F(3)];v=F(5,4);zp=[z[0],z[1]*v,z[2]]
        self.assertEqual(zp[2]/zp[0],z[2]/z[0])
        self.assertEqual(zp[1]**2/(zp[0]*zp[2]),v*v*z[1]**2/(z[0]*z[2]))

    def test_conditional_source_two_coordinates(self):
        z=[F(2),F(7),F(3)];p=[a/sum(z) for a in z]
        source=[F(-3),F(2),F(9)]
        a=(source[2]-source[0])/2;b=source[1]-(source[0]+source[2])/2
        q=p[2]-p[0];e=p[0]+p[2];s=sum(p[i]*source[i] for i in range(3))
        jq=sum(p[i]*(i-1)*source[i] for i in range(3))-q*s
        je=sum(p[i]*(i-1)**2*source[i] for i in range(3))-e*s
        self.assertEqual(jq,(e-q*q)*a-p[1]*q*b)
        self.assertEqual(je,p[1]*q*a-p[1]*e*b)

    def test_own_root_thermal_ratio(self):
        p=[F(2,5),F(1,5),F(2,5)];k=[F(2),F(7),F(10)]
        km=sum(x*y for x,y in zip(p,k))
        qx=sum(p[i]*(i-1)*(k[i]-km) for i in range(3))
        ex=sum(p[i]*(i-1)**2*(k[i]-km) for i in range(3))
        a=(k[2]-k[0])/2;b=k[1]-(k0:=k[0]+k[2])/2
        self.assertEqual(ex/qx,-p[1]*b/a)

    def test_one_proposal_floor_and_equality(self):
        p=[F(1,4),F(1,2),F(1,4)];x=[F(-2),F(1),F(5)]
        mean=sum(a*b for a,b in zip(p,x));ab=sum(a*abs(b) for a,b in zip(p,x))
        optimum=[a*abs(b)/ab for a,b in zip(p,x)]
        for proposal in (p,[F(1,3)]*3,optimum):
            second=sum(a*a*b*b/c for a,b,c in zip(p,x,proposal))
            self.assertGreaterEqual(second-mean**2,ab**2-mean**2)
        self.assertEqual(sum(a*a*b*b/c for a,b,c in zip(p,x,optimum)),ab**2)

    def test_root_brackets_have_correct_sign(self):
        for m,key in ((1,"Q1_comoving_split"),(64,"m64_star_oracle_budget")):
            laws=[M.FiniteLaw(rows,m) for rows in self.rows]
            root=self.report[key]["h"]
            for side,sgn in (("lower",-1),("upper",1)):
                x=F(root[side]);z=[M.homogeneous_eval(a.total,x) for a in laws]
                q=[M.homogeneous_eval(a.qnumer,x) for a in laws]
                self.assertGreaterEqual(sgn*(q[0]*z[1]+q[1]*z[0]),0)

    def test_q1_source_response_independent_centered_rows(self):
        with localcontext() as ctx:
            ctx.prec=90;entry=self.report["Q1_comoving_split"]
            h=dec(midpoint(entry["h"]))
            f,s=[direct_row_moments(rows,h,1) for rows in self.rows]
            D=(f["qx"]+s["qx"])/2;Yx=(f["Ex"]-s["Ex"])/dec(M.DELTA)
            Qxx=(f["qxx"]+s["qxx"])/2;Yxx=(f["Exx"]-s["Exx"])/dec(M.DELTA)
            jq=(f["Jq"]+s["Jq"])/2;jqx=(f["Jqx"]+s["Jqx"])/2
            jex=(f["JEx"]-s["JEx"])/dec(M.DELTA);w=jq/D
            v=(jex-Yxx*w)/D-Yx*(jqx-Qxx*w)/(D*D)
            self.assertLess(abs(v-dec(midpoint(entry["V_over_A"]))),Decimal('1e-32'))

    def test_m64_variance_budgets_independent_direct_squares(self):
        with localcontext() as ctx:
            ctx.prec=95;entry=self.report["m64_star_oracle_budget"]
            h=dec(midpoint(entry["h"]))
            f,s=[direct_row_moments(rows,h,64) for rows in self.rows]
            theta=f["full"]["mean"]-s["full"]["mean"]
            for block,k,out in (("full","iid_variance","iid_n_each"),
                ("full","one_proposal_floor","one_proposal_n_each_floor"),
                ("full","conditional_oracle_variance","rank1_conditional_n_each_oracle"),
                ("remainder","one_proposal_floor","after_shell_n_each_floor"),
                ("remainder","conditional_oracle_variance","after_shell_conditional_n_each_oracle")):
                val=9*(f[block][k]+s[block][k])/(theta*theta)
                target=dec(midpoint(entry[out]));self.assertLess(abs(val/target-1),Decimal('1e-22'))

    def test_new_claim_signs_and_scope(self):
        r=self.report;split=r["Q1_comoving_split"];budget=r["m64_star_oracle_budget"]
        self.assertGreater(F(split["bias_over_A"]["lower"]),0)
        self.assertLess(F(split["middle_over_A"]["upper"]),0)
        self.assertGreater(F(budget["one_proposal_n_each_floor"]["lower"]),10**10)
        self.assertLess(F(budget["after_shell_conditional_n_each_oracle"]["upper"]),700000)
        self.assertEqual(r["new_random_samples"],0)
        self.assertEqual(r["official_scorers_run"],0)
        self.assertIn("no practical sampler",r["claim_boundary"])

    def test_gauge_covariance_algebra(self):
        # Exact local jet identity; c0+c1*x is an arbitrary common clock.
        a,ax,at,atx,w,wx=map(F,(2,3,5,7,11,13));c0,c1=F(17),F(19)
        unchanged=(at-w*a,atx-wx*a-w*ax)
        shifted=(at+c0*a-(w+c0)*a,
            atx+c1*a+c0*ax-(wx+c1)*a-(w+c0)*ax)
        self.assertEqual(unchanged,shifted)

    def test_certified_sum_and_source_gauge_reconstruction(self):
        rows=self.rows;laws=[M.FiniteLaw(r,1) for r in rows]
        point=midpoint(self.report["Q1_comoving_split"]["h"])
        h=M.Interval.point(point)
        a=M.comoving_split(laws,h,0);b=M.comoving_split(laws,h,-1)
        for k in ("V_over_A","bias_over_A","middle_over_A"):
            self.assertTrue(a[k].overlaps(b[k]))
        self.assertTrue((a["bias_over_A"]+a["middle_over_A"]).overlaps(a["V_over_A"]))


if __name__ == "__main__":
    unittest.main()
