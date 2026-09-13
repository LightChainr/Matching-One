#!/usr/bin/env python3
"""Exact rational controls for the original-U influence and umbrella identities.
Synthetic finite probability laws only; not physical lattice enumeration.
"""
from fractions import Fraction as F
from itertools import product
import json
import unittest


def add(a,b): return (a[0]+b[0],a[1]+b[1])
def mul(a,b): return (a[0]*b[0],a[0]*b[1]+a[1]*b[0])
def div(a,b): return (a[0]/b[0],(a[1]*b[0]-a[0]*b[1])/b[0]**2)
def sm(xs):
    z=(F(0),F(0))
    for x in xs:z=add(z,x)
    return z

def raw_jet(rows,score,dh):
    ws=[(p,p*(a+dh*k)) for (k,q,p),a in zip(rows,score)]
    z=sm(ws)
    def avg(fn):return div(sm([mul(w,(F(fn(k,q)),F(0))) for (k,q,p),w in zip(rows,ws)]),z)
    K,Q,E=avg(lambda k,q:k),avg(lambda k,q:q),avg(lambda k,q:q*q)
    def covariance(fn,mean):
        kp=avg(lambda k,q:k*fn(k,q)); prod=mul(K,mean)
        return add(kp,(-prod[0],-prod[1]))
    return K,Q,E,covariance(lambda k,q:q,Q),covariance(lambda k,q:q*q,E)

F_ROWS=[(0,-1,F(3,20)),(1,-1,F(3,20)),(1,0,F(1,10)),(3,0,F(1,10)),(3,1,F(1,4)),(4,1,F(1,4))]
S_ROWS=[(0,-1,F(1,4)),(2,-1,F(1,4)),(1,0,F(1,10)),(2,0,F(1,10)),(3,1,F(3,20)),(4,1,F(3,20))]
ROWS=[F_ROWS,S_ROWS]


def full_formula():
    moment=[raw_jet(rows,[F(0)]*len(rows),F(0)) for rows in ROWS]
    D=(moment[0][3][0]+moment[1][3][0])/2
    B=moment[0][4][0]-moment[1][4][0]
    R=B/D
    thermal=[raw_jet(rows,[F(0)]*len(rows),F(1)) for rows in ROWS]
    Dh=(thermal[0][3][1]+thermal[1][3][1])/2
    Bh=thermal[0][4][1]-thermal[1][4][1]
    Rh=(Bh-R*Dh)/D
    phis=[]
    for ell,rows,ms in zip((1,-1),ROWS,moment):
        K,Q,E,Cq,Ce=[x[0] for x in ms]
        phis.append([(ell*((q*q-E)*(k-K)-Ce)-R*((q-Q)*(k-K)-Cq)/2-Rh*(q-Q)/2)/D for k,q,p in rows])
    return phis,D,R,Rh

class Controls(unittest.TestCase):
    def test_centering_and_common_thermal_null(self):
        phi,D,R,Rh=full_formula()
        self.assertEqual(sum(sum(p*x for (k,q,p),x in zip(rows,ph)) for rows,ph in zip(ROWS,phi)),0)
        for rows,ph in zip(ROWS,phi):self.assertEqual(sum(p*x for (k,q,p),x in zip(rows,ph)),0)
        self.assertEqual(sum(sum(p*k*x for (k,q,p),x in zip(rows,ph)) for rows,ph in zip(ROWS,phi)),0)

    def test_12_independent_weight_directions_via_dual_arithmetic(self):
        phi,D,R,Rh=full_formula()
        for which in range(2):
            for i in range(6):
                scores=[[F(0)]*6 for _ in range(2)];scores[which][i]=F(1)
                fixed=[raw_jet(rows,a,F(0)) for rows,a in zip(ROWS,scores)]
                dh=-(fixed[0][1][1]+fixed[1][1][1])/(2*D)
                moved=[raw_jet(rows,a,dh) for rows,a in zip(ROWS,scores)]
                self.assertEqual((moved[0][1][1]+moved[1][1][1])/2,0)
                B=add(moved[0][4],tuple(-x for x in moved[1][4]))
                den=mul(add(moved[0][3],moved[1][3]),(F(1,2),F(0)))
                direct=div(B,den)[1]
                predicted=ROWS[which][i][2]*phi[which][i]
                self.assertEqual(direct,predicted)

    def test_umbrella_normalization_and_exact_second_moment(self):
        phis,*_=full_formula();bias=[F(2),F(7),F(3)]
        for rows,phi in zip(ROWS,phis):
            normalizer=sum(p*bias[q+1] for k,q,p in rows)
            nu=[p*bias[q+1]/normalizer for k,q,p in rows]
            weights=[1/bias[q+1] for k,q,p in rows]
            z=sum(p*w for p,w in zip(nu,weights))
            for fn in (lambda k,q:k,lambda k,q:q,lambda k,q:q*q,lambda k,q:k*q,lambda k,q:k*k):
                actual=sum(v*w*fn(k,q) for (k,q,p),v,w in zip(rows,nu,weights))/z
                self.assertEqual(actual,sum(p*fn(k,q) for k,q,p in rows))
            actual=sum(p*(w/z*x)**2 for p,w,x in zip(nu,weights,phi))
            predicted=normalizer*sum(p*x*x/bias[q+1] for (k,q,p),x in zip(rows,phi))
            self.assertEqual(actual,predicted)

    def test_kantorovich_robustness_all_27_cases(self):
        p=[F(1,5),F(1,3),F(7,15)];rootM=[F(2),F(3),F(5)]
        optimum=sum(a*b for a,b in zip(p,rootM))**2
        for ratios in product((F(1,2),F(1),F(2)),repeat=3):
            alpha=[m*r for m,r in zip(rootM,ratios)]
            variance=sum(a*b for a,b in zip(p,alpha))*sum(a*m*m/b for a,m,b in zip(p,rootM,alpha))
            self.assertGreaterEqual(variance,optimum)
            self.assertLessEqual(variance,F(25,16)*optimum)

if __name__=='__main__':unittest.main(verbosity=2)
