from __future__ import annotations
import math,json
import mpmath as mp
mp.mp.dps=40
CUT=18

def pars(Q):
    Q=mp.mpf(Q);e0=2/mp.pi*mp.acos(mp.sqrt(Q)/2);return e0,4-2*e0

def probs(Q,tr=0,ti=1,cut=CUT):
    Q=mp.mpf(Q);e0,g=pars(Q);g0=g/4
    A=mp.mpf('0');R=mp.mpf('0')
    for m in range(-cut,cut+1):
      for n in range(-cut,cut+1):
        w=mp.e**(-mp.pi*g0*(m*m*ti*ti+(n-m*tr)**2)/ti)
        h=0 if (m==0 and n==0) else math.gcd(abs(m),abs(n))
        A += w*((-1)**h)
        if m or n:R += w*(mp.cos(mp.pi*e0*h)-((-1)**h))
    z0=A/2;z2=Q*z0;z=z0+R+z2
    return z0/z,R/z,z2/z

def diag(Q,tr=0,ti=1):
    p0,p1,p2=probs(Q,tr,ti)
    b=mp.log(p0/p2)/2;d=mp.log(p1/mp.sqrt(p0*p2));c=p1/(2*mp.sqrt(p0*p2))
    s=b;w0=p0*mp.e**(-s);w1=p1;w2=p2*mp.e**s;zz=w0+w1+w2
    return p0,p1,p2,b,d,c,(w0/zz,w1/zz,w2/zz)

def cqr(Q,r):return diag(Q,0,mp.mpf(r))[5]
def bisect_root(Q,lo=1,hi=3,it=45):
    lo=mp.mpf(lo);hi=mp.mpf(hi);flo=cqr(Q,lo)-1;fhi=cqr(Q,hi)-1
    assert flo<0<fhi
    for _ in range(it):
        mid=(lo+hi)/2;fm=cqr(Q,mid)-1
        if fm>0:hi=mid
        else:lo=mid
    return (lo+hi)/2

def ddlogQ(Q,h=mp.mpf('1e-5')):
    Q=mp.mpf(Q);qp=Q*mp.e**h;qm=Q*mp.e**(-h)
    return (diag(qp)[4]-diag(qm)[4])/(2*h)

rows={}
for Q in (1,2,3,4):
    p0,p1,p2,b,d,c,pbal=diag(Q);rstar=bisect_root(Q)
    der=None if Q==4 else ddlogQ(Q)
    a=pbal[0];r=pbal[1];odd=2*a*mp.mpf('.25');even=None if der is None else 2*a*r*der**2
    rows[str(Q)]={'P0':mp.nstr(p0,30),'P1':mp.nstr(p1,30),'P2':mp.nstr(p2,30),
      'b':mp.nstr(b,30),'minus_half_logQ':mp.nstr(-mp.log(Q)/2,30),'d':mp.nstr(d,30),'c':mp.nstr(c,30),
      'source_balanced_rank_law':[mp.nstr(x,30) for x in pbal],
      'd_dlogQ_even_shape':('endpoint_not_differentiated' if der is None else mp.nstr(der,30)),
      'fisher_odd_norm2':mp.nstr(odd,30),
      'fisher_even_norm2':('endpoint_not_differentiated' if even is None else mp.nstr(even,30)),
      'fisher_even_fraction':('endpoint_not_differentiated' if even is None else mp.nstr(even/(odd+even),30)),
      'collision_r':mp.nstr(rstar,30)}
print(json.dumps({'scope':'Arguin critical FK homology Q=1..4; cut=18 Gaussian sums, 40 dps; Q derivative via symmetric log-Q step 1e-5 for Q<4; collision root via bisection','rows':rows},indent=2))
