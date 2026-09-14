from __future__ import annotations
from collections import defaultdict
import mpmath as mp, json, math
from pathlib import Path
import argparse

DIRS=((1,0),(-1,0),(0,1),(0,-1))

def torus_rank(mask,L):
    N=L*L
    def idx(x,y): return (y%L)*L+(x%L)
    seen=set(); basis=[]
    def add_vec(v):
        nonlocal basis
        if v==(0,0): return
        if not basis: basis=[v]; return
        a=basis[0]
        if a[0]*v[1]-a[1]*v[0] !=0: basis=[a,v]
    for root in range(N):
        if not ((mask>>root)&1) or root in seen: continue
        lift={root:(0,0)}; stack=[root]; seen.add(root)
        while stack:
            u=stack.pop(); x=u%L; y=u//L; ux,uy=lift[u]
            for dx,dy in DIRS:
                v=idx(x+dx,y+dy)
                if not ((mask>>v)&1): continue
                prop=(ux+dx,uy+dy)
                if v not in lift:
                    lift[v]=prop; seen.add(v); stack.append(v)
                else:
                    dd=(prop[0]-lift[v][0],prop[1]-lift[v][1])
                    if dd!=(0,0):
                        assert dd[0]%L==0 and dd[1]%L==0
                        add_vec((dd[0]//L,dd[1]//L))
                        if len(basis)==2: return 2
    return len(basis)

def counts_by_rank_k(L):
    N=L*L
    c=[[0]*(N+1) for _ in range(3)]
    for m in range(1<<N): c[torus_rank(m,L)][m.bit_count()]+=1
    return c

def eval_mass(coeff,p):
    N=len(coeff)-1; q=1-p
    return mp.fsum(mp.mpf(c)*p**k*q**(N-k) for k,c in enumerate(coeff))

def solve(L):
    c=counts_by_rank_k(L); N=L*L
    f=lambda p: eval_mass(c[2],p)-eval_mass(c[0],p)
    root=mp.findroot(f,(mp.mpf('.55'),mp.mpf('.63')))
    P=[eval_mass(c[r],root) for r in range(3)]
    a=(P[0]+P[2])/2
    assert abs(P[0]-P[2])<mp.mpf('1e-45')
    theta=mp.acos(1-1/(2*a)) if a>=mp.mpf('0.25') else mp.pi
    def cg(s): return mp.log(P[0]*mp.e**(-s)+P[1]+P[2]*mp.e**s)
    cumul=[mp.diff(cg,0,n) for n in range(1,7)]
    kc=[]
    for r in (0,2):
        weights=[mp.mpf(c[r][k])*root**k*(1-root)**(N-k) for k in range(N+1)]
        mass=mp.fsum(weights)
        mu=mp.fsum(mp.mpf(k)*weights[k] for k in range(N+1))/mass
        var=mp.fsum((mp.mpf(k)-mu)**2*weights[k] for k in range(N+1))/mass
        k3=mp.fsum((mp.mpf(k)-mu)**3*weights[k] for k in range(N+1))/mass
        kc.append((mu,var,k3))
    g1=kc[1][0]-kc[0][0]; g2=kc[1][1]-kc[0][1]; g3=kc[1][2]-kc[0][2]
    z1=-2/g1
    z2=-4*g2/g1**3
    z3=8*(g1*g3-3*g2**2)/g1**5
    EX2=2*a
    EXK=P[2]*kc[1][0]-P[0]*kc[0][0]
    proj_delta=2*EXK/EX2
    assert abs(proj_delta-g1)<mp.mpf('1e-45')
    return {
      'L':L,'root':mp.nstr(root,50),
      'P0':mp.nstr(P[0],50),'P1':mp.nstr(P[1],50),'P2':mp.nstr(P[2],50),
      'a':mp.nstr(a,50),'nearest_zero_imag_if_pure_imag':mp.nstr(theta,50),
      'theta_over_pi':mp.nstr(theta/mp.pi,40),
      'cumulants_1_to_6':[mp.nstr(x,40) for x in cumul],
      'conditional_K_cumulants':{
        'rank0':{'mean':mp.nstr(kc[0][0],40),'variance':mp.nstr(kc[0][1],40),'third':mp.nstr(kc[0][2],40)},
        'rank2':{'mean':mp.nstr(kc[1][0],40),'variance':mp.nstr(kc[1][1],40),'third':mp.nstr(kc[1][2],40)},
        'delta_g1_g2_g3':[mp.nstr(g1,40),mp.nstr(g2,40),mp.nstr(g3,40)]},
      'topological_source_root_logit_derivatives_s1_s2_s3':[mp.nstr(z1,40),mp.nstr(z2,40),mp.nstr(z3,40)],
      'projector_deltaK_abs_error':mp.nstr(abs(proj_delta-g1),5),
      'rank_count_coefficients':c,
    }

def theta3(q):
    return mp.jtheta(3,0,q)

def eta_q(q):
    prod=mp.mpf(1); n=1
    while True:
        t=q**n; prod*=1-t
        if t < mp.mpf('1e-90'): break
        n+=1
    return q**(mp.mpf(1)/24)*prod

def continuum():
    num=(theta3(mp.e**(-3*mp.pi/8))*theta3(mp.e**(-8*mp.pi/3))
         -theta3(mp.e**(-3*mp.pi/2))*theta3(mp.e**(-2*mp.pi/3)))
    a=num/(2*eta_q(mp.e**(-2*mp.pi))**2)
    R_e=1-a
    theta=mp.acos(1-1/(2*a))
    return {
      'formula_source':'Newman-Ziff Eq. (11) from Pinson; a_cross=a_no_wrap=1-R_e',
      'R_e_any_wrap':mp.nstr(R_e,60),
      'a_cross_equals_no_wrap':mp.nstr(a,60),
      'theta':mp.nstr(theta,60),'theta_over_pi':mp.nstr(theta/mp.pi,60),
      'variance_X':mp.nstr(2*a,60),
      'raw_kurtosis':mp.nstr(1/(2*a),60),
      'excess_kurtosis':mp.nstr(1/(2*a)-3,60),
      'kappa4':mp.nstr(2*a*(1-6*a),60),
      'kappa6':mp.nstr(2*a*(1-30*a+120*a*a),60),
    }

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path('topological-rank-source-control.json'))
    args=parser.parse_args()
    mp.mp.dps=70
    out={'scope':'exact L=3,4 rank enumeration + analytic Pinson/Newman-Ziff square-torus control; no asymptotic fit',
         'finite':[solve(3),solve(4)],'continuum_from_Pinon_NewmanZiff':continuum()}
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
