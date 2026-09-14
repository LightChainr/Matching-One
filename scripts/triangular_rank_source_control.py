from __future__ import annotations
import argparse,json
from pathlib import Path
import mpmath as mp
TRI=((1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1))

def rank(mask,L):
    N=L*L
    def idx(x,y): return (y%L)*L+(x%L)
    seen=set(); basis=[]
    def add(v):
        nonlocal basis
        if v==(0,0): return
        if not basis: basis=[v]; return
        a=basis[0]
        if a[0]*v[1]-a[1]*v[0]!=0: basis=[a,v]
    for root in range(N):
        if not ((mask>>root)&1) or root in seen: continue
        lift={root:(0,0)}; st=[root]; seen.add(root)
        while st:
            u=st.pop(); x=u%L; y=u//L; ux,uy=lift[u]
            for dx,dy in TRI:
                v=idx(x+dx,y+dy)
                if not ((mask>>v)&1): continue
                prop=(ux+dx,uy+dy)
                if v not in lift:
                    lift[v]=prop; seen.add(v); st.append(v)
                else:
                    dd=(prop[0]-lift[v][0],prop[1]-lift[v][1])
                    if dd!=(0,0):
                        assert dd[0]%L==0 and dd[1]%L==0
                        add((dd[0]//L,dd[1]//L))
                        if len(basis)==2:return 2
    return len(basis)

def cumulants_from_counts(row):
    mass=mp.mpf(sum(row)); N=len(row)-1
    probs=[mp.mpf(v)/mass for v in row]
    def cg(t): return mp.log(mp.fsum(probs[k]*mp.e**(t*k) for k in range(N+1)))
    return [mp.diff(cg,0,n) for n in range(1,7)]

def run(L):
    N=L*L; full=(1<<N)-1; coeff=[[0]*(N+1) for _ in range(3)]; failures=[]
    for m in range(1<<N):
        r=rank(m,L); rc=rank(full^m,L)
        if r+rc!=2 and len(failures)<10: failures.append((m,r,rc))
        coeff[r][m.bit_count()]+=1
    assert not failures
    assert max(abs(coeff[2][k]-coeff[0][N-k]) for k in range(N+1))==0
    P=[mp.mpf(sum(row))/(2**N) for row in coeff]
    a=P[0]; c=P[1]/(2*a); theta=mp.acos(-c)
    kc=[cumulants_from_counts(row) for row in coeff]
    g=[kc[2][i]-kc[0][i] for i in range(6)]
    h=[kc[1][i]-(kc[0][i]+kc[2][i])/2 for i in range(6)]
    rho2=(a*g[0])**2/(2*a*N*mp.mpf('0.25'))
    curve_slope=-2*h[0]/g[0]
    curve_curv=4*(g[0]*h[1]-g[1]*h[0])/g[0]**3
    return {
      'L':L,'configurations':1<<N,'duality_failures':0,
      'rank_count_coefficients':coeff,'p_half_rank_counts':[sum(row) for row in coeff],
      'P0_P1_P2_at_half':[mp.nstr(x,50) for x in P],
      'c_rank_shape':mp.nstr(c,50),'theta_over_pi':mp.nstr(theta/mp.pi,40),
      'conditional_K_cumulants_1_to_6':[[mp.nstr(x,40) for x in row] for row in kc],
      'g_delta_rank2_minus_rank0':[mp.nstr(x,40) for x in g],
      'h_rank1_minus_extreme_average':[mp.nstr(x,40) for x in h],
      'diagnostics':{
        'g1_over_L_3over4':mp.nstr(g[0]/mp.mpf(L)**mp.mpf('.75'),40),
        'g3_over_g1_cubed':mp.nstr(g[2]/g[0]**3,40),
        'sqrtL_times_rho2_KX':mp.nstr(mp.sqrt(L)*rho2,40),
        'd_logc_db':mp.nstr(curve_slope,40),
        'd2_logc_db2':mp.nstr(curve_curv,40)
      }
    }

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();mp.mp.dps=80
    out={'scope':'exact standard six-neighbour triangular-site torus L=3,4; self-matching positive control, not a scaling fit',
         'systems':[run(3),run(4)]}
    a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
