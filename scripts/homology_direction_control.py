from __future__ import annotations
import math,json,argparse
from pathlib import Path
import mpmath as mp
SQ=((1,0),(-1,0),(0,1),(0,-1))
G0=mp.mpf(2)/3; E0=mp.mpf(2)/3

def canon(a,b):
    g=math.gcd(abs(a),abs(b)); a//=g;b//=g
    if a<0 or (a==0 and b<0): a=-a;b=-b
    return a,b

def rank_dir(mask,L):
    N=L*L
    def idx(x,y):return (y%L)*L+(x%L)
    seen=set(); basis=[]
    for root in range(N):
        if not(mask>>root&1) or root in seen:continue
        lift={root:(0,0)}; st=[root];seen.add(root)
        while st:
            u=st.pop();x=u%L;y=u//L;ux,uy=lift[u]
            for dx,dy in SQ:
                v=idx(x+dx,y+dy)
                if not(mask>>v&1):continue
                prop=(ux+dx,uy+dy)
                if v not in lift:
                    lift[v]=prop;seen.add(v);st.append(v)
                else:
                    dd=(prop[0]-lift[v][0],prop[1]-lift[v][1])
                    if dd!=(0,0):
                        assert dd[0]%L==0 and dd[1]%L==0
                        w=(dd[0]//L,dd[1]//L)
                        if not basis:basis=[w]
                        elif basis[0][0]*w[1]-basis[0][1]*w[0]!=0:return 2,None
    if not basis:return 0,None
    return 1,canon(*basis[0])

def finite(L):
    N=L*L; rows={}; rankrows=[[0]*(N+1) for _ in range(3)]
    for mask in range(1<<N):
        r,d=rank_dir(mask,L);k=mask.bit_count();rankrows[r][k]+=1
        if r==1: rows.setdefault(d,[0]*(N+1))[k]+=1
    def mass(row,p):return mp.fsum(mp.mpf(c)*p**k*(1-p)**(N-k) for k,c in enumerate(row))
    f=lambda p:mass(rankrows[2],p)-mass(rankrows[0],p)
    p=mp.findroot(f,(mp.mpf('.55'),mp.mpf('.63')))
    w={d:mass(row,p) for d,row in rows.items()}; p1=mp.fsum(w.values())
    probs={d:v/p1 for d,v in w.items()}
    w2=mp.fsum(v*(d[0]*d[0]+d[1]*d[1]) for d,v in probs.items())
    h10=mp.fsum(v*mp.cos(mp.pi*d[0]) for d,v in probs.items())
    h11=mp.fsum(v*mp.cos(mp.pi*(d[0]+d[1])) for d,v in probs.items())
    return {'L':L,'root':mp.nstr(p,50),'P1':mp.nstr(p1,45),
            'conditional_direction_probs':{f'{a},{b}':mp.nstr(v,40) for (a,b),v in sorted(probs.items())},
            'W2':mp.nstr(w2,40),'H_pi0_over_P1':mp.nstr(h10,40),'H_pipi_over_P1':mp.nstr(h11,40),
            'rank1_direction_count':len(rows)}

def dir_weight(a,b,tau_r=0,tau_i=1,kmax=80):
    s=mp.mpf('0')
    for k in range(-kmax,kmax+1):
        if not k:continue
        m=b*k;n=a*k;g=abs(k)
        ww=mp.e**(-mp.pi*G0*(m*m*tau_i*tau_i+(n-m*tau_r)**2)/tau_i)
        s+=ww*(mp.cos(mp.pi*E0*g)-((-1)**g))
    return s

def continuum_square(Amax=10):
    ds=set()
    for a in range(-Amax,Amax+1):
      for b in range(-Amax,Amax+1):
        if (a,b)==(0,0) or math.gcd(abs(a),abs(b))!=1:continue
        ds.add(canon(a,b))
    w={d:dir_weight(*d) for d in ds}; w={d:v for d,v in w.items() if v>mp.mpf('1e-35')}
    tot=mp.fsum(w.values()); probs={d:v/tot for d,v in w.items()}
    w2=mp.fsum(v*(d[0]*d[0]+d[1]*d[1]) for d,v in probs.items())
    h10=mp.fsum(v*mp.cos(mp.pi*d[0]) for d,v in probs.items())
    h11=mp.fsum(v*mp.cos(mp.pi*(d[0]+d[1])) for d,v in probs.items())
    top=sorted(probs.items(),key=lambda kv:kv[1],reverse=True)[:12]
    return {'tau':'i','conditional_direction_probs_top':{f'{a},{b}':mp.nstr(v,40) for (a,b),v in top},
            'W2':mp.nstr(w2,45),'H_pi0_over_P1':mp.nstr(h10,45),'H_pipi_over_P1':mp.nstr(h11,45),
            'primitive_cut':Amax,'qualification':'conditional rank-one direction law from Arguin Q=1 Gaussian homology sum; common partition prefactor cancels'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();mp.mp.dps=70
    out={'scope':'exact finite square NN L3,L4 balance-root direction law + continuum Q=1 square-torus Arguin/Pinson direction control; no rate fit',
         'finite':[finite(3),finite(4)],'continuum_square':continuum_square()}
    a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
