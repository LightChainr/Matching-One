from __future__ import annotations
import json, math
import numpy as np

SQ=((1,0),(-1,0),(0,1),(0,-1))
TRI=((1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1))
P_SQUARE={3:0.58651145511267563565455897660690173482430062489384,
          4:0.59067211233102829689590201143951286962111713272216}
G1={'square':{3:2.945015004456800709120992310384935214928,
              4:3.730045260798946011338924260846967120335},
    'triangular':{3:3.14457831325301204819277108433734939759,
                  4:3.918351810439168290431555407527239853451}}

def torus_rank(mask,L,dirs):
    N=L*L;seen=set();basis=[]
    def idx(x,y):return (y%L)*L+(x%L)
    def add(v):
        nonlocal basis
        if v==(0,0):return
        if not basis:basis=[v];return
        a=basis[0]
        if a[0]*v[1]-a[1]*v[0]!=0:basis=[a,v]
    for root in range(N):
        if not(mask>>root&1) or root in seen:continue
        lift={root:(0,0)};stack=[root];seen.add(root)
        while stack:
            u=stack.pop();x=u%L;y=u//L;ux,uy=lift[u]
            for dx,dy in dirs:
                v=idx(x+dx,y+dy)
                if not(mask>>v&1):continue
                prop=(ux+dx,uy+dy)
                if v not in lift:
                    lift[v]=prop;seen.add(v);stack.append(v)
                else:
                    dd=(prop[0]-lift[v][0],prop[1]-lift[v][1])
                    if dd!=(0,0):
                        assert dd[0]%L==0 and dd[1]%L==0
                        add((dd[0]//L,dd[1]//L))
                        if len(basis)==2:return 2
    return len(basis)

def transform(f,p,N):
    a=f.astype(float).copy();q=1-p;sq=math.sqrt(p*q)
    for i in range(N):
        step=1<<i;block=step<<1
        for start in range(0,len(a),block):
            f0=a[start:start+step].copy();f1=a[start+step:start+block].copy()
            a[start:start+step]=q*f0+p*f1
            a[start+step:start+block]=sq*(f1-f0)
    return a

def run(name,L,dirs,p):
    N=L*L
    X=np.array([torus_rank(m,L,dirs)-1 for m in range(1<<N)],dtype=float)
    c=transform(X,p,N)
    W=np.zeros(N+1)
    for mask,x in enumerate(c):W[mask.bit_count()]+=x*x
    var=float(W[1:].sum());mean=float(sum(k*W[k] for k in range(1,N+1))/var)
    m2=float(sum(k*k*W[k] for k in range(1,N+1))/var)
    m3=float(sum(k**3*W[k] for k in range(1,N+1))/var)
    even=float(sum(W[k] for k in range(2,N+1,2))/var)
    omega=2*mean/G1[name][L]-1
    lap={str(t):float(sum(W[k]/var*math.exp(-t*k/mean) for k in range(1,N+1))) for t in (.5,1,2,4)}
    return {'model':name,'L':L,'p':p,'configurations':1<<N,'mean_X':float(c[0]),'var_X':var,
            'level_weights':[float(x) for x in W],
            'spectral_size_mean':mean,
            'mean_over_L_3over4':mean/(L**.75),
            'spectral_size_cv2':m2/mean**2-1,
            'spectral_size_third_raw_normalized':m3/mean**3,
            'even_level_fraction':even,
            'jump2_share_from_spectrum':omega,
            'mean_normalized_laplace':lap}

out={'scope':'exact configuration rank values; floating p-biased tensor Fourier transform; square at exact finite balance roots, triangular at self-matching p=1/2',
     'systems':[]}
for L in (3,4):out['systems'].append(run('square',L,SQ,P_SQUARE[L]))
for L in (3,4):out['systems'].append(run('triangular',L,TRI,.5))
print(json.dumps(out,indent=2))
