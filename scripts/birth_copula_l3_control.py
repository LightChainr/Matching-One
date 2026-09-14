from __future__ import annotations
import itertools, math, json
from fractions import Fraction
import numpy as np
from scipy.optimize import linprog

SQ=((1,0),(-1,0),(0,1),(0,-1))
TRI=((1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1))

def torus_rank(mask,L,dirs):
    N=L*L
    def idx(x,y):return (y%L)*L+(x%L)
    seen=set();basis=[]
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

def static_coeff(rank,N):
    c=[[0]*(N+1) for _ in range(3)]
    for m,r in enumerate(rank):c[r][m.bit_count()]+=1
    return c

def marginals_from_static(c):
    N=len(c[0])-1
    q0=[Fraction(c[0][k],math.comb(N,k)) for k in range(N+1)]
    q2=[Fraction(c[2][k],math.comb(N,k)) for k in range(N+1)]
    p1=[q0[j-1]-q0[j] for j in range(1,N+1)]
    p2=[q2[j]-q2[j-1] for j in range(1,N+1)]
    return p1,p2

def lp_bounds(c):
    p1,p2=marginals_from_static(c);N=len(p1)
    V=[(i,j) for i in range(N) for j in range(i,N)];idx={v:k for k,v in enumerate(V)}
    A=[];b=[]
    for i in range(N):
        row=np.zeros(len(V));
        for j in range(i,N):row[idx[(i,j)]]=1
        A.append(row);b.append(float(p1[i]))
    for j in range(N):
        row=np.zeros(len(V));
        for i in range(j+1):row[idx[(i,j)]]=1
        A.append(row);b.append(float(p2[j]))
    A=np.array(A);b=np.array(b)
    diag=np.array([1.0 if i==j else 0.0 for i,j in V])
    d2=np.array([float((j-i)**2) for i,j in V])
    out={}
    for name,cost in [('p_equal',diag),('D2',d2)]:
        mn=linprog(cost,A_eq=A,b_eq=b,bounds=(0,None),method='highs')
        mx=linprog(-cost,A_eq=A,b_eq=b,bounds=(0,None),method='highs')
        assert mn.success and mx.success
        out[name]=[mn.fun,-mx.fun]
    return out

def run(name,dirs):
    L=3;N=9;rank=[torus_rank(m,L,dirs) for m in range(1<<N)];c=static_coeff(rank,N)
    hist={d:0 for d in range(N+1)};pairs={}
    for perm in itertools.permutations(range(N)):
        mask=0;j1=j2=None
        for step,v in enumerate(perm,1):
            mask|=1<<v;r=rank[mask]
            if j1 is None and r>=1:j1=step
            if r==2:j2=step;break
        assert j1 is not None and j2 is not None and j1<=j2
        d=j2-j1;hist[d]+=1;pairs[(j1,j2)]=pairs.get((j1,j2),0)+1
    tot=math.factorial(N)
    ED=Fraction(sum(d*n for d,n in hist.items()),tot)
    ED2=Fraction(sum(d*d*n for d,n in hist.items()),tot)
    PEQ=Fraction(hist[0],tot)
    static_ED=sum(Fraction(c[1][k],math.comb(N,k)) for k in range(N+1))
    p1,p2=marginals_from_static(c)
    return {'model':name,'L':3,'permutations':tot,'D_histogram':{str(k):v for k,v in hist.items() if v},
            'P_D0':str(PEQ),'E_D':str(ED),'E_D2':str(ED2),'Var_D':str(ED2-ED*ED),
            'static_E_D':str(static_ED),'static_ED_agrees':ED==static_ED,
            'J1_marginal':[str(x) for x in p1],'J2_marginal':[str(x) for x in p2],
            'marginal_only_LP_bounds':lp_bounds(c),
            'joint_nonzero_cells':len(pairs)}

out={'scope':'exact all 9! site permutations on L=3 square and triangular tori; LP bounds use only static birth marginals',
     'systems':[run('square',SQ),run('triangular',TRI)]}
print(json.dumps(out,indent=2))
