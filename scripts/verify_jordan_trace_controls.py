#!/usr/bin/env python3
"""Independent Fraction-arithmetic recheck; no SymPy or producer imports."""
from fractions import Fraction as F
from pathlib import Path
import argparse,json,time


def mm(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),F(0))
             for j in range(len(b[0]))] for i in range(len(a))]
def eye(n):return [[F(i==j) for j in range(n)] for i in range(n)]
def add(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def scale(a,c):return [[c*x for x in row] for row in a]
def trace(a):return sum(a[i][i] for i in range(len(a)))
def determinant(a):
    a=[row[:] for row in a];n=len(a);sgn=1;out=F(1)
    for c in range(n):
        piv=next((r for r in range(c,n) if a[r][c]),None)
        if piv is None:return F(0)
        if piv!=c:a[c],a[piv]=a[piv],a[c];sgn=-sgn
        q=a[c][c];out*=q
        for r in range(c+1,n):
            f=a[r][c]/q
            for j in range(c,n):a[r][j]-=f*a[c][j]
    return sgn*out


def verify(path):
    start=time.perf_counter();j=json.loads(Path(path).read_text());I=eye(3)
    P=[[F(1,3)]*3 for _ in range(3)]
    D=[[F(1,3),F(-2,3),F(1,3)],[F(-2,3),F(1,3),F(1,3)],[F(1,3),F(1,3),F(-2,3)]]
    N=[[F(1,6),F(1,6),F(-1,3)],[F(-1,6),F(-1,6),F(1,3)],[F(0)]*3]
    assert [[str(x) for x in row] for row in N]==j['N']
    matrices=0;moments=0
    for t in (F(-1,4),F(-1,8),F(0),F(1,8),F(1,4)):
        G0=add(add(P,scale(I,-1)),scale(D,t));G1=add(G0,scale(N,F(1,2)))
        for G in (G0,G1):
            assert all(sum(row)==0 for row in G)
            assert all(sum(G[i][k] for i in range(3))==0 for k in range(3))
            assert all(G[i][k]>0 for i in range(3) for k in range(3) if i!=k)
            assert determinant(G)==0
            matrices+=1
        T0=add(I,scale(G0,F(1,2)));T1=add(I,scale(G1,F(1,2)))
        A0=I;A1=I
        for m in range(21):
            expected=1+((1+t)/2)**m+((1-t)/2)**m
            assert trace(A0)==expected==trace(A1)
            if t==0:
                assert A0[0][2]==(1-F(1,2)**m)/3
                assert A1[0][2]==(1-F(1,2)**m)/3-F(m,6)*F(1,2)**m
            A0=mm(A0,T0);A1=mm(A1,T1);moments+=1
    yd=[F(x) for x in j['probability_readout']['semisimple_values']]
    yj=[F(x) for x in j['probability_readout']['Jordan_values']]
    a=determinant([[yd[i+k] for k in range(2)] for i in range(2)])
    b=determinant([[yj[i+k] for k in range(3)] for i in range(3)])
    assert [str(a),str(b)]==j['probability_readout']['hankel_determinants']
    assert a and b
    return {'schema':'matching-one.jordan-trace-fraction-verification.v1',
            'arithmetic':'stdlib Fraction; independently reconstructed matrices; Gaussian determinant',
            'generator_instances':matrices,'trace_moments':moments,
            'hankel_determinants':[str(a),str(b)],'success':True,
            'elapsed_seconds':time.perf_counter()-start}

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--out',type=Path)
    a=ap.parse_args();text=json.dumps(verify(a.input),indent=2)+'\n'
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        with a.out.open('x') as f:f.write(text)
    else:print(text,end='')
