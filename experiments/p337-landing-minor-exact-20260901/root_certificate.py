#!/usr/bin/env python3
"""Exact Sturm certificate for the L4 landing-minor at the matching root."""
from __future__ import annotations
import argparse,json,math
from fractions import Fraction as F
from pathlib import Path
import landing_minor as lm

def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return tuple(a or [F(0)])
def add(a,b):
    n=max(len(a),len(b));return trim([(a[i] if i<len(a) else F(0))+(b[i] if i<len(b) else F(0)) for i in range(n)])
def scale(a,c):return trim([c*x for x in a])
def sub(a,b):return add(a,scale(b,F(-1)))
def mul(a,b):
    z=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):z[i+j]+=x*y
    return trim(z)
def derivative(a):return trim([F(i)*a[i] for i in range(1,len(a))]) if len(a)>1 else (F(0),)
def evaluate(a,x):
    z=F(0)
    for c in reversed(a):z=z*x+c
    return z
def divmod_poly(a,b):
    a=list(trim(a));b=trim(b)
    if b==(F(0),):raise ZeroDivisionError
    if len(a)<len(b):return (F(0),),tuple(a)
    q=[F(0)]*(len(a)-len(b)+1)
    while len(a)>=len(b) and any(a):
        d=len(a)-len(b);c=a[-1]/b[-1];q[d]=c
        for i,x in enumerate(b):a[d+i]-=c*x
        a=list(trim(a))
    return trim(q),trim(a)
def sturm(a):
    s=[trim(a),derivative(a)]
    while s[-1]!=(F(0),):
        _,r=divmod_poly(s[-2],s[-1])
        if r==(F(0),):break
        s.append(scale(r,F(-1)))
    return tuple(s)
def variations(s,x):
    q=[(evaluate(a,x)>0)-(evaluate(a,x)<0) for a in s];q=[x for x in q if x]
    return sum(a!=b for a,b in zip(q,q[1:]))
def roots(a,l,r):
    if evaluate(a,l)==0 or evaluate(a,r)==0:raise ValueError('endpoint root')
    s=sturm(a);return variations(s,l)-variations(s,r)
def bernstein(layer,n):
    z=(F(0),)
    for k,v in enumerate(layer):
        if not v:continue
        z=add(z,mul((F(0),)*k+(v,),tuple(F(math.comb(n-k,j)*(-1)**j) for j in range(n-k+1))))
    return z

def build():
    t=lm.Torus(4);n=16;src=[t.source(m) for m in range(1<<n)];rank=[t.rank(m) for m in range(1<<n)]
    sl=[F(0)]*(n+1)
    for m,a in enumerate(src):sl[m.bit_count()]+=a
    mean=bernstein(sl,n);trans=((0,1),(1,2));hl={q:[F(0)]*n for q in trans};hal={q:[F(0)]*n for q in trans}
    for m in range(1<<n):
        if m&1:continue
        q=rank[m],rank[m|1];h=t.landing(m,1)[0]
        if q not in hl or not h:continue
        k=m.bit_count();hl[q][k]+=h;hal[q][k]+=h*(src[m]+src[m|1])/2
    matrix={}
    for q in trans:
        H=bernstein(hl[q],n-1);HK=bernstein([F(2*k+1,2)*hl[q][k] for k in range(n)],n-1)
        T=sub(HK,mul((F(0),F(n)),H));A=sub(bernstein(hal[q],n-1),mul(mean,H));matrix[q]=T,A
    T01,A01=matrix[(0,1)];T12,A12=matrix[(1,2)];det=sub(mul(T01,A12),mul(T12,A01));tsum=add(T01,T12)
    matching=(F(-1),F(0),F(0),F(0),F(8),F(0),F(32),F(-64),F(172),F(-704),F(1104),F(-608),F(-56),F(128),F(16),F(-32),F(6))
    left,right=F(59,100),F(3,5);lo,hi=left,right;sgn=(evaluate(matching,lo)>0)-(evaluate(matching,lo)<0)
    for _ in range(180):
        mid=(lo+hi)/2;s=(evaluate(matching,mid)>0)-(evaluate(matching,mid)<0)
        if s==sgn:lo=mid
        else:hi=mid
    p=(lo+hi)/2
    return {'schema':'matching-one/p337-landing-root-certificate/v1','interval':[str(left),str(right)],'matching_root_count':roots(matching,left,right),'determinant_root_count':roots(det,left,right),'thermal_sum_root_count':roots(tsum,left,right),'matching_signs':[str(evaluate(matching,left)),str(evaluate(matching,right))],'determinant_signs':[str(evaluate(det,left)),str(evaluate(det,right))],'thermal_sum_signs':[str(evaluate(tsum,left)),str(evaluate(tsum,right))],'root_midpoint_decimal':f'{float(p):.16g}','matrix_at_root_midpoint':{'T_01':f'{float(evaluate(T01,p)):.16g}','T_12':f'{float(evaluate(T12,p)):.16g}','A_01':f'{float(evaluate(A01,p)):.16g}','A_12':f'{float(evaluate(A12,p)):.16g}','determinant':f'{float(evaluate(det,p)):.16g}','thermal_sum':f'{float(evaluate(tsum,p)):.16g}','root_schur_even_residual':f'{float(2*evaluate(det,p)/evaluate(tsum,p)):.16g}'},'half_root_schur_even_residual':str(2*evaluate(det,F(1,2))/evaluate(tsum,F(1,2))),'determinant_coefficients':[str(x) for x in det]}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);a=ap.parse_args();z=build();assert z['matching_root_count']==1 and z['determinant_root_count']==0 and z['thermal_sum_root_count']==0 and all(F(x)<0 for x in z['determinant_signs']+z['thermal_sum_signs']);text=json.dumps(z,indent=2)+'\n';print(text,end='')
    if a.out:
        if a.out.exists():raise SystemExit('refusing overwrite')
        a.out.write_text(text)
if __name__=='__main__':main()
