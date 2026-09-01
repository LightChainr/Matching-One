#!/usr/bin/env python3
"""Exact finite falsifier for #537's pure-thermal landing rank-one gate."""
from __future__ import annotations
import argparse,json,math,time
from collections import deque
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from pathlib import Path

@lru_cache(None)
def parts(n):
    out=[(0,)]
    for _ in range(1,n): out=[p+(j,) for p in out for j in range(max(p)+2)]
    return tuple(out)
def canon(p):
    d={}; return tuple(d.setdefault(x,len(d)) for x in p)
def local4(p):
    z=0
    for s in (0,1):
        a,b,c,d=p[s:]+p[:s]
        if a==b or c==d: continue
        z+=int(a==c and b==d)+int(a==d and b==c)
        z+=int(a==c)+int(a==d)+int(b==c)+int(b==d)-4
    return z
@lru_cache(None)
def _g16(p):
    z=0
    for q in parts(max(p)+1):
        k=max(q)+1
        if k==1: continue
        c=tuple(q[x] for x in p)
        z+=(-1)**(k-2)*math.factorial(k-2)*local4(c[:4])*local4(c[4:])
    return z
def g16(p): return _g16(canon(tuple(p)))

class Torus:
    step=((0,1),(1,0),(0,-1),(-1,0)); mstep=step+((1,1),(1,-1),(-1,1),(-1,-1))
    def __init__(self,L):
        if L<4: raise ValueError('L>=4')
        self.L=L; self.N=L*L
        self.nb=[tuple(self.v(x+dx,y+dy) for dx,dy in self.step) for y in range(L) for x in range(L)]
        self.port=[(2*i,2*i+1,2*self.nb[i][2],2*self.nb[i][3]+1) for i in range(self.N)]
    def v(self,x,y): return x%self.L+self.L*(y%self.L)
    def xy(self,i): return i%self.L,i//self.L
    def mask(self,pts):
        m=0
        for x,y in pts:m|=1<<self.v(x,y)
        return m
    def rank(self,m):
        lift={}; cyc=[]
        for v in range(self.N):
            if not(m>>v)&1 or v in lift: continue
            lift[v]=(0,0); q=deque([v])
            while q:
                u=q.popleft(); ux,uy=lift[u]
                for w,(dx,dy) in zip(self.nb[u],self.step):
                    if not(m>>w)&1: continue
                    t=ux+dx,uy+dy
                    if w not in lift: lift[w]=t; q.append(w)
                    else:
                        c=t[0]-lift[w][0],t[1]-lift[w][1]
                        if c!=(0,0):cyc.append(c)
        if not cyc:return 0
        a,b=cyc[0]; return 2 if any(a*d-b*c for c,d in cyc[1:]) else 1
    def roots(self,m):
        p=list(range(2*self.N))
        def r(x):
            while p[x]!=x:p[x]=p[p[x]];x=p[x]
            return x
        for z,e in enumerate(self.port):
            if not(m>>z)&1: continue
            a=r(e[0])
            for x in e[1:]:
                b=r(x)
                if a!=b:p[b]=a
        return [r(i) for i in range(2*self.N)]
    def source(self,m):
        r=self.roots(m); lab=[tuple(r[e] for e in q) for q in self.port]; s=0
        for x,y in combinations(range(self.N),2):
            if (m>>x)&1 or (m>>y)&1: continue
            s+=g16(lab[x]+lab[y])
        return F(2*s,16*self.N*self.N)
    @staticmethod
    def sector(x,y,shift=0):return (math.floor((math.atan2(y,x)+math.pi/8)/(math.pi/4))-shift)%8
    @staticmethod
    def pair(ms,a,b):return any(i!=j and ms[i]&(1<<a) and ms[j]&(1<<b) for i in range(len(ms)) for j in range(len(ms)))
    def comps(self,m,R,matching,on,shift=0):
        pts=[(x,y) for y in range(-R,R+1) for x in range(-R,R+1) if (x,y)!=(0,0)]
        vs={q:self.v(*q) for q in pts}
        if len(set(vs.values()))!=len(vs):raise ValueError('noninjective landing box')
        st=self.mstep if matching else self.step; unseen={q for q,v in vs.items() if bool((m>>v)&1)==on}; out=[]
        while unseen:
            a=unseen.pop(); c={a}; q=[a]
            while q:
                x,y=q.pop()
                for dx,dy in st:
                    b=x+dx,y+dy
                    if b in unseen:unseen.remove(b);c.add(b);q.append(b)
            if not any(x in set(st) for x in c):continue
            z=0
            for x,y in c:
                if max(abs(x),abs(y))==R:z|=1<<self.sector(x,y,shift)
            if z:out.append(z)
        return out
    def landing(self,m,R,shift=0):
        o=self.comps(m,R,False,True,shift); c=self.comps(m,R,True,False,shift); P=self.pair
        ax=(P(o,0,4) and P(c,2,6)) or (P(o,2,6) and P(c,0,4))
        di=(P(o,1,5) and P(c,3,7)) or (P(o,3,7) and P(c,1,5))
        return int(ax)-int(di),tuple(o),tuple(c)

def family(L,R,p=F(1,2)):
    if R<1 or L<2*R+5:raise ValueError('R>=1,L>=2R+5')
    t=Torus(L); s=R+2; bit=1
    A=t.mask([(x,0) for x in range(1,L)]+[(x,s) for x in range(1,L)])
    B=t.mask([(0,y) for y in range(1,L)]+[(x,s) for x in range(L)])
    rows=[]
    for name,m,tr in [('A',A,(0,1)),('B',B,(1,2))]:
        h,o,c=t.landing(m,R); a0=t.source(m);a1=t.source(m|bit)
        rows.append(dict(name=name,transition=[t.rank(m),t.rank(m|bit)],expected=list(tr),h4=h,k=m.bit_count(),opened=o,closed=c,a0=str(a0),a1=str(a1),amid=str((a0+a1)/2)))
    k=2*L-2; w=p**k*(1-p)**((L-1)**2); S=F(2*k+1,2)-L*L*p; d=F(-2,L**4)
    return dict(L=L,R=R,p=str(p),states=rows,thermal=str(S),source_mid_B_minus_A=str(d),minor=str(w*w*S*d))

def l4():
    t=Torus(4); n=16; bit=1; src=[t.source(m) for m in range(1<<n)]; rank=[t.rank(m) for m in range(1<<n)]; mean=sum(src,F(0))/2**n
    M={(0,1):[F(0),F(0)],(1,2):[F(0),F(0)]}; cross=dict(pivotal=0,axis=0,diagonal=0,both=0,landed=0,h4=0); shiftbad=0; states={(0,1):0,(1,2):0}
    for m in range(1<<n):
        if m&bit:continue
        on=m|bit; tr=rank[m],rank[on]; h,o,c=t.landing(m,1); hs,_,_=t.landing(m,1,1)
        shiftbad+=hs!=-h
        piv=int(tr[0]<2 and tr[1]==2);cross['pivotal']+=piv
        if piv:
            ax=int(h!=0 or (t.pair(o,0,4) and t.pair(c,2,6)) or (t.pair(o,2,6) and t.pair(c,0,4)))
            # recover exact axis/diagonal from h and both by direct predicates
            A=(t.pair(o,0,4) and t.pair(c,2,6)) or (t.pair(o,2,6) and t.pair(c,0,4));D=(t.pair(o,1,5) and t.pair(c,3,7)) or (t.pair(o,3,7) and t.pair(c,1,5))
            cross['axis']+=int(A);cross['diagonal']+=int(D);cross['both']+=int(A and D);cross['landed']+=int(A or D);cross['h4']+=h
        if tr not in M or h==0:continue
        k=m.bit_count(); S=F(2*k+1,2)-8; amid=(src[m]+src[on])/2
        M[tr][0]+=F(h)*S/2**15;M[tr][1]+=F(h)*(amid-mean)/2**15;states[tr]+=1
    d=M[(0,1)][0]*M[(1,2)][1]-M[(1,2)][0]*M[(0,1)][1]
    return dict(cross=cross,expected=dict(pivotal=3121,axis=892,diagonal=474,both=88,landed=1278,h4=418),shift_violations=shiftbad,state_counts={'0_to_1':states[(0,1)],'1_to_2':states[(1,2)]},matrix={'T_01':str(M[(0,1)][0]),'T_12':str(M[(1,2)][0]),'A_01':str(M[(0,1)][1]),'A_12':str(M[(1,2)][1]),'det':str(d)})

def selftest():
    z=l4();assert z['cross']==z['expected'] and z['shift_violations']==0
    assert z['matrix']['det']=='-533831111/140737488355328'
    for L,R in ((7,1),(9,2),(11,3),(13,4)):
        q=family(L,R);A,B=q['states'];assert A['transition']==[0,1] and B['transition']==[1,2] and A['h4']==B['h4']==1
        assert F(B['amid'])-F(A['amid'])==F(-2,L**4) and F(q['minor'])!=0
    return z

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);a=ap.parse_args();t=time.perf_counter();z=selftest();payload={'schema':'matching-one/p337-landing-minor/v1','decision':'FINITE_PURE_THERMAL_RANK_ONE_REJECTED','boundary':'finite/source-Schur algebra only; no asymptotic lower bound or T_N rate','family':[family(7,1),family(9,2),family(11,3),family(13,4)],'l4':z,'wall_seconds':time.perf_counter()-t};text=json.dumps(payload,indent=2)+'\n';print(text,end='')
    if a.out:
        if a.out.exists():raise SystemExit('refusing overwrite')
        a.out.write_text(text)
if __name__=='__main__':main()
