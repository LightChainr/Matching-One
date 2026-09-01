#!/usr/bin/env python3
"""Independent bounded verifier for the exact L4 landing-minor certificate."""
from __future__ import annotations
import argparse,json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import landing_minor as primary

class PotentialUnionFind:
    def __init__(self,length:int):
        self.length=length;self.n=length*length;self.parent=list(range(self.n));self.size=[1]*self.n;self.dx=[0]*self.n;self.dy=[0]*self.n;self.cycles=[set() for _ in range(self.n)]
    def find(self,v):
        if self.parent[v]==v:return v,0,0
        r,px,py=self.find(self.parent[v]);dx=self.dx[v]+px;dy=self.dy[v]+py;self.parent[v]=r;self.dx[v]=dx;self.dy[v]=dy;return r,dx,dy
    def add(self,a,b,ex,ey):
        ra,ax,ay=self.find(a);rb,bx,by=self.find(b);rx,ry=ax+ex-bx,ay+ey-by
        if ra==rb:
            if rx%self.length or ry%self.length:raise AssertionError('nonperiod displacement')
            w=rx//self.length,ry//self.length
            if w!=(0,0):self.cycles[ra].add(w)
            return
        if self.size[ra]<self.size[rb]:ra,rb=rb,ra;rx,ry=-rx,-ry
        self.parent[rb]=ra;self.dx[rb]=rx;self.dy[rb]=ry;self.size[ra]+=self.size[rb];self.cycles[ra]|=self.cycles[rb]
    def rank(self,mask):
        L=self.length
        for y in range(L):
            for x in range(L):
                a=x+L*y
                if not(mask>>a)&1:continue
                for ex,ey in ((1,0),(0,1)):
                    b=((x+ex)%L)+L*((y+ey)%L)
                    if (mask>>b)&1:self.add(a,b,ex,ey)
        roots={self.find(v)[0] for v in range(self.n) if (mask>>v)&1};vec=[w for r in roots for w in self.cycles[r]]
        if not vec:return 0
        x,y=vec[0];return 2 if any(x*v-y*u for u,v in vec[1:]) else 1

def bfs_labels(t,mask):
    adj=[set() for _ in range(2*t.N)]
    for z,inc in enumerate(t.port):
        if not(mask>>z)&1:continue
        for a,b in combinations(inc,2):adj[a].add(b);adj[b].add(a)
    lab=[-1]*len(adj);c=0
    for s in range(len(adj)):
        if lab[s]>=0:continue
        lab[s]=c;q=[s]
        while q:
            u=q.pop()
            for v in adj[u]:
                if lab[v]<0:lab[v]=c;q.append(v)
        c+=1
    return lab

def source_bfs(t,mask):
    lab=bfs_labels(t,mask);site=[tuple(lab[e] for e in inc) for inc in t.port];s=0
    for x,y in combinations(range(t.N),2):
        if (mask>>x)&1 or (mask>>y)&1:continue
        s+=primary.g16(site[x]+site[y])
    return F(2*s,16*t.N*t.N)

def verify():
    t=primary.Torus(4);src=[];rank=[]
    for m in range(1<<16):
        a=t.source(m);b=source_bfs(t,m)
        if a!=b:raise AssertionError(('source',m,a,b))
        r=t.rank(m);q=PotentialUnionFind(4).rank(m)
        if r!=q:raise AssertionError(('rank',m,r,q))
        src.append(b);rank.append(q)
    mean=sum(src,F(0))/2**16;M={(0,1):[F(0),F(0)],(1,2):[F(0),F(0)]}
    for m in range(1<<16):
        if m&1:continue
        tr=rank[m],rank[m|1]
        if tr not in M:continue
        h=t.landing(m,1)[0]
        if not h:continue
        S=F(2*m.bit_count()+1,2)-8;amid=(src[m]+src[m|1])/2
        M[tr][0]+=F(h)*S/2**15;M[tr][1]+=F(h)*(amid-mean)/2**15
    d=M[(0,1)][0]*M[(1,2)][1]-M[(1,2)][0]*M[(0,1)][1]
    expected=F(-533831111,140737488355328)
    if d!=expected:raise AssertionError((d,expected))
    return {'states_checked':1<<16,'source_algorithm':'physical-edge adjacency BFS','rank_algorithm':'potential union-find','matrix':{'T_01':str(M[(0,1)][0]),'T_12':str(M[(1,2)][0]),'A_01':str(M[(0,1)][1]),'A_12':str(M[(1,2)][1]),'det':str(d)},'status':'verified'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);a=ap.parse_args();z=verify();text=json.dumps(z,indent=2)+'\n';print(text,end='')
    if a.out:
        if a.out.exists():raise SystemExit('refusing overwrite')
        a.out.write_text(text)
if __name__=='__main__':main()
