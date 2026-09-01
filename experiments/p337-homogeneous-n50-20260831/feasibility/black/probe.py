#!/usr/bin/env python3
"""Bounded exact rank-only lifted frontier; no Monte Carlo or source scoring."""
import argparse
from collections import defaultdict
import json
from math import gcd
from pathlib import Path
import resource
import signal
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'dual-color/scripts'))
from integer_period_torus import integer_torus_geometry, classify_configuration


def canonical(r, hx, hy, k, labels, px, py, retained):
    mapping, anchors, out = {}, {}, []
    for i in retained:
        c = labels[i]
        if c == 0:
            out.extend((0, 0, 0))
            continue
        if c not in mapping:
            mapping[c] = len(mapping) + 1
            anchors[c] = px[i], py[i]
        ax, ay = anchors[c]
        out.extend((mapping[c], px[i]-ax, py[i]-ay))
    return r, hx, hy, k, tuple(out)


def prepare(a, b):
    geo = integer_torus_geometry(((a, -b), (b, a)))
    n = geo.n
    prev = [[] for _ in range(n)]
    last = list(range(n))
    for e in geo.primal_edges:
        assert e.i != e.j
        u, v, dx, dy = e.i, e.j, e.dx, e.dy
        if u > v:
            u, v, dx, dy = v, u, -dx, -dy
        prev[v].append((u, dx, dy))
        last[u] = max(last[u], v)
    fs = [tuple(u for u in range(v) if last[u] >= v) for v in range(n+1)]
    return geo, prev, fs


def frontier(a, b, cap=100000, seconds=20, mib=600):
    geo, prev, fs = prepare(a,b)
    n = geo.n
    start = time.process_time()
    states = {(0,0,0,0,()): (1,2*n+1)}
    result = dict(geometry=[a,b],N=n,max_width=max(map(len,fs)),layers=[],complete=False)
    class Limit(Exception): pass
    def guard(*_):
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform=='darwin' else 1024)
        if time.process_time()-start > seconds or rss > mib:
            raise Limit('cpu_or_rss')
    signal.signal(signal.SIGALRM, guard)
    signal.setitimer(signal.ITIMER_REAL, .1, .1)
    try:
        for v in range(n):
            vertices = fs[v]+(v,)
            positions = {u:i for i,u in enumerate(vertices)}
            retained = [positions[u] for u in fs[v+1]]
            edges = [(positions[u],dx,dy) for u,dx,dy in prev[v]]
            following = {}
            for (oldr, oldhx, oldhy, k, packed), (count, sum_s) in states.items():
                basel = list(packed[0::3]); basex=list(packed[1::3]); basey=list(packed[2::3])
                for occupied in (0,1):
                    r,hx,hy=oldr,oldhx,oldhy
                    labels=basel+[max(basel,default=0)+1 if occupied else 0]
                    px=basex+[0]; py=basey+[0]
                    ds=-3*occupied
                    if occupied:
                        for u,ex,ey in edges:
                            if not labels[u]: continue
                            if r==1: ex,ey=hx*ey-hy*ex,0
                            if r==2: ex,ey=0,0
                            dx,dy=px[u]+ex-px[-1],py[u]+ey-py[-1]
                            cu,cv=labels[u],labels[-1]
                            if cu!=cv:
                                for j in range(len(labels)):
                                    if labels[j]==cv:
                                        labels[j]=cu;px[j]+=dx;py[j]+=dy
                            else:
                                ds+=2
                                if r==0 and (dx or dy):
                                    div=gcd(abs(dx),abs(dy));hx,hy=dx//div,dy//div
                                    if hx<0 or (hx==0 and hy<0):hx,hy=-hx,-hy
                                    px=[hx*y-hy*x for x,y in zip(px,py)];py=[0]*len(py)
                                    r=1;ds-=1
                                elif r==1 and dx:
                                    r=2;hx=hy=0;px=[0]*len(px);py=[0]*len(py);ds-=1
                    key=canonical(r,hx,hy,k+occupied,labels,px,py,retained)
                    val=following.get(key)
                    addition=sum_s+ds*count
                    following[key]=(count,addition) if val is None else (val[0]+count,val[1]+addition)
            states=following
            row=dict(layer=v+1,frontier=len(fs[v+1]),states=len(states),cpu_seconds=time.process_time()-start)
            result['layers'].append(row)
            print(json.dumps(row),flush=True)
            if len(states)>cap:raise Limit('state_cap')
        out=defaultdict(lambda:[0,0])
        for (r,hx,hy,k,packed),(c,s) in states.items():
            assert not packed
            out[(k,r-1)][0]+=c;out[(k,r-1)][1]+=s
        result['histogram']=[dict(K=k,q=q,count=c,sum_S=s) for (k,q),(c,s) in sorted(out.items())]
        assert sum(v[0] for v in out.values())==2**n
        result['complete']=True
    except Limit as e: result['stop']=str(e)
    finally:signal.setitimer(signal.ITIMER_REAL,0)
    result['cpu_seconds']=time.process_time()-start
    result['rss_mib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform=='darwin' else 1024)
    return result


def oracle(a,b):
    geo,_,_=prepare(a,b);n=geo.n
    out=defaultdict(lambda:[0,0])
    for mask in range(1<<n):
        active=[bool(mask>>i&1) for i in range(n)]
        _,components=classify_configuration(geo,active)
        vectors=[v for c in components for v in c.generators]
        r=0
        for x,y in vectors:
            if not (x or y):continue
            if r==0:r=1;ux,uy=x,y
            elif ux*y-uy*x:r=2;break
        k=mask.bit_count();edges=sum(active[e.i] and active[e.j] for e in geo.primal_edges)
        beta=edges-k+len(components)
        s=2*beta-r-3*k+2*n+1
        out[(k,r-1)][0]+=1;out[(k,r-1)][1]+=s
    return [dict(K=k,q=q,count=c,sum_S=s) for (k,q),(c,s) in sorted(out.items())]


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('a',type=int);ap.add_argument('b',type=int)
    ap.add_argument('--cap',type=int,default=100000);ap.add_argument('--seconds',type=int,default=20)
    ap.add_argument('--mib',type=int,default=600);ap.add_argument('--oracle',action='store_true');ap.add_argument('--output',required=True)
    args=ap.parse_args();resource.setrlimit(resource.RLIMIT_CPU,(args.seconds+30,args.seconds+30))
    result=frontier(args.a,args.b,args.cap,args.seconds,args.mib)
    if args.oracle:
        assert result['complete'];expected=oracle(args.a,args.b)
        assert result['histogram']==expected
        result['direct_oracle_equal']=True
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
