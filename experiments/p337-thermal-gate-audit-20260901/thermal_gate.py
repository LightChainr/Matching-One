#!/usr/bin/env python3
"""Exact finite controls for the canonical thermal-pivotal gate. Standard library.
No Monte Carlo, critical exponent fit, or infinite-volume inference is performed.
"""
from __future__ import annotations
import argparse, hashlib, json, math, time
from collections import Counter, deque
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from pathlib import Path

@lru_cache(None)
def partitions(n):
    if n < 1: raise ValueError('positive n required')
    out = [(0,)]
    for _ in range(1,n):
        out = [p+(j,) for p in out for j in range(max(p)+2)]
    return tuple(out)

def canon(p):
    names={}
    return tuple(names.setdefault(x,len(names)) for x in p)

def local4(p):
    """4*Kreg at Q1 from two C4 recouplings, canonical completion only."""
    total=0
    for s in (0,1):
        a,b,c,d=p[s:]+p[:s]
        if a==b or c==d: continue
        total += int(a==c and b==d)+int(a==d and b==c)
        total += int(a==c)+int(a==d)+int(b==c)+int(b==d)-4
    return total

@lru_cache(None)
def _g16(p):
    total=0
    for rho in partitions(max(p)+1):
        k=max(rho)+1
        if k==1: continue
        c=tuple(rho[x] for x in p)
        total += (-1)**(k-2)*math.factorial(k-2)*local4(c[:4])*local4(c[4:])
    return total

def g16(p): return _g16(canon(p))

def nshared(p): return len(set(p[:4]) & set(p[4:]))

def merge_pattern(p, selected):
    s=set(selected)
    if not s: return canon(p)
    v=min(s)
    return canon(tuple(v if x in s else x for x in p))

def join_audit():
    hist=Counter(); nonzero=Counter(); sign=Counter(); support_fail=0; checks=0
    maxabs=0; witness=None; changes_without_support_change=None
    for p in partitions(8):
        old=g16(p); so=nshared(p)
        for m in range(2,min(4,max(p)+1)+1):
            for S in combinations(range(max(p)+1),m):
                q=merge_pattern(p,S); new=g16(q); sn=nshared(q); d=new-old
                checks+=1; hist[(so,sn)]+=1
                if d:
                    nonzero[(so,sn)]+=1; sign['positive' if d>0 else 'negative']+=1
                    if so==0 or max(so,sn)<2: support_fail+=1
                    if abs(d)>maxabs: maxabs=abs(d);witness=(p,S,q,old,new,d)
                    if changes_without_support_change is None and so==sn and old*new<0:
                        changes_without_support_change=(p,S,q,old,new,d)
    assert support_fail==0
    return dict(checks=checks,nonzero=sum(nonzero.values()),sign=dict(sign),
        nonzero_by_shared_transition={f'{i}->{j}':v for (i,j),v in sorted(nonzero.items())},
        support_failures=support_fail,max_abs_delta16=maxabs,max_witness=witness,
        sign_flip_same_shared_count=changes_without_support_change)

class Torus:
    def __init__(self,a,b):
        self.a=a;self.b=b;self.n=a*a+b*b
        n=self.n
        key=lambda x,y: ((a*x+b*y)%n,(-b*x+a*y)%n)
        self.reps=[(0,0)]; idx={key(0,0):0}
        for x,y in self.reps:
            for dx,dy in ((1,0),(0,1)):
                k=key(x+dx,y+dy)
                if k not in idx: idx[k]=len(self.reps);self.reps.append((x+dx,y+dy))
        assert len(self.reps)==n
        self.step=((0,1),(1,0),(0,-1),(-1,0))
        self.nb=[tuple(idx[key(x+dx,y+dy)] for dx,dy in self.step) for x,y in self.reps]
        if any(len(set(v))<4 for v in self.nb): raise ValueError('aliased NN ports')
        self.ports=[(2*i,2*i+1,2*self.nb[i][2],2*self.nb[i][3]+1) for i in range(n)]
        self.white=[tuple(idx[key(x+dx,y+dy)] for dx in (-1,0,1) for dy in (-1,0,1) if dx or dy) for x,y in self.reps]
        self.faces=[(i,idx[key(x+1,y)],idx[key(x,y+1)],idx[key(x+1,y+1)]) for i,(x,y) in enumerate(self.reps)]
        self.pairs=list(combinations(range(n),2))
    def edge_roots(self,mask):
        par=list(range(2*self.n))
        def root(x):
            while par[x]!=x: par[x]=par[par[x]];x=par[x]
            return x
        for z,ed in enumerate(self.ports):
            if (mask>>z)&1:
                x=root(ed[0])
                for e in ed[1:]: par[root(e)]=x
        return [root(x) for x in range(2*self.n)]
    def rank(self,mask):
        lift={};cycles=[]
        for v in range(self.n):
            if not (mask>>v)&1 or v in lift: continue
            lift[v]=(0,0);queue=deque([v])
            while queue:
                x=queue.popleft();lx,ly=lift[x]
                for y,(dx,dy) in zip(self.nb[x],self.step):
                    if not (mask>>y)&1: continue
                    target=(lx+dx,ly+dy)
                    if y not in lift: lift[y]=target;queue.append(y)
                    else:
                        c=(target[0]-lift[y][0],target[1]-lift[y][1])
                        if c!=(0,0): cycles.append(c)
        if not cycles: return 0
        x,y=cycles[0]
        return 2 if any(x*v-y*u for u,v in cycles[1:]) else 1
    def components(self,mask,occupied,neighbors):
        seen=set();count=0
        for v in range(self.n):
            if bool((mask>>v)&1)!=occupied or v in seen: continue
            count+=1;seen.add(v);stack=[v]
            while stack:
                x=stack.pop()
                for y in neighbors[x]:
                    if bool((mask>>y)&1)==occupied and y not in seen:seen.add(y);stack.append(y)
        return count
    def dual_rank(self,mask):
        k=mask.bit_count()
        edges=sum(((mask>>i)&1)*sum((mask>>j)&1 for j in self.nb[i][:2]) for i in range(self.n))
        faces=sum(all((mask>>j)&1 for j in f) for f in self.faces)
        return self.components(mask,True,self.nb)-self.components(mask,False,self.white)-(k-edges+faces)+1
    def pair_values(self,mask,roots=None):
        roots=self.edge_roots(mask) if roots is None else roots
        labels=[tuple(roots[e] for e in ed) for ed in self.ports]
        return [0 if ((mask>>x)&1 or (mask>>y)&1) else g16(labels[x]+labels[y]) for x,y in self.pairs]

def fraction_entry(x):
    x=F(x);return {'exact':str(x),'decimal':float(x)}

def conditional_covariance_decomposition(n, obs, src16, p):
    """Derivative in p; src16 = 16*N^2*source (ordered pair integer sum)."""
    scale=16*n*n
    full=[p**k*(1-p)**(n-k) for k in range(n+1)]
    cond=[p**k*(1-p)**(n-1-k) for k in range(n)]
    meanO=sum((full[m.bit_count()]*obs[m] for m in range(1<<n)),F(0))
    meanA=sum((full[m.bit_count()]*src16[m] for m in range(1<<n)),F(0))/scale
    dA=F(0);dO=F(0);dOA=F(0);rearr=F(0);readout=F(0)
    signed=F(0);absolute=F(0);obs_only_abs=F(0);source_only_abs=F(0)
    jumps=Counter(); failures=0
    for z in range(n):
        bit=1<<z
        for m in range(1<<n):
            if m&bit:continue
            m1=m|bit;wt=cond[m.bit_count()]
            o0,o1=F(obs[m]),F(obs[m1]);a0,a1=F(src16[m],scale),F(src16[m1],scale)
            do=o1-o0;da=a1-a0
            dA+=wt*da;dO+=wt*do;dOA+=wt*(o1*a1-o0*a0)
            u=((o0+o1)/2-meanO)*da
            v=((a0+a1)/2-meanA)*do
            rearr+=wt*u;readout+=wt*v
            if (o1*a1-o0*a0)!=(o1+o0)*da/2+(a1+a0)*do/2:failures+=1
            signed+=wt*da;absolute+=wt*abs(da)
            if da==0:obs_only_abs+=wt*abs(v)
            if do==0:source_only_abs+=wt*abs(u)
            if do:jumps[str(do)]+=1
    # Independent Bernoulli score/cumulant calculation.
    score=F(0)
    for m in range(1<<n):
        k=m.bit_count(); a=F(src16[m],scale)
        score+=full[k]*(obs[m]-meanO)*(a-meanA)*(k-n*p)/(p*(1-p))
    direct=dOA-meanO*dA-meanA*dO
    assert failures==0 and direct==rearr+readout==score
    return dict(p=str(p),mean_O=fraction_entry(meanO),mean_source=fraction_entry(meanA),
        source_derivative=fraction_entry(dA),covariance_derivative=fraction_entry(direct),
        kernel_rearrangement=fraction_entry(rearr),observable_pivotal=fraction_entry(readout),
        omitted_readout_on_delta_source_zero_abs=fraction_entry(obs_only_abs),
        omitted_source_on_delta_observable_zero_abs=fraction_entry(source_only_abs),
        abs_source_Russo_mass=fraction_entry(absolute),
        source_signed_to_absolute=fraction_entry(abs(dA)/absolute if absolute else 0),
        score_residual=str(score-direct),leibniz_failures=failures)

def enumerate_case(a,b,full_pair_edge_check=True):
    start=time.perf_counter();t=Torus(a,b);n=t.n
    if n>13:raise ValueError('this reference enumerator is bounded to N<=13')
    ranks=[];sums=[];states=[];top_fail=0
    for m in range(1<<n):
        r=t.rank(m);top_fail+=r!=t.dual_rank(m);ranks.append(r)
        vals=t.pair_values(m);states.append(vals);sums.append(2*sum(vals))
    assert top_fail==0
    checks=0;failed=0;zero_shared_failed=0;transition_hist=Counter();same_s_sign=0
    witness=None; endpoint_hist=[0]*n; bulk_hist=[0]*n
    if full_pair_edge_check:
        for z in range(n):
            bit=1<<z
            for m in range(1<<n):
                if m&bit:continue
                roots=t.edge_roots(m);touched={roots[e] for e in t.ports[z]};newid=min(touched)
                for j,(x,y) in enumerate(t.pairs):
                    if ((m>>x)&1 or (m>>y)&1):continue
                    p=tuple(roots[e] for e in t.ports[x]+t.ports[y])
                    old=states[m][j];new=states[m|bit][j];delta=new-old
                    if z==x or z==y:
                        predicted=-old; endpoint_hist[m.bit_count()]+=2*delta
                    else:
                        bulk_hist[m.bit_count()]+=2*delta
                        merged=tuple(newid if v in touched else v for v in p)
                        predicted=g16(merged)-old
                        so,sn=nshared(p),nshared(merged)
                        if delta:
                            transition_hist[f'{so}->{sn}']+=1
                            if so==0:zero_shared_failed+=1
                            if so==sn and old*new<0:
                                same_s_sign+=1
                                if witness is None:witness=dict(mask=m,z=z,pair=[x,y],reps=t.reps,old16=old,new16=new,shared=so,rank_before=ranks[m],rank_after=ranks[m|bit])
                    checks+=1;failed+=predicted!=delta
        assert failed==zero_shared_failed==0
    rows=[]
    for p in (F(1,2),F(3,5),F(2,3)):
        for name,obs in [('q',[r-1 for r in ranks]),('E',[(r-1)**2 for r in ranks])]:
            r=conditional_covariance_decomposition(n,obs,sums,p)
            if full_pair_edge_check:
                endpoint=sum((p**k*(1-p)**(n-1-k)*v for k,v in enumerate(endpoint_hist)),F(0))/(16*n*n)
                bulk=sum((p**k*(1-p)**(n-1-k)*v for k,v in enumerate(bulk_hist)),F(0))/(16*n*n)
                assert endpoint == -2*F(r['mean_source']['exact'])/(1-p)
                assert endpoint+bulk == F(r['source_derivative']['exact'])
                r['raw_derivative_endpoint_dilution']=fraction_entry(endpoint)
                r['raw_derivative_bulk_reconnection']=fraction_entry(bulk)
            rows.append({'observable':name,**r})
    return dict(geometry=[a,b],N=n,states=1<<n,topology_failures=top_fail,
        state_edge_pair_checks=checks,pair_update_failures=failed,
        no_preexisting_shared_component_failures=zero_shared_failed,
        nonzero_delta_by_shared_transition=dict(transition_hist),
        sign_flip_same_shared_count=same_s_sign,lattice_witness=witness,
        source_checksum=hashlib.sha256(json.dumps(sums).encode()).hexdigest(),
        covariance_controls=rows,elapsed_seconds=time.perf_counter()-start)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--geometry',nargs=2,type=int)
    p.add_argument('--join-audit',action='store_true');p.add_argument('--skip-pair-edge',action='store_true')
    args=p.parse_args()
    if args.out.exists():raise SystemExit('refusing to overwrite output')
    start=time.perf_counter()
    result=join_audit() if args.join_audit else enumerate_case(*(args.geometry or (3,0)),full_pair_edge_check=not args.skip_pair_edge)
    result['wall_seconds']=time.perf_counter()-start
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('covariance_controls','lattice_witness')},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
