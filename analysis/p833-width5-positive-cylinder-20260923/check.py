#!/usr/bin/env python3
"""Width-parametric adaptation of PR708 lifted_boundary_rank.py.
Algorithm reused, no independent-author claim. Per-transition cache removed.
"""
from collections import deque,Counter
from itertools import product
import json,time,resource
class DSU:
 def __init__(self,n,h=False): self.pa=list(range(n));self.g=[0]*n;self.h=h
 def find(self,v):
  p=self.pa[v]
  if p==v:return v,0
  r,d=self.find(p);self.g[v]+=d;self.pa[v]=r
  return r,self.g[v]
 def edge(self,u,v,g):
  ru,du=self.find(u);rv,dv=self.find(v)
  if ru==rv:self.h |= du+g!=dv
  else:self.pa[rv]=ru;self.g[rv]=du+g-dv

def row(d,off,w,mask):
 for x in range(w):
  y=(x+1)%w
  if (mask>>x&1) and (mask>>y&1):d.edge(off+x,off+y,int(x==w-1))
def canonical(d,keep,w):
 e=[(-1,0)]*(2*w);groups={}
 for old,new in keep:
  r,g=d.find(old);groups.setdefault(r,[]).append((new,g))
 for group in groups.values():
  r,g0=min(group)
  for v,g in group:e[v]=(r,0 if d.h else g-g0)
 if not d.h:
  mixed=[(j,g) for j,(r,g) in enumerate(e) if j>=w and 0<=r<w]
  if mixed:
   c=-min(mixed)[1]
   e=[(r,g+c*(int(j>=w)-int(r>=w))) if r>=0 else (r,g) for j,(r,g) in enumerate(e)]
 return (int(d.h),tuple(e))
def initial(w,mask):
 d=DSU(2*w);keep=[]
 for x in range(w):
  if mask>>x&1:d.edge(x,w+x,0);keep.extend(((x,x),(w+x,w+x)))
 row(d,w,w,mask);return canonical(d,keep,w)
def advance(s,mask):
 h,e=s;w=len(e)//2;d=DSU(3*w,h)
 for v,(r,g) in enumerate(e):
  if r>=0 and v!=r:d.edge(r,v,g)
 keep=[(x,x) for x in range(w) if e[x][0]>=0]
 for x in range(w):
  if mask>>x&1:
   keep.append((2*w+x,w+x))
   if e[w+x][0]>=0:d.edge(w+x,2*w+x,0)
 row(d,2*w,w,mask);return canonical(d,keep,w)
def close(s):
 h,e=s;w=len(e)//2;adj={v:[] for v,(r,g) in enumerate(e) if r>=0}
 def edge(u,v,x,y):adj[u].append((v,x,y));adj[v].append((u,-x,-y))
 for v,(r,g) in enumerate(e):
  if r>=0 and v!=r:edge(r,v,g,0)
 for x in range(w):
  if x in adj and w+x in adj:edge(w+x,x,0,1)
 first=(1,0) if h else None;pos={}
 for root in adj:
  if root in pos:continue
  pos[root]=(0,0);stack=[root]
  while stack:
   u=stack.pop();px,py=pos[u]
   for v,dx,dy in adj[u]:
    pp=(px+dx,py+dy)
    if v not in pos:pos[v]=pp;stack.append(v)
    else:
     a,b=pp[0]-pos[v][0],pp[1]-pos[v][1]
     if a or b:
      if first is None:first=(a,b)
      elif first[0]*b-first[1]*a:return 2
 return int(first is not None)
def labels(keys):
 idx={};lab=[]
 for k in keys:
  if k not in idx:idx[k]=len(idx)
  lab.append(idx[k])
 return lab

def build(w,cap=200000):
 t=time.perf_counter();ss=list(dict.fromkeys(initial(w,m) for m in range(1<<w)));idx={s:i for i,s in enumerate(ss)};ini=[idx[initial(w,m)] for m in range(1<<w)]
 ts=[];i=0
 while i<len(ss):
  tr=[]
  for m in range(1<<w):
   s=advance(ss[i],m)
   if s not in idx:
    if len(ss)>=cap:raise RuntimeError(('cap',w,cap))
    idx[s]=len(ss);ss.append(s)
   tr.append(idx[s])
  ts.append(tr);i+=1
  if i%2000==0:print('BFS',w,i,len(ss),round(time.perf_counter()-t,3),flush=True)
 out=list(map(close,ss));lab=out[:];ref=[3]
 while True:
  new=labels((lab[i],tuple(lab[j] for j in tr)) for i,tr in enumerate(ts));ref.append(len(set(new)))
  if len(set(new))==len(set(lab)):lab=new;break
  lab=new
 reps=[lab.index(i) for i in range(len(set(lab)))];q=[[lab[j] for j in ts[i]] for i in reps];o=[out[i] for i in reps];qi=[lab[i] for i in ini]
 for i,tr in enumerate(ts):assert [lab[j] for j in tr]==q[lab[i]] and out[i]==o[lab[i]]
 result={'w':w,'states':ss,'transitions':ts,'initial':ini,'rank_output':out,'quotient_map':lab,'quotient_representatives':reps,'quotient_transitions':q,'quotient_initial':qi,'quotient_rank_output':o,'refinement_counts':ref}
 print('DONE',w,len(ss),len(q),ref,'time',time.perf_counter()-t,'RSS',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,flush=True)
 return result

def histories(q,ini):
 a=len(q[0]);words={}
 for x,y in product(range(a),repeat=2):words.setdefault(q[ini[x]][y],[x,y])
 todo=deque(words)
 while todo:
  s=todo.popleft()
  for b,t in enumerate(q[s]):
   if t not in words:words[t]=words[s]+[b];todo.append(t)
 assert len(words)==len(q)
 return [words[i] for i in range(len(q))]
def state_for(word,q,ini):
 s=ini[word[0]]
 for b in word[1:]:s=q[s][b]
 return s

def dihedral(cert,full_check=True):
 w=cert['w'];q=cert['quotient_transitions'];ini=cert['quotient_initial'];o=cert['quotient_rank_output'];hs=histories(q,ini);actions=[]
 for sg,sh in product((1,-1),range(w)):
  masks=[sum(((m>>j)&1)<<((sg*j+sh)%w) for j in range(w)) for m in range(1<<w)]
  action=[state_for([masks[m] for m in h],q,ini) for h in hs]
  assert len(set(action))==len(q)
  if full_check:
   for s in range(len(q)):
    assert o[action[s]]==o[s]
    for m in range(1<<w):assert action[q[s][m]]==q[action[s]][masks[m]]
  actions.append(action)
 lab=labels(min(a[s] for a in actions) for s in range(len(q)));reps=[lab.index(i) for i in range(len(set(lab)))];weight=[Counter((m.bit_count(),lab[q[s][m]]) for m in range(1<<w)) for s in reps]
 for s in range(len(q)):assert Counter((m.bit_count(),lab[q[s][m]]) for m in range(1<<w))==weight[lab[s]]
 return {'w':w,'weighted':[[[k,j,v] for (k,j),v in sorted(row.items())] for row in weight],'output':[o[s] for s in reps],'histories':[hs[s] for s in reps],'dlabel':lab,'dreps':reps}


"""Independent frontier partition engine for vertical site survival on a cylinder.
No winding/gain state; marked top-connected component, NN or matching NN+NNN.
"""
import json
from collections import Counter
from itertools import product
import sympy as s
x,z=s.symbols('x z')
def nextstate(st,mask,matching=False):
 w=len(st);pa=list(range(2*w+1));top=2*w
 def find(i):
  while pa[i]!=i:pa[i]=pa[pa[i]];i=pa[i]
  return i
 def join(i,j):pa[find(j)]=find(i)
 for i,(r,t) in enumerate(st):
  if r>=0:join(i,r)
  if t:join(i,top)
 for i in range(w):
  if mask>>i&1:
   if mask>>((i+1)%w)&1:join(w+i,w+(i+1)%w)
   for dx in (-1,0,1) if matching else (0,):
    j=(i+dx)%w
    if st[j][0]>=0:join(w+i,j)
 groups={}
 for i in range(w):
  if mask>>i&1:groups.setdefault(find(w+i),[]).append(i)
 if find(top) not in groups:return None
 out=[(-1,0)]*w
 for r,ii in groups.items():
  for i in ii:out[i]=(min(ii),int(r==find(top)))
 return tuple(out)
def build_survival(w,matching=False):
 ini=tuple((0,1) for _ in range(w));ss=[ini];idx={ini:0};ts=[]
 for st in ss:
  row=[]
  for mask in range(1<<w):
   new=nextstate(st,mask,matching)
   if new is None:row.append(-1);continue
   if new not in idx:idx[new]=len(ss);ss.append(new)
   row.append(idx[new])
  ts.append(row)
 # exact all-p strong lumping of survival probabilities includes implicit death.
 lab=[0]*len(ss);sizes=[]
 while True:
  keyidx={};new=[]
  for i,row in enumerate(ts):
   c=Counter((m.bit_count(),lab[j] if j>=0 else -1) for m,j in enumerate(row));key=(lab[i],tuple(sorted(c.items())))
   if key not in keyidx:keyidx[key]=len(keyidx)
   new.append(keyidx[key])
  sizes.append(len(keyidx))
  if len(keyidx)==len(set(lab)):lab=new;break
  lab=new
 reps=[lab.index(i) for i in range(len(set(lab)))];A=s.zeros(len(reps))
 for i,rep in enumerate(reps):
  for m,j in enumerate(ts[rep]):
   if j>=0:A[i,lab[j]]+=x**m.bit_count()
 return A,len(ss),sizes

# Numerical arrays below contain exact bounded integers, never floating ranks.
import argparse, hashlib, warnings
from pathlib import Path
from fractions import Fraction as Fr
from math import comb
import numpy as np
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import QQ, ZZ, GF
from sympy.polys.matrices import DomainMatrix
RR,XX=ring('x',QQ)


def oracle_advance(state,mask):
    """Independent graph traversal of the boundary-star graph, not gain-DSU."""
    h,entries=state;w=len(entries)//2
    active={i for i,(r,g) in enumerate(entries) if r>=0}
    active.update(2*w+i for i in range(w) if mask>>i&1)
    adj={i:[] for i in active}
    def edge(u,v,g):adj[u].append((v,g));adj[v].append((u,-g))
    for v,(r,g) in enumerate(entries):
        if r>=0 and v!=r:edge(r,v,g)
    for i in range(w):
        if mask>>i&1:
            if entries[w+i][0]>=0:edge(w+i,2*w+i,0)
            j=(i+1)%w
            if mask>>j&1:edge(2*w+i,2*w+j,int(i==w-1))
    seen={};groups=[]
    for root in sorted(active):
        if root in seen:continue
        seen[root]=0;stack=[root];group=[]
        while stack:
            u=stack.pop()
            if u<w:group.append((u,seen[u]))
            elif u>=2*w:group.append((u-w,seen[u]))
            for v,g in adj[u]:
                proposed=seen[u]+g
                if v not in seen:seen[v]=proposed;stack.append(v)
                elif proposed!=seen[v]:h=True
        if group:groups.append(group)
    e=[(-1,0)]*(2*w)
    for group in groups:
        r,g0=min(group)
        for v,g in group:e[v]=(r,0 if h else g-g0)
    mixed=[(j,g) for j,(r,g) in enumerate(e) if j>=w and 0<=r<w]
    if mixed and not h:
        c=-min(mixed)[1]
        e=[(r,g+c*(int(j>=w)-int(r>=w))) if r>=0 else (r,g) for j,(r,g) in enumerate(e)]
    return (int(h),tuple(e))


def graph_rank(w,rows):
    """Physical torus lifted DFS; no compressed boundary, no cut gains."""
    height=len(rows);pos={};first=None
    for yy in range(height):
        for xx in range(w):
            root=yy*w+xx
            if not(rows[yy]>>xx&1) or root in pos:continue
            pos[root]=(0,0);stack=[root]
            while stack:
                u=stack.pop();px,py=pos[u];ux,uy=u%w,u//w
                for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                    vx,vy=(ux+dx)%w,(uy+dy)%height
                    if not(rows[vy]>>vx&1):continue
                    v=vy*w+vx;proposed=(px+dx,py+dy)
                    if v not in pos:pos[v]=proposed;stack.append(v)
                    else:
                        a,b=proposed[0]-pos[v][0],proposed[1]-pos[v][1]
                        assert a%w==0 and b%height==0
                        if a or b:
                            if first is None:first=(a,b)
                            elif first[0]*b-first[1]*a:return 2
    return int(first is not None)


def geometry_checks(c,d):
    w=c['w'];ss=c['states'];ts=c['transitions']
    for i,row in enumerate(ts):
        for mask,j in enumerate(row):assert oracle_advance(ss[i],mask)==ss[j]
    q=c['quotient_transitions'];ini=c['quotient_initial'];out=c['quotient_rank_output'];ncheck=0
    for h in (2,3):
        for word in product(range(1<<w),repeat=h):
            assert graph_rank(w,word)==out[state_for(word,q,ini)];ncheck+=1
    # All quotient states prepared at one common physical length.
    reach={q[ini[a]][b]:[a,b] for a,b in product(range(1<<w),repeat=2)}
    for h in range(2,10):
        if len(reach)==len(q):break
        nxt={}
        for i,word in reach.items():
            for mask,j in enumerate(q[i]):nxt.setdefault(j,word+[mask])
        reach=nxt
    assert len(reach)==len(q)
    for i,word in reach.items():assert graph_rank(w,word)==out[i];ncheck+=1
    return {'graph_transition_checks':len(ss)*(1<<w),'physical_torus_checks':ncheck,
            'common_preparation_length':h,'all_deterministic_states_preparable':len(q)}


def strong_counts(d):
    lab=[int(r==1) for r in d['output']];counts=[2]
    while True:
        sig=[]
        for i,row in enumerate(d['weighted']):
            c=Counter()
            for k,j,a in row:c[k,lab[j]]+=a
            sig.append((lab[i],tuple(sorted(c.items()))))
        new=labels(sig);counts.append(len(set(new)))
        if len(set(new))==len(set(lab)):break
        lab=new
    return counts


def coefficient_arrays(d):
    n=len(d['output']);a=np.zeros((d['w']+1,n,n),dtype=np.int64)
    for i,row in enumerate(d['weighted']):
        for k,j,v in row:a[k,i,j]+=v
    for k in range(d['w']+1):assert np.all(a[k].sum(axis=1)==comb(d['w'],k))
    return a


def positive_model(d,circuits):
    n=len(d['output']);A=coefficient_arrays(d);N=np.zeros((len(circuits),n),dtype=np.int64)
    for i,row in enumerate(circuits):
        for j,c in row:N[i,j]=c
    assert np.all(N.sum(axis=1)==0)
    for k in range(6):assert np.all(N@A[k]==0)
    for row in circuits:assert len(set(d['output'][i] for i,c in row))==1
    # Disjoint circuits make independence and the dimension-13 upper certificate immediate.
    used=set();V=[];urows={};newout=[];groups=[]
    for circuit in circuits:
        ids=[i for i,c in circuit];assert used.isdisjoint(ids);used.update(ids)
        base=len(V);rank=d['output'][ids[0]]
        if len(circuit)==3:
            a,c=[i for i,v in circuit if v==1];b=next(i for i,v in circuit if v==-2)
            terms=[[(a,1)],[(c,1)]]
            urows[a]=[(base,2)];urows[b]=[(base,1),(base+1,1)];urows[c]=[(base+1,2)]
            groups.append({'type':'midpoint','physical':[a,b,c],'latent':[base,base+1]})
        else:
            a,dd=[i for i,v in circuit if v==1];b,c=[i for i,v in circuit if v==-1]
            terms=[[(dd,1)],[(a,1),(b,-1),(c,1)],[(a,1),(b,1),(c,-1)]]
            urows[a]=[(base+1,1),(base+2,1)];urows[b]=[(base,1),(base+2,1)]
            urows[c]=[(base,1),(base+1,1)];urows[dd]=[(base,2)]
            groups.append({'type':'parallelogram','physical':[a,b,c,dd],'latent':[base,base+1,base+2]})
        for term in terms:
            v=np.zeros(n,dtype=np.int64)
            for i,c in term:v[i]=c
            V.append(v);newout.append(rank)
    for i in range(n):
        if i not in used:
            base=len(V);v=np.zeros(n,dtype=np.int64);v[i]=1;V.append(v)
            urows[i]=[(base,2)];newout.append(d['output'][i])
    V=np.array(V);m=len(V);U2=np.zeros((n,m),dtype=np.int64)
    for i,term in urows.items():
        for j,c in term:U2[i,j]=c
    assert m==385 and np.all(U2>=0) and np.all(U2.sum(axis=1)==2) and np.all(V.sum(axis=1)==1)
    assert np.array_equal(V@U2,2*np.eye(m,dtype=np.int64))
    for r in range(3):assert np.array_equal((np.array(d['output'])==r)[:,None]*U2,U2*(np.array(newout)==r))
    T2=[]
    for k in range(6):
        VA=V@A[k];assert np.all(VA>=0)  # coefficientwise: works for every rowwise p_t
        assert np.array_equal(U2@VA,2*A[k])
        T=VA@U2;assert np.all(T>=0)
        assert np.array_equal(U2@T,2*A[k]@U2)
        assert np.all(T.sum(axis=1)==2*comb(5,k));T2.append(T)
    sparse=[[[k,int(j),int(T2[k][i,j])] for k in range(6) for j in np.flatnonzero(T2[k][i])] for i in range(m)]
    return {'states':m,'coefficient_denominator':2,'coefficients':sparse,'outputs':newout,'groups':groups,
            'U2':[[[int(j),int(c)] for j,c in urows[i]] for i in range(n)],
            'V':[[[int(i),int(v[i])] for i in np.flatnonzero(v)] for v in V],
            'nullity_certificate':13,'all_coefficients_nonnegative':True,'all_p_intertwining_exact':True}


def word_minor(d,rec,ternary=False):
    prime=rec['prime'];xx=rec['odds'];assert sp.isprime(prime)
    n=len(d['output']);out=np.array(d['output']) if ternary else np.array([int(r==1) for r in d['output']])
    A=np.zeros((n,n),dtype=np.int64)
    for i,row in enumerate(d['weighted']):
        for k,j,c in row:A[i,j]+=c*pow(xx,k,prime)
    A%=prime;cache={}
    def col(word):
        word=tuple(map(int,word))
        if word not in cache:
            cache[word]=(out==word[0]).astype(np.int64) if len(word)==1 else (out==word[0])*(A@col(word[1:])%prime)
        return cache[word]
    piv=rec['pivots'];mat=np.array([col(word)[piv] for word in rec['words']],dtype=np.int64).T
    assert mat.shape==(rec['dim'],rec['dim'])
    b=mat.copy();det=1
    for k in range(len(b)):
        inds=np.flatnonzero(b[k:,k]);assert len(inds)>0
        i=k+int(inds[0])
        if i!=k:b[[k,i]]=b[[i,k]];det=-det
        v=int(b[k,k]);det=det*v%prime;inv=pow(v,-1,prime)
        for i in range(k+1,len(b)):
            if b[i,k]:b[i,k:]=(b[i,k:]-int(b[i,k])*inv%prime*b[k,k:])%prime
    assert det%prime==rec['det']
    dm=DomainMatrix.from_list([[int(v) for v in row] for row in mat],GF(prime))
    assert int(dm.det())%prime==det%prime
    return {'dimension_lower_bound':rec['dim'],'prime':prime,'odds_mod':xx,
            'minor_determinant_mod':det%prime,'maximum_word_length':max(map(len,rec['words'])),
            'second_determinant_path_agrees':True}


def find_cores(d):
    cores=[];absorb=None
    for r,word in ((0,[0,0]),(2,[(1<<d['w'])-1]*2)):
        start=d['histories'].index(word);seen={start};todo=[start];exits=set()
        for i in todo:
            for k,j,c in d['weighted'][i]:
                if d['output'][j]!=r:exits.add(j)
                elif j not in seen:seen.add(j);todo.append(j)
        assert len(exits)==1
        aa=next(iter(exits));assert d['output'][aa]==1 and all(j==aa for k,j,c in d['weighted'][aa])
        assert absorb is None or absorb==aa;absorb=aa
        ss=sorted(seen)
        for start in ss:
            reach={start};queue=[start]
            for i in queue:
                for k,j,c in d['weighted'][i]:
                    if j in seen and j not in reach:reach.add(j);queue.append(j)
            assert reach==seen
        cores.append(ss)
    return cores,absorb


def polynomial_coeff(poly):
    return [str(poly.get((k,),0)) for k in range(poly.degree()+1)] if poly else []

def interval_mul(a,b):
    terms=[u*v for u in a for v in b];return min(terms),max(terms)

def interval_poly(cs,iv):
    a=(Fr(0),Fr(0))
    for c in reversed(cs):
        a=interval_mul(a,iv);c=Fr(c);a=(a[0]+c,a[1]+c)
    return a

def interval_ratio(a,b):
    assert b[0]>0 or b[1]<0
    vals=[u/v for u in a for v in b];return min(vals),max(vals)


def spectral_analysis(d):
    w=d['w'];cores,absorb=find_cores(d);mats=[];chars=[]
    for ids in cores:
        ind={j:i for i,j in enumerate(ids)};A=sp.zeros(len(ids))
        for i,u in enumerate(ids):
            for k,j,c in d['weighted'][u]:
                if j in ind:A[i,ind[j]]+=c*x**k
        mats.append(A);chars.append(A.charpoly(z).as_expr())
    # Independent no-gain marked-frontier implementation, both colours.
    companion=[]
    for match in (False,True):
        A,n,sizes=build_survival(w,match)
        if match:A=A.applyfunc(lambda a:sp.expand(x**w*a.subs(x,1/x)))
        target=chars[0 if match else 1]
        assert sp.expand(A.charpoly(z).as_expr()-target)==0
        perm=[0,2,1,3,4] if w==5 and match else list(range(A.rows))
        native=mats[0 if match else 1]
        assert all(sp.expand(A[i,j]-native[perm[i],perm[j]])==0 for i in range(A.rows) for j in range(A.cols))
        companion.append({'matching':match,'raw_states':n,'strong_refinements':sizes,'characteristic_match_exact':True,'exact_permutation_to_native_core':perm})
    res=sp.resultant(*chars,z);scale,factors=sp.factor_list(res,x)
    fs=[{'coefficients':[int(c) for c in reversed(sp.Poly(f,x).all_coeffs())],'multiplicity':int(m)} for f,m in factors]
    fc=max(fs,key=lambda t:len(t['coefficients']))['coefficients'];F=sp.Poly.from_list(fc[::-1],x,domain=ZZ)
    assert F.degree()==(17 if w==4 else 51)
    intervals=[ab for ab,m in F.intervals(eps=sp.Rational(1,10**30)) if ab[0]>0]
    assert len(intervals)==(1 if w==4 else 4)
    # Primitive linear subresultant, followed by adjugate rows modulo F.
    G=sum(QQ(c)*XX**i for i,c in enumerate(fc));sub=sp.subresultants(*chars,z)
    linear=sp.Poly(next(a for a in reversed(sub) if sp.degree(a,z)==1),z)
    nn=-RR.from_expr(linear.nth(0));dd=RR.from_expr(linear.nth(1));g=nn.gcd(dd)
    NUM=nn.exquo(g)%G;DEN=dd.exquo(g)%G;assert NUM and DEN
    NP=[RR.one];DP=[RR.one]
    for k in range(7):NP.append((NP[-1]*NUM)%G);DP.append((DP[-1]*DEN)%G)
    vecs=[]
    for A in mats:
        n=A.rows;AA=[[RR.from_expr(a) for a in row] for row in A.tolist()];cp=[RR.from_expr(a) for a in A.charpoly(z).all_coeffs()]
        h=[RR.one]+[RR.zero]*(n-1);v=[RR.zero]*n
        for j in range(n):
            f=(NP[n-1-j]*DP[j])%G;v=[(a+f*b)%G for a,b in zip(v,h)]
            if j<n-1:
                h=[sum((h[i]*AA[i][k] for i in range(n)),RR.zero) for k in range(n)];h[0]+=cp[j+1]
        assert all(v) and sum(v,RR.zero)%G
        assert all((sum((v[i]*AA[i][j] for i in range(n)),RR.zero)*DEN-v[j]*NUM)%G==0 for j in range(n))
        vecs.append(v)
    # Only the crossing in x in (7/5,3/2) is the common Perron mode.
    assert F.count_roots(sp.Rational(7,5),sp.Rational(3,2))==1
    a,b=sp.refine_root(F,sp.Rational(7,5),sp.Rational(3,2),eps=sp.Rational(1,10**320))
    iv=(Fr(str(a)),Fr(str(b)));signs=[]
    for v in vecs:
        sg=[]
        for e in v:
            lo,hi=interval_poly(polynomial_coeff(e),iv)
            sg.append(1 if lo>0 else -1 if hi<0 else 0)
        assert len(set(sg))==1 and 0 not in sg;signs.append(sg)
    common=interval_ratio(interval_poly(polynomial_coeff(NUM),iv),interval_poly(polynomial_coeff(DEN),iv))
    assert common[0]>0
    pos=[]
    for aa,bb in intervals:
        mid=(aa+bb)/2;is_perron=bool(sp.Rational(7,5)<mid<sp.Rational(3,2))
        if not is_perron:
            ai,bi=sp.refine_root(F,aa,bb,eps=sp.Rational(1,10**80));vi=(Fr(str(ai)),Fr(str(bi)))
            lam=interval_ratio(interval_poly(polynomial_coeff(NUM),vi),interval_poly(polynomial_coeff(DEN),vi))
            lower=min(interval_poly(polynomial_coeff(RR.from_expr(sum(mats[0].row(i)))),vi)[0] for i in range(mats[0].rows))
            assert lam[1]<lower  # Collatz lower bound excludes a Perron crossing.
        pos.append({'odds_interval':[str(aa),str(bb)],'p_decimal':str(sp.N(mid/(1+mid),30)),'perron':is_perron})
    mid=(a+b)/2
    rho=sp.N((NUM.as_expr()/DEN.as_expr()).subs(x,mid)/(1+mid)**w,60)
    pexpr=sum(c*sp.Symbol('p')**k*(1-sp.Symbol('p'))**(F.degree()-k) for k,c in enumerate(fc))
    return {'w':w,'cores':cores,'absorb':absorb,'core_histories':[[d['histories'][i] for i in core] for core in cores],
            'core_matrices':[[[str(a) for a in row] for row in A.tolist()] for A in mats],
            'characteristic_polynomials':[str(a) for a in chars],'resultant_factors':fs,'resultant_constant':str(scale),
            'odds_factor':fc,'p_minimal_polynomial_coefficients':[int(c) for c in reversed(sp.Poly(pexpr,sp.Symbol('p')).all_coeffs())],
            'positive_roots':pos,'p_interval':[str(a/(1+a)),str(b/(1+b))],
            'p_decimal':str(sp.N(mid/(1+mid),60)),'rho_decimal':str(rho),
            'positive_left_eigenvectors_certified':signs,'companion_vertical_survival':companion,
            'binary_invisible_contrast_support':sum(map(len,cores)),
            'eigenvalue_nonzero':True,'rank0_rank2_preparations_separate_under_full_rank':True}


def irreducibility_check(factor,records):
    F=sp.Poly.from_list(factor[::-1],x);possible=set(range(1,F.degree()))
    for rec in records:
        prime=rec['prime'];assert sp.isprime(prime) and int(F.LC())%prime
        ff=sp.Poly(F.as_expr(),x,modulus=prime);assert ff.gcd(ff.diff()).degree()==0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore');_,fac=ff.factor_list()
        ds=sorted(int(f.degree()) for f,m in fac for _ in range(m));assert ds==rec['degrees']
        subsets={0}
        for a in ds:subsets|={v+a for v in subsets}
        possible &= subsets
    assert not possible


def square_control(c):
    w=c['w'];q=np.array(c['quotient_transitions']);ini=c['quotient_initial'];out=np.array(c['quotient_rank_output']);n=len(q)
    assert w*w<=25  # all exact counts < 2^63
    a=np.zeros((n,w+1),dtype=np.int64)
    for mask,st in enumerate(ini):a[st,mask.bit_count()]+=1
    for h in range(2,w+1):
        nxt=np.zeros((n,w*h+1),dtype=np.int64)
        for mask in range(1<<w):
            k=mask.bit_count();np.add.at(nxt[:,k:k+a.shape[1]],q[:,mask],a)
        a=nxt
    count=np.array([a[out==r].sum(axis=0) for r in range(3)]);assert np.array_equal(count.sum(axis=0),[comb(w*w,k) for k in range(w*w+1)])
    pp=sp.Symbol('p');M=sp.Poly(sum(int(v)*pp**k*(1-pp)**(w*w-k) for k,v in enumerate(count[2]-count[0])),pp)
    intervals=[ab for ab,m in M.intervals(eps=sp.Rational(1,10**60)) if 0<ab[0]<ab[1]<1];assert len(intervals)==1
    a,b=intervals[0]
    return {'rank_by_occupation':count.tolist(),'polynomial_coefficients':[int(c) for c in reversed(M.all_coeffs())],
            'root_interval':[str(a),str(b)],'root_decimal':str(sp.N((a+b)/2,55))}


def run(input_path,seed):
    raw=input_path.read_bytes();assert hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()==seed['input_blob']
    old=json.loads(raw);result={'schema':'matching-one.833.width-five.v1','input_blob':seed['input_blob'],'widths':{},'scope':'current output included; iid sites within each row; uniform row parameter may vary by row'}
    certs={};ds={}
    for w in (4,5):
        c=build(w);d=dihedral(c);certs[w]=c;ds[w]=d
        if w==4:
            for k in old:
                if k in c:assert json.loads(json.dumps(c[k]))==old[k],k
        geom=geometry_checks(c,d);ref=strong_counts(d)
        result['widths'][str(w)]={'geometric_states':len(c['states']),'deterministic_classes':len(c['quotient_transitions']),
            'deterministic_refinements':c['refinement_counts'],'dihedral_classes':len(d['output']),
            'binary_all_p_strong_refinements':ref,'geometry_checks':geom,
            'dihedral_transition_checks':2*w*len(c['quotient_transitions'])*(1<<w)}
        print('geometry passed',w,result['widths'][str(w)],flush=True)
    d=ds[5];model=positive_model(d,seed['kernel_rows'])
    print('coefficientwise positive compression passed',flush=True)
    minor=word_minor(d,seed['generic_minor']);assert minor['dimension_lower_bound']==385
    result['positive_model']=model;result['generic_dimension_certificate']=minor
    result['generic_minimal_positive_and_linear_dimension']=385
    result['kernel_circuits']=[{'terms':row,'histories':[[d['histories'][i],c] for i,c in row]} for row in seed['kernel_rows']]
    print('385-state positive model and exact generic lower bound passed',flush=True)
    result['spectral']={}
    for w in (4,5):
        result['spectral'][str(w)]=spectral_analysis(ds[w])
        print('spectral result passed',w,flush=True)
    rec=seed['collision_minor'];factor=result['spectral']['5']['odds_factor'];assert factor==rec['factor']
    irreducibility_check(factor,rec['irreducibility'])
    result['degree_51_irreducibility']=rec['irreducibility']
    collision={}
    for key,tern in (('binary',False),('ternary',True)):
        r=rec[key];assert sum(c*pow(r['odds'],k,r['prime']) for k,c in enumerate(factor))%r['prime']==0
        collision[key]=word_minor(d,r,ternary=tern)
        print('collision lower bound passed',key,flush=True)
    assert collision['binary']['dimension_lower_bound']==384 and collision['ternary']['dimension_lower_bound']==385
    result['degree_51_roots_observer_dimensions']={'binary':384,'full_rank':385,'lower_certificates':collision}
    result['square_controls']={str(w):square_control(certs[w]) for w in (4,5)}
    comparison=[];pc=sp.Rational(seed['reference_pc'])
    for w in (4,5):
        pb=sp.Rational(result['spectral'][str(w)]['p_decimal']);sq=sp.Rational(result['square_controls'][str(w)]['root_decimal']);j=sp.Rational(seed['jacobsen_table_2'][str(w)])
        assert abs(pb-j)<sp.Rational(1,10**40)
        comparison.append({'w':w,'blindpoint_decimal':result['spectral'][str(w)]['p_decimal'],'jacobsen_table_2_decimal':seed['jacobsen_table_2'][str(w)],
            'difference_from_published_table_bound':'<1e-40','reference_pc':seed['reference_pc'],
            'blindpoint_absolute_error_diagnostic':str(sp.N(abs(pb-pc),20)),
            'square_balance_root':result['square_controls'][str(w)]['root_decimal'],
            'square_absolute_error_diagnostic':str(sp.N(abs(sq-pc),20))})
    result['comparison']=comparison
    result['not_claimed']=['complete width-five exceptional-parameter census','new percolation threshold sequence',
        'all-width convergence proof','positive minimum at isolated width-five collision parameters',
        'full repository CI','independent external review','publication novelty']
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser(description='Exact #833 width-five analysis; does not access network or start wider runs.')
    ap.add_argument('input',type=Path);ap.add_argument('--out',type=Path,default=Path('result.json'))
    ap.add_argument('--model-out',type=Path,help='Optional full exact 385-state model coefficients.')
    args=ap.parse_args();seed=json.loads(Path(__file__).with_name('certificate.json').read_text());start=time.perf_counter()
    result=run(args.input,seed)
    model=result['positive_model']
    if args.model_out:args.model_out.write_text(json.dumps(model,separators=(',',':'))+'\n')
    model['weighted_transition_entries']=sum(map(len,model['coefficients']))
    del model['coefficients']  # regenerated by --model-out, not silently omitted from an asserted raw file
    args.out.write_text(json.dumps(result,separators=(',',':'),ensure_ascii=False)+'\n')
    print('PASS elapsed_seconds=',time.perf_counter()-start,'peak_RSS_KiB=',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,flush=True)
