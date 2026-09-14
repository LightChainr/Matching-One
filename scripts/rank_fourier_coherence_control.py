from __future__ import annotations
import math,json,argparse
from pathlib import Path
import mpmath as mp
SQ=((1,0),(-1,0),(0,1),(0,-1)); TRI=((1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1))

def ranks_all(L,dirs):
 N=L*L
 def idx(x,y):return (y%L)*L+(x%L)
 out=[0]*(1<<N)
 for mask in range(1<<N):
  seen=set();basis=[];done=False
  for root in range(N):
   if not(mask>>root&1) or root in seen:continue
   lift={root:(0,0)};st=[root];seen.add(root)
   while st and not done:
    u=st.pop();x=u%L;y=u//L;ux,uy=lift[u]
    for dx,dy in dirs:
     v=idx(x+dx,y+dy)
     if not(mask>>v&1):continue
     prop=(ux+dx,uy+dy)
     if v not in lift:lift[v]=prop;seen.add(v);st.append(v)
     else:
      dd=(prop[0]-lift[v][0],prop[1]-lift[v][1])
      if dd!=(0,0):
       assert dd[0]%L==0 and dd[1]%L==0
       w=(dd[0]//L,dd[1]//L)
       if not basis:basis=[w]
       elif basis[0][0]*w[1]-basis[0][1]*w[0]!=0:done=True;break
   if done:break
  out[mask]=2 if done else (1 if basis else 0)
 return out

def biased_transform(vals,p,N):
 a=[mp.mpf(v) for v in vals];q=1-p;s=mp.sqrt(p*q)
 for bit in range(N):
  step=1<<bit
  for base in range(0,1<<N,2*step):
   for j in range(step):
    i0=base+j;i1=i0+step;x0=a[i0];x1=a[i1]
    a[i0]=q*x0+p*x1
    a[i1]=s*(x1-x0)
 return a

def balance(ranks,N,pair):
 def mass(p,r):
  q=1-p
  return mp.fsum((1 if rr==r else 0)*p**m.bit_count()*q**(N-m.bit_count()) for m,rr in enumerate(ranks))
 return mp.findroot(lambda p:mass(p,2)-mass(p,0),pair)

def run(name,L,dirs,p):
 N=L*L;r=ranks_all(L,dirs);X=[x-1 for x in r];c=biased_transform(X,p,N)
 mean=c[0]; powers=[mp.mpf('0')]*(N+1); sums=[mp.mpf('0')]*(N+1)
 for mask,v in enumerate(c):
  n=mask.bit_count();powers[n]+=v*v;sums[n]+=v
 var=mp.fsum(powers[1:]);rows=[]
 for n in range(1,N+1):
  if powers[n] < mp.mpf('1e-45'):continue
  coh=sums[n]*sums[n]/(mp.binomial(N,n)*powers[n])
  sym=sums[n]*sums[n]/mp.binomial(N,n)
  rows.append({'n':n,'power_fraction':mp.nstr(powers[n]/var,30),'coherence_Gamma':mp.nstr(coh,30),
               'K_symmetric_variance_fraction':mp.nstr(sym/var,30),'coherent_sum':mp.nstr(sums[n],30)})
 return {'model':name,'L':L,'p':mp.nstr(p,45),'meanX':mp.nstr(mean,30),'varX':mp.nstr(var,30),
         'K_explained_fraction':mp.nstr(mp.fsum(sums[n]**2/mp.binomial(N,n) for n in range(1,N+1))/var,30),
         'levels':rows}

if __name__=='__main__':
 mp.mp.dps=60; ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 datasets=[]
 for L in (3,4):
  r=ranks_all(L,SQ);p=balance(r,L*L,(mp.mpf('.55'),mp.mpf('.63')));datasets.append(run('square',L,SQ,p))
 for L in (3,4):datasets.append(run('triangular',L,TRI,mp.mpf('.5')))
 out={'scope':'complete p-biased Fourier-level power/coherence controls for topological rank X; exact finite enumeration, high-precision arithmetic, no scaling fit','systems':datasets}
 a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
