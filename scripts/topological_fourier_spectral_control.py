import json
from pathlib import Path
import mpmath as mp

SQ=((1,0),(-1,0),(0,1),(0,-1))
TRI=((1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1))

def ranks_all(L, dirs):
    N=L*L
    def idx(x,y):return (y%L)*L+(x%L)
    out=[0]*(1<<N)
    for mask in range(1<<N):
        seen=set(); basis=[]
        def add(v):
            nonlocal basis
            if v==(0,0):return
            if not basis:basis=[v];return
            a=basis[0]
            if a[0]*v[1]-a[1]*v[0]!=0:basis=[a,v]
        done=False
        for root in range(N):
            if not(mask>>root&1) or root in seen:continue
            lift={root:(0,0)}; st=[root];seen.add(root)
            while st and not done:
                u=st.pop();x=u%L;y=u//L;ux,uy=lift[u]
                for dx,dy in dirs:
                    v=idx(x+dx,y+dy)
                    if not(mask>>v&1):continue
                    prop=(ux+dx,uy+dy)
                    if v not in lift:
                        lift[v]=prop;seen.add(v);st.append(v)
                    else:
                        dd=(prop[0]-lift[v][0],prop[1]-lift[v][1])
                        if dd!=(0,0):
                            assert dd[0]%L==0 and dd[1]%L==0
                            add((dd[0]//L,dd[1]//L))
                            if len(basis)==2:done=True;break
            if done:break
        out[mask]=len(basis)
    return out

def eval_rank_probs(ranks,p,N):
    P=[mp.mpf('0')]*3
    for m,r in enumerate(ranks):
        k=m.bit_count(); P[r]+=p**k*(1-p)**(N-k)
    return P

def balance_root(ranks,N):
    f=lambda p: eval_rank_probs(ranks,p,N)[2]-eval_rank_probs(ranks,p,N)[0]
    return mp.findroot(f,(mp.mpf('.55'),mp.mpf('.63')))

def spectral_moments(ranks,L,p):
    N=L*L; q=1-p
    P=eval_rank_probs(ranks,p,N)
    EX=P[2]-P[0]; EX2=P[0]+P[2]; var=EX2-EX*EX
    e1=mp.mpf('0'); bit0=1
    for m in range(1<<N):
        if m&bit0:continue
        k=m.bit_count(); d=ranks[m|bit0]-ranks[m]
        e1 += (d*d)*p**k*q**(N-1-k)
    total_pair=mp.mpf('0')
    for i in range(N):
      bi=1<<i
      for j in range(i+1,N):
        bj=1<<j; e=mp.mpf('0'); skip=bi|bj
        for m in range(1<<N):
            if m&skip:continue
            k=m.bit_count()
            dd=ranks[m|bi|bj]-ranks[m|bi]-ranks[m|bj]+ranks[m]
            e += (dd*dd)*p**k*q**(N-2-k)
        total_pair += 2*e
    mean_size = p*q*N*e1/var
    fac2=(p*q)**2*total_pair/var
    var_size=fac2+mean_size-mean_size**2
    return dict(p=p,P=P,varX=var,e_delta2=e1,mean_size=mean_size,fac2=fac2,var_size=var_size)

def walsh_level_energy(ranks,L):
    N=L*L; n=1<<N; a=[r-1 for r in ranks]; h=1
    while h<n:
      for i in range(0,n,2*h):
        for j in range(i,i+h):
          x=a[j];y=a[j+h];a[j]=x+y;a[j+h]=x-y
      h*=2
    levels=[0]*(N+1); max_even=0
    for mask,v in enumerate(a):
      levels[mask.bit_count()]+=v*v
      if mask and mask.bit_count()%2==0:max_even=max(max_even,abs(v))
    return levels,max_even

if __name__=='__main__':
    mp.mp.dps=50; out={}
    for name,dirs in [('square',SQ),('triangular',TRI)]:
      rows=[]
      for L in (3,4):
        ranks=ranks_all(L,dirs); N=L*L
        p=mp.mpf('.5') if name=='triangular' else balance_root(ranks,N)
        sp=spectral_moments(ranks,L,p)
        row={k:(mp.nstr(v,30) if isinstance(v,mp.mpf) else v) for k,v in sp.items() if k!='P'}
        row['P']=[mp.nstr(x,30) for x in sp['P']]
        if name=='triangular':
          lev,maxeven=walsh_level_energy(ranks,L)
          row['walsh_level_energy_integer_numerators']=lev
          row['max_abs_even_nonconstant_walsh_numerator']=maxeven
        rows.append(row)
      out[name]=rows
    target=Path('topological-fourier-spectral-control.json')
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
