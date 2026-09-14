from __future__ import annotations
import mpmath as mp, json
mp.mp.dps=60

SYSTEMS={
 'square':{
  3:[[1,9,36,78,90,45,0,0,0,0],[0,0,0,6,36,72,48,0,0,0],[0,0,0,0,0,9,36,36,9,1]],
  4:[[1,16,120,560,1812,4272,7448,9440,8082,3984,792,32,0,0,0,0,0],[0,0,0,0,8,96,560,1984,4580,6368,4704,1472,160,0,0,0,0],[0,0,0,0,0,0,0,16,208,1088,2512,2864,1660,560,120,16,1]]},
 'triangular':{
  3:[[1,9,36,75,45,0,0,0,0,0],[0,0,0,9,81,81,9,0,0,0],[0,0,0,0,0,45,75,36,9,1]],
  4:[[1,16,120,560,1808,4128,6304,5616,2160,304,0,0,0,0,0,0,0],[0,0,0,0,12,240,1704,5520,8550,5520,1704,240,12,0,0,0,0],[0,0,0,0,0,0,0,304,2160,5616,6304,4128,1808,560,120,16,1]]}}

def pd(counts,p):
    N=len(counts[0])-1;q=1-p
    P=[];D=[]
    for cs in counts:
        pp=mp.mpf('0');dd=mp.mpf('0')
        for k,c in enumerate(cs):
            if not c: continue
            w=mp.mpf(c)*p**k*q**(N-k)
            pp+=w
            dd+=w*(k/p-(N-k)/q)
        P.append(pp);D.append(dd)
    return P,D

def root(counts):
    N=len(counts[0])-1
    f=lambda p:sum(mp.mpf(counts[2][k]-counts[0][k])*p**k*(1-p)**(N-k) for k in range(N+1))
    return mp.findroot(f,(mp.mpf('.45'),mp.mpf('.65')))

def cond_mom(counts,p):
    N=len(counts[0])-1;q=1-p;out=[]
    for cs in counts:
        ww=[mp.mpf(c)*p**k*q**(N-k) for k,c in enumerate(cs)]
        Z=sum(ww); mu=sum(k*w for k,w in enumerate(ww))/Z
        var=sum((k-mu)**2*w for k,w in enumerate(ww))/Z
        out.append((Z,mu,var))
    return out

def speed(counts,p):
    P,D=pd(counts,p)
    return mp.sqrt(sum(dd*dd/pp for pp,dd in zip(P,D) if pp>0))

def fisher_length(counts):
    e=mp.mpf('1e-10')
    return mp.quad(lambda x:speed(counts,x),[e,.001,.01,.05,.2,.5,.8,.95,.99,.999,1-e])

def sphere_area(counts):
    def f(p):
        P,D=pd(counts,p); P0,P1,P2=P;D0,D1,D2=D
        bp=mp.mpf('.5')*(D0/P0-D2/P2)
        b=mp.mpf('.5')*mp.log(P0/P2)
        phip=-bp/(2*mp.cosh(b))
        return 4*mp.sqrt(P1)*phip
    e=mp.mpf('1e-10')
    return mp.quad(f,[e,.001,.01,.05,.2,.5,.8,.95,.99,.999,1-e])

def mutual_info(counts,p):
    N=len(counts[0])-1;q=1-p
    pk=[mp.binomial(N,k)*p**k*q**(N-k) for k in range(N+1)]
    pj=[];joint=[]
    for cs in counts:
        row=[mp.mpf(c)*p**k*q**(N-k) for k,c in enumerate(cs)]
        pj.append(sum(row));joint.append(row)
    I=mp.mpf('0')
    for j in range(3):
        for k,x in enumerate(joint[j]):
            if x:I+=x*mp.log(x/(pj[j]*pk[k]))
    return I

def local_geometry(counts,p):
    stats=cond_mom(counts,p); a=stats[0][0]; r=stats[1][0]
    mu=[x[1] for x in stats];var=[x[2] for x in stats]
    g1=mu[2]-mu[0];g2=var[2]-var[0]
    h1=mu[1]-(mu[0]+mu[2])/2;h2=var[1]-(var[0]+var[2])/2
    s=-2*h1/g1
    D2=4*(g1*h2-g2*h1)/g1**3
    kg=mp.sqrt(r/(2*a))*(D2-mp.mpf('.5')+mp.mpf('.5')*s*s)/(1+r*s*s)**mp.mpf('1.5')
    odd=a*g1*g1/2
    even=2*a*r*h1*h1
    full=(len(counts[0])-1)*p*(1-p)
    return dict(a=a,P1=r,g1=g1,h1=h1,Dprime=s,Dsecond=D2,geodesic_curvature=kg,
                fisher_fraction_odd=odd/full,fisher_fraction_even=even/full,
                fisher_fraction_total=(odd+even)/full)

res={'scope':'exact coefficient inputs L=3,4; high-precision quadrature for global Fisher invariants; no scaling fit','systems':{}}
for name,sys in SYSTEMS.items():
    res['systems'][name]={}
    for L,counts in sys.items():
        p=root(counts); loc=local_geometry(counts,p)
        row={'balance_root':mp.nstr(p,50),
             'local':{k:mp.nstr(v,40) for k,v in loc.items()},
             'fisher_length':mp.nstr(fisher_length(counts),40),
             'excess_over_endpoint_geodesic_pi':mp.nstr(fisher_length(counts)-mp.pi,40),
             'fisher_sphere_area_to_P1_zero_edge':mp.nstr(sphere_area(counts),40),
             'mutual_information_rank_K_nats':mp.nstr(mutual_info(counts,p),40),
             'sqrtL_times_mutual_information':mp.nstr(mp.sqrt(L)*mutual_info(counts,p),40)}
        res['systems'][name][str(L)]=row
print(json.dumps(res,indent=2))
