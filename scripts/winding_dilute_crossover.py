#!/usr/bin/env python3
"""Dilute NN-site winding density: exact combinatorics and analytic controls.

The all-width theorem is in dilute-winding-crossover.md. No large-width
stationary density or Monte Carlo is calculated here. Small-width rational
functions below are retained inputs from the preceding intensity delivery.
Python stdlib for combinatorial tests; mpmath for report numerics.
"""
from __future__ import annotations
import argparse
from collections import Counter, deque
from fractions import Fraction
from itertools import combinations, product
from math import comb
from pathlib import Path
import json

SMALL_INTENSITIES = {'2': {'expression': 'p**2*(p - 1)**2*(p**2 + p + 1)/(p**2 - p + 1)', 'numerator_descending': [1, -1, 0, -1, 1, 0, 0], 'denominator_descending': [1, -1, 1], 'low_p_through_8': 'p**2 - p**4 - 2*p**5 + 2*p**7 + 2*p**8 + O(p**9)'}, '3': {'expression': '-p**3*(p - 1)**3*(p**6 + p**3 + 2*p**2 + 2*p + 1)/(p**6 - 3*p**5 + 3*p**4 + p**3 - p**2 - p + 1)', 'numerator_descending': [-1, 3, -3, 0, 1, 1, 0, -1, -1, 1, 0, 0, 0], 'denominator_descending': [1, -3, 3, 1, -1, -1, 1], 'low_p_through_8': 'p**3 - p**6 - 3*p**7 + O(p**9)'}, '4': {'expression': 'p**4*(p - 1)**4*(p**19 - 5*p**18 + 10*p**17 - 8*p**16 - 3*p**15 + 12*p**14 - 14*p**13 + 14*p**12 - 8*p**11 - 5*p**10 + 3*p**9 + 3*p**8 + 3*p**7 + 4*p**6 - 3*p**5 - 3*p**4 - 9*p**3 - 7*p**2 - 3*p - 1)/((p**2 - p - 1)*(p**2 - p + 1)*(p**15 - 7*p**14 + 21*p**13 - 33*p**12 + 25*p**11 - 2*p**10 - 8*p**9 + 3*p**8 - p**7 - 2*p**6 + 6*p**5 - 4*p**4 + 4*p**3 - 2*p**2 - p + 1))', 'numerator_descending': [1, -9, 36, -82, 110, -69, -38, 146, -199, 179, -95, 7, 21, -7, -10, 24, -28, 27, -20, 2, 5, -1, 1, -1, 0, 0, 0, 0], 'denominator_descending': [1, -9, 36, -82, 111, -78, 0, 50, -40, 5, 17, -21, 19, -12, 1, 5, -7, 3, 1, -1], 'low_p_through_8': 'p**4 + 4*p**6 - 8*p**7 + 7*p**8 + O(p**9)'}}


def neighbours(v: tuple[int, int], w: int):
    x,y=v
    return [((x+1)%w,y), ((x-1)%w,y), (x,y+1), (x,y-1)]


def winding_components(vertices: set[tuple[int,int]], w: int) -> int:
    """Independent physical lifted-coordinate BFS on a free-height cylinder."""
    seen: dict[tuple[int,int], tuple[int,int]]={}
    count=0
    for root in sorted(vertices):
        if root in seen:
            continue
        seen[root]=(0,0)
        todo=[root]; winds=False
        while todo:
            v=todo.pop(); x,y=v; vx,vy=seen[v]
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                t=((x+dx)%w,y+dy)
                if t not in vertices:
                    continue
                value=(vx+dx,vy+dy)
                if t not in seen:
                    seen[t]=value;todo.append(t)
                else:
                    error=(value[0]-seen[t][0],value[1]-seen[t][1])
                    if error[0]:
                        assert error[0]%w==0 and error[1]==0
                        winds=True
        count += int(winds)
    return count


def gap_subsets(w: int, k: int):
    """Subsets with cyclic distance at least three between selected columns."""
    if w<6 or k<0:
        raise ValueError('require w>=6 and k>=0')
    for c in combinations(range(w),k):
        if k<=1 or all((c[(i+1)%k]-c[i])%w>=3 for i in range(k)):
            yield c


def gap_count(w: int, k: int) -> int:
    if k==0:
        return 1
    if w<3*k:
        return 0
    return w*comb(w-2*k,k)//(w-2*k)


def separated_cycles(w: int, r: int):
    """Each positive winding cycle is addressed by its incoming height at x=0."""
    for cols in gap_subsets(w,2*r):
        for ups in combinations(cols,r):
            up=set(ups); events={x:(1 if x in up else -1) for x in cols}
            y=0; vset=set()
            for x in range(w):
                vset.add((x,y))
                if x in events:
                    y += events[x]
                    vset.add((x,y))
            assert y==0
            yield vset


def short_external_contacts(cycle: set[tuple[int,int]], w: int) -> int:
    """Off-cycle vertices touching >=2 cycle vertices (one-site returns)."""
    contacts=Counter()
    for v in cycle:
        for t in neighbours(v,w):
            if t not in cycle:
                contacts[t]+=1
    return sum(n>=2 for n in contacts.values())


def minimal_nonrow_census(max_width: int=8) -> list[dict]:
    """Only fixed occupation numbers in TWO free rows, not a full-size scan."""
    results=[]
    for w in range(3,max_width+1):
        points=[(x,y) for y in (0,1) for x in range(w)]
        tested=0; hits={}
        for k in (w,w+1,w+2):
            n=0
            for subset in combinations(points,k):
                tested+=1; s=set(subset)
                if any(all((x,y) in s for x in range(w)) for y in (0,1)):
                    continue
                n += int(winding_components(s,w)>0)
            hits[k]=n
        assert hits[w]==0 and hits[w+1]==0
        assert hits[w+2]==w*(w-3)
        results.append({'width':w,'tested_fixed_size_sets':tested,
                        'nonrow_winding_counts':hits,'predicted_first_count':w*(w-3)})
    return results


def cycle_geometry_checks(max_width: int=16) -> dict:
    totals={'cycles':0,'gap_subset_checks':0,'max_contact_over_r':0}
    for w in range(6,max_width+1):
        for r in range(0,min(2,w//6)+1):
            k=2*r
            subsets=list(gap_subsets(w,k))
            assert len(subsets)==gap_count(w,k)
            # The probability that k labelled iid columns are separated:
            # union bound over pairs at cyclic distance 0,1,2.
            from math import factorial
            assert Fraction(len(subsets)*factorial(k),w**k)>=1-Fraction(5*k*(k-1),2*w)
            totals['gap_subset_checks']+=1
            n=0
            for c in separated_cycles(w,r):
                assert len(c)==w+2*r
                assert winding_components(c,w)==1
                # An induced simple cycle: no occupied chords.
                assert all(sum(t in c for t in neighbours(v,w))==2 for v in c)
                z=short_external_contacts(c,w)
                assert z<=16*r
                if r:
                    totals['max_contact_over_r']=max(totals['max_contact_over_r'],z/r)
                n+=1
            assert n==gap_count(w,k)*comb(k,r)
            totals['cycles']+=n
    return totals


def evaluate_small_nu(w: int, p: Fraction) -> Fraction:
    s=SMALL_INTENSITIES[str(w)]
    def horner(a):
        v=Fraction(0)
        for x in a: v=v*p+x
        return v
    return horner(s['numerator_descending'])/horner(s['denominator_descending'])


def rational_series(w: int, degree: int=10) -> list[Fraction]:
    s=SMALL_INTENSITIES[str(w)]
    num=list(reversed(s['numerator_descending']))
    den=list(reversed(s['denominator_descending']))
    out=[]
    for n in range(degree+1):
        value=Fraction(num[n] if n<len(num) else 0)
        value -= sum(den[k]*out[n-k] for k in range(1,min(n,len(den)-1)+1))
        out.append(value/den[0])
    return out


def universal_lower_factor(w: int, p: Fraction) -> Fraction:
    if w<6 or not 0<p<=Fraction(1,8):
        raise ValueError('proved bound requires w>=6, 0<p<=1/8')
    return 1-128*w*p*p-64*w*(3*p)**(w-1)


def finite_cycle_lower(w: int, p: Fraction) -> Fraction:
    """Rigorous finite sum / p**w. May be weak, never a point estimate."""
    answer=Fraction(0)
    for r in range(w//6+1):
        n=w+2*r
        failure=16*r*p+12*n*p*p/(1-3*p)+16*n*(3*p)**(w-1)
        answer+=gap_count(w,2*r)*comb(2*r,r)*p**(2*r)*max(Fraction(0),1-failure)
    return answer


def central_trinomial(w: int) -> int:
    return sum(comb(w,2*k)*comb(2*k,k) for k in range(w//2+1))


def mp_number(x,mp):
    if isinstance(x,Fraction): return mp.mpf(x.numerator)/x.denominator
    return mp.mpf(x)


def walk_upper_normalized(w: int, p, mp):
    """Upper bound U_w/p^w; integral is diagnostic, not interval quadrature."""
    p=mp_number(p,mp)
    def value(theta):
        b=1-2*p*mp.cos(theta)
        t_over_p=2/(b+mp.sqrt(b*b-4*p*p))
        return mp.exp(w*mp.log(t_over_p))
    return mp.quad(value,[0,mp.pi/2,mp.pi])/mp.pi


def bessel_contrast(lam,mp):
    lam=mp.mpf(lam)
    return (mp.log(mp.besseli(0,2*lam))+mp.log(mp.besseli(0,6*lam))
            -2*mp.log(mp.besseli(0,4*lam)))/mp.log(mp.mpf(4)/3)


def matrix_loop_coefficients(max_n: int) -> list[Fraction]:
    """Exact 2^n L_n for the correlated two-state renewal control.

    det(I-A(2z,y))=1+(2+y+y^-1)*(-6z-5z²+2z³+z⁴)/32.
    Scale coefficients by 32^n to use integers throughout.
    """
    a={1:-6,2:-5,3:2,4:1}
    polys=[{}]
    ans=[Fraction(0)]
    for n in range(1,max_n+1):
        c={}
        if n<=4:
            base=-n*a[n]*32**(n-1)
            c={-1:base,0:2*base,1:base}
        for k in range(1,min(n-1,4)+1):
            scale=-a[k]*32**(k-1)
            for j,v in polys[n-k].items():
                for d,b in ((-1,1),(0,2),(1,1)):
                    c[j+d]=c.get(j+d,0)+scale*b*v
        c={k:v for k,v in c.items() if v}
        assert all(v>=0 for v in c.values())
        polys.append(c)
        ans.append(Fraction(c.get(0,0),32**n))
    return ans


def direct_matrix_loops(max_n: int=10) -> list[Fraction]:
    """Independent matrix-of-Laurent-polynomials trace expansion through n.
    Kernel below is the R-tilted kernel; return 2^w L_w.
    """
    P=((Fraction(3,4),Fraction(1,4)),(Fraction(1,4),Fraction(3,4)))
    A={}
    for i,j,x,dy in product(range(2),range(2),(1,2),(0,1)):
        y=(1 if j==0 else -1)*dy
        A[(i,j,x,y)]=P[i][j]/4
    power={(i,i,0,0):Fraction(1) for i in range(2)}
    result=[Fraction(0) for _ in range(max_n+1)]
    for k in range(1,max_n+1):
        nxt={}
        for (i,j,x,y),v in power.items():
            for (a,b,dx,dy),u in A.items():
                if j!=a or x+dx>max_n:continue
                key=(i,b,x+dx,y+dy)
                nxt[key]=nxt.get(key,Fraction(0))+v*u
        power=nxt
        for (i,j,x,y),v in power.items():
            if i==j and y==0:
                result[x]+=Fraction(x,k)*v
    return result


def numerical_report(dps: int=70):
    import mpmath as mp
    mp.mp.dps=dps
    fmt=lambda x:mp.nstr(x,26)
    # Exact retained densities, not newly computed large-width site data.
    small=[]
    for w in (2,3,4):
        for p in (Fraction(1,32),Fraction(1,64),Fraction(1,128)):
            n=evaluate_small_nu(w,p)
            lam=mp_number(p,mp)*w
            val=mp_number(n/p**w,mp)
            small.append({'width':w,'p':str(p),'nu_exact':str(n),
                'nu_over_p_to_w':fmt(val),'bessel':fmt(mp.besseli(0,2*lam)),
                'ratio':fmt(val/mp.besseli(0,2*lam))})
    curves=[]
    for lam in ('0.05','0.1','0.25','0.5','1','2','4','10','30'):
        curves.append({'lambda':lam,'beta_crossover':fmt(bessel_contrast(lam,mp))})
    bounds=[]
    for w,p in ((64,Fraction(1,256)),(256,Fraction(1,256)),
                (1024,Fraction(1,1024)),(1024,Fraction(1,256))):
        lam=mp_number(p,mp)*w
        b=mp.besseli(0,2*lam)
        lower=finite_cycle_lower(w,p)
        up=walk_upper_normalized(w,p,mp)
        uf=universal_lower_factor(w,p)
        assert mp_number(lower,mp)<=up
        assert up/b<=mp.exp(6*w*mp_number(p,mp)**2)
        assert mp_number(lower,mp)/b>=mp_number(uf,mp)
        bounds.append({'width':w,'p':str(p),'lambda':fmt(lam),
            'finite_lower_over_bessel':fmt(mp_number(lower,mp)/b),
            'walk_upper_over_bessel':fmt(up/b),
            'universal_lower_factor':fmt(mp_number(uf,mp)),
            'universal_upper_factor':fmt(mp.exp(6*w*mp_number(p,mp)**2))})
    coeff=matrix_loop_coefficients(384)
    matrix=[]
    D=mp.mpf(2)/3
    for w in (8,16,32,64,128):
        v=mp_number(coeff[w],mp)
        r=coeff[w]*coeff[3*w]/coeff[2*w]**2
        matrix.append({'width':w,'normalized_gaussian_amplitude':fmt(v*mp.sqrt(2*mp.pi*D*w)),
                       'beta_effective':fmt(mp.log(mp_number(r,mp))/mp.log(mp.mpf(4)/3))})
    # Different diffusion constants for the identical one-step distribution.
    c4,c8,c12=(central_trinomial(w) for w in (4,8,12))
    mr=Fraction(c4*c12,c8*c8)
    return {'small_width_retained_site_controls':small,'bessel_contrast_curve':curves,
       'matching_fixed_width_low_p_limit_not_finite_p_data':{
           'c4':c4,'c8':c8,'c12':c12,'R4_exact':str(mr),
           'beta_effective_limit':fmt(mp.log(mp_number(mr,mp))/mp.log(mp.mpf(4)/3))},
       'proved_bound_numerical_controls_not_density_estimates':bounds,
       'matrix_renewal_controls_not_site_model':{
          'mean_forward_length':'3/2','single_step_transverse_variance':'1/2',
          'asymptotic_variance_per_renewal':'1','D':'2/3',
          'naive_D_ignoring_correlations':'1/3','sequence':matrix},
       'log_nu_absolute_error_to_beta_error_multiplier':fmt(4/mp.log(mp.mpf(4)/3))}


def build_report(dps: int=70):
    geometry=cycle_geometry_checks()
    census=minimal_nonrow_census()
    a=matrix_loop_coefficients(10);b=direct_matrix_loops(10)
    assert a==b
    for w in (3,4):
        series=rational_series(w,10)
        assert series[w]==1 and series[w+1]==0
        assert series[w+2]==w*(w-3)
    return {'schema':'matching-one/dilute-prefactor/v1',
       'scope':'NN SITE dilute double limit plus an explicitly separate matrix-renewal control; no fixed-p OZ theorem',
       'geometry_checks':geometry,'minimal_nonrow_census':census,
       'independent_matrix_trace_coefficients':[str(v) for v in a[1:]],
       'numerics':numerical_report(dps)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--dps',type=int,default=70)
    args=parser.parse_args()
    if args.output.exists(): parser.error('refusing to overwrite an existing result')
    if args.dps<40: parser.error('use at least 40 decimal digits')
    result=build_report(args.dps)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(args.output)

if __name__=='__main__': main()
