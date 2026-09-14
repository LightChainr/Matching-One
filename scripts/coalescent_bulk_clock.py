#!/usr/bin/env python3
"""Exact finite checks for circular cut genealogies and a shared-label noise clock.

The cycle enumeration is an exact finite combinatorial calculation, NOT a
new percolation simulation. The conditional Gaussian/Poisson formulas below
belong to the explicitly specified limiting process. No asymptotic claim is
inferred by fitting these controls.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, permutations
import json
from math import comb, factorial, log, sqrt
from pathlib import Path


def partition_key(blocks):
    return tuple(sorted(tuple(sorted(b)) for b in blocks))


def cycle_history(cycle, deletion_order):
    """Vertices are terminal cells; edges are the n distinct surviving cuts.

    n=2 has two parallel edges; deleting the last edge changes topology but
    not the one-block partition. Rotations are removed by fixing vertex 0.
    """
    n = len(cycle)
    blocks = [{j} for j in range(n)]
    out = [partition_key(blocks)]
    for edge in deletion_order:
        x, y = cycle[edge], cycle[(edge + 1) % n]
        i = next(i for i, b in enumerate(blocks) if x in b)
        j = next(i for i, b in enumerate(blocks) if y in b)
        if i != j:
            blocks[i] |= blocks[j]
            blocks.pop(j)
            out.append(partition_key(blocks))
    assert len(out) == n
    return tuple(out)


@lru_cache(None)
def enumerate_histories(n):
    if not 2 <= n <= 6:
        raise ValueError('exact history enumeration supports 2 <= n <= 6')
    histories = Counter()
    orders = tuple(permutations(range(n)))
    for tail in permutations(range(1, n)):
        cycle = (0,) + tail
        for order in orders:
            histories[cycle_history(cycle, order)] += 1
    total = factorial(n - 1) * factorial(n)
    assert sum(histories.values()) == total
    expected_histories = 1
    for k in range(2, n + 1):
        expected_histories *= comb(k, 2)
    assert len(histories) == expected_histories
    assert set(histories.values()) == {2 ** (n - 1)}
    levels = [Counter() for _ in range(n)]
    prefixes = {}
    for hist, count in histories.items():
        for depth, part in enumerate(hist):
            levels[depth][part] += count
        for depth in range(n - 1):
            prefix = hist[:depth + 1]
            prefixes.setdefault(prefix, Counter())[hist[depth + 1]] += count
    transition_checks = 0
    for prefix, future in prefixes.items():
        k = len(prefix[-1])
        assert len(future) == comb(k, 2)
        assert len(set(future.values())) == 1
        transition_checks += 1
    eppf_checks = 0
    mass_checks = 0
    for depth, law in enumerate(levels):
        k = n - depth
        second_mass = F(0)
        for part, count in law.items():
            prob = F(factorial(n-k) * factorial(k) * factorial(k-1),
                     factorial(n) * factorial(n-1))
            for block in part:
                prob *= factorial(len(block))
            assert F(count, total) == prob
            eppf_checks += 1
            second_mass += F(count, total) * F(
                sum(len(block)*(len(block)+1) for block in part), n*(n+1))
        assert second_mass == F(2, k + 1)
        mass_checks += 1
    return {'n': n, 'cyclic_orders': factorial(n-1),
            'edge_orders': factorial(n), 'realizations': total,
            'ranked_histories': len(histories),
            'multiplicity_per_history': 2**(n-1),
            'full_history_transition_checks': transition_checks,
            'eppf_equalities': eppf_checks,
            'dirichlet_second_mass_equalities': mass_checks}


def pair_separated_by_edges(n, survival):
    """Two uniformly chosen terminal cell labels; average cyclic separation."""
    survival = F(survival)
    if not 0 <= survival <= 1 or n < 2:
        raise ValueError('bad parameters')
    gone = 1 - survival
    return sum((1-gone**j)*(1-gone**(n-j)) for j in range(1,n))/F(n-1)


def pair_separated_by_levels(n, survival):
    survival = F(survival)
    ans = F(0)
    for k in range(2,n+1):
        prob = comb(n,k) * survival**k * (1-survival)**(n-k)
        same = F(2*(n-k), (k+1)*(n-1))
        ans += prob * (1-same)
    return ans


# Small bivariate rational jets. Keys are powers of the two Fourier variables.
DEGREE = 4


def add(*xs):
    ans = {}
    for x in xs:
        for m,v in x.items():
            ans[m] = ans.get(m,F(0)) + v
    return {m:v for m,v in ans.items() if v}


def scale(x,a):
    return {m:v*F(a) for m,v in x.items() if v*F(a)}


def mul(x,y):
    ans = {}
    for (i,j),a in x.items():
        for (k,l),b in y.items():
            if i+j+k+l <= DEGREE:
                key=(i+k,j+l)
                ans[key]=ans.get(key,F(0))+a*b
    return {m:v for m,v in ans.items() if v}


def inv(x):
    c=x.get((0,0),F(0))
    if not c:
        raise ZeroDivisionError('zero jet constant')
    tail=scale(add(x,{(0,0):-c}),-1/c)
    ans={(0,0):F(1)}; power=ans.copy()
    for _ in range(DEGREE):
        power=mul(power,tail); ans=add(ans,power)
    return scale(ans,1/c)


def score_cf_jet(k,a):
    k,a=F(k),F(a)
    A={(0,0):F(1),(2,0):F(1,2)}
    numerator=add(A,{(0,0):k-1})
    denominator2=add(numerator,{(0,2):k/2,(1,1):k*a})
    one=mul(numerator,mul(inv(A),inv(denominator2)))
    return mul(one,one)


def score_moments(k,a):
    k,a=F(k),F(a)
    if k < 1 or a*a*k > 1:
        raise ValueError('conditional Gaussian covariance is not admissible')
    jet=score_cf_jet(k,a)
    cov=-jet.get((1,1),F(0))
    var1=-2*jet.get((2,0),F(0)); var2=-2*jet.get((0,2),F(0))
    fourth1=24*jet.get((4,0),F(0)); fourth2=24*jet.get((0,4),F(0))
    mixed=4*jet.get((2,2),F(0))
    assert cov==2*a and var1==var2==2
    assert fourth1==fourth2==18
    assert mixed==4+F(2,k)+12*a*a
    return {'k':k,'a':a,'variance':var1,'covariance':cov,
            'fourth':fourth1,'mixed_fourth':mixed,
            'raw_score_correlation':a,
            'squared_score_correlation':(F(1,k)+6*a*a)/7}


def six_variable_transform(k,a,r,s,t,z,v,u1,u2):
    """Laplace(Y1,Y2), PGF(H1,H2), Fourier(Q1,Q2), all scalar inputs.

    Formula is real even when a cross term makes B negative. The complete
    quadratic form in its denominator is positive semidefinite for k*a^2<=1.
    """
    k,a,r,s,t,z,v,u1,u2=map(F,(k,a,r,s,t,z,v,u1,u2))
    if k<1 or a*a*k>1 or min(r,s,t)<0 or not (0<=z<=1 and 0<=v<=1):
        raise ValueError('parameters outside stated transform domain')
    A=1+s+r*(1-z)+u1*u1/2
    B=k*(t+u2*u2/2+a*u1*u2)+r*(1-v)*(z+k-1)
    assert A>0 and A+k-1+B>0
    return ((A+k-1)/(A*(A+k-1+B)))**2


def retention_mean(k):
    k=float(k)
    if k<1:
        raise ValueError('k must be at least one')
    if k==1:
        return 1.0
    return 2*((k-1)-log(k))/(k-1)**2


def retention_density(x,k):
    a=float(k)-1
    if not 0<x<1 or a<=0:
        return 0.0
    return 2*a*x*(3*k-1-a*x)/(1+a*x)**4


def retention_moment_mixture(k, power, terms=10000):
    """Independent NB/Beta mixture; floating control, not an interval bound."""
    k=float(k); prob=1/k**2; ans=prob; beta_moment=1.0
    for n in range(1,terms):
        prob *= (n+1)/n * (k-1)/k
        beta_moment *= (n+1)/(n+1+power)
        ans += prob*beta_moment
        if prob<1e-17 and n>50:
            break
    return ans


def simpson(f,n=20000):
    if n%2:
        raise ValueError('even n required')
    total=f(0)+f(1)
    for j in range(1,n):
        total+=(4 if j%2 else 2)*f(j/n)
    return total/(3*n)


def label_covariance_checks():
    ans=[]
    for e in (F(1,100),F(1,1000),F(1,10000)):
        for t1,t2 in ((F(1),F(2)),(F(1),F(4)),(F(2),F(5))):
            u,v=e*t1,e*t2
            outcomes=((0,0,1-v),(0,1,v-u),(1,1,u))
            covariance=sum(prob*(x-u)*(y-v) for x,y,prob in outcomes)
            assert covariance==u*(1-v)
            # Covariance of the rescaled empirical sheet per unit length.
            assert covariance/e==t1*(1-e*t2)
            ans.append({'epsilon':e,'t1':t1,'t2':t2,
                        'scaled_covariance':covariance/e})
    return ans


def gaussian_even_moment(m):
    if m % 2:
        return 0
    ans = 1
    for j in range(1,m,2):
        ans *= j
    return ans


def xy_moment(a,b,h):
    """E[X^a Y^b (X+Y)^h] for independent rate-one exponentials."""
    return sum(comb(h,j)*factorial(a+j)*factorial(b+h-j) for j in range(h+1))


def closure_stationarity_checks():
    """Exact invariant-generator identities for 4-coordinate closure.

    Test functions: x_minus^a*x_plus^b*Q^m*(H falling n).
    eta is the temporal label-refresh coefficient; dilute models have eta=1/w.
    All expectations below use the CONDITIONAL Gaussian and Poisson marks.
    """
    checked=0
    for eta in (F(0),F(1,8),F(1,4),F(1)):
      for a in range(4):
       for b in range(4):
        for m in (0,2,4):
         for n in range(4):
          r=F(3,7)
          d=m//2; h=d+n
          gaussian=gaussian_even_moment(m)
          common=gaussian*r**n
          mean=F(xy_moment(a,b,h))*common
          drift=(a+b)*mean + F(m,2)*(1-eta)*mean
          if m:
              diffusion=F(m*(m-1),2)*eta*gaussian_even_moment(m-2)*r**n*xy_moment(a,b,h)
          else:
              diffusion=F(0)
          immigration=n*mean
          cuts=F(0)
          for j in range(h+1):
              coef=comb(h,j)*common
              left=factorial(a+j+1)*factorial(b+h-j)
              right=factorial(a+j)*factorial(b+h-j+1)
              cuts += coef*(F(left,a+j+1)-left+F(right,b+h-j+1)-right)
          assert drift+diffusion+immigration+cuts==0
          checked+=1
    return checked


def report():
    genealogy=[enumerate_histories(n) for n in range(2,7)]
    pair_checks=0
    for n in range(2,13):
        for p in (F(0),F(1,7),F(1,3),F(1,2),F(4,5),F(1)):
            assert pair_separated_by_edges(n,p)==pair_separated_by_levels(n,p)
            pair_checks+=1
    scores=[]
    for b in (4,8):
        for root in (1,2,3):
            ratio=F(root*root)
            k=ratio**b
            a=F(1,root**(b+1))
            d=score_moments(k,a)
            d.update({'barrier_order':b,'parameter_ratio':ratio,
                      'tagged_length_correlation':1/k,
                      'self_normalized_squared_correlation':float(1/ratio)*retention_mean(k)})
            scores.append(d)
    retention=[]
    for k in (2,3,4,8):
        integral=simpson(lambda x: retention_density(x,k)) + 1/k**2
        first=simpson(lambda x: x*retention_density(x,k)) + 1/k**2
        # Density is continuous to the endpoint; Simpson uses its limiting values.
        # Our density function excludes endpoints; correct their quadrature weights.
        endpoint=2*(k-1)*(2*k)/(k**4)
        integral += endpoint/(3*20000)
        first += endpoint/(3*20000)
        mixture=retention_moment_mixture(k,1)
        assert abs(integral-1)<1e-10
        assert abs(first-retention_mean(k))<1e-10
        assert abs(mixture-retention_mean(k))<1e-12
        retention.append({'clock_ratio':k,'no_split_atom':1/k**2,
                          'mean_retention':retention_mean(k),
                          'NB_Beta_mean_check':mixture,
                          'quadrature_mass':integral})
    return {'scope':'Exact finite cycle combinatorics and limiting-process identities; no new site simulation.',
            'genealogy':genealogy,'pair_survival_equalities':pair_checks,
            'nested_label_covariances':label_covariance_checks(),
            'four_coordinate_generator_equalities':closure_stationarity_checks(),
            'score_process_moments':scores,'retention_law_controls':retention,
            'conditioning':'Terminal cells uniformly relabelled; cyclic order forgotten. Not fixed spatial tags.',
            'literature':'Bertoin--Goldschmidt math/0408128v1, Proposition 1 and section 2.3.'}


def encode(x):
    if isinstance(x,F):
        return {'fraction':str(x),'decimal':float(x)}
    raise TypeError(type(x).__name__)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path)
    args=p.parse_args()
    text=json.dumps(report(),ensure_ascii=False,sort_keys=True,indent=2,default=encode)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text,encoding='utf-8')
    else:
        print(text,end='')

if __name__=='__main__':
    main()
