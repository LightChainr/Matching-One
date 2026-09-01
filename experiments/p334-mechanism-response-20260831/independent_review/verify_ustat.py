#!/usr/bin/env python3
"""Independent exact audit of determinant and determinant-square U statistics.

No external dependencies, cloud calls, research imports, or data writes except
the verification report alongside this file. Fractions make all comparisons
exact; production implementations can use the displayed constant-size moments.
"""
from __future__ import annotations

import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path
import random


def falling(n, k):
    if n < k:
        raise ValueError(f'Need n >= {k}, got {n}')
    return math.prod(range(n-k+1, n+1))


def set_partitions(items):
    """Each set partition exactly once, independent of any research scorer."""
    if not items:
        yield ()
        return
    first, *rest = items
    for partition in set_partitions(rest):
        yield ((first,),) + partition
        for j, block in enumerate(partition):
            yield partition[:j] + ((first,) + block,) + partition[j+1:]


def mobius_distinct_sum(rows, columns):
    """Sum of products on ordered distinct sample indices in O(n), fixed k.

    positions are labelled even when columns repeat, e.g. (a,a,d,d).
    p_B=sum_i product_{r in B} X[i,columns[r]].
    coefficient of partition pi = product_B (-1)^(|B|-1)(|B|-1)!.
    """
    k = len(columns)
    positions = tuple(range(k))
    subsets = [s for size in range(1, k+1)
               for s in itertools.combinations(positions, size)]
    moments = {subset: sum((math.prod(row[columns[j]] for j in subset)
                            for row in rows), F(0)) for subset in subsets}
    total = F(0)
    for partition in set_partitions(list(positions)):
        coefficient = math.prod((-1)**(len(block)-1)*math.factorial(len(block)-1)
                                for block in partition)
        total += coefficient*math.prod(moments[tuple(sorted(block))] for block in partition)
    return total


def explicit_distinct_sum(rows, columns):
    return sum((math.prod(rows[i][j] for i,j in zip(indices, columns))
                for indices in itertools.permutations(range(len(rows)), len(columns))), F(0))


def repeated_22(rows, x, y):
    """Expanded O(n) numerator for (mu_x)^2 (mu_y)^2."""
    sx=sum((r[x] for r in rows), F(0));sy=sum((r[y] for r in rows), F(0))
    s20=sum((r[x]**2 for r in rows), F(0));s02=sum((r[y]**2 for r in rows), F(0))
    s11=sum((r[x]*r[y] for r in rows), F(0))
    s21=sum((r[x]**2*r[y] for r in rows), F(0))
    s12=sum((r[x]*r[y]**2 for r in rows), F(0))
    s22=sum((r[x]**2*r[y]**2 for r in rows), F(0))
    return (sx*sx*sy*sy-s20*sy*sy-s02*sx*sx-4*s11*sx*sy
            +s20*s02+2*s11*s11+4*s21*sy+4*s12*sx-6*s22)


def four_different(rows):
    """Independent explicit 15-partition expansion for a,b,c,d."""
    def s(*columns):
        return sum((math.prod(row[j] for j in columns) for row in rows), F(0))
    a,b,c,d=(s(j) for j in range(4))
    return (a*b*c*d
            -s(0,1)*c*d-s(0,2)*b*d-s(0,3)*b*c
            -s(1,2)*a*d-s(1,3)*a*c-s(2,3)*a*b
            +s(0,1)*s(2,3)+s(0,2)*s(1,3)+s(0,3)*s(1,2)
            +2*(s(0,1,2)*d+s(0,1,3)*c+s(0,2,3)*b+s(1,2,3)*a)
            -6*s(0,1,2,3))


def det_estimators(rows):
    n=len(rows)
    a,b,c,d=(sum((row[j] for row in rows), F(0)) for j in range(4))
    diagonal_ad=sum((r[0]*r[3] for r in rows), F(0))
    diagonal_bc=sum((r[1]*r[2] for r in rows), F(0))
    u2=(a*d-diagonal_ad-b*c+diagonal_bc)/falling(n,2)
    u4=(repeated_22(rows,0,3)+repeated_22(rows,1,2)-2*four_different(rows))/falling(n,4)
    return u2,u4


def h2(x,y):
    return (x[0]*y[3]+y[0]*x[3]-x[1]*y[2]-y[1]*x[2])/2


def unordered_pairing_estimator(rows):
    n=len(rows)
    total=F(0)
    for i,j,k,l in itertools.combinations(range(n),4):
        total += (h2(rows[i],rows[j])*h2(rows[k],rows[l])
                  +h2(rows[i],rows[k])*h2(rows[j],rows[l])
                  +h2(rows[i],rows[l])*h2(rows[j],rows[k]))/3
    return total/math.comb(n,4)


def finite_distribution_expectation(support, probabilities, n=4):
    mean=tuple(sum((p*x[j] for p,x in zip(probabilities,support)),F(0)) for j in range(4))
    target=mean[0]*mean[3]-mean[1]*mean[2]
    e2=e4=clipped=naive=F(0)
    for choices in itertools.product(range(len(support)),repeat=n):
        rows=[support[j] for j in choices]
        probability=math.prod(probabilities[j] for j in choices)
        u2,u4=det_estimators(rows)
        e2+=probability*u2;e4+=probability*u4
        clipped+=probability*max(F(0),u4)
        naive+=probability*u2*u2
    assert e2==target and e4==target*target
    return {'n':n,'support_size':len(support),'mean_matrix_entries':list(map(str,mean)),
            'target_det_mean':str(target),'target_det_mean_squared':str(target*target),
            'exact_E_U2':str(e2),'exact_E_U4':str(e4),
            'exact_E_clipped_U4':str(clipped),'exact_E_U2_squared':str(naive)}


def main():
    rng=random.Random(20260831)
    check_count=0
    for n in range(4,9):
        for case in range(5):
            rows=[tuple(F(rng.randrange(-5,6),rng.randrange(1,5)) for _ in range(4)) for _ in range(n)]
            if case==0:
                rows=[(r[0],r[1],r[1],r[3]) for r in rows]  # symmetric case included
            for columns in ((0,3),(1,2),(0,0,3,3),(1,1,2,2),(0,1,2,3)):
                assert mobius_distinct_sum(rows,columns)==explicit_distinct_sum(rows,columns)
                check_count+=1
            assert repeated_22(rows,0,3)==mobius_distinct_sum(rows,(0,0,3,3))
            assert repeated_22(rows,1,2)==mobius_distinct_sum(rows,(1,1,2,2))
            assert four_different(rows)==mobius_distinct_sum(rows,(0,1,2,3))
            u2,u4=det_estimators(rows)
            assert u4==unordered_pairing_estimator(rows)
            assert u2==sum((h2(rows[i],rows[j]) for i,j in itertools.combinations(range(n),2)),F(0))/math.comb(n,2)
    support=[tuple(map(F,x)) for x in ((1,0,0,2),(0,1,2,0),(2,-1,1,1))]
    nonnull=finite_distribution_expectation(support,[F(1,2),F(1,3),F(1,6)])
    null_support=[tuple(map(F,x)) for x in ((1,1,1,0),(1,1,1,2))]
    null=finite_distribution_expectation(null_support,[F(1,2),F(1,2)])
    balanced=null_support*2
    assert det_estimators(balanced)==(F(0),F(-1,3))
    assert null['exact_E_clipped_U4']=='1/8' and null['exact_E_U2_squared']=='1/4'
    large=null_support*36
    assert det_estimators(large)==(F(0),F(-1,71))
    assert falling(72,4)==72*71*70*69
    for n in (4,5,8,72):
        constant=[(F(2),F(3),F(4),F(5))]*n
        assert det_estimators(constant)==(F(-2),F(4))
    result={'status':'all_exact_fraction_assertions_passed','explicit_ordered_sum_checks':check_count,
            'random_rational_cases':25,'n_explicit_enumeration':[4,5,6,7,8],
            'independent_unordered_pairing_checks':25,'bell_partitions_order4':len(list(set_partitions([0,1,2,3]))),
            'finite_nonnull_distribution':nonnull,'finite_null_distribution':null,
            'negative_sample_n4':{'U2':'0','U4':'-1/3'},
            'negative_sample_n72':{'U2':'0','U4':'-1/71'},
            'n72_falling_denominator_order2':falling(72,2),'n72_falling_denominator_order4':falling(72,4),
            'arithmetic':'Exact fractions throughout; no tolerance comparisons, no data-dependent clipping',
            'scientific_boundary':'Conditional iid quartet matrices; arbitrary dependence between entries within a matrix is allowed. Old8 and new64 can pool only if they form disjoint conditionally iid quartets under the same prefix law.'}
    destination=Path(__file__).resolve().parent/'verification.json'
    destination.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
