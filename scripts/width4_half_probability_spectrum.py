#!/usr/bin/env python3
"""Exact observable spectral identity at p=1/2, certified on the rank automaton.

Integer sequences count configurations; s_m = 16**m M_(4,m)(1/2).
The identity is verified by a finite-dimensional annihilator certificate,
NOT extrapolated from a fitted finite list. No eigenvalue numerics required.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from collections import Counter
from lifted_boundary_rank import reachable_states, minimize_rank_machine

# Monic characteristic factors, coefficients in descending order.
FACTORS = [
    [-2, [1, -2]],
    [1, [1, 1]],
    [-1, [1, -15, 4]],
    [2, [1, -3, -4, 2]],
    [1, [1, -3, 2, -2]],
    [1, [1, -11, -9, 3, -10, -2]],
]
RECURRENCE = [1,-33,353,-1301,697,4327,-5785,-551,5534,
              -8814,872,2716,-3992,1480,80,-64]


def polynomial_product(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j] += x*y
    return out


def power_sums(coefficients, maximum):
    degree = len(coefficients)-1
    out = [degree]
    for n in range(1,maximum+1):
        if n <= degree:
            value = -sum(coefficients[i]*out[n-i] for i in range(1,n))
            value -= n*coefficients[n]
        else:
            value = -sum(coefficients[i]*out[n-i] for i in range(1,degree+1))
        out.append(value)
    return out


def integer_determinant(matrix):
    """Fraction-free Bareiss elimination, including row-pivot signs."""
    a = [list(row) for row in matrix]
    n, sign, previous = len(a), 1, 1
    for k in range(n-1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k+1,n) if a[i][k]), None)
            if pivot is None: return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value = pivot*a[i][j]-a[i][k]*a[k][j]
                if value % previous:
                    raise AssertionError('nonexact Bareiss division')
                a[i][j] = value//previous
            a[i][k] = 0
        previous = pivot
    return sign*a[-1][-1]


def squarefree_monic(coefficients):
    from fractions import Fraction
    def trim(a):
        while a and a[0]==0:a.pop(0)
        return a
    def remainder(a,b):
        a=list(a)
        while len(a)>=len(b):
            c=a[0]/b[0]
            for i,x in enumerate(b):a[i]-=c*x
            trim(a)
        return a
    a=list(map(Fraction,coefficients));n=len(a)-1
    b=[(n-i)*x for i,x in enumerate(a[:-1])]
    while b:a,b=b,remainder(a,b)
    return len(a)==1


def report():
    states, rows, initials = reachable_states(4)
    labels, reps, A, output, refinements = minimize_rank_machine(states,rows)
    b = Counter(labels[s] for s in initials)
    c = [r-1 for r in output]
    size = len(A)
    def times(v):
        return [sum(v[j] for j in row) for row in A]
    def read(v):
        return sum(weight*v[i] for i,weight in b.items())
    sequence = []
    v = c[:]
    for _ in range(2*15+1):
        sequence.append(read(v)); v = times(v)
    poly = [1]
    for _, factor in FACTORS:
        poly = polynomial_product(poly,factor)
    if not squarefree_monic(poly):
        raise AssertionError('observable annihilator is not squarefree')
    if poly != RECURRENCE:
        raise AssertionError('factor product != claimed recurrence')
    powers = [c]
    for _ in range(15): powers.append(times(powers[-1]))
    remainder = [sum(RECURRENCE[i]*powers[15-i][j] for i in range(16))
                 for j in range(size)]
    nonzeros = sum(x!=0 for x in remainder)
    # The polynomial need NOT annihilate c. Only its scalar observable image.
    # For a size-d matrix, d consecutive zero moments imply all via C-H.
    for k in range(size):
        if read(remainder):
            raise AssertionError(f'nonzero annihilator moment at {k}')
        remainder = times(remainder)
    moments = [(weight,power_sums(factor,31)) for weight,factor in FACTORS]
    for m in range(1,32):
        if sum(w*values[m] for w,values in moments) != sequence[m-1]:
            raise AssertionError('trace-power identity initial values mismatch')
    # Start at physical length m=2, not the convenient aliased m=1 term.
    hankel = [[sequence[i+j+1] for j in range(15)] for i in range(15)]
    det = integer_determinant(hankel)
    if not det: raise AssertionError('order-15 minimality not established')
    return {
        'schema':'matching-one.width4-half-probability-spectrum.v1',
        'scope':'M_(4,m)(1/2), all integer m>=2; no other p or width asserted',
        'integer_sequence_definition':'s_m=16^m*M_(4,m)(1/2)=b*A^(m-1)*c',
        'trace_formula':'s_m = sum_j weight_j * sum_{f_j(lambda)=0} lambda^m (multiplicity included)',
        'factors_descending_and_weights':FACTORS,
        'minimal_scalar_recurrence_order':15,
        'recurrence_descending':RECURRENCE,
        'recurrence_is_squarefree_over_Q':True,
        'hankel_determinant_15_start_m2':det,
        'finite_dimensional_certificate_dimension':size,
        'annihilator_zero_scalar_moments_verified':size,
        'pA_c_nonzero_entries':nonzeros,
        'initial_trace_power_identities_verified':31,
        'first_31_s_m':sequence,
        'normalization':'s_m counts configurations; divide by 16^m, not 16^(m-1)',
        'limitations':['one-parameter point only','not a global generator spectrum',
                       'not a continuum Jordan test','no all-width prefactor theorem'],
    }


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);args=ap.parse_args()
    data=report();text=json.dumps(data,indent=2,allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open('x',encoding='utf-8') as f:f.write(text)
    else:print(text,end='')
