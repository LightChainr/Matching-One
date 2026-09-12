#!/usr/bin/env python3
"""Exact all-fugacity sector identities from PR #708's rank automaton.

Standard library only. No Monte Carlo, numerical rank, or interpolation guess.
The large-radix evaluation is injective by the coefficient bound checked below.
All physical claims concern width four and lengths m >= 2.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from math import factorial
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results/research-control-20260912'
DEFINITION = RESULTS / 'width4-parametric-definition.json'
SOURCE = RESULTS / 'width4-rank-closure-certificate.json'


def poly_add(a, b):
    c = [0] * max(len(a), len(b))
    for i, v in enumerate(a): c[i] += v
    for i, v in enumerate(b): c[i] += v
    return c


def poly_mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, v in enumerate(a):
        for j, w in enumerate(b): c[i+j] += v*w
    return c


def polynomial_value(a, t):
    value = 0
    for c in reversed(a): value = value*t+c
    return value


def x_poly_product(a, b):
    """Descending x coefficients, each an ascending integer t polynomial."""
    out = [[0] for _ in range(len(a)+len(b)-1)]
    for i, p in enumerate(a):
        for j, q in enumerate(b): out[i+j] = poly_add(out[i+j], poly_mul(p,q))
    return out


def packed_multiply(value, coefficients, bits):
    return sum(c*(value << (bits*k)) for k,c in enumerate(coefficients) if c)


def trace_sequence(coefficients, maximum, *, t=None, bits=None):
    """Newton identities; exactly one of integer t or radix bit size is used."""
    degree = len(coefficients)-1
    if (t is None) == (bits is None): raise ValueError('supply t or bits')
    co = [polynomial_value(a,t) for a in coefficients] if bits is None else None
    def multiply(v,i):
        return v*co[i] if bits is None else packed_multiply(v,coefficients[i],bits)
    result = [degree]
    for m in range(1,maximum+1):
        if m <= degree:
            result.append(-sum(multiply(result[m-i],i) for i in range(1,m))
                          -m*multiply(1,m))
        else:
            result.append(-sum(multiply(result[m-i],i) for i in range(1,degree+1)))
    return result


def weighted_quotient(source):
    """Common strong lumping for every coefficient A_k in A(t)=sum t^k A_k."""
    rows=source['quotient_transitions']; output=source['quotient_rank_output']
    labels=output[:]; counts=[len(set(labels))]
    while True:
        lookup={}; new=[]
        for i,row in enumerate(rows):
            totals=Counter((mask.bit_count(),labels[j]) for mask,j in enumerate(row))
            key=(labels[i],tuple(sorted(totals.items())))
            new.append(lookup.setdefault(key,len(lookup)))
        counts.append(len(lookup))
        if len(lookup)==len(set(labels)):
            labels=new; break
        labels=new
    reps=[labels.index(i) for i in range(len(set(labels)))]
    weighted=[]
    for i in reps:
        counts_row=Counter((mask.bit_count(),labels[j]) for mask,j in enumerate(rows[i]))
        weighted.append([(k,j,n) for (k,j),n in sorted(counts_row.items())])
    # Check ALL source rows against their claimed lumped transitions.
    for i,row in enumerate(rows):
        expected=Counter((mask.bit_count(),labels[j]) for mask,j in enumerate(row))
        actual={(k,j):n for k,j,n in weighted[labels[i]]}
        if dict(expected)!=actual: raise AssertionError('not a common lumping')
    return dict(labels=labels,reps=reps,rows=weighted,
                output=[output[i] for i in reps],
                initial=[labels[i] for i in source['quotient_initial']],counts=counts)


def block_rows(source, ids):
    pos={v:i for i,v in enumerate(ids)}; rows=[]
    for i in ids:
        counts=Counter((mask.bit_count(),pos[j]) for mask,j in
                       enumerate(source['quotient_transitions'][i]) if j in pos)
        rows.append([(k,j,n) for (k,j),n in sorted(counts.items())])
    return rows


def strongly_connected(rows):
    n=len(rows)
    for reverse in (False,True):
        edges=[set() for _ in rows]
        for i,row in enumerate(rows):
            for _,j,c in row:
                if c: edges[j if reverse else i].add(i if reverse else j)
        visited={0}; todo=[0]
        while todo:
            for j in edges[todo.pop()]:
                if j not in visited: visited.add(j);todo.append(j)
        if len(visited)!=n:return False
    return True


def direct_matrix_traces(rows, maximum, bits):
    """Independent multiplication of the SMALL raw block, not Newton recurrence."""
    n=len(rows); powers=[[int(i==j) for j in range(n)] for i in range(n)]
    traces=[n]
    for _ in range(maximum):
        powers=[[sum(c*(powers[j][column]<<(bits*k)) for k,j,c in row)
                 for column in range(n)] for row in rows]
        traces.append(sum(powers[i][i] for i in range(n)))
    return traces


def perron_quotient(rows):
    """Positive equitable block quotient; the Perron eigenvalue is retained."""
    labels=[0]*len(rows)
    while True:
        lookup={}; new=[]
        for i,row in enumerate(rows):
            totals=Counter()
            for k,j,n in row: totals[k,labels[j]]+=n
            key=(labels[i],tuple(sorted(totals.items())))
            new.append(lookup.setdefault(key,len(lookup)))
        if len(lookup)==len(set(labels)):labels=new;break
        labels=new
    reps=[labels.index(i) for i in range(len(set(labels)))];out=[]
    for i in reps:
        totals=Counter()
        for k,j,n in rows[i]:totals[k,labels[j]]+=n
        out.append([(k,j,n) for (k,j),n in sorted(totals.items())])
    for i,row in enumerate(rows):
        totals=Counter()
        for k,j,n in row:totals[k,labels[j]]+=n
        if dict(totals)!={(k,j):n for k,j,n in out[labels[i]]}:
            raise AssertionError('equitable block identity failed')
    return {'labels':labels,'rows':out}


def scalar_sequences(quotient, maximum, t):
    v=[0]*len(quotient['rows'])
    for mask,j in enumerate(quotient['initial']):v[j]+=t**mask.bit_count()
    p0=[];p2=[]
    for _ in range(maximum):
        p0.append(sum(vv for vv,r in zip(v,quotient['output']) if r==0))
        p2.append(sum(vv for vv,r in zip(v,quotient['output']) if r==2))
        new=[0]*len(v)
        for i,row in enumerate(quotient['rows']):
            for k,j,c in row:new[j]+=c*t**k*v[i]
        v=new
    return p0,p2,[b-a for a,b in zip(p0,p2)]


def bareiss(matrix):
    a=[list(row) for row in matrix];n=len(a);previous=1;sign=1
    for k in range(n-1):
        if not a[k][k]:
            j=next((j for j in range(k+1,n) if a[j][k]),None)
            if j is None:return 0
            a[k],a[j]=a[j],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=pivot*a[i][j]-a[i][k]*a[k][j]
                if numerator%previous:raise AssertionError('Bareiss division not exact')
                a[i][j]=numerator//previous
            a[i][k]=0
        previous=pivot
    return sign*a[-1][-1]


def squarefree(coefficients):
    a=list(map(Fraction,coefficients));n=len(a)-1
    b=[(n-i)*a[i] for i in range(n)]
    def remainder(a,b):
        a=a[:]
        while len(a)>=len(b):
            c=a[0]/b[0]
            for i,x in enumerate(b):a[i]-=c*x
            while a and not a[0]:a.pop(0)
        return a
    while b:a,b=b,remainder(a,b)
    return len(a)==1


def factor_weights(definition, which):
    terms=Counter()
    def add_block(size,sgn):
        for name,multiplicity in definition['blocks'][str(size)]['factorization']:
            terms[name]+=sgn*multiplicity
    if which in ('P0','M'):add_block(5,1 if which=='P0' else -1)
    if which in ('P2','M'):add_block(15,1)
    if which!='M':add_block(16,-1);terms['positive_square']+=2
    return {name:w for name,w in terms.items() if w}


def verify(source_path=SOURCE,definition_path=DEFINITION):
    started=time.perf_counter();raw=Path(source_path).read_bytes();source=json.loads(raw)
    definition=json.loads(Path(definition_path).read_text())
    blob=hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    if blob!=definition['source_blob']:raise ValueError('wrong upstream certificate bytes')
    q=weighted_quotient(source);dimension=len(q['rows'])
    maximum=dimension+15+16+1
    bits=4*maximum+8;radix=1<<bits
    if 34*(16**maximum)>=radix:raise AssertionError('insufficient radix for injectivity')
    blocks={};traces={}
    for size_text,entry in definition['blocks'].items():
        n=int(size_text); rows=block_rows(source,entry['quotient_state_indices'])
        if not strongly_connected(rows):raise AssertionError('claimed block not irreducible')
        if not any(j==i for i,row in enumerate(rows) for _,j,_ in row):
            raise AssertionError('aperiodicity loop missing')
        characteristic=[[1]]
        for name,multiplicity in entry['factorization']:
            for _ in range(multiplicity):
                characteristic=x_poly_product(characteristic,
                  definition['factors'][name]['coefficients_desc_x_ascending_t'])
        # Raw matrix traces through n determine its characteristic polynomial.
        actual=direct_matrix_traces(rows,n,bits)
        predicted=trace_sequence(characteristic,n,bits=bits)
        if actual!=predicted:raise AssertionError('factorization trace check failed')
        # Factor coefficients and true characteristic coefficients have this
        # conservative l1 bound, so evaluation cannot hide coefficient errors.
        candidate_bound=sum(sum(abs(c) for c in a) for a in characteristic)
        actual_bound=factorial(n)*32**n
        if candidate_bound+actual_bound>=radix:raise AssertionError('factor bound too large')
        traces[n]=trace_sequence(characteristic,maximum,bits=bits)
        pfq=perron_quotient(rows)
        pfname={5:'closed_2',15:'open_5',16:'shared_4a'}[n]
        pfco=definition['factors'][pfname]['coefficients_desc_x_ascending_t']
        if direct_matrix_traces(pfq['rows'],len(pfq['rows']),bits)!=trace_sequence(pfco,len(pfq['rows']),bits=bits):
            raise AssertionError('Perron quotient characteristic mismatch')
        pfq['factor_name']=pfname
        blocks[size_text]=dict(state_indices=entry['quotient_state_indices'],rows=rows,perron_quotient=pfq,
          characteristic_coefficients_desc_x_ascending_t=characteristic,
          direct_trace_checks=n,irreducible_and_aperiodic_for_positive_t=True)
    v=[0]*dimension
    for mask,j in enumerate(q['initial']):v[j]+=1<<(bits*mask.bit_count())
    small=[]
    for m in range(1,maximum+1):
        p0=sum(vv for vv,r in zip(v,q['output']) if r==0)
        p2=sum(vv for vv,r in zip(v,q['output']) if r==2)
        monomial=2<<(2*m*bits)
        if p0!=traces[5][m]-traces[16][m]+monomial:raise AssertionError(('P0',m))
        if p2!=traces[15][m]-traces[16][m]+monomial:raise AssertionError(('P2',m))
        if m<=6:
            small.append({'length':m,'rank0_coefficients':[(p0>>(bits*k))&(radix-1) for k in range(4*m+1)],
                          'rank2_coefficients':[(p2>>(bits*k))&(radix-1) for k in range(4*m+1)]})
        if m<maximum:
            new=[0]*dimension
            for i,row in enumerate(q['rows']):
                for k,j,c in row:new[j]+=c*(v[i]<<(bits*k))
            v=new
    # Full-polynomial resultant certificate, using an exact Sylvester determinant.
    fc=definition['factors']['closed_2']['coefficients_desc_x_ascending_t']
    fo=definition['factors']['open_5']['coefficients_desc_x_ascending_t']
    c=[polynomial_value(a,radix) for a in fc];o=[polynomial_value(a,radix) for a in fo]
    sylvester=[]
    for k in range(2):sylvester.append([0]*k+o+[0]*(1-k))
    for k in range(5):sylvester.append([0]*k+c+[0]*(4-k))
    result=bareiss(sylvester)
    expected=radix**10*(radix+1)**2*polynomial_value(definition['crossing_polynomial_ascending_t'],radix)
    if result!=expected:raise AssertionError('resultant does not factor as claimed')
    norm_c=sum(sum(abs(v) for v in a) for a in fc);norm_o=sum(sum(abs(v) for v in a) for a in fo)
    bound=norm_o**2*norm_c**5+4*sum(abs(v) for v in definition['crossing_polynomial_ascending_t'])
    if bound>=radix:raise AssertionError('resultant coefficient bound too large')
    minimality=[]
    # Characteristic factors give upper bounds. These exact minors give matching
    # lower bounds, using the physical tail beginning at m=2, never a tolerance.
    for t,channel,order in [(1,'M',15),(2,'M',14),(3,'M',16),(3,'P0',17),(3,'P2',23)]:
        seq=scalar_sequences(q,2*order+1,t)[{'P0':0,'P2':1,'M':2}[channel]]
        det=bareiss([[seq[1+i+j] for j in range(order)] for i in range(order)])
        if det==0:raise AssertionError(('zero Hankel minor',t,channel,order))
        minimality.append(dict(fugacity=t,probability=str(Fraction(t,1+t)),channel=channel,
          minimum_scalar_order=order,hankel_start_length=2,determinant=str(det)))
    generic_factors=[definition['factors'][name]['coefficients_desc_x_ascending_t']
                     for name in factor_weights(definition,'M')]
    rec=[[1]]
    for f in generic_factors:rec=x_poly_product(rec,f)
    if len(rec)-1!=16 or not squarefree([polynomial_value(c,3) for c in rec]):
        raise AssertionError('generic squarefree recurrence proof failed')
    return {'schema':'matching-one.width4-parametric-traces.v1','source_head':definition['source_head'],
      'source_blob':blob,'weighted_quotient':q,'blocks':blocks,
      'identities':{'P0':'tr(B5^m)-tr(B16^m)+2*t^(2m)',
                    'P2':'tr(B15^m)-tr(B16^m)+2*t^(2m)',
                    'M':'tr(B15^m)-tr(B5^m)'},
      'normalization':'divide by (1+t)^(4*m); t=p/(1-p); endpoints by continuity',
      'trace_factor_weights':{ch:factor_weights(definition,ch) for ch in ('P0','P2','M')},
      'certificate':{'common_weighted_lumping_counts':q['counts'],
        'all_509_rows_checked':True,'zero_polynomials_checked_each':maximum,
        'P0_sufficient_moments':dimension+5+16+1,'P2_sufficient_moments':maximum,
        'radix_bits':bits,'difference_coefficient_bound':'34*16^m < radix',
        'method':'exact Kronecker evaluation with proved coefficient bound, then finite-dimensional Cayley-Hamilton',
        'small_block_raw_matrix_trace_checks':36,'resultant_sylvester_identity':True},
      'generic_minimum_orders':{'P0':17,'P2':23,'M':16},'minimality':minimality,
      'generic_M_recurrence_coefficients_desc_x_ascending_t':rec,
      'generic_recurrence_squarefree_over_Q_t':True,'small_integer_sector_polynomials':small,
      'seconds':time.perf_counter()-started,
      'limits':['width four only, physical m>=2','small-block traces are NOT individual probabilities',
        'generic rank is not uniform at exceptional parameters','not a continuum operator identification or all-width pTL intertwiner']}


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',type=Path,default=SOURCE)
    ap.add_argument('--definition',type=Path,default=DEFINITION);ap.add_argument('--out',type=Path)
    args=ap.parse_args();report=verify(args.source,args.definition)
    text=json.dumps(report,indent=2,allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open('x',encoding='utf-8') as f:f.write(text)
    else:print(text,end='')
if __name__=='__main__':main()
