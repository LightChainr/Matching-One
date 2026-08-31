#!/usr/bin/env python3
"""Exact finite Q1 activation of the prescribed antisymmetric colour trace.

Reads immutable N25 coefficients; performs no sampling or root search.
The root/sign decisions use Fraction intervals, never binary floats.
This standalone contribution does not import an unmerged research branch.
"""
from __future__ import annotations
import argparse
import csv
from fractions import Fraction as F
import hashlib
import io
import json
from math import comb, factorial
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'experiments/p337-antisymmetric-trace-20260831'
N = 25
DELTA = F(1152, 625)

class Interval:
    def __init__(self, lo, hi=None):
        self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError('reversed interval')
    @staticmethod
    def cast(x):
        return x if isinstance(x, Interval) else Interval(x)
    def __add__(self, other):
        b = self.cast(other)
        return Interval(self.lo+b.lo, self.hi+b.hi)
    __radd__ = __add__
    def __neg__(self):
        return Interval(-self.hi, -self.lo)
    def __sub__(self, other):
        return self + -self.cast(other)
    def __rsub__(self, other):
        return self.cast(other) + -self
    def __mul__(self, other):
        b = self.cast(other)
        z = (self.lo*b.lo, self.lo*b.hi, self.hi*b.lo, self.hi*b.hi)
        return Interval(min(z), max(z))
    __rmul__ = __mul__
    def reciprocal(self):
        if self.lo <= 0 <= self.hi:
            raise ZeroDivisionError('interval contains zero')
        return Interval(1/self.hi, 1/self.lo)
    def __truediv__(self, other):
        return self*self.cast(other).reciprocal()
    def __rtruediv__(self, other):
        return self.cast(other)*self.reciprocal()
    def __pow__(self, n):
        if not isinstance(n, int):
            raise TypeError('integer power required')
        if n < 0:
            return self.reciprocal()**(-n)
        r = Interval(1)
        for _ in range(n):
            r = r*self
        return r
    def compact(self, bits=128):
        scale = 1 << bits
        a = self.lo.numerator*scale//self.lo.denominator
        b = -((-self.hi.numerator*scale)//self.hi.denominator)
        return Interval(F(a, scale), F(b, scale))
    def record(self):
        a = self.compact()
        return {'lo': str(a.lo), 'hi': str(a.hi),
                'sign': 'positive' if a.lo > 0 else 'negative' if a.hi < 0 else 'unresolved'}

def blob_sha(data):
    return hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()

def polynomial(count_coefficients):
    """Convert sum c[k] p^k(1-p)^(N-k); c already includes multiplicity."""
    if len(count_coefficients) != N+1:
        raise ValueError('wrong Bernstein degree')
    return [sum(count_coefficients[k]*(-1)**(j-k)*comb(N-k, j-k)
                for k in range(j+1)) for j in range(N+1)]

def derivative(c):
    return [j*c[j] for j in range(1, len(c))]

def evaluate(c, x):
    y = 0*x
    for a in reversed(c):
        y = y*x+a
    return y

def jet(c, x):
    return [evaluate(c,x), evaluate(derivative(c),x),
            evaluate(derivative(derivative(c)),x)]

def integer_partitions(n, least=1):
    if n == 0:
        yield ()
    for k in range(least, n+1):
        for tail in integer_partitions(n-k, k):
            yield (k,)+tail

def character_coefficients(Q):
    """Exact conjugacy-class sums for exterior-square standard character."""
    out = {k:F(0) for k in ('constant','one_port','A','B','norm')}
    for parts in integer_partitions(Q):
        counts = {k:parts.count(k) for k in set(parts)}
        z = 1
        for k,c in counts.items():
            z *= k**c*factorial(c)
        x1, x2 = counts.get(1,0), counts.get(2,0)
        chi = F((x1-1)*(x1-2),2)-x2
        functions = {'constant':1,'one_port':x1,'A':x1*x1,
                     'B':x1+2*x2,'norm':chi}
        for name,f in functions.items():
            out[name] += chi*f/z
    return out

def load_inputs(data=DATA):
    source = json.loads((data/'SOURCES.json').read_text())
    for name, expected in source['input_sha256'].items():
        if hashlib.sha256((data/'inputs'/name).read_bytes()).hexdigest() != expected:
            raise ValueError('input digest mismatch: '+name)
    ordinary = {}
    for geometry in ('axis','tilted'):
        b = (data/'inputs'/('ordinary_'+geometry+'.csv')).read_bytes()
        if blob_sha(b) != source['ordinary_blobs'][geometry]:
            raise ValueError('ordinary source Git blob mismatch')
        rows = [{k:int(v) for k,v in r.items()} for r in csv.DictReader(io.StringIO(b.decode()))]
        if [r['k'] for r in rows] != list(range(N+1)):
            raise ValueError('missing occupation layer')
        if any(r['count'] != comb(N,r['k']) for r in rows):
            raise ValueError('incorrect finite-population multiplicity')
        mq = [F(r['sum_q'],r['count']) for r in rows]
        if not all(y>=x for x,y in zip(mq,mq[1:])) or mq[-1]<=mq[0]:
            raise ValueError('unique-root Bernstein monotonicity failed')
        ordinary[geometry] = rows
    with (data/'inputs/trace_types.csv').open(encoding='utf-8', newline='') as handle:
        trace = list(csv.DictReader(handle))
    for row in trace:
        for k in ('k','g','count'):
            row[k] = int(row[k])
        if row['geometry'] not in ordinary or row['type'] not in ('A','B') or row['count']<=0:
            raise ValueError('invalid trace row')
    if len({(r['geometry'],r['type'],r['k'],r['g']) for r in trace}) != len(trace):
        raise ValueError('duplicate compact trace row')
    for geometry in ordinary:
        for kind in ('A','B'):
            if sum(r['count'] for r in trace if r['geometry']==geometry and r['type']==kind) != source['type_totals'][geometry][kind]:
                raise ValueError('trace support total mismatch')
        for k in range(N+1):
            support = sum(r['count'] for r in trace if r['geometry']==geometry and r['k']==k)
            if support > ordinary[geometry][k]['count']-ordinary[geometry][k]['sum_e']:
                raise ValueError('trace exceeds rank-one support')
    bounds = json.loads((data/'inputs/root_bracket.json').read_text())['root_interval']
    return ordinary, trace, Interval(bounds['lo'],bounds['hi']), source

def check_source_checkout(checkout, data=DATA):
    """Optional exact re-extraction; needs the pinned source objects in Git."""
    _, stored, _, source = load_inputs(data)
    actual = []
    for geometry, path in source['source_histogram_paths'].items():
        b = subprocess.run(['git','-C',str(checkout),'show',source['science_commit']+':'+path],
                           check=True, capture_output=True).stdout
        if blob_sha(b) != source['source_histogram_blobs'][geometry]:
            raise ValueError('full histogram blob mismatch')
        for row in csv.DictReader(io.StringIO(b.decode())):
            r = {k:int(v) for k,v in row.items()}
            for kind, rule in source['type_rule'].items():
                if all(r[k]==v for k,v in rule.items()):
                    actual.append(dict(geometry=geometry,type=kind,k=r['k'],g=r['g'],count=r['count']))
    key = lambda r:(r['geometry'],r['type'],r['k'],r['g'],r['count'])
    if sorted(actual,key=key) != sorted(stored,key=key):
        raise ValueError('compact source extraction differs')
    return True

def build_result(data=DATA):
    ordinary, trace, p, source = load_inputs(data)
    qpoly = [polynomial([r['sum_q'] for r in ordinary[g]]) for g in ('axis','tilted')]
    epoly = [polynomial([r['sum_e'] for r in ordinary[g]]) for g in ('axis','tilted')]
    Q = [(x+y)/F(2) for x,y in zip(*qpoly)]
    if not evaluate(Q,p.lo)<0<evaluate(Q,p.hi):
        raise ValueError('inherited root bracket fails exact signs')
    q = [jet(c,p) for c in qpoly]; e = [jet(c,p) for c in epoly]
    D = (q[0][1]+q[1][1])/2
    T = (q[0][2]+q[1][2])/2
    B = (e[0][1]-e[1][1])/DELTA
    H = (e[0][2]-e[1][2])/DELTA
    if D.lo<=0:
        raise ValueError('pooled slope is not certified positive')
    sources = {'type_A':{'A':F(1),'B':F(0)}, 'type_B':{'A':F(0),'B':F(1)},
               'symmetric_control':{'A':F(-1),'B':F(-1)},
               'antisymmetric_activation':{'A':F(-1,2),'B':F(1,2)}}
    results = {}
    for name,coeff in sources.items():
        J = []; fractions = []
        for i,g in enumerate(('axis','tilted')):
            c = [F(0)]*(N+1)
            for row in trace:
                if row['geometry']==g:
                    c[row['k']] += coeff[row['type']]*row['count']
            b,bp,_ = jet(polynomial(c),p)
            m,mp,_ = q[i]; ee,ep,_ = e[i]
            J.append((-m*b,-mp*b-m*bp,-ee*b,-ep*b-ee*bp))
            fractions.append({'mean':b.record(),'p_derivative':bp.record()})
        jq = (J[0][0]+J[1][0])/2
        jqp = (J[0][1]+J[1][1])/2
        jep = (J[0][3]-J[1][3])/DELTA
        terms = [jep/D, -H*jq/D**2, -B*jqp/D**2, B*T*jq/D**3]
        v = sum(terms, Interval(0))
        mid = (v.lo+v.hi)/2
        results[name] = {'V_over_A':v.record(),
                         'V_display_only':float(mid)*(N**(13/8)/2),
                         'four_terms_over_A':[t.record() for t in terms],
                         'source_fractions':fractions}
    control = results['symmetric_control']['V_display_only']
    if abs(control-(-0.001904836180602413))>1e-15:
        raise ValueError('published symmetric control failed')
    chars = {}
    for Qn in range(4,10):
        answer = character_coefficients(Qn)
        expected = {'constant':0,'one_port':0,'A':1,'B':-1,'norm':1}
        if answer != expected:
            raise ValueError('character identity failed')
        chars[str(Qn)] = {k:str(v) for k,v in answer.items()}
    return {'schema':'matching-one/antisymmetric-trace-review/v1',
            'scope':'N25 fixed first deck seam, declared rational continuation of [Q-2,1,1]; not a continuum-field identification',
            'source_commit':source['science_commit'],
            'derivative_coordinate':'d/dlogQ at Q=1 of the affine closed-trace epsilon response; log(m) derivative is twice this',
            'root_interval':{'lo':str(p.lo),'hi':str(p.hi)},
            'pooled_slope':D.record(),'new_root_search':False,
            'new_samples':0,'new_enumerations':0,
            'arithmetic':'Fraction outward intervals; floating values are presentation only',
            'character_checks':chars,'responses':results,
            'decision':'nonzero_finite_activation' if results['antisymmetric_activation']['V_over_A']['sign']!='unresolved' else 'unresolved',
            'input_boundary':'Exact result conditional on pinned finite-population coefficients; no independent configuration re-enumeration'}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--verify',type=Path)
    parser.add_argument('--check-source-checkout',type=Path)
    args = parser.parse_args()
    if args.check_source_checkout:
        check_source_checkout(args.check_source_checkout)
    result = build_result()
    if args.verify and json.loads(args.verify.read_text()) != result:
        raise SystemExit('stored certificate differs from exact reconstruction')
    text = json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open('x') as f:
            f.write(text)
    elif not args.verify:
        print(text,end='')
    else:
        print('exact certificate reproduced')

if __name__=='__main__':
    main()
