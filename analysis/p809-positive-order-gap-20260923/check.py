#!/usr/bin/env python3
"""Positive order 94 versus linear order 93 at the degree-11 blind point.

A proof certificate for the whole family of width-four NN-site preparations,
including the current bit B0=1{rank=1}. No optimization is part of verification.
The D4/oracle and finite-field helpers are reused from PR831's check.py, pinned
at e099524bdc7e9c8844c56d88e5d7b411225f85c9, not a new rank implementation.
The new tests are the forbidden-language face, exact gap and positive left bound.

Run: python check.py input/width4-rank-closure-certificate.json --out reproduced.json
"""
from __future__ import annotations
import argparse, json, time
from collections import Counter, deque
from fractions import Fraction
from hashlib import sha1
from itertools import product
from math import isqrt
from pathlib import Path
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ, QQ, GF
from sympy.polys.matrices import DomainMatrix
x = sp.Symbol('x')

def require(test, message):
    if not test:
        raise AssertionError(message)


def labels(keys):
    index, out = {}, []
    for key in keys:
        key = tuple(key)
        if key not in index:
            index[key] = len(index)
        out.append(index[key])
    return out


def state_for(word, transitions, initial):
    s = initial[word[0]]
    for b in word[1:]:
        s = transitions[s][b]
    return s


def physical_histories(q, initial):
    words = {}
    for a, b in product(range(16), repeat=2):
        words.setdefault(q[initial[a]][b], [a, b])
    todo = deque(words)
    while todo:
        s = todo.popleft()
        for b, t in enumerate(q[s]):
            if t not in words:
                words[t] = words[s] + [b]
                todo.append(t)
    require(len(words) == len(q), 'physical reachability')
    return [words[s] for s in range(len(q))]


def d4_quotient(q, output, initial):
    histories = physical_histories(q, initial)
    actions = []
    for sign, shift in product((1, -1), range(4)):
        masks = [sum(((b >> x) & 1) << ((sign*x+shift) % 4)
                     for x in range(4)) for b in range(16)]
        action = [state_for([masks[b] for b in word], q, initial)
                  for word in histories]
        require(len(set(action)) == len(q), 'D4 not a permutation')
        for s in range(len(q)):
            require(output[action[s]] == output[s], 'D4 output')
            for b in range(16):
                require(action[q[s][b]] == q[action[s]][masks[b]], 'D4 transition')
        actions.append(action)
    lab = labels((min(a[s] for a in actions),) for s in range(len(q)))
    representatives = [lab.index(i) for i in range(len(set(lab)))]
    weighted = [Counter((b.bit_count(), lab[q[s][b]]) for b in range(16))
                for s in representatives]
    for s in range(len(q)):
        require(Counter((b.bit_count(), lab[q[s][b]]) for b in range(16)) == weighted[lab[s]],
                'all-p strong quotient')
    return lab, representatives, weighted, [output[s] for s in representatives], histories


def graph_rank(rows):
    """Independent physical lifted-graph DFS; not certificate transitions."""
    w,h,positions,first = 4,len(rows),{},None
    for y0 in range(h):
        for x0 in range(w):
            root = y0*w+x0
            if not(rows[y0]>>x0&1) or root in positions:
                continue
            positions[root]=(0,0);stack=[root]
            while stack:
                u=stack.pop();x,y=u%w,u//w;px,py=positions[u]
                for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                    nx,ny=(x+dx)%w,(y+dy)%h
                    if not(rows[ny]>>nx&1):
                        continue
                    v=ny*w+nx;proposed=(px+dx,py+dy)
                    if v not in positions:
                        positions[v]=proposed;stack.append(v)
                    else:
                        a,b=proposed[0]-positions[v][0],proposed[1]-positions[v][1]
                        require(a%w==0 and b%h==0, 'physical displacement')
                        if a or b:
                            if first is None:
                                first=(a,b)
                            elif first[0]*b-first[1]*a:
                                return 2
    return int(first is not None)


def prime_check(m):
    return m >= 2 and all(m % k for k in range(2,isqrt(m)+1))


def eval_mod(cs, t, m):
    a = 0
    for c in reversed(cs):
        a = (a*t+int(c)) % m
    return a


def determinant_mod(a, m):
    a = [[int(v)%m for v in row] for row in a]
    d = 1
    for k in range(len(a)):
        i = next((i for i in range(k,len(a)) if a[i][k]), None)
        if i is None:
            return 0
        if i != k:
            a[i],a[k] = a[k],a[i]; d = -d
        pivot = a[k][k]; d = d*pivot % m; inv = pow(pivot,-1,m)
        for i in range(k+1,len(a)):
            c = a[i][k]*inv % m
            if c:
                for j in range(k+1,len(a)):
                    a[i][j] = (a[i][j]-c*a[k][j]) % m
            a[i][k] = 0
    return d % m


def observer_mod(weighted, output, m, t):
    """Exact finite-field lower bound; return a nonzero actual word minor."""
    n = len(output)
    rows = [[(j,mult*pow(t,k,m)%m) for (k,j),mult in row.items()] for row in weighted]
    echelon, pivots, vectors, words = [], [], [], []
    def add(v, word):
        u = v[:]
        for i,b in zip(pivots,echelon):
            c = u[i]
            if c:
                u = [(z-c*y)%m for z,y in zip(u,b)]
        i = next((i for i,z in enumerate(u) if z), None)
        if i is None:
            return
        inv = pow(u[i],-1,m)
        echelon.append([z*inv%m for z in u]); pivots.append(i)
        vectors.append(v); words.append(word)
    for b in sorted(set(output)):
        add([int(z==b) for z in output], [b])
    k = 0
    while k < len(vectors) and len(vectors) < n:
        av = [sum(c*vectors[k][j] for j,c in row)%m for row in rows]
        for b in sorted(set(output)):
            add([z if output[i]==b else 0 for i,z in enumerate(av)], [b]+words[k])
        k += 1
    minor = [[v[i] for v in vectors] for i in pivots]
    d = determinant_mod(minor,m)
    require(d != 0, 'word-minor lower bound failed')
    return {'dimension':len(vectors), 'prime':m, 'root':t, 'minor_determinant_mod':d,
            'maximum_word_length':max(map(len,words)), 'minor_rows':pivots,
            'words':words}


def algebraic_field(cs):
    f = sp.Poly.from_list(cs[::-1],x,domain=ZZ)
    field = QQ.algebraic_field(sp.CRootOf(f.as_expr(),0))
    return field, field.from_sympy(field.ext.as_expr())


def field_value(cs, t, field):
    result = field.zero
    for c in reversed(cs): result = result*t+field.convert(QQ(c))
    return result


def fraction_value(cs, t):
    out = Fraction(0)
    for c in reversed(cs):
        out = out*t + Fraction(c)
    return out


def interval_multiply(a, b):
    values = [u*v for u in a for v in b]
    return min(values), max(values)


def variation(signs):
    nonzero = [a for a in signs if a]
    return sum(a != b for a, b in zip(nonzero, nonzero[1:]))


def run(path, seed):
    started = time.monotonic()
    raw = path.read_bytes()
    blob = sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    require(blob == seed['input_blob'], 'input is not the pinned PR708 certificate')
    source = json.loads(raw)
    q = source['quotient_transitions']; initial = source['quotient_initial']
    lab, reps, weighted, ranks, histories = d4_quotient(
        q, source['quotient_rank_output'], initial)
    require(len(ranks) == 94, 'D4 quotient changed')
    bits = [int(r == 1) for r in ranks]
    record = seed['invisible']; n = len(ranks)
    f = sp.Poly.from_list(record['factor'][::-1], x, domain=ZZ)

    # Exact irreducibility makes the good finite-field minors characteristic-zero
    # lower bounds at every root of f, not merely numerical ranks at one prime.
    possible = set(range(1, f.degree()))
    for item in seed['irreducibility']:
        m = item['prime']; require(prime_check(m), 'not a prime')
        fm = sp.Poly(f.as_expr(), x, modulus=m)
        require(fm.degree() == f.degree() and fm.gcd(fm.diff()).degree() == 0,
                'bad or repeated-root modular reduction')
        _, factors = fm.factor_list()
        degrees = sorted(g.degree() for g, e in factors for _ in range(e))
        require(degrees == sorted(item['degrees']), 'irreducibility witness changed')
        sums = {0}
        for d in degrees:
            sums |= {a+d for a in sums}
        possible &= sums
    require(not possible, 'rational irreducibility not proved')
    require(f.count_roots(-sp.oo, 1) == 0 and f.count_roots(1, 2) == 1
            and f.count_roots(2, sp.oo) == 0, 'wrong real algebraic embedding')
    lo, hi = Fraction(1), Fraction(2)
    require(fraction_value(record['factor'], lo) < 0
            < fraction_value(record['factor'], hi), 'root is not bracketed')
    for _ in range(seed['root_bisections']):
        mid = (lo+hi)/2
        if fraction_value(record['factor'], mid) > 0:
            hi = mid
        else:
            lo = mid
    field, t = algebraic_field(record['factor'])
    den = (1+t)**4

    def element(c):
        c = Fraction(c)
        return field.convert(QQ(c.numerator, c.denominator))

    def sign(a):
        """Prove signs by rational interval Horner evaluation at the isolated root."""
        if a == field.zero:
            return 0
        interval = (Fraction(0), Fraction(0))
        for c in a.to_list():  # ANP coefficient order is descending.
            interval = interval_multiply(interval, (lo, hi))
            cc = Fraction(str(c))
            interval = (interval[0]+cc, interval[1]+cc)
        if interval[0] > 0:
            return 1
        if interval[1] < 0:
            return -1
        raise AssertionError('algebraic sign unresolved by the declared rational interval')

    def matrix(I, J):
        return [[sum((multiplicity*t**k for (k, target), multiplicity in weighted[i].items()
                       if target == j), field.zero)/den for j in J] for i in I]

    v = [field.zero]*n
    for i, cs in zip(record['support'], record['vectors']):
        v[i] = field_value(cs, t, field)
    lam = field_value(record['eigenvalue'], t, field)/den
    require(sum(v, field.zero) == field.zero and any(v), 'invalid invisible contrast')
    require(all(bits[i] == 0 for i in record['support']), 'contrast is not B=0 supported')
    image = [field.zero]*n
    for i in record['support']:
        for (k, j), multiplicity in weighted[i].items():
            image[j] += v[i]*multiplicity*t**k/den
    require(all(a == lam*b for a, b in zip(image, v)), 'vK=lambda v failed')
    require(sign(-lam-element(seed['mode_magnitude_lower'])) == 1,
            'the invisible eigenvalue is not below the negative cutoff')

    gm = seed['global_minor']; m, a = gm['prime'], gm['odds']
    require(prime_check(m) and eval_mod(record['factor'], a, m) == 0
            and (a+1) % m != 0, 'invalid global minor reduction')
    global_lower = observer_mod(weighted, bits, m, a)
    require(global_lower['dimension'] == 93, 'global 93-dimensional lower bound changed')
    # The exact nonzero invisible vector is the matching upper bound.
    print('PASS exact linear dimension 93 and one-dimensional invisible space', flush=True)

    G0 = seed['good_zero']; absorber = seed['absorber']; bridge = seed['bridge_one']
    G = G0+[absorber]; J = seed['resonant_zero']
    require(len(set(G0)) == 7 and len(set(J)) == 10 and not set(G) & set(J),
            'incorrect face/core index sets')
    require(all(bits[i] == 0 for i in G0+J) and bits[absorber] == bits[bridge] == 1,
            'output typing failed')
    require(all(j in G for i in G0 for k, j in weighted[i]), 'good physical family not closed')
    require(all(j == absorber for k, j in weighted[absorber]), 'rank-one state not absorbing')
    require(all(j in set(J)|{absorber, bridge} for i in J for k, j in weighted[i]),
            'resonant core has an unaccounted exit')
    require(set(record['support']) <= set(G0)|set(J)
            and any(v[i] for i in G0) and any(v[i] for i in J),
            'contrast no longer straddles the two cores')

    # 85 forbidden-word columns. A word is forbidden precisely when it contains 10.
    fm = seed['forbidden_minor']; m, a = fm['prime'], fm['odds']
    require(prime_check(m) and eval_mod(record['factor'], a, m) == 0
            and (a+1) % m != 0, 'invalid forbidden minor reduction')
    rows = [[(j, multiplicity*pow(a,k,m) % m) for (k,j),multiplicity in row.items()]
            for row in weighted]
    cache = {}
    def column(word):
        word = tuple(word)
        if word not in cache:
            if len(word) == 1:
                out = [int(b == word[0]) for b in bits]
            else:
                tail = column(word[1:])
                out = [sum(w*tail[j] for j,w in row) % m if bits[i] == word[0] else 0
                       for i,row in enumerate(rows)]
            cache[word] = out
        return cache[word]
    words = fm['bad_words']; pivots = fm['minor_rows']
    require(len(words) == len(pivots) == 85 and len(set(pivots)) == 85,
            'not an 85-by-85 forbidden minor')
    require(all(any(a == 1 and b == 0 for a,b in zip(w,w[1:])) for w in words),
            'one alleged forbidden word contains no 10')
    vectors = [column(w) for w in words]
    require(all(vec[i] == 0 for vec in vectors for i in G), 'forbidden event on the good face')
    minor = [[vec[i] for vec in vectors] for i in pivots]
    determinant = determinant_mod(minor, m)
    require(determinant == fm['minor_determinant'] != 0, 'forbidden minor is singular')
    # A second arithmetic implementation checks both nonzero finite-field determinants.
    ff = GF(m)
    other_bad = int(DomainMatrix([[ff(z) for z in row] for row in minor], (85,85), ff).det())
    require(other_bad == determinant, 'independent forbidden determinant disagrees')
    require(gm['prime'] == m and gm['odds'] == a, 'global cross-check uses a different field')
    gvectors = [column(w) for w in global_lower['words']]
    gmatrix = [[ff(vec[i]) for vec in gvectors] for i in global_lower['minor_rows']]
    other_global = int(DomainMatrix(gmatrix, (93,93), ff).det())
    require(other_global == global_lower['minor_determinant_mod'],
            'independent global determinant disagrees')
    # Eight good rows are zero, and v restricted to their complement is another
    # exact nonzero dependence. The upper bound is therefore 94-8-1=85.
    print('PASS forbidden-language dimension 85; eight independent good preparations', flush=True)

    Rmat = matrix(G0, G0); Bmat = matrix(J, J)
    vg = [v[i] for i in G0]
    require(all(sum((vg[i]*Rmat[i][j] for i in range(7)), field.zero) == lam*vg[j]
                for j in range(7)), 'good-face restriction is not a lambda mode')
    weights = seed['left_supervector']
    require(len(weights) == 10 and all(isinstance(c,int) and c > 0 for c in weights),
            'left supervector is not strictly positive')
    upper = element(seed['upper_cut'])
    margins = [upper*weights[j]-sum((weights[i]*Bmat[i][j] for i in range(10)), field.zero)
               for j in range(10)]
    require(all(sign(z) == 1 for z in margins), 'positive left growth bound failed')

    coefficients = DomainMatrix(Rmat, (7,7), field).charpoly()
    ring_z, z = ring('z', field)
    characteristic = sum((c*z**(7-k) for k,c in enumerate(coefficients)), ring_z.zero)
    sturm = [characteristic, characteristic.diff(z)]
    while sturm[-1]:
        remainder = -(sturm[-2] % sturm[-1])
        if not remainder:
            break
        sturm.append(remainder)
    signs = [[sign(pol.evaluate(z, element(cut))) for pol in sturm]
             for cut in (seed['low_cut'], seed['upper_cut'])]
    require(signs == seed['expected_sturm_signs'] and all(row[0] for row in signs),
            'Sturm endpoint signs changed or an endpoint is a root')
    counts = [variation(row) for row in signs]
    require(counts[0] == counts[1], 'good-face spectrum has a real root in the gap')
    print('PASS exact Perron obstruction: no real root in [1/50,19/25], '
          'resonant growth <19/25, |lambda|>1/20', flush=True)

    # Physical check proportional to the reused geometry: the two named cores,
    # their bridge and absorber, each with all one-row suffixes.
    physical = sorted(set(G+J+[bridge])); graph_checks = 0
    for i in physical:
        h = histories[reps[i]]
        for suffix in ([], *[[b] for b in range(16)]):
            word = h+suffix
            s = state_for(word, q, initial)
            require(graph_rank(word) == ranks[lab[s]], 'independent physical graph mismatch')
            graph_checks += 1
    prepared = {s:[b] for b,s in enumerate(initial)}
    for _ in range(4):
        nxt = {}
        for s,h in prepared.items():
            for b,j in enumerate(q[s]):
                nxt.setdefault(j,h+[b])
        prepared = nxt
    same_height = {}
    for s,h in prepared.items():
        same_height.setdefault(lab[s], h)
    require(len(same_height) == 94, 'not all preparations are available at height five')

    result = {
        'schema':'matching-one.p809-positive-gap.v1',
        'input_blob':blob,
        'parameter_odds_polynomial':record['factor'],
        'root_interval_odds':[str(lo),str(hi)],
        'parameter_decimal_diagnostic':str(sp.N(sp.Rational(lo.numerator,lo.denominator)/
                                      (1+sp.Rational(lo.numerator,lo.denominator)),30)),
        'contract':'all 94 history preparations and arbitrary probability mixtures; '
                   'time-homogeneous nonnegative edge-emitting model; current B0 included',
        'linear_dimension':93,
        'minimal_positive_states_binary':94,
        'minimal_positive_states_full_rank':94,
        'global_lower_bound_minor':global_lower,
        'independent_arithmetic_minors_mod':{'prime':m,'global':other_global,'forbidden':other_bad},
        'forbidden_pattern':'10',
        'forbidden_word_space_dimension':85,
        'forbidden_minor':{'prime':m,'odds_root':a,'determinant_mod':determinant,
                           'max_word_length':max(map(len,words))},
        'good_zero_states':G0,'absorber':absorber,'resonant_zero_states':J,'bridge_one':bridge,
        'good_histories':{str(i):histories[reps[i]] for i in G},
        'resonant_histories':{str(i):histories[reps[i]] for i in J},
        'left_supervector':weights,
        'left_bound':seed['upper_cut'],'left_bound_margin_signs':[sign(a) for a in margins],
        'spectral_gap_endpoints':[seed['low_cut'],seed['upper_cut']],
        'sturm_signs':signs,'sturm_variations':counts,
        'mode_magnitude_lower':seed['mode_magnitude_lower'],
        'characteristic_polynomial_coefficients_algebraic_descending':
            [[str(a) for a in c.to_list()] for c in coefficients],
        'invisible_eigenvalue_coefficients_algebraic_descending':[str(a) for a in lam.to_list()],
        'physical_graph_checks':graph_checks,
        'same_height_five_preparations':[same_height[i] for i in range(94)],
        'proof_status':'author proof with exact finite certificate; not independently externally reviewed',
        'not_claimed':['one fixed natural preparation also needs 94 states',
                       'delayed recording at p_a needs 94 states',
                       'positive order at the degree-17 point',
                       'new general positive-realization theory',
                       'approximate-model state lower bound'],
        'not_run':['full repository CI','new width','Monte Carlo','paid compute'],
        'sympy_version':sp.__version__,
    }
    print('PASS all checks; graph comparisons',graph_checks,
          '; elapsed seconds',time.monotonic()-started,flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--certificate', type=Path,
                        default=Path(__file__).with_name('certificate.json'))
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.input, json.loads(args.certificate.read_text()))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False)+'\n')
