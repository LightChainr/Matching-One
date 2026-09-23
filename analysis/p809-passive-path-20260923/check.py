#!/usr/bin/env python3
"""Passive rank-trajectory observability on the pinned PR708 width-four machine.
Python 3.10+ standard library. Reads the original certificate; no new width engine.
Usage: python check.py CERTIFICATE --out result.json
"""
from __future__ import annotations
import argparse
from collections import Counter, deque
from fractions import Fraction
from hashlib import sha1
from itertools import combinations, product
import json
from math import comb, gcd, lcm
from pathlib import Path
import sys
import time

INPUT_BLOB = '50b7297deefe7c50215aea2ed534ca5810461af3'


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def labels(keys):
    ids, out = {}, []
    for key in keys:
        key = tuple(key)
        if key not in ids:
            ids[key] = len(ids)
        out.append(ids[key])
    return out


def state_for(word, q, initial):
    s = initial[word[0]]
    for b in word[1:]:
        s = q[s][b]
    return s


def histories(q, initial):
    found = {}
    for a, b in product(range(16), repeat=2):
        found.setdefault(q[initial[a]][b], [a, b])
    todo = deque(found)
    while todo:
        s = todo.popleft()
        for b, t in enumerate(q[s]):
            if t not in found:
                found[t] = found[s] + [b]
                todo.append(t)
    require(len(found) == len(q), 'physical history missing')
    return [found[s] for s in range(len(q))]


def strong_lump(q, rank):
    old = rank[:]
    while True:
        new = labels((old[s], tuple(sorted(Counter(old[t] for t in row).items())))
                     for s, row in enumerate(q))
        if len(set(new)) == len(set(old)):
            return new
        old = new


def apply(q, vectors):
    return [[sum(vectors[t][j] for t in row) for j in range(len(vectors[0]))]
            for row in q]


def two_time(q, rank, a, b):
    require(0 <= a < b, 'invalid ordered snapshot times')
    c = [[int(r == j) for j in range(3)] for r in rank]
    for _ in range(b - a):
        c = apply(q, c)
    c = [[int(rank[s] == i) * c[s][j] for i in range(3) for j in range(3)]
         for s in range(len(q))]
    for _ in range(a):
        c = apply(q, c)
    require(all(sum(row) == 16**b for row in c), 'joint-law normalization')
    return c


def full_paths(q, rank, horizon):
    # Integer counts for words (Y_0,...,Y_h), lexicographic base-three ordering.
    v = [[int(r == j) for j in range(3)] for r in rank]
    table = []
    for h in range(horizon + 1):
        table.append(len(set(labels(v))))
        if h == horizon:
            return v, table
        av = apply(q, v)
        v = [[x if rank[s] == j else 0 for j in range(3) for x in av[s]]
             for s in range(len(q))]


class IntegerBasis:
    """Fraction-free exact Q-span; primitive integer pivots, never tolerance rank."""
    def __init__(self):
        self.pivots = {}

    def add(self, vector):
        v = vector[:]
        for i, b in self.pivots.items():
            if v[i]:
                c, d = v[i], b[i]
                v = [x*d - y*c for x, y in zip(v, b)]
                g = gcd(*v)
                if g > 1:
                    v = [x//g for x in v]
        nonzero = [i for i, x in enumerate(v) if x]
        if not nonzero:
            return False
        i = nonzero[0]
        g = gcd(*v) * (1 if v[i] > 0 else -1)
        self.pivots[i] = [x//g for x in v]
        return True


def path_basis(q, rank):
    basis, columns, frontier = IntegerBasis(), [], []
    for j in range(3):
        v = [int(r == j) for r in rank]
        require(basis.add(v), 'initial output not independent')
        columns.append(((j,), v))
        frontier.append(((j,), v))
    dimensions = [len(columns)]
    # Process only newly independent columns. Descendants of dependent columns
    # lie in descendants of this basis, by linearity. Old descendants were handled.
    for h in range(1, len(q) + 1):
        new = []
        for word, v in frontier:
            av = [sum(v[t] for t in row) for row in q]
            for j in range(3):
                u = [x*int(r == j) for x, r in zip(av, rank)]
                if basis.add(u):
                    new.append(((j,) + word, u))
                    columns.append(((j,) + word, u))
        dimensions.append(len(columns))
        if not new or len(columns) == len(q):
            return columns, dimensions
        frontier = new
    raise AssertionError('observable-space iteration did not stop')


def det_mod(rows, prime=65521):
    # A nonzero residue of a square integer determinant certifies full rational rank.
    a = [[x % prime for x in row] for row in rows]
    d, n = 1, len(a)
    for c in range(n):
        p = next((r for r in range(c, n) if a[r][c]), None)
        if p is None:
            return 0
        if p != c:
            a[c], a[p] = a[p], a[c]
            d = -d
        pivot = a[c][c]
        d = d*pivot % prime
        inv = pow(pivot, -1, prime)
        a[c] = [x*inv % prime for x in a[c]]
        for r in range(c + 1, n):
            v = a[r][c]
            if v:
                a[r] = [(x-v*y) % prime for x, y in zip(a[r], a[c])]
    return d % prime


def one_nullvector(rows):
    # Used only for the 93x94 short-horizon annihilator. Exact integer forward
    # elimination followed by one rational back-substitution; no SymPy dependency.
    basis = IntegerBasis()
    for row in rows:
        basis.add(row)
    n = len(rows[0])
    free = [i for i in range(n) if i not in basis.pivots]
    require(len(free) == 1, 'expected a one-dimensional nullspace')
    v = [Fraction(0) for _ in range(n)]
    v[free[0]] = Fraction(1)
    for i in sorted(basis.pivots, reverse=True):
        row = basis.pivots[i]
        v[i] = -sum((Fraction(row[j])*v[j] for j in range(i+1, n)), Fraction(0))/row[i]
    d = lcm(*(x.denominator for x in v))
    z = [int(x*d) for x in v]
    g = gcd(*z)
    z = [x//g for x in z]
    if next(x for x in z if x) > 0:
        z = [-x for x in z]
    require(all(sum(x*y for x, y in zip(row, z)) == 0 for row in rows), 'nullspace residual')
    return z


def graph_rank(rows, width=4):
    """Independent physical lifted-graph traversal, not a certificate transition."""
    height, seen, first = len(rows), {}, None
    for y in range(height):
        for x in range(width):
            root = width*y+x
            if not (rows[y] >> x & 1) or root in seen:
                continue
            seen[root] = (0, 0)
            stack = [root]
            while stack:
                u = stack.pop()
                x0, y0 = u % width, u // width
                px, py = seen[u]
                for dx, dy in ((1,0),(-1,0),(0,1),(0,-1)):
                    xx, yy = (x0+dx) % width, (y0+dy) % height
                    if not (rows[yy] >> xx & 1):
                        continue
                    v, at = width*yy+xx, (px+dx, py+dy)
                    if v not in seen:
                        seen[v] = at
                        stack.append(v)
                    else:
                        a, b = at[0]-seen[v][0], at[1]-seen[v][1]
                        require(a % width == 0 and b % height == 0, 'invalid deck displacement')
                        if a or b:
                            if first is None:
                                first = (a, b)
                            elif first[0]*b-first[1]*a:
                                return 2
    return int(first is not None)


def bernstein_power(hist, sites):
    out = [0]*(sites+1)
    for k, c in enumerate(hist):
        for j in range(sites-k+1):
            out[k+j] += c*comb(sites-k, j)*(-1)**j
    return out


def physical_controls(q, rank, initial):
    checks, pair12, pair13, hist12 = 0, [], [], []
    for history in ([0,7], [7,0]):
        s0 = state_for(history, q, initial)
        cache = {}
        for depth in (1, 2, 3):
            for future in product(range(16), repeat=depth):
                s = s0
                for b in future:
                    s = q[s][b]
                r = graph_rank(history + list(future))
                require(r == rank[s], 'physical prefix/certificate discrepancy')
                cache[future] = r
                checks += 1
        p12, p13 = [0]*9, [0]*9
        hist = [[0]*9 for _ in range(9)]
        for a, b in product(range(16), repeat=2):
            idx = 3*cache[(a,)]+cache[(a,b)]
            p12[idx] += 1
            hist[idx][a.bit_count()+b.bit_count()] += 1
        for a,b,c in product(range(16), repeat=3):
            p13[3*cache[(a,)]+cache[(a,b,c)]] += 1
        pair12.append(p12); pair13.append(p13); hist12.append(hist)
    dpoly = [0,0,0,1,-1,-1,1,0,0]  # p^3 (1-p)^2 (1+p)
    signs = [1,-1,0,-1,1,0,0,0,0]
    for j, sign in enumerate(signs):
        diff = [a-b for a,b in zip(hist12[0][j], hist12[1][j])]
        require(bernstein_power(diff, 8) == [sign*x for x in dpoly], 'passive all-p difference')
    for k, history in enumerate(([0,7],[7,0])):
        s = state_for(history, q, initial)
        require(pair12[k] == two_time(q,rank,1,2)[s], 'joint 1,2 mismatch')
        require(pair13[k] == two_time(q,rank,1,3)[s], 'joint 1,3 mismatch')
    return checks, pair12, pair13


def run_memory(q, rank, initial):
    s = state_for([7,0], q, initial)
    t = state_for([7,0,13], q, initial)
    absorb = state_for([0,15], q, initial)
    require(rank[s] == 0 and rank[t] == rank[absorb] == 1, 'run-state outputs')
    require(all(v == absorb for v in q[absorb]), 'not absorbing')
    first = {b:q[s][b] for b in range(16) if rank[q[s][b]] == 1}
    require(first == {13:t,15:absorb}, 'entry into all-one run')
    transient = [b for b in range(16) if q[t][b] == t]
    absorbed = [b for b in range(16) if q[t][b] == absorb]
    exits = [b for b in range(16) if rank[q[t][b]] == 0]
    require(transient == [5,13] and absorbed == [7,15] and len(exits) == 12, 'run transition masks')
    require(set(transient+absorbed+exits) == set(range(16)), 'unaccounted run transition')
    for masks, expected in ((transient,[0,0,1,-1,0]),(absorbed,[0,0,0,1,0]),(exits,[1,0,-1,0,0])):
        hist = [sum(b.bit_count() == k for b in masks) for k in range(5)]
        require(bernstein_power(hist, 4) == expected, 'all-p run transition polynomial')
    checks = 0
    for word in ([7,0],[7,0,13],[0,15]):
        u = state_for(word, q, initial)
        for b in range(16):
            require(graph_rank(word+[b]) == rank[q[u][b]], 'physical run transition')
            checks += 1
    examples = []
    for p in (Fraction(1,3),Fraction(1,2),Fraction(2,3)):
        weights = [p**b.bit_count()*(1-p)**(4-b.bit_count()) for b in range(16)]
        v = Counter()
        for b,u in first.items():
            v[u] += weights[b]
        rho = p*p*(1-p)
        hazards = []
        for h in range(1,13):
            den = sum(v.values())
            num = sum(value*sum(weights[b] for b,u in enumerate(q[state]) if rank[u] == 0)
                      for state,value in v.items())
            expected_den = (p**4+p**3*(1-p)*(1-p*p)*rho**(h-1))/(1-rho)
            expected_hazard = (1-rho)*(1-p)*(1-p*p)*rho**(h-1)/(p+(1-p)*(1-p*p)*rho**(h-1))
            require(den == expected_den and num/den == expected_hazard, 'run formula mismatch')
            require(num/den <= Fraction(4,27)**(h-1), 'uniform run-hazard upper bound')
            hazards.append(str(num/den))
            new = Counter()
            for state,value in v.items():
                for b,u in enumerate(q[state]):
                    if rank[u] == 1:
                        new[u] += value*weights[b]
            v = new
        examples.append({'p':str(p),'hazards_n1_to12':hazards})
    return checks, {'history':[7,0], 'transient_history':[7,0,13], 'absorbing_history':[0,15],
                   'entry_transient_mask':13,'entry_absorb_mask':15,
                   'transient_stay_masks':transient,'transient_to_absorb_masks':absorbed,
                   'rho':'p^2*(1-p)','exit_probability':'1-p^2','absorb_probability':'p^3',
                   'survival':'[p^4+p^3*(1-p)*(1-p^2)*rho^(n-1)]/(1-rho)',
                   'hazard':'(1-rho)*(1-p)*(1-p^2)*rho^(n-1)/[p+(1-p)*(1-p^2)*rho^(n-1)]',
                   'uniform_hazard_bound':'(4/27)^(n-1), n>=1', 'examples':examples}


def run(path):
    raw = path.read_bytes()
    blob = sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    require(blob == INPUT_BLOB, 'not the pinned original certificate')
    data = json.loads(raw)
    q, rank, initial = data['quotient_transitions'],data['quotient_rank_output'],data['quotient_initial']
    require(len(q) == 509 and all(len(row) == 16 for row in q), 'machine shape')
    words = histories(q,initial)
    common = {state:[b] for b,state in enumerate(initial)}
    common_length = 1
    while len(common) < len(q):
        new_common = {}
        for state,word in common.items():
            for b,target in enumerate(q[state]):
                new_common.setdefault(target,word+[b])
        common = new_common
        common_length += 1
        require(common_length <= 6, 'no bounded common-length preparation')
    require(common_length == 5, 'common-length preparation changed')
    require(all(len(w) == 5 and state_for(w,q,initial) == s for s,w in common.items()),
            'common-length representative mismatch')
    lump = strong_lump(q,rank)
    require(len(set(lump)) == 94, 'parent strong lump')
    actions=[]
    for sign,shift in product((1,-1),range(4)):
        bm = [sum(1 << ((sign*x+shift)%4) for x in range(4) if b>>x&1) for b in range(16)]
        perm = [state_for([bm[b] for b in w],q,initial) for w in words]
        require(sorted(perm) == list(range(509)), 'D4 action not bijective')
        require(all(rank[perm[s]] == rank[s] and perm[q[s][b]] == q[perm[s]][bm[b]]
                    for s in range(509) for b in range(16)), 'D4 output-transition covariance')
        actions.append(perm)
    d4 = labels((min(a[s] for a in actions),) for s in range(509))
    require(labels((x,) for x in lump) == d4, 'strong lump not D4')
    reps = [lump.index(i) for i in range(94)]
    qq = [[lump[t] for t in q[s]] for s in reps]
    rr = [rank[s] for s in reps]
    # All 509 states, not merely the selected representative, must have identical
    # successor aggregate counts in their class. This is the probability quotient.
    for s,row in enumerate(q):
        require(rank[s] == rr[lump[s]] and Counter(lump[t] for t in row) == Counter(qq[lump[s]]),
                'probability quotient invalid')
    paths2, classes2 = full_paths(q,rank,2)
    paths3, classes3 = full_paths(q,rank,3)
    require(classes3 == [3,35,91,94], 'rank-path class sequence')
    pair_data=[]
    for a,b in ((1,2),(1,3),(2,3)):
        table = two_time(q,rank,a,b)
        count = len(set(labels(table)))
        if b == 3:
            require(labels(table) == d4, 'passive pair does not separate D4')
        min_sum = min(sum(abs(x-y) for x,y in zip(table[reps[s]],table[reps[t]]))
                      for s,t in combinations(range(94),2))
        pair_data.append({'times':[a,b],'classes':count,'minimum_TV_between_D4_classes':str(Fraction(min_sum,2*16**b))})
    # Earliest possible horizon: these two states agree on the entire rank path to 2.
    delay_words = [[3,15,5],[5,15,3]]
    sa,sb = [state_for(w,q,initial) for w in delay_words]
    require(lump[sa] != lump[sb] and paths2[sa] == paths2[sb] and paths3[sa] != paths3[sb],
            'minimal-horizon witness failed')
    ep_columns,ep_dimensions = endpoint_basis(qq,rr)
    require(ep_dimensions == [3,5,7,9,11,13,15,17,19,21,23,25,27,29,30,30], 'marginal space dimension')
    conditioned = IntegerBasis()
    for v in ep_columns:
        for j in range(3):
            conditioned.add([int(r == j)*x for r,x in zip(rr,v)])
    require(len(conditioned.pivots) == 53, 'conditional marginal space dimension')
    pc,pe,pd = pair_basis(qq,rr,10)
    require(pd == [3,9,19,33,50,69,78,86,92,93,94], 'two-time space dimension')
    pair_determinant = det_mod(pc)
    require(pair_determinant != 0, 'two-time full-rank minor')
    small_count,small = small_mixture(q,rank,initial,lump,reps,ep_columns,common)
    columns,dimensions = path_basis(qq,rr)
    require(dimensions == [3,9,25,53,76,86,90,91,92,93,94], 'exact observable dimension')
    require(len(columns) == 94, 'not full observable space')
    mat = [col for word,col in columns]
    determinant = det_mod(mat)
    require(determinant != 0, 'modular determinant certificate failed')
    z = one_nullvector([col for word,col in columns if len(word) <= 10])
    require(sum(z) == 0, 'signed witness not a difference of measures')
    positive = sum(x for x in z if x>0)
    last_word,last_col = columns[-1]
    dot = sum(x*y for x,y in zip(z,last_col))
    require(dot != 0 and len(last_word) == 11, 'delayed mixture not separated at 10')
    tv_curve,contrasts = signed_path_control(qq,rr,z,10)
    gap = Fraction(abs(dot),positive*16**10)
    require(all(x == 0 for x in tv_curve[:10]) and tv_curve[10] == gap, 'forward full-path TV')
    require(len(contrasts) == 2, 'unexpected full-path contrast support')
    required = Fraction(1,2)/gap
    sample_lower = (required.numerator+required.denominator-1)//required.denominator
    mixture = {'positive_total':positive,'negative_total':-sum(x for x in z if x<0),
               'common_preparation_length':common_length,
               'signed_histories':[{'history':words[reps[i]],'common_length_history':common[reps[i]],'weight':x}
                                   for i,x in enumerate(z) if x],
               'all_rank_paths_through_horizon_9_equal':True,
               'separating_rank_word':list(last_word),
               'difference_on_separating_word':str(Fraction(dot,positive*16**10)),
               'full_path_TV_h0_to10':[str(x) for x in tv_curve],
               'equal_prior_error_at_most_one_quarter_necessary_independent_trials':sample_lower,
               'sample_bound_scope':'fixed horizon 10; an independent fresh mixed preparation per trial; these two hypotheses only'}
    physical_count,p12,p13 = physical_controls(q,rank,initial)
    require(Fraction(sum(abs(a-b) for a,b in zip(p12[0],p12[1])),512) == Fraction(3,32), 'passive witness TV')
    memory_count,memory = run_memory(q,rank,initial)
    absorb = state_for([0,15],q,initial)
    for a,b in ((0,15),(15,0)):
        require(all(q[q[s][a]][b] == absorb for s in range(509)), 'reset pair failed')
    # Validate the terminal-forgetting bound at p=1/2 directly on all 509 states.
    nonone = [[int(r != 1)] for r in rank]
    for h in range(101):
        rhs = 2*15**h-14**h
        require(max(x[0] for x in nonone) <= rhs, 'uniform terminal bound failed')
        nonone = apply(q,nonone)
    cutoff = next(h for h in range(1000) if Fraction(2*15**h-14**h,16**h) <= Fraction(1,100))
    return {'schema':'matching-one.p809-passive-path.v1','input_git_blob':blob,
            'scope':'width-four NN square-site; reclosed-prefix ranks, not a monotone site-insertion process',
            'D4_classes':94,'full_path_classes_h0_to3':classes3,'two_snapshot_designs':pair_data,
            'minimum_snapshot_count_for_D4_separation':2,'minimum_maximum_horizon_for_D4_separation':3,
            'minimal_horizon_witness':delay_words,
            'passive_witness':{'histories':[[0,7],[7,0]],'joint_1_2_counts':p12,'denominator_1_2':256,
                               'joint_1_3_counts':p13,'denominator_1_3':4096,
                               'joint_difference_at_general_p':'p^3*(1-p)^2*(1+p) * [[1,-1,0],[-1,1,0],[0,0,0]]',
                               'TV_joint_1_2_at_half':'3/32','shared_prediction_error_at_least':'3/64'},
            'marginal_dimensions_h0_to15':ep_dimensions,
            'marginal_space_after_rank_projection_dimension':len(conditioned.pivots),
            'two_time_dimensions_h0_to10':pd,
            'two_time_basis_events_a_b_i_j':pe,
            'two_time_integer_minor_determinant_mod_65521':pair_determinant,
            'small_marginal_mixture':small,
            'observable_dimensions_h0_to10':dimensions,
            'linear_dimension_including_constant':94,'affine_probability_dimension':93,
            'common_preparation_length':common_length,
            'basis_rank_words':[list(w) for w,v in columns],
            'integer_minor_determinant_mod_65521':determinant,'delayed_mixture':mixture,
            'infinite_markov_order_witness':memory,
            'terminal_forgetting':{'bound':'(1-p^w)^h+(1-(1-p)^w)^h-(1-p^w-(1-p)^w)^h',
                                   'uniform_over_initial_histories':True,'half_width4_one_percent_horizon':cutoff},
            'checks':{'D4_covariance_relations':8*509*16,'physical_graph_ranks':physical_count+memory_count+small_count,
                      'integer_basis':'primitive-integer Gaussian elimination over Q; no numerical rank tolerance',
                      'independent_full_rank_check':'two nonzero 94x94 integer determinants modulo 65521 (path and two-time bases)',
                      'mixture_annihilator':'93 exact zero residuals; independent forward signed-path propagation gives TV=0 through 9 and the reported TV at 10',
                      'run_formula_fraction_checks':36,'terminal_bound_horizons_checked':101},
            'boundaries':['author-level proof and executed finite certificate; no external independent review',
                          'no new width, Monte Carlo, fit, or full-repository CI',
                          '94 realization lower bound concerns the family of arbitrary initial preparations, not one fixed initial law',
                          'joint-law separation is not certain classification from one observed pair',
                          'exact linear rank is not a numerical-conditioning or statistical-efficiency guarantee',
                          'finite-order Markov exclusion is time-homogeneous and for the specified fixed preparation',
                          'run-hazard bound is local to the specified all-one conditioning, not a global approximate-Markov theorem']}


def endpoint_basis(q, rank):
    basis, columns, frontier = IntegerBasis(), [], []
    for j in range(3):
        v = [int(r == j) for r in rank]
        require(basis.add(v), 'initial marginal outputs')
        columns.append(v); frontier.append(v)
    dimensions = [len(columns)]
    while frontier:
        new = []
        for v in frontier:
            av = [sum(v[t] for t in row) for row in q]
            if basis.add(av):
                new.append(av); columns.append(av)
        dimensions.append(len(columns))
        frontier = new
    # Terminal equality of dimensions here is an A-invariant-span certificate,
    # not an extrapolation from equal successive class counts.
    return columns, dimensions


def pair_basis(q, rank, max_horizon):
    basis, events, columns = IntegerBasis(), [], []
    ep = [[[int(r == j) for j in range(3)] for r in rank]]
    for j in range(3):
        v = [int(r == j) for r in rank]
        require(basis.add(v), 'initial pair outputs')
        columns.append(v); events.append([0,0,j,j])
    dimensions = [len(columns)]
    for h in range(1, max_horizon+1):
        ep.append(apply(q,ep[-1]))
        for a in range(h):
            f = [[int(rank[s] == i)*ep[h-a][s][j] for i in range(3) for j in range(3)]
                 for s in range(len(q))]
            for _ in range(a):
                f = apply(q,f)
            for i,j in product(range(3), repeat=2):
                v = [row[3*i+j] for row in f]
                if basis.add(v):
                    columns.append(v); events.append([a,h,i,j])
        dimensions.append(len(columns))
    return columns, events, dimensions


def small_mixture(q, rank, initial, lump, reps, endpoint_columns, common):
    words = [[1,15],[1,7,13],[5,7,13],[5,15,5]]
    weights = [2,-2,-1,1]
    states = [state_for(w,q,initial) for w in words]
    # Exact A-invariant endpoint space gives all-horizon equality at p=1/2.
    require(all(sum(a*col[lump[s]] for a,s in zip(weights,states)) == 0
                for col in endpoint_columns), 'small mixture not marginal-equivalent')
    signatures = [tuple(col[lump[s]] for col in endpoint_columns) for s in states]
    require(len(set(signatures)) == 4, 'mixture collision merely repeats a pure marginal class')
    joint = two_time(q,rank,2,3)
    difference = [sum(a*joint[s][j] for a,s in zip(weights,states)) for j in range(9)]
    require(difference == [0,0,0,0,-36,36,0,36,-36], 'small mixture pair response')
    earlier = two_time(q,rank,1,3)
    require(all(sum(a*earlier[s][j] for a,s in zip(weights,states)) == 0 for j in range(9)),
            'small mixture early-pair collision')
    physical = []
    checks = 0
    for word,s in zip(words,states):
        c2 = {}
        for a,b in product(range(16),repeat=2):
            r = graph_rank(word+[a,b]); checks += 1
            require(r == rank[q[q[s][a]][b]], 'small mixture physical depth2')
            c2[a,b] = r
        hist = [0]*9
        for a,b,c in product(range(16),repeat=3):
            r = graph_rank(word+[a,b,c]); checks += 1
            require(r == rank[q[q[q[s][a]][b]][c]], 'small mixture physical depth3')
            hist[3*c2[a,b]+r] += 1
        require(hist == joint[s], 'small mixture joint histogram')
        physical.append(hist)
    tv = Fraction(sum(abs(x) for x in difference),6*16**3)
    require(tv == Fraction(3,512), 'small mixture TV')
    return checks, {'histories':words,'signed_weights':weights,'normalizer':3,
                    'common_length_histories':[common[s] for s in states],
                    'all_single_time_marginals_equal':True,'joint_1_3_equal':True,
                    'joint_2_3_difference_numerators':difference,'joint_denominator':3*16**3,
                    'joint_2_3_TV':str(tv)}


def signed_path_control(q, rank, weights, horizon):
    """Forward propagation of signed preparations: independent of column-span code."""
    current = {}
    for s,value in enumerate(weights):
        if value:
            current.setdefault((rank[s],),Counter())[s] = value
    transitions = [Counter(row) for row in q]
    normalizer = sum(x for x in weights if x > 0)
    tv_by_horizon = []
    contrasts = {}
    for h in range(horizon+1):
        contrasts = {word:sum(v.values()) for word,v in current.items() if sum(v.values())}
        tv_by_horizon.append(Fraction(sum(abs(v) for v in contrasts.values()),2*normalizer*16**h))
        if h == horizon:
            break
        new = {}
        for word,vector in current.items():
            arrival = Counter()
            for s,value in vector.items():
                for t,count in transitions[s].items():
                    arrival[t] += value*count
            for t,value in arrival.items():
                if value:
                    new.setdefault(word+(rank[t],),Counter())[t] = value
        current = new
    return tv_by_horizon, contrasts


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('certificate',type=Path)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    start=time.perf_counter()
    result=run(args.certificate)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    text = '{\n' + ',\n'.join('  '+json.dumps(k)+': '+json.dumps(v,ensure_ascii=False,separators=(',',':')) for k,v in result.items()) + '\n}\n'
    args.out.write_text(text,encoding='utf-8')
    print(json.dumps({'output':str(args.out),'elapsed_seconds':time.perf_counter()-start,'python':sys.version.split()[0]}))

if __name__=='__main__':
    main()
