#!/usr/bin/env python3
"""Exact #809 endpoint-reversal controls on the unmodified PR708 certificate.
Python 3.10+, standard library only. No sampling, new width, or transfer build.
Usage: python check.py /path/to/width4-rank-closure-certificate.json --out result.json
"""
from __future__ import annotations
import argparse
from collections import Counter, deque
from fractions import Fraction
from hashlib import sha1, sha256
from itertools import product
import json
from math import comb
from pathlib import Path
import sys
import time

INPUT_BLOB = '50b7297deefe7c50215aea2ed534ca5810461af3'
INPUT_SHA256 = '508640f6a6be1461ca4738acdcf3226790f485dbaa3c43b8a00966c143f3f890'


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def labels(keys) -> list[int]:
    index, result = {}, []
    for key in keys:
        key = tuple(key)
        if key not in index:
            index[key] = len(index)
        result.append(index[key])
    return result


def partition_equal(a, b) -> bool:
    return labels((x,) for x in a) == labels((x,) for x in b)


def physical_histories(q, initial) -> list[list[int]]:
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
    require(len(words) == len(q), 'unreachable state at physical length >=2')
    return [words[s] for s in range(len(q))]


def state_for(word, q, initial) -> int:
    s = initial[word[0]]
    for b in word[1:]:
        s = q[s][b]
    return s


def deterministic_refinement(base, q):
    old = list(base)
    counts = [len(set(old))]
    while True:
        new = labels((old[s], *(old[t] for t in row)) for s, row in enumerate(q))
        counts.append(len(set(new)))
        if partition_equal(old, new):
            return new, counts
        old = new


def lump_refinement(base, q, all_p=False):
    old = list(base)
    counts = [len(set(old))]
    while True:
        keys = []
        for s, row in enumerate(q):
            rates = Counter((b.bit_count() if all_p else 0, old[t])
                            for b, t in enumerate(row))
            keys.append((old[s], tuple(sorted(rates.items()))))
        new = labels(keys)
        counts.append(len(set(new)))
        if partition_equal(old, new):
            return new, counts
        old = new


def suffix_signatures(q, rank, mode, maximum=2):
    """Exact integer Bernstein counts, not samples of a probability grid.

    homogeneous: one common p for every future site;
    row_schedule: one independent p_t per future row;
    stationary_columns: four column probabilities repeated at every row.
    """
    n = len(q)
    key0 = (0,) if mode == 'homogeneous' else ((0, 0, 0, 0) if mode == 'stationary_columns' else ())
    coefficients = {key0: [[int(rank[s] == j) for j in range(3)] for s in range(n)]}
    cumulative = [0] * n
    table = []
    for depth in range(maximum + 1):
        keys = sorted(coefficients)
        signatures = [tuple(v for key in keys for v in coefficients[key][s]) for s in range(n)]
        single = labels(signatures)
        cumulative = labels((cumulative[s], *signatures[s]) for s in range(n))
        table.append({'horizon': depth, 'single': len(set(single)), 'through_horizon': len(set(cumulative))})
        if depth == maximum:
            break
        new = {}
        for key, values in coefficients.items():
            for b in range(16):
                if mode == 'homogeneous':
                    new_key = (key[0] + b.bit_count(),)
                elif mode == 'row_schedule':
                    new_key = (b.bit_count(),) + key
                elif mode == 'stationary_columns':
                    new_key = tuple(key[x] + ((b >> x) & 1) for x in range(4))
                else:
                    raise ValueError(mode)
                target = new.setdefault(new_key, [[0, 0, 0] for _ in range(n)])
                for s in range(n):
                    value = values[q[s][b]]
                    for j in range(3):
                        target[s][j] += value[j]
        coefficients = new
    return cumulative, table


def half_signatures(q, rank, maximum=3):
    v = [[int(r == j) for j in range(3)] for r in rank]
    cumulative = [0] * len(q)
    table = []
    for h in range(maximum + 1):
        cumulative = labels((cumulative[s], *v[s]) for s in range(len(q)))
        table.append({'horizon': h, 'single': len(set(labels(v))), 'through_horizon': len(set(cumulative))})
        v = [[sum(v[t][j] for t in row) for j in range(3)] for row in q]
    return cumulative, table


def graph_rank(rows, width=4):
    """Physical lifted-graph DFS, independent of certificate transitions."""
    height, positions, first = len(rows), {}, None
    for y0 in range(height):
        for x0 in range(width):
            root = y0 * width + x0
            if not (rows[y0] >> x0 & 1) or root in positions:
                continue
            positions[root] = (0, 0)
            stack = [root]
            while stack:
                u = stack.pop()
                x, y = u % width, u // width
                px, py = positions[u]
                for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    nx, ny = (x + dx) % width, (y + dy) % height
                    if not (rows[ny] >> nx & 1):
                        continue
                    v, proposed = ny * width + nx, (px + dx, py + dy)
                    if v not in positions:
                        positions[v] = proposed
                        stack.append(v)
                    else:
                        a, b = proposed[0] - positions[v][0], proposed[1] - positions[v][1]
                        require(a % width == 0 and b % height == 0, 'nonperiodic cycle')
                        if a or b:
                            if first is None:
                                first = (a, b)
                            elif first[0] * b - first[1] * a:
                                return 2
    return int(first is not None)


# Small exact bivariate power-polynomial algebra, (x degree, y degree) -> integer.
def padd(*polys):
    out = Counter()
    for poly in polys:
        for key, value in poly.items():
            out[key] += value
    return {key: value for key, value in out.items() if value}


def pmul(a, b):
    out = Counter()
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            out[i + k, j + l] += x * y
    return {key: value for key, value in out.items() if value}


def ppow(a, k):
    out = {(0, 0): 1}
    for _ in range(k):
        out = pmul(out, a)
    return out


def pscale(a, k):
    return {key: k * v for key, v in a.items() if k * v}


def bernstein2(hist):
    out = Counter()
    for k in range(5):
        for l in range(5):
            for a in range(5 - k):
                for b in range(5 - l):
                    out[k + a, l + b] += hist[k][l] * comb(4 - k, a) * comb(4 - l, b) * (-1) ** (a + b)
    return {key: value for key, value in out.items() if value}


def dipole_polynomial(poly):
    # Substitute x=p+d and y=p-d; the returned coordinates are (p degree,d degree).
    plus, minus = {(1, 0): 1, (0, 1): 1}, {(1, 0): 1, (0, 1): -1}
    return padd(*(pscale(pmul(ppow(plus, i), ppow(minus, j)), v) for (i, j), v in poly.items()))


def quadratic_remainder(coefficients):
    # Divide by p^2+2p-1, monic, over Z[p].
    a = dict(coefficients)
    while a and max(a) >= 2:
        k = max(a)
        v = a.pop(k)
        a[k - 1] = a.get(k - 1, 0) - 2 * v
        a[k - 2] = a.get(k - 2, 0) + v
        a = {j: x for j, x in a.items() if x}
    return a


def pair_control(histories, q, rank, initial):
    counts, laws = [], []
    checks = 0
    for history in histories:
        s = state_for(history, q, initial)
        hist = [[[0] * 5 for _ in range(5)] for _ in range(3)]
        law = [Fraction(0) for _ in range(3)]
        for a, b in product(range(16), repeat=2):
            output = graph_rank(history + [a, b])
            require(output == rank[q[q[s][a]][b]], 'physical graph/certificate mismatch')
            checks += 1
            hist[output][a.bit_count()][b.bit_count()] += 1
            law[output] += Fraction(2 ** (4 - a.bit_count() + b.bit_count()), 3 ** 8)
        require(sum(law) == 1, 'probability normalization')
        counts.append(hist)
        laws.append(law)
    differences = []
    for j in range(3):
        differences.append(bernstein2([[counts[0][j][k][l] - counts[1][j][k][l]
                                       for l in range(5)] for k in range(5)]))
    return differences, laws, checks


def extended_ch(q, rank, lump, expected):
    """Independent all-horizon upper check: powers 0..93 over Q[p]."""
    reps = [lump.index(i) for i in range(len(set(lump)))]
    grouped = [Counter((b.bit_count(), lump[t]) for b, t in enumerate(q[s])) for s in reps]
    rank = [rank[s] for s in reps]
    expected = [expected[s] for s in reps]
    values = [[[int(rank[s] == j)] for j in range(3)] for s in range(len(reps))]
    first = {}
    for s, a in enumerate(expected):
        first.setdefault(a, s)
    for h in range(len(reps)):
        require(all(values[s] == values[first[a]] for s, a in enumerate(expected)), 'CH moment difference')
        if h + 1 == len(reps):
            break
        size = 4 * h + 1
        new = [[[0] * (size + 4) for _ in range(3)] for _ in reps]
        for s, edges in enumerate(grouped):
            for (k, t), weight in edges.items():
                for j in range(3):
                    src, dst = values[t][j], new[s][j]
                    for i, value in enumerate(src):
                        dst[k + i] += weight * value
        values = new
    return len(reps)


def run(path: Path, ch: bool):
    raw = path.read_bytes()
    blob = sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
    require(blob == INPUT_BLOB and sha256(raw).hexdigest() == INPUT_SHA256, 'wrong or corrupted input certificate')
    cert = json.loads(raw)
    q, rank, initial = cert['quotient_transitions'], cert['quotient_rank_output'], cert['quotient_initial']
    n = len(q)
    require(n == 509 and all(len(row) == 16 for row in q), 'wrong machine shape')
    words = physical_histories(q, initial)
    rev = [state_for(word[::-1], q, initial) for word in words]
    require(all(rev[rev[s]] == s and rank[rev[s]] == rank[s] for s in range(n)), 'reversal involution')
    gwords = physical_histories(cert['transitions'], cert['initial'])
    require(all(state_for(word[::-1], q, initial) == rev[cert['quotient_map'][s]]
                for s, word in enumerate(gwords)), 'reversal differs on geometric representatives')
    actions = []
    for sign, shift in product((1, -1), range(4)):
        maps = [sum(1 << ((sign * x + shift) % 4) for x in range(4) if b >> x & 1) for b in range(16)]
        perm = [state_for([maps[b] for b in word], q, initial) for word in words]
        require(sorted(perm) == list(range(n)), 'D4 not bijective')
        require(all(perm[rev[s]] == rev[perm[s]] and rank[perm[s]] == rank[s] for s in range(n)), 'D4/reversal output')
        require(all(perm[q[s][b]] == q[perm[s]][maps[b]] for s in range(n) for b in range(16)), 'D4 transition covariance')
        actions.append(perm)
    d4 = labels((min(perm[s] for perm in actions),) for s in range(n))
    reversal = labels((min(s, rev[s]),) for s in range(n))
    both = labels((min([perm[s] for perm in actions] + [rev[perm[s]] for perm in actions]),) for s in range(n))
    require([len(set(v)) for v in (d4, reversal, both)] == [94, 274, 62], 'orbit counts')
    d4_reps = [d4.index(i) for i in range(94)]
    reversed_d4 = [d4[rev[s]] for s in d4_reps]
    fixed_d4 = sum(i == j for i, j in enumerate(reversed_d4))
    require(fixed_d4 == 30, 'self-reversing D4 orbit count')
    homogeneous, ht = suffix_signatures(q, rank, 'homogeneous')
    scheduled, st = suffix_signatures(q, rank, 'row_schedule')
    columns, ct = suffix_signatures(q, rank, 'stationary_columns')
    half, ft = half_signatures(q, rank)
    require(partition_equal(homogeneous, both) and partition_equal(half, both), 'homogeneous orbit separation')
    require(partition_equal(scheduled, d4) and partition_equal(columns, reversal), 'source-family orbit separation')
    det, dt = deterministic_refinement(rank, q)
    rec, rt = deterministic_refinement(both, q)
    fixed_lump, lt = lump_refinement(rank, q)
    joint_lump, jt = lump_refinement(rank, q, all_p=True)
    require(len(set(det)) == len(set(rec)) == 509, 'recursive closure')
    require(partition_equal(fixed_lump, d4) and partition_equal(joint_lump, d4), 'strong lumping')
    # First antisymmetric two-row pulse: [K_prime,K] C at p=1/2, common denominator128.
    C = [[int(r == j) for j in range(3)] for r in rank]
    A = [[sum(C[t][j] for t in row) for j in range(3)] for row in q]
    D = [[sum((2 * b.bit_count() - 4) * C[t][j] for b, t in enumerate(row)) for j in range(3)] for row in q]
    pulse = [[sum((2 * b.bit_count() - 4) * A[t][j] - D[t][j] for b, t in enumerate(row)) for j in range(3)] for row in q]
    pulse_labels = labels((both[s], *pulse[s]) for s in range(n))
    require(partition_equal(pulse_labels, d4), 'first pulse fails to recover D4 quotient')
    require(all(pulse[rev[s]] == [-x for x in pulse[s]] for s in range(n)), 'pulse not reversal-odd')
    x, y, one = {(1, 0): 1}, {(0, 1): 1}, {(0, 0): 1}
    ux, uy, xy = padd(one, pscale(x, -1)), padd(one, pscale(y, -1)), pmul(x, y)
    diffxy = padd(x, pscale(y, -1))
    polynomial = padd(pmul(ppow(x, 2), y), ppow(x, 2), pmul(x, ppow(y, 2)), xy, ppow(y, 2))
    expected_simple = pscale(pmul(pmul(pmul(ux, uy), diffxy), polynomial), -1)
    simple, laws, graph_checks = pair_control([[0, 7], [7, 0]], q, rank, initial)
    require(simple == [expected_simple, pscale(expected_simple, -1), {}], 'simple witness polynomial')
    tv = sum(abs(a - b) for a, b in zip(*laws)) / 2
    require(tv == Fraction(2, 27), 'simple witness TV')
    exceptional, _, extra_checks = pair_control([[3, 15, 5], [5, 15, 3]], q, rank, initial)
    expected_exception = pmul(pmul(pmul(pmul(xy, ux), uy), diffxy), padd(xy, x, y, pscale(one, -1)))
    require(exceptional == [{}, expected_exception, pscale(expected_exception, -1)], 'exception polynomial')
    dipole = dipole_polynomial(expected_exception)
    linear = {i: v for (i, j), v in dipole.items() if j == 1}
    cubic = {i: v for (i, j), v in dipole.items() if j == 3}
    require(quadratic_remainder(linear) == {}, 'linear root exception')
    require(quadratic_remainder(cubic) == {0: -20, 1: 48}, 'cubic response exception')
    # A same-input recursive counterexample, with one future uniform row.
    s, t = [state_for(word, q, initial) for word in ([0, 7], [0, 11])]
    require(both[s] == both[t], 'recursive witness not initially equivalent')
    delta = [A[q[s][5]][j] - A[q[t][5]][j] for j in range(3)]
    require(delta == [-1, 1, 0], 'recursive witness contrast /16')
    for word in ([0, 7, 5], [0, 11, 5]):
        for b in range(16):
            require(graph_rank(list(word) + [b]) == rank[q[state_for(word, q, initial)][b]], 'recursive physical witness')
    s, t = [state_for(word, q, initial) for word in ([0, 7], [7, 0])]
    pa, pb = Counter(both[u] for u in q[s]), Counter(both[u] for u in q[t])
    fail = {str(k): pa[k] - pb[k] for k in sorted(set(pa) | set(pb)) if pa[k] != pb[k]}
    require(bool(fail), '62 quotient unexpectedly lumpable')
    return {
        'schema': 'matching-one.p809-endpoint-reversal.v1',
        'input_git_blob': blob, 'input_sha256': INPUT_SHA256,
        'deterministic_states': n, 'D4_orbits': 94, 'reversal_orbits': 274, 'D4_reversal_orbits': 62,
        'reversal_fixed_states': sum(rev[s] == s for s in range(n)),
        'D4_orbit_reversal_pairs': (94 - fixed_d4) // 2, 'D4_orbit_reversal_fixed': fixed_d4,
        'homogeneous_all_p': ht, 'row_scheduled_p': st, 'stationary_column_sources': ct, 'fixed_p_half': ft,
        'deterministic_refinement': dt, 'recursive_closure_of_62': rt,
        'strong_lumping_half': lt, 'strong_lumping_all_p': jt,
        'first_antisymmetric_pulse_at_half_classes': len(set(pulse_labels)),
        'simple_witness': {'histories': [[0, 7], [7, 0]], 'future_probabilities': ['1/3', '2/3'],
                           'rank_laws': [[str(z) for z in law] for law in laws], 'TV': str(tv),
                           'unavoidable_shared_prediction_TV_error_at_least': str(tv / 2)},
        'recursive_witness': {'histories': [[0, 7], [0, 11]], 'common_appended_mask': 5,
                              'next_uniform_row_law_difference': ['-1/16', '1/16', '0']},
        'not_lumpable_witness': {'histories': [[0, 7], [7, 0]], 'class_probability_difference_numerators': fail, 'denominator': 16},
        'linear_exception': {'histories': [[3, 15, 5], [5, 15, 3]], 'p': 'sqrt(2)-1',
                             'rank1_linear_coefficient': '2*p^2*(p-1)^2*(p^2+2*p-1)',
                             'rank1_cubic_coefficient_at_exception': '48*p-20 = 48*sqrt(2)-68'},
        'checks': {'D4_transition_covariance': 8 * n * 16,
                   'reversal_geometric_representatives': len(gwords),
                   'physical_graph_configurations': graph_checks + extra_checks + 32,
                   'witness_polynomial_identities': 'exact coefficientwise',
                   'optional_CH_powers': extended_ch(q, rank, joint_lump, both) if ch else None},
        'boundaries': ['author proof plus exact finite controls, not independent review',
                       'no width above four; no Monte Carlo; no transfer-matrix production',
                       'no full-repository CI; no field/Jordan/critical-exponent identification',
                       'class counts are not minimal linear or positive realization dimensions'],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--extended-ch', action='store_true')
    args = parser.parse_args()
    start = time.perf_counter()
    result = run(args.certificate, args.extended_ch)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps({'output': str(args.out), 'elapsed_seconds': time.perf_counter() - start,
                      'python': sys.version.split()[0], 'result_sha256': sha256(args.out.read_bytes()).hexdigest()}))


if __name__ == '__main__':
    main()
