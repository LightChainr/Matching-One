#!/usr/bin/env python3
"""Exact second-order observability on the existing P398 process, widths 4/5.

Uses only the Python standard library. This is an independent implementation of
P398's declared noncrossing join/detach model and original three readouts/four
sources, not a new percolation model. All claim-bearing arithmetic is integer
or Fraction. A nonzero finite-field minor certifies a rational rank lower bound;
parity/ambient dimension supplies the matching upper bound.

H = join@0 - join@(w-2) = 2 * the repository's H_odd. G+epsilon*H is a
Markov generator for |epsilon|<=1 (irreducible for |epsilon|<1 here).
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
from typing import Sequence

PRIME = 1_000_000_007
Matrix = list[list[int]]


def canonical(values: Sequence[int]) -> tuple[int, ...]:
    names: dict[int, int] = {}
    return tuple(names.setdefault(x, len(names)) for x in values)


def noncrossing_states(width: int) -> list[tuple[int, ...]]:
    if width not in (4, 5):
        raise ValueError('This bounded analysis is declared only at widths 4 and 5.')
    quadruples = list(combinations(range(width), 4))
    def grow(prefix):
        if len(prefix) == width:
            if not any(prefix[a] == prefix[c] and prefix[b] == prefix[d]
                       and prefix[a] != prefix[b] for a, b, c, d in quadruples):
                yield prefix
            return
        for label in range(max(prefix, default=-1) + 2):
            yield from grow(prefix + (label,))
    return list(grow(()))


def zeros(n: int, m: int) -> Matrix:
    return [[0] * m for _ in range(n)]


def transpose(a: Matrix) -> Matrix:
    return [list(c) for c in zip(*a)]


def matvec(a: Matrix, x: Sequence[int]) -> list[int]:
    return [sum(t * z for t, z in zip(row, x)) for row in a]


def multiply(a: Matrix, b: Matrix) -> Matrix:
    bt = transpose(b)
    return [[sum(x * y for x, y in zip(row, col)) for col in bt] for row in a]


def build_model(width: int) -> dict:
    states = noncrossing_states(width)
    index = {state: i for i, state in enumerate(states)}
    n = len(states)
    joins, detaches = [], []
    for point in range(width):
        j, d = zeros(n, n), zeros(n, n)
        for i, state in enumerate(states):
            image = canonical([state[point] if x == state[(point + 1) % width] else x
                               for x in state])
            j[i][index[image]] += 1
            j[i][i] -= 1
            q = list(state)
            q[point] = width
            d[i][index[canonical(q)]] += 1
            d[i][i] -= 1
        joins.append(j)
        detaches.append(d)
    g = [[sum(a[i][j] for a in joins + detaches) for j in range(n)] for i in range(n)]
    h = [[joins[0][i][j] - joins[width - 2][i][j] for j in range(n)] for i in range(n)]
    reflect = [index[canonical(state[::-1])] for state in states]
    f = [[max(state) + 1, sum(state.count(x) == 1 for x in set(state)),
          int(state[0] == state[-1])] for state in states]
    # n times each ORIGINAL probability source: no floating-point uniform law.
    src = zeros(4, n)
    initial = [tuple(range(width)), (0,) * width,
               tuple([0] + list(range(1, width - 1)) + [0])]
    for k, state in enumerate(initial):
        src[k][index[state]] = n
    src[3] = [1] * n
    for i in range(n):
        assert sum(g[i]) == sum(h[i]) == 0
        assert f[i] == f[reflect[i]]
        for j in range(n):
            assert g[reflect[i]][reflect[j]] == g[i][j]
            assert h[reflect[i]][reflect[j]] == -h[i][j]
            if i != j:
                assert g[i][j] >= abs(h[i][j])
    assert all(row[i] == row[reflect[i]] for row in src for i in range(n))
    for matrix in (g, transpose(g)):
        reached, queue = {0}, [0]
        for i in queue:
            for j, rate in enumerate(matrix[i]):
                if j != i and rate > 0 and j not in reached:
                    reached.add(j)
                    queue.append(j)
        assert len(reached) == n
    return dict(width=width, states=states, G=g, H=h, reflection=reflect,
                F=f, sources_scaled=src, source_denominator=n)


def parity_basis(model: dict, parity: int) -> tuple[Matrix, list[int]]:
    pi = model['reflection']
    representatives = [i for i, j in enumerate(pi) if i < j or (parity == 1 and i == j)]
    u = zeros(len(pi), len(representatives))
    for k, i in enumerate(representatives):
        u[i][k] = 1
        if pi[i] != i:
            u[pi[i]][k] = parity
    return u, representatives


class ModularBasis:
    def __init__(self, prime: int = PRIME):
        self.rows: list[tuple[int, list[int]]] = []
        self.prime = prime

    def add(self, vector: Sequence[int]) -> bool:
        p = self.prime
        v = [x % p for x in vector]
        for pivot, row in self.rows:
            factor = v[pivot]
            if factor:
                v = [(x - factor * y) % p for x, y in zip(v, row)]
        if not any(v):
            return False
        pivot = next(i for i, x in enumerate(v) if x)
        inverse = pow(v[pivot], p - 2, p)
        self.rows.append((pivot, [x * inverse % p for x in v]))
        return True


def determinant_mod(a: Matrix, prime: int = PRIME) -> int:
    n = len(a)
    if any(len(row) != n for row in a):
        raise ValueError('determinant expects a square matrix')
    b, result = [[x % prime for x in row] for row in a], 1
    for k in range(n):
        pivot = next((i for i in range(k, n) if b[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            b[k], b[pivot] = b[pivot], b[k]
            result = -result
        value = b[k][k]
        result = result * value % prime
        inverse = pow(value, prime - 2, prime)
        for i in range(k + 1, n):
            factor = b[i][k] * inverse % prime
            for j in range(k + 1, n):
                b[i][j] = (b[i][j] - factor * b[k][j]) % prime
    return result % prime


def word_span(operators: list[tuple[str, Matrix]], seeds: Matrix, side: str) -> dict:
    """Exact integer vectors selected by independent modular pivots.

    At a certified full ambient bound these span over Q as well. When the
    observed rank is below a proven bound, this routine reports only a lower
    bound; a modular stabilization alone is NOT a rational upper-bound proof.
    """
    basis, selected, frontier = ModularBasis(), [], []
    for j, vector in enumerate(transpose(seeds)):
        if basis.add(vector):
            item = dict(seed=j, word='', vector=vector)
            selected.append(item)
            frontier.append(item)
    levels = [len(selected)]
    while frontier:
        new = []
        for name, a in operators:
            for item in frontier:
                vector = matvec(a, item['vector'])
                if basis.add(vector):
                    word = name + item['word'] if side == 'column' else item['word'] + name
                    nxt = dict(seed=item['seed'], word=word, vector=vector)
                    selected.append(nxt)
                    new.append(nxt)
        levels.append(len(selected))
        frontier = new
    return dict(selected=selected, rank_lower_bound=len(selected), levels=levels)


def realization_certificate(operators: list[tuple[str, Matrix]], f: Matrix,
                            src: Matrix, upper_bound: int) -> dict:
    r = word_span(operators, f, 'column')
    o = word_span([(name, transpose(a)) for name, a in operators], transpose(src), 'row')
    gram = [[sum(x * y for x, y in zip(row['vector'], col['vector']))
             for col in r['selected']] for row in o['selected']]
    basis, ri = ModularBasis(), []
    for i, row in enumerate(gram):
        if basis.add(row):
            ri.append(i)
    basis, ci = ModularBasis(), []
    for j, col in enumerate(transpose([gram[i] for i in ri])):
        if basis.add(col):
            ci.append(j)
    minor = [[gram[i][j] for j in ci] for i in ri]
    det = determinant_mod(minor)
    rank = len(ri)
    if rank != upper_bound or not det:
        raise ValueError(f'Lower bound {rank} does not attain asserted upper bound {upper_bound}')
    return dict(minimal_order=rank, ambient_upper_bound=upper_bound,
                reach_levels=r['levels'], observe_levels=o['levels'], prime=PRIME,
                hankel_minor_determinant_mod_prime=det,
                row_words=[{k: v for k, v in o['selected'][i].items() if k != 'vector'} for i in ri],
                column_words=[{k: v for k, v in r['selected'][j].items() if k != 'vector'} for j in ci],
                integer_hankel_minor=minor)


def charpoly_integer(a: Matrix) -> list[int]:
    """Faddeev-LeVerrier, integer division checked; coefficients high first."""
    n = len(a)
    b = [[int(i == j) for j in range(n)] for i in range(n)]
    coefficients = [1]
    for k in range(1, n + 1):
        ab = multiply(a, b)
        trace = sum(ab[i][i] for i in range(n))
        assert trace % k == 0
        c = -trace // k
        coefficients.append(c)
        b = [[ab[i][j] + (c if i == j else 0) for j in range(n)] for i in range(n)]
    assert not any(x for row in b for x in row)
    return coefficients


def pulse_data(model: dict) -> dict:
    g, h, f, src = (model[k] for k in ('G', 'H', 'F', 'sources_scaled'))
    u, reps = parity_basis(model, -1)
    a = [row for i, row in enumerate(multiply(g, u)) if i in reps]
    b = [row for i, row in enumerate(multiply(h, f)) if i in reps]
    c = multiply(multiply(src, h), u)
    assert multiply(g, u) == multiply(u, a)
    assert multiply(h, f) == multiply(u, b)
    d = len(reps)
    certificate = realization_certificate([('G', a)], b, c, d)
    current, moments = b, []
    for _ in range(2 * d + 1):
        moments.append(multiply(c, current))
        current = multiply(a, current)
    # Scalar witnesses chosen from the ORIGINAL dictionary. No new readout.
    source, output = (2, 2) if model['width'] == 4 else (0, 2)
    scalar = [Fraction(m[source][output], model['source_denominator']) for m in moments]
    assert all(x.denominator == 1 for x in scalar)
    scalar = [int(x) for x in scalar]
    hankel = [[scalar[i + j] for j in range(d)] for i in range(d)]
    determinant = determinant_mod(hankel)
    assert determinant != 0
    polynomial = charpoly_integer(a)
    for k in range(len(scalar) - d):
        assert sum(polynomial[j] * scalar[k + d - j] for j in range(d + 1)) == 0
    # Numerator of C(zI-A)^(-1)B from the first d moments, descending order.
    numerator = [sum(polynomial[j] * scalar[k-j] for j in range(k + 1)) for k in range(d)]
    while len(numerator) > 1 and numerator[0] == 0:
        numerator.pop(0)
    return dict(odd_dimension=d, G_odd=a, B_odd=b, C_odd_scaled=c,
                source_denominator=model['source_denominator'],
                odd_characteristic_polynomial_high_first=polynomial,
                rank_certificate=certificate, markov_moments_scaled=moments,
                scalar_witness=dict(source_index=source, readout_index=output,
                                    moments=scalar, hankel_minor=hankel,
                                    hankel_determinant_mod_prime=determinant,
                                    laplace_numerator_high_first=numerator,
                                    laplace_denominator_high_first=polynomial))


def analyse(width: int) -> dict:
    model = build_model(width)
    n = len(model['states'])
    g, h, f, src = (model[k] for k in ('G', 'H', 'F', 'sources_scaled'))
    up, reps = parity_basis(model, 1)
    gp = [row for i, row in enumerate(multiply(g, up)) if i in reps]
    fp = [f[i] for i in reps]
    sp = multiply(src, up)
    assert multiply(up, fp) == f
    assert multiply(g, up) == multiply(up, gp)
    baseline = realization_certificate([('G', gp)], fp, sp, len(reps))
    controlled = realization_certificate([('G', g), ('H', h)], f, src, n)
    contrast = [[src[i][j] - src[3][j] for j in range(n)] for i in range(3)]
    # Constants are invariant under G and killed by H and the source contrasts.
    # Therefore n-1 is a mathematical upper bound, not a numerical rank guess.
    contrasted = realization_certificate([('G', g), ('H', h)], f, contrast, n - 1)
    return dict(width=width, states=model['states'], reflection=model['reflection'],
                G=g, H=h, readouts=f, sources_scaled=src, source_denominator=n,
                baseline=baseline, double_pulse=pulse_data(model),
                arbitrary_controlled_word=controlled, controlled_source_contrasts=contrasted)


def report() -> dict:
    return dict(schema='matching-one.p398-double-pulse-exact.v1',
                model='P398 noncrossing join/detach calibration process; NOT percolation',
                source_ref='eb89e9422791d9e3c3a78f0e65d56912b815a7bd',
                source_files={'scripts/p398_intervention_transport.py':'ac13194751b8ba727a45499900e3f12ce1b90461',
                              'scripts/planar_state_operations.py':'c5f57dcb8606dbb62306cde839edefae680f863d'},
                H_convention='join@0 - join@(w-2); twice the archived H_odd',
                source_order=['delta_all_singletons','delta_single_block','delta_wrapped_pair','uniform'],
                readout_order=['blocks','singletons','wrap'],
                pulse_kernel='K(tau)=S H exp(tau G) H F; S is the probability-source matrix',
                claim_boundary=['Exact homogeneous linear/bilinear realization orders include the constant mode.',
                                'A finite-field nonzero minor is a rational lower bound; parity and dimension give the upper bound.',
                                'No new sampling, asymptotic width law, positivity-minimality, CFT or threshold claim.',
                                'The short-pulse kernel is a limit of derivatives of physical nonnegative-rate windows, not exp(epsilon H).'],
                widths=[analyse(4), analyse(5)])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    text = json.dumps(report(), indent=2, allow_nan=False) + '\n'
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open('x', encoding='utf-8') as handle:
            handle.write(text)
    else:
        print(text, end='')
