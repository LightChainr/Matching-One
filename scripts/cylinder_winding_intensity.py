#!/usr/bin/env python3
"""Exact once-per-component winding rewards on fixed-width site cylinders.

The production transfer retains one frontier, integer horizontal cut gains,
and a winding flag per active component.  It is NOT the torus-rank automaton.
An independent graph-potential BFS supplies finite controls.

Normal report: Python stdlib + SymPy (only for symbolic rational functions).
No simulation. Width 2--4 is exhaustively controlled. Larger widths are a
bounded capacity probe, not an asserted cheap asymptotic computation.
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from math import comb, log, pi, sqrt
from pathlib import Path
from time import perf_counter
from typing import NamedTuple


class State(NamedTuple):
    labels: tuple[int, ...]
    gains: tuple[int, ...]
    winding: tuple[int, ...]


class GainDSU:
    """Potential(b)-potential(a)=gain on each oriented edge a->b."""
    def __init__(self, n: int) -> None:
        self.parent = list(range(n))
        self.delta = [0] * n
        self.winding = [False] * n

    def find(self, a: int) -> tuple[int, int]:
        if self.parent[a] != a:
            root, shift = self.find(self.parent[a])
            self.delta[a] += shift
            self.parent[a] = root
        return self.parent[a], self.delta[a]

    def join(self, a: int, b: int, gain: int) -> None:
        ra, da = self.find(a)
        rb, db = self.find(b)
        if ra == rb:
            self.winding[ra] |= db - da != gain
        else:
            self.parent[rb] = ra
            self.delta[rb] = gain + da - db
            self.winding[ra] |= self.winding[rb]


def empty_state(width: int) -> State:
    if width < 2:
        raise ValueError('width must be at least two; retain lifted parallel edges')
    return State((-1,) * width, (0,) * width, ())


def advance(state: State, mask: int, matching: bool = False) -> tuple[State, int]:
    """Append one spatial row; reward every retired winding component once."""
    width = len(state.labels)
    if not 0 <= mask < 1 << width:
        raise ValueError('row mask outside width')
    dsu = GainDSU(2 * width)
    old = [i for i, k in enumerate(state.labels) if k >= 0]
    new = [i for i in range(width) if mask >> i & 1]
    representatives: dict[int, int] = {}
    for i in old:
        k = state.labels[i]
        if k in representatives:
            dsu.join(representatives[k], i, state.gains[i])
        else:
            representatives[k] = i
    for k, i in representatives.items():
        dsu.winding[dsu.find(i)[0]] = bool(state.winding[k])
    for i in new:
        j = (i + 1) % width
        if mask >> j & 1:
            dsu.join(width + i, width + j, (i + 1) // width)
        for dx in ((-1, 0, 1) if matching else (0,)):
            j = (i + dx) % width
            if state.labels[j] >= 0:
                dsu.join(width + i, j, (i + dx) // width)
    all_roots = {dsu.find(i)[0] for i in old + [width + i for i in new]}
    kept_roots = {dsu.find(width + i)[0] for i in new}
    reward = sum(dsu.winding[r] for r in all_roots - kept_roots)
    labels = [-1] * width
    gains = [0] * width
    flags: list[int] = []
    root_map: dict[int, tuple[int, int]] = {}
    for i in new:
        root, potential = dsu.find(width + i)
        if root not in root_map:
            root_map[root] = (len(flags), potential)
            flags.append(int(dsu.winding[root]))
        k, origin = root_map[root]
        labels[i] = k
        gains[i] = 0 if dsu.winding[root] else potential - origin
    return State(tuple(labels), tuple(gains), tuple(flags)), int(reward)


def build_transfer(width: int, matching: bool = False,
                   state_cap: int = 5000) -> tuple[list[State], list[list[tuple[int, int]]]]:
    """BFS exhausts every row successor; fail explicitly at the chosen cap."""
    if width > 10:
        raise ValueError('reference Python builder limited to width <=10')
    states = [empty_state(width)]
    index = {states[0]: 0}
    transfer: list[list[tuple[int, int]]] = []
    for state in states:
        row = []
        for mask in range(1 << width):
            nxt, reward = advance(state, mask, matching)
            if nxt not in index:
                if len(states) >= state_cap:
                    raise RuntimeError('state cap reached; no incomplete closure returned')
                index[nxt] = len(states)
                states.append(nxt)
            row.append((index[nxt], reward))
        transfer.append(row)
    return states, transfer


def reward_lump(transfer: list[list[tuple[int, int]]]) -> tuple[list[list[tuple[int, int]]], list[int]]:
    """Common all-p stochastic lumping retaining the joint next-class/reward law.

    This is not a claimed minimal positive/linear realization. Coefficients are
    grouped by the new row's number of occupied sites, not numeric p samples.
    """
    blocks = [0] * len(transfer)
    while True:
        classes: dict[tuple, int] = {}
        refined = []
        for row in transfer:
            counts = Counter((mask.bit_count(), reward, blocks[j])
                             for mask, (j, reward) in enumerate(row))
            signature = tuple(sorted(counts.items()))
            if signature not in classes:
                classes[signature] = len(classes)
            refined.append(classes[signature])
        if refined == blocks:
            break
        blocks = refined
    reduced = [[(blocks[j], reward) for j, reward in transfer[blocks.index(k)]]
               for k in range(max(blocks) + 1)]
    return reduced, blocks


def row_weights(width: int, p: Fraction) -> list[Fraction]:
    if not 0 <= p <= 1:
        raise ValueError('p outside [0,1]')
    return [p ** mask.bit_count() * (1 - p) ** (width - mask.bit_count())
            for mask in range(1 << width)]


def solve_fraction(a: list[list[Fraction]], b: list[Fraction]) -> list[Fraction]:
    """Exact Gauss-Jordan solve with explicit singularity checking."""
    n = len(b)
    if len(a) != n or any(len(row) != n for row in a):
        raise ValueError('matrix must be square')
    mat = [[Fraction(x) for x in row] + [Fraction(rhs)] for row, rhs in zip(a, b)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if mat[i][j]), None)
        if pivot is None:
            raise ValueError('singular exact linear system')
        mat[j], mat[pivot] = mat[pivot], mat[j]
        v = mat[j][j]
        mat[j] = [x / v for x in mat[j]]
        for i in range(n):
            if i != j and mat[i][j]:
                v = mat[i][j]
                mat[i] = [x - v * y for x, y in zip(mat[i], mat[j])]
    return [row[-1] for row in mat]


def stationary_reward(transfer: list[list[tuple[int, int]]], p: Fraction) -> dict:
    if not 0 < p < 1:
        raise ValueError('stationary pressure calculation requires 0<p<1')
    n = len(transfer)
    width = (len(transfer[0]) - 1).bit_length()
    weights = row_weights(width, p)
    matrices = [[[Fraction(0) for _ in range(n)] for _ in range(n)] for _ in range(3)]
    for i, row in enumerate(transfer):
        for probability, (j, reward) in zip(weights, row):
            for k in range(3):
                matrices[k][i][j] += probability * reward ** k
    k0, k1, k2 = matrices
    if any(sum(row) != 1 for row in k0):
        raise AssertionError('not stochastic')
    eq = [[k0[j][i] - int(i == j) for j in range(n)] for i in range(n)]
    eq[-1] = [Fraction(1)] * n
    stationary = solve_fraction(eq, [Fraction(0)] * (n - 1) + [Fraction(1)])
    if min(stationary) < 0 or sum(stationary) != 1:
        raise AssertionError('invalid stationary distribution')
    for j in range(n):
        assert sum(stationary[i] * k0[i][j] for i in range(n)) == stationary[j]
    g = [sum(row) for row in k1]
    mean = sum(v * r for v, r in zip(stationary, g))
    # The derivative of the Perron right eigenvector has stationary mean zero.
    h_eq = [[int(i == j) - k0[i][j] + stationary[j] for j in range(n)] for i in range(n)]
    h = solve_fraction(h_eq, [v - mean for v in g])
    assert sum(a * b for a, b in zip(stationary, h)) == 0
    variance = sum(stationary[i] * sum(k2[i]) for i in range(n)) - mean * mean
    variance += 2 * sum(stationary[i] * k1[i][j] * h[j] for i in range(n) for j in range(n))
    assert variance >= 0
    return {'mean': mean, 'variance_rate': variance, 'stationary': stationary}


def stationary_certificate(transfer: list[list[tuple[int, int]]], p: Fraction,
                           candidate: list[Fraction]) -> dict:
    """A rigorous forward-error bound, via the common empty-row reset.

    Candidate must be a probability vector. All arithmetic here is rational;
    floating residuals supplied without outward rounding are not certificates.
    """
    n = len(transfer)
    if len(candidate) != n or min(candidate) < 0 or sum(candidate) != 1:
        raise ValueError('candidate must be a normalized nonnegative vector')
    if not 0 < p < 1:
        raise ValueError('requires interior probability')
    if any(row[0][0] != 0 for row in transfer):
        raise ValueError('empty-row reset must lead to state zero')
    width = (len(transfer[0])-1).bit_length()
    weights = row_weights(width, p)
    pushed = [Fraction(0)]*n
    rewards = [Fraction(0)]*n
    for i, row in enumerate(transfer):
        for weight, (j, reward) in zip(weights, row):
            pushed[j] += candidate[i]*weight
            rewards[i] += weight*reward
    residual = sum(abs(a-b) for a,b in zip(pushed, candidate))
    reset = (1-p)**width
    estimate = sum(a*b for a,b in zip(candidate, rewards))
    error = max(rewards)*residual/reset
    return {'estimate': estimate, 'absolute_error_bound': error,
            'lower': max(Fraction(0), estimate-error), 'upper': estimate+error,
            'stationarity_l1_residual': residual, 'empty_row_reset': reset}


def symbolic_intensity(transfer: list[list[tuple[int, int]]]) -> dict:
    import sympy as sp
    from sympy.polys.matrices import DomainMatrix
    p = sp.Symbol('p')
    width = (len(transfer[0]) - 1).bit_length()
    n = len(transfer)
    k = sp.zeros(n)
    g = sp.zeros(n, 1)
    for i, row in enumerate(transfer):
        for mask, (j, reward) in enumerate(row):
            probability = p ** mask.bit_count() * (1-p) ** (width-mask.bit_count())
            k[i, j] += probability
            g[i] += probability * reward
    eq = (k.T - sp.eye(n)).applyfunc(sp.expand)
    eq[-1, :] = sp.ones(1, n)
    rhs = sp.zeros(n, 1)
    rhs[-1] = 1
    dm = DomainMatrix.from_Matrix(eq).to_field()
    solution = dm.inv().matmul(DomainMatrix.from_Matrix(rhs).convert_to(dm.domain)).to_Matrix()
    expression = sp.factor((solution.T * g)[0])
    num, den = sp.fraction(expression)
    return {'expression': str(expression),
            'numerator_descending': [int(x) for x in sp.Poly(num, p).all_coeffs()],
            'denominator_descending': [int(x) for x in sp.Poly(den, p).all_coeffs()],
            'low_p_through_8': str(sp.series(expression, p, 0, 9))}


def evaluate_polynomial(coefficients: list[int], p: Fraction) -> Fraction:
    out = Fraction(0)
    for x in coefficients:
        out = out * p + x
    return out


def graph_components(mask: int, width: int, height: int,
                     matching: bool = False, torus: bool = False) -> tuple[int, int]:
    """Independent BFS, raw (dx,dy) potentials. Return rank and essential count.

    No production DSU, frontier state, reward, or reduction is used here.
    In a free cylinder only horizontal winding can occur. Parallel lifted
    edges on width two are deliberately retained.
    """
    steps = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    if matching:
        steps += [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    potential: dict[int, tuple[int, int]] = {}
    generators: list[tuple[int, int]] = []
    count = 0
    for start in range(width * height):
        if not (mask >> start & 1) or start in potential:
            continue
        potential[start] = (0, 0)
        todo = [start]
        essential = False
        while todo:
            v = todo.pop()
            x, y = v % width, v // width
            vx, vy = potential[v]
            for dx, dy in steps:
                ny = y + dy
                if not torus and not 0 <= ny < height:
                    continue
                j = (ny % height) * width + (x + dx) % width
                if not (mask >> j & 1):
                    continue
                proposed = (vx + dx, vy + dy)
                if j not in potential:
                    potential[j] = proposed
                    todo.append(j)
                else:
                    gx, gy = proposed[0] - potential[j][0], proposed[1] - potential[j][1]
                    if gx or gy:
                        essential = True
                        generators.append((gx, gy))
        count += essential
    if not generators:
        return 0, int(count)
    x, y = generators[0]
    rank = 2 if any(x*b-y*a for a, b in generators) else 1
    return rank, int(count)


def reward_count(mask: int, width: int, height: int,
                 transfer: list[list[tuple[int, int]]]) -> int:
    state = 0
    count = 0
    for y in range(height):
        row_mask = mask >> (width*y) & ((1 << width)-1)
        state, reward = transfer[state][row_mask]
        count += reward
    _, final_reward = transfer[state][0]  # a deliberate EMPTY closing row
    return count + final_reward


def reward_histogram(width: int, height: int, transfer: list[list[tuple[int, int]]]) -> Counter:
    """Count law by total occupancy, not a false configurationwise lumping claim."""
    current = Counter({(0, 0, 0): 1})
    for _ in range(height):
        nxt = Counter()
        for (state, occupied, count), coefficient in current.items():
            for mask, (j, reward) in enumerate(transfer[state]):
                nxt[j, occupied+mask.bit_count(), count+reward] += coefficient
        current = nxt
    histogram = Counter()
    for (state, occupied, count), coefficient in current.items():
        reward = transfer[state][0][1]
        histogram[occupied, count+reward] += coefficient
    return histogram


def central_trinomial(n: int) -> int:
    return sum(comb(n, k) * comb(n-k, k) for k in range(n//2+1))


def renewal_controls() -> dict:
    """Exactly solvable directed renewal LOOP; not a site-cluster theorem."""
    controls = []
    for w in (8, 16, 32, 64, 128):
        c = [central_trinomial(k*w) for k in (1, 2, 3)]
        ratio = Fraction(c[0]*c[2], c[1]**2)  # 3 powers and exponential fugacity cancel
        beta_eff = log(float(ratio)) / log(4/3)
        amplitude_ratio = float(Fraction(c[0], 3**w)) * sqrt(w) / sqrt(3/(4*pi))
        controls.append({'width': w, 'exact_tripling_ratio': str(ratio),
                         'effective_beta': beta_eff,
                         'normalized_amplitude': amplitude_ratio})
    return {'model': 'X=1, Y=-1,0,1 with equal probabilities; distinct from site percolation',
            'predicted_beta': '1/2', 'transverse_diffusion': '2/3',
            'predicted_amplitude_squared_times_pi': '3/4', 'controls': controls}


def centering_counterexample() -> list[dict]:
    """C1/semiconcave mass is insufficient for an o(1/w) linearized centre.

    Exact model on h>0: kappa(a+h)=d-h-h/log(e/h).  The intensity is
    w^{-1} exp(-w*kappa) and log m=dw. At intensity level one, solve
    h+h/log(e/h)=log(w)/w. Then w*(h-log(w)/w) -> -1, not zero.
    Numbers below are Decimal diagnostics; the limit has a direct proof.
    """
    from decimal import Decimal, localcontext
    result=[]
    with localcontext() as ctx:
        ctx.prec=80
        for exponent in (4,8,16,32):
            w=Decimal(10)**exponent
            linear=w.ln()/w
            lo,hi=Decimal(0),linear
            for _ in range(400):
                h=(lo+hi)/2
                if h+h/(1-h.ln())<linear:lo=h
                else:hi=h
            exact=(lo+hi)/2
            result.append({'width':'10^'+str(exponent),
                           'scaled_linearization_error':str(w*(exact-linear))})
    return result


def run_report() -> dict:
    import sympy as sp
    start = perf_counter()
    tables = {}
    generators = []
    physical_checks = 0
    symbolic_checks = 0
    for width in (2, 3, 4):
        for matching in (False, True):
            states, transfer = build_transfer(width, matching)
            reduced, blocks = reward_lump(transfer)
            tables[width, matching] = transfer
            expression = symbolic_intensity(reduced)
            point_controls = []
            for p in (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4)):
                exact = stationary_reward(reduced, p)
                rational = evaluate_polynomial(expression['numerator_descending'], p) / evaluate_polynomial(expression['denominator_descending'], p)
                assert exact['mean'] == rational
                point_controls.append({'p': str(p), 'intensity': str(exact['mean']),
                                       'variance_per_row': str(exact['variance_rate']),
                                       'fano_rate': str(exact['variance_rate']/exact['mean'])})
                truncated = [Fraction(int(x*(1 << 20)), 1 << 20) for x in exact['stationary']]
                truncated[0] += 1-sum(truncated)
                certificate = stationary_certificate(reduced, p, truncated)
                assert certificate['lower'] <= exact['mean'] <= certificate['upper']
                point_controls[-1]['rounded_stationary_certificate'] = {k: str(v) for k,v in certificate.items()}
                symbolic_checks += 1
            height = 4
            histogram = Counter()
            for mask in range(1 << (width*height)):
                _, observed = graph_components(mask, width, height, matching, torus=False)
                assert reward_count(mask, width, height, transfer) == observed
                histogram[mask.bit_count(), observed] += 1
                physical_checks += 1
            assert reward_histogram(width, height, reduced) == histogram
            generators.append({'width': width, 'graph': 'matching' if matching else 'NN',
                               'frontier_states': len(states), 'reward_lumps': len(reduced),
                               'max_absolute_gain': max(abs(v) for state in states for v in state.gains),
                               'transition_entries': len(states)*(1 << width),
                               'symbolic_intensity': expression,
                               'point_controls': point_controls,
                               'states': [list(state) for state in states],
                               'transition_table': transfer, 'reward_lump_map': blocks,
                               'reduced_table': reduced})
    p = sp.Symbol('p')
    for w in (2, 3, 4):
        g4 = next(g for g in generators if g['width']==w and g['graph']=='NN')
        g8 = next(g for g in generators if g['width']==w and g['graph']=='matching')
        e4 = sp.sympify(g4['symbolic_intensity']['expression'])
        e8 = sp.sympify(g8['symbolic_intensity']['expression'])
        assert sp.cancel(e4-e8.subs(p,1-p)) == 0
        for point in g4['point_controls']:
            paired = next(v for v in g8['point_controls'] if Fraction(v['p'])==1-Fraction(point['p']))
            assert paired['intensity'] == point['intensity']
            assert paired['variance_per_row'] == point['variance_per_row']
    duality_checks = 0
    joint_support = set()
    for width, height in ((2,2), (3,3), (4,4)):
        full = (1 << (width*height)) - 1
        for mask in range(full+1):
            rank4, n4 = graph_components(mask, width, height, False, torus=True)
            rank8, n8 = graph_components(full^mask, width, height, True, torus=True)
            assert rank4+rank8 == 2
            assert n4-n8 == rank4-1
            assert (n4,n8) in ((1,0),(0,1)) or n4==n8>=1
            joint_support.add((n4,n8))
            duality_checks += 1
    return {'schema': 'matching-one.cylinder-winding-intensity.v1',
            'scope': 'two graphs; exact fixed-width component rewards, not a large-width prefactor proof',
            'models': generators,
            'checks': {'open_cylinder_graph_configurations': physical_checks,
                       'torus_complement_pairs': duality_checks,
                       'symbolic_point_controls': symbolic_checks,
                       'all_p_duality_identities': 3,
                       'joint_torus_support_observed': sorted(joint_support)},
            'renewal': renewal_controls(),
            'centering_counterexample': centering_counterexample(),
            'elapsed_seconds': perf_counter()-start}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--probe-width', type=int)
    parser.add_argument('--state-cap', type=int, default=5000)
    args = parser.parse_args()
    if args.probe_width is not None:
        start = perf_counter()
        states, transfer = build_transfer(args.probe_width, state_cap=args.state_cap)
        print(json.dumps({'width': args.probe_width, 'states': len(states),
                          'transition_entries': sum(map(len,transfer)),
                          'elapsed_seconds': perf_counter()-start}, indent=2))
        return
    if args.output is None:
        parser.error('--output is required unless --probe-width is used')
    if args.output.exists():
        raise FileExistsError(f'refusing to overwrite {args.output}')
    result = run_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result['checks'], indent=2))
    print(f'elapsed_seconds={result["elapsed_seconds"]:.3f}')


if __name__ == '__main__':
    main()
