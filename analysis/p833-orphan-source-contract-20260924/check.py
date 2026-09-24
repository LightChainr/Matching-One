#!/usr/bin/env python3
"""Exact controls for the local orphan-component identity and source contracts.

The physical input and input/p833_check.py are unmodified artifacts of PR834.
That module is used explicitly for the already established boundary/oracle
algorithms, not silently installed on main. This file adds no width engine.
All mathematical checks use bounded integers, GF(2), a finite prime field, or
Fraction. Floating time/RSS are telemetry only. Author proof, not external review.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
import importlib.util
from itertools import product
import json
from pathlib import Path
import resource
import sys
import time
import numpy as np


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location('p833_reference', path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f'Cannot load reference module: {path}')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def label(keys):
    index, result = {}, []
    for key in keys:
        if key not in index:
            index[key] = len(index)
        result.append(index[key])
    return result


def canonical_relation(rel):
    rel = sorted((int(i), int(v)) for i, v in rel if v)
    if rel[0][1] < 0:
        rel = [(i, -v) for i, v in rel]
    return tuple(rel)


def gf2_bits(a):
    return int.from_bytes(np.packbits(a, bitorder='little').tobytes(), 'little')


def gf2_rank(vectors, low_pivot=False):
    basis = {}
    for x in vectors:
        while x:
            pivot = ((x & -x).bit_length()-1) if low_pivot else x.bit_length()-1
            if pivot not in basis:
                basis[pivot] = x
                break
            x ^= basis[pivot]
    return len(basis)


def rank_mod(matrix, prime):
    a = np.array(matrix, dtype=np.int64).copy() % prime
    row = 0
    for col in range(a.shape[1]):
        candidates = np.flatnonzero(a[row:, col])
        if not len(candidates):
            continue
        k = row + int(candidates[0])
        a[[row, k]] = a[[k, row]]
        a[row] = a[row] * pow(int(a[row, col]), -1, prime) % prime
        for j in range(row+1, len(a)):
            a[j] = (a[j] - int(a[j, col])*a[row]) % prime
        row += 1
        if row == len(a):
            break
    return row


def discover_raw_relations(c):
    """Sparse modular discovery, followed elsewhere by exact integer checks.

    A row signature has one nonzero entry for each deterministic input mask.
    Work separately in the three current-output blocks. Every returned row is
    an exact signed square, not a modular dependence relabeled as rational.
    """
    q, out = c['quotient_transitions'], c['quotient_rank_output']
    n, prime = len(q), 1000003
    relations = []
    for rank in range(3):
        basis = {}
        for i in (j for j in range(n) if out[j] == rank):
            row = {b*n+t: 1 for b, t in enumerate(q[i])}
            coeff = {i: 1}
            while row:
                pivot, scale = min(row), row[min(row)]
                if pivot not in basis:
                    inv = pow(scale, -1, prime)
                    row = {j: a*inv % prime for j, a in row.items()}
                    coeff = {j: a*inv % prime for j, a in coeff.items()}
                    basis[pivot] = (row, coeff)
                    break
                br, bc = basis[pivot]
                for dest, source in ((row, br), (coeff, bc)):
                    for j, a in source.items():
                        value = (dest.get(j, 0)-scale*a) % prime
                        if value:
                            dest[j] = value
                        else:
                            dest.pop(j, None)
            if not row:
                rel = [[j, a if a <= prime//2 else a-prime] for j, a in coeff.items()]
                check = Counter()
                for j, a in rel:
                    for b, t in enumerate(q[j]):
                        check[b*n+t] += a
                assert all(a == 0 for a in check.values())
                relations.append(rel)
    return relations


def discover_controlled_words(c, target):
    """Generate actual Boolean controlled experiments, not reduced fake words."""
    q = np.array(c['quotient_transitions'], dtype=np.int64)
    out = (np.array(c['quotient_rank_output']) == 1).astype(np.uint8)
    n, alphabet = q.shape
    basis, cols, dag, depths, pivots = {}, [], [], [], []
    def add(col, description, depth):
        original = gf2_bits(col)
        value = original
        while value:
            pivot = value.bit_length()-1
            if pivot not in basis:
                basis[pivot] = value
                cols.append(original); dag.append(description)
                depths.append(depth); pivots.append(pivot)
                return
            value ^= basis[pivot]
    for r in (0, 1):
        add((out == r).astype(np.uint8), [r, -1, -1], 0)
    cursor = 0
    while cursor < len(cols) and len(cols) < target:
        raw = np.frombuffer(cols[cursor].to_bytes((n+7)//8, 'little'), dtype=np.uint8)
        f = np.unpackbits(raw, bitorder='little')[:n]
        successors = f[q]
        for mask in range(alphabet):
            for r in (0, 1):
                add(successors[:, mask]*(out == r), [r, mask, cursor], depths[cursor]+1)
        cursor += 1
    assert len(cols) == target
    return {'rank': target, 'field': 2, 'witness_dag': dag,
            'depths': depths, 'pivots': pivots}


def actions(c, hs, src):
    w = c['w']
    q, ini = c['quotient_transitions'], c['quotient_initial']
    mask_maps, state_maps = [], []
    for sign, shift in product((1, -1), range(w)):
        mm = [sum(((m >> j) & 1) << ((sign*j+shift) % w)
                  for j in range(w)) for m in range(1 << w)]
        aa = [src.state_for([mm[m] for m in h], q, ini) for h in hs]
        assert len(set(aa)) == len(q)
        mask_maps.append(mm)
        state_maps.append(aa)
    return np.array(mask_maps), state_maps


def check_local_relations(c, d, seed, src, hs, state_maps):
    q = c['quotient_transitions']
    ranks = c['quotient_rank_output']
    reps = [c['states'][i] for i in c['quotient_representatives']]
    relations = seed['raw_relations']
    w, n = c['w'], len(q)
    assert len(relations) == 110
    assert len({max(i for i, a in rel) for rel in relations}) == 110
    counts, local_records, geometric_updates = Counter(), [], 0
    for rel in relations:
        assert sorted(a for i, a in rel) == [-1, -1, 1, 1]
        assert len({ranks[i] for i, a in rel}) == 1
        counts[ranks[rel[0][0]]] += 1
        states = [reps[i] for i, a in rel]
        active = [{t for t, (root, gain) in enumerate(s[1]) if root >= 0}
                  for s in states]
        common, union = set.intersection(*active), set.union(*active)
        optional = sorted(union-common)
        assert len(optional) == 2 and min(optional) >= w
        i, j = [t-w for t in optional]
        if (i+1) % w != j:
            i, j = j, i
        assert (i+1) % w == j
        arc = [(i-1) % w, i, j, (j+1) % w]
        endpoints = [w+arc[0], w+arc[-1]]
        assert set(endpoints) <= common
        assert len({s[0] for s in states}) == 1
        occupancy_patterns = {tuple(int(t in a) for t in optional) for a in active}
        assert occupancy_patterns == {(0, 0), (0, 1), (1, 0), (1, 1)}
        for idx, s in enumerate(states):
            entries = s[1]
            root = entries[endpoints[0]][0]
            block = {t for t, (r, gain) in enumerate(entries) if r == root}
            assert block == active[idx] & {w+t for t in arc}
            assert min(block) >= w  # the component does not meet the pinned seam
            if not s[0]:
                shift, expected = 0, {arc[0]: 0}
                for a, b in zip(arc, arc[1:]):
                    shift += int(a == w-1 and b == 0)
                    expected[b] = shift
                for t in block:
                    assert entries[t][1]-entries[endpoints[0]][1] == expected[t-w]
            # Removing the two optional terminals leaves the same old skeleton.
            ds, keep = src.DSU(2*w, s[0]), []
            for r in {r for r, g in entries if r >= 0}:
                vv = [t for t, (rr, g) in enumerate(entries)
                      if rr == r and t not in optional]
                if vv:
                    first = vv[0]
                    for t in vv:
                        ds.edge(first, t, entries[t][1]-entries[first][1])
                        keep.append((t, t))
            reduced = src.canonical(ds, keep, w)
            if idx == 0:
                template = reduced
            assert reduced == template
        # Whole next boundary states, not just equal terminal-rank numbers.
        for mask in range(1 << w):
            native, graph = Counter(), Counter()
            for state, coefficient in rel:
                native[q[state][mask]] += coefficient
                graph[src.oracle_advance(reps[state], mask)] += coefficient
                geometric_updates += 1
            assert all(v == 0 for v in native.values())
            assert all(v == 0 for v in graph.values())
        local_records.append({'states': [i for i, a in rel],
                              'optional_columns': [i, j], 'arc': arc})
    # Independent integer relations with distinct final columns; their span
    # annihilates every controlled word since each row is in one output block.
    lookup = {canonical_relation(rel): i for i, rel in enumerate(relations)}
    seen, orbit_records = set(), []
    for r, rel in enumerate(relations):
        if r in seen:
            continue
        orbit = {lookup[canonical_relation((g[i], a) for i, a in rel)]
                 for g in state_maps}
        seen.update(orbit)
        projected = Counter()
        for i, a in rel:
            projected[d['dlabel'][i]] += a
        projected = canonical_relation(projected.items())
        match = [k for k, old in enumerate(seed['parent_kernel_rows'])
                 if canonical_relation(old) == projected]
        assert len(match) == 1
        orbit_records.append({'members': sorted(orbit), 'size': len(orbit),
                              'parent_relation': match[0], 'projected': projected})
    assert len(seen) == 110 and len(orbit_records) == 13
    assert Counter(x['size'] for x in orbit_records) == {5: 4, 10: 9}
    assert sorted(x['parent_relation'] for x in orbit_records) == list(range(13))
    return {'independent_relations': 110, 'output_support_counts': dict(counts),
            'whole_boundary_update_checks': geometric_updates,
            'orbit_sizes': sorted(x['size'] for x in orbit_records),
            'local_records': local_records, 'relation_orbits': orbit_records}


def binary_refinement(c):
    q = c['quotient_transitions']
    lab = [int(r == 1) for r in c['quotient_rank_output']]
    counts = [len(set(lab))]
    while True:
        nxt = label((lab[i], tuple(lab[j] for j in row)) for i, row in enumerate(q))
        counts.append(len(set(nxt)))
        if len(set(nxt)) == len(set(lab)):
            break
        lab = nxt
    assert counts == [2, 808, 3318, 3438, 3438]
    return counts


def controlled_observability(c, record):
    q = np.array(c['quotient_transitions'], dtype=np.int64)
    out = (np.array(c['quotient_rank_output']) == 1).astype(np.uint8)
    dag = record['witness_dag']
    columns, depths = [], []
    for index, (r, mask, tail) in enumerate(dag):
        assert r in (0, 1)
        if tail < 0:
            assert mask == -1
            value, depth = (out == r).astype(np.uint8), 0
        else:
            assert 0 <= tail < index and 0 <= mask < 32
            value = (out == r)*columns[tail][q[:, mask]]
            depth = depths[tail]+1
        columns.append(value)
        depths.append(depth)
    assert len(columns) == 3328 and max(depths) == 5
    assert depths == record['depths']
    rank = gf2_rank(gf2_bits(v) for v in columns)
    assert rank == 3328
    # Opposite orientation and opposite pivot order on the claimed square minor.
    matrix = np.array(columns, dtype=np.uint8).T[record['pivots'], :]
    assert matrix.shape == (3328, 3328)
    assert gf2_rank((gf2_bits(row) for row in matrix), low_pivot=True) == 3328
    return {'dimension': rank, 'field': 2, 'actual_controlled_word_minor': 3328,
            'minor_determinant_mod_2': 1, 'maximum_future_rows': max(depths),
            'both_row_and_column_elimination_pass': True}


def uniform_minor(d, record):
    """Replay the actual 385-word integer minor; no floating-point ranks."""
    prime, xx = record['prime'], record['odds']
    n = len(d['output'])
    out = np.array([int(r == 1) for r in d['output']])
    A = np.zeros((n, n), dtype=np.int64)
    for i, row in enumerate(d['weighted']):
        for k, j, count in row:
            A[i, j] += count*pow(xx, k, prime)
    A %= prime
    cache = {}
    def column(word):
        word = tuple(map(int, word))
        if word not in cache:
            cache[word] = ((out == word[0]).astype(np.int64) if len(word) == 1
                           else (out == word[0])*(A@column(word[1:]) % prime))
        return cache[word]
    matrix = np.array([column(word)[record['pivots']] for word in record['words']], dtype=np.int64).T
    assert matrix.shape == (385, 385)
    det = 1
    for k in range(385):
        choices = np.flatnonzero(matrix[k:, k]); assert len(choices)
        row = k+int(choices[0])
        if row != k:
            matrix[[row,k]] = matrix[[k,row]]; det = -det
        pivot = int(matrix[k,k]); det = det*pivot % prime
        factors = matrix[k+1:,k]*pow(pivot,-1,prime) % prime
        matrix[k+1:,k:] = (matrix[k+1:,k:] - factors[:,None]*matrix[k,k:]) % prime
    assert det == record['det'] == 821642
    return {'dimension_lower_bound': 385, 'prime':prime, 'odds_mod':xx,
            'minor_determinant_mod':int(det), 'maximum_word_length':max(map(len,record['words'])),
            'parent_witness_actual_columns_recomputed':True}


def orbit_positive_model(c, d, seed, src, mask_maps):
    original = src.positive_model(d, seed['parent_kernel_rows'])
    n, m = len(d['output']), original['states']
    U2, V = np.zeros((n, m), dtype=np.int64), np.zeros((m, n), dtype=np.int64)
    for i, row in enumerate(original['U2']):
        for j, a in row:
            U2[i, j] = a
    for i, row in enumerate(original['V']):
        for j, a in row:
            V[i, j] = a
    q = np.array(c['quotient_transitions'])
    dl, reps = np.array(d['dlabel']), np.array(d['dreps'])
    tr = dl[q[reps]]
    orbits = sorted({tuple(sorted(set(map(int, mask_maps[:, b])))) for b in range(32)})
    assert len(orbits) == 8
    N = np.zeros((13, n), dtype=np.int64)
    for row, terms in enumerate(seed['parent_kernel_rows']):
        for j, a in terms:
            N[row, j] = a
    orbit_records, compressed = [], []
    for orbit in orbits:
        C = np.zeros((n, n), dtype=np.int64)
        for i in range(n):
            for b in orbit:
                C[i, tr[i, b]] += 1
        assert np.all(C.sum(axis=1) == len(orbit))
        assert not np.any(N@C)
        VC = V@C
        assert np.all(VC >= 0) and np.array_equal(U2@VC, 2*C)
        T2 = VC@U2
        assert np.all(T2 >= 0)
        assert np.array_equal(U2@T2, 2*C@U2)
        assert np.all(T2.sum(axis=1) == 2*len(orbit))
        row = orbit[0]
        adjacent = sum(((row >> i) & 1)*((row >> ((i+1) % 5)) & 1) for i in range(5))
        orbit_records.append({'masks': orbit, 'size': len(orbit),
                              'occupied_sites': row.bit_count(),
                              'occupied_neighbor_pairs': adjacent,
                              'minimum_VC_entry': int(VC.min())})
        compressed.append({'masks': orbit, 'denominator': 2*len(orbit),
                           'rows': [[[int(j), int(T2[i, j])]
                                     for j in np.flatnonzero(T2[i])] for i in range(m)]})
    # The 24 necessary constraints use only histories already mapped to the
    # same D5 class, hence to the same fixed F=orbit_map*U encoding.
    rows = []
    for s, t, latent in seed['source_law_constraints']['witnesses']:
        assert dl[s] == dl[t]
        rows.append([int(U2[dl[q[s, b]], latent]-U2[dl[q[t, b]], latent])
                     for b in range(32)])
    rows = np.array(rows, dtype=np.int64)
    assert rank_mod(rows, seed['source_law_constraints']['prime']) == 24
    for orbit in orbits:
        assert np.all(rows[:, orbit].sum(axis=1) == 0)
    # Parent's real word minor supplies a lower bound for the universal family
    # since uniform Bernoulli(1/2) is one member. Recomputed, not just a flag.
    lower = uniform_minor(d, seed['parent_generic_minor'])
    model = {'schema': 'matching-one.width5.orbit-source-positive.v1',
             'states': m, 'outputs': original['outputs'], 'U2': original['U2'],
             'source_orbit_operators': compressed}
    return {'state_count': m, 'orbit_count': 8, 'orbit_records': orbit_records,
            'all_orbit_coefficients_nonnegative': True,
            'all_orbit_intertwinings_exact': True,
            'fixed_encoding_source_constraint_rank': 24,
            'allowed_unnormalized_source_dimension': 8,
            'allowed_probability_simplex_dimension': 7,
            'uniform_source_dimension_lower_bound': lower}, model


def physical_examples(c, src):
    q, ini = c['quotient_transitions'], c['quotient_initial']
    out = c['quotient_rank_output']
    square = [[0, 15], [0, 7, 13], [0, 14, 11], [0, 15, 9]]
    checks = 0
    # Equal-length representatives used for the physical square witness.
    square_equal = [[0]+h if len(h) == 2 else h for h in square]
    ids = [src.state_for(h, q, ini) for h in square]
    assert [src.state_for(h, q, ini) for h in square_equal] == ids
    assert all(out[i] == 0 for i in ids)
    for word in product(range(32), repeat=2):
        transcripts = []
        for h, state in zip(square_equal, ids):
            rword = []
            for step in range(1, 3):
                got = src.graph_rank(5, h+list(word[:step]))
                exp = out[src.state_for(h+list(word[:step]), q, ini)]
                assert got == exp
                rword.append(got)
                checks += 1
            transcripts.append(tuple(rword))
        assert sorted([transcripts[0], transcripts[3]]) == sorted([transcripts[1], transcripts[2]])
    # One independently addressed column breaks the pre-existing orientation quotient.
    h1, h2 = [1, 1], [2, 2]
    coeff = [0, 0, 0]
    for b in range(32):
        sign = 1 if b & 1 else -1
        r1 = src.graph_rank(5, h1+[b]); r2 = src.graph_rank(5, h2+[b]); checks += 2
        coeff[r1] += sign
        coeff[r2] -= sign
    assert coeff == [-16, 16, 0]
    return {'physical_graph_rank_checks': checks,
            'raw_square_equal_length_histories': square_equal,
            'square_signed_coefficients': [1, -1, -1, 1],
            'single_column_histories': [h1, h2],
            'single_column_probability_response': ['-delta', 'delta', '0']}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input', type=Path, help='PR834 generated width5-rank-closure.json')
    ap.add_argument('--reference', type=Path, default=Path(__file__).parent/'input/p833_check.py')
    ap.add_argument('--parent-seed', type=Path, default=Path(__file__).parent/'input/p833_certificate.json')
    ap.add_argument('--certificate-out', type=Path, help='Export generated exact witness DAG and local relations')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--model-out', type=Path)
    ap.add_argument('--full-input-check', action='store_true',
                    help='Rerun PR834 independent whole-boundary and physical checks')
    args = ap.parse_args()
    started = time.perf_counter()
    for p in (args.input, args.reference, args.parent_seed):
        if not p.is_file():
            ap.error(f'Missing file: {p}')
    src = load_module(args.reference)
    c = json.loads(args.input.read_text()); parent = json.loads(args.parent_seed.read_text())
    c['states'] = [(int(h), tuple(tuple(e) for e in entries)) for h, entries in c['states']]
    assert c['w'] == 5 and len(c['states']) == 10766
    assert len(c['quotient_transitions']) == 3438
    # The source cert is read, not rebuilt and mislabeled as new production.
    ss = c['states']; ts = c['transitions']; qm = c['quotient_map']
    for i, row in enumerate(ts):
        assert len(row) == 32 and all(0 <= j < len(ss) for j in row)
        assert [qm[j] for j in row] == c['quotient_transitions'][qm[i]]
        assert c['rank_output'][i] == c['quotient_rank_output'][qm[i]]
    d = src.dihedral(c, full_check=True)
    assert len(d['output']) == 398
    geometry = src.geometry_checks(c, d) if args.full_input_check else {'not_rerun': True}
    print('input and D5 checks passed', flush=True)
    seed = {'parent_commit': '8f40247491af737c03aa3b1156db5918a7024093',
            'parent_kernel_rows': parent['kernel_rows'],
            'parent_generic_minor': parent['generic_minor'],
            'raw_relations': discover_raw_relations(c),
            'controlled_binary_witnesses': discover_controlled_words(c, 3328),
            'source_law_constraints': {'prime': 1000003, 'witnesses':
                [[2,1,0],[2,1,35],[2,1,36],[2,1,43],[2,1,46],[2,1,47],
                 [4,1,0],[4,1,35],[4,1,36],[4,1,43],[4,1,46],[4,1,47],
                 [6,3,35],[8,1,35],[8,1,36],[9,5,35],[9,5,36],
                 [10,5,35],[10,5,36],[13,11,35],[16,1,35],[18,5,35],[18,5,36],[20,5,35]]},
            'epsilon_example': [1,1000000000]}
    hs = src.histories(c['quotient_transitions'], c['quotient_initial'])
    mm, aa = actions(c, hs, src)
    local = check_local_relations(c, d, seed, src, hs, aa)
    print('110 local relations and 13 orbits passed', flush=True)
    refinement = binary_refinement(c)
    obs = controlled_observability(c, seed['controlled_binary_witnesses'])
    print('controlled observable dimension 3328 passed', flush=True)
    positive, model = orbit_positive_model(c, d, seed, src, mm)
    print('all 8 orbit operators and maximal source law passed', flush=True)
    examples = physical_examples(c, src)
    n, width, depth = 3438, 5, 3
    denominator = width*depth*n*(n-1)
    eps = Fraction(*seed['epsilon_example'])
    overlap_budget = denominator*eps
    assert denominator == 177246090 and overlap_budget < 1
    # Algebraic inverse check for the two-bit noise matrix. Its fifth tensor
    # power is invertible for epsilon != 1/2, so linear dimensions agree.
    det_noise = 1-2*eps
    assert det_noise != 0
    a,b = 1-eps,eps
    assert (a*a-b*b)/det_noise == 1
    result = {'schema': 'matching-one.orphan-source-contract.result.v1',
              'date': '2026-09-24', 'parent_pr': 834,
              'parent_commit': seed['parent_commit'],
              'scope': 'width5 NN; controlled row extension; current output included; all pure histories and mixtures; not an autonomous near-critical gap',
              'input_checks': geometry, 'dihedral_transition_relations': 3438*32*10,
              'local_mechanism': {k: v for k,v in local.items() if k not in ['local_records','relation_orbits']}, 'controlled_binary_refinement': refinement,
              'controlled_linear_dimension': obs,
              'controlled_positive_minimum': 3438,
              'controlled_positive_linear_gap': 110,
              'positive_invariant_sources': positive,
              'physical_examples': examples,
              'strictly_interior_control_bound': {
                  'epsilon_interval': ['0', f'1/{denominator}'],
                  'example_epsilon': str(eps), 'commands': 32,
                  'sites_per_distinguishing_experiment_at_most': width*depth,
                  'pairwise_initial_overlap_bound': '30 epsilon',
                  'total_overlap_budget_at_example': str(overlap_budget),
                  'positive_minimum_on_interval': 3438,
                  'linear_dimension_on_interval': 3328},
              'not_claimed': ['gap for autonomous homogeneous p near pc',
                              'minimum for arbitrary intermediate command noise',
                              'all-width completeness of local generators',
                              'independent external review or publication novelty',
                              'full repository CI or change to original-U']}
    if args.certificate_out:
        args.certificate_out.parent.mkdir(parents=True, exist_ok=True)
        seed['local_geometry'] = local
        args.certificate_out.write_text(json.dumps(seed, ensure_ascii=False, separators=(',',':'))+'\n')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    if args.model_out:
        args.model_out.parent.mkdir(parents=True, exist_ok=True)
        args.model_out.write_text(json.dumps(model, ensure_ascii=False, separators=(',',':'))+'\n')
    print(json.dumps({'seconds': time.perf_counter()-started,
                      'peak_RSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      'python': sys.version.split()[0], 'numpy': np.__version__}), flush=True)


if __name__ == '__main__':
    main()
