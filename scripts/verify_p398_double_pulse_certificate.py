#!/usr/bin/env python3
"""Check the stored Hankel rank witnesses by integer, fraction-free elimination.

Rebuild each selected word coefficient before computing its determinant. This
checks the claim-bearing ranks, not document wording or a repository registry.
Requires only the standard library. Output paths use exclusive creation.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from p398_double_pulse_exact import (build_model, matvec, multiply, parity_basis, transpose)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / 'results/research-control-20260912/p398-double-pulse-exact.json'


def bareiss(matrix):
    m = [row[:] for row in matrix]
    sign, previous, n = 1, 1, len(m)
    if not n:
        return 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if m[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            m[k], m[pivot] = m[pivot], m[k]
            sign = -sign
        value = m[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = value * m[i][j] - m[i][k] * m[k][j]
                if numerator % previous:
                    raise ValueError('Nonexact Bareiss division')
                m[i][j] = numerator // previous
            m[i][k] = 0
        previous = value
    return sign * m[-1][-1]


def word_minor(operators, f, src, cert):
    rows, columns = [], []
    for item in cert['row_words']:
        vector = src[item['seed']][:]
        for name in item['word']:
            vector = matvec(transpose(operators[name]), vector)
        rows.append(vector)
    for item in cert['column_words']:
        vector = [row[item['seed']] for row in f]
        for name in reversed(item['word']):
            vector = matvec(operators[name], vector)
        columns.append(vector)
    return [[sum(x*y for x, y in zip(row, col)) for col in columns] for row in rows]


def verify(payload):
    checks = []
    for row in payload['widths']:
        width = row['width']
        model = build_model(width)
        for key, other in [('G','G'), ('H','H'), ('sources_scaled','sources_scaled'), ('readouts','F')]:
            if row[key] != model[other]:
                raise ValueError(f'Model mismatch at width {width}: {key}')
        g, h, f, src = (model[k] for k in ('G','H','F','sources_scaled'))
        u, reps = parity_basis(model, 1)
        gp = [multiply(g, u)[i] for i in reps]
        fp, sp = [f[i] for i in reps], multiply(src, u)
        odd = row['double_pulse']
        u, reps = parity_basis(model, -1)
        gm = [multiply(g, u)[i] for i in reps]
        bm, cm = [multiply(h, f)[i] for i in reps], multiply(multiply(src,h),u)
        for got, expected in [(odd['G_odd'],gm),(odd['B_odd'],bm),(odd['C_odd_scaled'],cm)]:
            if got != expected:
                raise ValueError('Odd sector mismatch')
        contrast = [[src[i][j]-src[3][j] for j in range(len(g))] for i in range(3)]
        cases = [
            ('baseline', {'G':gp}, fp, sp, row['baseline']),
            ('double_pulse', {'G':gm}, bm, cm, odd['rank_certificate']),
            ('controlled', {'G':g,'H':h}, f, src, row['arbitrary_controlled_word']),
            ('contrasts', {'G':g,'H':h}, f, contrast, row['controlled_source_contrasts']),
        ]
        for name, operators, inputs, outputs, cert in cases:
            matrix = word_minor(operators, inputs, outputs, cert)
            if matrix != cert['integer_hankel_minor']:
                raise ValueError(f'Word coefficient mismatch: {width}, {name}')
            det = bareiss(matrix)
            if not det or det % cert['prime'] != cert['hankel_minor_determinant_mod_prime']:
                raise ValueError(f'Rank witness failed: {width}, {name}')
            checks.append(dict(width=width, object=name, dimension=cert['minimal_order'],
                               exact_integer_determinant=str(det), modular_residue_matches=True))
        scalar = odd['scalar_witness']
        d = len(scalar['hankel_minor'])
        matrix = [[scalar['moments'][i+j] for j in range(d)] for i in range(d)]
        if matrix != scalar['hankel_minor']:
            raise ValueError('Scalar Hankel mismatch')
        det = bareiss(matrix)
        if not det or det % 1_000_000_007 != scalar['hankel_determinant_mod_prime']:
            raise ValueError('Scalar rank witness failed')
        checks.append(dict(width=width, object='scalar_pulse', dimension=d,
                           exact_integer_determinant=str(det), modular_residue_matches=True))
    return dict(schema='matching-one.p398-double-pulse-bareiss-check.v1',
                method='Independent fraction-free Bareiss elimination of stored integer Hankel matrices, with exact division at every step.',
                checks=checks)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DEFAULT)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    text = json.dumps(verify(json.loads(args.input.read_text(encoding='utf-8'))), indent=2) + '\n'
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open('x', encoding='utf-8') as handle:
            handle.write(text)
    else:
        print(text, end='')
