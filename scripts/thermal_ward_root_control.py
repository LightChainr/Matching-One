#!/usr/bin/env python3
"""Conditional thermal Ward identities and existing-data root diagnostics.

No new percolation production is performed. Exact algebra uses Fraction;
mpmath evaluates modular series and roots and is not interval arithmetic.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from functools import lru_cache
from math import comb
import json
from pathlib import Path

H = F(5, 8)
C = F(0)
Vector = dict[tuple[int, ...], F]


def add(*vectors: Vector) -> Vector:
    out: Vector = {}
    for vec in vectors:
        for word, coeff in vec.items():
            out[word] = out.get(word, F(0)) + coeff
    return {w: c for w, c in out.items() if c}


def scale(a: F, vec: Vector) -> Vector:
    return {w: a*c for w, c in vec.items() if a*c}


@lru_cache(None)
def normal(word: tuple[int, ...]) -> Vector:
    """PBW basis uses decreasing positive indices for negative modes."""
    for j in range(len(word)-1):
        a, b = word[j:j+2]
        if a < b:
            swapped = word[:j] + (b, a) + word[j+2:]
            merged = word[:j] + (a+b,) + word[j+2:]
            return add(normal(swapped), scale(F(b-a), normal(merged)))
    return {word: F(1)}


@lru_cache(None)
def act_word(m: int, word: tuple[int, ...], h: F = H, c: F = C) -> Vector:
    if m < 0:
        return normal((-m,) + word)
    if m == 0:
        return {word: h + sum(word)}
    if not word:
        return {}
    n, rest = word[0], word[1:]
    first: Vector = {}
    for w, a in act_word(m, rest, h, c).items():
        first = add(first, scale(a, normal((n,) + w)))
    second = scale(F(m+n), act_word(m-n, rest, h, c))
    central = {rest: c*F(m**3-m, 12)} if m == n and c else {}
    return add(first, second, central)


def act(m: int, vec: Vector, h: F = H, c: F = C) -> Vector:
    out: Vector = {}
    for word, coeff in vec.items():
        out = add(out, scale(coeff, act_word(m, word, h, c)))
    return out


def descendants() -> tuple[Vector, Vector, Vector]:
    chi = {(2,): F(1), (1, 1): -F(2, 3)}
    quotient_q = {(4,): F(1), (3, 1): F(20, 11), (1, 1, 1, 1): -F(160, 561)}
    qhat = scale(-F(11, 29), add(quotient_q, scale(-F(80, 33), act(-2, chi))))
    return chi, quotient_q, qhat


def partitions(n: int, top: int | None = None):
    if n == 0:
        yield ()
        return
    for k in range(min(n, n if top is None else top), 0, -1):
        for tail in partitions(n-k, k):
            yield (k,) + tail


def rational_rank(vectors: list[Vector], basis: list[tuple[int, ...]]) -> int:
    rows = [[v.get(w, F(0)) for w in basis] for v in vectors]
    rank = 0
    for col in range(len(basis)):
        pivot = next((j for j in range(rank, len(rows)) if rows[j][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        p = rows[rank][col]
        rows[rank] = [x/p for x in rows[rank]]
        for j in range(len(rows)):
            if j != rank and rows[j][col]:
                t = rows[j][col]
                rows[j] = [x-t*y for x, y in zip(rows[j], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def null_derivative_subspace() -> list[Vector]:
    chi, _, _ = descendants()
    return [act(-1, {w: F(1)}) for w in partitions(3)] + [
        act(-2, chi), act(-1, act(-1, chi))]


def sigma(n: int, power: int) -> int:
    return sum(d**power for d in range(1, n+1) if n % d == 0)


def eta_product_coeffs(alpha: F, degree: int) -> list[F]:
    """Coefficients of prod_(n>=1)(1-q^n)^alpha, without q^(alpha/24)."""
    coeffs = [F(1)] + [F(0)]*degree
    for n in range(1, degree+1):
        fac = [F(1)]
        for k in range(1, degree//n+1):
            fac.append(-fac[-1]*(alpha-k+1)/k)
        new = [F(0)]*(degree+1)
        for i in range(degree+1):
            for k, a in enumerate(fac):
                if i+n*k <= degree:
                    new[i+n*k] += coeffs[i]*a
        coeffs = new
    return coeffs


def eta_ward_coeffs(h: F, degree: int) -> list[F]:
    """(q d/dq - h E2/12) q^(h/12) sum a_n q^n = 0."""
    a = [F(1)]
    for n in range(1, degree+1):
        a.append(-2*h*sum(F(sigma(k, 1))*a[n-k] for k in range(1, n+1))/n)
    return a


def eisenstein(tau, weight: int = 4, terms: int = 140):
    import mpmath as mp
    if mp.im(tau) <= 0 or weight not in (2, 4, 6):
        raise ValueError('Require Im(tau)>0 and weight in {2,4,6}.')
    q = mp.exp(2j*mp.pi*tau)
    factors = {2: -24, 4: 240, 6: -504}
    return 1 + factors[weight]*mp.fsum(sigma(n, weight-1)*q**n for n in range(1, terms+1))


RANK_INPUT = {
  5: {
    'rank0': [1,25,300,2300,12650,53120,176900,478700,1068575,1982350,3054880,3869650,3931075,3067350,1723100,639850,141575,15900,550,0,0,0,0,0,0,0],
    'rank1': [0,0,0,0,0,10,200,2000,13000,60600,213280,581350,1229100,1970750,2301450,1861060,994350,340550,72700,9000,520,0,0,0,0,0],
    'rank2': [0,0,0,0,0,0,0,0,0,25,600,6400,40125,162200,432850,767850,907050,725125,407450,168100,52610,12650,2300,300,25,1],
    'cylinder_root': '0.5922358232050258',
  },
  6: {
    'rank0': [1,36,630,7140,58905,376992,1947780,8347320,30254904,94088944,253787076,598517172,1241138298,2270804400,3669278472,5226303348,6525511038,7070366736,6544591544,5064164712,3186853893,1579450836,595366164,164526840,31987068,4124196,313290,10656,72,0,0,0,0,0,0,0,0],
    'rank1': [0,0,0,0,0,0,12,360,5436,54336,399780,2288088,10538070,39962304,126774180,339782112,772416936,1485756360,2397393616,3199280040,3466778256,2990089320,2015358768,1046966400,415851744,125773992,28780308,4889208,589680,45576,1704,0,0,0,0,0,0],
    'rank2': [0,0,0,0,0,0,0,0,0,0,0,36,1332,22896,244548,1817100,9944136,41373504,133150140,334051848,654239961,998362404,1185572268,1099296360,803838888,470907108,225093258,89243416,29670588,8302104,1946088,376992,58905,7140,630,36,1],
    'cylinder_root': '0.592507356205638',
  },
}


def rank_probability(counts: list[int], p):
    import mpmath as mp
    n = len(counts)-1
    return mp.fsum(c*p**k*(1-p)**(n-k) for k, c in enumerate(counts) if c)


def root_diagnostics() -> list[dict]:
    import mpmath as mp
    pc = mp.mpf('0.59274605079210')
    target = mp.re(eisenstein(1j))
    records = []
    for size, data in RANK_INPUT.items():
        n = size**2
        good = all(sum(data[f'rank{j}'][k] for j in range(3)) == comb(n,k) for k in range(n+1))
        if not good:
            raise ValueError(f'Cardinality normalization failed for L={size}')
        f = lambda p: rank_probability(data['rank2'], p)-rank_probability(data['rank0'], p)
        root = mp.findroot(f, (mp.mpf('.58'), mp.mpf('.61')))
        cylinder = mp.mpf(data['cylinder_root'])
        ratio = (root-pc)/(cylinder-pc)
        records.append({'L': size, 'cardinality_sum_check': good,
                        'torus_root': mp.nstr(root, 36),
                        'cylinder_root_input': data['cylinder_root'],
                        'shift_ratio': mp.nstr(ratio, 24),
                        'relative_deviation_from_E4_i': mp.nstr(ratio/target-1, 24),
                        'ratio_pc_derivative': mp.nstr((root-cylinder)/(cylinder-pc)**2, 24),
                        'normalization_only_not_rank_assignment_audit': True})
    return records


def report(dps: int = 70) -> dict:
    import mpmath as mp
    if dps < 45:
        raise ValueError('Use at least 45 digits for reproducible controls.')
    mp.mp.dps = dps
    chi, q, qhat = descendants()
    basis = list(partitions(4))
    sub = null_derivative_subspace()
    reduced = add(qhat, {(4,): -F(1)})
    assert not act(1, chi) and not act(2, chi)
    assert not act(1, qhat)
    assert add(act(1,q), scale(-F(80,11),act(-1,chi))) == {}
    assert rational_rank(sub,basis) == rational_rank(sub+[reduced],basis) == 4
    tau = mp.mpc('.21','1.17')
    e = eisenstein(tau)
    rho = mp.mpc(F(1,2).numerator/F(1,2).denominator,mp.sqrt(3)/2)
    serialize = lambda v: {','.join(map(str,w)): str(a) for w,a in sorted(v.items())}
    return {
        'schema': 'thermal-null-ward-modular-root-v1',
        'scope': 'Exact conditional Virasoro algebra; numerical modular functions; reanalysis of supplied rank counts. No site scaling theorem or new large enumeration.',
        'repo_snapshot': '8e1282f6d4397f03c80a2883b917c2d31a4b93a3',
        'exact': {
            'thermal_h': str(H), 'chi': serialize(chi),
            'normalized_quasiprimary': serialize(qhat),
            'L1_quasiprimary': serialize(act(1,qhat)),
            'quotient_and_derivative_rank_at_level4': rational_rank(sub,basis),
            'normalized_class': 'L_-4 modulo null descendants and total derivatives',
            'retained_null_correction_coefficient': '80/87',
            'retained_null_form': 'U4 = L_-4 epsilon + (80/87)L_-2 chi modulo total derivatives',
            'null_diagonal_magnetic_weight': str(F(2,3)*H*(H+1)-H),
            'cylinder_Lminus2_coefficient_at_magnetic_weight': str(F(5,96)-H/12),
            'cylinder_chiral_Lminus4_coefficient': str(H/240),
            'torus_chiral_G4_coefficient': str(3*H),
            'torus_chiral_pi4_E4_coefficient': str(H/15),
            'eta_without_q_leading_coeffs': [str(x) for x in eta_product_coeffs(2*H,10)],
            'eta_product_matches_null_ward_ODE': eta_product_coeffs(2*H,10) == eta_ward_coeffs(H,10),
            'first_trace_coefficient_from_primary_3point': str((2*F(5,96)+H*(H-1))/(2*F(5,96))),
        },
        'modular_numerical': {
            'E4_i': mp.nstr(mp.re(eisenstein(1j)),36),
            'E4_2i': mp.nstr(mp.re(eisenstein(2j)),36),
            'hex_E4_abs': mp.nstr(abs(eisenstein(rho)),8),
            'hex_E4_derivative': mp.nstr(-(2j*mp.pi/3)*eisenstein(rho,6),30),
            'S_transform_error': mp.nstr(abs(eisenstein(-1/tau)-tau**4*e),8),
            'T_transform_error': mp.nstr(abs(eisenstein(tau+1)-e),8),
            'arithmetic_is_not_outward_interval': True,
        },
        'existing_data_diagnostics': root_diagnostics(),
        'existing_data_provenance': {
            'rank5': 'results/geometric-consistency/rank-sector-C-L5-20260914.json; blob e3b03a7583b5346083907a4531ba51dc4f3d8819',
            'rank6': 'results/geometric-consistency/rank-sector-C-L6-20260914.json; blob 4146a7bd51ca8f2e5c964b5c157ce1b2c4bd1f2e',
            'cylinders': 'results/geometric-consistency/fixed-width-charge-spectrum-derivatives-w4-w8-20260914.json; printed root values; no independent Perron rerun',
            'pc_reference_not_certificate': '0.59274605079210',
        },
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path)
    ap.add_argument('--dps', type=int, default=70)
    args = ap.parse_args()
    text = json.dumps(report(args.dps), ensure_ascii=False, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text,encoding='utf-8')
    else:
        print(text,end='')

if __name__ == '__main__':
    main()
