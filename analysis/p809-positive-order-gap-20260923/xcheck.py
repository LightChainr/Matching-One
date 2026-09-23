#!/usr/bin/env python3
"""Control-side independent cross-check of the #809 positive-gap note.

Does not import or copy check.py. Rebuilds the D4 quotient from the pinned PR708
certificate with its own orbit code, computes exact observability ranks modulo
primes that were NOT used by check.py, and re-derives the invisible vector,
the core structure and the two spectral facts numerically at 50 digits.
Numerics here are corroboration only; the proof certificate is check.py.

Run: python xcheck.py input/width4-rank-closure-certificate.json certificate.json
"""
import json, sys
from collections import deque
import mpmath as mp

mp.mp.dps = 50
cert = json.load(open(sys.argv[1])); seed = json.load(open(sys.argv[2]))
q, init, rk = cert['quotient_transitions'], cert['quotient_initial'], cert['quotient_rank_output']
n0 = len(q)

# 1. D4 orbits from histories (own implementation)
hist = {}
for a in range(16):
    for b in range(16):
        hist.setdefault(q[init[a]][b], [a, b])
todo = deque(hist)
while todo:
    s = todo.popleft()
    for b in range(16):
        t = q[s][b]
        if t not in hist:
            hist[t] = hist[s] + [b]; todo.append(t)
assert len(hist) == n0

def run(word):
    s = init[word[0]]
    for b in word[1:]:
        s = q[s][b]
    return s

def act(mask, refl, rot):
    return sum(1 << (((-c if refl else c) + rot) % 4) for c in range(4) if mask >> c & 1)

orbit_min = []
for s in range(n0):
    imgs = [run([act(m, r, t) for m in hist[s]]) for r in (0, 1) for t in range(4)]
    assert all(rk[i] == rk[s] for i in imgs)
    orbit_min.append(min(imgs))
order = list(dict.fromkeys(orbit_min))
cls = [order.index(o) for o in orbit_min]
N = len(order); rep = [cls.index(c) for c in range(N)]
bits = [int(rk[rep[c]] == 1) for c in range(N)]

def profile(s):
    d = {}
    for b in range(16):
        key = (bin(b).count('1'), cls[q[s][b]]); d[key] = d.get(key, 0) + 1
    return d
assert all(profile(s) == profile(rep[cls[s]]) for s in range(n0))
prof = [profile(rep[c]) for c in range(N)]
print(f'D4 classes: {N}; strongly lumpable for every p')

# 2. exact observability rank modulo primes (roots of f, and generic points)
f = seed['invisible']['factor']
def Kmod(a, m):
    inv = pow(pow(1 + a, 4, m), -1, m); K = [[0] * N for _ in range(N)]
    for i in range(N):
        for (k, j), mult in prof[i].items():
            K[i][j] = (K[i][j] + mult * pow(a, k, m) * inv) % m
    return K
def obs_rank(K, m):
    basis = {}
    def add(v):
        v = v[:]
        for p, r in basis.items():
            if v[p]:
                c = v[p]; v = [(x - c * y) % m for x, y in zip(v, r)]
        p = next((i for i, x in enumerate(v) if x), None)
        if p is None:
            return False
        inv = pow(v[p], -1, m); v = [x * inv % m for x in v]
        for pp in basis:
            if basis[pp][p]:
                c = basis[pp][p]; basis[pp] = [(x - c * y) % m for x, y in zip(basis[pp], v)]
        basis[p] = v; return True
    work = deque()
    for b in (0, 1):
        v = [int(bits[i] == b) for i in range(N)]
        if add(v): work.append(v)
    while work:
        u = work.popleft(); Ku = [sum(K[i][j] * u[j] for j in range(N)) % m for i in range(N)]
        for b in (0, 1):
            v = [Ku[i] if bits[i] == b else 0 for i in range(N)]
            if add(v): work.append(v)
    return len(basis)
found = 0
for m in range(1500, 4000):
    if any(m % d == 0 for d in range(2, int(m ** .5) + 1)):
        continue
    roots = [a for a in range(m) if sum(c * pow(a, k, m) for k, c in enumerate(f)) % m == 0 and (a + 1) % m]
    if roots:
        print(f'rank at root of f mod {m} (x={roots[0]}): {obs_rank(Kmod(roots[0], m), m)}'); found += 1
    if found == 3:
        break
for m, a in ((2003, 5), (3001, 777)):
    print(f'rank at generic x={a} mod {m}: {obs_rank(Kmod(a, m), m)}')

# 3. numerics at p_a
xa = mp.findroot(lambda z: sum(c * z ** k for k, c in enumerate(f)), mp.mpf('1.1557')); pa = xa / (1 + xa)
print('p_a =', mp.nstr(pa, 25))
K = mp.matrix(N, N)
for i in range(N):
    for (k, j), mult in prof[i].items():
        K[i, j] += mult * pa ** k * (1 - pa) ** (4 - k)
basis = []
def gs(v):
    v = mp.matrix(v)
    for _ in range(2):
        for bv in basis:
            v = v - (bv.T * v)[0] * bv
    nv = mp.norm(v)
    return None if nv < mp.mpf(10) ** -35 else v / nv
work = []
for b in (0, 1):
    u = gs([int(bits[i] == b) for i in range(N)])
    if u is not None: basis.append(u); work.append(u)
while work:
    u = work.pop(0); Ku = K * u
    for b in (0, 1):
        w = gs([Ku[i] if bits[i] == b else 0 for i in range(N)])
        if w is not None: basis.append(w); work.append(w)
print('numeric Krylov dimension (tol 1e-35):', len(basis))
rec = seed['invisible']; v = [mp.mpf(0)] * N
for i, cs in zip(rec['support'], rec['vectors']):
    v[i] = sum(mp.mpf(c) * xa ** k for k, c in enumerate(cs))
k0 = max(range(N), key=lambda i: abs(v[i]))
e = mp.matrix(N, 1); e[k0] = 1
Bm = mp.matrix(N, len(basis))
for c, bv in enumerate(basis):
    for i in range(N): Bm[i, c] = bv[i]
u = e - Bm * (Bm.T * e)
print('independent null vector vs certificate v: max|u - c v| =',
      mp.nstr(max(abs(u[i] - u[k0] / v[k0] * v[i]) for i in range(N)), 5))
lam = sum(mp.mpf(c) * xa ** k for k, c in enumerate(rec['eigenvalue'])) / (1 + xa) ** 4
vK = mp.matrix(v).T * K
print('lambda =', mp.nstr(lam, 15), ' |vK - lambda v| =', mp.nstr(max(abs(vK[0, j] - lam * v[j]) for j in range(N)), 5))
G0, A, Br, J = seed['good_zero'], seed['absorber'], seed['bridge_one'], seed['resonant_zero']; G = G0 + [A]
nz = lambda i: {j for j in range(N) if K[i, j] != 0}
print('support(v) in G0uJ:', {i for i in range(N) if abs(v[i]) > mp.mpf(10) ** -40} <= set(G0) | set(J),
      '| G closed:', all(nz(i) <= set(G) for i in G), '| absorber:', nz(A) == {A},
      '| J exits in J+{absorber,bridge}:', all(nz(i) <= set(J) | {A, Br} for i in J))
R = mp.matrix([[K[i, j] for j in G0] for i in G0]); Q = mp.matrix([[K[i, j] for j in J] for i in J])
eR = mp.eig(R, left=False, right=False)
print('eig(R):', sorted(mp.nstr(z.real, 8) + ('' if abs(z.imag) < 1e-40 else '%+.3gi' % float(z.imag)) for z in eR))
print('real eig(R) in [1/50,19/25]:', [mp.nstr(z.real, 8) for z in eR if abs(z.imag) < 1e-40 and 0.02 <= z.real <= 0.76])
print('Perron(Q) =', mp.nstr(max(abs(z) for z in mp.eig(Q, left=False, right=False)), 10))
ell = seed['left_supervector']
print('min margin (19/25)l - lQ =', mp.nstr(min(mp.mpf(19) / 25 * ell[j] - sum(ell[i] * Q[i, j] for i in range(10)) for j in range(10)), 6))
