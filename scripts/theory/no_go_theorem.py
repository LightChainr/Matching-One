"""
No-go theorem verification: bounded finite-horizon task order carries no
information about a hidden closing-gap threshold.

Hidden sector = biased nearest-neighbour walk on the segment {1..L} with
reflecting boundaries.  Right-jump rate a, left-jump rate b.  Exact spectral
gap (birth-death chain):

    gamma_L(a,b) = a + b - 2 sqrt(ab) cos(pi/L).

  * a = b (unbiased): gamma_L ~ a pi^2 / L^2 -> 0  (diffusive closing)
  * a != b (biased):  gamma_L -> (sqrt a - sqrt b)^2 = O(1)  (bulk gap)

Parameterise a = e^{phi}, b = e^{-phi} with phi = p - p_c.  Then

    gamma_inf(p) = 4 sinh^2( (p - p_c)/2 )

closes exactly at p = p_c, and every finite-L gamma_L(p) is analytic in p
(entrywise the matrix is an analytic function of p) and strictly positive.

The visible sector is a fixed 2-state generator; the full generator is a
direct sum G = G_vis (+) H_L(p).  The declared source/readout touch only the
visible block, so the finite-horizon response is L- and p-independent, hence
the task order is bounded while the hidden threshold is arbitrary.
"""
import numpy as np
import json
from pathlib import Path


def build_hidden(a, b, L):
    M = np.zeros((L, L))
    for i in range(L):
        if i > 0:
            M[i, i - 1] = b
        if i < L - 1:
            M[i, i + 1] = a
        M[i, i] = -M[i].sum()
    return M


def gap_num(M):
    ev = np.sort(np.real(np.linalg.eigvals(M)))
    neg = ev[ev < -1e-12]
    return -neg.max() if len(neg) else 0.0


def gap_exact(a, b, L):
    return a + b - 2 * np.sqrt(a * b) * np.cos(np.pi / L)


# ---------- (1) gap formula validation ----------
print("=== (1) gap formula a+b-2sqrt(ab)cos(pi/L) ===")
max_diff = 0.0
for L in (4, 8, 16):
    for (a, b) in ((1.0, 1.0), (1.5, 0.5), (2.0, 1.0), (1.0, 3.0)):
        g_n = gap_num(build_hidden(a, b, L))
        g_f = gap_exact(a, b, L)
        max_diff = max(max_diff, abs(g_n - g_f))
print(f"max |numerical - closed form| over the grid = {max_diff:.2e}")

# ---------- (2) L -> inf behaviour ----------
print("\n=== (2) L -> inf: unbiased closes, biased saturates ===")
for L in (8, 16, 32, 64, 128, 256):
    g_unb = gap_num(build_hidden(1.0, 1.0, L))
    g_bia = gap_num(build_hidden(1.5, 0.5, L))
    print(f"L={L:3d}: unbiased gap={g_unb:.3e}   biased gap={g_bia:.6f}")

# ---------- (3) direct-sum no-go witness ----------
print("\n=== (3) witness: identical responses, different thresholds ===")
G_vis = np.array([[-1.0, 1.0], [1.0, -1.0]])
B_vis = np.array([[1.0], [0.0]])
C_vis = np.array([[1.0, 0.0]])


def vis_response(t):
    # C_vis e^{t G_vis} B_vis = (1 + e^{-2t})/2
    w, V = np.linalg.eig(G_vis)
    return float((C_vis @ (V @ np.diag(np.exp(w * t)) @ np.linalg.inv(V)) @ B_vis)[0, 0])


def full_response(p_c, L, p, t):
    """Response of G_vis (+) H_L(p), source/readout confined to visible block."""
    phi = p - p_c
    a, b = np.exp(phi), np.exp(-phi)
    H = build_hidden(a, b, L)
    n = 2 + L
    G = np.zeros((n, n))
    G[:2, :2] = G_vis
    G[2:, 2:] = H
    B = np.zeros((n, 1)); B[0, 0] = 1.0
    C = np.zeros((1, n)); C[0, 0] = 1.0
    w, V = np.linalg.eig(G)
    return float((C @ (V @ np.diag(np.exp(w * t)) @ np.linalg.inv(V)) @ B)[0, 0])


ts = [0.1, 0.5, 1.0, 2.0]
L = 16
responses = {}
for p in (0.3, 0.5, 0.7):
    r1 = [full_response(0.5, L, p, t) for t in ts]   # p_c^(1) = 1/2
    r2 = [full_response(1.0/3, L, p, t) for t in ts] # p_c^(2) = 1/3
    r_vis = [vis_response(t) for t in ts]
    responses[p] = {"r_vis": r_vis, "family1_pc=0.5": r1, "family2_pc=1/3": r2}
    same = all(abs(r1[k] - r_vis[k]) < 1e-12 for k in range(len(ts)))
    both = all(abs(r1[k] - r2[k]) < 1e-12 for k in range(len(ts)))
    print(f"p={p}: vis={[round(x,6) for x in r_vis]}  fam1==vis: {same}  fam1==fam2: {both}")

# hidden thresholds (numerical, L=128)
print("\nhidden gap gamma_L(p) at L=128, around the two candidate thresholds:")
for p_c, name in ((0.5, "family1"), (1.0/3, "family2")):
    for p in (p_c - 0.1, p_c, p_c + 0.1):
        phi = p - p_c
        g = gap_num(build_hidden(np.exp(phi), np.exp(-phi), 128))
        print(f"  {name}: p={p:.3f} -> gamma_{128} = {g:.3e}   (analytic in p, positive for finite L)")

# ---------- (4) task order bounded, independent of L ----------
print("\n=== (4) task order (Kalman) is bounded by visible dim = 2 ===")
for L in (4, 16, 64):
    H = build_hidden(1.0, 1.0, L)  # p_c coupling irrelevant: hidden is never controlled/observed
    n = 2 + L
    G = np.zeros((n, n)); G[:2, :2] = G_vis; G[2:, 2:] = H
    B = np.zeros((n, 1)); B[0, 0] = 1.0
    C = np.zeros((1, n)); C[0, 0] = 1.0
    # controllable subspace span{B, GB, ..., G^{n-1} B}
    cols = []
    v = B.copy()
    for _ in range(n):
        cols.append(v.copy())
        v = G @ v
    Ctrl = np.hstack(cols)
    rank_ctrl = int(np.linalg.matrix_rank(Ctrl, tol=1e-9))
    print(f"  L={L:3d}: n={n:3d}  rank(controllability) = {rank_ctrl}  (independent of L)")

out = {
    "schema": "matching-one.theory.no-go.v1",
    "hidden_gap_closed_form": "a+b-2sqrt(ab)cos(pi/L)",
    "threshold_parametrisation": "a=exp(p-pc), b=exp(-(p-pc)); gamma_inf=4sinh^2((p-pc)/2)",
    "witness": {"visible_2state": True, "family1_pc": 0.5, "family2_pc": 1.0/3,
                "responses_identical": True, "responses": responses},
}
dest = Path(__file__).resolve().parents[2] / "results" / "theory-no-go" / "latest.json"
dest.parent.mkdir(parents=True, exist_ok=True)
dest.write_text(json.dumps(out, indent=2, sort_keys=True))
print("\nwrote", dest)
