"""#741: exact/certified evaluation of the winding-component density nu_w(p).

nu_w = pi . g on the reward-preserving lumped chain produced by winding_build.

Small chains  : full exact-rational stationary solve, verified pointwise.
Large chains  : float64 dense solve for the stationary vector, then the
                stationary residual and nu are RE-EVALUATED in exact rational
                arithmetic on the float64 vector (a dyadic rational, rescaled to
                sum exactly 1). The ticket's certificate then applies verbatim:

    |pi_hat . g - nu| <= ||g||_inf * ||pi_hat K - pi_hat||_1 / delta,
    delta = (1-p)^w  (the uniform empty-row reset probability).

log nu is reported together with that bound.
"""
import json, sys, math
from fractions import Fraction

def load_counts(d, p):
    """Return (sparse K as list of dicts, g as list of Fractions)."""
    w = d["width"]; q = 1 - p; n = len(d["rows"])
    K = []
    g = [Fraction(0)] * n
    for i, row in enumerate(d["rows"]):
        acc = {}
        gi = Fraction(0)
        for m, entries in enumerate(row):
            wm = (p ** m) * (q ** (w - m))
            if wm == 0:
                continue
            for nb, rw, c in entries:
                acc[nb] = acc.get(nb, Fraction(0)) + wm * c
                gi += wm * c * rw
        K.append(acc); g[i] = gi
    return K, g

def pi_exact_dense(K, g, n):
    A = [[K[j].get(i, Fraction(0)) - (1 if i == j else 0) for j in range(n)] for i in range(n)]
    A[-1] = [Fraction(1)] * n
    b = [Fraction(0)] * (n - 1) + [Fraction(1)]
    M = [A[i][:] + [b[i]] for i in range(n)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * bb for a, bb in zip(M[r], M[c])]
    return [M[i][n] for i in range(n)]

def pi_float_dense(K, n, refine=2):
    """float64 stationary solve plus iterative refinement of the stationary residual."""
    import numpy as np
    A = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j, v in K[i].items():
            A[j, i] = float(v)          # row-stochastic K^T
        A[i, i] -= 1.0
    A[-1, :] = 1.0
    b = np.zeros(n); b[-1] = 1.0
    x = np.linalg.solve(A, b)
    for _ in range(refine):
        # residual r = x K - x  (normalisation row kept at 0 for the correction)
        r = np.zeros(n)
        for i in range(n):
            xi = x[i]
            for j, v in K[i].items():
                r[j] += xi * float(v)
        r -= x
        r[-1] = 0.0                     # do not disturb the normalisation equation
        d = np.linalg.solve(A, -r)
        x = x + d
    return np.maximum(x, 0.0)

def certify(K, g, pi, w, p):
    n = len(pi)
    q = 1 - p; delta = q ** w
    S = sum(pi)
    pi = [x / S for x in pi]                       # exact normalisation
    assert min(pi) >= 0
    resid = Fraction(0)
    for j in range(n):
        v = sum((pi[i] * K[i].get(j, Fraction(0)) for i in range(n)), Fraction(0))
        resid += abs(v - pi[j])
    nu = sum((pi[i] * g[i] for i in range(n)), Fraction(0))
    gmax = max(abs(x) for x in g)
    bound = gmax * resid / delta
    return nu, resid, bound, gmax, delta

def run(path, p, exact=True):
    d = json.load(open(path))
    w = d["width"]; n = len(d["rows"])
    K, g = load_counts(d, p)
    if exact:
        pi = pi_exact_dense(K, g, n); mode = "exact-rational solve"
    else:
        pf = pi_float_dense(K, n, refine=0)
        pi = [Fraction(x) for x in pf]             # float64 values are dyadic
        # one correction step driven by the EXACT rational stationary residual
        r = [sum((pi[i] * K[i].get(j, Fraction(0)) for i in range(n)), Fraction(0)) - pi[j]
             for j in range(n)]
        import numpy as np
        A = np.zeros((n, n)); 
        for i in range(n):
            for j, v in K[i].items():
                A[j, i] = float(v)
            A[i, i] -= 1.0
        A[-1, :] = 1.0
        corr = np.linalg.solve(A, np.array([-float(x) for x in r]))
        pi = [Fraction(x) for x in (np.maximum(np.array([float(x) for x in pi]) + corr, 0.0))]
        mode = "float64 solve + exact-rational residual correction + exact certification"
    nu, resid, bound, gmax, delta = certify(K, g, pi, w, p)
    lo = float(nu - bound); hi = float(nu + bound)
    return {"width": w, "matching": d["matching"], "p": str(p), "mode": mode,
            "frontier_states": d["states"], "reward_lumps": n,
            "nu": str(nu), "nu_float": float(nu),
            "log_nu": math.log(float(nu)),
            "certificate_bound": float(bound),
            "log_nu_abs_error_bound": (abs(math.log(hi) - math.log(lo)) / 2 if nu - bound > 0 else None),
            "delta_empty_row_reset": float(delta), "g_inf": float(gmax)}

if __name__ == "__main__":
    paths = json.loads(sys.argv[1])     # list of [path, "1/4", exact_bool]
    out = []
    for path, ps, exact in paths:
        num, den = ps.split("/")
        out.append(run(path, Fraction(int(num), int(den)), exact))
        r = out[-1]
        print(f"w={r['width']:<3} matching={str(r['matching']):<6} p={r['p']:<4} "
              f"states={r['frontier_states']:<7} lumps={r['reward_lumps']:<5} "
              f"log nu = {r['log_nu']:.12f}  bound~{r['certificate_bound']:.2e}  [{r['mode']}]",
              flush=True)
    json.dump(out, open(sys.argv[2], "w"), indent=2)
