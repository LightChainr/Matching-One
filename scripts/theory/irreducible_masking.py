"""
Irreducible exact-masking: on a connected, positive, nearest-neighbour walk the
slowest mode can be made EXACTLY invisible by choosing the readout C and source
B orthogonal to the slow right/left eigenvectors.  This is fine-tuning, not
generic — the partial-fraction coefficient of e^{lambda1 t} in C e^{tG} B is
(c^T v)(w^T b), which vanishes iff c^T v = 0 or w^T b = 0.
"""
import numpy as np
import json
from pathlib import Path


def build(a, b, L):
    M = np.zeros((L, L))
    for i in range(L):
        if i > 0:
            M[i, i - 1] = b
        if i < L - 1:
            M[i, i + 1] = a
        M[i, i] = -M[i].sum()
    return M


def main():
    L, a, b = 4, 1.5, 0.5
    G = build(a, b, L)
    assert np.allclose(G.sum(axis=1), 0)

    evals, V = np.linalg.eig(G)
    order = np.argsort(-np.real(evals))     # 0 first, then slowest, ...
    evals, V = evals[order], V[:, order]
    lam1 = evals[1]
    v = V[:, 1]                              # right slow eigenvector
    _, W = np.linalg.eig(G.T)
    W = W[:, np.argsort(-np.real(np.linalg.eigvals(G.T)))]
    w = W[:, 1]
    w = w / (w @ v)                          # dual normalisation w^T v = 1

    def slow_coeff(C, B):
        return float((C @ v) * (w @ B))

    B_gen = np.array([1., 0, 0, 0])
    C_gen = np.array([1., 0, 0, 0])
    C_ft = np.array([1.0, -v[0] / v[1], 0.0, 0.0])
    C_ft /= np.linalg.norm(C_ft)
    B_ft = np.array([1.0, -w[0] / w[1], 0.0, 0.0])

    out = {
        "schema": "matching-one.theory.irreducible-masking.v1",
        "eigenvalues": [float(np.real(x)) for x in evals],
        "slowest_mode": float(np.real(lam1)),
        "partial_fraction_coeff": {
            "generic": slow_coeff(C_gen, B_gen),
            "masked_readout_Cv0": slow_coeff(C_ft, B_gen),
            "masked_source_wTB0": slow_coeff(C_gen, B_ft),
            "masked_both": slow_coeff(C_ft, B_ft),
        },
        "claim": ("on a connected positive walk the slowest mode is exactly masked "
                  "iff the readout kills the right slow eigvec or the source kills "
                  "the left slow eigvec; generic readout/source see it (0.48)."),
    }
    dest = Path(__file__).resolve().parents[2] / "results" / "theory-irreducible-masking" / "latest.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2, sort_keys=True))
    print(json.dumps(out, indent=2))
    print("wrote", dest)
    return out


if __name__ == "__main__":
    main()
