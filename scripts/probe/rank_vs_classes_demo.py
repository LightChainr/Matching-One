#!/usr/bin/env python3
"""Predictive class count vs response-matrix rank: two elementary families.

Family 1 (the #549 obstruction): states a = 0..k; every single-round
shared-update fork test has success probability affine in a.  The response
matrix for ANY finite set of such tests has rank <= 2 while the class count is
k+1.

Family 2 (Vandermonde): if a language contains tests whose probabilities are
theta^j for j = 0..k, the response matrix on k+1 distinct states is
Vandermonde and has rank k+1.

Nothing here is Monte Carlo.  Run: python3 scripts/probe/rank_vs_classes_demo.py
"""

import numpy as np


def fork_probability(k: int, a: int) -> float:
    """#549 closed form: [343k^3 - 182k^2 + 25k + 4a] / [8k(8k-1)^2]."""
    return (343 * k ** 3 - 182 * k ** 2 + 25 * k + 4 * a) / (8 * k * (8 * k - 1) ** 2)


def affine_family(k: int, n_tests: int, seed: int = 0) -> np.ndarray:
    """k+1 states (a = 0..k) x n_tests arbitrary affine-in-a probabilities."""
    rng = np.random.default_rng(seed)
    slopes = rng.uniform(-0.2, 0.2, size=n_tests)
    intercepts = rng.uniform(0.2, 0.6, size=n_tests)
    a = np.arange(k + 1, dtype=float)
    return intercepts[None, :] + slopes[None, :] * a[:, None]


def vandermonde_family(k: int) -> np.ndarray:
    """k+1 distinct states theta_j, tests giving theta^j (j = 0..k)."""
    theta = np.linspace(0.5, 2.0, k + 1)
    return np.vander(theta, N=k + 1, increasing=True)


def main() -> None:
    for k in (1, 4, 9, 24):
        M = affine_family(k, n_tests=7)
        r = np.linalg.matrix_rank(M, tol=1e-9)
        print(f"k={k:>2}  classes={k+1:>3}  affine-family rank = {r}   (predicted <= 2)")
    print()
    for k in (1, 4, 9, 24):
        M = vandermonde_family(k)
        r = np.linalg.matrix_rank(M, tol=1e-9)
        print(f"k={k:>2}  classes={k+1:>3}  Vandermonde rank = {r}   (predicted = k+1)")
    # exact reproduction of the #549 fork probabilities themselves
    print()
    for k in (1, 4, 9):
        vals = [fork_probability(k, a) for a in range(k + 1)]
        diffs = np.diff(vals)
        print(f"k={k}: F_k,a = {np.round(vals, 6)}, constant increment = {np.allclose(diffs, diffs[0])} "
              f"(value {diffs[0]:.6e})")
    print()
    print("CONCLUSION: class count and response rank are different coordinates;")
    print("affine-in-label fork languages cap rank at 2 regardless of k.")


if __name__ == "__main__":
    main()
