# Two essential homology births in the threshold-rank filtration

Status: tiny exact semantics gate for Issue 269. No production stream or continuum claim is changed.

For a site permutation, let `B_k` contain its first `k` sites and let

```text
R_k = rank im[H1(B_k) -> H1(T^2)].
```

The oracle recomputes `R_k` for every `k`, rather than trusting the names of the historical endpoint
variables. It separately evaluates the old forward primal rank-two birth and reverse matching
rank-two birth, including the conversion `K_minus=N-r+1`.

## Exact gates

Every permutation of axis `L=2` and Gaussian `(2,1)` is exhausted: `4!+5!=144` paths. Each path must
satisfy

```text
K_minus = min{k : R_k >= 1},
K_plus  = min{k : R_k = 2},
R_k = 1[k >= K_minus] + 1[k >= K_plus].
```

After swapping the primal and matching graphs and reversing the permutation, the two endpoint gates
are

```text
K_minus^G(pi) + K_plus^Ghat(reverse(pi)) = N+1,
K_plus^G(pi)  + K_minus^Ghat(reverse(pi)) = N+1.
```

Thus the existing `K_minus/K_plus` convention is directly connected to the first and second ambient
homology births on these exact controls, including simultaneous births when the rank-one plateau is
empty.

## Projective line and integral index

At every step with `R_k=1`, the oracle gathers the actual lifted winding generators. It reduces them
to a canonical primitive line in `P^1(Q)` and separately computes the gcd of their integral
coefficients, namely the saturation index of the image subgroup inside that line.

The black NN and complementary white matching lines must agree, have zero determinant, and remain
constant throughout one rank-one plateau. The integral black/white indices are recorded as a pair at
every step; rational-line equality is never promoted silently to equality of integer subgroups.

## Boundary

This exhaustive result covers only the two declared tiny quotients. The merged surface proof supplies
the regular-torus theorem, but production sufficient statistics, general-period regression, and the
self-identifying short-period boundary remain separate work. The filtration variables are exact
topological coordinates, not by themselves local CFT fields.
