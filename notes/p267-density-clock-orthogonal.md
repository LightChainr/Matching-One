# P267 density-clock orthogonalization

## Why this decomposition is identifiable from the stored paths

For a uniform K-site subset of a translation-regular square quotient, the
external Euler coordinate `O_ext=V-E_NN+F0` has exact conditional mean

\[
\mu_N(k)=k-2N(k)_2/(N)_2+N(k)_4/(N)_4.
\]

The path archive records, for every `k/batch/orientation`, `O_ext`, `O_ext^2`,
both complex `J_D/J_S`, both `O_ext*J` products, `J_D*conj(J_S)`, and
`|J_S|^2`.  Therefore it supports all of the following without returning to
the sparse per-replica stream:

\[
\begin{aligned}
O_0 &= O_{ext}-\mu_N(K),\\
D_0 &= J_D-E[J_D\mid K],\\
S_0 &= J_S-E[J_S\mid K],\\
\beta_0 &= \frac{\operatorname{Re}E[D_0\overline{S_0}]}{E|S_0|^2}.
\end{aligned}
\]

At fixed K the cross product can be evaluated in the frozen order,

\[
E[(O-\mu_K)(J-E[J\mid K])\mid K]
=E[OJ\mid K]-E[O\mid K]E[J\mid K].
\]

The explicit `mu_K` cancels algebraically after the source conditional mean is
also removed.  This is not a failure of the clock definition: it is exactly the
within-layer term in the law of total covariance, while the removed piece is
the between-layer density clock.  The implementation records the expanded and
simplified forms as an equality audit.

## Frozen execution

- retrospective protocol: `fe26a8f`;
- pre-reveal scorer/exact test: `285229b`;
- inputs: Target 1 2M and two-observer 2M at N325/N425;
- all four seed/counter domains are disjoint;
- full and every delete-one estimate recompute the intrinsic root and both
  source projections;
- observer restricted to `O_ext`; `O_far` is not scored.

The fixed-K residual retains 0.321–0.323 of the raw magnitude in all four
independent size/block cells.  Its N425/N325 complex transfer is reproduced:
amplitude 1.05183 versus 1.05233 and phase 1.37369 versus 1.37273 radians.
Independent-block compatibility is p=0.943 at N325 and p=0.995 at N425.

This closes the “pure density clock” mechanism: it is quantitatively dominant
but incomplete.  The surviving one-third response is a stable within-density
coupling and is now the higher-information object for a subsequent field
selector.

