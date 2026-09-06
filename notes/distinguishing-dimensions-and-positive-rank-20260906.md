# Programs F + E: distinguishing dimensions, and why ordinary rank does not separate here

## Program F — minimal distinguishing experiment sets (exact on #549)

Three genuinely different "witness dimensions" for the `k+1` predictive
classes `a = 0..k` of the #549 family, measured exactly:

| notion | definition | value |
|---|---|---|
| separating dimension | fewest experiments whose response *values* distinguish all classes | **1** (one fork: `F_{k,a}` is strictly increasing in `a`) |
| affine-separating dimension | dimension of the affine span of the classes' response vectors | **2** (`{1, a}`) |
| linear-witness dimension | fewest experiments making the response vectors *linearly independent* | **k+1** (the grouped-fork columns of the rank theorem) |

**COUNTEREXAMPLE to the one-witness fallacy.**  One experiment assigns `k+1`
distinct real values (`m_distinguish = 1`) yet the linear rank is only 2; you
need `k+1` experiments for linear independence.  So "one witness per class ⇒
rank ≥ #classes" is false, and the probe now has an exact canonical instance
of how far apart the three numbers sit: `(1, 2, k+1)`.

The correct general inequality is only

```
separating <= affine-separating <= linear-witness = rank of the response space,
```

with the #549 family attaining the maximal gap between the first two and the
third.

## Program E — positive vs signed realization on the #549 response space

For the grouped-fork response matrix `M` of depth `<= d` (columns spanning
`{1,a,...,a^d}`), the matrix is `(k+1) x (d+1)`, full column rank `d+1`, and
has exactly `d+1` columns — so its nonnegative rank satisfies

```
d+1 = rank(M) <= r_+(M) <= #columns = d+1  ==>  r_+(M) = d+1.
```

**Observation (no separation in this family).**  In the #549 grouped-fork
algebra, nonnegative rank equals ordinary rank at every depth.  The richness
is in the *language depth*, not in a positive-vs-signed gap.  This is the
opposite of P398 (where signed/balanced order is small and exact positive
quotients are large): the two families exhibit different phenomena, exactly as
the program brief warned against assuming.

**Abstract separation exists but is not repository-native yet (OPEN).**  The
classical slack-matrix construction (e.g. the facet–vertex slack matrix of a
regular `n`-gon) has ordinary rank <= 3 while nonnegative rank is `Theta(n)`
(Fiorini–Rothvoß–Tiwari, extension complexity; the polygon bound via rectangle
covers).  Whether a *branching-reliability* family inside the cut-network
algebra can reproduce an analogous ordinary-rank-bounded / nonnegative-rank-
unbounded separation is the open target: it needs a family whose response
entries form a slack-type matrix (nonnegative, with a hidden low-rank signed
realisation but no low-rank nonnegative one).

Files: values above follow from `scripts/rank_lower_bound.py` and
`results/rank-vs-depth/latest.json`.
