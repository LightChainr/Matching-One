# P250/P249/P255 bivariate common-state freeze

The scalar x/y character score failed because the axis-difference row is not
closed by the axis-average two-root basis.  The next acquisition therefore
measures a two-dimensional displacement table rather than another scalar
phase.

## Exact observable

Record the full C4-closed Manhattan-radius-four diamond, 41 displacements and
328 real coordinates across plus/minus hands and Z5 charges 1/2.  The stencil
contains `xy`, `x^2 y`, and `x y^2`, plus five total-degree-four holdouts.
All coordinates share one fresh random block and a complete 400-batch
covariance.

The old `(parent,fiber)` section is not C4-equivariant: rotation acts affinely
as `j' = k j + s(x)` with `k=3` for the plus child and `k=2` for the minus
child.  Before collection, the exact gate solves

```text
t(Rx) = s(x) + k t(x) mod 5
```

and records `zeta^(r t(x)) O_r(x)`.  This removes the position-dependent
cocycle.  Testing a constant rotation matrix without this gauge would test a
coordinate artifact.

## Frozen state test

For ranks 1, 2, and 3, fit separate shared complex x/y recurrences using only

```text
(0,0), (1,0), (2,0), (3,0), (0,1), (0,2), (0,3).
```

Pair the x/y roots by the minimum axis-fit residual and retain channel-specific
complex amplitudes.  The commuting common-state hypothesis must first predict
`(1,1),(2,1),(1,2)`, which were not used in the axis fit.  It must then predict
the total-degree-four points `(4,0),(3,1),(2,2),(1,3),(0,4)`.  The first rank
with both covariance scores at `p>=0.01` is the minimal identifiable commuting
dimension.  Failure through rank three is an allowed result.

For `Ty=R Tx R^-1`, score equality of the two characteristic polynomials as a
necessary similarity invariant.  Separately score the signed C4 displacement
orbits under the exact realified channel map, whose fourth power is identity.
Neutral two-point data do not uniquely determine R inside the state
centralizer, so no stronger claim is frozen.  `D^5` is explicitly not
identifiable because the pair cancels simultaneous deck phase.

## Fresh run

Use 80,000 samples, 400 batches, seed `25050510120261130`, counters
`[0,80000)`, and 16 workers on XPk2PZ.  The existing axis-only N505 stream is
used only for this sample-size choice: degree at most three has minimum complex
resolution `z=7.19`; the degree-four boundary is a collective heldout rather
than an individual-row detection gate.
