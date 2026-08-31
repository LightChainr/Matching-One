# When a descendant preserves Jordan rank

Issue #216 uses the nonzero thermal `Q4` descendant to inherit a rank-two
Jordan pair.  This certificate isolates the exact linear-algebra condition
behind that step and records a negative control that prevents a stronger,
incorrect claim.

## Exact criterion

Write a rank-`r` source chain in bottom-to-top order as

```text
N e_0 = 0,        N e_j = e_(j-1).
```

For a homogeneous level-`ell` map `A`, the descendant relation is

```text
D_target A = A D_source + ell A,
```

or equivalently `N_target A = A N_source` after subtracting the shifted
eigenvalues.  Applying `N_target^(r-1)` to the image of the top vector gives

```text
N_target^(r-1) A e_(r-1) = A e_0.
```

Therefore an equal-rank target contains the full image chain exactly when
the bottom image `A e_0` is nonzero.  Intertwining by itself is not enough.

## Thermal Q4 instance

For the rank-two thermal pair, `D_source=5/4 I+N` and homogeneous left level
four shifts the descendant eigenvalue to `21/4`.  The common `Q4` action is
the identity on the two-dimensional Jordan label, so the intertwining
identity is exact and the bottom survives.  The independently derived
ordinary Gram norm `4930` supplies the required non-null check.  The image
therefore remains rank two.

## Commuting collapse control

On a rank-three chain choose `A=N`.  It commutes with `N` and satisfies the
same shifted intertwining identity, but it kills `e_0`.  The image of the top
has only the two nonzero steps `e_1,e_0`; the image rank is two, not three.
This explicitly rules out “commuting descendants always preserve rank.”

## Reproduction

```text
python3 scripts/descendant_jordan_rank_survival.py
python3 -m unittest tests/test_descendant_jordan_rank_survival.py -v
```

## Boundary

This exact finite-chain statement does not establish lattice overlap, fix a
logarithmic coefficient, derive the logarithmic torus Ward response, or show
that `P4[S']` reads out the top component.  Those parts of Issue #216 remain
open.
