# Digital Alexander duality for the square-site matching pair

Status: proof for honest periodic square-cell tori.  Short-period quotient
degeneracies remain covered by the separate finite oracle.

## 1. The two embedded objects

Let `S=T^2` carry its periodic square-cell decomposition.  For a black vertex
set `B`, let `G_B` be the induced nearest-neighbour graph and let `U` be a
closed regular neighbourhood of `G_B`.  Then `U` deformation retracts to
`G_B`, so both have the same image in `H_1(S;Q)`.

The white matching graph contains white NN edges and both possible diagonals
of every square.  As an abstract graph it need not be embedded: when all four
corners are white its diagonals cross.  For ambient homology, however, it has
a canonical embedded reduction.

In each face retain a white diagonal only when its endpoints are exactly the
two white corners.  In every other case where a matching diagonal is present,
its endpoints are connected by a white NN path along the boundary of that
same face.  Replacing the diagonal by this path changes a cycle by a chain
contained in one contractible face and therefore does not change its class in
`H_1(S;Q)`.

The executable certificate checks all 16 face patterns.  Exactly the two
opposite-white-pair patterns retain one diagonal; six redundant active
diagonals receive explicit boundary-path replacements; no pattern retains
crossing diagonals.

Call the resulting embedded graph `G_W`.  Cell by cell, `G_W` is a 1-skeleton
for `V=closure(S\U)`.  The inclusion of a CW 1-skeleton surjects onto the
first homology of the CW complex, while the extra local cycles of the full
matching graph bound chains inside individual faces.  Consequently

```text
im[H_1(white matching graph) -> H_1(S)]
= im[H_1(G_W) -> H_1(S)]
= im[H_1(V) -> H_1(S)].
```

This is the precise role of the complementary 4/8 adjacency convention.

## 2. Complementary-subsurface duality lemma

Let `U,V` be complementary compact subsurfaces of a closed oriented surface
`S`, meeting on their common boundary.  Over `Q`, put

```text
A = im[H_1(U) -> H_1(S)],
C = im[H_1(V) -> H_1(S)].
```

Then `C=A^perp` for the nondegenerate intersection pairing on `H_1(S)`.

Proof.  Poincare duality identifies the annihilator of `A` with

```text
ker[H^1(S) -> H^1(U)].
```

The long exact sequence of `(S,U)` identifies this kernel with the image of
`H^1(S,U)`.  Excision replaces the latter by `H^1(V,boundary V)`, and
Poincare-Lefschetz duality identifies that relative cohomology group with
`H_1(V)`.  Naturality of the duality and excision maps identifies its map into
`H^1(S)` with Poincare duality applied to inclusion `H_1(V)->H_1(S)`.  Hence
the image is exactly `C`, proving `C=A^perp`.

For `S=T^2`, the intersection form has dimension two.  Therefore

```text
rank A + rank C = 2.
```

## 3. Matching-rank theorem

Applying the lattice bridge to the subsurface lemma gives, configuration by
configuration,

```text
r_black + r_white = 2.
```

Thus the only possible rank pairs are

```text
(0,2), (1,1), (2,0),
```

exactly as observed in the finite oracle.  With

```text
q = 1[r_black>0] - 1[r_white>0],
```

we obtain

```text
q = r_black-1 = 1-r_white,
2q = r_black-r_white.
```

Hence the finite matching function is the homological balance observable

```text
M_N(p) = E[q]
       = 1/2 (E[r_black]-E[r_white]).
```

## 4. Scope and relation to the literature

The proof uses the ambient-homology image that Duncan--Kahle--Schweinhart
(`arXiv:2011.11903`) call a giant-cycle observable.  The 4/8 local bridge is
consistent with the modern digital-connectivity treatment of Cote and
Uzcategui-Aylwin (`arXiv:2503.17861`).  The resulting `q` is the wrapping
quantity in the finite periodic matching relation of Mertens and Ziff
(`arXiv:1603.07289`).

Those works provide the surrounding definitions; the repository-specific
proof obligation is discharged by the regular-neighbourhood lemma and the
16-pattern local certificate above.

The theorem scope assumes each periodic unit square is an embedded cell with
four distinct corners.  Tiny quotients with local identifications are not
silently covered by this cellwise proof, although the separate exact oracle
has already verified the declared small controls directly.

This homological identification does not by itself select a CFT field, prove
the absence of the lower `V_(2,2)` spin-4 contribution, or determine the
square-site critical threshold.
