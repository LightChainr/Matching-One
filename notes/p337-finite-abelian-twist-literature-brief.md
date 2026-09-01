# P337 finite-abelian twist tomography: literature boundary and next bridge

## Repository result

The branch-only result
`theory/p337-finite-abelian-twist-tomography-20260830@4b17955946207f39762f836156aa063f64fbf67d`
is an exact finite-volume cohomology transform, conditional on the unmerged
integral-saturation parent `c1a72e5`.

For a finite abelian group `A` of order `n`, let `alpha` range over
`Hom(Z^2,A)` and let `T_alpha` be the probability that `alpha` annihilates the
occupied ambient homology image. Then

```text
S_n = sum_alpha T_alpha = n^2 P0 + n P1 + P2.
```

Together with `S_1=1`, the two values `S_2,S_3` invert the complete unmarked
rank source `(P0,P1,P2)`. For `A=F_q`, each nonzero twist has a projective-line
kernel and returns one primitive winding-line bin modulo `q` after subtracting
`P0`. Aggregate twists and individual twists therefore answer different
questions: source tomography versus modular projective-line tomography.

The N65 archive reanalysis at
`analysis/p337-p334-flat-twist-n65-20260830@a7cb19a1b0f1115a7e5dfd23421a5dd25892fe78`
checks the transform without new samples. F3 is the smallest basis in this
geometry that separates the axes pair, the diagonals pair and their two
reflection-odd directions. Its balanced H4 contrast is statistically tied
with raw `chi4`; the sharper exploratory `F3_diagonal_odd` row is a different
projective reflection-odd sector and reuses the same 20k block.

## What the literature supports

Twisted torus partition functions are a well-established route to resolving
state spaces and symmetry sectors. Jacobsen, Ribault and Saleur use insertion
of a global-symmetry group element along a torus cycle to decompose Potts and
O(n) state spaces into irreducible characters:

- https://arxiv.org/abs/2208.14298

Critical Potts topological defects and their torus/cylinder partition
functions can also be realized as deformable lattice MPO defects, with fusion
and conformal-character extraction:

- https://arxiv.org/abs/2107.11177

Integrable seams provide another explicit lattice-to-conformal realization of
twisted boundary conditions in minimal models, including the three-state
Potts model:

- https://arxiv.org/abs/hep-th/0106182

These papers support the strategy of comparing a vector of twist-sector
partition functions with continuum characters. They do **not** identify the
P337 observable automatically. P337's `alpha` is a cohomology functional that
annihilates the occupied ambient homology image; it is not yet shown to equal
an `S_Q` group-element insertion, an MPO topological defect or an integrable
seam.

## Missing bridge

A continuum/module claim needs three additional pieces:

1. a finite-lattice operator or transfer-matrix identity relating the
   cohomology constraint to a declared seam/defect insertion;
2. the exact basis convention and induced modular/D4 action on the complete
   twist vector, not only one H4-like projection;
3. a fresh multiscale score with fixed normalization showing which twist
   characters carry stable amplitudes.

The most informative next acquisition is therefore a prereveal N130 child
with the complete F3 twist vector and the conditional birth/exit composition
in one aligned covariance block. Primary rows should separate the ordinary
balanced axis-minus-diagonal character from the reflection-odd diagonal
character. Agreement or disagreement across the child map then distinguishes
an ordinary H4 line-sorting continuation, a modular projective sector and a
purely microscopic N65 sorting effect.

## Boundary

- A finite list of primes resolves an integral winding line only modulo those
  primes unless an external size bound is supplied.
- Integral saturation makes finite-field rank equal integral rank; there is no
  additional Smith-defect tomography in `r_q` on the theorem-supported
  carriers.
- The N65 conditional-flux and F3 rows are post-hoc views of one exploratory
  20k block and must remain one dependency group.
- No local CFT field, affine-TL module, universal amplitude or exponent is
  identified by the current exact transform or smoke result.
