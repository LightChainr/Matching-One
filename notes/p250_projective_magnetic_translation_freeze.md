# P250 projective/magnetic-translation freeze

## Question

The fresh 80k bivariate stream rejected a common commuting state of rank at
most three while preserving the exact signed C4 channel covariance.  This
existing-data reanalysis asks the next narrower question: can the mixed rows be
the matrix elements of the minimal five-dimensional projective translation
algebra over Z5?

This is a post-result hypothesis designed after commit `4c1f5d8`; it is not a
preregistered interpretation of the 80k acquisition.  The scorer and all model
choices are frozen before the new residual scores are computed.

## Exact boundary first

The actual plus and minus fivefold covers have zero unit-plaquette fiber
curvature at every one of the 101 parent vertices.  Therefore the microscopic
translation bundle has `m=0`.  A fitted nonzero `m` below would be an effective
projective representation of the observed state/observable projection, not a
discovery of literal deck magnetic flux.

For each `m=1,2,3,4`, the canonical matrices obey

```text
Tx = X, Ty = Z_m,
Tx Ty = omega^(-m) Ty Tx,
R Tx R^-1 = Ty,
R Ty R^-1 = Tx^-1,
R^4 = I,
D = Tx Ty Tx^-1 Ty^-1 = omega^(-m) I,
D^5 = I.
```

The plus/minus hands receive conjugate center charges `(m,-m)`.  `R`, the
clock/shift matrices, and the center are exact and are never fitted.

## Frozen score

For a hand and displacement `(a,b)`, fit

```text
G(a,b) = rho^(|a|+|b|) v* W_m(a,b) v,
```

with one complex five-vector per hand, one shared real `rho`, Hermitian
source/sink, and the exact observed C4 map fixing the charge-2 rows.  Only the
seven axis points through degree three train the model.  The mixed points
`(1,1),(2,1),(1,2)` and all five degree-four points are scored separately with
the full 400-batch delete-one covariance.  A Weyl model survives only if both
scores have `p>=0.01`.

The control is deliberately favorable to the commuting explanation.  Free
commuting ranks four and five share their eigenpairs across channels but have
channel-specific complex amplitudes.  Because four axis moments cannot
identify those ranks, they may train on every first-quadrant point of total
degree at most three, including the mixed points, and face only degree four as
held out.  Thus a commuting survivor diagnoses insufficient dimension in the
earlier rank-1/2/3 gate; a Weyl survivor accompanied by commuting-rank-5
failure supports noncommuting state algebra.

Failure of both is still informative but narrow: it rejects the canonical
one-vector Weyl realization and a free diagonal rank-five realization, not all
five-dimensional projective, Jordan, or periodic-image models.

The numerical optimizer is a local analysis dependency, recorded as NumPy
2.4.6 and SciPy 1.18.0.  The immutable input batch SHA-256 is
`dfbd83d680080939bcae5cad0090af95d82748ed780ce73807d4b002d6a9ca72`.
