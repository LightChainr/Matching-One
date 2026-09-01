# P250 augmented old-plus-degree-five joint-operator freeze

## Decision target

The radius-four joint-annihilation gate and the radius-five line-direction gate
left disjoint marginal survivor sets.  This zero-sample analysis does not cast
another vote between R2 and R3.  It asks whether one fixed projective
coefficient vector satisfies the old annihilation equations and the already
acquired degree-five extension equations simultaneously.

The five primary candidates are exactly the fixed-map family present before the
fresh radius-five reveal: identity plus conjugation and Alexander reflection
composed with `R0..R3` plus conjugation.  However, this augmented joint gate was
designed only after the marginal R2/R3 conflict was seen.  It is therefore a
**primary retrospective zero-sample synthesis**, not a prospective
confirmation.  The other eleven maps already frozen by the later
joint-annihilation manifest are scoreable from the same raw coordinates, but
remain secondary retrospective existing-data views.

## One operator, not two p-values

For each map, stack:

1. the existing `24 x 6` complex candidate-mapped radius-four matrix; and
2. sixteen complex degree-three-shift rows, eight per hand.

Each extension row has one degree-three column and two degree-four columns from
the old 80k stream, followed by three degree-five columns from the independent
1.2M shell.  The resulting `40 x 6` matrix is tested for rank at most five.

Reuse the candidate-specific old `5 x 5` pivot from `99d23a7`.  The first
nineteen complex Schur coordinates must therefore replay the old result
exactly; the sixteen new coordinates are additional equations for the same
projective `q`.  The pivot is never reselected after seeing the fresh stream or
inside a delete-one recomputation.

## Independence and covariance

Within the old stream, one delete-one removes the same aligned batch from both
hands, both charges and all radius-four displacements.  Within the fresh
stream, one delete-one removes the same aligned batch from both hands, both
charges and all twenty shell points.

The streams have distinct seeds.  Their batch labels `0..399` are not paired:

- recompute 400 old-delete residual vectors with the fresh full mean fixed;
- recompute 400 fresh-delete residual vectors with the old full mean fixed;
- form complete old and fresh jackknife covariance matrices separately;
- add them once, retaining their separate matrices in the result.

The old covariance includes the cross-covariance between the replayed old rows
and the mixed extension rows, because the latter contain old degree-three and
degree-four moments and inherit uncertainty in `q`.  The fresh contribution to
the old 38-real-coordinate block must be zero up to numerical roundoff.

The project-compatible finite-batch Hotelling calibration uses
`B=min(400,400)` and is saved alongside the asymptotic chi-square value.  It is
an approximation for a sum of two independent jackknife covariance estimates,
not an exact two-sample pivot law.

Both response JSONs, both exact gates, and the pinned runner sources are part
of the input contract.  The scorer checks common `p=0.59274605079`, N505
geometry, hands `plus/minus`, charges `1/2`, the affine-fiber C4 gauge, and the
shared uniform translated-origin convention.  The translation salts are
intentionally distinct.  Because the fresh response does not repeat
`child_order`, its N505 provenance is checked through its declared runner,
which imports `CHILD_ORDER=505` and both hand contexts from the pinned shared
geometry runner.  Radius-five names are matched by an explicit alias table;
names are never inferred by string similarity.

The human-readable score is a small summary.  Residual vectors and complete
old/fresh covariance components live in a compressed NPZ payload; total
covariance is derived exactly as `old + fresh`.  One canonical plus-hand block
and eight minus-hand geometry blocks are stored once and referenced by all
candidate summaries.

## Interpretation gate

Also score each `20 x 6` hand-specific augmented matrix with its frozen old
pivot.  A rejection of the joint map is a bridge-only result only if both
hand-specific rank-at-most-five gates survive.  If either hand fails, the live
conclusion is that the full degree-five five-state extension is incomplete;
one must not promote the fallback to a general `5+5` quotient.

The observed result does not support the current general `5+5` fallback,
because the plus-hand full-Schur gate also fails.  It also does not name rank
six.  Survival would mean compatibility with a declared truncated chart, not
exact rank five, a completed `Tx,Ty` algebra, a microscopic quotient
isomorphism or a continuum field identity.
