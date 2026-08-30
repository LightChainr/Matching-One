# P250 candidate-constrained joint annihilation freeze

## Question

The model-free score at `a770ac9` keeps separate plus/minus rank-five charts
compatible but rejects one raw shared chart.  The bridge score at `a46ed63`
compares the two minimum-singular directions after each hand has been fitted
separately.  It does not retain the magnitude of either annihilation residual.

This zero-sample score asks the stronger question for every previously declared
map:

```text
H_plus q_plus = 0,
H_minus(T B2) q_minus = 0,
q_minus = q_plus or conjugate(q_plus).
```

For a linear map, vertically stack `H_plus` and the transformed minus Hankel
matrix.  For a conjugating map, vertically stack `conjugate(H_plus)` and the
transformed minus matrix.  In both cases the candidate null is exactly that
the resulting complex `24 x 6` matrix has rank at most five.  The spatial map
acts on both the left shift and the coefficient monomial; a column-only map is
not the frozen observable.

## Projective chart and covariance

Choose the maximum-volume `5 x 5` pivot of the full candidate mean.  Hold its
rows and columns fixed through all 400 delete-one calculations.  Set the one
nonpivot coefficient to one, solve the pivot equations for the other five
complex coefficients, and score the complete `19`-complex-coordinate Schur
residual.  Its `38` real coordinates are the generic codimension
`48 - dim_R(CP5) = 38` of the candidate null.

Deleting batch `b` deletes its plus and minus rows, both charges and every
displacement together.  Each candidate keeps its full `38 x 38` covariance.
The six decision-relevant maps also preserve their `228 x 228` cross-candidate
covariance, but no combined vote is formed.

## Mandatory identity replay

`orientation_preserving_R0_linear` must reproduce the `a770ac9` shared-block
rank-five calculation before any other candidate is reported.  The matrix and
algorithm are identical, so the pivot, residual, covariance, resolved modes,
statistic and finite-batch probability must agree to numerical roundoff.  This
is an implementation regression, not scientific evidence.

`orientation_preserving_R0_conjugate` is a different scientific control.  It
must not be substituted for the linear replay.

## Radius-five chronology

The later branch result `11130ae` retains only
`Alexander_R2_conjugation` at the frozen radius-five line-direction gate.  It
still does not impose the two annihilation equations jointly.  R2 is therefore
the post-radius-five primary candidate here, while all sixteen maps declared
at `a46ed63` remain visible and the four Alexander-plus-conjugation conventions
retain their union decision.

The 1.2M radius-five stream uses an independent seed, but its published score
also reuses this 80k archive.  Its p-value and this score are not independent
votes.  This analysis tests only a degree-two radius-four truncated map; it
does not establish exact rank five, a closed transfer algebra, ordered words,
a finite-quotient graph isomorphism or a continuum operator identity.
