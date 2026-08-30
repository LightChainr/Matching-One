# P250 augmented annihilation-extension map freeze

## Decision question

The radius-four joint score at `99d23a7` and radius-five line-extension score
at `11130ae` reuse the same 80k archive.  Their R2/R3 survivor conflict is
therefore not two independent votes.  This zero-sample analysis replaces those
separate gates by one candidate-constrained operator:

```text
conjugate(H_plus_old+degree-five-extension) q = 0,
H_minus_candidate(old+degree-five-extension) q = 0.
```

For each candidate this is a single complex `40 x 6` matrix.  The null is rank
at most five.  A maximum-volume `5 x 5` pivot is chosen from the full mean and
held fixed; the complete Schur residual has 35 complex (70 real) coordinates.
The D4 map acts on both the coefficient basis and the degree-three left shifts.

## Frozen candidates and de-duplication

The candidates are exactly the five maps frozen before the radius-five data:

- identity plus conjugation;
- Alexander reflection plus conjugation at `R^0`, `R^1`, `R^2`, `R^3`.

This is the intersection of the radius-five candidate list with the older
radius-four declaration.  No linear or new post-reveal map is introduced.
Every candidate receives exactly one augmented p-value.  In particular, the
old R2/R3 annihilation p-values and the radius-five R2/R3 direction p-values
are provenance, not additional decision gates.  We neither multiply them nor
select the candidate with the largest p-value.

## Two-source influence covariance

The old and fresh streams have independent seeds but different roles in the
same nonlinear Schur residual.  Their batch indices must not be paired.

1. Delete one old 80k batch while holding the fresh mean fixed.  Recompute the
   old Hankel block and all degree-three/four entries in the extension rows.
2. Delete one fresh 1.2M batch while holding the old mean fixed.  Recompute only
   the degree-five extension entries.
3. Compute the two jackknife influence covariance matrices separately and add
   them.  The same rule applies to the saved cross-candidate covariance.

The primary statistic is the correlation-normalized pseudoinverse quadratic
form at relative eigenvalue cutoff `1e-10`.  The finite-batch probability uses
the existing conservative project convention with `min(400,400)=400` batches.
The decision threshold is `alpha=0.01`.

## Interpretation

The result will report all five frozen candidates and the Alexander union.  It
tests whether one fixed truncated map accounts for annihilation and extension
simultaneously through degree five.  It is not an exact rank-five or flatness
certificate, does not identify physical state dimension, and does not observe
ordered `TxTy` versus `TyTx`.  Radius-six and rank-eight results are not reopened.

No new simulation is authorized.  Joint residuals must not be computed before
this freeze is committed.
