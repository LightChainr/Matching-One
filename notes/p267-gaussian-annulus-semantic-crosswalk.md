# Target 3: a real Gaussian x annulus context rectangle

Status: post-reveal semantic audit and missing-cell freeze. No new production
has been authorized or started.

## The existing blocks do not form a rectangle

PR277 and P253 both discuss low-dimensional transfer and both use Gaussian
quotients, but their numerical rows are not the same observables.

PR277 starts from complete threshold-rank curves. Its source is a global
Bernoulli-parameter derivative. `U` is a global orientation/rank-projector
ratio; the `r2..r6` coordinates are width-normalized Hermite--Krawtchouk
derivatives. The transfer coordinate is Gaussian cover generation.

P253 holds `p` fixed, toggles one root, conditions on primal/matching pivotal
mass, and reads a local landing-shell H4 mark. `A_plus` is matching-even and
`A_minus` is matching-odd. The transfer coordinate is annulus radius.

These differ in source, observer, normalization, and transfer coordinate.
Dividing the PR277 row by a pivotal mass after reveal or multiplying P253 by a
power of `N` does not change the sigma-algebra. The exact crosswalk therefore
contains zero eligible cross-context pairs, and no old-data numerical 2x2 is
reported.

## Rows selected for the future rectangle

The same local row is already present in both parity sectors of P253, so it is
the cheaper semantics to transport:

```text
ordinary row      A_plus  = (H4_primal+H4_matching)/pivotal_mass,
matching-odd row  A_minus = (H4_primal-H4_matching)/pivotal_mass.
```

Using the full 16-dimensional P253 covariance, a post-reveal R8 prediction
from R2/R4/R7 was compared for the frozen `lambda=0,1/2,1` confluent bases.
The N425 pair maximizes the minimum adjacent-model Mahalanobis distance, but
the value is only `0.04332` (`N325: 0.0003704`). This selects N425 as the
annulus context; it does not identify a lambda.

PR277 does contain more sensitive internal coordinates (`N85_U` and
`N65_r5` are its best individual rows), but neither is synonymous with a local
landing-H4 row. They remain design diagnostics, not rectangle cells.

## Minimal missing acquisition

Fill only the Gaussian-cover row of the matrix by running the P253 fixed-p
root-toggle observable at `R=2` along the two PR277 cover lineages. Both
`A_plus` and `A_minus` are emitted from every sample and share one covariance
block.

The existing cyclic multiradius runner was compiled and preflighted. It
accepts all eight primitive parent designs at R2, but correctly refuses the
nonprimitive norm-4 children (`Smith=(2,N/2)`) because its vertex encoding
requires `gcd(a,b)=1`. Therefore production remains unauthorized. The minimal
implementation is an adapter: reuse P253 landing/pivotal/statistics verbatim
on the already-existing general integer-period geometry backend. No generic
new FK or topology framework is needed.

## Frozen held-out score

For each row, the Gaussian context uses the four-generation recurrence

```text
x3-2*x2+x1-lambda*(x2-2*x1+x0)=0.
```

The annulus context predicts R8 from R2/R4/R7 in the corresponding confluent
basis. The shared model fixes one `lambda` across both contexts and rows. The
minimal context-enriched adversary fixes one lambda per context, still shared
between `A_plus/A_minus`. All three diagonal and all nine context-pair scores
are reported; a 100,000-draw frozen Gaussian bootstrap calibrates the discrete
minimum improvement. A context gain demonstrates context dependence, not by
itself path/state memory.

## Scientific card

1. SEMANTIC RESULT: PR277 global thermal derivatives and P253 local root-toggle amplitudes are not the same rows.
2. NO PSEUDO-RECTANGLE: zero old cross-context pairs pass source, observer, normalization and batch gates.
3. ROW FREEZE: N425 `A_plus` ordinary and `A_minus` matching-odd are the common future readouts.
4. MISSING CELLS: fixed-p R2 local pivotal rows on the two norm-4 Gaussian cover lineages.
5. DECISION: shared versus context-specific lambda is scored only after a general-period adapter passes exact preflight.

