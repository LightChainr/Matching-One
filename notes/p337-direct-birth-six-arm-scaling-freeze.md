# P337 direct-birth four-generation scaling freeze

## Question and chronology

The conditional collision law

```text
D_N = A N^(-5/6)
```

was frozen before the P334 N325/N425 collision counts at `652bf15`.  This
analysis transports that same high-risk hypothesis to the exact P337
same-lineage doubling sequence `N=85,170,340,680`.  It does not choose an
exponent from the first three approximate rates and adds no simulation.

The four projective-birth archives already exist.  N680 remains a heldout point
in the present analysis chronology: its direct-birth count has not been used
to select the models, transformation, covariance, or threshold below.

## Event and paired-orientation covariance

`DIRECT_RANK2` is the typed `0 -> 2` birth with no intermediate primitive line,
equivalently the diagonal event `K1=K2`.  For each size, orientation and batch,
compute its count divided by the recorded number of paths.  The lineage scalar
is the equal-orientation mean

```text
Dbar_N = (D_first + D_second)/2.
```

The two orientations share a counter-keyed permutation.  Delete one common
batch index from both and preserve the resulting `2 x 2` covariance; do not
treat them as separate evidence.  The four sizes use independent seeds and
counter domains, so their covariance blocks are combined block-diagonally.

## Frozen model ladder

All fits use `log(Dbar)` with its delete-one delta covariance.

1. Primary conditional line: `log Dbar = log A - (5/6) log N`, one fitted
   amplitude and three goodness-of-fit degrees.
2. Effective power diagnostic: `log Dbar = log A - beta log N`, two fitted
   parameters and two residual degrees.
3. Minimal curvature diagnostic:
   `log Dbar = a - beta log N + kappa*(log2(N/85)-1.5)^2`, with one residual
   degree.  `kappa` is the only new curvature coordinate and cannot replace the
   fixed primary decision.

Two transparent secondary views are frozen.  First, jointly score all three
adjacent doubling contrasts against `2^(-5/6)`, retaining the covariance caused
by shared intermediate sizes.  Second, estimate the fixed-exponent amplitude
from N85/N170/N340 only and score N680 using both source-amplitude and heldout
variance.

The decision threshold is `alpha=0.01`.  Every model is reported; none is
selected by largest p-value.

## External lineage and boundary

The N325/N425 production values from `2e99533` are shown only as an external
lineage comparison.  Their shapes are not the P337 doubling genealogy and they
do not enter the four-size GLS fit or evidence count.

Even exact compatibility with `5/6` would establish only a **conditional
six-arm scaling line for this direct-birth observable**.  It would not prove an
arm-event correspondence, a universal amplitude, or universality across
geometries.  Failure could reflect corrections or lineage amplitude drift;
it is not evidence for or against thermal Q4.
