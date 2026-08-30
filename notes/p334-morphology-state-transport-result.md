# P334: morphology-state transport is partial and boundary-rotating

## Question and coordinates

This zero-new-sample analysis asks whether the P334 current-state morphology
contains a low-dimensional predictor that survives a change of size or
orientation.  It does not refit eleven free shape coefficients at every target.

The raw counts are first placed in physical dimensionless coordinates:

- mass and core counts are divided by area `N`;
- frontiers, boundary cut/contact counts, articulations, and bridges are divided
  by linear size `L=sqrt(N)`;
- the number of extra essential carriers is already dimensionless.

Within each line stratum, the source block supplies a canonical quotient of the
one exact redundant morphology direction.  The source environment coefficients
then define rank-one and rank-two predictive subspaces.  Target outcomes may fit
amplitudes only inside that fixed subspace.  Thus “held out” refers to the
subspace, not to a fully specified point forecast.

## Held-out transport

The rank-one result, precision-anchored to the disjoint 2M production age
slopes, is:

| source -> target | source coefficient energy | target retentions | residual joint p |
|---|---:|---:|---:|
| N325 -> N425 | `81.2%` | `45.1%, 46.5%` | `0.900` |
| N425 -> N325 | `95.3%` | `43.2%, 53.0%` | `0.0166` |
| first -> second | `83.5%` | `63.9%, 48.6%` | `0.0986` |
| second -> first | `94.4%` | `59.1%, 57.6%` | `0.692` |

The transported direction is real enough to remove roughly half the point
signal, but it is not a sufficient common state: every route retains `43--64%`
of the production association.  The large residual p-values in some routes
reflect source-subspace uncertainty carried through the full delete-one
covariance; they are not evidence that the retained effect is exactly zero.

Rank two does not repair the mismatch.  Its retained magnitudes span
`43.7--71.8%`, sometimes worse than rank one.  Both ranks were reported as
fixed descriptive levels; no model was selected by the largest target p-value.

## Where the state rotates

Across all four environments, the first coefficient direction carries `86.6%`
of squared coefficient energy (delete-one SE `10.8%`), and the first two carry
`93.8%` (SE `5.5%`).  That apparent low rank coexists with substantial
direction change:

```text
same orientation, N325 vs N425: cos(first)=+0.445, cos(second)=-0.866
same size, first vs second:      cos(N325)=-0.570, cos(N425)=+0.903
```

The stable localization is by physical sector rather than by individual angle.
`85.7% +/- 3.5%` of the cross-scale coefficient-rotation energy lies in the
frontier/boundary/contact coordinates.  Only `2.3%` lies in carrier-count,
articulation, and bridge coordinates; mass/core contributes `12.0%`.
Orientation rotation is also boundary-dominated (`79.6%`), though less precise.

So the previous failure should not be answered by adding more articulation or
bridge detail.  The next useful current-state coordinate would need to resolve
the *organization* of the active boundary/frontier, not merely its total length
or number of bottlenecks.

## Scientific card

- Mechanism changed: a partial transferable rank-one morphology direction
  exists, but neither rank one nor rank two closes birth age across both size
  and orientation.
- Rotation localized: scale mismatch is overwhelmingly a boundary/frontier/
  contact phenomenon, not a global mass/core or articulation/bridge effect.
- Not proved: intrinsic memory, a universal transfer state, or absence of a
  targeted boundary-shape sufficient statistic.
- Dependency: this is post-reveal reuse of the same fresh 20k blocks plus the
  disjoint production precision anchor; it is not an independent vote.
