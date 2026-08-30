# P334 result: one-pass bottleneck shape does not commonly absorb birth age

## Result

The fresh 20,000-path/size pilot added only current-state quantities obtainable
in one carrier traversal: boundary cut/contact counts, 2-core size, articulation
vertices, and bridges.  No full configuration was saved.  Two independent deck
cycle lengths were deliberately excluded because the observed state has ambient
rank one; defining both directions would condition on the future rank-two birth.

The frozen production-anchored decision is negative for common absorption:

| size | first retention | second retention | joint p |
|---|---:|---:|---:|
| N325 | `75.6%` | `67.0%` | `2.62e-5` |
| N425 | `30.4%` | `27.2%` | `0.435` |

The rule required both size-level joint p-values to be at least `0.01` and all
four retained magnitudes to be at most `25%`.  The new vector substantially
attenuates the N425 point estimates, but it leaves a large, jointly nonzero
N325 production-anchored association.  It therefore fails as a common
all-size/all-orientation bottleneck state.

The frozen N325 geometry-only map transferred to N425 leaves residual age
slopes `-0.0583` and `-0.0772`.  Their joint test is underpowered (`p=0.452`),
but the effect sizes remain comparable to the production age signal rather
than closing toward zero.

## Identifiability and calibration

The twelve nominal `Mshape` columns have rank eleven in every row.  The scorer
therefore uses a Moore--Penrose least-squares solution on the complete frozen
column span, without selecting or dropping a proxy.  Age itself adds one
independent rank, so its coefficient remains identifiable.  Individual
coefficients inside the redundant proxy span are coordinate descriptions, not
separately identified mechanisms.

H2 still gives the expected stronger attenuation, with fresh-pilot retained
magnitudes from `5.2%` to `37.0%`.  This is only the calibration ceiling:
`H2/(N-k0)` is the conditional one-step completion hazard by definition.

## Scientific card

- Mechanism space changed: current boundary cut, contact multiplicity, 2-core,
  articulation, and bridge summaries are removed as a common explanation of
  the production birth-age association.
- Not proved: intrinsic temporal memory, complete-state non-Markovianity, or
  the absence of a more targeted current bottleneck coordinate.
- Observer/source/geometry: rank-one state at N325 `k0=193` and N425 `k0=252`,
  two paired norm-five quotient orientations, fresh disjoint counter domains.
- Dependency: orientations within a size share random paths; the independent
  2M archive supplies the precision anchor and is not counted as another
  geometry vote.
- Next lift: a useful next coordinate must be more targeted than global
  carrier boundary/core summaries; repeating larger samples of this vector has
  low expected information value.

## Provenance

Raw files were locked at `d1be2e0` before any proxy distribution or score was
read.  Two scorer invocations stopped before producing output: first on an
exact redundant proxy column, then on an age-identifiability guard incorrectly
applied to the geometry-only transfer.  Commits `2caedb3` and `1540b28` made
only those fail-closed interface corrections; data, columns, hypotheses, and
thresholds did not change.  The completed frozen score was then produced once.
