# P250 result: no frozen fixed map jointly annihilates and extends

## Result

The frozen augmented operator combines the old radius-four joint annihilation
rows and the already acquired radius-five degree-three shifts in one `40 x 6`
complex matrix per candidate.  At `alpha=0.01`, all five parameter-free maps
are rejected:

| frozen candidate | chi-square / df | finite-batch p | decision |
|---|---:|---:|---|
| identity + conjugation | `290.595 / 70` | `3.67e-14` | rejected |
| Alexander R0 + conjugation | `266.284 / 70` | `2.46e-12` | rejected |
| Alexander R1 + conjugation | `189.993 / 70` | `1.03e-6` | rejected |
| Alexander R2 + conjugation | `210.366 / 70` | `3.55e-8` | rejected |
| Alexander R3 + conjugation | `216.611 / 70` | `1.24e-8` | rejected |

No candidate is selected by relative p-value.  The Alexander union is rejected
because all four preregistered members fail the same augmented gate.

## What the joint score resolves

The old radius-four gate retained R3 while rejecting R2; the later direction
extension retained R2 while placing R3 just below threshold.  Those were not
independent votes because both reused the 80k archive.  The augmented score
removes the ambiguity: **neither R2 nor R3 supplies one fixed rank-five map
that satisfies annihilation and extension simultaneously.**  The apparent
survivor swap was a consequence of asking two projections of one constrained
problem, not evidence for two viable physical intertwiners.

The result also sharpens the module interpretation.  Separate hand-specific
degree-five extensions can be compatible while every declared cross-hand
identification fails.  Thus the missing structure is in the inter-sector map,
not merely another vote on a scalar line.  This agrees with the later need for
higher-rank/common-state or morphism-enriched descriptions, without using
those later data as part of this decision.

## Covariance and de-duplication

The 80k and 1.2M streams were kept as separate statistical sources.  For every
candidate, the old delete-one replicate changes both the radius-four block and
the degree-three/four entries in extension rows; the fresh delete-one replicate
changes only degree-five entries.  Their centered influence arrays have shape
`5 x 400 x 70` for each source.  Each source covariance is formed separately
and the two are added.  Batch numbers from independent streams are never
paired or pooled.  The maximum numerical covariance-addition discrepancy is
`6.62e-24`.

The published old/fresh R2 and R3 p-values are retained only as chronology.
They are not multiplied, combined, or counted as separate evidence.  The
machine-readable influence artifact reconstructs every candidate covariance
and the saved full cross-candidate covariance.

## Boundary and next object

This is a Level-S elimination of five fixed, parameter-free, rank-five maps on
moments through total degree five.  It is not an exact rank or flat-extension
certificate, a physical state-count proof, an ordered `TxTy/TyTx` observation,
or a microscopic graph isomorphism.  General hand-specific maps,
modulus-dependent maps, higher-rank common states, noncommutative translations,
and context-dependent morphisms remain outside this gate.

The informative next object is therefore not another R2/R3 line vote.  It is a
genuinely morphism-enriched/common-state observable with ordered or contextual
rows, consistent with the later #249/#255 direction.

## Scientific card

- Mechanism space changed: all five fixed parameter-free degree-five rank-five
  cross-hand maps are eliminated; the R2/R3 survivor conflict is closed.
- Not proved: exact rank, physical state dimension, noncommutation, or path
  memory.
- Observer/sector/source/geometry: norm-505 projective-leg Z5 charged two-point
  moments, plus/minus hands, charges 1/2, radius-four archive plus radius-five
  shell.
- Dependency group: one 80k old block and one independent 1.2M fresh block;
  five candidate views share both blocks.
- Next lift: ordered/common-state morphism rows, not a renewed scalar-map vote.
