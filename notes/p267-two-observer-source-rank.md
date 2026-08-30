# P267 Target 1: two-observer x two-source projective rank lane

## Mechanism decision

The previous Target 1 stream established a large bulk coupling of `O_far` to
`JD_perp`, but one observer cannot distinguish a genuine two-source plane from
two lattice sources transported in the same projective direction.  The next
discriminant is therefore the complex coupling matrix

```
             JD_perp        JS
O_far        C_11           C_12
O_sep4       C_21           C_22
```

and not any one entry.  Its determinant is the primary observable.  A zero
determinant is the rank-one/common-projective-lane hypothesis; a resolved
nonzero determinant shows that the two observers see a two-dimensional source
plane.  It still does not name either direction as a continuum field.

## New external observer

At each pre-insertion site, one counter-derived C4 turn is shared across the
primal/matching traces and both torus orientations.  At physical distance
`R=6` the runner evaluates a local arm/landing mark on an axial anchor and a
diagonal anchor.  It retains the ordered pair as
`axis_landing + i*diagonal_landing`; the spatial H4 projection used in the
matrix is `O_sep4=axis_landing-diagonal_landing` (twice the normalized
projection).  The internal axis/diagonal arm types remain separate controls.

This is not a root contact observable: the exact quotient gate proves that the
source root lies outside every local ring in all four production geometries.
The two direction orbits are distinct and the scalar/H4 response matrix is
`[[1,1],[1,-1]]`, with determinant `-2`.  Thus the one-orbit alias of `83e98fc`
is removed before data collection.

## Exact and conditional layers

- **Exact:** quotient anchors/rings are injective and source-separated at R=6;
  axis and diagonal orbits are disjoint; the response map has rank two; reverse
  complement uses the identical C4 turn; the runner retains all same-batch
  products needed for the coupling matrix.
- **Frozen statistical lane:** `JD_perp=JD-beta*JS` uses per-orientation
  `beta=Re<JD,JS>/<|JS|^2>` and recomputes beta and the intrinsic root in every
  delete-one replicate.  The primary null is the complex determinant zero.
- **Mechanism inference only:** a resolved rank-two matrix would exclude the
  simplest single-projective-source explanation of the old D/S phase lock.  It
  would not by itself identify thermal Q4 epsilon, H4, or a scaling exponent.

## 20k smoke and production choice

Both N325 and N425 passed all complement/mapping audits.  Their 20k determinant
estimates were `-1.263-2.246i` and `1.660-3.051i`; the imaginary parts agree in
sign despite independent counter ranges.  The corresponding two-dimensional
chi-squares were only `1.88` and `2.64`, as expected for a smoke run.  Scaling
the standard errors to 2M gives ample power even if only roughly one fifth of
the smoke amplitude survives.  The frozen production is therefore 2M samples
per size, 100 batches, on independent Huawei hosts.

The manifest `experiments/p267_two_observer_source_rank_20260830.yaml` freezes
the matrices, radius, seeds, counters, matrix order, joint chi-square decision,
and claim boundary before production.
