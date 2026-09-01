# Issue #275 forward identifiability on existing K1/K2 production

## Decision

The two mechanism classes are **PARTIALLY_IDENTIFIABLE, but not uniformly identifiable**.
A fixed or spectrally separated semisimple second mode has a different three-generation
image from a Jordan shear.  However, after all amplitudes are retained and no amplitude
bound or spectral gap is invented, the `kappa -> 1` closure of the semisimple image
contains the full Jordan image.  The present result therefore does not identify a Jordan
operator.

No Monte Carlo sample is generated.  The score consumes the complete selected submatrix
of the pinned cross-size aligned-delete-one covariance.

## Existing-data scores

| model | kappa | chi-square / df | survival p | reading |
|:--|--:|--:|--:|:--|
| fixed ordinary semisimple `q2` | 0.5 | 11.792466 / 4 | 0.0189631 | excluded at .05 |
| Jordan affine log | 1 | 6.432527 / 4 | 0.169092 | not excluded |
| free kappa, unconstrained | 2.114984508 | 2.156631 / 3 | 0.540541 | descriptive optimum |
| physical decaying semisimple | `0<kappa<1` | infimum 6.432527 at `kappa -> 1-` | diagnostic df3 p=0.0923618 | collision boundary, no interior winner |

The unconstrained optimum is `kappa=2.114984508`, equivalent to
`delta=-log2(kappa)=-1.080647096`.  Its negative
relative exponent means that the second mode grows relative to the proposed leading
`N^-13/8` response; it is not a more-irrelevant bulk singlet.

## Raw coordinates and units

For activation `i=1,2`,

```text
d_i(N) = [F_i(first,p_bar)-F_i(second,p_bar)] / Delta cos(4 theta),
u_i(N) = N^(13/8) d_i(N).
```

`d_i` is a dimensionless probability response per unit `Delta cos(4 theta)` at the
same pooled moving matching root as the source archive.  `A_top=d1+d2`,
`E_top=d2-d1`, and the linearized root shifts are deterministic transforms and are not
counted as extra evidence coordinates.

## Exact forward maps

For each lineage and activation, the fixed-`kappa` semisimple design and Jordan design are

```text
X_S(kappa) = [[1,1], [1,kappa], [1,kappa^2]],
X_J        = [[1,0], [1,1],     [1,2]].
```

Thus, from the first two generations,

```text
semisimple: d_hat_i(4N) = (4N)^(-13/8) *
             [(1+kappa)(2N)^(13/8)d_i(2N) - kappa N^(13/8)d_i(N)],
Jordan:     d_hat_i(4N) = (4N)^(-13/8) *
             [2(2N)^(13/8)d_i(2N) - N^(13/8)d_i(N)].
```

There are 4 blocks.  At fixed
`kappa != 1`, both complete designs have rank 8;
their combined rank is 12 and
their image intersection has dimension 4.
Per block the intersection is exactly the pure leading constant line `(1,1,1)`.

At a single geometry the change of basis from `(A_top,E_top,B_bulk)` to
`(F1,F2,B_bulk)` has determinant `1/2`.  Both the semisimple direct sum and an
`E_top/B_bulk` Jordan shear are rank three and have image `R^3`; adding an arbitrary
bulk column without a transport law is therefore exactly nonidentifying.

## Why separation is non-uniform

Write `kappa=exp(-epsilon)`, `b=-c/epsilon`, and
`a_semisimple=a+c/epsilon`.  For generation `g=0,1,2`,

```text
a_semisimple + b kappa^g  ->  a + c g.
```

Therefore the closure of the unbounded semisimple image contains the Jordan plane.
The fixed `q2` rejection removes only `kappa=1/2`; it does not eliminate every distinct
bulk singlet whose transfer eigenvalue is allowed to approach the topological one.

## Unique missing physical input

The missing input is **restricted_trace_modulus_or_phase_transport**: one
semantics-matched transfer relation on rank-0/rank-2 restricted traces for the
same `B` source, original `q/E` observables, and pooled-root physical
normalizer.
The preferred route is a theory-derived restricted-trace vector for the independent
singlet and Jordan top component, scored on existing rho-child or P43/P57 assets.  If
that relation cannot be derived, the only new acquisition justified by this audit is one
phase-calibrated second physical rotation of the same `B` column—not another untyped
topological coordinate.

## Provenance and boundary

- Input: `results/norm4-two-activation-h4/latest.json` (`f29ce76fa5be92abb2a233c7efb6e4d94f37236242656d0dcf17cbfc3fd1e462`).
- Manifest: `analysis/p275_forward_identifiability_manifest.yaml` (`80f027e73ce10da8f8afbce193bb81d448aa1726e62da5c57901c624ff0bc386`).
- All four residual rows are correlated views of the registered dependency groups.
- Scores are post-reveal existing-data diagnostics, not prospective validation.
- Jordan compatibility is not continuum-field, lattice-overlap, or normalizer identification.
