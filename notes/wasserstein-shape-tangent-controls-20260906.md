# The affine tangent of #582 is the ray statistic of #579, with two columns

**Date:** 2026-09-06
**Tickets:** #582 (the ask), #579 (the machinery)
**Claim level:** C0 — synthetic controls only, no lattice claim, no data
**Artifact:** `results/wasserstein-shape-tangent-controls/latest.json`
**Script:** `scripts/wasserstein_shape_tangent.py`

#582 asks for a basis-free finite-size velocity of a whole threshold law, and for
the part of it that survives after the best possible center/width adjustment. Its
no-new-production first pass wants exact affine-null controls and a
covariance-aware affine tangent projection. Both are delivered here, and the
second one turns out not to need building.

## The observation

In one dimension the quantile function is an isometric coordinate for `W2`, so
the finite-size displacement is the function `v(u) = Q_M(u) − Q_N(u)`. An
infinitesimal affine map `p → a p + b` sends `Q` to `aQ + b`, so the affine
tangent at `Q` is

```text
span{ 1, Q }  ⊂  L²(0,1),
```

a **two-dimensional subspace**. Asking for the covariance-weighted distance from
`v` to that subspace is

```text
D = min_{α,β} (v − α·1 − β·Q)ᵀ S⁺ (v − α·1 − β·Q),
```

which is `projective_inference.subspace_residual` with a two-column basis — the
same routine that scores a model ray with a one-column basis. #582's shape flow
and #579's model exclusion are the same statistic at different `dim(V)`.

That is not a coincidence worth admiring; it is a reason not to write a second
implementation, and it means the χ² calibration, the pseudo-inverse rank
reporting and the degrees-of-freedom bookkeeping already tested under #579 apply
unchanged. Nine quantile levels minus two tangent columns leaves 7 df.

## What the weighting buys, which #582's own formulation does not have

#582 proposes `‖v_shape‖_L2`. That norm treats every quantile coordinate alike.
The same issue's tail section then notes that a quantile's uncertainty carries a
factor `1/f(Q(u))`, so coordinates in a thin tail are simultaneously the noisiest
and, in an unweighted norm, potentially the loudest — which is why that
formulation has to freeze a central window `u_min` and warns against choosing it
by whichever value makes the rank smallest.

Weighting by `S⁺` handles this without a window. The control:

```text
displacement = affine, plus a ONE-SIGMA excursion in the noisiest quantile
               (that coordinate's standard deviation is 400x the others')

unweighted ‖v_shape‖_L2  =  310.2          "large shape flow"
weighted D               =  1.00 on 7 df   0.006 sigma
```

The unweighted norm reports shape. The weighted statistic reports a residual
consistent with its degrees of freedom, because it knows that coordinate was
worth 1σ. The window does not stop being worth declaring — rank pathologies are
real — but it stops being where the conditioning lives.

## The controls, all exact

| control | result |
|---|---|
| pure translation | `D = 0` exactly |
| pure scale | `D = 6.5e−35` |
| translation and scale | `D = 1.6e−33` |
| negative scale | `D = 9.3e−34` |
| one declared shape generator | `0.310` on the affine basis, `1.0e−32` once the generator joins it |
| two declared shape generators | `0.402` alone, **`0.326` with only the first**, `1.2e−32` with both |

The second-to-last row is the one that matters for #582's decision table. A
projector that annihilated everything would pass the affine-null rows and fail
here; a rank-1 method applied to a rank-2 flow would close the last row with one
generator and report every flow as rank 1. Neither happens.

## What this does not establish

- Nothing about any threshold law. Every displacement here is synthetic.
- Not that a weighted shape norm removes the need to declare a quantile window.
  It removes the need for the window to carry the conditioning, which is a
  smaller claim.
- Not that the shape subspace of any real lineage is low-dimensional. That is
  the measurement #582 asks for, and it needs the delete-one covariance of a
  reconstructed quantile grid, which this pass does not touch.

## What comes next in #582's own order

The remaining first-pass items are data-side: a typed threshold-CDF → quantile
reconstruction utility, the exact tiny percolation control, and the aligned
delete-one covariance of the quantile grid for one lineage. The statistic is
ready for them; the covariance is the work.

One warning carried over from #579 and worth repeating here. The covariance must
be of one random object across all quantile coordinates, so the lineage's
transitions have to share a stream and the delete-one has to delete a batch from
every level at once. A quantile grid whose levels were estimated independently
has a diagonal `S` by construction, and a diagonal `S` on a hundred correlated
coordinates of one distribution will make almost any displacement look
significant.
