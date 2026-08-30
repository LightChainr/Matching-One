# Exact aligned-batch common-ray wedge control

This control addresses the parameter-free diagnostic in Issue 439 without
using production data. Four generations carry aligned batch coordinates
`(A_M, A_K)`. For each adjacent pair it computes

```text
D_N = A_M(2N) A_K(N) - A_K(2N) A_M(N).
```

Every delete-one-batch replicate first recomputes all four coordinate means and
then recomputes the three products. The multivariate jackknife covariance is
formed from those nonlinear replicates; products are never treated as
independent measurements.

Two exact four-batch controls are frozen:

- a common-ray family with loading `A_K=2 A_M` at every generation, for which
  all full-sample and delete-one wedges and the covariance vanish exactly;
- a loading-drift family with generation loadings `(2,5/2,3,7/2)`, whose
  full wedges are `(-25/16,-25/64,-25/256)` and whose exact jackknife
  covariance has rank one.

The oracle checks every principal covariance minor exactly and fails closed on
floats, incomplete generation sets, duplicate batches, invalid fields, or a
non-doubling size lineage.

## Boundary

This is a statistical-algebra certificate only. It does not import the
N=85/170/340/680 archives, estimate empirical covariance, fit or reject a
common ray, compare transfer laws, decompose direct rank-2 events, forecast
N=1360, or support a physics conclusion. Issue 439 remains open.
