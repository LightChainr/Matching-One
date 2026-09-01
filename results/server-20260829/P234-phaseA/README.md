# P234 Phase A: fixed-delta parent-pair result

Status: completed production result for `delta=1/(8 sqrt(2))` at
`L=64,96,128,192`, with 100,000 configurations and 100 batches per size.

The raw archive contains the exact cluster-sign-integrated `(LL,LD,DD)`
sufficient statistics and their same-stream covariance.  The batch columns
`sum_D1` and `sum_D2` also retain the bilocal endpoint connection probability.
At one fixed physical cutoff this probability is proportional to `pi_a^2`, so
the relative Camia--Feng normalization can be recovered without a second run.

## Chronology

Commit `9252137` froze the transform and primary score while `L=192` was still
running.  The three-size preview had `chi2=5.339/4`, `p=0.254`.  The subsequently
revealed `L=192` block is therefore a held-out fourth point for that score.

## Held-out L=192 block

```text
(LL,LD,DD)
  = (-2.65159e-7, 2.410107e-5, 3.634710e-4)
SE
  = (4.67728e-6, 2.79549e-6, 7.07998e-6)
J = -1.1659 +/- 2.9328
```

The connection-normalized four-size vectors in `L` order are

```text
64   (+0.03167, 0.23189, 1.81512)
96   (-0.04088, 0.21866, 1.59275)
128  (+0.02498, 0.21772, 1.61755)
192  (-0.00490, 0.22740, 1.75092)
```

## Frozen primary score

The primary model was

```text
Y_LL = c_LL/L,
Y_LD = C_LD + c_LD/L,
Y_DD = C_DD + c_DD/L.
```

The held-out result rejects this single monotone finite-mesh correction:

```text
chi2 = 22.95795 / 7 df
p    = 0.0017336
```

The rejection is localized rather than generic.  Separate coordinate fits give

```text
LL zero-limit:       chi2=2.358/3, p=0.501
LD finite limit:     chi2=0.201/2, p=0.904
DD finite limit:     chi2=21.632/2, p=2.01e-5
```

Thus the bottom-null and nonzero mixed-pairing pieces survive, while the
top--top coordinate detects the signed nearest-vertex cutoff displacement.

## Post-reveal mechanism diagnostic

After the frozen score was revealed, the natural realized-cutoff prefactor was
used and the deterministic signed error
`delta_realized-delta_declared` was added beside `1/L`.  This is explicitly a
post-reveal mechanism analysis, not a replacement preregistered score:

```text
chi2 = 2.36234 / 4 df
p    = 0.66944
```

The inferred same-gauge cutoff shear is

```text
kappa_proxy = -delta (d DD/d delta)/(2 LD)
            = 2.788 +/- 1.341.
```

This is evidence for cutoff/basis sensitivity concentrated in the top partner,
not yet a universal logarithmic coupling measurement.  The planned
`delta=1/(12 sqrt(2))` and `1/(16 sqrt(2))` lines provide the independent test
of the derivative and remove reliance on four rounding offsets.

Primary machine: Huawei DevEnvC_HZsCM6, aarch64, 16 workers.
