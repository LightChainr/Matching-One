# P253 N365 held-out annulus recurrence

Status: completed frozen third-geometry acquisition.

The old N325/N425 readouts determine one saturated dyadic recurrence

```text
g8 - T g4 + D g2 = 0
```

per plus/minus channel.  N365, using the primitive orientation pair
`(14,13)` versus `(19,2)`, was acquired from an independent replica-counter
block and enters only as a held-out third readout.

## Result

```text
channel  residual at R8     SE          z       p(1 df)
plus     +0.00168196         0.0084103  +0.200   0.841
minus    -0.00016686         0.0009644  -0.173   0.863

joint chi2 = 0.07103 / 2 df
joint p    = 0.96511
```

The off-grid `R=7` propagation check, using the old point class and only the
new N365 `R=2,4` amplitudes, is also compatible:

```text
plus positive-real branch: z=+0.348, p=0.728
minus principal-complex:    z=-0.166, p=0.869
```

This is a genuine third-readout confirmation that the old two-state radial
recurrence was not merely an algebraic saturation of two geometry outputs.
It does not distinguish a repeated root from the exact gap-one ordinary pair
or a complex pair.  The limiting uncertainty is inherited from the old
`T,D` covariance, especially in the plus channel; increasing N365 statistics
alone is therefore low information per CPU.

Run: 200,000 samples, 200 batches, 16 Huawei aarch64 threads, radii
`2,4,7,8`; elapsed production time was about two seconds after compilation.
