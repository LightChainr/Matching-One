# Rank-gap thermal-window score

The paired threshold-rank gap is resolved extremely precisely, but the six
available sizes do **not** yet lie on a constant-amplitude `N^(5/8)` law.

The orientation-pooled means and fixed-exponent scaled views are:

| N | E[G] | delete-one SE | N^(-5/8) E[G] | delete-one SE |
|---:|---:|---:|---:|---:|
| 65 | 5.548175210 | 0.000302968 | 0.408393969 | 0.000022301 |
| 85 | 6.607616795 | 0.000383379 | 0.411299347 | 0.000023864 |
| 130 | 8.699175715 | 0.000501935 | 0.415205892 | 0.000023957 |
| 170 | 10.338020465 | 0.000532105 | 0.417260056 | 0.000021477 |
| 185 | 10.911188908 | 0.000297775 | 0.417724305 | 0.000011400 |
| 265 | 13.725339051 | 0.000374766 | 0.419753581 | 0.000011461 |

The fixed-exponent, one-amplitude GLS score is

```text
common scaled amplitude = 0.4167447324743813
chi-square / df          = 271122.20 / 5
```

This rejects a correction-free amplitude plateau, not the thermal exponent as
an asymptotic mechanism.  The scaled gap rises monotonically by about 2.8%
from N=65 to N=265.  The two common-random-number doubling diagnostics have the
same correction direction:

```text
(G_130/G_65) / 2^(5/8) = 1.0166797851 +/- 0.0000841746
(G_170/G_85) / 2^(5/8) = 1.0144923860 +/- 0.0000783133
```

Thus the new joint observable supplies a sharp finite-size correction target:
the leading thermal-window bridge is plausible in direction but is not a
precision asymptotic law on the current sizes.

The joint-shape views also drift coherently.  The gap coefficient of variation
increases from `0.76360` to `0.81195`, while `Corr(K_minus,K_plus)` increases
from `0.43529` to `0.49202`.  These statistics are secondary: they show that
the paired topological window has not reached a size-independent standardized
shape, and they are not counted as independent confirmations of the exponent.

All input moment identities passed exactly.  N=65/85/130/170 use aligned
delete-one covariance; N=185 and N=265 are independent because their counter
intervals are disjoint.
