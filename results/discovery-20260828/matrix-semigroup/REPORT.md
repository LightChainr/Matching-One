# Matrix Gaussian-semigroup discovery score

**Classification:** post-reveal discovery.  N=185/265 are used in the fit, so
none of the model comparisons below is a preregistered claim.

## Main result

The four P48 channels admit a compact **three-dimensional** dynamic state:
one rank-1 identity/matching-even state read out by `P4[S]` and `P4[D']/Mbar'`,
plus one rank-2 thermal/matching-odd state read out by `P4[D]` and
`P4[S']/Mbar'`.

Using the intrinsic coordinate `u=Mbar(p)` does **not** remove the S-prime
drift.  The reconstructed metric amplitude

```text
Mbar' / N^(3/8): 1.751289 (N=65) -> 1.745139 (N=265)
```

changes by only 0.35%, whereas `N^(5/4) P4[S']` changes by roughly 43% on the
same clean sequence.  The pure four-channel score becomes slightly worse:

| coordinate | chi-square / dof | AICc |
|---|---:|---:|
| original p derivatives | 284.440 / 20 | 294.545 |
| intrinsic u derivatives | 289.248 / 20 | 299.353 |

Therefore the S-prime correction is not a simple thermal-metric artifact.

## Minimal matrix state score

The score uses the full aligned delete-one-batch covariance across the clean
N=65/85/130/170 Gaussian-doubling curves and independent 500M N=185/265
covariance blocks.

| thermal rank-2 action | chi-square / dof | AICc | BIC |
|---|---:|---:|---:|
| none (four pure constants) | 289.248 / 20 | 299.353 | 301.960 |
| ordinary `A+B/N` in S-prime only | 41.350 / 19 | 54.683 | 57.240 |
| Jordan `A+B log N` in S-prime only | **23.274 / 19** | **36.608** | **39.165** |
| ordinary correction in D and S-prime | 41.346 / 18 | 58.288 | 60.415 |
| Jordan correction in D and S-prime | 23.233 / 18 | 40.174 | 42.301 |

The smallest useful model is consequently a rank-2 thermal state whose second
component is visible in the `S'` readout but almost annihilated by the central
`D` readout.  Adding the companion to `D` improves chi-square by only 0.041 in
the Jordan model and costs one parameter.  This is the matrix mechanism that
unifies the previously successful `D` channel with the drifting `S'` channel.

The discovery fit is

```text
U_T2(N) = N^(13/8) P4[S']/Mbar'
         = -0.409331(102) + 0.368618(23) log N.
```

Equal-parameter comparison favors this Jordan action over the ordinary action
by `Delta AICc = 18.08` and `Delta BIC = 18.08`.  This is a discovery ranking,
not a replacement for the frozen Issue #72 chronology.

## Exact semigroup and event-involution blocks

The exact `cross/either` complementarity is the affine involution

```text
C = [[1,0],[1,-1]] on (1,S), C^2=1.
```

After centering at `S-1/2`, it is `diag(1,-1)`; differentiation gives the same
minus sign for `S'`.  The candidate norm actions are

```text
ordinary: R_q = diag(1,q^-1) on (1,N^-1)
Jordan:   J_q = [[1,0],[log(q),1]] on (1,log N)
```

and obey `T(q1)T(q2)=T(q1*q2)` exactly.  Both commute with the event involution;
the nilpotent Jordan direction stays inside one involution-parity block.

## New numerical predictions

For norm-2 children of the new geometries,

| child geometry | Jordan `U_T2` | ordinary `U_T2` | model gap |
|---|---:|---:|---:|
| N=370, `(17,9)/(19,3)` | **1.77049 +/- 0.03465** | 1.58270 +/- 0.02474 | 0.18779 |
| N=530, `(19,13)/(23,1)` | **1.90297 +/- 0.04242** | 1.61261 +/- 0.02644 | 0.29035 |

Here `U_T2=N^(13/8)P4[S']/Mbar'`.  Using the observed asymptotic slope scale
`Mbar'/N^(3/8) ~= 1.744`, these correspond approximately to

```text
N=370 P4[S']: Jordan 0.0019035, ordinary 0.0017016
N=530 P4[S']: Jordan 0.0013052, ordinary 0.0011061
```

These sizes and numerical targets were not part of the earlier P48/P64
artifacts.

## One-new-run rank-1 versus Jordan discriminator

For

```text
U(N)=N^(13/8)P4[S']/Mbar',
C_R=3U(N)-8U(2N)+5U(5N),
```

`C_R=0` exactly under `A+B/N`, but
`C_R=B log(3125/256)` under a Jordan log.  Existing clean parent/double points
therefore make either norm-5 child a one-new-run discriminator:

| closure | ordinary completion | Jordan completion | Jordan `C_R` |
|---|---:|---:|---:|
| N=65,130 -> **325** | U(325)=1.53875 | U(325)=1.71631 | 0.88779 |
| N=85,170 -> **425** | U(425)=1.86013 | U(425)=2.15558 | 1.47727 |

The N=425 lineage has the larger separation and is the more informative
single run.  The dual combination

```text
log(2/5)U(N) + log(5)U(2N) - log(2)U(5N)
```

annihilates the Jordan model exactly and is nonzero for `A+B/N`.

## Files

- `discovery.json`: covariance-aware fits, exact matrices, predictions, and
  discriminator values.
- `state_points.csv`: reconstructed intrinsic-coordinate state at every size.
- `scripts/discover_matrix_semigroup.py`: full replay from committed raw
  histograms.

