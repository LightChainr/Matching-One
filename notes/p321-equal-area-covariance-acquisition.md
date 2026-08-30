# P321 equal-area rectangle covariance acquisition

Status: acquisition and scoring implementation plus a local variance smoke.
The smoke is not a test of the thermal-Q4 `E4` curve.

## Acquisition contract

At fixed area `N`, run the square as the first matrix against each of
`rho=16/9,9/4,4,9`.  All four invocations use exactly the same seed, replica
interval and batch boundaries.  The runner records a campaign manifest and
the scorer fails closed unless:

1. every metadata field in the CRN contract agrees;
2. the first and second period matrices equal the frozen design;
3. the repeated square histogram rows are byte-identical;
4. the repeated square joint-moment rows are byte-identical;
5. every histogram passes the existing exact moment audit.

Aligned delete-one-batch roots are then assembled in the fixed order

```text
(rho=1, 16/9, 9/4, 4, 9)
```

and used to estimate the complete `5 x 5` root covariance.  Rectangle-minus-
square contrasts and their `4 x 4` covariance are derived by the exact linear
map, not by treating the four pair runs as independent.

## Frozen scale scorer

Once all three equal-area scales are present, the only fitted scale law is

```text
p(N,rho) = pc + C_N(rho) N^-2 + D_N(rho) N^-3.
```

It is one GLS block with a common `pc`, five `C_N` values, five `D_N` values,
and blockwise full cross-rho covariance.  No exponent is fitted.  With
`N=144,576,1296`, this has 15 observations, 11 parameters and 4 residual
degrees of freedom.

The conversion to the transverse-width convention is applied coefficient by
coefficient and to its covariance:

```text
C_width(rho) = C_N(rho) / rho^2.
```

Only after that conversion does the scorer form the conditional ordinary
thermal-Q4 residual

```text
C_width(rho) - [E4(i rho)/E4(i)] C_width(1).
```

The primary conditional score uses `rho=16/9,9/4,4`.  `rho=9` remains an
endpoint diagnostic.  The E4 values come only from the already frozen oracle;
there is no smoke-driven model selection or amplitude fit.

## Local N=144 variance smoke

The local Apple-M4 smoke used 20,000 replicas per shape, 20 aligned batches,
one thread, seed `32114420260830` and counter interval `[0,20000)`.  The four
engine invocations took 1.162 seconds in total; covariance scoring took 12.74
seconds.  All four repeated square histograms and moment rows were byte-
identical.

| rho | root | root SE | rectangle-square | contrast SE |
|---:|---:|---:|---:|---:|
| 1 | .59253605 | 3.32e-4 | -- | -- |
| 16/9 | .59301893 | 3.85e-4 | +4.83e-4 | 4.76e-4 |
| 9/4 | .59251633 | 3.91e-4 | -1.97e-5 | 5.91e-4 |
| 4 | .59207124 | 3.30e-4 | -4.65e-4 | 5.65e-4 |
| 9 | .59158127 | 7.05e-4 | -9.55e-4 | 8.93e-4 |

The contrast correlations range from `0.306` to `0.808`.  This is the useful
smoke result: diagonal scoring would materially misprice the eventual shape
curve.  The numerical root differences themselves have order-one smoke z
scores and are not interpreted.

The scale scorer correctly reports `insufficient_scales` on this N144-only
archive; consequently no E4 score is emitted.

## Execution boundary

The N144 smoke demonstrates the entire covariance path.  N576/N1296 are not
authorized by this branch.  A later acquisition should preserve one raw
randomness coordinator and distinct counter domains across scales while
retaining aligned batches within each fixed-N shape family.

