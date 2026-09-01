# Norm-5 local variance pilot

**Result boundary:** variance/power calibration only.  The sampled point estimates in this directory are deliberately not used as scientific evidence because 1M--10M replicas are far below the frozen production design.

## Outcome

The Apple-Silicon pilot supports the frozen 500M-replica allocation for the primary fixed-`p` matching-odd angular test:

| child | orientations | 1M SE of `Delta M` | projected 500M SE | H4 target | H12 target |
|---:|---|---:|---:|---:|---:|
| 325 | `(17,6)` vs `(18,1)` | `8.07305e-4` | `3.61012e-5` | `-4.98626e-5` | `+1.33951e-4` |
| 425 | `(16,13)` vs `(19,8)` | `8.34571e-4` | `3.73253e-5` | `-3.77114e-5` | `+1.01308e-4` |

Using only the child variances, the signed H4/H12 target separation is about `5.09 sigma` at N=325 and `3.73 sigma` at N=425.  Their independent-child quadrature is about `6.31 sigma`; parent uncertainty and the final covariance-aware score will reduce this.  This is consistent with the frozen maximin design estimate of about `5.60 sigma` for the joint H4/H12 comparison.

The same 500M curves are much less decisive for the secondary `P4[S']` q=2-versus-Jordan question.  From the 10M full-curve pilot:

| child | 10M SE of `P4[S']` | projected 500M SE | projected 500M SE of `Y=N^(5/4)P4[S']` | frozen q2/Jordan `Y` gap |
|---:|---:|---:|---:|---:|
| 325 | `2.16571e-3` | `3.06278e-4` | `0.42264` | `0.53296` |
| 425 | `2.11261e-3` | `2.98767e-4` | `0.57653` | `0.74010` |

Thus 500M gives only about `1.26--1.28 sigma` child-only separation per lineage for that derivative correction.  It remains a useful frozen secondary readout, but it is not a free decisive Jordan test.  Roughly 2B replicas per child would halve these child errors and move the two-lineage child-only separation to about `3.6 sigma` before parent uncertainty.

## Runtime

The two 10M jobs used five OpenMP threads each and ran concurrently.  N=325 took `116.30 s`; N=425 took `150.06 s`.  These timings are machine-specific and are not used in the statistical score.

## Reproducibility

- source commit: `bcbb135f55cdd94e46a614144b9f089c920cc3ca`
- source SHA-256: `7893216c66802b28eb67eb27ac61976835291c4ad734f94a0d255a3e6d7e179a`
- local binary SHA-256: `464b0398bad42c5edb910c17090884726d70c434e7af8d65011b95b64e0f8282`
- RNG: counter-derived SplitMix64 plus unbiased Fisher--Yates
- 10M counter intervals: `[8002000000,8012000000)` and `[8012000000,8022000000)`
- analysis precision: 60 decimal digits for the P48 derivative reconstruction

Raw histograms, batch moments, metadata and derived tables are retained beside this report.  See `commands.txt`, `environment.txt` and `checksums.sha256`.
