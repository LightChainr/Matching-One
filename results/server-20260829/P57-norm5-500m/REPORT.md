# Issue 57 norm-5 production reveal

Two independent 500M-replica threshold-rank runs complete the frozen
N=325/N=425 Gaussian norm-5 experiment.  The primary typed score strongly
selects the H4 multiplier over its H12 adversarial alias.

## Primary fixed-p harmonic score

| model | chi-square / df | survival | delta chi-square from H4 |
|---|---:|---:|---:|
| H4 | 0.4163 / 2 | 0.8121 | 0 |
| H12 | 35.1931 / 2 | 2.280e-8 | 34.7768 |
| H8 | 16.0120 / 2 | 3.335e-4 | 15.5957 |
| zero effect | 1.7764 / 2 | 0.4114 | 1.3600 |

The two H4 residuals are `+0.645 sigma` and `-0.021 sigma`.  The H12
residuals are `-4.612 sigma` and `-3.770 sigma`.  This is a decisive norm-5
H4-versus-H12 discrimination, but the two target points alone do not reject a
zero child effect.  The result supports the frozen H4 alias choice; it is not
an independent discovery of a nonzero H4 amplitude.

## Frozen full-curve secondary score

The intrinsic functional cocycle remains unresolved:

| model | chi-square / effective df | survival |
|---|---:|---:|
| q=2 analytic (`c=8/5`) | 10.6482 / 6 | 0.0999 |
| rank-2 Jordan (`c=log(5)/log(2)`) | 9.0201 / 6 | 0.1725 |

Both fixed models survive.  Jordan improves chi-square by only 1.628 and is
not selected by this run.  These diagnostics reuse the primary raw curves and
are not an additional independent evidence block.

## Production and provenance

- environment: Huawei Cloud DevEnv `DevEnvC_ZyTrST`, id
  `f415a4bcbd9a438b85f5f29e4a507ea4`;
- platform: Huawei Cloud EulerOS 2.0, aarch64, 16 CPUs, 30 GiB visible RAM;
- production source commit: `fd6cfb647049bdaeb6bf1bb7a5f6a8a996507107`;
- engine: counter-derived SplitMix64, unbiased Fisher--Yates, same-permutation
  coupling within each orientation pair;
- N=325: counter interval `[10000000000,10500000000)`, 3192.15 s;
- N=425: counter interval `[10500000000,11000000000)`, 4103.82 s;
- both runs: 500M samples per orientation, 100 batches, 8 OpenMP threads.

Raw histograms, paired moments, metadata, stdout/stderr, typed score JSON and
SHA-256 checksums are retained in this directory.
