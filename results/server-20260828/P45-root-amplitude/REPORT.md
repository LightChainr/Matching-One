# P45: high-stat angular-normalized root amplitude

## Frozen test

The preregistered prediction is

```text
A_p = -N^2 DeltaRoot / DeltaCos4 = A_M / B
    = 0.4510066187069702 +/- 0.02013371335254959.
```

The prediction and primary sizes `N=65,85` were frozen before this run. The
score uses the same-batch cross-size jackknife covariance and adds the common
source-amplitude uncertainty to the residual covariance.

| N | A_p | jackknife SE | A_M | B | closure C |
|---:|---:|---:|---:|---:|---:|
| 65 | 0.4203381 | 0.0215723 | 0.736082 +/- 0.037778 | 1.751303 +/- 0.000045 | 1.0000781 +/- 0.0000217 |
| 85 | 0.3949486 | 0.0307810 | 0.691043 +/- 0.053854 | 1.749549 +/- 0.000041 | 0.9999122 +/- 0.0000218 |

Primary frozen score:

```text
chi-square = 2.42667 / 2 dof
```

The two sizes are also mutually compatible. A free common amplitude gives
`0.413007 +/- 0.019242`, with internal `chi-square=0.56497/1`. It is only
about `1.37` combined standard errors below the frozen prediction. The zero
effect is excluded with `chi-square=461.28/2`.

## Mechanism decomposition

The direct and linearized amplitudes agree:

```text
N=65: A_p direct=0.4203381, A_M/B=0.4203053
N=85: A_p direct=0.3949486, A_M/B=0.3949832
```

Thus the lower free-common point estimate is amplitude-level variation, not a
failure of root linearization. `B` varies by less than 0.2 percent from the
frozen reference, while the measured `A_M` values are moderately below the
P31 pooled fixed-p amplitude. Neither deviation makes the frozen joint score
fail.

The previous 10-million pilot had `A_p=0.2650 +/- 0.0751` and
`0.3575 +/- 0.0990`; its frozen score was `5.97/2`. The fresh 100-million run
shows why that pilot was not decisive.

## Execution and provenance

- Clean source commit: `6d2d68a62b433e337b97eadaf1870cb58e2f7666`.
- Source SHA-256: `251d1e7a74cb778ad38a11c6d5b79e86dc79701c0fff077af7f7bbf5109deefa`.
- Binary SHA-256: `1f918e8cd51f3ff8a5ec88eee6aac3a358d70da98086d18803bc583064901cbf`.
- Compiler/flags: GCC 10.3.1, `-O3 -std=c++17 -fopenmp`.
- RNG policy: deliberate common stream across `N`, aligned batches, full
  measured covariance.
- Seed/counters: `2026104501`, `[5000000000,5100000000)`.
- Sampling: 100,000,000 paired permutations per size, 100 batches.
- Wall time: 136.4 s (`N=65`) and 170.4 s (`N=85`) in concurrent 8-thread jobs.
- Exact tiny self-test passed. A separate 1-thread/4-thread 1,000-sample run
  produced byte-identical histogram and moment files.
- All raw artifact checksums were verified after download.

## Interpretation

The frozen angular-normalized radial amplitude test passes at the two primary
sizes and independently confirms a nonzero root-moving orientation sector.
Together with the parameter-free doubling result, this strengthens the
spin-4/`N^-13/8` mechanism. It does not establish a unique asymptotic exponent
or the `x=21/4` LCFT operator: secondary sizes and the corrected full thermal
projector remain required for that stronger claim.
