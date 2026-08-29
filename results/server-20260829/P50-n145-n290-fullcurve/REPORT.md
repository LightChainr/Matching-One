# N=145 -> 290 held-out full-curve score

## Outcome

The third Gaussian-doubling lineage validates the frozen scalar-plus-H4
finite-size correction to the center slope, but rejects closure of the whole
three-level thermal-even `DeltaM` curve under one fixed multiplier.

| frozen score | result | two-sided / survival p |
|---|---:|---:|
| joint `DeltaM`, `u=0,0.025,0.05` | chi-square = 9.35200 / 2 | 0.009316 |
| raw slope, `2^(3/8)` | z = -22.6903 | 5.58e-114 |
| scalar+H4 corrected slope | z = -0.66608 | 0.50536 |
| raw root ratio, `-1/4` | z = 1.38768 | 0.16523 |
| frozen induced root ratio | z = 1.38520 | 0.16599 |

The slope result is the clean positive finding: the two-parameter correction
trained only on `65 -> 130` and `85 -> 170` predicts the held-out
`145 -> 290` slope, while the asymptotic multiplier alone misses it by more
than 22 standard errors.  The root statistic is compatible with both frozen
targets at this precision and therefore does not discriminate them.

## Where the joint DeltaM failure lives

Both histograms preserve the declared genealogy order.  The observed central
contrasts are

```text
X_even(145,0) =  0.0005268040321
X_even(290,0) = -0.0002784936192
observed ratio = -0.52864747
frozen ratio   = -0.32420989
```

The central residual by itself is only `-1.385` standard errors.  The joint
rejection appears because the three intrinsic thermal levels are extremely
correlated and resolve a small shape deformation.

The residual correlation eigenvalues are

```text
1.23785e-13, 1.76308e-6, 2.99999824.
```

At the frozen relative cutoff `1e-10`, the numerical rank is two.  The common
amplitude-like mode contributes about `1.9243` to chi-square, while the active
shape mode contributes about `7.4277`.  The numerical null mode is excluded.
Thus the result does not overturn the earlier central/fixed-p H4 evidence; it
shows that a single multiplier does not transport the local thermal-even curve
shape.

## Correlated P4 diagnostics

These use the same raw block and are not additional evidence rows.

| channel | frozen multiplier | z | two-sided p |
|---|---:|---:|---:|
| `P4_S` | 0.5 | 0.67515 | 0.49958 |
| `P4_D` | 0.3242099 | 1.38514 | 0.16601 |
| `P4_S_prime` | 0.4204482 | 2.69536 | 0.00703 |
| `P4_D_prime` | 0.6484198 | -0.00909 | 0.99275 |

The near-exact `D_prime` transfer together with the failed `S_prime` transfer
is a sharper mechanism boundary than a generic scalar-correction failure.
It supports treating the derivative response as a small non-diagonal or
Jordan thermal block: one projected direction is already an excellent
eigenchannel, while its companion is not governed by the same scalar law.

The preregistered conjugation phase-node diagnostic in Issue #138 is not yet
identified with one of these real projected scalars.  It should be formed
from the two genealogy-resolved transfer responses with full covariance,
without selecting jet modes from the N290 result.  No new Monte Carlo block is
needed.

## Covariance and sign protocol

N145 and N290 use independent seeds.  Each size is jackknifed internally and
the transfer covariance is

```text
Cov(residual) = Cov_N290 + ratio^2 Cov_N145.
```

Both engines already emit `first - second` in frozen genealogy order.  The
negative norm-2 character belongs to the frozen target ratio and must not be
implemented by a second child sign flip.

The three-level covariance is effectively rank deficient.  The scorer uses a
correlation eigendecomposition and declares active eigenvalues relative to the
largest one; direct full-rank inversion would give a numerically meaningless
larger chi-square.

## Production provenance

- Huawei Cloud ARM64 DevEnv, 16 vCPU / 32 GiB
- raw execution commit: `d99a679dda985f5c65469cb25dbe46054f6e8b3c`
- N145: 100,000,000 replicas, 100 batches, seed `2026105003`, counters
  `[7000000000,7100000000)`, 8 threads, 290.71 seconds
- N290: 100,000,000 replicas, 100 batches, seed `2026105004`, the same counter
  interval in an independent stream, 8 threads, 705.83 seconds
- corrected scorer: PR #195

Machine-readable score: `analysis/score.json`.
