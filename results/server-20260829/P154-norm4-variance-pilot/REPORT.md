# Issue #154 norm-4 variance pilot and production allocation

## Decision

Promote the prospective norm-4 closure with a covariance-balanced allocation:

- extend each aligned source size `N=65,130,85,170` from 100M to 2B;
- run each independent target `N=260,340` at 1B;
- preserve the q2-first, Jordan-second score order and the complete source plus
  target covariance frozen on PR #202.

This is the lowest-CPU point in the tested grid expected to reach about three
sigma in both directional fixed-model tests. The N260/N340 pilot means were not
used to choose the allocation and are not interpreted as model evidence.

## Huawei pilot

The general-period engine ran on Huawei DevEnv
`f415a4bcbd9a438b85f5f29e4a507ea4` (16 vCPU, 32 GiB). Both sizes used 100
paired batches and independent RNG streams.

| target | samples | wall seconds | CPU sec / million | jackknife SE of U |
|---|---:|---:|---:|---:|
| N260 | 10M | 67.4057 | 53.9246 | 0.493102 |
| N340 | 10M | 85.9319 | 68.7455 | 0.920741 |

All four target stderr files (1M and 10M) are empty. The initial 1M run exposed
an analysis-floor problem: the generic double-precision full-curve solver
evaluated `p=0.9`, where the N340 binomial recurrence underflows. The dedicated
variance analyzer now computes only the required stable intrinsic-center slope
and `P4[S_prime]`; it does not import unrelated full-curve coordinates.

## Why target-only scaling was insufficient

With the existing 100M common-random-number source block, source uncertainty
dominates once target samples become large. Even 2B samples per target would
give expected square-root noncentralities only 2.019 (q2 null if Jordan is true)
and 1.567 (Jordan null if q2 is true).

The frozen 4x4 source covariance changes the useful allocation:

| total/source size | samples/target | q2-null sqrt(lambda) | Jordan-null sqrt(lambda) | incremental CPU sec |
|---:|---:|---:|---:|---:|
| 100M | 2B | 2.019 | 1.567 | 245340 |
| 1B | 2B | 3.891 | 3.473 | 310444 |
| **2B** | **1B** | **3.155** | **3.025** | **260111** |

The selected row spends more effort on the cheaper source sizes and less on the
expensive noncyclic targets. Its scientific justification is also independent:
it is the first prospective noncyclic `T4=T2^2` closure and the frozen
conjugation-even phase-node test, not merely a substitute for more norm-5 data.

## Reproducibility

- 1M generation commit: `b243d3a930ec637add9ff4f9f92de97af7b97250`
- 10M generation commit: `4e0561d4f93451d1e5a6c196ff0bfb77f20e1f8b`
- N260 pilot seed: `2026105401`
- N340 pilot seed: `2026105402`
- 10M disjoint counter interval: `[8101000000,8111000000)`
- source covariance block: seed `2026104501`, counters
  `[5000000000,5100000000)`, 100M per source size
- machine-readable result: `analysis/variance_forecast.json`
- frozen execution plan: `experiments/norm4_variance_pilot_20260829.yaml`

Production uses the disjoint counter intervals and sample counts in the frozen
execution plan. Pilot raw histograms, moments, metadata, stdout and stderr are
preserved under `raw/`.
