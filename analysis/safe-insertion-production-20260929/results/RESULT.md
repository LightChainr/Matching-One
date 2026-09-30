# Prospective square-L512 safe-insertion score

Status: **scored**. Fixed 14 batches x 10000 prefixes; q=4.

All differences are early minus late. The point estimates are pooled plug-in estimates.

| Metric | Estimate | Aligned batch SE |
|---|---:|---:|
| deltaY | 0.006681826608 | 0.01208705518 |
| deltaE | 356.3977798 | 645.0460174 |
| deltaZ2 | -6.254031427e-08 | 1.131920089e-07 |
| secondary_deltaNu | -0.03986511616 | 0.5698779455 |
| weighted_earlyY | 0.09563676962 | 0.007712358034 |
| weighted_lateY | 0.08895494301 | 0.01060802794 |

Rows: 140000; rank counts: {'rank0': 42329, 'rank1': 55028, 'rank2': 42643}.
Rank-one prefixes: 55028; safe: 55028; dead (nu=m): 0.
Valid probes: 220112; nonzero probes: 1233; prefixes with any nonzero probe: 1215.

| Comparison / cohort | Eligible | Supported | Coverage of eligible |
|---|---:|---:|---:|
| primary / early | 24201 | 24109 | 0.9961985042 |
| primary / late | 30827 | 30354 | 0.9846563078 |
| secondary / early | 24201 | 24201 | 1 |
| secondary / late | 30827 | 30821 | 0.9998053654 |

Primary excludes dead prefixes because a safe-step mean is undefined there; secondary includes every rank-one prefix, including dead prefixes.

Within each shared exact (D,nu) cell: deltaE=(m-nu)*deltaY/2 and deltaZ2=-(m-nu)*deltaY/[m(m-1)]. The same normalized overlap weights combine all primary quantities. The exact-D secondary recomputes its own weights.

The JSON records exact integer sufficient statistics, all cells, support coverage, all 14 aligned deletion vectors, the full 6x6 covariance and source/input manifests. Probes are averaged within prefix; metrics are not independent and covariance is necessarily singular. No covariance inverse or independent-metric test is used.

- The primary tests one safe-step completion-creation mean at fixed exact (D,nu_b).
- Derived edge and two-step-survival contrasts use cell-specific s before aggregation.
- Secondary exact-D nu replication is independent of the old block, not of this primary.
- A null weighted moment can hide cell cancellation and does not establish recursive closure.
- Finite square L512 insertion clock only; no limiting, label-clock, or universality claim.
- No historical data pooling, subgroup search, fitted exponent, or acquisition is performed.
