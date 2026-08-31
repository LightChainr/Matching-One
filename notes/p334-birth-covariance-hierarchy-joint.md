# Most joint-birth response is transport among fixed prefixes

The exact-score hierarchy localizes 96–98% of the normalized-rank covariance
response to relationships among prefix-conditional means. The prefix sampling
law itself is unchanged. The smaller within-prefix contribution is a response
of conditional covariance, not a fraction of the baseline noise variance.

The following signed response shares use the supplied twenty original-batch
LOO values. Errors are shared-batch SE, in percentage points.

| N / control-output | between prefixes / total | within prefix / total | between rank cells / total |
|---|---:|---:|---:|
| 325 / plus->S | 97.909 ± 1.341% | 2.091 ± 1.341% | 89.612 ± 2.176% |
| 425 / plus->S | 97.063 ± 1.830% | 2.937 ± 1.830% | 89.498 ± 2.759% |
| 325 / minus->D | 96.202 ± 1.645% | 3.798 ± 1.645% | 86.009 ± 2.458% |
| 425 / minus->D | 95.996 ± 1.870% | 4.004 ± 1.870% | 85.061 ± 2.172% |

The rank-cell split is consumed from `2bc35294`, without rerunning its
calculation. Even after rank cell is fixed, the between-prefix contribution
remains resolved: plus->S `3.69856e-8 ± 9.35027e-9` /
`3.19177e-8 ± 8.33132e-9`, and minus->D
`−1.33437e-7 ± 2.75886e-8` / `−1.06430e-7 ± 2.12545e-8`.
Thus rank-cell transport is the largest layer, with a smaller measurable
response of conditional means among prefixes sharing the same rank cell.
The supplied lifetime-variance cancellation coordinates are retained in the
complete factor; no additional cell readout is performed here.

## Same-source estimator change

Both estimators target the covariance derivative of
`X=K1/(N+1), Y=K2/(N+1)`. The old intrinsic result has removed uniform-clock
covariance. The new exact-score calculation uses the complete label census,
distinct-quartet products inside a prefix, and distinct-prefix global products.
The paired difference is computed from matching old/new LOO rows:

| N / direction | exact-score total minus old matched-mask intrinsic +/- paired SE |
|---|---:|
| 325 / plus->S | −1.33669e-8 ± 6.15837e-8 |
| 425 / plus->S | −1.12883e-7 ± 5.79992e-8 |
| 325 / minus->D | −1.14141e-7 ± 2.98280e-7 |
| 425 / minus->D | +1.15491e-7 ± 1.62115e-7 |

The largest paired shift is about 1.95 SE. These are estimator changes on the
same random block, not a new geometry, physical intervention or replication.
The response shares are signed ratios of derivatives, not probabilities;
outside this observed regime such shares need not lie between zero and one.

## Compact handoff

Source hierarchy: `44dc9e3396e39105cae85a29d04b39d0afc82d84`, using exact-score
moments `375cd3a1`. Previous shared covariance: `e2ef9983`. Rank-cell transport:
`2bc3529468fbcba589182acaf98fa4855eb0a85e`. All derive from the original
`e32a8593` fork and `959a7fa2` contact block.

`scripts/p334_birth_covariance_hierarchy_joint.py` reads only the supplied
committed summaries. It appends all 510 raw batch coordinates, the supplied
hierarchy/rank-cell derived LOO columns, and the eight named ratio/difference
LOO columns to the previous factor. Full covariance is `factor.T@factor`, with
rank at most 19 and no inverse. Outputs are under
`results/p334-birth-covariance-hierarchy-joint/`.

Scientific card: the joint-birth covariance response is now localized primarily
to transport of conditional means across existing rank cells and prefixes;
within-prefix response is smaller. This is one dependency block and one shared
covariance ledger. No new raw-path reads, MC, DP, tests, models or independent
evidence were added by this handoff.
