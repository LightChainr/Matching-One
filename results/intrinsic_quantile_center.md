# Intrinsic quantile-center spectroscopy

Source: `scripts/intrinsic_quantile_center.py`.
Claim level: C0 definition freeze, C1 N=10 oracle. Not a P43/#57 target.

Frozen levels `u={0.025, 0.05}` and

```text
Q_N = c_{0.05} - c_{0.025}  ~  C N^{-3/4}
Q_{2N}/Q_N = 2^{-3/4}     on a true doubling lineage
```

Numeric doubling ratio `0.5946035575013605`.
Do not add quantile levels after looking at outcomes.

## N=10 Beta(3,3) oracle

The C4 self-matching control has `M(p)=2 I_p(3,3)-1`, which is odd about
`p=1/2`. Every intrinsic midpoint is therefore `1/2` and `Q_10=0` exactly.

```text
Q_10 = 0.0
c_u = 1/2 for both frozen u
```

## Descriptive P49 N=130/170

These two sizes are children of different doubling lineages, not a
doubling pair. Scaled `Q_N N^{3/4}` is reported as a development
diagnostic only. Do not score P43 against these numbers.

| N | `c_{0.025}` | `c_{0.05}` | `Q_N` | `Q_N N^{3/4}` | `w_{0.025} N^{3/8}` | `w_{0.05} N^{3/8}` |
|---:|---:|---:|---:|---:|---:|---:|
| 130 | 0.5927467678 | 0.5927427326 | -4.035242e-06 | -1.553557e-04 | 0.014309 | 0.028632 |
| 170 | 0.5927420089 | 0.5927386535 | -3.355397e-06 | -1.579722e-04 | 0.014319 | 0.028652 |

## Boundary

- Retrospective P49 numbers are development only.
- A claim-bearing score must recompute `p_±^u` inside each delete-one
  replicate and must be frozen before the target coordinates are read.
- P43 N=185/265 and Issue #57 norm-5 are not targets of this freeze.
