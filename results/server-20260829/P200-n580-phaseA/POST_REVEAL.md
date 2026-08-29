# N580 Phase A post-reveal interpretation

## Frozen joint score

The target state was

```text
(I_S, I_Du, T_D, T_Su)
= (-0.0122960443, -0.0162140264, -0.0859457053, 2.0354062688).
```

Using the full target-plus-prediction covariance frozen before acquisition:

| model | chi-square | df | survival |
|---|---:|---:|---:|
| ordinary q2 | 3.34640 | 4 | 0.50161 |
| rank-2 Jordan | 3.27859 | 4 | 0.51233 |

Both models are locally compatible with this N580 noncyclic Smith target, and
their joint scores are practically indistinguishable.  The largest correlated
diagnostic is the thermal-difference coordinate: `-1.46 sigma` for q2 and
`-1.55 sigma` for Jordan.  No coordinate reaches a two-sigma marginal deviation.

## Scientific reading

This result does not restore the ordinary-q2 model after the independent N170
norm-4 result rejected its scalar and thermal-jet predictions.  It says instead
that the single N145-to-N580 step has weak power to separate q2 from a first
Jordan correction.  Jordan remains compatible with both geometries and is the
relative surviving one-generator description, but N580 alone does not identify
it.

The informative discrepancy is therefore between paths, not between the two
N580 rows.  A universal factor-additive scalar law can look adequate on N580
while failing on N170.  The next observable should expose interaction between
the C2 and C5 quotient directions at fixed N650, rather than add another scale
point.  The typed mixed-join connected residual developed in #200 is designed
for exactly this distinction.

## Acquisition

```text
N=580
samples=100,000,000 per orientation
batches=100
threads=8
elapsed=23:31.34
peak RSS=9,444 kB
target commit=d9d6bb1a9bd9230139fef616a16280cad34d28f6
seed=2026102001
replica counters=[12000000000,12100000000)
```

The scorer commit is `510ed09f02e6d885f6d965b07e770b208b84ff7a`.

