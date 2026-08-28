# P31 independent same-N orientation confirmation

The Huawei ARM server completed the frozen five-size confirmation with
100,000,000 paired replicas per size, 100 batches, seed `2026093001`, counters
`[1000000000,1100000000)`, 16 CPU threads, and
`p_ref=0.592746050790`.  Total wall time was 525.8 seconds for 500 million
paired replicas.  Orientation order is first minus second throughout.

| N | pair | Delta M | batch SE | z | scaled A4 |
|---:|---|---:|---:|---:|---:|
| 65 | `(8,1)` / `(7,4)` | +1.24948e-3 | 7.80e-5 | 16.03 | 0.8093 |
| 85 | `(9,2)` / `(7,6)` | +1.01189e-3 | 9.01e-5 | 11.23 | 0.8666 |
| 130 | `(11,3)` / `(9,7)` | +4.67000e-4 | 8.95e-5 | 5.22 | 0.9330 |
| 145 | `(12,1)` / `(9,8)` | +4.42250e-4 | 8.38e-5 | 5.27 | 0.7501 |
| 170 | `(13,1)` / `(11,7)` | +2.37640e-4 | 9.20e-5 | 2.58 | 0.6277 |

Every sign agrees with `Delta cos(4 theta)`.  In particular, the independent
seed reproduces the N=65 and N=85 signs at 16.0 and 11.2 standard errors.
Against the earlier 30-million seed, the seed-to-seed differences at
N=65,85,145 are only 1.32, 1.38, and 0.53 combined standard errors.

Pooling seeds within size gives scaled amplitudes

```text
N=65   0.7812 +/- 0.0458
N=85   0.8139 +/- 0.0671
N=130  0.9330 +/- 0.1788
N=145  0.7201 +/- 0.1307
N=170  0.6277 +/- 0.2430
```

The common-amplitude summary is `A4=0.7885 +/- 0.0352` with chi-square
`1.53` for 4 degrees of freedom.  Thus all five sizes are statistically
compatible, while the visible drift is retained for the frozen model
challenge rather than interpreted away.

The run used source SHA-256
`dc54f6a2ae8ecf9e9cca630f33b183248e750b56e191d6b7c96d611eba0ee548`.
The metadata git field records the pre-commit integration point `34b1824`;
repository history in this PR contains the identical source file.  Raw batch
counts, full sector analysis, covariance/CRN diagnostics, and the P31/P32
interface CSV are retained in this directory.
