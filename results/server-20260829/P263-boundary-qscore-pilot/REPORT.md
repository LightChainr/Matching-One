# Issue #263 frozen boundary Q-score pilot

## Outcome

Both acquisitions frozen in
`experiments/p263_boundary_qscore_phaseD_20260829.yaml` completed on the
16-core Huawei DevEnv while the existing N520 and N680 jobs retained their
eight cores.  The P263 runner used one core.  There was no predeclared
level-1 stop/go threshold: level 1 and level 2 are separate frozen scores and
are not combined as independent continuum evidence.

| frozen score | 14|23 event counts at lambda=(1/4,1/3,2/3,3/4) | active residual | chi2 / df | survival p |
|---|---|---|---:|---:|
| level 1, 200k/geometry | (129, 205, 1211, 1996) | (0.706182, -0.225028, -0.533449) | 1.097746 / 3 | 0.777618 |
| level 2, 500k/geometry | (147, 244, 1241, 2212) | (1.582687, 0.133978, -1.930542) | 5.872172 / 3 | 0.117998 |

Neither level rejects the frozen Cai high-branch tangent shape.  This is not
positive continuum evidence: the residual standard errors are
`(1.67084,1.10361,1.01192)` at level 1 and
`(2.53629,1.88263,1.67143)` at level 2.  The doubled-span run is less precise
despite 2.5 times as many samples.

## Exact acquisition audit

The executable is the dedicated square-bond runner at source commit
`fbd017408d0841d022aae425253799f778c38b09`, with binary SHA-256
`a08cd26a1763ebf803a9013b524c074d3f43073725bcd62a7db8beda169ebc55`.
No geometry, score formula, seed, or counter domain changed.
Preflight inspection found no existing P263 process/artifact and no reference
to seeds `2026102631` or `2026102632` anywhere under `/workspace`; the active
N520/N680 jobs used different seeds and counter domains.

| level | spans | seed | rows | batches/geometry | samples/geometry | sample counter domain | wall time | peak RSS |
|---|---|---:|---:|---:|---:|---|---:|---:|
| 1 | (15,14,15,14) | 2026102631 | 400 | 100 | 200,000 | [0,200000) | 2:07.86 | 2.9 MiB |
| 2 | (30,28,30,28) | 2026102632 | 400 | 100 | 500,000 | [0,500000) | 21:07.95 | 3.2 MiB |

All three Cai event totals are retained:

| level | lambda | 1234 | 12|34 | 14|23 |
|---|---:|---:|---:|---:|
| 1 | 1/4 | 3819 | 1882 | 129 |
| 1 | 1/3 | 4127 | 1417 | 205 |
| 1 | 2/3 | 4316 | 323 | 1211 |
| 1 | 3/4 | 5071 | 226 | 1996 |
| 2 | 1/4 | 3864 | 1879 | 147 |
| 2 | 1/3 | 4100 | 1430 | 244 |
| 2 | 2/3 | 4344 | 290 | 1241 |
| 2 | 3/4 | 5208 | 215 | 2212 |

The remote Python 3.9 runtime cannot execute the unchanged parent-oracle use
of `int.bit_count`.  The raw streams were therefore checksum-verified after
transfer and scored without code changes under local Python 3.13.7 with
mpmath 1.4.1.  This affected only where scoring ran, not the acquisition or
score definition.

## Post-reveal conformal-geometry secondary

Span doubling refines the fixed normalized rectangle
`[-2,4] x [0,4]`; it does not send the top and side boundaries to infinity.
The Schwarz--Christoffel elliptic map to the upper half-plane is

```text
z = sn((K/3)(w-1) | m),
K'(m)/K(m) = 4/3,
m = 0.21549970429193269707024331706636106943156632191124...
```

It sends the four declared Euclidean cross ratios to

| declared lambda | effective rectangle-to-UHP lambda |
|---:|---:|
| 1/4 | 0.23161561099460535818 |
| 1/3 | 0.31280861631947778502 |
| 2/3 | 0.64998175749754660991 |
| 3/4 | 0.73654763309312352114 |

The full secondary also applies the boundary-primary Jacobian and Cai
prefactor correction

```text
h' [sum_i log|f'(x_i)| + 2 log K(f(x_i)) - 2 log K(x_i)],
h' = sqrt(3)/(3 pi).
```

| acquisition | original frozen chi2 / 3 | effective-lambda only chi2 / 3 | full rectangle-conformal chi2 / 3 |
|---|---:|---:|---:|
| level 1 | 1.097746 | 1.163058 | 1.128359 (p=0.770233) |
| level 2 | 5.872172 | 5.945728 | 5.906372 (p=0.116255) |

Thus the coordinate correction is real but small at these marked points.  It
does not explain the scale-to-scale change.  This is explicitly a post-reveal
secondary and does not replace the frozen primary score.

## Mechanism-first interpretation

At Q=1, Cai's boundary weight is exactly `h=1/3`, so a doubled span predicts

```text
P(level 2) / P(level 1) = 2^(-4/3) = 0.3968502629920499.
```

The pooled observed probability ratios are 0.404223 for `1234`, 0.396466 for
`12|34`, and 0.434228 for `14|23`.  The first two closely reproduce the
ordinary boundary-field scale law.  Because the level-2 sample multiplier is
2.5, the predicted event-count multiplier is only 0.992126; the larger run
was designed to yield almost the same effective number of connectivity
events.  The primary `14|23` channel is 9.4% above the leading pooled scaling
prediction, a plausible channel-specific finite-box/lattice correction rather
than a resolved new exponent.

The tangent estimator is harder: it divides a covariance with the extensive
global score `J/2` by a rare-event probability.  Its covariance eigenvalues
grow from `(0.13385,1.02969,3.87010)` at level 1 to
`(0.52623,2.59096,9.65354)` at level 2.  The information-gain lesson is that
larger fixed-shape boxes plus proportional raw samples do not sharpen this
global-score tangent.  A future experiment should prioritize a declared
variance reduction/local score decomposition and a geometry whose conformal
cross ratios are frozen after mapping; this pilot itself cannot choose among
those post-reveal designs.

## Claim layers

- Exact finite/run facts: integer sufficient statistics, event counts,
  seeds/counters, resource profiles, checksums, and the elliptic rectangle map.
- Mechanism inference: ordinary probabilities follow the `h=1/3` scale law,
  while the global Q-score loses precision with span.
- Exploratory conjecture: the excess scaling ratio in `14|23` may be a
  channel-specific finite-box correction; the present covariance is too broad
  to identify its form.

The raw streams, primary scores, and post-reveal secondary are under `raw/`
and `analysis/`.  `checksums.sha256` covers every committed result artifact.
