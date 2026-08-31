# P334: the N425 lifetime signal includes both lengthening and broadening

The full-population N425 `W²` direction is not explained by a change of mean
lifetime alone. Its squared-mean and variance terms point in the same direction:
the second orientation has a longer **and** more variable rank-one lifetime.
At N325 these orientation differences remain unresolved. This is a post-reveal
explanation of the already observed direction on the original archive, not a
fresh test, an inferred field count or a scaling law.

## Named decomposition on the complete original population

Set `W=K2-K1`, `C=(K1+K2)/2`. The operator below is exactly the original
orientation difference divided by `delta_cos4` (negative for these pairs).
At N425:

| Quantity | first | second | H4 difference +/- shared-batch SE |
|---|---:|---:|---:|
| E W | 18.29220 | 18.67795 | 0.432019 +/- 0.146872 |
| (E W)² | 334.60458 | 348.86582 | 15.97182 +/- 5.44701 |
| Var W | 225.55962 | 234.96563 | 10.53423 +/- 4.39429 |
| E W² | 560.16420 | 583.83145 | 26.50606 +/- 8.83952 |

The unnormalized mean-lifetime difference is `-0.38575 +/- 0.13114` ranks.
The exact population identity is

```
E W² = (E W)² + Var W.
```

The point contributions to the W² contrast are about 60% squared mean and
40% variance. These are descriptive shares with a common uncertain denominator,
not independently established fractions. Through
`J1_width=-E W²/[4(N+1)(N+2)]`, the two components are
`-2.19511e-5 +/- 7.48618e-6` and
`-1.44779e-5 +/- 6.03936e-6`; their covariance yields the previously read
total `-3.64290e-5 +/- 1.21487e-5`.

This mean-lifetime direction also has a simpler full-observer interpretation:
`integral E_top = 1-E W/(N+1)`. Thus its N425 H4 contrast is about
`-0.00101413 +/- 0.000344771` (the same approximately 2.94-SE direction,
not another independent observation). It does not require the first thermal
moment to become visible.

## Center spread, true ensemble spread, and the connected birth direction

The N425 H4 clock-center readouts are

```
E C:       -0.162420 +/- 0.187388
Var C:      0.439528 +/- 3.85172
Cov(C,W):   2.728043 +/- 2.58165.
```

The point center-spread difference does not oppose the lifetime contribution.
Its uncertainty is substantial. For the equally weighted normalized birth-rank
mixture `Y in {K1/(N+1),K2/(N+1)}`, the **ensemble** identity is

```
Var Y = [Var C + E W²/4]/(N+1)².
```

Its N425 H4 contrast is `3.89365e-5 +/- 2.58576e-5`: the total ensemble-spread
direction is not resolved simply because a component is. This variance is
distinct from energy around a fixed `p_ref`; the latter additionally contains
the squared displacement of the mixture mean from that reference.

There is one specifically connected, joint-birth coordinate:

```
Cov(K1,K2) = Var C - Var W/4.
```

At N425 it is `119.88467` versus `117.92562`, with H4 difference
`-2.19403 +/- 3.79450`. A coupling change is not resolved. Neither the lifetime
signal alone nor a fixed-reference total shape energy establishes a copula
change: the latter depends only on the two birth marginals. No further moment
families or best-window searches are introduced here.

## Source and reproducible covariance

- Full birth source: `9c495ab13e65f2bc93dc0849ee3b73f88724c4b1`, all 20,000
  paired paths per N, including R0/R1/R2. No conditional replacement is used
  in this population-moment decomposition.
- H4 normalization source: `3edc785a`.
- Script: `scripts/p334_lifetime_square_mechanism.py`.
- Full numerical result: `results/p334-lifetime-square-mechanism/score.json`.

The script first pools raw moments over all 20 batches. Each LOO replicate
drops one original batch of 1,000 paired paths and **recenters the remaining
19,000 paths** before recomputing squared means, variances and covariances.
It never averages batch-centered variances or batch products. The score saves
the ten raw moments per batch, all derived first/second/raw-difference/H4
coordinates, and their complete common LOO covariance. These empirical
plug-in descriptors use `ddof=0` to preserve the exact finite-archive identity.
No covariance inverse, new Monte Carlo, DP, replay or validation suite is used.

The analysis uses the dedicated local research Python environment:

```
/Users/lc/python-envs/research-py311/bin/python scripts/p334_lifetime_square_mechanism.py
```

The next handoff is the same-source join to the fixed-reference rank/canonical
shape energies and source-connected first-birth debt, keeping their shared
covariance rather than upgrading these reused coordinates to fresh evidence.
