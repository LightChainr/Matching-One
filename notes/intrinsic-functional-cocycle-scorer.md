# Intrinsic functional-cocycle scorer

This implements the pre-target full-curve discriminator proposed in Issues
#101, #119, #125 and #138.  It is deliberately narrower than a new scaling
fit: the unknown leading function and the unknown first correction function
are eliminated algebraically at each of the already-used intrinsic levels
`u={0,.025,.05}`.

For each size, the scorer solves `Mbar(p)=+/-u` inside every delete-one-batch
replicate, reconstructs

```text
T_N(u) = D_even(u,N) + S_odd(u,N),
Z_N(u) = N^(13/8) T_N(u),
```

and scores

```text
R_c = Z_(5N) - c Z_(2N) + (c-1) Z_N
```

for `c=8/5` (ordinary relative-q=2 correction) before
`c=log(5)/log(2)` (rank-2 Jordan cocycle).

The production covariance boundary is explicit.  Existing N65/85/130/170
curves share counters and form one synchronized jackknife group.  The Huawei
N325 and N425 runs use disjoint counters, so each is an independent group.
The scorer sums the three group covariance contributions rather than creating
spurious covariance by aligning unrelated batch numbers.

Example after the two target runs are available:

```bash
python3 scripts/score_intrinsic_functional_cocycle.py \
  --histograms n65.hist.csv n85.hist.csv n130.hist.csv n170.hist.csv \
               n325.hist.csv n425.hist.csv \
  --covariance-groups 65,85,130,170 325 425 \
  --json functional-cocycle-score.json
```

The same run also reports the intrinsic center diagnostics
`J=P4[S']/Mbar'`, `N^(13/8)J`, and `Xi=J/P4[D]` with their joint covariance.
Those diagnostics test whether the bare-`p` thermal metric explains the
center `S'` drift; they do not alter the frozen functional score order.
