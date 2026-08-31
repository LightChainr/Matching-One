# P154 prospective transmission: a one-coefficient map, currently a no-launch design

**Decision:** do not run another source-existence experiment or merely change
the old replay seed. A concrete one-coefficient transmission hypothesis can
be written down, but its nominal original-U effect is far below a bounded
fresh block's resolution. Keep this design as support unless a different,
explicit angular transmission map predicts a substantially larger effect.

## Current evidence and two distinct sources

The completed package `04743caf:experiments/p154-temporal-source-20260831/`
uses `L=K-ceil(sqrt(N))`. Its N260 source has entry response
`−.228771±.000879` in the first geometry, early-rank1 exit response
`−.047207±.000350`, but original-U derivative `+.84255±4.88218`.
It establishes temporal transmission to rank, not transmission to U.

The newer planning folder in the results-first-synthesis worktree instead
consumes **lag one**, committed as `4daae57e:results/norm4-lagged-source/`.
For that different source, N260/N340 total U derivatives are
`+.37897±.72219 / +1.27653±1.18227`. The planning-only N260 split is
entry `+.56268±.57048`, completion `−.18372±.55442`; neither sign is resolved.
They cannot supply honest opposite-sign production predictions.
Both sources reuse the old2.4M permutations; neither is fresh production.

## Exact source → births → U interface

Write F1/F2 for cumulative first/second birth, so
`q=F1+F2−1`, `E=1−F1+F2`. If their source derivatives are e/x, then
`Jq=e+x`, `JE=x−e`. For lag one specifically,
`e=T01+T02`, `x=T12+T02`; direct0→2 cancels from JE but not Jq.

All following quantities are at the same pooled matching root. Let

```text
D = mean(q'), B=P4[E'], H=P4[E''], T=mean(q''),
a = −mean(e+x)/D,
dD = mean(e'+x') + a*T,
U = N^(13/8) B/(2D),
v = N^(13/8)/(2D) * {P4[x'−e'] + a*H − (B/D)*dD}.
```

The entry/completion readout split uses this **same** a and dD in both
components. Pooled entry/exit magnitudes alone do not determine the angular
p-derivatives in this formula. No hazard-to-U proportionality is automatic.

## A minimal prospective hypothesis: two rigid birth shifts

Assume, as a new finite-window mechanism rather than an identity,

`e_g(p)=−tau1*F1_g'(p)`, `x_g(p)=−tau2*F2_g'(p)`

with the same two constants for both geometries and throughout the root
window. Define m=(tau1+tau2)/2 and delta=tau2−tau1. The common shift m
cancels **exactly** from root-comoving U. Only the single relative-shift
coefficient delta needs calibration for the U prediction.

Let `c=mean(E')/D`, `Q2=P4[q'']`, `E2=mean(E'')`. Then

```text
R = d_t mean(P_rank1) along the root
  = delta*D*(1−c²)/2,

K_rel = N^(13/8)/(4D) * [c*H−Q2 − (B/D)*(c*T−E2)],
delta = 2R/[D*(1−c²)],
v_clock = delta*K_rel.
```

Positive birth slopes imply `|c|<1`; the measured negative R therefore
implies a negative relative shift **under this hypothesis**. The sign of
v additionally depends on K_rel; entry suppression alone does not fix it.

At most two competing point predictions are proposed:

1. **Angular screening:** v=0 despite resolved pooled rank transmission.
   This is not the already rejected early-rank-only sufficiency claim.
2. **Rigid two-clock transmission:** v=delta*K_rel with delta calibrated
   from R and no fitted angular gain, offset or exponent.

These are distinct when delta*K_rel is nonzero. Neither is asserted to be
a complete physical theory. A future result outside both rejects both
operational predictions rather than justifying a new fit on the target.

## Existing-number budget: why this does not yet warrant production

The following are algebraic planning values, using the **source-root**
baseline jets already published in
`7da1eeb0:results/norm4-source-endpoint-1m/latest.json`, `by_N[N].source`.
The separate high-statistics `anchor` is at a different root and is not
silently spliced in. No histogram or old path was reprocessed.

| Source / N | delta | K_rel | Nominal v_clock | Old SE(v), 1M | Optimistic count for3-SE separation from0 |
|---|---:|---:|---:|---:|---:|
| lag1 /260 | −.0005042894 | −19.43667 | +.00980171 | .722186 | 4.89×10^10 |
| lag1 /340 | −.0003826928 | −45.79846 | +.01752674 | 1.182274 | 4.10×10^10 |
| ceil(sqrtN) /260 | −.0021455794 | −19.43667 | +.04170293 | 4.882175 | 1.23×10^11 |
| ceil(sqrtN) /340 | −.0018831736 | −45.79846 | +.08624645 | 9.921738 | 1.19×10^11 |

Budget formula is `M=10^6*(3*SE_1M/|v_clock|)^2`. It assumes unchanged
estimator efficiency and ignores training/prediction uncertainty: an
optimistic planning extrapolation, not a rigorous sample-size certificate.
An8M lag-one block has a3-SE resolution of about .766 at N260 and1.254
at N340. It does not distinguish these two predictions near .01–.02.

K_rel and delta have **not** received a new joint uncertainty calculation
here. The table is not a numerically frozen predictive interval, and its
positive point signs are not established mechanism signs. The remaining
scientific requirement is an independently declared, measurable map for
the angular p-jets—e.g. a contact-to-T01/T12 directional gain—not another
demonstration of large pooled source response. Dividing by weak U, arbitrary
source rescaling, or choosing the largest noisy entry/completion sign does
not remove this gap.

## Exact producer status and a conditional fresh protocol

The mature lag-one kernel and actual existing CLI are:

```text
4daae57e:src/norm4_lagged_source_replay.cpp
  norm4-lagged-source-replay N output.csv
4daae57e:scripts/replay_norm4_lagged_source.py
  python scripts/replay_norm4_lagged_source.py --workers 4
```

**These commands replay old data; do not launch them as a fresh experiment.**
Seeds, counters, sample counts, nested1000+9000 batches and sizes are fixed
in C++, and the driver exposes only workers. The longer-lag package also
hardcodes old counters and its run.py checks exact old-profile agreement.
The new planning folder currently contains an archive scorer, not a fresh
producer or frozen CONTRACT. No runnable new-seed CLI exists there yet.

If an alternative map passes the above power gate, the smallest honest
fresh implementation would reuse the lag-one event kernel and only:

- add explicit new seed/counter/count arguments and fresh output metadata;
  keep N260 and its exact periods `(16,−2;2,16)` / `(14,−8;8,14)`, both
  directions sharing each fresh permutation;
- collect source conditional means and unmarked q/E profiles on that same
  fresh block, or freeze a training centering table and retain its training
  uncertainty; never treat old profiles as the new sample;
- use fresh root/slope estimates and100 paired batch omissions for the
  original four-term v, with training coefficient uncertainty kept separate.

Required files are the C++ observer, immutable
`bfab0330:src/threshold_rank_integer_period_mc.cpp`, the prospective contract,
the frozen transmission coefficient/map and centering convention, and a
fresh driver/scorer. No new topology engine or general framework is needed.
But implementing that thin interface is not the current scientific blocker:
the tested nominal map has too little U separation for a bounded block.

This task modified no other team's files/processes, ran no sampling/replay,
and did not expand the archived covariance. The completed temporal result
was read with SHA256 `93f753a9a447773a1367fb728304e4a1028d83fc4d3c685bce506ded6280534c`;
the lag-one and longer-lag protocols were kept distinct throughout.
