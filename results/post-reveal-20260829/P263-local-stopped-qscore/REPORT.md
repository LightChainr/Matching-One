# Issue #263 — localized stopped-transcript Q-score

## Outcome

This secondary replaces the extensive global `J/2` sample score by an unbiased
coupled estimate of its event-conditioned Doob projection. It reveals only the
clusters reached from the four marked boundary vertices, stops as soon as a
forbidden cross-group merge proves `14|23` false, and runs independent coupled
completions only when `14|23` is true.

The exact tiny oracle proves that irrelevant far-field score noise cancels. The
frozen 20k rectangle pilot verifies the executable estimator and records a
20.3%–22.7% mean revealed-edge fraction, but it does **not** yet demonstrate a
finite-rectangle variance reduction: its covariance trace is 31.2504 versus
29.8515 in the earlier, different-seed 20k global-score smoke. That comparison
is only a sensitivity check because the two runs are not paired and their rare
event counts differ.

## Exact identity and local increments

Let `X` be the outer critical bond configuration, `I` the `14|23` indicator,
and `T` an adaptive edge-reveal transcript that determines `I`. Write

`S = J/2 - E[J/2]`, where `J=2k+b`.

For an independent fair-bond configuration `C`, let `C<-T` overwrite only the
edges revealed by `T`, with their transcript values. Conditional on `T`,
`C<-T` has the law of `X|T`, while `C` retains the unconditional law. Therefore

`E[(J(C<-T)-J(C))/2 | T] = E[S|T]`,

and, because `I` is measurable from `T`,

`E[I(T)*(J(C<-T)-J(C))/2] = E[I E[S|T]] = Cov(I,J/2)`.

This is an exact unbiased identity for the covariance numerator. The reported
`d_Q log P` still divides by an empirical event probability and is therefore a
finite-sample ratio estimator.

The endpoint difference is also a transcript-edge telescoping sum. If
`C^(t)` overwrites the first `t` transcript edges, then

`(J(C<-T)-J(C))/2 = sum_t (J(C^(t))-J(C^(t-1)))/2`.

Each term is exactly `+1/2` or `-1/2`: adding an edge changes `J` by `+1` when
its endpoints were already connected and by `-1` when it merges two clusters;
removal gives the inverse bridge/non-bridge increment. Thus a transcript of
`r` edges obeys `|delta J|/2 <= r/2`, and every untouched far-field contribution
cancels sample by sample.

Equivalently, for the reveal filtration `F_t`, `M_t=E[S|F_t]` is the Doob
martingale and `E[S|T]=M_tau-M_0`. Coupled completions provide an unbiased Monte
Carlo estimate of this stopped martingale value without ever inserting the
unrevealed global score into the statistic.

## Event-determining exploration

The runner starts a FIFO exploration simultaneously from all four marked
terminals. It reveals only edges incident to vertices in terminal-connected
clusters. A connection between `{x1,x4}` and `{x2,x3}` is permanent, so it stops
immediately and returns false. If no such merge occurs, it exhausts all marked
cluster frontiers and accepts exactly when `x1~x4` and `x2~x3` are two distinct
clusters. Hence both early false outcomes and the final decision are measurable
from the stored transcript.

The current reference implementation computes the endpoint `J` difference by
building full DSUs for `C` and `C<-T` on positive samples. Its **statistical
support** is localized, but positive-sample computational work is not yet purely
arm-local; dynamic bridge-aware connectivity would be needed for that further
optimization.

## Tiny exact variance oracle

The oracle uses the four-cycle `14|23` event and appends disjoint fair spectator
edges invisible to that event. The exact target is `-1/256` for every row.

| spectator edges | global variance | ideal stopped variance | local K=1 / global | local K=8 / global |
|---:|---:|---:|---:|---:|
| 0 | 15/65536 | 15/65536 | 51.1333 | 7.2667 |
| 1 | 271/65536 | 15/65536 | 2.8303 | 0.4022 |
| 2 | 527/65536 | 15/65536 | 1.4554 | 0.2068 |
| 3 | 783/65536 | 15/65536 | 0.9796 | 0.1392 |
| 4 | 1039/65536 | 15/65536 | 0.7382 | 0.1049 |

The ideal stopped variance is exactly constant while the global variance is
`(15+256m)/65536`. One completion draw contributes exactly `47/4096` of
conditional noise, divided by `K`. This gives both the mechanism and its honest
boundary: with no irrelevant bulk (`m=0`), finite-`K` coupling is substantially
worse; the gain appears only after canceled extensive noise exceeds completion
noise.

## Frozen 20k mechanism pilot

The manifest froze level 1, 20,000 outer samples per geometry, 20 synchronized
batches, `K=8`, outer seed `2026102633`, and independent completion seed
`2026102634` before acquisition.

| lambda | events | dQ log P estimate | mean revealed / all edges | mean revealed on event | completion share of outer variance |
|---:|---:|---:|---:|---:|---:|
| 1/4 | 12 | -3.671875 | 0.2030 | 6541.8 | 9.9% |
| 1/3 | 21 | 2.428571 | 0.2083 | 5838.2 | 16.2% |
| 2/3 | 113 | -1.605088 | 0.2268 | 6109.3 | 15.4% |
| 3/4 | 209 | -2.025718 | 0.2272 | 5384.8 | 17.0% |

The unchanged amplitude-gauge target gives `chi2_3=6.5865424325`, with residual
covariance trace `31.2504005315`. The earlier global-score 20k smoke had
`chi2_3=6.7421914932`, trace `29.8515153094`, and event counts
`[19,19,133,210]`; it used a different RNG stream, so the trace ratio 1.0469 is
not a matched variance estimate. At this scale, rare-event count remains the
dominant limitation and `K=8` completion noise is secondary but non-negligible.

## Claim boundary and next mechanism test

- Exact: the covariance identity, transcript-edge telescoping, local bound,
  Rao–Blackwell dominance of the ideal conditional estimator, and `1/K`
  completion-noise law.
- Demonstrated on the tiny oracle: exact cancellation of spectator-area noise.
- Measured in the rectangle pilot: transcript fractions and completion-noise
  decomposition; no observed variance advantage over the unpaired 20k smoke.
- Conjecture: as rectangles grow, variance can change from bulk-area scale to
  the marked-arm/transcript scale. This requires measuring transcript growth
  across levels and can fail if marked explorations themselves become
  extensive.

The minimal next test is a paired runner that keeps the same outer transcript
and records both the stopped conditional estimator and its matched global
score solely for variance comparison. It should vary `K` and resolution while
leaving geometry, event, and tangent target fixed. It must remain a secondary
mechanism experiment, not a rescore of Phase E.
