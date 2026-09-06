# Freezing the state once: #580's intervention-transport control on P398

**Date:** 2026-09-06
**Issue:** LightChainr/Matching-One#580
**Script:** `scripts/p398_intervention_transport.py`
**Artifact:** `results/p398-intervention-transport/latest.json`
**Tests:** `tests/test_p398_intervention_transport.py` (16)
**Compute:** deterministic. No Monte Carlo, no width increase, no mark search.

## The question, and why P398 can answer it exactly

Same-environment low rank is cheap. #580 asks the harder thing: freeze the state
functions **once**, at `eta = 0`, and predict a *different* environment by
changing only the generator. P398 is the right place to ask because the
intervention is not estimated — it is written down.

On the canonical noncrossing connectivity states of width `w` (issue #11's
frozen codec) each boundary point carries two competing moves. Tilting their
rates by one scalar gives

```text
G_eta = (1+eta) * sum(join_cyclic_adjacent) + (1-eta) * sum(detach) - exit diagonal
      = G_0 + eta H
```

with `H` exact. State counts are the Catalan numbers: 14, 42, 132, 429, 1430 for
widths 4 to 8.

## What was frozen before any number was computed

| frozen at `eta = 0` | allowed to move with `eta` |
| --- | --- |
| the state functions `Phi` | the microscopic generator `G_eta` |
| the readouts and their coordinates `Phi^T f` | the stationary law `pi_eta` (reported, never used to choose a function) |
| the sources and their coordinates `Phi^T mu` | |
| the counting-measure inner product | |
| the reduced tangent `B = Phi^T H Phi` | |

Also declared before the run: the `eta` ladder `{0, ±1/8, ±1/4, ±1/2}`, the lag
grid `{0.25, 0.5, 1, 1.5, 2, 3, 4}`, the split of the eight readouts into three
that seed the span (`blocks`, `singletons`, `wrap`, plus the constant function)
and five that never touch it (`max_block`, `linked_pairs`, `halves_linked`,
`covering_depth`, `boundary_span`), the contrast reference source, the 0.10 pass
threshold, the 3x control margin, and the three intervention directions.

## Two things worth stating before the numbers

**Galerkin compression of an affine family is automatically affine.** For a
frozen span, `Phi^T G_eta Phi = A_0 + eta B` identically — measured drift
`8.2e-14` across every block in the run. So Contract II adds nothing to
Contract I *unless* `A_eta` is allowed to see the target. That is what separates
the two rungs here: `M_fixed_span` is scored at its best possible in-span
prediction at `eta`, `M_transport` at `A_0 + eta B` with nothing refit.

**`eta = 0` belongs in the ladder.** At `eta = 0` the transported model *is* the
baseline Galerkin model, so its error there is pure rank truncation. Every
finite-`eta` number has to be read as an excess over it, or the truncation gets
reported as a transport failure. Adding `eta = 0` is what turned this reading
from "6.8% error" into "0.6% excess".

## Result 1 — transport is nearly free, and it is not a pencil artifact

Width 8 (1430 states), Krylov span, declared intervention, relative Frobenius
error on the source-contrast response tensor over the three span-building
readouts:

| rank | `eta=0` | `+1/4` | excess | `-1/4` | `+1/2` | `E_refit` | `E_fixed_span` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 0.1893 | 0.1955 | +0.0062 | 0.1706 | 0.1976 | 0.1856 | 0.1856 |
| 4 | 0.1560 | 0.1581 | +0.0021 | 0.1387 | 0.1550 | 0.1471 | 0.1471 |
| 6 | 0.0625 | 0.0683 | +0.0057 | 0.0506 | 0.0696 | 0.0510 | 0.0510 |
| 8 | 0.0262 | 0.0261 | -0.0001 | 0.0296 | 0.0283 | 0.0198 | 0.0198 |
| 12 | 0.0060 | 0.0174 | +0.0113 | 0.0116 | 0.0318 | 0.0054 | 0.0097 |

The excess is stable in system size. At rank 6, worst `|eta| = 1/4`:

| width | states | `eta=0` | worst `1/4` | excess | refit span drift |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 4 | 14 | 0.0246 | 0.0289 | +0.0043 | 7.0e-16 |
| 5 | 42 | 0.0348 | 0.0398 | +0.0050 | 4.6e-16 |
| 6 | 132 | 0.0446 | 0.0497 | +0.0051 | 1.4e-15 |
| 7 | 429 | 0.0539 | 0.0592 | +0.0053 | 2.2e-15 |
| 8 | 1430 | 0.0625 | 0.0683 | +0.0057 | 1.2e-14 |

A random-direction span of the same rank sits at 0.997 at width 8. The
realizable least-squares generator fit, given its best shot over five
pseudo-inverse cut-offs, lands at 0.0535 — it does not beat the frozen
transported model by enough to matter.

**The pencil control comes back negative, which is the interesting part.** The
baseline generator is `J + D` and the declared intervention is `J - D`, so the
intervened generator never leaves the two-dimensional pencil `span{J, D}` that
built the state space. Measured directly: one step of that intervention lands
the declared readouts back inside the Krylov flag to `1.8e-14`, and the refit
span does not move at all (`1.2e-14`). A transport pass under such an
intervention would be an algebraic tautology.

So a third direction was scored: tilting the join at **one** boundary point.
That one drags the declared readouts out of the flag by `0.71` and does move the
refit span (`9.0e-02`). It transports anyway — excess `+0.0016` at rank 6, if
anything *better* than the in-pencil directions. Transport here is not an
artifact of where the intervention sits relative to the algebra.

## Result 2 — the state is a state of its readouts, not of the system

The same frozen realization, scored on the five readouts it never saw:

| readout | in span? | transported error, width 8, rank 6, `eta=+1/4` |
| --- | --- | ---: |
| `singletons` | yes | 0.042 |
| `blocks` | yes | 0.078 |
| `max_block` | no | 0.124 |
| `wrap` | yes | 0.230 |
| `linked_pairs` | no | 0.267 |
| `boundary_span` | no | 0.314 |
| `halves_linked` | no | 0.522 |
| `covering_depth` | no | 0.726 |

This is a **baseline representation** failure, not a transport failure: the
held-out block sits at 0.240 already at `eta = 0` and the transport excess on it
is only `+0.025`. Rank 12 brings it to 0.141, improving slowly.

The ordering is interpretable. Additive local counts of the partition are
captured; nesting depth and long-range boundary linking are not. `wrap` is the
warning inside the good news: it is a declared readout, and it fails at 0.230
while the pooled declared-block score reads 0.068, because the pooled relative
Frobenius norm is magnitude-weighted and `blocks`/`singletons` are `O(w)` while
`wrap` is `O(1)`. The artifact therefore carries a second, explicitly post-hoc
verdict computed one readout at a time at the same threshold; by that reading
nothing passes.

## Result 3 — the slow-mode subspace is not the state

The other frozen span family is the dominant invariant subspace of the
uniformized baseline generator — the span that closes *exactly* at `eta = 0`
(measured closure `4.6e-04` to `7.6e-03` at width 8 after 400 orthogonal
iterations, and `2.6e-16` at width 4 where it converges fully). It fails at
everything: transported error 0.20 to 0.28 at width 8 across all ranks, and
0.14 to 0.22 even when the span is **refit at each `eta`**. A "slow-mode state" is simply not the object
these observables need.

## Exact controls

All passed, at every width where they are affordable and for all three
interventions:

- rows of `G_eta` sum to exactly zero in rationals; no negative off-diagonal
  entry for `|eta| <= 1`; `G_eta = G_0 + eta H` exactly; chain strongly
  connected (widths 4 and 5, all `eta` in `{-1, -1/4, 0, 1/4, 1}`);
- `exp(t G) 1 = 1` to `4.4e-16`;
- uniformization matches an independent dense Taylor exponential to `9.8e-15`;
- the stationary power iteration matches an exact rational Gaussian elimination
  to `4.1e-14`;
- a full-rank span reproduces the intervened generator exactly (drift `0`).

Moves that do not change the state are dropped rather than recorded: keeping
them would cancel in the generator but not in the uniformization rate, shifting
every fitted timescale.

## The decision

Against #580's table, at width 8, rank 6, `eta = ±1/4`:

- **Contract I (fixed-function closure): passes.** `E_fixed_span` equals
  `E_refit` to the digit for both in-pencil directions and to `6e-4` for the
  out-of-pencil one. Low rank here is *not* environment-specific.
- **Contract II (tangent transport): passes.** Excess over the `eta = 0`
  truncation is `+0.0057` at `±1/4` and `+0.0071` at `+1/2`.
- **Contract III (fixed realization on held-out readouts): fails.** 0.264
  against a 0.10 threshold, and the failure is already present at `eta = 0`.

**Verdict: `STATE_TRANSPORTS_ON_ITS_OWN_READOUTS_ONLY`.** This branch is not in
#580's table, which assumes the declared and held-out readouts stand or fall
together. They do not, and the branch was added rather than folding the case
into the refinement branch — the one added intervention-generated direction
lowers the error by 0.003 and closes nothing.

## What this means for the wider programme

The finite low-dimensional state is real, transportable and cheap to transport —
**relative to a declared observable dictionary**. It is not a state of the
system. Change the dictionary from additive counts to nesting or long-range
linking and the same frozen span stops representing the dynamics at all, before
any intervention is applied.

That is a direct input to #275. A "mechanism class" defined by a low-rank state
is identified only up to its observable dictionary, so two mechanism classes
that share a dictionary will share an image whatever their microscopic content.
The relevant question for #275 is therefore not "does a low-dimensional state
exist" — here it does, and it survives a controlled mechanism change — but
"which dictionary is the claim about". Two models scored on the same dictionary
returning the same image is the expected outcome, not a surprise, and should be
read as `UNIDENTIFIABLE_WITH_CURRENT_ASSETS` unless the dictionaries themselves
separate.

It also sharpens #581's bulk-versus-topology split with an independent piece of
evidence: on this exactly solvable process, the readouts that a compact
low-rank state captures are exactly the additive bulk-like ones, and the ones it
misses are exactly the topological ones (`covering_depth`, `halves_linked`,
`wrap`).

## Boundary

Deterministic P398 calculation only. Intervention invariance is a falsification
gate, not latent-state identification (Bing et al., arXiv:2312.03580), and this
is a statement about the noncrossing planar process, not about square-site
Matching One. No square-site application is opened by this note: that requires a
declared intervention with common observables before and after, and a supplied
map between microscopic state spaces.
