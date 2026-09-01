# Production-first model elimination

Date: 2026-08-30. Integration status: **open PR #267**. Acquisition: **existing-data reuse only**.

## Late branch update: production now leads the queue

The first production adapter and the first adaptive intervention pilot are no longer missing. Branch `f5779b9` certifies eight real A/E rows and eliminates three fixed lines; branch `f54fb8c` gives positive P250 adaptive response at N325/N425. Conversely, `54b3e8` proves integrated Bernoulli energy is exactly E_top/2, `b239871` stops the radius-1 one-point singlet pilot, and `36bf1fe` leaves an A+C clock plane. The current executable order is therefore P250 N505 adaptive x spectrum, P334 shared-update nested forks, then a theory-led Phase-E mixed-correlator acquisition. See [`production-priority-pivot-20260830.md`](production-priority-pivot-20260830.md).

## Decision

Matching One has enough mature observables, covariance bundles and exact maps to
make production scoring the default next action. A tool is scientifically mature
when it can be applied to one already archived dependency block with a declared
model image and claim boundary. At that point another synthetic fixture or
semantic wrapper is optional infrastructure, not a prerequisite for analysis.

The working loop is:

1. use exact identities, selection zeros and rank/minor relations first;
2. score the surviving declared models on the complete production covariance;
3. reconstruct an outward interval or algebraic certificate when the stronger
   proof object would change the decision;
4. escalate only the still-live bounded classes to SOS or new acquisition.

This is an attention rule, not a permission or locking system. Exploratory models
may run in parallel, and an inconclusive certificate leaves a model open.

## What has actually reached production

The current boundary is uneven and should be read literally:

| requested layer | present state | exact boundary | next mechanism-changing output |
|---|---|---|---|
| Issue #275 Phase D, canonical `A_top/E_top` coordinates | **production-scored** on ten archived aligned direction pairs, including the P57 norm-5 N325/N425 blocks | full delete-one covariance and finite-model Mahalanobis exclusion; not a finite-sample theorem or a field identification | localize the surviving even response by activation, geometry and state current |
| P154 norm-4 production | **Phase-D production-scored** at N65/N85/N130/N170/N260/N340 | immutable PR #273 blobs supply `K1/K2`, canonical `A_top/E_top` and full aligned-delete-one covariance for the declared 1.9B/1B blocks; all six are K1-dominant and reinforcing, but `E_top=1-P1` is topology, not the energy operator | treat this as the topological arm of Phase E; do not repeat scalar/common-line scoring or call it an energy projection |
| Issue #275 Phase E, local energy versus clock plane | **integrated alias answered; radius-1 mean pilot stopped; mixed production still missing** | P154 has no local row or `B2,I0*B,I2*B`; `J_bulk,integrated=E_top/2` is not independent, and A+C currently survives | acquire a theory-led local B with mixed moments and ask whether it replaces C; do not extend the stopped mean-H4 row |
| Issue #370 exact certificate framework | **main framework plus first real branch-only production adapter `f5779b9`** | eight high-statistics rows eliminate `E=0,-A,+A`; a free ray survives. Exact margins remain conditional on the declared Gaussian outer set | reuse on the next new scientific covariance problem; do not build another fixture or adapter in front of acquisition |
| Issue #370 production elimination | **statistical Level-S applications exist** for E_top and P250, plus an observer-wise E_top outer audit | the outer audit is floating Gaussian Bonferroni/Fieller, not directed rational interval arithmetic, Positivstellensatz or SOS | escalate one bounded survivor only when the stronger proof object would change the survivor set |

In particular, the norm-5 request is no longer waiting for a synthetic harness:
its N325/N425 production rows already participate in the canonical E_top score.
The former immediate data-use gap in norm-4 is now closed by
[`results/norm4-two-activation-h4/latest.md`](../results/norm4-two-activation-h4/latest.md)
and
[`results/norm4-etop-production-elimination/latest.md`](../results/norm4-etop-production-elimination/latest.md).
Phase E remains open for a sharper reason than “more fitting is needed”:
canonical `E_top` is an Alexander-even rank-plane coordinate, not an energy
operator. The `J_bulk` alternative requires a third, separately observed
ordinary-energy/even-singlet row.
The common-line and fixed-power distances live entirely in the A/E plane and
cannot substitute for that missing observer.

### Issue #370 exact synthetic ladder

The former toy-fixture gap is now closed on main:

| integrated PR / main commit | exact witness | boundary |
|---|---|---|
| #411 / `7ecfaa3` | moments `(2,3,5)` give scalar Hankel minor `1`, excluding M1; `diag(1,2)` realizes M2d | one generator and three declared moments |
| #410 / `de2f44e` | `f_n=n+1` has a nonzero rank-two minor and a unique repeated-root recurrence, excluding minimal reachable/observable M2d; a 2x2 Jordan block realizes M2j | not all product-only or multigenerator models |
| #414 / `73ec55f` | endpoint moments remain M2j-compatible, but one typed morphism row gives stacked determinant `-1` and forces common predictive rank at least three | generic three-state realization, not physical deck/Smith/projector M3c |
| #413 / `9a94388` | `z=0`, `7z-1=0` has a primitive Bezout contradiction before optimization | one synthetic semantic-zero row |
| #412 / `58e01a4` | reachable-source minor `1/1024` exposes exact gauge amplification `1024` while responses, trace and determinant remain invariant | the zero-minor boundary and full chart atlas remain uncovered |

These are main-integrated exact controls, not five independent scientific results.
Main commits `dafc5e1 -> 2193eb7 -> cf4c4b5 -> c9b069a -> cc3c204` now add a
fail-closed envelope, a finite order-two recurrence compiler, exact Hankel-minor
enumeration, verification of a supplied univariate Bezout witness and verification
of a supplied rational realization. These are reusable arithmetic primitives,
but their canonical adapters are still tied to the synthetic fixtures: they do
not build arbitrary typed model problems, search witnesses, consume Level-S
covariance regions or constitute an implementation-independent verifier. The
next output is the production-confidence adapter for one real frozen block.
A separately implemented final verifier, directed intervals/SOS and a physical
M3c construction remain survivor-dependent follow-ups.

Later main through `5ac456d` adds exact Pell, mixed-curvature, product-Walsh,
Johnson-slice, acquisition-semantic and K-centered Euler controls. They are
design oracles, not new production rows. The semantic gate is concrete: P267
marked birth is Palm-like/path-adaptive; only retrospective fixed-K `O_ext`
supports Johnson, while `O_far` and the six-level high-pass are not scoreable
from current aggregates. PR #437 constructs the exact population high-pass but
does not create those coupled observations or a low-variance estimator.
The direct N16 branching counterexample `d9f813b` also shows that a complete
unbranched survival law can miss a common-update-then-clone gap. Apply these
facts to sharpen acquisition; do not use them to postpone `J_bulk`, P334
branching production or adaptive P250.

## First production application: this archive requires an even topological coordinate

The exact response coordinates are

```text
A_top = Delta4 F1 + Delta4 F2,
E_top = Delta4 F2 - Delta4 F1.
```

[`results/etop-production-elimination/latest.md`](../results/etop-production-elimination/latest.md)
uses all ten archived sizes and their complete aligned-delete-one covariance. It
does not generate Monte Carlo samples. The declared model distances are:

| declared model | result |
|---|---:|
| pure Alexander-odd response, `E_top=0` | `445.618/10`, `p=1.80e-89` |
| second-activation directional response zero, `Delta4 F2=0` | `182.905/10`, `p=5.84e-34` |
| first-activation directional response zero, `Delta4 F1=0` | `1041.049/10`, `p=2.68e-217` |
| one common projective line, `E_top=lambda A_top` | `28.593/9`, `p=7.59e-4` |
| one uncorrected fixed power, `E_top=c N^(-13/8)` | `37.482/9`, `p=2.16e-5` |

This requires an Alexander-even directional component and eliminates one common
rank-plane line across the declared archive. It does not identify a continuum
field, exclude corrected or multi-field H4 mechanisms, or prove that the finite
state dimension is two.

A deliberately more conservative observer-specific audit in
[`results/etop-production-model-certificate/latest.md`](../results/etop-production-model-certificate/latest.md)
uses separate 99% familywise Gaussian outer sets on eight high-statistics rows
and does not pool cross-dataset p-values. It still excludes `E_top=0`, `F1=0`
and `F2=0`, but leaves a common ratio in
`[-0.93056,-0.45518]`. The difference is useful calibration: the full-covariance
archive rejects the common-line parameterization, while the more conservative
no-pooling outer-set statement does not. They use different declared evidence
sets and certificate strengths, so they are calibration layers—not independent
votes or contradictory observers.

### Norm-4 Phase-D extension

The same exact coordinate transform now reads the six P154 archives directly
from immutable PR #273 Git blobs. N65/N85/N130/N170 share one declared aligned
counter group; N260 and N340 are separate groups. All six direction pairs are
scoreable. K1 and K2 reinforce at every size, all six are K1-dominant, and the
K2 point estimates at N260/N340 remain individually unresolved.

On the full saved covariance, the declared distances are:

| norm-4 model image | result |
|---|---:|
| `E_top=0` | `5324.015/6`, `p=2.85e-1150` |
| `Delta4 F2=0` | `3761.101/6`, `p=3.43e-811` |
| `Delta4 F1=0` | `21141.019/6`, `p=1.08e-4583` |
| one common `E_top=lambda A_top` line | `106.665/5`, `p=2.08e-21` |
| one uncorrected `E_top=c N^-13/8` amplitude | `177.528/5`, `p=1.80e-36` |

These results close the canonical Phase-D reuse on P154. They do not decide
Phase E. With only activation-derived A/E rows, `J_bulk` has no independently
measured ordinary-even coordinate. A valid comparison must add that observer
on the same aligned batches, retain its joint covariance with A/E and freeze
both common-transfer model images before scoring.

## Second production application: P250 is no longer a bridge election

The retrospective zero-sample augmented P250 gate combines the old radius-four
rows with the already acquired degree-five shifts while retaining separate old
and fresh delete-one influence covariance. All five pre-existing fixed bridge
maps and all eleven secondary maps reject in the same augmented Schur problem.

More importantly, the declared plus-hand `20 x 6` old-plus-degree-five
rank-at-most-five chart rejects at `p=6.01e-7`. The conflict therefore cannot be
assigned only to a wrong cross-hand bridge, and the earlier direction-only
extension cannot be promoted to general `5+5` support.

This closes the current R2-versus-R3 map-voting route. Exact endpoint audit
`af7dd01` proves that every present endpoint observable factors through word
abelianization, and `3128e3e` proves fixed-site delete/add overwrites commute or
absorb. Open PR #416 supplies the correct endpoint baseline: the complete
uniform-anchor N505 autocorrelation has at least 100 nonzero spatial Fourier
modes. Branch score `dbeb29c` then applies that structure to the archived
radius-4/5/6 production blocks; all four hand-charge channels remain compatible
with the nonnegative 101-frequency cone (`p=.163--.378`). Endpoint rank growth
therefore cannot be used as a hidden-state count.

The order-sensitive escape is already in production. Exact branch `6fbbe5e`
defines covariant state-dependent cut/join supports, and branch `f54fb8c`
finds positive response at N325/N425 with all 1,503 defined rectangles positive.
The next experiment is N505 adaptive x spatial-spectrum joint covariance. It
must report unconditional `I_defined*R_minus`, the conditional response,
defined/tie rates, all four selected supports, complement partners and joint
periodogram covariance. Another endpoint rank or shell vote has lower value.

A later `branch_only` locked radius-six Level-S certificate at `33c557b`
completes the candidate-independent rank ladder through eight. Endpoint-Hankel
`rank<=5,6,7` rejects separately in plus and frozen minus-R2 coordinates;
`rank<=8` is first compatible (`p=.1978/.1214`). This is a statistical lower
bound of eight in one observer/dependency convention, not exact state dimension,
rank-eight flatness, cross-hand closure, noncommutation or a field count. The R2
bridge was not reached because its rank-five support prerequisite failed.

## Third production application: the Euler response is only partly an occupancy clock

The N325/N425 external-Euler archives retain aligned batch-by-occupation
aggregates, so the global density nuisance can be removed exactly without new
simulation. With `n_occ` denoting the pre-insertion occupation count—not K1/K2—

```text
mu_ext(k) = k - 2N (k)_2/(N)_2 + N (k)_4/(N)_4.
```

The frozen radius-two local Euler nuisance has a corresponding exact conditional
mean on the declared locally injective period quotients. The scorer then centers
`J_D` and `J_S` within each `n_occ`, recomputes the same-next-site Gram
coefficient inside every aligned delete-one batch, and projects `J_S` only after
that centering.

[`results/euler-occupancy-clock/latest.md`](../results/euler-occupancy-clock/latest.md)
shows that the clock accounts for 53.14% and 56.16% of the raw far-D complex
magnitude at N325 and N425. A coherent 46.64% and 43.65% remains after
fixed-occupation centering and the JS projection; the joint zero-residual score
is `chi2=390067.1/4`, `log10 p=-84696.7`. The JS step itself retains more than
99.5% of the fixed-occupation residual.

This eliminates the declared “global occupation clock plus direct JS”
explanation. It does not identify a second microscopic source, Q4 field,
continuum exponent or independent evidence block: every row is a coordinate of
the same `external-euler-n325-n425` production group. The next discriminator is
K1 age/K2 completion/`DIRECT_RANK2` conditioning where retained fields permit,
then one genuinely new `n12` or winding/seam source—not another density transform
or reversible Gram rotation.

## Consequences for the next portfolio

1. **P250 adaptive morphism:** run N505 adaptive x 101-frequency spectrum on one
   covariance, with supports, typed intermediates, ties and complement IDs.
2. **P334/P337 continuation:** four proxy families fail. Move to checkpoint `b2`,
   one common update and two clones; keep risk composition, conditional hazard
   and `mark12_h4` separate.
3. **Norm-4 Phase E:** integrated energy is aliased and the radius-1 mean pilot
   stopped. Acquire theory-led `B,B2,I0*B,I2*B` only to test whether B replaces C.
4. **P333 typed radical:** charge one and charge two both fail at width 4. Retain
   rooted/landing connectivity after emission or couple irreps.
5. **Issue #370 production lift:** the first real adapter is complete. Reuse it
   only when a new observation creates a survivor-changing certification need.
6. **Continuum passports:** normalization-free three-point and modulus ratios
   remain high-value orthogonal discriminators after a lattice observable map is
   explicit. They should not delay production scoring of already typed finite
   observables.
7. **Threshold track:** main's exact degree-one height-100 exclusion removes
   low-height rational roots only. Degree 2--4 and constant-basis searches stay
   open, separate from operator identification.

## Literature boundary

- [Noncommutative polynomial optimization hierarchies](https://arxiv.org/abs/0903.4368)
  justify SOS escalation for a bounded, gauge-complete class; they do not turn a
  numerical relaxation or incomplete gauge chart into a certificate.
- [Emerging Jordan blocks in Potts and loop models](https://arxiv.org/abs/2403.19830)
  explains why finite-size diagonalizability is not a Jordan no-go.
- [The percolation energy field and its logarithmic partner](https://arxiv.org/abs/2508.16047),
  [exact loop-model three-point functions](https://arxiv.org/abs/2604.05503),
  [torus one-point functions](https://arxiv.org/abs/2604.24491) and
  [Temperley--Lieb local operators](https://arxiv.org/abs/2602.15742) motivate
  independent fusion, modulus and finite-module coordinates. None supplies the
  Matching-One lattice-to-field dictionary by itself.
- [Anchored random clusters and SLE excursions](https://arxiv.org/abs/2605.04395)
  supplies bulk-boundary anchored and pivotal observables in a half-plane SLE
  setting; it motivates boundary-conditioned readouts but does not identify the
  torus occupancy-clock residual or its birth-age kernel.
