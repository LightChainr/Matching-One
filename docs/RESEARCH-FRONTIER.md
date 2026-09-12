# Research frontier: one active paper, two visible reserves

**2026-09-13, owner-delegated reset reconciled with the latest #650.**
The single active deliverable is a probability-theorem manuscript and independent
proof audit: #735's arbitrary-period root theorem, with #613/#736's sharp axial
full-law boundary. Return one acceptance package to #735. The other directions
below remain scientifically visible but are not parallel active assignments.

[ROADMAP](ROADMAP.md) controls execution; live Issue updates override this dated
snapshot. [STATUS](STATUS.md) distinguishes current scope from the
[verbatim historical ledger](history/STATUS-before-20260913.md). The baseline
read was main at `eb89e9422791d9e3c3a78f0e65d56912b815a7bd` (#702--#705 merged).
Navigation was integrated by #738; the research PRs below remain unmerged.

## 1. Active probability paper: correct root does not imply a concentrated law

For independent NN square-site percolation on an honest torus,
`M=P2-P0`, `F=E[r]/2=(1+M)/2`, `Q=F^{-1}`. The balance root is Q(1/2).
The conditional rare-sector function `H=P2/(P0+P2)` is not F.

#735 supplies arbitrary-integer-period root consistency as the genuine shortest
period ell grows, with no aspect or area-versus-systole restriction. At fixed
subcritical p it compares rates: `P2/P0 <= exp[-kappa(p) N/ell]`.
#736 adds the exact axial contrast: for 2<=w<=m and wm->infinity,
all fixed interior Q(u) converge to p_c **iff** `log(m)/w -> 0`.
These are proofs under named probability/topology inputs, not a new numerical
critical point, a fitted root-shift exponent or an established priority claim.

The new contribution in this reset is #736's necessity argument. Critical
square-site RSW and finite-product continuity give a seam-closed occupied ring
with probability >=exp(-eta*w) at a FIXED p_eta<p_c. Independent transverse
bands amplify it. On a subsequence with log(m)>=d*w this keeps every fixed
lower quantile away from p_c, while the median can remain consistent. Failure
of a union upper bound alone would never have proved this implication.

**Current acceptance work:** consolidate the theorem statement and imported
inputs, independently audit #735's local support/injectivity, oblique entry
count, backtracking-path extraction and disjoint-band independence; audit #736's
seam and quantifier order; compare the exact statements with their closest
strip-percolation, RSW and homological-percolation predecessors. Deliver one
manuscript, not another growing list of review tickets. Finite code checks and
verbatim quotations are not substitutes for the all-size argument.

**Potential strengthening, not a new prerequisite or second dispatch:** build
uniform fixed-subcritical winding corridors around an arbitrary shortest
integer period u, with probability >=exp(-eta*|u|), then pack disjoint corridors.
This could establish necessity of `log N/ell -> 0` beyond axial tori. It is NOT
proved here. Orientation and ambient nonprimitivity must be handled on the
physical NN lattice; changing a period basis does not rotate its interaction.
A counterexample could instead reveal a finer geometric invariant. Neither
outcome invalidates the already supplied narrower results by implication.

#736's executed controls cover 21,760 exhaustive strip configurations and 2,365
deterministic uneven-cell masks, with exact Fraction inequalities, three local
mathematical tests and compilation. Full repository CI was not run. The theorem
rests on its argument and imported inputs, not on extrapolating these checks.

## 2. Representation reserve: sources determine which states can be merged

#708 gives lifted rank-future states and 509 deterministic continuation classes
at width four. #710 gives the all-p scalar spectrum and exact cancellation of
common P0/P2 modes in M. #733 gives source-compatible quotients, physical two-row
site responses and exact final-rank-conditioned backward sampling. Two adjacent
independently addressed columns require all 509 classes within the investigated
common-strong-lumping class. These completed calculations are not assignments
to repeat merely because their PRs are absent from main.

The distinction worth pursuing after reactivation is an actual all-width
closure-weight theorem or a source-faithful response separating concrete
candidate representations. Deterministic classes, strong lumpings, scalar
Hankel order, positive order and noise-limited effective order are different.
A full-operator spectral gap need not be the gap visible to M. P398 is the
known periodic IC O(1) TL calibration process (#718), not square-site percolation.

### Concurrent #737: a concrete bridge to spatial structure factors

The spatial-Hessian and invisible-mark notes were read at the pinned head below.
For a mean-zero logit source h, H=sum_i h_i n_i and occupation count K, the
finite logit-root curvature is

    z_*''(h) = -[Var(H|r=2)-Var(H|r=0)] / [E(K|r=2)-E(K|r=0)].

This identifies a difference of rank-conditioned spatial structure factors,
not another uniform thermal derivative. Additive probability and additive logit
fields have different second-order corrections. The reported 4x4 additive-root
Hessian has both signs; no infinite-volume disorder relevance follows.

The invisible-mark example preserves the entire uniform (r,K) law but separates
under spatial sources. It uses positive DEPENDENT marked measures, not an
ambiguity about the already specified Bernoulli law or the two original-U fields.
A separate orbit calculation in this reset reproduced sizes 16/128/64, contrast
matrix `[[-512,-128],[0,-256]]` and determinant 131072. This checks that finite
algebra, not the full Hessian or complete CI. #737 discloses three missing
delivery artifacts; resolve actual `delivery/` paths and reproduction inputs
before treating archive execution logs as a fresh-checkout installation.

### Exact class growth need not be robust at finite precision

#549's k+1 branching classes have identical full unbranched survival. In its
specified single-fork readout, adjacent separation is `1/[2k(8k-1)^2]`, while
the whole family's span is `1/[2(8k-1)^2]`. More exact classes coexist with a
shrinking observable range. This is not a bound for every branching experiment.
The reserve question is whether budgeted admissible interventions amplify the
separation, not whether larger k gives a larger exact class count.

## 3. Theory-blocked flagship: original U, not a substitute observer

#275 remains P1 with UNIDENTIFIABLE_WITH_CURRENT_ASSETS. This does not disprove
Jordan or equate two physical theories. The missing input is two candidate maps
through the SAME source, six-coordinate thermal jet, physical normalizer,
rank-one denominator and pooled moving-root counterterm. Establish baseline
representability, then compare nuisance-profiled prediction images using existing
covariance. More precision cannot supply unspecified theoretical columns.

#735's positive irreducible stochastic controls share all ordinary trace jets
while differing in Jordan structure. Actual visibility depends on
`C N_lambda^k P_lambda B`. A derivative-jet lift can have its own nilpotent part;
its repeated pole does not identify the physical operator. Spin, a power and a
logarithm are not substitutes for the source/readout map. #737 does not supply
the two missing continuum candidates.

## Pinned unmerged assets

| PR | Head | Role |
|---|---|---|
| [#708](https://github.com/LightChainr/Matching-One/pull/708) | `f782061c1a592ed2f9fd0e9dabaa45f0e54bc4e7` | Finite rank closure |
| [#710](https://github.com/LightChainr/Matching-One/pull/710) | `d543ba1afe052a36682d2c1723c8a07281c28950` | Parametric spectrum and geometric onsets |
| [#718](https://github.com/LightChainr/Matching-One/pull/718) | `72d2ee6b8e5d25bd96aaef586c113bc94201f64a` | Axial root theorem and IC TL identification |
| [#733](https://github.com/LightChainr/Matching-One/pull/733) | `525c8e98d99c44d4c76280a1f0f83064df9bfa72` | Physical sources, sampler and dictionary correction |
| [#734](https://github.com/LightChainr/Matching-One/pull/734) | `bb41df7da8021d922a48ad353709f35912977c37` | Prior-art sweep, not novelty certification |
| [#735](https://github.com/LightChainr/Matching-One/pull/735) | `9d29d014df28af7c635e6859d98a95ffe2b34d06` | Arbitrary-period root proof and Jordan controls |
| [#736](https://github.com/LightChainr/Matching-One/pull/736) | `64d809b4404f80ff3f9adf9713337cc76008e92d` | New axial full-law iff proof |
| [#737](https://github.com/LightChainr/Matching-One/pull/737) | `7b459c4809baee9dfc2b910091e00a3a8d508d1c` | Spatial-source Hessian and invisible marks |

## Corrections that must travel with the result

| Old inference | Current reading |
|---|---|
| #628 bond duality fails | Code/convention defects, corrected #631/#646/#653 |
| Raw M(p)+M(1-p) diagnoses normalized shape | Wrong centre/gauge; corrected anchored quantiles #702/#706 |
| #675 rules out any single-operator representation | Unrestricted no-go withdrawn; efficiency and existence differ |
| #715 P398 is unnamed | Periodic IC TL via the explicit map #718/#729 |
| #717 finite 2D/0D categorically differs from rank2/rank0 | Same-site event dictionary #733; no all-width intertwiner follows |
| #724/#731 linear split implies semisimplicity | Counterexample #735; m*lambda^m alone is not either diagnostic |
| Primitive Gaussian C3 H8 label identifies local spin | Later unit-rotation/H0 correction in #275 supersedes the near-alias |
| N580 compatibility identifies bare-aspect scaling | Same-block nominal compatibility after #703/#704, not a physical law |

N is site count; square-period linear size is sqrt(N). Do not compare N=425
with a published L=425 lattice. Site/bond, rank-two cross/rank-one spiral,
fixed-p/moving-root and physical/jet operators must remain typed.

Global H4 finite-size evidence, #537's proof obligations, #622's completed #706
analysis, cut networks, finite terminal algebra and publication units remain in
the [atlas](RESEARCH-ATLAS.md). No new Monte Carlo, GPU campaign, width or
angle ladder, generic venue survey, or census expansion is started by exposing
a reserve. Preserve their historical data and failures. Allocation should reduce
open claims, not automatically generate a fresh branch of tasks from every result.
