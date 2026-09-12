# Research frontier: three claims, not a catalogue of tasks

**Owner-delegated reset, 2026-09-13.** Baseline read: main at
`eb89e9422791d9e3c3a78f0e65d56912b815a7bd`, including #702--#705.
This document changes research allocation, not historical experimental verdicts.
[ROADMAP](ROADMAP.md) owns the current work choices; live Issue updates override
this dated snapshot. [STATUS](STATUS.md) separates current adjudications from the
[verbatim older ledger](history/STATUS-before-20260913.md).

## 1. Primary: geometry can separate a correct root from an incorrect full law

The most immediately tractable theorem programme is now #613, not another
attempt to assign a continuum field to a fitted exponent.

For independent NN square-site percolation on an honest torus, let
`M=P2-P0`, `F=E[r]/2=(1+M)/2`, and `Q=F^{-1}`. The balance root is Q(1/2).
The conditional rare-sector function `H=P2/(P0+P2)` is a different object.

**Supplied results.** #735 proves balance-root consistency on arbitrary integer
period lattices as the genuine shortest period ell grows, without an aspect or
area-versus-systole condition. Its fixed-subcritical bound is
`P2/P0 <= exp[-kappa(p) N/ell]`. #736 supplies a sharp contrast on axial w-by-m
tori: all fixed interior Q(u) converge to p_c **iff** `log(m)/w -> 0`.
Both statements depend on their named probability/topology inputs. Neither
claims a numerical threshold, an L^-4 shift, or publication priority.

**What was done in this reset.** #736 adds necessity to #613's sufficiency.
Critical square-site RSW and finite-product continuity produce a seam-closed
occupied winding ring with probability at least exp(-eta*w) at a fixed
p_eta<p_c. Independent transverse bands amplify it. If log(m)>=d*w along a
subsequence, every fixed lower quantile is bounded away from p_c, even though
the median can remain consistent. This closes the logarithmic-width case,
not only the earlier fixed-width or w=o(log m) examples.

**Next missing lemma, not another simulation.** Extend the fixed-subcritical
winding-ring construction uniformly around an arbitrary shortest integer
period u, with cost exp(-eta*|u|) and bounded-thickness disjoint corridors.
Oblique orientation and ambient nonprimitivity must be handled on the physical
NN lattice. Changing a period basis does not rotate its interaction. This is
the specific route towards necessity of `log N/ell -> 0` beyond axial tori;
that general necessity is still open in this repository.

A second possible outcome is a geometric obstruction showing that a finer
quantity than the shortest period is necessary. Either would change the paper.
Preparing a paper also requires one independent examination of the supplied
proof and a systematic comparison with strip-percolation/RSW and homological
percolation literature. Repeated negative web searches are not a novelty proof.

## 2. Parallel: the minimal state is determined by the allowed sources

#708 supplies rank-future sufficient lifted states and 509 deterministic
continuation classes at width four. #710 supplies the all-p scalar spectrum:
common modes in P0 and P2 cancel in M. #733 supplies physical site interventions,
source-compatible quotients and an exact final-rank-conditioned backward sampler.
Two independently addressed adjacent columns already require all 509 classes
within the investigated common-strong-lumping class.

The scientific claim is not "509 fields" or "a new transfer-matrix technique".
It is that state compression can be exact for one observer/source language and
fail for another, with an explicit topology-preserving construction and witnesses.
The two-row square-site source is the preferred physical laboratory: its linear
response vanishes while a mixed response is nonzero. P398 is the known periodic
identified-connectivity O(1) TL calibration chain (#718), not the site process.

**Next distinction:** determine the actual all-width topological closure weights
and their observable cancellations, or separate two concrete candidate closures
using the same physical source. A full-operator spectral gap need not be the gap
seen by M. Fixed-width exponential convergence is not a width-uniform theorem.
Reuse the width-four certificate; do not fund a duplicate state enumeration.

### Exact complexity is not finite-noise complexity

#549's cut-network family has k+1 distinct branching-predictive classes within
one identical full unbranched-survival class. In its specified single-fork
readout the adjacent gap is `1/[2k(8k-1)^2]` and the entire family's span is
`1/[2(8k-1)^2]`. Thus this readout becomes less separating as exact class count
grows. This direct consequence of the existing formula is not a lower bound
for every experiment in the full branching language.

A worthwhile reserve question is whether a specified budget of admissible
interventions can amplify those separations. More values of k alone do not
answer it, and an exact discrete class-count theorem is not a Euclidean or
noise-robust dimension theorem.

## 3. High-upside flagship, theory-blocked: original-U identification

#275 is retained as P1 rather than the default P0 expenditure. Its current-asset
closure is UNIDENTIFIABLE_WITH_CURRENT_ASSETS, not "Jordan disproved" and not
"the two physical theories are identical".

The next input remains two candidate-specific maps through the SAME source,
`s0,s1,s2,partial_p s0,partial_p s1,partial_p s2`, physical normalizer, rank-one
denominator and pooled moving-root counterterm. Check each coordinate's baseline
representability, then the nuisance-profiled prediction images on existing data.
No new source or statistical coordinate should be substituted merely because it
separates better without those physical maps.

#735 strengthens the diagnostic warning: positive irreducible stochastic matrices
can share every ordinary trace and thermal derivative while differing in Jordan
structure. For a specified response it is `C N_lambda^k P_lambda B` that determines
visibility. A repeated pole of a derivative-jet lift does not certify a Jordan
block of the physical operator. Spin, a power, and a logarithm are not substitutes
for this source/readout map.

## Pinned research assets and integration boundaries

All entries below were open/unmerged when this reset was prepared. The navigation
commit does not merge them or certify their complete test suites.

| PR | Pinned head | Evidence role |
|---|---|---|
| [#708](https://github.com/LightChainr/Matching-One/pull/708) | `f782061c1a592ed2f9fd0e9dabaa45f0e54bc4e7` | Finite lifted closure and deterministic/scalar rank certificates |
| [#710](https://github.com/LightChainr/Matching-One/pull/710) | `d543ba1afe052a36682d2c1723c8a07281c28950` | Parametric spectrum; geometric minimal-winding configurations |
| [#718](https://github.com/LightChainr/Matching-One/pull/718) | `72d2ee6b8e5d25bd96aaef586c113bc94201f64a` | Axial root consistency; P398 identification with published IC TL |
| [#733](https://github.com/LightChainr/Matching-One/pull/733) | `525c8e98d99c44d4c76280a1f0f83064df9bfa72` | Physical site sources, conditional sampler, event dictionary correction |
| [#734](https://github.com/LightChainr/Matching-One/pull/734) | `bb41df7da8021d922a48ad353709f35912977c37` | Bounded prior-art survey, not an originality certificate |
| [#735](https://github.com/LightChainr/Matching-One/pull/735) | `9d29d014df28af7c635e6859d98a95ffe2b34d06` | Arbitrary-period root proof, marked robustness, Jordan controls |
| [#736](https://github.com/LightChainr/Matching-One/pull/736) | `64d809b4404f80ff3f9adf9713337cc76008e92d` | New axial full-law iff proof and exact finite gluing controls |

#736 reports 21,760 exhaustive strip configurations, 2,365 deterministic uneven-
cell controls, exact Fraction inequalities, three local mathematical tests and
compilation. Full repository CI was not run. The all-size theorem rests on the
written argument and imported inputs, not on those finite checks.

## Corrections that must travel with the result

| Old inference | Current reading |
|---|---|
| #628 bond duality fails | Implementation/convention errors; corrected by #631/#646/#653 |
| Raw M(p)+M(1-p) diagnoses normalized shape | Wrong centre/gauge; use the corrected anchored-quantile quantity (#702/#706) |
| #675 excludes any one-operator representation | Unrestricted no-go withdrawn; distinguish efficient local closure from existence |
| #715 P398 has no named process | Periodic IC TL under the explicit fattening map (#718/#729) |
| #717 finite 2D/0D weights categorically differ from rank2/rank0 | Same finite site ensemble has the event dictionary in #733; no all-width intertwiner follows |
| #724/#731 linear eigenvalue split implies semisimplicity | Explicit counterexample in #735; m*lambda^m alone is not either diagnostic |
| Primitive Gaussian C3's old H8 label identifies local spin | Later unit-rotation/H0 result in #275 supersedes that near-alias |
| N580 compatibility identifies bare-aspect scaling | Same-block nominal compatibility after #703/#704; no physical-law identification |

N is site count. For square-period geometries the linear size is sqrt(N).
Do not compare N=425 with a published L=425 lattice, or treat site, bond,
rank-two cross and rank-one spiral as interchangeable observables.

## What is intentionally not expanded

No default new Monte Carlo, GPU campaign, width ladder, descriptor ladder,
angle grid, generic venue survey or polynomial-height census follows from this
reset. A genuinely discriminating idea may change allocation. The test is the
mathematical distinction it resolves, not its conformity to an ever-growing
checklist. #537's proof obligations, #622's completed shape analysis, the finite
terminal algebra and P2/P3/P4 manuscript work remain visible through the
[atlas](RESEARCH-ATLAS.md); their historical data and failures are retained.
