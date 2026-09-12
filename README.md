# Matching One

**Exact topology, geometric threshold laws, and the information needed to identify a physical response.**

Matching One studies square-lattice site percolation. Its strongest current
questions are no longer "which exponent fits?" or "which familiar field name
looks plausible?" They are:

1. When does a topological balance root locate the infinite-volume threshold,
   even when the whole finite-volume threshold law does not concentrate?
2. Which finite states preserve a declared observable under future updates and
   physical source interventions?
3. Which additional prediction or measurement can actually distinguish the
   remaining mechanisms of the original normalized response U?

**Start with [the research frontier](docs/RESEARCH-FRONTIER.md), then
[the execution roadmap](docs/ROADMAP.md).** This entry point was reset under the
owner's research delegation dated 2026-09-13. Historical data, failed tests of
physical hypotheses, and unmerged branches are retained. Main integration,
mathematical correctness, independent evidence and publication novelty are
separate questions.

## The exact observable

For an occupied configuration on a torus,

    r = rank im[H1(occupied complex) -> H1(torus)] in {0,1,2},
    X = r-1,
    M = E[X] = P2-P0,       E_top = E[X^2] = P2+P0,
    F(p) = E_p[r]/2 = (1+M(p))/2.

The balance root is Q(1/2), where Q=F^{-1}. In contrast,
H=P2/(P0+P2) conditions on the rare non-rank-one sectors; it is NOT F.
The identity X^3=X defines a two-dimensional nonconstant observable algebra,
not a two-state dynamical model or a two-field continuum theory.
These observables sit in the established matching/wrapping/homological
percolation literature; their names alone are not a novelty claim.

## Three research lanes

| Lane | Current asset | Next result worth obtaining |
|---|---|---|
| **Geometric separation of root and full-law consistency** | [#735](https://github.com/LightChainr/Matching-One/pull/735): arbitrary-period balance roots; [#736](https://github.com/LightChainr/Matching-One/pull/736): sharp axial full-law criterion | Check the supplied proof and novelty once; prove or refute the uniform oblique winding-corridor lemma |
| **Observable- and source-dependent minimal state** | [#708](https://github.com/LightChainr/Matching-One/pull/708), [#710](https://github.com/LightChainr/Matching-One/pull/710), [#733](https://github.com/LightChainr/Matching-One/pull/733): exact rank closure, visible spectrum, physical site sources | An all-width structural closure theorem or one source-faithful separating response, not another width count |
| **Original-U physical identifiability** | [#275](https://github.com/LightChainr/Matching-One/issues/275): current assets do not identify the named candidates | Two same-source, same-normalizer forward maps; nuisance-profiled separation before further acquisition |

The linked research PRs are **open/unmerged at this snapshot**; these links do
not silently install their code on main. The [frontier](docs/RESEARCH-FRONTIER.md)
pins their commits and distinguishes results from remaining hypotheses.

For axial w-by-m tori with m>=w and w->infinity, the new organizing contrast is

    balance-root consistency:       no aspect-ratio restriction (#718/#735),
    convergence of every fixed Q(u): log(m)/w -> 0 iff (#736 + #613 inputs).

The iff proof uses critical square-site RSW, finite-product continuity, a
seam-closed crossing ring and independent transverse bands. It is not derived
by fitting sizes. Its arbitrary-oblique necessity extension remains open,
and no near-critical rate or new numerical p_c is claimed.

## What survives from the numerical program

The global orientation-sensitive matching-odd signal and several frozen H4
transfer comparisons remain substantive finite-size evidence. A successful
harmonic comparison is not a standalone nonzero detection and not a field
identification. Single-multiplier full-curve closure and several scalar
correction models have failed; their raw blocks are not discarded.

The primitive-C3 Gaussian observer has a later exact unit-rotation explanation;
its old H8 near-alias must not be exported to square-site U. N580's recovered
covariance and bounded-H8 reanalysis do not establish a bare-aspect physical law.
See [current adjudications](docs/STATUS.md), not the historical opening text of
an issue, for these distinctions. N denotes site count; on square-period
geometries the linear scale is sqrt(N), not N.

## Other valuable assets, without automatic expansion

The cut-network/branching program [#549](https://github.com/LightChainr/Matching-One/pull/549)
separates complete unbranched survival from branching prediction. Its exact
class-count lower bound is not a noise-robust effective-dimension theorem.
Bounded algebraic exclusion and projective inference remain publication units,
not reasons to grow a polynomial census or claim a new statistical principle.
The [research atlas](docs/RESEARCH-ATLAS.md) preserves less prominent work and
closed-unmerged assets; its old execution priorities are historical.

## Reproduction and working rules

[REPRODUCIBILITY.md](REPRODUCIBILITY.md) describes the existing engines and data.
[GOVERNANCE.md](GOVERNANCE.md) keeps exploration lightweight; raw data and
chronology are preserved. [AGENTS.md](AGENTS.md) routes new work to the current
scientific question instead of stale acquisition instructions.

The old README, roadmap, map and status ledger are preserved byte-for-byte in
[docs/history](docs/history/README.md). This is a navigation and allocation
reset, not a claim that older experiments were never run.
