# Translation-orbit Hall and the missing transverse reservoir

## Semantic correction

Every output face must have a lower base that is itself rank one on the same fixed projective line. Enforcing this omitted gate changes the N=6 two-carrier result from the invalid `1152/1152` to `588/1152`. The mark-only refutation remains valid; the claimed base-only repair is withdrawn.

## Exact orbit theorem for arbitrary HNF

Simultaneous translation acts freely on every ordered face token: a stabilizer must fix its first marked site, hence is the identity. Every source and cover orbit therefore has size `N`, independently of whether the Smith group is cyclic. Translating the first marked site to zero gives a unique orbit representative.

A matching in the normalized orbit graph lifts uniformly to a fractional matching of the full graph. Bipartite matching integrality then produces a collision-free full matching. Conversely a full matching averages to orbit flow. Thus full Hall and orbit Hall are exactly equivalent, with an exact factor-`N` compression.

## Minimal obstruction and exact extra reservoir

The exhaustive gate has no hard row below `N=6` and 4 minimal rows at `N=6`. On each, base-only two-carrier transport matches `588/1152` and one transverse-mark release with fixed bases matches `768/1152`. Neither resource is sufficient.

Their conjunction is sufficient: exchange one occupied base site for one vacant site in either carrier, and independently release exactly one of the four crossed output marks to a transverse quotient site. The compatibility graph matches `1152/1152` on every minimal row. This is the smallest successful move in the tested two-axis lattice (base transport, transverse release).

The interpretation is precise: the base exchange carries the translated Alexander birth-square configuration; the released mark supplies the transverse line absent from the original D/F four-mark data.

## First Smith gate

At `N=8`, matrix `[[2, 0], [0, 4]]` with Smith invariants `[2, 4]`, matching carrier, and middle layer `k=4`, base-only orbit matching is `2496/5760`. The corrected reservoir saturates `5760/5760`. The raw source set has 46,080 tokens; orbit normalization reduces it exactly to 5,760.

The remaining general statement is now sharp: prove orbit-Hall saturation of this corrected graph from digital Alexander complement for every HNF, or locate its next minimal obstruction. There is no reason to revisit mark-only or base-only switching.

## Scientific card

- **Correction:** fixed-line output-base semantics invalidates the earlier base-only repair.
- **Obstruction:** N=6 base-only deficiency `564`; transverse-only deficiency `384`.
- **Reservoir:** one carrier base exchange plus one free transverse output mark.
- **Certificate:** `1152/1152` on all minimal rows and `5760/5760` on the first Smith-(2,4) middle-layer gate.
- **General step:** orbit Hall is exactly equivalent to full Hall for arbitrary HNF; only uniform Alexander saturation remains open.
