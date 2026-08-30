# The D x F cross-switch obstruction and its minimal repair

## Exact refutation

The hoped-for universal mark-only injection is false already at `N=6`. The exhaustive minimal gate finds no hard row below `N=6` and four matching-carrier rows at `N=6`, with matrices [[2, 0], [0, 3]], [[2, 1], [0, 3]], [[3, 0], [0, 2]], [[6, 3], [0, 1]]. Each has `(D,M,Y,F)=(12,48,24,24)` and therefore `4DF=1152` hard tokens.

Forcing one flat mark to align with one coexit mark produces no positive target on any of the 4608 hard tokens. Allowing every quotient phase does not cure collisions: each source has exactly one admissible phase, but the unmarked output has only 120 distinct images and maximum fiber 24, so maximum matching is 120 and Hall deficiency is 1032.

The exact lost labels are visible. Decorating by translation phase gives 720 images and leaves fourfold fibers. Decorating also by the source replica gives 1152 singleton images. That is a reconstruction certificate, not an unmarked TM injection: the decorations would add target capacity.

## Obstruction invariant

Let `u_i=1_S(i)+1_T(i)` be the sitewise occupation multiplicity of the two lower configurations after translating the flat face. All balanced configuration crossovers preserve `u`. Exhausting every such crossover still leaves Hall deficiencies [840, 852]. Thus the obstruction is not a bad choice of pairing; the entire union fiber is too small.

Breaking the union invariant in only one carrier is still insufficient: the exact Hall deficiencies are [276, 336]. The signed-union covariance certificate at `a9f7d28` is orthogonal because it concerns local embedding overlaps, not the fixed-line global-rank union fiber.

## Minimal corrected theorem

Allow exactly one occupied-to-vacant replacement independently in each lower configuration, then cross the ordered marks. The resulting two-carrier transport graph has a collision-free matching of `1152/1152` on every one of the four minimal rows. Zero transports, all union-preserving crossovers, and one-carrier transport have already failed, so two carrier transports are minimal within this explicit move hierarchy.

For an arbitrary fixed-line HNF row, Hall saturation of this two-carrier transport graph is sufficient for aggregate TM. The implication is exact; the remaining open step is a topology proof of Hall saturation, naturally organized by the Alexander-dual birth square. The original universal mark-only injection should be considered refuted, not merely unfinished.

## Scientific card

- **Question:** Can `D x F` be removed by a universal configuration-level cross-switch?
- **Answer:** Not by a mark switch or a sitewise-union-preserving base crossover; `N=6` is an exact counterexample.
- **Obstruction:** Translation phase and source replica are collapsed inside a fixed two-base union fiber.
- **Minimal repair:** One occupied-to-vacant transport in both carriers yields full matching on every minimal row.
- **Next theorem:** Prove Hall for that two-carrier transport graph from Alexander complement; do not return to delta-local or one-carrier switching.
