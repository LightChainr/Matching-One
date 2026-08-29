# Exhaustive short-period quotient frontier

Status: finite exact boundary analysis for Issue 269.

The regular square-cell proof assumes that each unit face embeds with four distinct corners. Small
period lattices may identify corners, turn edges into loops, or create repeated local incidences. Rather
than silently importing the surface proof into those cases, this oracle exhausts the first quotient
frontier directly.

## Enumeration

Every index-`n` sublattice of `Z^2` has a two-dimensional Hermite-normal-form representative. The
oracle enumerates every such representative for `2<=n<=7`, constructs the square NN and matching
NN+diagonal quotient graphs, and evaluates all `n!` site permutations. This gives 40 quotient
geometries and 49,878 complete filtrations.

Each quotient is classified by whether the four corners `(0,0),(1,0),(0,1),(1,1)` remain distinct.
Results for self-identifying faces are kept separate from the honest-cell cases.

## Exact gates

For every path the oracle independently checks:

- direct first/second ambient-rank births against historical `K_minus/K_plus`;
- reconstruction of every `R_k` from the two endpoints;
- configurationwise `r_black+r_white=2`;
- both primal/matching reversed endpoint identities with sum `N+1`;
- equality and plateau constancy of the black/white saturated primitive winding line.

Integral saturation indices remain separate from rational lines. Any failure is counted by geometry
and gate, with the first permutation counterexample archived rather than suppressed.

## Boundary

This is an exhaustive finite frontier through index 7 only. A null search does not prove all degenerate
quotients; a discovered counterexample would delimit the regular-cell theorem rather than contradict it.
No production stream or continuum interpretation is changed.
