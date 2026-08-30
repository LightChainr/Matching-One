# Exhaustive short-period quotient frontier

Status: finite exact boundary analysis for Issue 269.

The regular square-cell proof assumes that each unit face embeds with four distinct corners. Small
period lattices may identify corners, turn edges into loops, or create repeated local incidences. Rather
than silently importing the surface proof into those cases, this oracle exhausts the first quotient
frontier directly.

## Enumeration

Every index-`n` sublattice of `Z^2` has a two-dimensional Hermite-normal-form representative. The
oracle enumerates every such representative for `2<=n<=12` and constructs the square NN and matching
NN+diagonal quotient graphs. This gives 126 quotient geometries and 13,961,736,918 complete filtrations.

For each geometry the oracle first computes the exact primal and complementary matching marks for
all `2^n` occupied-site subsets. Every site permutation is a maximal chain in this Boolean subset
lattice. Exact dynamic programming counts those chains, while a prefix set `S` receives the exact
weight `|S|!(n-|S|)!` when plateau-step statistics are accumulated. A direct cached-versus-uncached
test exhausts all paths on the axis `L=2` control, and an independent regression compares the new
subset-lattice summary field-for-field with explicit permutation enumeration on every HNF quotient
through index 6. The acceleration therefore preserves the existing path-level observable semantics.

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

This is an exhaustive finite frontier through index 12 only. A null search does not prove all degenerate
quotients; a discovered counterexample would delimit the regular-cell theorem rather than contradict it.
No production stream or continuum interpretation is changed.
