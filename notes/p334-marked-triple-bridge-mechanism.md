# P334: the missing third-order state is a marked interior bridge

The two real N425 prefixes with identical `(H2,b2,safe-degree-square sum)`
have **one versus six positive interior bridge sites**, carrying respectively
5 versus 19 genuine minimal triple triggers. This locates the first physical
clock difference in actual quotient-site geometry, without a triple census or
another reliability solve.

Source clocks and networks are fixed at `87b6ca5b`: counters 43042508631 (A)
and 43042514803 (B), N425 second orientation, k0=252, age10,
ell=(12,-19). Both have H2=15, b2=12,397 and safe-degree-square sum 3,890,796.
Their genuine triple counts had been recovered algebraically as 5 and 19.
This follow-up reads only those two saved physical networks.

## Exact graph lemma

Start with one ordinary two-terminal factor of the already-certified physical
birth event. Its terminals s,t and occupied-component nodes are always present;
its other vertices are the actual initially vacant sites. Eliminate each
nonterminal fixed node by making its current neighbors a clique, then removing
that node. This preserves s–t connectivity for **every subset of random sites**:
the new edge replaces a path through an always-present vertex, and every
previous path through that vertex is represented by a clique edge. Iteration
therefore preserves the complete monotone event. Keep both terminals.

Let D be random vertices adjacent to both s and t. These are precisely the
original one-site triggers, so no inclusion-minimal trigger of size three can
contain them. Delete D. Partition the remaining sites into

- L: adjacent only to s;
- R: adjacent only to t;
- I: adjacent to neither terminal.

**Lemma.** A three-site set is an inclusion-minimal s–t trigger iff it is
uniquely `{x,y,z}` with x∈L, y∈I, z∈R, edges xy,yz present, and edge xz absent.

Proof: take a simple s–t path using only the three selected sites. Minimality
requires it to use all three, so it is s–x–y–z–t. Its endpoint sites lie in L
and R. If y touches s or t, a two-site subpath already connects the terminals;
thus y∈I. Likewise, edge xz would give a two-site trigger. Conversely, the
stated path connects s,t, while none of its singleton or pair subsets does.
The L/I/R types make its middle site and its two endpoint sites unique. ∎

Consequently the exact marked count is

\[
m_3(y)=|N(y)\cap L|\,|N(y)\cap R|
       -e\bigl(N(y)\cap L,N(y)\cap R\bigr),\qquad
g_3=\sum_{y\in I}m_3(y).
\]

Multiple physical two-port factors are disjoint in their random sites. A
minimal triple must belong wholly to one factor: otherwise any factor that
triggers uses a proper subset. Thus their g3 and marked counts simply add.
Terminal reversal exchanges L/R and leaves m3 unchanged. This lemma uses only
the exact two-terminal representation, not a continuum or exponent assumption.

## Actual middle-site decomposition

Here IDs are the **original vacant-site IDs**, not newly invented features.
For the HNF [[425,268],[0,1]], site v has the exact quotient representative
(v,0). Factor numbering is zero-based in the saved two-port order.

| Prefix | Factor | Middle y | N(y)∩L | N(y)∩R | Product minus existing pair edges | m3(y) |
| --- | ---: | ---: | --- | --- | ---: | ---: |
| A: 43042508631 | 1 | **198** | 40,199,309,355 | 196,354 | 4×2−3 | **5** |
| B: 43042514803 | 1 | 24 | 27 | 71,180,292,338 | 1×4−3 | 1 |
| B | 1 | 25 | 27 | 71,136,180,292,338 | 1×5−3 | 2 |
| B | 1 | **184** | 27,185,341 | 71,180,338 | 3×3−3 | **6** |
| B | 1 | **340** | 27,72,341 | 71,180,338 | 3×3−4 | **5** |
| B | 2 | 94 | 95,362 | 93 | 2×1−0 | 2 |
| B | 2 | 361 | 204,362 | 93,360 | 2×2−1 | 3 |

For A, the three excluded L–R pair edges at y=198 are
(40,196), (40,354), (355,354). The five genuine triples are therefore

```
(199,198,196), (199,198,354), (309,198,196),
(309,198,354), (355,198,196).
```

For B, sites 24,25,184,340 contribute 14 triples in factor1; sites 94,361
contribute another 5 in factor2. The two largest middle sites, 184 and340,
carry 11/19 of B's minimal triples. A's entire triple channel is concentrated
at one site, whereas B has a distributed set of interior bridges in two
different physical connection channels. The numerical split 14+5 does not
identify B's five-triple factor with A's factor; these are different prefixes.

The typed port-address pairs are A factor1=(-109,316), B factor1=(133,558),
and B factor2=(5,430); each separation is425. Synthetic clique edges can pass
through occupied components, so these are paths in the exact occupied-component
reduction, not necessarily three adjacent microscopic sites. The artifact
retains every fixed-node elimination, excluded pair edge, ordered triple, and
expanded contracted-node path. The source maps retain all occupied-component
members and the original prefix/counter for physical replay.

## What this explains

The same scalar state had six minimal-pair edges and five trigger wedges in
both prefixes. It does not constrain which **interior sites simultaneously
bridge still-unconnected L/R pairs**. That additional marked geometric state
is exactly what causes their first difference:

\[
f_3^A-f_3^B=g_3^B-g_3^A=14,
\quad P_B(T=3)-P_A(T=3)=\frac{14}{\binom{173}{3}}
=\frac7{424023}.
\]

For a uniform first-three-site subset, `m3(y)/C(173,3)` is its exact genuine
triple contribution to third-step birth, resolved by **structural middle** y.
This mark must not be confused with the actual last insertion site: within
one minimal triple any of its three sites can arrive last. Nor does m3 alone
explain the complete later clock, whose higher-order terms can differ.

This yields a concrete candidate marked-state observable for subsequent
research: the interior-bridge profile `y → m3(y)`, retaining its factor and
port type. It is an exact, named finite-lattice statistic beyond pair counts
and pair-degree second moments, not a CFT field identification or path-memory
claim.

## Artifact

`scripts/p334_marked_triple_bridges.py` reads the two existing maps and writes
`results/p334-marked-triple-bridges/middle_site_bridges.json`. It took0.11s in
the existing local research environment. No global C(173,3) traversal, new
Monte Carlo, network DP, remote connection, or full-repository test suite was
run. This is a posthoc mechanism readout of the same source block as the147
clock analysis, not new independent production evidence.
