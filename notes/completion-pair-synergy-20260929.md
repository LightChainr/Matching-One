# Completion-pair synergy and the successor law of completion count

2026-09-29. Bounded author derivation for the finite uniform-insertion model.
This continues the [hazard identity](completion-pivotal-hazard-kernel-20260929.md)
and [directional closure note](directional-markov-completion-20260929.md).
It gives no literature-novelty assessment, empirical validation or scaling claim.
The next discriminating quantity is creation of completion sites after a safe
insertion: its distribution, not just the already fixed immediate exit hazard.

## 1. The completion graph gives the exact next-state row

Let V contain N sites and let `r:2^V -> {0,1,2}` be monotone, with
`r(empty)=0`, `r(V)=2`. Fix a rank-one prefix A of size k and put `m=N-k`.
Conditional on the revealed prefix, remaining sites arrive in uniform order.
Write `A+v=A union {v}` and define

    C(A) = {v in V\A : r(A+v)=2},       c=nu_2(A)=|C(A)|,
    S(A) = (V\A)\C(A),                 s=|S(A)|=m-c.

Define the simple undirected graph G_A on S(A) by

    {v,w} in E(G_A) iff v != w and r(A+v+w)=2.

Its edges are pairs that complete together although neither completes alone;
they need not be adjacent sites of the physical lattice. Put `e=|E(G_A)|`.
For a safe v, every old completion site remains vacant and completing by
monotonicity. Among the other safe sites, exactly its graph neighbours become
completion sites. Thus the stronger set identity, and its count version, are

    C(A+v)=C(A) union N_G(v),           nu_2(A+v)=c+deg_G(v).       (1)

Let `a_j(A)=#{v in S(A):deg_G(v)=j}`. Then `sum_j a_j=s` and
`sum_j j*a_j=2e`. Assume D is a persistent rank-one mark: inclusion between
rank-one sets preserves D. For ambient homology its unoriented rational line
is persistent because homology images are nested; arbitrary marks need not be.
Observe `X_k=0`, `(1,D(A_k),nu_2(A_k))`, or `2`, according to rank.
For `m>0`, the complete next-state row from this fixed A is

    P(X_(k+1)=2 | A)=c/m,
    P(X_(k+1)=(1,D(A),c+j) | A)=a_j(A)/m.                        (2)

If `s>0`, survival has probability s/m and the conditional completion-count
law is `P(nu_2(A_(k+1))=c+j | A,survival)=a_j/s`. Direction is unchanged.
If `s=0`, exit is certain and this conditional law is undefined, not zero.

## 2. Two drifts and an exact two-insertion survival formula

Extend nu_2 by zero after rank two, solely as a killed observable. Starting
from A let `Z'=1{r(A_(k+1))=1}nu_2(A_(k+1))`. Summing (1) over safe sites,

    E[Z' | A]=(s*c+2e)/m,
    E[Z'-c | A]=(2e-c^2)/m.                                     (3)

This killed drift includes a loss of c on each completing insertion, whose
probability is c/m. It can be negative. In contrast, the ordinary count
increment conditioned on remaining rank one satisfies, for `s>0`,

    E[nu_2(A_(k+1))-c | A,survival]=2e/s >= 0.                   (4)

For `m>=2` and `s>0`, the conditional next exit probability after that safe step is
`(c+2e/s)/(m-1)`. Its difference from c/m is
`c/[m(m-1)]+2e/[s(m-1)]`: removal of a vacancy and creation of completion
sites both matter. This conditional increase is not the killed drift (3).

For two insertions, survival means choosing an ordered pair of safe sites
that is not an edge. There are `s(s-1)-2e` such pairs out of `m(m-1)`:

    P(J2>k+2 | A)=[s(s-1)-2e]/[m(m-1)],                         (5)
    P(J2=k+2 | A)=[s*c+2e]/[m(m-1)].                            (6)

Here `J2=min{j:r(A_j)=2}`. Equation (6) concerns first completion at the
second insertion; first-insertion completion has probability c/m. These
three probabilities sum to one. Equivalently, surviving the first step
and then using its mean completion count gives (5), with no Markov assumption.
The undivided formulas remain valid when s=0 or s=1: two-step survival is
zero in either case. For s=0 there is no surviving first-step cohort; for
s=1, e=0 and its sole safe insertion makes every remaining site completing.
When m=1 only the one-step formulas apply; a second insertion does not exist.

## 3. Exact contrast between matched birth cohorts

At the same count k, match the exact direction d and integer count c. For
a fixed cutoff `a<k`, let E and L be the positive-probability cells

    E={J1<=a, X_k=(1,d,c)},       L={a<J1<=k, X_k=(1,d,c)},
    J1=min{j:r(A_j)>=1},         Delta_e=E[e(A_k)|E]-E[e(A_k)|L].

The future remains uniform conditional on each microscopic prefix, so (5)
can be averaged directly. For `m>=2`, the early-minus-late contrast is

    P(J2>k+2 | E)-P(J2>k+2 | L)=-2*Delta_e/[m(m-1)].             (7)

Both first-step survivals are exactly s/m. If `s>0`, this probability is
constant across configurations in the cell, so conditioning on first-step
survival does not reweight its distribution of A_k. Consequently

    E[nu_(k+1)|E,J2>k+1]-E[nu_(k+1)|L,J2>k+1]=2*Delta_e/s,
    P(J2>k+2|E,J2>k+1)-P(J2>k+2|L,J2>k+1)
        =-2*Delta_e/[s(m-1)].                                  (8)

Here nu_(k+1) denotes nu_2(A_(k+1)); the second line requires m>=2.
The killed successor-mean difference is `2*Delta_e/m`, by (3).
For s=0, e=0 in both cohorts and (7) is zero, while (8) is undefined.
Absent cohort cells have no conditional contrast. Exact matching matters:
mixing different c also mixes the baseline `s(s-1)` and survival weights.
For common weights w_(d,c) summing to one at fixed k, the standardized
contrast is exactly `-2 sum_(d,c) w_(d,c) Delta_e(d,c)/[m(m-1)]`.
Cell effects can cancel; a zero weighted contrast need not mean zero per cell.

For `s>0`, the full next-count contrast, conditional on survival, is
`[E[a_j|E]-E[a_j|L]]/s` at value c+j. Equal mean e tests only one moment
of this law. Equation (7) supplies an exact two-step consumer for that moment.

## 4. What closure would actually require

**Strong rank-one closure** means one transition kernel works for every
microscopic preparation having the same `(k,D,c)`. By (2), this holds iff
the entire vector `(a_j(A))_j` is identical across every such rank-one cell.
Necessity follows by preparing each A separately; sufficiency follows by
averaging the common row. Equal e alone gives only a common two-step survival.
Different degree vectors refute this strong claim, even if their means agree.

For a **full strong closure** including arbitrary rank-zero preparations,
one must additionally have equal rank-zero counts of insertion sites leading
to each next observed state, including each `(1,d,c')` and direct rank two.
The rank-one graph criterion says nothing about those entry rows. Rank two
is absorbing. Direct `0 -> 2` jumps have no intervening rank-one direction.

For **weak closure from the uniform empty start**, let
`F_k^X=sigma(X_0,...,X_k)` under that particular source law. The exact
rank-one condition is, for every k and every j on positive-risk cells,

    E[a_j(A_k) | F_k^X]=E[a_j(A_k) | X_k].                       (9)

The fixed external k is implicit on the right. The exit row c/m is already
fixed by X_k; (9) supplies all the other probabilities in (2). On rank zero
the observed past is all zero, so its averaged transition row automatically
depends only on k and current state. Uniform empty preparation makes A_k
uniform over rank-zero sets of size k there. With absorption in rank two,
(9) is necessary and sufficient for the full weak Markov property. In
particular, strong rank-one rows suffice for this weak full process without
strong rank-zero lumpability, as in the [L3 construction](birth-completion-geometry-20260929.md).

The observed past now contains **the entire nu history**, as well as J1 and
D. Unlike the direction-only observer, that past is not encoded by entry
time alone. A population discrepancy in the matched birth-cohort successor
law, or in (7), refutes weak Markovness because birth cohort is observed past.
Equality for one split does not prove it; even equality of full successor
laws for all birth cutoffs need not establish (9) after conditioning on
different nu histories. Conversely, distinct microscopic degree vectors
need not refute weak closure: the source-conditioned history mixtures may
still satisfy (9). Reachability of two configurations alone is insufficient.

The parent's separate [exact L4 companion](completion-pair-synergy-L4-20260929.md)
reports matched birth-cohort discrepancies for square and triangular models.
These establish actual weak failure from the uniform empty start, as well as
strong failure; the finite witnesses and census evidence belong to that note.

## 5. Why the present graph is not automatically a recursive state

After a safe insertion v, (1) gives

    S(A+v)=S(A)\({v} union N_G(v)).

Edges of G_A induced on this new safe set persist by monotonicity. Additional
edges `{w,z}` can appear if `r(A+v+w+z)=2` although `r(A+w+z)=1`.
Thus graph incidence can matter when vertices and their neighbours are
removed, and triple-completion information can matter for newly created
edges. A degree histogram does not specify the former; even the whole
current graph does not in general specify the latter under monotonicity alone.

For example, three currently safe sites may have no completing pair, while
their triple completes; after the first safe insertion the other two form
an edge. Monotonicity also permits their triple to remain rank one with the
same initial pair data. These are abstract local possibilities, not asserted
lattice witnesses. They explain why even adding e or G_A is no automatic
closure proof. The bounded question here is the next-count law (2) and its
two-step implication (7), not an unending hierarchy of proposed descriptors.

## 6. A physical prediction and the next scientific action

For occupied-graph homology, a safe insertion can join components or extend
a path without adding a second homology direction. Its graph neighbours are
then newly completing vacancies: together with that insertion they supply
the second direction. The increment `nu_2(A+v)-c` counts these newly enabled
opportunities. This gives a concrete interpretation of pair synergy as
creation of future completion sites while rank and D remain unchanged.

A prospective comparison can fix the lattice, count k, birth cutoff and
common-support weights, then retain independent naturally sampled rank-one
prefixes matched on exact `(D,c)`. In each prefix with s>0, choose a uniformly
random safe site, insert it in a copy, and measure `Y=nu_2(A+v)-c`.
Conditional on that prefix, `E[Y]=2e/s`; at matched c, this samples the actual
next-step law conditioned on survival. Within each cell, the early-minus-late
two-step survival contrast is exactly `-s/[m(m-1)]` times its mean Y contrast.
Forced safe insertion is a conditional continuation experiment, not the
unconditional process; its survival factor s/m restores the latter prediction.

Multiple continuations of one prefix share its geometry and birth history.
They can reduce continuation noise but are not independent configurations;
average within prefix and retain prefix/batch dependence in cohort uncertainty.
The full Y distribution can address (2), while its mean has the fixed consumer
(7). The reported L4 failure already redirects finite-model work from accepting
`(D,nu)` closure to how safe insertions create completion sites. A prospective
comparison would address that specific mechanism at its declared target size.
Equality leaves recursive sufficiency open and does not justify a feature scan.

The parent's L4 census is not recomputed here. Existing L512 evidence is
reported only in its [source note](completion-hazard-independent-block-20260929.md)
and [retrospective companion](directional-hazard-contrast-20260929.md); this
derivation neither reanalyses nor upgrades it. All formulas here use the count
clock; hidden-count label-time closure requires its own argument. No tests,
enumeration, cloud work, web research, or data analysis were performed here.
