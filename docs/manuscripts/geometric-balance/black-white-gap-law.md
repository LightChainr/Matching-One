# Rare black barriers and giant white components on a cylinder

2026-09-14. A continuation of the mathematical question in #764/#739.

The point of this note is a connection, not a publication plan: a rare, short
component of one colour determines the longitudinal scale of a very large
component of the other colour. All colour probabilities below use the SAME
Bernoulli labels. No independent recolouring is introduced.

## 1. Objects and results

Use the cylinder `(Z/wZ) x Z`, with its physical horizontal edge gains. Black
sites have probability p and NN connectivity; white sites have probability
q=1-p and matching (NN+NNN) connectivity. Take w>=3. The width-two controls
retain the two lifted parallel horizontal edges explicitly.

An essential component contains a cycle of nonzero horizontal winding. Its
complete span is `L=max y-min y+1`. A component is counted ONCE, with anchor
its lowest occupied row (the least column in that row is only a tie-breaking
vertex label). Write nu_B,nu_W for the numbers of essential components per
vertical row, and E_B,E_W for their component-Palm averages. Neither average
is a random-vertex average. At any fixed w and 0<p<1 all components of both
colours are finite: full white rows cut black connectivity, full black rows
cut white connectivity, and both occur infinitely often in both directions.

### Exact finite-width conclusions

The complete essential components form a bi-infinite alternating chain

    ... < B_i < W_i < B_(i+1) < W_(i+1) < ... .               (1.1)

Their lowest rows and their highest rows are strictly increasing in this
order. In particular there is at most one essential anchor of either colour
in a given row. The translation-covariant pairing gives

    nu_B=nu_W=:nu.                                          (1.2)

Let a_i=min y(B_i), D_i=a_(i+1)-a_i, and let L_i be B_i's span.
For the intervening WHITE complete component,

    D_i-L_i <= L(W_i) <= D_i+L_(i+1).                        (1.3)

The stationary gap identity `nu E_B D_i=1` then gives

    |nu E_W L - 1| <= nu E_B L.                             (1.4)

These statements do not require independence of gaps, a dilute parameter,
a correlation-length prefactor, or an asymptotic limit. They are exact for
the matching pair above.

### Fixed-subcritical, increasing-width conclusions

For each FIXED `0<p<pc_NN`, put `nu=nu_w(p)`. Then the following hold as
w->infinity:

    nu E_B L^k = O(nu w^k) -> 0,             every fixed k;
    nu L_white => Exp(1),  with all fixed positive moments;
    nu^k E_W L^k -> Gamma(k+1),              k>0;
    nu E_W L -> 1,   Var_W(L)/(E_W L)^2 -> 1.                (1.5)

Here `Gamma(k+1)` is a moment, not a Gamma-distributed component length.
The white parameter is q=1-p, on the supercritical side of the plane
matching graph. A Brownian-span assertion for a FIXED SUBCRITICAL white
parameter would address a different experiment.

Under a uniformly chosen longitudinal position, the containing interval,
and asymptotically the unique white essential component spanning that row,
has instead

    nu L_seen_from_row => Gamma(shape=2,rate=1).             (1.6)

This is length bias. It is not a second physical phase transition.

The full black and white ANCHOR point processes, scaled vertically by nu,
converge jointly to the SAME rate-one Poisson point process. They are not two
independent Poisson processes. Consecutive white spans under component Palm
converge jointly to independent exponential spacings; this independence is
an asymptotic consequence, not a finite-width assumption.

The proofs below use only elementary cylinder topology, the plane site
volume tail [AV], the exploration coupling, and local-dependence Poisson
approximation [AGG]. They do NOT use the earlier mass-inversion theorem,
Gumbel affine-window theorem, p-differentiability, OZ sewing, or a unit residue.

## 2. The alternating chain is a topological fact

Use the usual complementary regular-neighbourhood representation of the
4/8 digital matching pair. One construction thickens occupied NN vertices
and edges and fills elementary all-black faces. Complementary regions have
the connectivity and essential topology of the white matching graph. It is
enough to check the sixteen face patterns. Crossing white diagonals in an
all-white face are connected, not two crossing distinct components. This is
the same local boundary-connectivity fact underlying [MZ, Section II].

A compact connected essential neighbourhood on the open annulus has exactly
two essential boundary curves. Additional boundary curves enclose discs.
The region immediately on its upper boundary belongs to ONE opposite-colour
component, because the boundary is connected and the local matching rule
connects its surrounding opposite-colour sites. This opposite component is
itself essential. The same argument applies below. Two disjoint connected
essential components are ordered by which end of the annulus lies on their
side of an essential separating curve. A contractible component cannot be
inserted into this order: it does not separate the two ends.

Starting from any essential component and following its upper neighbour
therefore alternates the colours. There is no omitted essential component
between neighbours, since they meet along the same separating boundary.
The full monochromatic rows ensure that the chain extends to both ends and
exhausts all essential components. Local finiteness excludes accumulation
inside a bounded set of rows. This proves (1.1).

### Height order and contact bounds

If D is above C, choose a simple essential curve inside C. A lattice vertex
of D at or below the lowest row of C could be joined downwards to the bottom
end without intersecting this curve: below that row there are no vertices
or edge segments of C, and the proposed vertical ray at the vertex's column
cannot pass through C's distinct lattice vertex. Such a point belongs to
the lower, not upper, region. Hence `min y(D)>min y(C)`. The analogous upward
ray gives `max y(D)>max y(C)`. The argument works also for the planarized
matching curve: its added face points stay between their endpoint heights.

Neighbouring opposite-colour components have boundary-adjacent sites whose
vertical coordinates differ by at most one. Thus, for B_i,W_i,B_(i+1),

    a_i < min W_i < a_(i+1),
    min W_i <= max B_i+1 = a_i+L_i,
    a_(i+1)-1 <= max W_i < max B_(i+1).

Writing `L(W_i)=max W_i-min W_i+1` gives (1.3). Branches, holes and multiple
essential cycles in a component do not change this argument. The height is
the complete component's height, not a chosen loop's height.

There is another useful version. Let N_B(y),N_W(y) count essential components
whose vertical projections include row y. Their active indices form a
contiguous segment of the alternating chain, because both endpoints of the
spans are ordered. This segment is nonempty: neighbouring spans have no
integer-row gap by the contact bound. Consequently

    N_B(y)+N_W(y)>=1,       |N_B(y)-N_W(y)|<=1.               (2.1)

In particular, if no black essential component spans y, exactly one white
essential component does. Campbell counting gives `E N_B(0)=nu E_B L`.

## 3. Stationarity, gaps and the mean law

Pair B_i with its upper neighbour W_i. This pairing is bijective and
commutes with vertical translations. Counting matched components in a long
interval, or applying the mass-transport identity, proves (1.2). No renewal
property is needed.

For the black anchor process, assign integer rows
`a_i,...,a_(i+1)-1` to anchor a_i. They partition Z. Therefore

    nu E_B D_i=1.                                          (3.1)

Moving to the next black anchor preserves the component-Palm law, so
`E_B L_(i+1)=E_B L_i`. The bijection also transports the W_i mark to white
component Palm without size bias. Averaging (1.3) proves (1.4).

Equivalently, from (2.1),

    1-nu E_B L <= nu E_W L <= 1+nu E_B L.                    (3.2)

Thus the entire mean relation is topological plus stationary. Poisson
approximation is only needed to determine the DISTRIBUTION around this mean.

## 4. Uniform black volume control and rarity

[AV, Theorems 2--3] gives, for fixed p<pc_NN, constants C,c>0 such that the
plane occupied root cluster has `Pr(|C(0)|>=n)<=C exp(-cn)`.

To apply it uniformly in the cylinder width, use the following coupling.
Explore a cylinder cluster using a fixed queue. At the first query of a
quotient vertex, choose a neighbouring lift of its already assigned parent
and read that as-yet unqueried plane label. Distinct quotient vertices have
distinct lifts. A repeated quotient vertex uses its old answer and does not
query a new copy. The answers thus have exactly the independent cylinder
law, and the discovered occupied TREE injects into the independent plane
root cluster. Cycles need not lift as cycles; size, which is all this step
uses, is preserved. Hence the same C,c give the cylinder volume bound.

A winding component contains at least w distinct vertices. An isolated full
black row, with white rows immediately above and below, is one complete
essential component. Thus

    [p(1-p)^2]^w <= nu <= C w exp(-cw).                      (4.1)

For component Palm, selecting its anchored lowest-row site gives

    nu Pr_B(K>=n) <= C w exp(-cn).                          (4.2)

Split the moment integral at a sufficiently large constant times w and use
the lower bound in (4.1). For every fixed r>0 this proves

    E_B K^r=O(w^r),      E_B L^r=O(w^r).                    (4.3)

In particular epsilon_w=nu E_B L tends to zero exponentially up to a
polynomial. The mean conclusion `nu E_W L->1` already follows, without any
point-process limit. Equation (2.1) also proves that the probability a fixed
row is not spanned by exactly one white essential component tends to zero.
This is row COVERAGE, not probability that a uniformly chosen site belongs
to the white component.

## 5. A self-contained rare-anchor Poisson argument

This section shortens the needed probabilistic input. It does not assume
that the cylinder intensity's exponent has already been identified with the
plane inverse correlation length.

Let H=w^2. Let I_y^H indicate a COMPLETE black essential component with
lowest row y and span at most H. There is at most one such component.
The event can be decided from the rows y-1,...,y+H, including the rows that
certify no connection outside. Its probability is nu_H and

    0<=nu-nu_H<=Cw exp(-cH)=:eta_H.                          (5.1)

Consider n<=T/nu anchor rows, T fixed. Indicators at distances greater than
H+1 have disjoint supports. Thus, in the notation of [AGG], b3=0 and

    b1<=n(2H+3)nu^2.                                       (5.2)

Here is the essential b2 step. Two distinct complete anchored components
have vertex-disjoint occupied winding witnesses. Each witness is contained
in a common slab of at most 3H+4 rows. If E is the increasing event of ANY
occupied essential cycle in this slab, then

    {I_y^H=I_z^H=1} subset E disjoint-occurrence E.

We use site BK on E, NOT on the nonmonotone complete-anchor events.
An occurrence of E either belongs to a full component of span <=H, whose
anchor lies in an enlarged interval of <=4H+4 rows, or a vertex in the slab
belongs to a component of size >H. The intensity and volume bounds give

    Pr(E)<= (4H+4)nu + Cw(3H+4) exp(-cH)=:a_H.

Hence

    b2<=n(2H+3)a_H^2.                                      (5.3)

Equations (4.1),(5.1)--(5.3) give b1+b2=o(1): the leading term is
`O_T(H^3 nu)`, and every remaining term contains exp(-c w^2) divided by at
most an exponential in w. The probability of omitting any true anchor in
the observation interval is at most n eta_H=o(1).

[AGG, Theorem 2], with our total variation convention `sup_A|P(A)-Q(A)|`,
then gives the joint count approximation on any finite partition. Their
printed norm is twice this convention; the bound we use is 2(b1+b2) when
b3=0. On rescaling row coordinates by nu, the independent Poisson intensity
measure converges to Lebesgue measure. The black anchors therefore converge
to a stationary rate-one Poisson point process on R.

### Uniform gap tails, not just bounded-window convergence

Choose a block of length `ceil(1/nu)+2H+2`. Restrict local anchors to the
interior so their H+2 supporting rows are wholly inside the block. By the
just proved Poisson limit, the probability of at least one such anchor is
bounded below by a fixed a>0 for all sufficiently large w. These events on
disjoint blocks are independent. Therefore the stationary void probability
V_w(n) satisfies, with constants independent of sufficiently large w,

    V_w(n)<=C0 exp(-c0 nu n).                               (5.4)

The event used here is a complete local anchor, not just a path belonging
to a possible component rooted outside the block.

## 6. From voids to component-Palm spacings

For a stationary integer anchor process with gaps D, direct counting yields

    V_w(n):=Pr(no anchor in rows 1,...,n)
           =nu E_B[(D-n)_+].                              (6.1)

Each gap of length D contains exactly `(D-n)_+` starting positions with no
anchor among the next n rows. This identity is valid without independent
gaps. At n=0 it includes (3.1).

Linearly interpolate (6.1) and put X_w=nu D. The Poisson limit gives

    E[(X_w-t)_+] -> exp(-t),       t>=0.                    (6.2)

These are convex functions of t, with a differentiable limit. Their
one-sided derivatives bracket `-Pr(X_w>t)`. Difference quotients and then
a shrinking increment imply

    X_w => Exp(1).

There is no unjustified differentiation of an arbitrary asymptotic error.
For t>=1, `Pr(X_w>t)<=E[(X_w-(t-1))_+]`; (5.4) gives a uniform exponential
tail. All fixed positive moments of X_w therefore converge.

The same point-process limit also gives finite vectors of successive Palm
gaps: the limit is a vector of independent Exp(1) variables. One way to
justify the Palm step explicitly is Campbell's formula on a bounded
rescaled interval. The summand is bounded by its anchor count. The latter
has uniformly bounded second moment: distant truncated indicators are
independent, nearby pair contributions are bounded by (5.3), and the full
versus truncated second-moment discrepancy is bounded by a polynomial in
n times exp(-cH). This supplies uniform integrability. Local continuous
Palm tests pass to the Poisson Palm law; truncate very long gaps using
(5.4) before removing the truncation.

## 7. The white law, its moments, and the coupled point processes

Under the black-to-white pairing, (1.3) implies

    |L(W_i)-D_i| <= max(L_i,L_(i+1)).                       (7.1)

Together with (4.3), nu times this difference tends to zero in every fixed
L^r space. Equations (6.2),(7.1) prove (1.5), including all fixed moments and
finite vectors of consecutive white spans. In particular the limiting
white span coefficient of variation squared is ONE, not the Brownian-range
constant pi/3-1 applicable to a different, subcritical component law.

For the lower anchor c_i of W_i the sharper geometric relation is

    1<=c_i-a_i<=L_i.                                       (7.2)

On a bounded rescaled interval, the expected number of paired anchors whose
rescaled displacement exceeds epsilon is bounded by a constant times
`Pr_B(nu L_i>epsilon)`, plus negligible endpoint terms. This tends to zero.
Thus the white anchor process and black anchor process collapse onto the
SAME Poisson process. For disjoint bounded intervals J_l,

    (N_B(J_l/nu),N_W(J_l/nu))_l => (Z_l,Z_l)_l,
    Z_l independently Poisson(|J_l|).                     (7.3)

This does not violate the torus prohibition on having no wrapping component
of either colour. Here N counts ANCHORS IN A WINDOW. Both counts can vanish
while one white component anchored outside spans the entire window.

### A row-selected component is length biased

Let D^dagger be the black gap containing a uniformly chosen row. Exactly,

    E f(nu D^dagger)=nu E_B[D f(nu D)].                      (7.4)

Consequently `nu D^dagger=>Gamma(2,1)`, with all fixed moments. Its mean
is asymptotically 2, whereas the component-Palm gap mean is 1.

There is also an exact incidence-Palm identity for white spans: selecting
one pair (row, white component whose projection includes that row) gives

    E_inc f(nu L)= E_W[L f(nu L)]/E_W L.                    (7.5)

It has the same Gamma(2,1) limit, by the already proved white moments.
The symmetric difference of the paired intervals
`[a_i,a_(i+1)-1]` and `[min W_i,max W_i]` has length at most
`L_i+L_(i+1)`. The proportion of rows affected is at most 2 epsilon_w.
Moreover (2.1) says there is a unique white projected component with
probability at least 1-epsilon_w. Thus ordinary row observation and (7.5)
agree asymptotically. Higher moment convergence follows by the incidence
moment formula and uniform integrability at one higher order.

## 8. Closing the long cylinder into a torus changes the zero-count cell

Take an axial torus with `m nu_w(p)->t in (0,infinity)`, p fixed subcritical.
The same truncated-anchor proof works with periodic dependency neighbourhoods.
The chance of any component with size >H is at most `wm C exp(-cH)=o(1)`.
Hence vertical wrapping is absent with probability tending to one, and the
number K of black horizontal components tends to Poisson(t).

The finite matching topology [MZ] then determines the WHITE global count:

    (W_black,W_white) => (K, K+1{K=0}).                    (8.1)

If K=0, there is one white rank-two component. If K>=1, there are K white
rank-one components. Consequently, for 0<=z,y<=1,

    E[z^W_black y^W_white]
      -> exp(t(zy-1))+exp(-t)(y-1).                        (8.2)

The limit is NOT independent Poisson and is NOT simply two identical
Poisson counts. The one exceptional zero-count term is the global closure.
The black rank probabilities tend to `(exp(-t),1-exp(-t),0)`.

Conditionally on K=k>=1, the cyclic macroscopic gaps between anchors are
times t a Dirichlet(1,...,1) spacing vector, after a uniformly selected anchor
fixes the cyclic origin. White components occupy those gaps up to vanishing
scaled boundary spans. At k=0 the white cross is a different topology, not
an artificial additional anchor from an independent process.

## 9. Exact all-height controls and a new bulk-volume question

The attached script reuses UNMODIFIED tagged_winding_span.py, Git blob
`52f3611990ce2b1331d9e5296e0262f5e402e0d7`. It computes both colours at the
complementary parameters p=1/4 and q=3/4, using exact rational arithmetic.
It does not use a cutoff histogram to estimate a long white mean.

| w | nu | E_B L | E_W L | nu E_W L | CV_W^2 |
|---|---:|---:|---:|---:|---:|
| 2 | .05679086538 | 2.07733008 | 16.46886447 | .9352810651 | .8865951229 |
| 4 | .004429300955 | 3.13710294 | 224.88167088 | .9960685995 | .9846477960 |
| 6 | .0004219779152 | 3.93801838 | 2369.27391382 | .9997812668 | .9978072200 |

The JSON stores fractions, not only these rounded displays. The errors in
`nu E_W L-1` satisfy (1.4) exactly. It also evaluates the rational transform
`E[(1+s nu)^(-L)]`; its limiting target is 1/(1+s), without rounding exp(-s nu).

### Exact volume identities

Let K be the number of sites in a COMPLETE white essential component. Define

    theta_w=Pr(a uniform site belongs to an essential white component),
    chi_w=E[|C(0)| 1{C(0) is white and essential}].

Two direct mass transports give

    theta_w=(nu/w) E_W K,
    chi_w=(nu/w) E_W K^2.                                  (9.1)

Thus

    nu chi_w/[2 w theta_w^2] = (1+CV_W(K)^2)/2.              (9.2)

The independent direct activity transfer verifies these finite identities.
For w=2,3,4 at q=3/4 it gives `Corr_W(L,K)^2` approximately
`.99362833, .99849621, .99964572`. These are finite exact-rational results,
not proof of a large-width correlation limit.

### Bulk-filling conjecture (NOT proved here)

Write theta_8(q)=Pr_plane(0 lies in the infinite matching cluster). The new
question is whether a giant white slab has asymptotically deterministic bulk
density theta_8(q):

    (nu L, nu K/w) => (E,theta_8(q) E),    E~Exp(1),         (9.3)

with at least second moments. If true, it gives

    theta_w -> theta_8(q),
    chi_w ~ 2 w theta_8(q)^2/nu,
    Corr_W(L,K)^2 -> 1.                                    (9.4)

This is a supercritical-dual statement. It does not contradict the separate
subcritical hypothesis of decorrelation between a diffusive shape and a
centred additive occupation fluctuation. The missing input is a bulk-density
law in the randomly delimited, exponentially long slab, under COMPLETE
component Palm rather than vertex Palm. Exponentially rare boundary selection
must not be divided out without a conditional argument.

## 10. What this does not yet say about the ultimate white pole

It is tempting to infer that the fixed-width far-tail exponent gamma_white
satisfies gamma_white/nu->1. Equation (1.5), even with convergence of every
fixed moment, does NOT by itself establish that claim: the tail limit takes
height to infinity before width.

An explicit inference counterexample uses nu_n=1/n and a mixture of positive
geometric variables: weight 1-2^(-n^2) at success parameter 1/n, weight
2^(-n^2) at parameter 1/n^2. Then L/n tends to Exp(1) with all fixed moments,
but the ultimate tail rate is `-log(1-1/n^2)`, whose ratio to nu_n tends to
zero. This is not asserted to occur in site percolation.

A useful further conjecture is the genuine site metastable relation
`gamma_white/nu_black->1`. Proving it requires uniform long-time control,
not another confirmation of the first few moments. It is recorded here as
a separate question, not silently included in (1.5).

## 11. Executed scope

The script checks 41,984 finite annular configurations with black endpoint
rings, using lifted BFS and an independently coded displacement DSU for both
colours. It checks alternation, strict height order, contact bounds, gap
sandwiches and row coverage. Those are finite geometric controls, not a
substitute for Sections 2--7.

Every nonempty periodic point pattern of periods 2 through 11 is used to
verify the discrete stationary gap identities, including the length bias.
All-height complementary moments use widths 2 through 6; joint activity
controls use widths 2 through 4. No Monte Carlo, extra heavy-width build,
external machine, or full repository suite is run. The mathematical proofs
and the explicit bulk-filling conjecture are kept separate.

## Sources actually used

[MZ] S. Mertens and R. M. Ziff, *Percolation in Finite Matching Lattices*,
arXiv:1603.07289v2, Section II (surrounding-component connectivity and the
single/spiral/cross count classification). HTML relevant section read.
https://arxiv.org/html/1603.07289v2

[AV] T. Antunovic and I. Veselic, *Sharpness of the phase transition and
exponential decay of the subcritical cluster size for percolation on
quasi-transitive graphs*, arXiv:0707.1089v3, Theorems 2--3 and Section 3
site/BK conventions. Relevant theorem statements and model definitions read.
https://arxiv.org/html/0707.1089v3

[AGG] R. Arratia, L. Goldstein and L. Gordon, *Two Moments Suffice for Poisson
Approximations: The Chen--Stein Method*, Ann. Probab. 17 (1989), 9--25,
Theorem 2, printed p.11. The PDF statement and total-variation convention
were visually read, not reconstructed from the parsed formula alone.
https://dornsife.usc.edu/larry-goldstein/wp-content/uploads/sites/221/2023/06/AGG-1.pdf

[Palm context] G. Nieuwenhuis, *Bridging the gap between a stationary point
process and its Palm distribution*, Statistica Neerlandica (1994).
Abstract only; not used as the proof of an inversion formula. Equations
(3.1),(6.1),(7.4),(7.5) are derived by counting in this note.
https://doi.org/10.1111/j.1467-9574.1994.tb01430.x
