# Winding-cluster intensities, Poisson windows, and the two birth fluctuations

2026-09-13. Continuation of the SAME geometric-balance paper, PR #739.
Author-supplied proof, not an independent referee acceptance or priority claim.
The earlier axial mass/centre argument is used explicitly. The new step is a
local, once-per-component intensity and a two-colour Poisson approximation.
No fixed-width continuum identification or unproved Ornstein--Zernike
prefactor is used.

## 1. What is settled, and what is conditional

Let G4 be the NN square graph and G8 its eight-neighbour matching graph. All
occupations are independent SITES. On the axial w-by-m torus suppose

    w -> infinity,       log(m)/w -> d in (0,infinity).

Use the earlier manuscript's mass

    kappa_G(p) = lim_n -log Pr_p^G(0 <-> n e1)/n,

and its continuous strictly decreasing inverse on the subcritical interval.
Write a=kappa_4^{-1}(d), c=kappa_8^{-1}(d), b=1-c. The earlier result gives
T1 -> a and T2 -> b, with 0<a<pc<b<1. None of these facts is re-labelled as an
independent result of the finite controls below.

This continuation proves, for the actual site models:

* Near a (or near c for G8), the number of winding COMPONENTS is asymptotically
  Poisson with mean m nu_w^G(p), where nu is a spatial component density on the
  infinite cylinder C_w x Z. It counts neither paths nor winding vertices.
* nu has the same exponential rate as the planar mass:
  -log nu_w^G(p)/w -> kappa_G(p). A window of w^2 interior rows approximates
  it to an error exp[-c_I w^2], up to polynomial factors, uniformly on a fixed
  compact subcritical parameter interval I.
* The two colour counts, using the SAME uniform labels at the two separated
  birth windows, converge jointly to INDEPENDENT Poisson laws.
* Consequently, birth CDFs at specified intensity levels have opposite Gumbel
  forms for every d>0. A separate uniform cluster-volume and semiconvexity
  argument upgrades this to median-centred affine 1/w Gumbel limits for all
  d outside an at-most-countable exceptional set. This does not require an OZ
  prefactor. At exceptional d we give a convex-log-intensity subsequential
  classification rather than pretending that a unique affine law is proved.
* Even in the actual site model, d alone does not fix the CDF at p=a(d).
  By changing only the subexponential length factor, any boundary CDF value in
  [0,1] is possible.

The exceptional set refers to possible nondifferentiability of the mass; it
is NOT asserted to be nonempty. All limits here fix d>0 and send w to infinity,
with m growing exponentially. The fixed-width 2/3/4 laboratory is not this limit.

A later section derives the displacement FROM THE INFINITE centre, including
a possible log(w)/w term, UNDER a separate prefactor hypothesis. That
hypothesis is not proved for G4/G8 here. The affine fluctuation theorem around
the TRUE FINITE MEDIAN is proved at regular d without this hypothesis. The
distinction between locating a centre and resolving fluctuations around it is
part of the mathematical conclusion.

## 2. External tools and reused inputs

[AV] supplies uniform subcritical one-arm decay on compact p intervals and
states Harris/FKG and BK for both site and bond product spaces. We use its
site statement, not a square-bond numerical critical probability.

[AGG], Theorem 2, supplies Poisson PROCESS approximation for locally dependent
indicators. In this note d_TV is sup_A |P(A)-Q(A)| (half the L1 convention).
For indicators I_i with means pi_i, dependency neighbourhoods B_i including i,
put

    b1 = sum_i sum_{j in B_i} pi_i pi_j,
    b2 = sum_i sum_{j in B_i,j != i} E(I_i I_j).

When I_i is independent of the joint family outside B_i (b3=0), the published
process bound implies

    d_TV(Law((I_i)), product_i Poi(pi_i)) <= 2(b1+b2).             (2.1)

Its contractions give the same bound for sums or finitely many typed sums.
This factor of two has been checked against the paper's doubled-TV convention,
not inferred from OCR of a formula. Zero-mean coordinates may be discarded.

The preceding `exponential-birth-centres.md` proves two model-specific inputs:
for a cylinder strip of t rows the positive-winding event has upper bound

    2p w^3 (t+2w+4) exp[-(w-1) kappa_G(p)],                      (2.2)

with an inessential enlargement of the boundary count; and a fixed-height,
seam-closed ring can be built with probability >=exp[-(kappa_G(p)+eps)w]
for each fixed subcritical p and eps>0. The first-span proof uses only planar
edges of an injecting w-by-w vertex square. The matching diagonals do not
invalidate it. We recall where these enter below.

Digital Alexander duality is the earlier finite identity

    r_4(omega)+r_8(omega^c)=2.                                   (2.3)

No correlation-length differentiability, OZ amplitude, Gumbel hypothesis, or
independence of the finite birth times is included among these inputs.

## 3. One anchor per full winding component

Work first on the infinite cylinder C_w x Z, periodic horizontally and free
vertically. Because a whole empty row has probability (1-p)^w>0 and separated
rows are independent, every occupied component is vertically bounded almost
surely for fixed w and p<1, on BOTH graphs (all steps have vertical increment
at most one).

For any component with nonzero horizontal winding, let j be its lowest row.
Choose as its anchor (j,x), where x is the smallest label in {0,...,w-1}
among its vertices in row j. This arbitrary horizontal tie-break is only a
counting convention. In particular, individual anchor probabilities need NOT
be equal in x. The sum over x is independent of the choice of tie-break.
Define

    nu_w^G(p) = E[number of winding-component anchors in row 0]. (3.1)

This is NOT an event probability when more than one anchor is possible.
It is a mean count per unit vertical length and is at most w.

For an integer H>=1, retain only components occupying at most H consecutive
rows. Their intensity is nu_{w,H}. The event that (j,x) is such an anchor can
be decided on the window

    C_w x {j-1,j,...,j+H}                                       (3.2)

alone. Compute the component of (j,x) in this window; require that it meets
neither guard row j-1 nor j+H, that its lowest row is j, that it has nonzero
horizontal winding, and that x is its bottom-row tie-break.

The two guard rows are essential: without them a locally winding cluster might
join a larger component outside the window and get counted more than once.
Since vertical jumps are at most one, not meeting a guard row proves that the
computed component is the FULL infinite-cylinder component.

Thus nu_{w,H}(p) is a finite polynomial in the product measure on w(H+2) sites.
Its events include closed sites and are generally NOT increasing. No use of
Harris or BK below is made directly on these anchor events.

On a vertically periodic torus with m>4H+4, the same local test defines
indicators I_{j,x}. A component confined to at most H consecutive cyclic rows
has a unique bottom after the complementary gap, so it contributes exactly
one anchor. Translation in j, not an incorrect equality of all x marginals,
gives

    E Z_{w,m,H} = m nu_{w,H},      Z=sum_{j,x} I_{j,x}.            (3.3)

## 4. Uniform vertical localization on the cylinder

Fix a compact subcritical interval I. By domination at sup I, the local
one-arm probability a_R is bounded by A exp(-cR), uniformly on I. Set

    R=w/64,      h=floor(w/8),

and take w sufficiently large. A path crossing a cylinder band of h rows
must run from its first to its last row and hence produce an R-arm from one
of the w possible entry sites. A first-exit arm uses only the radius
R+sqrt(2) neighbourhood, which injects into C_w x Z. Hence

    Pr(cross a specified h-row band) <= w a_R <= exp(-c0 w)       (4.1)

for some c0>0 and all large w, uniformly on I.

To cross a vertical distance H-1, a path must cross floor(H/h) disjoint h-row
bands. Take the last entrance before the first exit in each band to handle
backtracking. Band crossing events use disjoint SITE sets, so are independent.
For H>=w, floor(H/h)>=4H/w, after an immaterial adjustment for endpoints.
Absorb that adjustment in constants. There are C_I,c_I>0 such that

    Pr(a specified H-row strip has a bottom/top crossing)
           <= C_I exp(-c_I H).                                 (4.2)

This is a cylinder estimate derived from planar local arms, not an unjustified
application of an infinite-plane cluster law to a periodic graph.

If a torus component cannot fit inside H consecutive cyclic rows, a lifted
path contains such a strip crossing. A component with any vertical homology
also does so. Union over the m strip positions gives

    Pr(Bad_H) <= C_I m exp(-c_I H).                             (4.3)

On Bad_H^c all nonzero homology is horizontal and every winding component is
counted by Z. In particular

    {r_G=0} = {Z=0} on Bad_H^c.                                 (4.4)

For the infinite-cylinder anchor at row zero, height greater than H forces
the specified strip crossing. At most w anchors can lie in a row, so

    0 <= nu_w - nu_{w,H} <= C_I w exp(-c_I H).                  (4.5)

We henceforth set H=w^2. There is no unknown fitted cutoff constant. The
window has O(w^3) site variables but length polynomial rather than exponential
in w. Exact enumeration of it is NOT claimed to be computationally cheap.

For every fixed w, local-polynomial exhaustion and uniform empty-row tails on
compact subsets of p<1 also prove continuity of nu_w(p). We do not assume
that this component density is globally monotone in p.

## 5. Density and planar mass have the same exponential rate

The upper bound follows from (2.2) applied to an H-row window:

    nu_{w,H} <= poly(w) exp[-(w-1) kappa_G(p)].                   (5.1)

Equation (4.5), with H=w^2, is superexponentially smaller. For the lower bound,
use the earlier fixed-height D ring with probability
exp[-(kappa_G(p)+eps)w]. Unless its full component has vertical height greater
than H, a component anchor lies in one of at most H+D nearby rows. The chance
that the component containing the ring travels a distance H/2 is at most
poly(H+D) exp(-c_I H/2), by (4.2). Therefore

    (H+D) nu_w(p)
       >= exp[-(kappa_G(p)+eps)w] - poly(H+D) exp(-c_I H/2).      (5.2)

Together these prove

    -log nu_w^G(p)/w -> kappa_G(p).                             (5.3)

The same conclusion holds for p_w -> p in a compact subcritical interval.
For the lower bound use a fixed seed at p-delta and monotonicity of the seed
event, subtracting a uniform long-component tail at p+delta. For the upper
bound use monotonicity of the enclosing winding event. Let delta decrease
after the size limit, using the already proved continuity of kappa. This
avoids falsely treating nu itself as an increasing event.

## 6. Poisson approximation with a vanishing explicit error

For a fixed (j,x), the dependency neighbourhood consists of all anchors whose
windows (3.2) overlap. Its size is at most

    D_w = w(2H+3).                                               (6.1)

Disjoint windows are functions of disjoint product variables, so b3=0 exactly.
If two windows overlap, their union lies in a cylinder band of at most 2H+3
rows (including an appropriate lift across the vertical seam). Let

    B_w(p) = 2p w^3 (4H+4w+4) exp[-(w-1) kappa_G(p)].             (6.2)

By first-span counting this bounds the probability of any winding in that
union; it also bounds each anchor probability. If I_i=I_j=1 for distinct
anchors, they are DIFFERENT full components. Each contains a winding witness,
and these occupied witness sets are disjoint. Let E be the increasing event
that the enlarged band contains a horizontal winding. Then

    {I_i=I_j=1} subset E square E,
    E(I_i I_j) <= Pr(E square E) <= Pr(E)^2 <= B_w(p)^2.          (6.3)

The middle inequality is the site BK inequality. Anchors themselves can be
positively correlated because they share CLOSED guard sites; treating them
as increasing would be wrong. The executable controls exhibit that effect.

There are mw anchor indices. Thus b1,b2 are each at most mw D_w B_w^2.
Equation (2.1), contraction to the count, localization, and a coupling of
Poisson laws with nearby means give

    d_TV(Law(actual winding-component count), Poi(m nu_w))
       <= C_I m(w+1) exp(-c_I w^2)
          + 4 m w^2(2w^2+3) B_w(p)^2.                          (6.4)

For the count outside Bad_H, define it as the number of components with any
nonzero ambient homology; its value on Bad_H changes the comparison by at
most Pr(Bad_H). This makes (6.4) an assertion about the true finite graph.

Choose a compact interval I around c_G(d)=kappa_G^{-1}(d) so small that

    2 inf_{p in I} kappa_G(p) > d.                               (6.5)

If log m/w -> d, the RHS of (6.4) tends to zero exponentially in w, uniformly
on I. Polynomial factors in (6.4) are harmless; H=w^2 and the first term is
superexponentially small. In particular

    Pr_p(r_G=0) = exp[-m nu_w^G(p)] + o(1), uniformly on I.       (6.6)

No limit of m nu_w is needed for this absolute-error assertion. When
m nu_w(p_w)->lambda in (0,infinity), the count is Poi(lambda) asymptotically.
When the mean tends to zero/infinity, the corresponding void probabilities
follow as well. A finite real-rate numerical constant is not supplied here.

The PROCESS version of (2.1) yields more: after projecting anchor row j to
j/m on the unit circle, the winding-component positions converge to a
homogeneous Poisson process of intensity lambda. For a fixed arc the mean is
its number of discrete rows times nu, tending to lambda times its length;
disjoint arcs have independent Poisson limits. No uniform distribution of the
horizontal tie-break x is asserted.

## 7. Why the two actual births decouple in these windows

Use the SAME independent uniform site labels U_v. At a lower parameter p_1
call U_v<=p_1 black. At an upper parameter p_2 call U_v>p_2 white and use G8.
Choose p_1 near a and p_2 near b; then p_1<p_2 for large w. The two colours
are NOT independent at a site. There is a middle category.

Apply the anchor construction to both graphs, using the same cutoff H=w^2.
The combined dependency graph still consists of overlapping windows, now
with a colour tag. Same-colour pairs obey (6.3). A black enclosing winding
event E4 is decreasing in the U variables, while a white enclosing winding
event E8 is increasing in them. Product positive association therefore gives

    Pr(E4 intersect E8) <= Pr(E4) Pr(E8).                        (7.1)

Again, apply this to the enclosing increasing-in-colour events, NOT to the
nonmonotone anchors. Consequently all cross-colour neighbouring pairs have
bound B4 B8. For the combined indicators,

    b1,b2 <= m w^2(2H+3)(B4+B8)^2.                              (7.2)

Choose the two compact intervals so both minimum masses exceed d/2. The
joint Poisson process error tends to zero. If

    m nu_w^4(p_1,w) -> lambda_1,
    m nu_w^8(1-p_2,w) -> lambda_2,                              (7.3)

the two winding-component counts converge to independent Poi(lambda_1) and
Poi(lambda_2), despite being coupled through the same U_v.

Duality identifies T1>p1 with no black winding, and T2<=p2 with no white
winding. Hence

    Pr(T1<=p1, T2<=p2)
       -> (1-exp(-lambda_1)) exp(-lambda_2).                    (7.4)

This is an asymptotic independence result for the two separated transition
windows, not a statement that finite birth times are independent, not a
process-in-p convergence theorem, and not a result uniform as d decreases to
zero and the two windows merge.

## 8. The intensity-clock Gumbel law, without guessing a prefactor

Put Lambda_w^G(p)=m nu_w^G(p). For each fixed lambda>0 define p_w^G(lambda)
as the FIRST parameter in a fixed small interval around c_G(d) where
Lambda_w^G reaches lambda. It exists for large w: the mass rate makes Lambda
tend to zero at the lower endpoint and infinity at the upper endpoint.
Continuity makes its value exactly lambda. Its first-hit definition is
monotone in lambda even if nu is not globally monotone. Bounds using the
monotone enclosing winding event show p_w^G(lambda)->c_G(d).

Define deterministic level coordinates

    l_w(x)=p_w^4(exp(x)),
    u_w(y)=1-p_w^8(exp(-y)).                                    (8.1)

Both increase with their displayed argument. Equations (6.6)--(7.4) give

    Pr(T1<=l_w(x)) -> 1-exp[-exp(x)],
    Pr(T2<=u_w(y)) -> exp[-exp(-y)],
    Pr(T1<=l_w(x), T2<=u_w(y))
          -> (1-exp[-exp(x)]) exp[-exp(-y)].                    (8.2)

Thus the intensity-level CDFs have opposite Gumbel forms and factorize. The
clock here is a cylinder component density determined in polynomial-height
windows, not -log of the full-torus survival CDF defined tautologically.

Equivalent fixed quantile calibrations are

    first-birth u-quantile:       Lambda_4 = -log(1-u),
    second-birth u-quantile:      Lambda_8 = -log(u),
    either marginal median:      Lambda = log(2).               (8.3)

For the equal mixture, at each fixed u != 1/2 the local approximation is

    u<1/2:   Lambda_4 = -log(1-2u),
    u>1/2:   Lambda_8 = -log(2u-1).                             (8.4)

In the lower window the whole rank law is

    P0=exp(-Lambda_4)+o(1),  P1=1-exp(-Lambda_4)+o(1),  P2=o(1),

because localization excludes vertical homology. In the upper window duality
instead gives

    P0=o(1),  P1=1-exp(-Lambda_8)+o(1),  P2=exp(-Lambda_8)+o(1).

These identify both edges of the rank-one plateau; they are not a replacement
for the rare-odds analysis at the matching median in its interior.

Here 'quantile calibration' means the CDF at the specified intensity inverse
converges to u. A quantitative p-error requires control of the local inverse
clock; (8.3) is NOT silently called an affine quantile expansion.

The matching median lies inside the rank-one plateau, not in either window.
It remains a rare-odds balance problem and is not determined by averaging
(8.3). This paper's root theorem is still a separate input to its pc limit.

## 9. A genuine non-universality at the boundary

Fix a subcritical p0 and d=kappa_G(p0). For any t>0 choose

    m_w(t)=floor(t/nu_w^G(p0)).                                 (9.1)

By (5.3), log m_w(t)/w -> d and m_w(t)>=w eventually. By (6.6),

    Pr_{p0}(r_G>0) -> 1-exp(-t).                                (9.2)

Thus every boundary probability in (0,1) occurs in the ACTUAL site model
while keeping the same exponential geometry rate d. Taking t_w=exp(-sqrt w)
or exp(sqrt w) similarly gives 0 or 1, respectively. The absolute Poisson
error still vanishes since log m/w -> d.

This is not a toy with a different percolation rule. It proves that d by
itself cannot select the boundary CDF or a universal finite-centre correction.
It does not show that any pair of black and white boundary probabilities can
be prescribed simultaneously by the one common length m; their density ratio
would additionally matter.

At a first-birth median a_w,m in I, (6.6) does prove

    m nu_w^4(a_w,m) -> log(2),                                  (9.3)

and similarly for the white intensity at the second-birth median. This is a
sharper and correctly normalized finite-centre equation than kappa=d alone.

## 10. Affine fluctuations at all but countably many d, without OZ

We can go further than a conditional statement about p-scaling. The missing
regularity can be obtained at every differentiability point of kappa using a
uniform cylinder cluster-volume bound. This section supplies that argument.

### 10.1 Uniform subcritical cluster-volume tails on finite-width cylinders

**Lemma.** On every fixed compact I subset (0,pc(G)), there are C,c>0 and w0
such that for all w>=w0, p in I and cylinder vertices v,

    Pr_p^{C_w x Z}(|C_v|>=n) <= C exp(-cn).                     (10.1)

This does not follow by relabelling the plane graph. Here is an explicit
coarse-block proof using only the local one-arm bound already assumed.

Choose a fixed integer scale s>=4. Partition the horizontal circle into cells
of integer widths between s and 2s, and the vertical line into cells of s rows.
For w>=32s this can be done with at least sixteen horizontal cells. Each block
contains at most 2s^2 vertices. Call a block bad if some occupied vertex in it
has an occupied arm of Euclidean length at least s. The event is decided in
the block enlarged by s+sqrt(2), and its probability is at most

    delta_s <= 2s^2 A exp(-c0 s),                              (10.2)

uniformly on I. These supports inject into the cylinder. Their dependency
graph has a bounded degree independent of w and s; a conservative bound of
D=289 including each block itself suffices (only blocks within eight coarse
steps in each coordinate can overlap).

The set of blocks visited by a connected site cluster is connected in a graph
of maximum degree eight. If it contains more than 81 blocks, it cannot be
contained in the Euclidean (s+sqrt(2))-enlargement of ANY visited block: such an enlargement meets
fewer than 81 blocks. Thus every visited block is bad. There are at most
64^k connected k-block sets through a specified block: encode a canonical
spanning tree by a depth-first walk of length 2(k-1), with at most eight
choices per step. Any k-set contains at least k/D mutually independent bad
block events (greedy selection in the dependency graph). If the full visited
set is larger than k, extract a connected k-set through the root block, all
of whose blocks are still bad. Therefore

    Pr(cluster visits at least k blocks) <= 64^k delta_s^(k/D).

Take s once and for all so delta_s<=128^(-D); then this is at most 2^(-k).
Since a block holds at most 2s^2 vertices, absorb the finitely many k<=81
cases in C to prove (10.1). The constants may be extremely conservative and
are not a numerical algorithm. Horizontal periodicity, diagonal matching
steps, uneven cells, and arbitrary vertical length are all covered.

### 10.2 Cluster activities cancel the irrelevant bulk sites

A horizontally winding finite cylinder component C with lowest row zero
contributes

    p^{n(C)} (1-p)^{b(C)},
    n(C)=|C|,   b(C)=|external vertex boundary of C|.

The boundary has distinct sites; no edge multiplicity is substituted for
b(C). Sum over such connected sets C of height at most H to get exactly
nu_{w,H}. This is a sum of component probabilities, NOT a probability that
only one component exists. Each full component is represented once.

Under the probability measure on component shapes proportional to these
activities, denote expectation by E_*. Then with q=1-p,

    S_C=n(C)/p-b(C)/q,
    (log nu_{w,H})'=E_* S_C,
    (log nu_{w,H})''=Var_*(S_C)-E_*[n(C)/p^2+b(C)/q^2].           (10.3)

These identities are exact. For logit parameter z they become

    partial_z log nu=E_*[(1-p)n-pb],
    partial_z^2 log nu=Var_*[(1-p)n-pb]-p(1-p)E_*[n+b].           (10.4)

They do not involve an artificial O(wH) independent bulk occupation count.

The isolated full horizontal row has activity p^w(1-p)^{2w}; it belongs to
nu_{w,H} for every H>=1. Hence nu_{w,H}>=exp(-K_I w) uniformly on I.
Meanwhile (10.1) gives

    Pr_*(n(C)>=n) <= C w exp(K_I w-cn),

because at most w bottom-row sites can represent such a component. Splitting
the tail sum at a sufficiently large multiple of w yields

    E_* n(C) = O_I(w),       E_* b(C) = O_I(w),                 (10.5)

uniformly in H and large w; use b(C)<=8n(C). In particular for H=w^2,

    F_w(p)=(1/w)log nu_{w,w^2}(p)
    satisfies F_w''(p)>=-C_I on I.                             (10.6)

Thus the scaled log-intensities are uniformly SEMICONVEX. Discarding the
nonnegative variance in (10.3) proves the lower curvature bound; we do not
need a variance asymptotic or cluster renewal theorem.

### 10.3 Consequences for the mass and its slopes

Section 5 gives F_w(p)->-kappa_G(p). Equation (10.6) says that
F_w(p)+C_I p^2/2 is convex. The finite pointwise limit is convex too, so
-kappa is locally semiconvex (equivalently kappa is locally semiconcave).
Convex secant bounds imply local uniform convergence. The same bounds imply:
at every differentiability point a of kappa, for EVERY sequence p_w->a,

    F_w'(p_w) -> -kappa_G'(a).                                 (10.7)

For completeness, bound the derivative of the convex function at p_w between
its secant slopes with endpoints p_w+-delta. Use local uniform convergence,
then let delta decrease to zero. This works for moving p_w, not merely fixed p.
In one dimension a convex function has at most countably many derivative
jumps. Hence kappa is differentiable outside an at-most-countable set of p.

The slope does not vanish at a differentiability point. The earlier FK
argument actually gives the quantitative comparison

    kappa(p)-kappa(q) >= (q-p) kappa(q)/rho,  p<q<pc,            (10.8)

where rho is the universal constant in the stated Friedgut--Kalai theorem.
To prove it, suppose the inequality fails. Choose d'<kappa(q) sufficiently
close to kappa(q) that rho(kappa(p)-d')/d'<q-p. The earlier winding rate on
m=ceil(exp(d'w)), together with FK, forces the q-winding probability to one;
the first-span upper bound forces it to zero. This is a contradiction.
Dividing (10.8) by q-p proves

    v_G(a):=-kappa_G'(a) >= kappa_G(a)/rho >0                    (10.9)

where the derivative exists. No numerical rho is calibrated. These arguments
do not assert analyticity of kappa or rule out every possible corner.

### 10.4 Actual median-centred Gumbel theorem

Let a_w,m and b_w,m be the TRUE medians of T1 and T2. For a regular d, meaning
that both kappa_4 at a(d) and kappa_8 at c(d)=1-b(d) are differentiable, put
v4=-kappa_4'(a)>0 and v8=-kappa_8'(c)>0. Then

    X_w=v4 w (T1-a_w,m),       Y_w=v8 w (T2-b_w,m)

converge jointly to INDEPENDENT variables with CDFs

    Pr(X<=x)=1-2^(-exp(x)),
    Pr(Y<=y)=2^(-exp(-y)).                                     (10.10)

Regular d excludes at most a countable subset of (0,infinity), since each
mass is a bijection and its exceptional p set is countable.

**Proof.** The centres converge to a and b by the earlier mass argument.
The Poisson formula at an exact marginal median implies
m nu_{w,w^2}(a_w,m)->log 2, and likewise on the matching side at 1-b_w,m.
For fixed x, apply (10.7) throughout the shrinking interval from a_w,m to
p_w=a_w,m+x/(v4 w), and integrate:

    log[nu_{w,w^2}(p_w)/nu_{w,w^2}(a_w,m)] -> x.

Apply (6.6). For the white side the displacement is negative in its own
occupation probability, giving -y. Apply (7.4) for joint factorization. This
proves (10.10) without a prefactor expansion or a formula for the medians'
distance from a,b. A finite-centre fluctuation theorem and a deterministic
centering correction are genuinely different results.

### 10.5 Moments and an archive-facing consequence

The earlier FK concentration about true medians has exponential tails on
scale 1/log N; here log N~dw. Consequently the w-scaled variables have
uniformly bounded moments of every fixed order. Joint convergence therefore
also gives convergence of means, variances, and products. Set

    h0=EulerGamma+log(log 2) = 0.210702744319868533594073... .

Then at regular d,

    E T1=a_w,m-h0/(v4 w)+o(1/w),
    E T2=b_w,m+h0/(v8 w)+o(1/w),
    Var(T1)=pi^2/(6v4^2 w^2)+o(w^-2),
    Var(T2)=pi^2/(6v8^2 w^2)+o(w^-2),
    Cov(T1,T2)=o(w^-2).                                        (10.11)

The quarter and three-quarter quantiles of the mixture lie within o(1/w) of
a_w,m and b_w,m: the other birth is a fixed positive distance away and the
limiting CDF (10.10) crosses its own median strictly. Hence, for the already
available rank-gap observable G=E(T2-T1),

    G-IQR(F) = [h0/w](1/v4+1/v8)+o(1/w).                        (10.12)

This refines the earlier O(1/log N) comparison. It requires regular d, the
correct ambient-rank birth statistics, and the actual mass slopes. It is not
a free numerical prediction at width two or four.


### 10.6 A scale-free consequence that does not require numerical mass slopes

Write J1=IQR(Law(T1)), J2=IQR(Law(T2)), and J=IQR(F) for the mixture.
Let

    D0=log[log(4)/log(4/3)] = 1.57253358368551918078557... .

The marginal quantile functions of (10.10) are

    x(u)=log[-log(1-u)/log 2],
    y(u)=log[log 2/(-log u)].

Both marginal IQRs in the scaled coordinate equal D0. Consequently, at every
regular fixed d,

    J1=D0/(v4 w)+o(1/w),    J2=D0/(v8 w)+o(1/w),
    [G-J]/(J1+J2) -> h0/D0 = 0.133989344651100063543034... ,             (10.13)
    Var(Tj)/Jj^2 -> pi^2/(6 D0^2)
                   = 0.665194479964364626613157... .       (10.14)

Thus a properly typed exponential-aspect birth archive can test the shape
without knowing A, beta, or the numerical mass slopes. These are statements
about population quantiles/moments, not finite-sample unbiasedness claims.
They are not predictions for fixed-aspect square sequences (d=0) or the
fixed-width 2/3/4 controls. No new data collection is commissioned here.

## 11. Exceptional d: a constrained crossover, not an arbitrary profile

The intensity-clock statement (8.2) holds for EVERY d>0. We can also constrain
all subsequential p-scaled profiles at a corner of the mass.

Let a=c_G(d), let a_w,m be the first positive-rank median for G, and define

    psi_w(x)=log[m nu_{w,w^2}(a_w,m+x/w)].

Here psi_w(0)->log(log 2), and (10.6) gives psi_w''>=-C_I/w. Secant bounds
from local uniform convergence to -kappa bound psi_w' on compact x intervals
between the limiting one-sided slopes

    v_-=-kappa'_-(a),        v_+=-kappa'_+(a),
    0<v_-<=v_+<infinity.

Every subsequence has a further subsequence on which psi_w converges locally
uniformly to an increasing convex function psi with

    psi(0)=log(log 2),       v_-<=one-sided psi'<=v_+ .          (11.1)

Its limiting first-birth CDF is

    1-exp[-exp(psi(x))].                                       (11.2)

For the second birth use the matching clock at -y. Joint limits still
factorize by (7.4). When v_-=v_+, (11.1) forces a line and recovers the
Gumbel theorem. At a genuine corner, the exact crossover function and its
uniqueness are NOT determined by this argument. In particular we do not
claim that all d have already been proved regular.

Uniform FK moment bounds plus joint factorization also give
Cov(T1,T2)=o(w^-2) for every fixed d, even if the two marginal scaled laws
require subsequences. There is no contradiction with the finite positive
covariances measured in the controls. No statement is uniform in d down to
zero, where the windows merge and different critical scaling may intervene.

## 12. The remaining deterministic centering correction

To locate the TRUE FINITE medians relative to a(d),b(d), not merely to resolve
fluctuations about them, a prefactor theorem remains necessary. Its consequence
can be stated exactly, without pretending to have computed that prefactor.

Suppose near a=c_G(d), uniformly for |p-a|<=C log(w)/w,

    nu_w(p)=A(p) w^(-beta) exp[-w kappa(p)] (1+o(1)),             (12.1)

where A(a)>0, A is continuous, beta is fixed, kappa is twice differentiable
near a, and v=-kappa'(a)>0. Suppose

    log m=dw+gamma log w+c0+o(1).                               (12.2)

Then a first-birth u-quantile satisfies

    Q_G(u)=a+[(beta-gamma)log w-c0-log A(a)+log(-log(1-u))]/(vw)
                +o(1/w).                                     (12.3)

For the marginal median use log(log 2). Apply the analogous formula to the
white site graph and reflect p=1-q for the second birth. At a location with
log[m nu_w(p)]=x the law is 1-exp[-exp(x)], while around the true median it
is (10.10). The distinct constants are not interchangeable.

No A or beta value for square-site winding COMPONENT density has been derived
here. In particular, a plane two-point OZ factor w^-1/2 cannot simply be
copied into (12.1): periodic seam closure and once-per-component counting
change the normalization. The existence of a pure polynomial prefactor is itself left unproved, not
only the numerical values of A and beta. This unresolved centering problem does NOT
invalidate the proved regular-d finite-median fluctuation law.

## 13. What the calculation actually checked

The accompanying script has NO dependency on earlier state certificates.
It uses physical lifted graph traversal, retaining parallel periodic edges,
and a second, window-restricted implementation of component anchoring.

* 75,776 graph/configuration pairs across 2x5, 2x6, 3x5 and both graphs.
  Every eligible full component gets exactly one local anchor. Void mismatch
  is contained in the explicit localization-failure event.
* Local window polynomials independently reproduce the expected number of
  anchors divided by m at two rational probabilities. Marginals, neighbour
  pairs, count laws, and b1/b2 are exact fractions.
* All 3^8=6,561 low/middle/high assignments on 2x4, with probabilities
  1/4,1/2,1/4, check the common-label coupling and typed count law. The small
  system's count covariance is NONZERO (4831/1048576); no finite independence
  is inferred. Opposite association is checked for the enclosing winding
  events, not falsely asserted for the anchor counts.
* 8,192 deterministic full/empty-row masks on 3x12 check seams, cutoffs and
  disjoint-window neighbourhoods without random sampling.
* Component-shape activity sums are checked against full-window occupancy
  sums at widths/cutoffs (2,1), (2,2), (3,2), on both graphs. That is 168
  interior candidate shapes and 8,832 surrounding-window configurations.
  At three rational p values, 18 exact checks compare first and second
  log-derivatives via distinct external boundary counts against derivatives
  of the full-window polynomial. Logit chain rules and the curvature lower
  bound agree as well. These checks test (10.3)--(10.4), not the asymptotic
  uniform cluster-volume proof.
* For one-row isolated loops the exact local intensity is
  nu_{w,1}(p)=p^w(1-p)^(2w). Two loops two rows apart have joint probability
  p^(2w)(1-p)^(3w), strictly GREATER than the product of their anchor marginals.
  This is the closed-guard countercontrol for an invalid direct BK argument.

These tiny cutoffs are NOT H=w^2 and their void/true-rank discrepancy need not
be small. The output records that discrepancy rather than claiming a tiny
system verifies the asymptotic localization constants. Total variation against
Poisson uses floating exponentials as a diagnostic; all site probabilities
and the published-theorem RHS are stored exactly. No finite check proves the
asymptotic result, its source inputs, or publication novelty.

## 14. Position within the same paper

The present closure is:

    geometric consistency -> two mass-defined centres -> local component
    intensity -> joint Poisson windows -> natural-clock laws at every d ->
    uniform component tails/semiconvexity -> median-centred affine laws at
    regular d -> explicit leading fluctuation moments.

No new width ladder is required for this proof. At regular d, the affine
finite-median law is now proved. The remaining numerical centre displacement
is a local prefactor expansion of the winding-component density; the possible
countable exceptional d need separate regularity information. Merely fitting
a Gumbel curve or copying a plane two-point prefactor would not supply either.
The unconditional and conditional statements above are kept separate.

This is not presented as a novel Chen--Stein method. [AGG] already explicitly
explains declumping, extremes, and process approximation. [DL] is nearby
bond-wedge/rectangle work with inverse-correlation and Poisson arguments.
The specific model work here is physical SITE winding-component anchoring,
vertical localization on a periodic cylinder, and the opposite-colour
common-label joint limit. No systematic priority certification is claimed.

## References and actual reading

[AGG] R. Arratia, L. Goldstein, L. Gordon, *Two Moments Suffice for Poisson
Approximations: The Chen--Stein Method*, Annals of Probability 17 (1989), 9--25.
Author-hosted PDF, Section 2, Theorems 1--2; printed pages 10--11 rendered and
checked. Uses the doubled total-variation convention; (2.1) uses half of it.
https://dornsife.usc.edu/larry-goldstein/wp-content/uploads/sites/221/2023/06/AGG-1.pdf
DOI: 10.1214/aop/1176991491.

[AV] T. Antunovic, I. Veselic, *Sharpness of the phase transition and exponential
decay of the subcritical cluster size for percolation on quasi-transitive
graphs*, J. Stat. Phys. 130 (2008), 983--1009. Primary HTML, Theorems 2--3 and
Section 3 Fundamental Tools, especially product-site Harris and BK.
https://arxiv.org/html/0707.1089v3

[FK/DKS] P. Duncan, M. Kahle, B. Schweinhart, *Homological percolation on a
torus: plaquettes and permutohedra*, Theorem 6 reproduces Friedgut--Kalai.
Used for the earlier two-birth concentration, the quantitative mass-slope
inequality, and uniform integrability of median-centred fluctuations. It is
not the model-specific Poisson proof.
https://arxiv.org/html/2011.11903v4

[DL] M. Damron, W.-K. Lam, *Asymptotics for first passage percolation on
logarithmic subgraphs of Z^2*, arXiv:2502.18235v1. Sections 1--2 and the Poisson
comparison in Section 5 read for scope. Bond/open-boundary context, not a
substitute for site-periodic closure.
https://arxiv.org/html/2502.18235v1

[Previous working proof] `exponential-birth-centres.md`, owner handoff on
2026-09-13, recorded in #739 comment 5650571466. It supplies the axial mass
rate and its continuity/inversion with named inputs. The local control code
in this delivery does not depend on that earlier script.
