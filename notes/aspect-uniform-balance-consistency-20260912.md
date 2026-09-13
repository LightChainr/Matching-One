# Rectangular homological balance is consistent without an aspect-ratio bound

2026-09-12. Direct continuation of #276/#613/#716. This is a mathematical
consequence of established subcritical sharpness and the repository's digital
Alexander identity. It is NOT a new numerical value of p_c, a claimed novel
sharpness theorem, or an inference from the width-2/3/4 spectra.

## 1. The distinction and theorem

Let G_(w,m) be the square-site NN graph on the axis torus with periods (w,0),
(0,m), integers m>=w>=2. Let r be its ambient rational H1 rank, P_j=Pr_p(r=j),
M=P_2-P_0 and F=(1+M)/2. Parallel lifted edges at period two are retained.
The matching complement graph is the NN+NNN graph on the same sites.

The following established inputs are used:

1. For site percolation on each infinite graph, every fixed subcritical p has
   an exponential one-arm bound a_R(p)<=C(p) exp[-c(p)R]. Both graphs are
   locally finite transitive finite-range graphs. See Duminil-Copin--Tassion,
   arXiv:1502.03050v3, Theorem 1.1(3) and the explicit site adaptation in §1.2.
   The main theorem there is printed in bond language; §1.2 is essential here.
2. p_c(NN)+p_c(NN+NNN)=1. Grimmett--Li, RSA 65 (2024), 832--856,
   DOI 10.1002/rsa.21226, introduction Eqs. (1.1),(1.3) and the amenable
   p_u=p_c discussion, supply the matching-pair relation with its provenance.
3. The repository's honest-torus digital-Alexander identity is configurationwise
   r_NN(omega)+r_matching(omega^c)=2. This is not inferred from the new census.

**Theorem.** The unique zero p_(w,m) of M satisfies

    lim_(w->infinity) sup_(m>=w) |p_(w,m)-p_c(NN)| = 0.                 (1)

There is NO upper bound on m/w, and NO condition w/log(wm)->infinity in (1).
This improves the *root* conclusion on axis rectangles, not the full threshold-
distribution conclusion of #613. It does not extend the geometry to arbitrary
skew quotients without an additional strip construction.

In fact define the conditional rank profile

    H_(w,m)(p)=P_2(p)/(P_0(p)+P_2(p)) = Pr(r=2 | r != 1).             (2)

This is strictly increasing from zero to one. Every compact set of its interior
quantiles converges uniformly to p_c, uniformly over all m>=w. H is NOT the
birth-time mixture F, and it is not proposed as a new independent observation.

## 2. A local event with an embedded finite support

Fix R>=1 with 2R+2<=w. For either infinite graph, let a_R be the unconditional
probability that the occupied origin is connected, inside the sup-norm box
[-R,R]^2, to its boundary. The root being occupied is included.

For each torus vertex v, define A_v by the corresponding local lifted box. Its
vertices are distinct and its relevant edges coincide with the infinite graph:
2R+2<=w<=m rules out identification of the box and spurious boundary adjacency.
Therefore Pr(A_v)=a_R. If an occupied cycle has a nonzero winding class, its
lift cannot remain in that box. The initial segment up to first exit witnesses
A_v for a vertex on the cycle. In particular

    intersection_v A_v^c  is contained in {r=0}.                    (3)

All events A_v^c are decreasing. Harris positive association for product sites
gives, without pretending these overlapping boxes are independent,

    P_0 >= Pr(intersection_v A_v^c) >= (1-a_R)^(wm).                 (4)

For completeness, association on finitely many independent Bernoulli variables
follows by induction: condition on the last variable; the conditional covariance
is nonnegative by induction, and the covariance of the two conditional means is
nonnegative because they are monotone in the same direction. Products of
nonnegative decreasing indicators remain decreasing, so iterate this result.

This exponential LOWER bound is the step absent from the elementary union-bound
proof of full-law concentration. It remains meaningful even if wm a_R is huge.

## 3. Disjoint slab crossings give the other exponential rate

Take L=floor(m/(R+1)) disjoint slabs, with vertex rows

    y=k(R+1),...,k(R+1)+R,  k=0,...,L-1.

Let B_k be an occupied bottom-to-top crossing entirely in slab k, with the
horizontal coordinate periodic. The slabs share no vertex, so their B_k events
are independent under site product measure. Inter-slab edges are not part of
these events.

If r=2, some occupied cycle has nonzero vertical winding. A vertically lifted
copy of it traverses each chosen slab. Between the last visit to the bottom
before a first visit to the top lies a path contained in that slab. The NN and
matching steps both change y by at most one. Thus

    {r=2} is contained in intersection_(k=0)^(L-1) B_k.             (5)

A slab crossing starts at one of w bottom vertices. Before a lift reaches the
top row at distance R, it reaches the boundary of an embedded radius-R box
about that starting vertex (possibly through a horizontal side first). This
is precisely an A_v event, so a union bound gives Pr(B_k)<=w a_R. No assumed
independence between possible starting vertices is used.
Consequently, whenever w a_R<1,

    P_2 <= (w a_R)^floor(m/(R+1)).                                  (6)

Since m>=w>=2R+2, floor(m/(R+1))>=m/[2(R+1)]. Combining (4),(6),

    P_2/P_0 <= exp[-m gamma_w],
    gamma_w = -log(w a_R)/[2(R+1)] + w log(1-a_R).                   (7)

The second term is NEGATIVE. It must not be dropped in a finite sign test.
Given a certified upper bound alpha>=a_R, one fully rational sufficient test is

    w alpha < (1-alpha)^[2w(R+1)].                                 (8)

The executable checks compare these rationals exactly. A positive result is a
finite-width sign certificate for all m>=w, conditional on the validity of the
supplied arm upper bound. It does not itself determine a new numerical p_c.

## 4. Proof of the aspect-uniform statement

For fixed p<p_c(NN), choose R=floor(w/8) for w>=8. Exponential decay gives

    -log(w a_R)/[2(R+1)] >= [c R-log(Cw)]/[2(R+1)] -> c/2,
    w log(1-a_R) -> 0.

Hence gamma_w>=c/4>0 for all sufficiently large w. Equation (7) proves an
exponentially small ratio P_2/P_0 uniformly in m>=w, even if both probabilities
are themselves extremely small. In particular M(p)<0.

For fixed p>p_c(NN), set p*=1-p<p_c(NN+NNN) and apply exactly the same argument
to the matching graph. The configurationwise rank-sum identity interchanges
P_0 and P_2 under complementation. Thus P_0/P_2<=exp[-c* m/4] and M(p)>0.
This uses actual site sharpness on the degree-eight graph, not a bond result
or a planar-drawing assertion about its crossing diagonals.

Finite rank is increasing under occupied-set inclusion and is nonconstant. On
some empty-to-full chain an insertion changes the relevant increasing event;
that pivotal configuration has positive product probability at every interior
p. The finite Russo derivative is therefore strictly positive for Pr(r>0)
and for Pr(r=2). Thus P_0'<0, P_2'>0, and M'>0 on (0,1). Endpoints are -1,+1.
The unique M zero lies between p_c-epsilon and p_c+epsilon for all sufficiently
large w and EVERY m>=w. This proves (1).

Moreover, H<=P_2/P_0 below p_c and 1-H<=P_0/P_2 above p_c. The same inequalities
uniformly trap all H quantiles in [delta,1-delta] for fixed delta>0. For a fixed
width the result does not assert that an inner m-limit of the matching roots
exists or equals a specific pTL eigenvalue crossing. If such limits exist, or
one chooses any accumulation point of the actual finite roots at each width,
all those choices converge to p_c as width tends to infinity.

## 5. Why this does not contradict the thin-torus theorem

Take w_j=j and m_j=ceil(exp(j^2)). The previous full/empty-row bounds imply
for every fixed p in (0,1): r->1, F(p)->1/2, and the equal mixture of the two
birth-time laws tends to (delta_0+delta_1)/2. In the SAME sequence, (1) gives

    p_(w_j,m_j)->p_c,
    H_(w_j,m_j)(p)->1{p>p_c},  p!=p_c.                              (9)

Thus the physical matching median can consistently target p_c while most of
its threshold mixture escapes to the endpoints. Sign/order of two rare sectors
is not equivalent to absolute concentration of the full law.
#613's original full-law theorem remains valid as stated. This result neither
shows its geometric condition necessary nor supplies a sharp condition for the
entire birth distribution. It supplies a less restrictive root theorem.

## 6. Robustness to an intrinsic source and information accounting

For a source s independent of p, zero tilted mean of X=r-1 is equivalent to

    ell(p)+2s=0,  ell=log(P_2/P_0).                                 (10)

For any rectangular sequence w->infinity,m>=w, if |s|/m->0 then the exponential
sign margins on either side of p_c dominate 2s. Its source-balanced root also
converges to p_c. More generally H levels with |logit(u)|=o(m) obey the same
trapping along that sequence. This is an asymptotic statement, not permission
to change a frozen source/observer contract to get a preferred score.

Conditioning is not a free precision gain. Write E=P_0+P_2. For a single
independent rank snapshot the Fisher information for p is exactly

    I_rank = (E')^2/[E(1-E)] + E (H')^2/[H(1-H)].                   (11)

Derive it directly from P_0=E(1-H), P_2=EH, P_1=1-E; the cross terms cancel.
The second term is E times the information in a conditionally retained sample.
At balance, M'=2E H', and that term equals (M')^2/E. A sharper-looking conditional
CDF does not remove the cost of seeing r!=1. This identity is for independent
one-p snapshots; it is not a lower bound for whole-permutation, conditional,
importance-sampling or exact-transfer estimators.

## 7. Executed checks, and scope

`scripts/rectangular_rank_odds.py` verifies the local-arm implication and each
required slab crossing on all 65,536 configurations of a 4x4 torus for BOTH NN
and NN+NNN graphs: 131,072 graph/configuration evaluations. Complementary rank
sum is checked on all 65,536 pairings. Exact probability versions of (4),(6),
and (11) are checked at p=1/10,1/3,1/2. Another 236 fixed structural graph cases
check larger embedded radii; they are not a size census or Monte Carlo run.
Three sign certificates use the elementary simple-path upper bound

    a_R <= d (d-1)^(R-1) p^(R+1),

at (d,p,w,R)=(4,1/10,16,2),(8,1/20,16,2),(4,1/4,256,32).
Their approximate positive gamma margins are .0818820,.2524822,.0635379;
the signs are checked with rationals, not these decimals.

These checks do not prove an asymptotic theorem by extrapolation. Sections 2--4
do the mathematical work. No quantitative rate in w or certified numerical
interval for p_c is obtained because near-critical arm constants are not supplied.
No all-width spectral representation, CFT exponent, or ordinary/Jordan model
identification is needed. Novelty relative to published homological-estimator
proofs has not been established.

Sources (primary text sections read):
- https://arxiv.org/html/1502.03050v3, Thm 1.1(3), §1.2 site adaptation.
- https://onlinelibrary.wiley.com/doi/10.1002/rsa.21226, Eqs. (1.1),(1.3), amenable case.
- Repository #276/#613 and #702 digital-Alexander scope; #716 supplies the
  separate thin-geometry full-law contrast, not this relative-probability proof.
