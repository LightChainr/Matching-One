# A sharp geometry criterion for the full axial birth-time law

2026-09-13. New proof continuation of #613/#716; separate from #735's
arbitrary-period balance-root theorem. No Monte Carlo, fitted exponent, new
numerical critical point, or claim of priority in the literature.

## Result

Let T_n be the NN square-site torus with axial periods (w_n,0),(0,m_n),
2 <= w_n <= m_n and w_n m_n -> infinity. At parameter p let r_n be the
ambient rational first-homology rank and put

    F_n(p) = E_p[r_n]/2,       Q_n = F_n^{-1}.

Under the standard inputs listed below, the following are equivalent:

1. Q_n(u) -> p_c for every fixed u in (0,1).
2. The convergence is uniform on each compact subinterval of (0,1).
3. F_n(p) -> 0 for every p < p_c, and -> 1 for every p > p_c.
4. log(m_n)/w_n -> 0 (equivalently log(w_n m_n)/w_n -> 0).

The new step is necessity, including logarithmic-width subsequences. Merely
observing that an upper union bound stops vanishing would NOT prove necessity.
We instead construct winding events with arbitrarily small exponential cost
at a fixed p below p_c, and repeat them in independent transverse bands.

This is an axial theorem. The corresponding necessity statement for arbitrary
integer-period lattices is NOT proved here. #735's arbitrary-period result
concerns the balance root, not the whole F_n law.

## Imported inputs, not results of this note

A. Critical box crossing for NN square-site percolation: for each fixed aspect
ratio the occupied crossing probability is bounded away from zero at p_c,
uniformly in scale. Zeng, arXiv:1309.2273, Theorem 1.1, states precisely this
site result. The primary HTML theorem and gluing discussion were read on
2026-09-13. It is a preprint; no journal publication is asserted here.

B. Site subcritical exponential one-arm decay on NN and on NN+NNN, and the
matching critical-point relation p_c(NN)+p_c(NN+NNN)=1. These are the existing
#613 inputs. Duminil-Copin--Tassion, arXiv:1502.03050, Theorem 1.1(3) is
printed for bonds; section 1.2 explicitly discusses the site adaptation.
Grimmett--Li's matching relation is proved in the companion Hyperbolic site
percolation, arXiv:2203.00981, combined here with amenable p_u=p_c. The
companion's theorem body was not independently re-audited in this delivery.

C. On honest tori, the repository's deterministic digital-Alexander identity
r_NN(omega)+r_matching(omega^c)=2. This is an input, not inferred from a small
census. For these axial tori w,m >= 2 are in the stated honest scope.

Harris positive association and finite-product continuity are also used. No
conformal invariance, value of a critical exponent, correlation-length power,
or asymptotic surface-tension formula is assumed.

## 1. A deterministic ring of crossings, including the seam

Fix an integer s >= 1 and w >= 4s. Write k=floor(w/s) and divide the horizontal
circle into k consecutive integer cells with boundaries

    0=x_0 < x_1 < ... < x_k=w,
    s <= x_{i+1}-x_i <= 2s.

For example, divide w by k and distribute the remainder one unit per cell.
Work in the cylinder C_w x {0,...,s}; there is NO vertical periodic edge in this
band. Extend boundaries by x_{i+k}=x_i+w. In the planar lift require:

- V_i: an occupied vertical crossing of [x_i,x_{i+1}] x [0,s];
- H_i: an occupied horizontal crossing of [x_i,x_{i+2}] x [0,s].

There are 2k events. Every rectangle injects into the cylinder, including the
one straddling the seam. A horizontal crossing H_i meets a vertical crossing
V_i and V_{i+1}: restrict the horizontal path to an appropriate left-to-right
subpath in each cell, then use planar NN crossing intersection.

Choose these paths once on the cylinder and lift them periodically. Their
union connects V_0 successively to V_1,...,V_k=V_0+(w,0). By travelling within
the first and last vertical crossings, some vertex z is connected to z+(w,0)
in the lifted occupied graph. Projection gives a closed walk of nonzero
horizontal homology. A closed walk decomposes into cycles, so at least one
cycle has nonzero ambient image.

Thus

    intersection_i (V_i intersection H_i)
        subset {there is a horizontal essential occupied cycle in the band}. (1)

This is stronger than a path merely joining the two cut sides: the last
vertical connector explicitly closes the periodic seam. It is a sufficient
event, not an equality or a claim that every winding cluster looks this way.

## 2. Arbitrarily small exponential cost below criticality

Critical RSW supplies a constant c in (0,1), independent of s, bounding below
all the above crossing probabilities at p_c. To see uniformity without a
variable-aspect theorem, a 4s-by-s horizontal crossing implies each shorter
horizontal crossing, and an s-by-s vertical crossing fits inside each cell.
Use the minimum of these two fixed-aspect RSW constants.

Set a=c/2. For each FIXED s, only finitely many integer rectangle sizes occur.
Their probabilities are finite polynomials in p. Continuity therefore gives
one p_s in (0,p_c) at which every such crossing probability is at least a.
This p_s is independent of w and m. We do not claim a quantitative lower bound
on p_c-p_s, or a uniform-in-s neighborhood of criticality.

Harris association, not independence, applies to the overlapping rectangles.
The ring event G_{w,s} in (1) satisfies

    P_{p_s}(G_{w,s}) >= a^(2 floor(w/s))
                      >= exp[-2 log(1/a) w/s].                  (2)

Consequently for every eta > 0 there are a FIXED s and a FIXED p_eta < p_c
such that, for every w >= 4s,

    P_{p_eta}(G_{w,s}) >= exp(-eta w).                           (3)

The quantifier order is essential: choose eta, then s, then p_eta, and only
then let the torus sizes diverge. This is not a critical RSW bound silently
applied at a sequence-dependent subcritical parameter.

## 3. Independent repetition across the long direction

In the w-by-m torus take bands with vertex rows

    j(s+1), ..., j(s+1)+s,      0 <= j < b=floor(m/(s+1)).

Their vertex sets are disjoint. There can be physical edges between bands,
but each G_{w,s} uses only sites and paths inside its own band. Its probability
and independence are not affected by unused outside edges. Thus

    P_0^{w,m}(p_eta)
       <= (1-exp(-eta w))^floor(m/(s+1))
       <= exp[-floor(m/(s+1)) exp(-eta w)].                     (4)

This is an upper bound for P_0: just one successful band already forces r>0.
It does not say that a single band forces r=2.

## 4. Necessity of subexponential aspect growth

Suppose log(m_n)/w_n does not tend to zero. Extract a subsequence with
log(m_n) >= d w_n for some d>0.

If w_n -> infinity along a further subsequence, choose eta=d/2 in (3).
Then floor(m_n/(s+1)) exp(-eta w_n) -> infinity, so (4) gives P_0(p_eta)->0.
Since

    F_n(p) = (P_1+2P_2)/2 >= (1-P_0)/2,

we have liminf F_n(p_eta) >= 1/2 at a FIXED p_eta<p_c. For every fixed
u<1/2, eventually Q_n(u) <= p_eta < p_c. In particular Q_n(1/4) cannot tend
to p_c.

If instead w_n stays bounded along a further subsequence, m_n -> infinity.
At any fixed p in (0,p_c), fully occupied rows give independent winding events
with probability p^{w_n} bounded below. Hence P_0 <= (1-p^{w_n})^{m_n}->0,
and the same lower-quantile obstruction applies.

Every subsequence witnessing failure has one of these two further subsequences.
This proves 1 => 4. It does NOT prove that the median fails: Q_n(1/2) is the
balance root and has its own, stronger consistency theorem.

## 5. Sufficiency and the quantile equivalences

If log(m_n)/w_n->0 then w_n->infinity and log(w_n m_n)/w_n->0. At fixed p<p_c,
a nonzero winding cycle must escape an embedded ball of radius w_n/4-O(1).
Site exponential decay and a union bound over w_n m_n sites give

    P_p(r_n>0) <= C(p) w_n m_n exp[-c(p) w_n] -> 0.

Above p_c apply the same argument to the matching complement and use C.
This is the already established #613 sufficiency, not new work here. It yields
3. Monotonicity traps every Q_n(u), uniformly for u in [delta,1-delta], between
p_c-epsilon and p_c+epsilon for sufficiently large n; hence 3 => 2 => 1.
Together with section 4 this proves the four-way equivalence.

## 6. What the new result changes

For axial sequences with w->infinity there are now two different geometric
requirements:

    balance-root consistency:    w -> infinity            (#718/#735);
    entire birth-law convergence: log(m)/w -> 0            (this note + #613).

For m=ceil(exp(d*w)), d>0, the root remains consistent while at least every
fixed lower quantile u<1/2 is bounded away from p_c along the sequence by the
proof above. This addresses the logarithmic-width boundary, not only #716's
w=o(log m) regime. No value of the separated limiting quantiles is computed.
For m polynomial in w, all fixed quantiles converge. For m=ceil(exp(sqrt(w))),
all fixed quantiles still converge; bounded aspect ratio was never necessary.

Interpretation: local convergence of the underlying graphs to Z^2 does not
by itself make a global two-birth statistic concentrate. The number of chances
to create a short winding cycle must be compared with its probability cost.
This does not identify an irrelevant field or provide an L^-4 root-shift law.

## 7. Finite checks actually executed

`scripts/axial_ring_gluing.py` independently implements rectangle crossing BFS
and a graph-potential winding detector. Exhaustive cylinder strips at
(w,s)=(4,1),(5,1),(6,1),(7,1) cover 21,760 configurations, with zero instances
of the ring event without nonzero horizontal winding. It also checks the
product-of-marginals <= ring-event <= winding-event inequalities with Fraction
arithmetic at p=1/3,1/2,2/3. These are not binomially weighted twice.

At w=4,s=1 there are 17 ring-event configurations but 35 winding configurations,
so at p=1/2 the two probabilities are 17/256 and 35/256. This guards against
promoting the sufficient construction into an event dictionary identity.

Another 2,365 deterministic masks at (9,2),(11,2),(13,3) check uneven cells
and seam closure; 2,345 satisfy the ring event, with no implication failures.
These masks are controls, NOT samples estimating any physical probability.
Three mathematical unit tests and Python compilation passed locally.
Full repository CI was NOT run. Finite checks support the implementation;
sections 1--5 and the imported inputs, not enumeration, support the theorem.

## 8. Prior art and the next genuinely new question

Primary theorem text read this delivery:
- Zeng, https://arxiv.org/html/1309.2273, Theorem 1.1 and crossing/gluing setup.
- Kohler-Schindler--Tassion, https://arxiv.org/html/2011.04618, Theorem 1 and
  Comment 1 (site extension). Its general RSW statement is background, not
  automatically an identification of the critical point of this site model.
- Duminil-Copin--Tassion, https://arxiv.org/html/1502.03050, Theorem 1.1(3)
  together with the site-adaptation paragraph in section 1.2.

Context checked, not used as a proof of necessity:
- Easo, Sharpness and Locality for Percolation on Finite Transitive Graphs,
  https://doi.org/10.1007/s00039-025-00726-w (2025), publisher introduction and
  theorem discussion: giant bond clusters, not this rank-birth distribution.
- Duncan--Kahle--Schweinhart, https://arxiv.org/abs/2011.11903: abstract and
  metadata this delivery; the existing #613 source review handles its scope.
- Grimmett--Li, https://arxiv.org/abs/2203.00981: abstract/metadata this delivery;
  see B above for the reused theorem provenance.

A bounded search did not locate this exact axial quantile iff formulation.
That is NOT a novelty certificate: strip-percolation/finite-size-criterion and
homological-percolation citation chains still require systematic comparison.

Next target: for an arbitrary integer-period lattice with shortest period ell,
construct a bounded-thickness occupied ring around that actual period with
probability at least exp(-eta ell) for some fixed p_eta<p_c, uniformly in its
orientation and ambient primitivity, and pack disjoint such corridors. The
axial proof does not rotate the physical NN lattice to claim this for free.
Do not replace this missing geometric lemma with another size census.
