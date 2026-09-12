# Thin tori: two endpoint births, a flat rank CDF, and the logarithmic geometry boundary

2026-09-12. New completed probability analysis connecting #276/#337/#582/#613
with the finite-width spectra of #705/#707/#708 and the preceding parametric
delivery. No Monte Carlo, new width ladder, or continuum identification.
No literature novelty claim is made. The proofs below are self-contained once
the repository's digital-Alexander rank complement theorem is supplied.

## 1. Objects and a distinction that must not be lost

Take the square-site NN torus with periods `(w,0),(0,m)`, w,m >= 2. Distinct
lifted periodic edges are retained at period two. Its matching graph has the
NN and diagonal NN+NNN edges. Let r(p) in {0,1,2} be the rational ambient rank,
P_j(p)=P(r(p)=j), M=P_2-P_0, and F=E[r]/2=(1+M)/2.

Attach independent uniform labels U_v in (0,1), occupying v when U_v <= p.
Let T_1 and T_2 be the first p at which r reaches 1 and 2. Simultaneous births
are allowed in this definition. Since rank is monotone,

    r(p) = 1{T_1 <= p} + 1{T_2 <= p}.

Consequently F is the CDF of a random variable T chosen by selecting either
birth with independent probability 1/2. The distribution of r at a fixed p and
the distribution of the random threshold T are DIFFERENT probability laws.
Strict positivity of the finite Bernoulli rank derivative gives a unique
inverse Q_m(u) for every u in (0,1), at every finite m.

We first hold w FIXED and let m -> infinity. This is not the growing Gaussian
square-torus geometry used by the N145/N290/N725 productions.

## 2. A finite bound: most configurations have rank one, not zero or two

A full occupied row forces a transverse cycle, so rank zero requires that no
row be fully occupied. An empty row prohibits longitudinal winding, so rank
two requires that no row be empty. Rows are independent. Thus exactly,

    P_0(p) <= (1-p^w)^m,
    P_2(p) <= [1-(1-p)^w]^m,
    P(r != 1) <= (1-p^w)^m + [1-(1-p)^w]^m.                 (1)

For every fixed 0<p<1, both bounds vanish exponentially at fixed w. Therefore

    r(p) -> 1 in probability,
    M_m(p) -> 0,  F_m(p) -> 1/2,
    law(T) => (delta_0 + delta_1)/2.                        (2)

Convergence in p is uniform on every compact subset of (0,1). Weak convergence
of the threshold law does not assert CDF convergence at the atom at zero.

This already prevents confusing the cylinder limit with the square-torus step
at p_c. In particular, a median may converge to an interior spectral crossing
while ALL fixed quantiles below 1/2 go to zero and ALL those above 1/2 go to one.
The inverse of the limiting flat CDF does not determine the limit of its finite
medians. There is no contradiction with strict finite-volume monotonicity.

## 3. Shortest matching circuits and their exact multiplicity

Define the central trinomial coefficient

    c_w = [z^0](z^-1+1+z)^w
        = sum_{k=0}^{floor(w/2)} w!/[k! k! (w-2k)!].        (3)

For w=2,3,4 this is 3,7,19. For m>2w, the matching graph has exactly m*c_w
distinct vertex sets of size w supporting a nonzero ambient cycle.

Proof. A nonzero vertical winding requires at least m steps. A nonzero
horizontal winding requires at least w steps, because |Delta x|<=1. At size w
choose a simple cycle with nonzero gain (any nonzero closed walk decomposes
into simple graph cycles). Equality forces every horizontal step to be +1
or every step to be -1. Orient it +1, and start at its unique column-zero
vertex. At column x its height increment lies in {-1,0,1}; all w increments
sum to zero because m>w. There are c_w such words and m starting heights.
The vertex set recovers the height in each column and hence the word, so there
is no division by w and no duplicate vertex-set count. The w=2 case uses the
two distinct lifted edges between endpoint pairs. Conversely every such word
forms the required cycle. This proves the classification, not just the count.

For the primary NN graph the analogous size-w sets are exactly the m full
straight rows: no vertical step is possible in a w-step transverse cycle.

## 4. Joint limiting law of the two births

For each FIXED integer w>=2,

    (m^(1/w) T_1, m^(1/w)(1-T_2)) => (A_w,B_w),            (4)

where A_w and B_w are INDEPENDENT, with survival functions

    P(A_w>x)=exp(-x^w),
    P(B_w>y)=exp(-c_w*y^w),          x,y>=0.                (5)

Independence is a limit statement. Finite birth times are not independent.

### Step 1: nonminimal circuits are negligible on these scales

On a graph of maximum degree Delta, the number of simple cycles of length ell
is bounded above by N*Delta*(Delta-1)^(ell-1), with overcounting harmless.
Set p=x*m^(-1/w). The probability of a nonzero cycle using more than w sites
is at most a constant (depending on w) times

    m * sum_{ell>=w+1} [(Delta-1)*p]^ell = O_w(m^(-1/w)).

Eventually (Delta-1)*p<1; the same argument applies to sparse matching sites
at q=y*m^(-1/w), with Delta=8. All vertical cycles are included in the
negligible sum once m>w. Thus the first primary birth is, with asymptotically
unit probability at the tested scale, the first full-row pattern. The first
reverse matching birth is similarly the first size-w motif of section 3.

### Step 2: local rare-pattern counts converge to Poisson variables

Full primary rows have disjoint supports, so their count is exactly binomial
with parameters (m,p^w), converging to Poisson(x^w).

There are m*c_w matching motifs. Each has probability q^w. A motif intersects
only a bounded number of other motifs, the bound depending on fixed w, not m.
Distinct intersecting motifs have at least w+1 sites in their union.

For any fixed factorial moment, split selected motifs into connected components
of their support-overlap graph. All-singleton selections have disjoint supports;
their contribution tends to the corresponding power of c_w*y^w. A component
with two or more motifs has only O_w(m) placements, and costs at least q^(w+1)
instead of the q^w of a single motif. Combining it with other components gives
an O_w(m*q^(w+1)) contribution, which vanishes. This proves all fixed factorial
moments of the Poisson limit. The same argument applies to mixed factorial
moments with marked low-label full rows and high-label matching motifs.

Specifically, low motifs require U<=p and high motifs require U>=1-q.
For large m, p+q<1. Overlapping opposite-type motifs are incompatible, while
disjoint supports are independent. The excluded overlap tuples are O(m^(k-1))
instead of O(m^k); same-type overlap components vanish as above. Mixed factorial
moments therefore tend to (x^w)^a*(c_w*y^w)^b. The locally dependent counts have
uniformly bounded moments of each fixed order by this same component bound;
the moment-determinate independent Poisson limit follows. This is a direct
rare-pattern proof, not an invocation that pairwise independence suffices.

The exact complement-rank identity identifies the first reverse matching birth
with 1-T_2 (threshold endpoint conventions change no continuous probabilities).
The joint zero-count event then has limit exp(-x^w-c_w*y^w), proving (4)-(5).

No uniform growing-w Poisson estimate is asserted. In particular c_w and the
local-dependence constants grow with w.

## 5. Quantiles, expectations, and the archived integer birth ranks

For fixed u away from 0,1/2,1, inversion of the continuous strict limiting
endpoint CDFs yields

    m^(1/w) Q_m(u) -> [-log(1-2u)]^(1/w),       0<u<1/2,
    m^(1/w)[1-Q_m(u)] -> [-log(2u-1)/c_w]^(1/w), 1/2<u<1. (6)

These convergences are uniform on compact subintervals of the displayed ranges.
The signs and factors two matter: F selects either birth with mass 1/2.

All positive joint moments also converge. Indeed T_1 is no larger than the
first complete primary row time, giving the uniform tail bound

    P(m^(1/w) T_1>x) <= exp(-x^w).

The reverse matching birth is no later than a complete reverse row, giving
exactly the same dominating tail for m^(1/w)(1-T_2). Uniform integrability
of arbitrary products follows, for example, by Cauchy-Schwarz and higher moments.
In particular, with g_w=Gamma(1+1/w),

    E[T] = 1/2 + (g_w/2)*(1-c_w^(-1/w))*m^(-1/w)
                       + o(m^(-1/w)),
    E[T_2-T_1] = 1-g_w*(1+c_w^(-1/w))*m^(-1/w)
                       + o(m^(-1/w)).                  (7)

The mean of the threshold law approaches 1/2, not the interior median limit q_w
established by the width-2/3/4 spectral results.

For the repository's permutation clocks let N=w*m, and K_1,K_2 be the first
occupation ranks at which ambient rank reaches 1,2. Conditional on the random
permutation, T_i is its K_i-th uniform order statistic. The permutation and
ordered uniform values are independent. Thus E[T_i|K_i]=K_i/(N+1), and the
conditional order-statistic variance is K_i(N+1-K_i)/[(N+1)^2(N+2)]. At the
lower endpoint, E[K_1/N]=O(m^(-1/w)), so after multiplying by m^(1/w) the
mean-square discrepancy between T_1 and K_1/N vanishes; the power is
m^(1/w-1), which tends to zero even for w=2. The reverse birth rank is
N-K_2+1. Therefore (4) also holds with

    (m^(1/w) K_1/N, m^(1/w)(N-K_2+1)/N).

This does not assume uniform order-statistic noise is absent; it controls it
at the endpoint scale. Nor does it mix the Gaussian square lineage with strips.

## 6. Even a diverging shortest period is not sufficient

The finite row bound (1), unlike the preceding Poisson estimate, permits w to
vary. If w=w_m>=2 and w_m=o(log m), then for every fixed 0<p<1,

    m*p^(w_m) -> infinity,  m*(1-p)^(w_m) -> infinity.

Hence (2) continues to hold. For example w_j=j and m_j=ceil(exp(j^2)) has BOTH
periods tending to infinity, yet F_j(p)->1/2 at every interior p and the threshold
law splits between zero and one. This is a concrete counterexample to replacing
#613's geometry condition by only 'shortest period tends to infinity'.

The #613 hypothesis ell_N/log N -> infinity is not contradicted: for these
rectangles ell_N=w_j and ell_N/log N~1/j ->0. We do not claim its exact necessity
or a complete phase diagram for w comparable to log m.

One elementary partial statement at logarithmic width is available. If
w_m/log m -> a with 0<a<1/log 2, then (1) forces rank one throughout

    exp(-1/a) < p < 1-exp(-1/a).                         (8)

Outside this interval the row argument is inconclusive, not a contrary theorem.

## 7. Consequence for the intrinsic topological source

The existing rank-source definition (#337) is

    Z_top(p,s)=P_0 e^-s+P_1+P_2 e^s
              =1+M*sinh(s)+E*(cosh(s)-1), E=P_0+P_2.

On the thin sequences of section 2 or 6, E->0 and M->0 throughout every interior
compact. Thus Z_top(p,s)->1 uniformly for bounded s, and every rank-source
cumulant vanishes. Source reflection symmetry is restored TRIVIALLY by
concentration on X=r-1=0, at a whole interval of p values. It is not by itself
a marker of the two-dimensional critical point.

This is an actual square-site counterexample, not a P398 analogy or a hidden
Markov gadget. It does not invalidate the source algebra, the finite matching
root, or the growing-square critical-point bridge.

## 8. Executed controls and sources

`thin_torus_extremes.py` independently traverses minimal white subsets on 2x5,
3x7,4x9: 45+1330+58905=60,280 subsets. It finds exactly 15,49,171 minimal
motifs, with no unclassified set. Full/empty-row controls check 64+4096=4160
complete configurations. A subset DP records exact joint K_1/K_2 distributions
on 2x3 and 3x3; an independent loop over all 720 permutations checks the 2x3 DP.
Its binomial reconstruction agrees exactly with physical graph probabilities
at p=1/3,1/2,2/3. Finite continuous-birth covariances are 1123/58800 and 1/56,
explicitly NOT zero. These controls do not replace the all-fixed-width proof.

`thin_torus_spectral_controls.py` uses the already proved finite spectra to
check endpoint constants exactly by integer polynomial expansion: the leading
survival eigenvalues are 1-p^w+O(p^(w+1)) and
1-c_w*(1-p)^w+O((1-p)^(w+1)) at widths 2/3/4. High-precision values at large m
are small-matrix powers, not enumeration of a large system.

Sources actually consulted in this continuation:
- Existing issues #276, #337, #582, #613, #321 and draft #708.
- Mertens & Ziff, PRE 94 (2016) 062152, arXiv:1603.07289v2, section II,
  Eqs. (20)-(21), and section IV Eq. (31): PRIMARY_TEXT_READ at
  https://arxiv.org/html/1603.07289 . Their growing-square limit is not (2).
- Newman & Ziff, PRE 64 (2001) 016706, arXiv:cond-mat/0101295v2,
  Eqs. (1)-(2) and section II: PRIMARY_TEXT_READ at
  https://arxiv.org/html/cond-mat/0101295 . This supports the canonical/microcanonical
  reconstruction distinction; it is not cited as a proof of the new limits.
- Jacobsen, J. Phys. A 48 (2015) 454003, arXiv:1507.03027v1, fixed-width
  cylinder construction and section 6.1 Table 2: PRIMARY_TEXT_READ at
  https://arxiv.org/html/1507.03027v1 . The numerical q_2,q_3,q_4 were published.

A bounded search for primary local-dependence Poisson references encountered
blocked full text. The factorial-moment argument above is supplied rather than
asserting an unread theorem's exact hypotheses. No broad novelty search was done.
