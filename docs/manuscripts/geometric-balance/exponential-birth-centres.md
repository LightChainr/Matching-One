# Locating the two births in exponential rectangles

2026-09-13. Continuation of the same geometric-balance manuscript, #739.
This note replaces the previously unspecified axial birth centres by two
well-defined inverse-correlation-length equations. It is a proof supplied in
this analysis, using the published inputs below; it is not an independent
acceptance of the older manuscript or a claim of literature priority.

## 1. Main statement

Let G4 be the nearest-neighbour square lattice and G8 its matching lattice,
with steps (+-1,0), (0,+-1), (+-1,+-1). All random variables are independent
**site** occupations. Write pc(G) for the infinite-volume site threshold.
For either graph define the horizontal inverse correlation length (mass)

\[
 \kappa_G(p)=\lim_{n\to\infty}-\frac1n
  \log\Pr_p^G(0\leftrightarrow(n,0)),\quad 0<p<p_c(G).                 \tag{1}
\]

Existence, positivity, continuity, strict decrease, and the endpoint limits
needed for this definition are established in Sections 2--5, from explicitly
listed inputs. The mass is microscopic and axial. It is NOT a critical
exponent, a continuum conformal weight, or a fixed-width transfer gap.

Consider honest axial tori with periods (w,0),(0,m), where

\[
 w\to\infty,\qquad m\ge w,\qquad \frac{\log m}{w}\longrightarrow d\in(0,\infty).
                                                                        \tag{2}
\]

Let T1,T2 be the first times that the occupied G4 graph acquires ambient
homology ranks at least one and two in the uniform-label coupling. Let a_w,m
and b_w,m be their medians. There are unique numbers

\[
 \boxed{a(d)=\kappa_{G4}^{-1}(d),\qquad
 b(d)=1-\kappa_{G8}^{-1}(d).}                                         \tag{3}
\]

Then

\[
 0<a(d)<p_c(G4)<b(d)<1,\quad
 (T_1,T_2)\xrightarrow{\Pr}(a(d),b(d)),\quad
 a_{w,m}\to a(d),\quad b_{w,m}\to b(d).                              \tag{4}
\]

The equal mixture of the two birth laws converges to

\[
 \boxed{\tfrac12\delta_{a(d)}+\tfrac12\delta_{b(d)}.}                 \tag{5}
\]

Every fixed mixture quantile below 1/2 tends to a(d); every fixed quantile
above 1/2 tends to b(d). Formula (5) by itself does not locate the finite
mixture median. The separate balance-root theorem of #735/#739 still gives
that median's limit pc(G4). It need not be the midpoint of (3).

The functions a and b are continuous and strictly decreasing/increasing,
respectively. As d decreases to zero, both tend to pc(G4). As d tends to
infinity, a tends to zero and b to one. Thus the previously established two
extreme geometric regimes are the endpoints of this axial phase diagram.

In particular the whole sequence of centres converges for each d, not only
selected subsequences. The equal mixture, its integrated rank-one mass, and
its interquartile range satisfy

\[
 G_{w,m}:=\int_0^1 P_1(p)\,dp\to b(d)-a(d),\quad
 \operatorname{IQR}(F_{w,m})\to b(d)-a(d),\quad
 \operatorname{Var}(T)\to\tfrac14[b(d)-a(d)]^2.                      \tag{6}
\]

No independence of the two finite birth times is assumed or concluded.

### Published inputs and what is proved here

* **AV:** subcritical cluster-volume exponential decay, and divergence of
  susceptibility as p increases to pc, for independent site percolation on
  the two transitive, locally finite graphs. Antunovic--Veselic [AV],
  Theorems 2--3 and Proposition 5, explicitly include site percolation; see
  also their Section 6. We do not copy a square-bond theorem onto sites.
* **FK:** the Friedgut--Kalai transitive sharp-threshold inequality, in the
  formulation reproduced as Theorem 6 of Duncan--Kahle--Schweinhart [DKS]:
  for an increasing event invariant under a transitive action on N independent
  Bernoulli variables, a universal rho>0 satisfies
  \[
   \mu_p(A)>\epsilon,\quad
   q-p\ge\rho\frac{\log(1/(2\epsilon))}{\log N}
       \quad\Longrightarrow\quad\mu_q(A)>1-\epsilon.                \tag{7}
  \]
  Parameters must lie in [0,1] and 0<epsilon<1/2.
* **D:** the existing manuscript's finite identity
  r_G4(omega)+r_G8(omega^c)=2, and the infinite matching relation
  pc(G4)+pc(G8)=1. The present note reuses these stated inputs rather than
  claiming a new proof of digital duality.

The seam-closed finite-seed construction, the first-span upper bound,
comparison of their exponential rates, and the deduction of (3) are proved
below. No Ornstein--Zernike asymptotic, surface-tension differentiability,
RSW estimate or conjectured 4/3 exponent is needed for these deductions.
Sharpness and Harris correlation are established theory, not new machinery.

## 2. A normalized two-point function and a reflection bound

Let tau_p(x)=Pr_p(0<->x), including occupation of both endpoints, and put
s_p(x)=tau_p(x)/p. In particular s_p(0)=1. Conditioning a common endpoint to
be occupied leaves a product measure. Harris correlation on that measure gives

\[
 \tau_p(x+y)\ge \frac{\tau_p(x)\tau_p(y)}p,\qquad
 s_p(x+y)\ge s_p(x)s_p(y).                                           \tag{8}
\]

This argument is valid even when the connection events overlap away from the
common endpoint. It does not incorrectly treat them as independent.
Fekete's lemma along e1 gives

\[
 \kappa(p)=\inf_{n\ge1}-\frac1n\log s_p(ne_1),\qquad
 s_p(ne_1)\le e^{-n\kappa(p)}.                                      \tag{9}
\]

The limit in (1) is the same because log p/n tends to zero. A straight occupied
path gives kappa(p)<=-log p. AV supplies kappa(p)>0 for p<pc: a connection to
(n,0) needs at least n+1 occupied vertices, since every step has horizontal
increment at most one.

Both graphs are invariant under the reflections of the square. Compose a
connection to (n,y) and its reflection/translate into a connection to (2n,0).
Equation (8) implies

\[
 \frac{\tau_p(n,y)^2}{p}\le\tau_p(2n,0)
      \le p e^{-2n\kappa(p)}.
\]

Consequently

\[
 \boxed{\tau_p(x,y)\le p\exp[-\kappa(p)\max(|x|,|y|)].}              \tag{10}
\]

Only reflection and quarter-turn symmetry are used. There is no assumed
rotational invariance of the microscopic lattice, and no unproved assumption
that the most likely connecting path is straight.

## 3. A sharp exponential upper bound for any torus winding

Let f_G(w,m;p)=Pr_p(r_G>0) on the axial torus. For m>=w>=2,

\[
 \boxed{f_G(w,m;p)\le 2p\,m w^3\,e^{-(w-1)\kappa_G(p)}.}            \tag{11}
\]

**First-span proof.** A nonzero-homology closed walk has a planar lift whose
endpoint differs by (aw,bm), with (a,b) not both zero. Follow the lift until
its x-range or y-range first reaches w-1. Each physical step changes either
coordinate by at most one, so the entire retained prefix has both ranges at
most w-1. It lies in a translate of the square of vertices {0,...,w-1}^2.
This square injects into the w-by-m torus. The walk within it connects two
opposite sides of the square in one coordinate.

The relevant edges are the planar edges of this cut square. Extra periodic
edges joining its opposite sides are NOT admitted to this event. The product
law on its vertices is exactly the planar finite-box law. There are at most
N=wm translations, two crossing directions, and w^2 pairs of endpoints.
Equation (10), followed by the union bound, proves (11).

This argument retains the full distance w-1. An embedded radius-w/2 arm bound
would lose a factor two in the exponential rate and could not identify the
same centre. Backtracking, vertical winding and winding with both projections
nonzero are included. Planarity of the matching graph is not required: the
walk and its lifted displacement remain well defined despite crossing edges.

## 4. Closing finite connecting seeds into a winding ring

For fixed p in (0,pc), every epsilon>0 has a finite integer height D and w0
such that, for every w>=w0 and m>=D,

\[
 \boxed{1-f_G(w,m;p)\le
  \exp\{-\lfloor m/D\rfloor e^{-(\kappa_G(p)+\epsilon)w}\}.}         \tag{12}
\]

Here D and w0 may depend on p and epsilon, but not on w or m.

**Selecting a finite seed.** By (9), choose an integer L such that
-log s_p(Le1)/L < kappa(p)+epsilon/4. Connections in finite boxes
[-R,L+R] x [-R,R] increase to the full-plane connection as R grows. Choose a
fixed R for which the conditional finite-box probability

\[
 q=\Pr_p(0\leftrightarrow Le_1\text{ in the box}\mid0\text{ occupied})
       >e^{-(\kappa(p)+\epsilon/2)L}.                              \tag{13}
\]

No asymptotic theorem about the shape of a connecting cluster is needed for
this exhaustion step. Set D=2R+1, enlarging R if necessary to at least one.

**Closing the seam.** Write w=kL+r, 0<=r<L. Translate the seed k times along
the x-axis. On the cylinder of circumference w, each individual box injects
once w>L+2R, and all seeds use sites in the same D-row band. Require their
connection events and, when r>0, the remaining straight occupied path from
kL to w. A connected walk then runs from 0 to w e1 in the lift, including the
last endpoint, which is the translate of the first. Its projection has
nonzero homology. This is a closed ring, not an open cut-side crossing.

The first seed has probability p q. If the first j seeds occur, the common
endpoint is occupied. Applying conditional Harris as in (8) to the next seed
shows inductively that the first k seeds have probability at least p q^k.
The remainder path contributes at worst p^r by the same common-endpoint
argument. Therefore the ring probability is at least

\[
 p^{r+1}q^k\ge e^{-(\kappa(p)+\epsilon)w}                           \tag{14}
\]

for all sufficiently large w. For r=0 there is no remainder requirement;
retaining the factor p is a harmless conservative lower bound. Possible
additional overlaps at the final seam only increase the Harris lower bound.

Pack floor(m/D) disjoint D-row site bands. Their events are independent
because their SITE supports are disjoint. Unused edges between bands are
irrelevant. One successful ring already forces r_G>0, proving (12).

### The rate statement, before making any inverse in p

If log m/w -> d with d>=0, (11)--(12) imply for every fixed p<pc(G)

\[
 \boxed{\lim_{w\to\infty}\frac1w\log f_G(w,m;p)
       =-\max\{\kappa_G(p)-d,0\}.}                                \tag{15}
\]

For kappa(p)>d the event vanishes exponentially. For kappa(p)<d it occurs
with probability tending to one, and (12) bounds its absence by an
exponential of a negative exponential in w. At kappa(p)=d, (15) asserts only
that log f/w tends to zero. It does NOT specify f's limit at that parameter.

For the lower rate use 1-exp(-x)>=min{x/2,1-exp(-1)}. Then let epsilon in
(12) decrease to zero AFTER taking liminf. The upper rate follows from (11)
and f<=1. This proves (15) without continuity, strict monotonicity, or an
assumed limit of the birth medians.

## 5. Why the mass equation has exactly one solution

The following arguments avoid silently importing a bond-only
inverse-correlation-length theorem.

### 5.1 Continuity in the subcritical interval

Monotone coupling makes kappa nonincreasing. Moreover it is the infimum over
L,R of the continuous functions -log q_(L,R)(p)/L, with q the conditional
finite seed probability above. Hence kappa is upper semicontinuous; together
with monotonicity this proves left continuity.

For right continuity fix p<p0<pc. AV gives
Pr_(p0)(|C0|>n)<=C0 exp(-c0 n); the prefactor may be taken one in the cited
formulation but is immaterial. For any p<=q<=p0 and connected finite set A
containing 0,

\[
 \frac{\Pr_q(C_0=A)}{\Pr_p(C_0=A)}
  =(q/p)^{|A|}[(1-q)/(1-p)]^{|\partial A|}\le(q/p)^{|A|}.             \tag{16}
\]

Restrict the connection 0<->ne1 to clusters of size at most A0*n. The
remaining probability is bounded by the p0 tail. Using (9),

\[
 \tau_q(ne_1)\le (q/p)^{\lceil A_0n\rceil}\tau_p(ne_1)
                   +C_0e^{-c_0 A_0n+O(1)}.
\]

Taking exponential rates yields

\[
 \kappa(q)\ge\min\{\kappa(p)-A_0\log(q/p),\ c_0 A_0\}.             \tag{17}
\]

Choose A0 large enough that c0 A0>kappa(p), and then q close to p. This proves
right continuity. The use of a CLUSTER-VOLUME exponential tail here is
essential; a one-arm radius tail alone does not justify this cutoff.

### 5.2 Strict monotonicity from transitive sharp thresholds

Suppose 0<p<q<pc and kappa(p)=kappa(q)=A>0. Choose d in (0,A) so close to A
that rho*(A-d)/d<q-p, where rho is from (7). Choose zeta>0 small enough that
rho*(A-d+zeta)/d<q-p, and take m=ceil(exp(d*w)). By (15),

\[
 f_G(w,m;p)>\epsilon_w:=e^{-(A-d+\zeta)w}
\]

for large w. The event r_G>0 is increasing and invariant under all torus
translations, which act transitively on the N=wm sites. Also

\[
 \rho\frac{\log(1/(2\epsilon_w))}{\log(wm)}
 \longrightarrow\rho\frac{A-d+\zeta}{d}<q-p.
\]

Thus FK implies f_G(w,m;q)>1-epsilon_w ->1. But (11) and kappa(q)=A>d force
f_G(w,m;q)->0. This contradiction proves strict decrease. Notice the logical
order: (15) was proved without strict monotonicity, so there is no circular
use of the desired inverse-centre formula.

### 5.3 Range of the mass

For p sufficiently small, counting nonbacktracking paths of length at least n
on the degree-z graph gives kappa(p)>=-log((z-1)p). Together with the straight
path bound,

\[
 \max\{0,-\log((z-1)p)\}\le\kappa(p)\le-\log p.                    \tag{18}
\]

Thus kappa(p)->infinity as p decreases to zero. If kappa were bounded below
by c>0 as p increases to pc, (10) would imply the uniform susceptibility bound

\[
 \chi(p)=\sum_x\tau_p(x)
 \le p+8p\sum_{n\ge1}n e^{-cn}<\infty.
\]

This contradicts AV Proposition 5 and Theorem 2, which give
chi(p)->infinity as p increases to pc. Therefore kappa(p)->0.
Together with continuity and strict decrease, this makes kappa a bijection
from (0,pc) to (0,infinity), with a continuous strictly decreasing inverse.

## 6. Deducing the two centres and the rank-one phase

Let c_G(d)=kappa_G^{-1}(d). For p<c_G(d), (11) gives f_G->0. For
c_G(d)<p<pc(G), (12) gives f_G->1. Monotonicity extends the latter conclusion
to all larger p. These pointwise limits imply concentration of the first
positive-rank birth at c_G(d) and convergence of its median there.

Apply this first to G4. For G8, use its OWN site parameter t, its OWN
critical threshold, and its OWN mass kappa_G8. Finite digital duality gives

\[
 P_2^{G4}(p)=1-f_{G8}(w,m;1-p).                                    \tag{19}
\]

Consequently the second G4 birth concentrates at 1-c_G8(d), not at
1-c_G4(d). The matching critical relation places pc(G4) strictly between the
two limits. The three off-boundary phases are

\[
\begin{array}{c|c}
 p<a(d)&P_0\to1\\
 a(d)<p<b(d)&P_1\to1\\
 p>b(d)&P_2\to1.
\end{array}                                                        \tag{20}
\]

Joint concentration follows from the union bound on the two actual coupled
births. Mixture weak convergence, fixed nonmedian quantiles and (6) follow
without a finite-birth independence assumption. Since all times lie in [0,1],
convergence of the moments follows by boundedness.

The plateau in (20) contains pc. The balance root can select pc inside it
while the two births remain separated. Its location is governed by rare-sector
odds that disappear from the limiting unscaled CDF. No average of the two
centres is substituted for that root.

## 7. A small exact computation improves the old quantitative brackets

For d=log 4, the previous nonbacktracking/full-row bounds gave

    a in [1/12,1/4],      b in [3/4,27/28].

A 3-by-3 planar seed already improves them, without estimating kappa itself.
Condition its left-middle vertex occupied and ask for a connection to the
right-middle vertex using only planar edges of the box. Denote this
conditional probability by q_G(p). By (9),

\[
 \kappa_G(p)\le-\tfrac12\log q_G(p).                                \tag{21}
\]

The seed probabilities have elementary independent derivations:

\[
\begin{aligned}
 q_{G4}(p)
 &=p\{p+(1-p)(2p^3-p^6)\}
   =p^2+2p^4-2p^5-p^7+p^8,\\
 q_{G8}(p)&=p[1-(1-p)^3]=3p^2-3p^3+p^4.                            \tag{22}
\end{aligned}
\]

For G4, if the centre vertex is closed, either the complete top or complete
bottom three-site detour must be occupied. For G8, any occupied vertex in the
middle column connects the two occupied endpoints. These derivations and an
independent enumeration of the eight remaining bits give identical polynomials.

Let alpha and beta be the unique roots of q_G4(alpha)=1/16 and
q_G8(beta)=1/16. Exact rational bisection gives

\[
 \alpha=0.239805566880062101\ldots,\qquad
 \beta=0.156394361442176291\ldots.
\]

Then (21) and the strict monotonicity of kappa give the rigorous bounds

\[
 \boxed{\frac1{12}\le a(\log4)\le\alpha,\qquad
 1-\beta\le b(\log4)\le\frac{27}{28}.}                              \tag{23}
\]

In particular 1-beta=0.843605638557823708... . The JSON stores rational
outward endpoints rather than relying on rounded decimals. Alpha and beta are
roots of FINITE SEED polynomials; they are NOT computed values or point
estimates of the two infinite centre locations. The gain here is a certified
bound: actual centres still require the actual planar mass functions.

A finite seed is a one-sided certificate. Larger boxes can improve that upper
bound, but this note does not open a box/width sweep or claim an efficient
algorithm for arbitrary precision at a centre.

## 8. Relation to prior work and remaining question

The competition between a crossing's exponential cost and the number of
attempts is not a new mechanism. Grimmett's wedge/sponge programme, as recalled
and sharpened by Damron--Lam [DL], already connects a logarithmic geometry to
inverse correlation length. Their Section 1.1.2 records an inverse-
correlation-length formula for a BOND percolation wedge threshold; Section 2
studies open-boundary rectangle crossings and uses stronger two-point
asymptotics. Those are close antecedents, not a site-torus theorem that can be
copied without examining the event and boundary conditions.

Here the distinction is two periodic, ambient-rank births in the actual
square SITE / matching SITE pair. Equations (11)--(14) explicitly close the
seam and provide matching exponential rates using elementary product
arguments. We proved the mass properties needed for inversion with the
published AV and FK inputs, rather than asserting unverified site extensions
of Ornstein--Zernike results. No claim is made that inverse correlation length,
its qualitative properties, or this general entropy-cost mechanism are new.
A systematic literature priority conclusion has not been established.

What remains at fixed d is NOT whether these two axial limiting centres
exist or how to characterize them: (3) settles that in this argument. What is
not determined is the finite-width displacement and fluctuation law at
kappa(p)=d. At equality, polynomial prefactors in a true winding probability,
periodic closure weights, and subexponential factors in m can matter. Neither
(15) nor an uncalibrated FK constant supplies a Gumbel law, a 1/w shift or its
coefficient. Changing to genuinely oblique growing periods may also require
a directional mass, not merely substituting the Euclidean systole in (3).
These are possible extensions of the same paper, not new automatic queues.

## 9. Finite controls actually executed

`scripts/winding_rate_centres.py` uses only the Python standard library.
It does not import an older automaton or transfer certificate.

* Both 3-by-3 conditional seed polynomials are enumerated exactly; their roots
  against 1/16 are isolated to dyadic intervals of width 2^-64.
* An independent graph-potential traversal on 3x3, 3x4 and 4x4, for BOTH
  adjacencies, checks 140,288 graph/configuration pairs. For each of the 91,668
  nonzero-winding configurations it produces a lifted closed walk, extracts
  the first-span prefix, and independently checks a crossing on the cut planar
  square. Extra seam edges are excluded from that square.
* At 3x3 the finite seed-ring event, winding probability and conditional-Harris
  lower bound are compared at p=1/100,1/4,1/2 with Fraction arithmetic.
  The ring event is only sufficient, not all windings.
* The cluster-size truncation comparison used in right continuity is checked
  directly on the finite 3x3 product law, at p=1/5 and q=1/4.
* Tests separately compare (22) to enumeration, check strict dyadic endpoint
  signs, seam witnesses, probability inequalities and invalid inputs.

These finite checks do not prove the infinite AV/FK inputs, the asymptotic
inversion or the all-size geometry. Those claims rest on the arguments above.
Full Matching-One repository CI was not run. No Monte Carlo, new numerical pc,
critical exponent fit, or external compute was used.

## References and exact use

[AV] T. Antunovic and I. Veselic, *Sharpness of the phase transition and
exponential decay of the subcritical cluster size for percolation on
quasi-transitive graphs*, Journal of Statistical Physics 130 (2008), 983--1009.
https://arxiv.org/html/0707.1089v3
Theorems 2--3, Proposition 5, model definitions and Section 6 read in primary
HTML. Supplies subcritical site cluster-volume tails and susceptibility
divergence; no matrix or torus-centre formula is attributed to this source.

[FK] E. Friedgut and G. Kalai, *Every monotone graph property has a sharp
threshold*, Proceedings of the AMS 124 (1996), 2993--3002.
https://www.ams.org/journals/proc/1996-124-10/S0002-9939-96-03732-X/
The precise statement used here is the published theorem as reproduced in
[DKS] Theorem 6. No new direct reading of the original PDF is claimed.

[DKS] P. Duncan, M. Kahle and B. Schweinhart, *Homological percolation on a
torus: plaquettes and permutohedra*.
https://arxiv.org/html/2011.11903v4
Section 1.3, Theorems 5--6, read in primary HTML. Used for (7) and the Harris
formulation, NOT to assert that their model-specific torus theorem is the
square-site statement (3).

[DL] M. Damron and W.-K. Lam, *Asymptotics for first passage percolation on
logarithmic subgraphs of Z^2*, arXiv:2502.18235v1 (2025).
https://arxiv.org/html/2502.18235v1
Section 1.1.2 and Section 2 read in primary HTML. Context and closest mechanism:
bond wedges, correlation length and rectangle crossings, not a claimed proof
of our two site-rank birth locations. Original older wedge sources were not
independently read in this delivery.

[D] Matching critical relation and finite digital-Alexander identity are
specified, sourced and used in the parent geometric-balance manuscript at
PR739 head 758800f92fd84ee036e9b10e91a7facd54d3712c.
