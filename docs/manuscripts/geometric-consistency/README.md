# When the matching root locates criticality but the threshold law does not

**Consolidated research manuscript — 13 September 2026.**

This document replaces a dispersed reading path with one mathematical argument.
It consolidates the probability results of #613, #735 and #736, and supplies a
new oblique winding construction for the necessity direction of the full-law
criterion. It does not merge those PRs or import their unrelated representation,
Jordan, source-response or sampling results. The construction and proofs below
are this delivery's author analysis; the finite controls are not an independent
referee's proof acceptance. Priority in the literature is not established.

## Abstract

We study independent square-lattice site percolation on finite tori obtained
from arbitrary full-rank integer period lattices. Let N be the number of sites,
ell the length of the shortest nonzero period, and r the rank of the occupied
ambient first-homology image. The matching root is the unique solution of
P(r=2)=P(r=0). Its convergence to the infinite-lattice critical probability is
uniform over all period lattices with ell tending to infinity, without aspect
or shear restrictions. In contrast, the mixture of the two ambient-homology
birth times converges to a point mass at the critical probability if and only
if log(N)/ell tends to zero.

The necessity statement, previously obtained only for axial rectangles in
#736, follows from an explicit necklace of axis-aligned occupied circuits and
connecting crossings around an arbitrary shortest period. Critical box crossing
and finite-product continuity make its exponential probability cost arbitrarily
small at a fixed subcritical parameter. Lattice translations pack order N/ell
independent necklaces. The construction never rotates the square-lattice
interaction and applies when the shortest period is nonprimitive in Z^2.

The inputs are site sharpness on the nearest-neighbour and matching graphs,
their critical-point relation, critical square-site box crossing, and finite
matching duality. No critical exponent, conformal limit, transfer-matrix
representation, new numerical critical point, or near-critical rate is assumed.

## 1. Model, observables and statements

Let Lambda be a rank-two sublattice of Z^2. On the flat torus

    T_Lambda = R^2/Lambda,       V_Lambda = Z^2/Lambda,
    N = [Z^2:Lambda],           ell = min_{u in Lambda, u != 0} |u|_2,

occupy each vertex independently with probability p. Use the physical nearest-
neighbour square-grid edges; the ambient embedding, including edge lifts, is
part of the model. All tori considered are honest square-cell tori: each unit
square is embedded with four distinct corners. Every sufficiently large ell
has this property. Parallel lifted edges are not identified merely because
they have the same endpoints.

For the induced occupied graph G_omega, define, over Q,

    r(omega) = rank im[H_1(G_omega) -> H_1(T_Lambda)] in {0,1,2},
    P_j(p) = Pr_p(r=j),
    M(p) = P_2(p)-P_0(p),        F(p) = E_p[r]/2 = (1+M(p))/2.

A rank-one diagonal/spiral class remains rank one. Nonzero projection on two
coordinate axes is not rank two.

Give the sites independent uniform [0,1] labels, and let T_j be the first p
at which r reaches j, for j=1,2. A single insertion can create more than one
rank, so T_1=T_2 is allowed. Pathwise,

    r(p) = 1_{T_1 <= p} + 1_{T_2 <= p}.

Thus F is the distribution function of T_J, where J is independent and
uniform on {1,2}. It is not the law of r at a fixed p. Let Q=F^{-1} on (0,1).

Both {r>0} and {r=2} are increasing nonconstant events. A positive pivotal
configuration on an empty-to-full chain, together with the finite Bernoulli
Russo formula, makes their derivatives positive at each interior p. Therefore
F and M are strictly increasing there. The endpoint values M(0)=-1, M(1)=1
give a unique root q_Lambda, and q_Lambda=Q(1/2).

### Theorem A — root consistency, consolidated from #735

For the actual nearest-neighbour square-site model,

    lim_{L -> infinity} sup_{Lambda: ell(Lambda) >= L}
           |q_Lambda - p_c^site(Z^2)| = 0.                      (A)

More precisely, for each fixed p<p_c, there are c_p>0 and L_p<infinity,
independent of N, orientation and shear, such that

    P_2^Lambda(p)/P_0^Lambda(p) <= exp[-c_p N/ell], ell>=L_p.   (A1)

At each fixed p>p_c, the reciprocal ratio satisfies the corresponding bound.

### Theorem B — sharp geometry for the whole law

Let Lambda_n be any sequence of honest integer-period tori with N_n tending
to infinity. The following are equivalent:

1. Q_n(u) tends to p_c for every fixed u in (0,1).
2. This convergence is uniform on each compact subinterval of (0,1).
3. F_n(p) tends to 0 for p<p_c and to 1 for p>p_c.
4. The mixture law of T_1,T_2 converges weakly to delta_{p_c}.
5. log(N_n)/ell_n tends to zero.                              (B)

The result includes arbitrary shear, orientation and Smith class. It does
not assert that ell tending to infinity is *necessary* for a balance-root
sequence to converge: Theorem A is a uniform sufficient statement for roots.
The necessity in Theorem B is for the complete set of fixed interior
quantiles, not for its median alone.

The new work in this document is the general-period necessity in Theorem B.
Its key estimate is given in Proposition 6 below.

## 2. Probability and topology inputs

We use the following established inputs with their actual model types.

**S: Site sharpness.** At fixed subcritical p, the probability that an occupied
origin connects to Euclidean distance R is at most C(p) exp[-c(p) R], on both
NN Z^2 and the eight-neighbour matching graph. Duminil-Copin–Tassion [DT],
Theorem 1.1(3), is printed in bond language; their Section 1.2 explicitly gives
the adaptation to site percolation on transitive graphs and refers to
Aizenman–Barsky. We use that site result, not square-bond p_c=1/2.

**D: Matching criticality.** p_c^site(NN)+p_c^site(NN+NNN)=1. The source chain
is Grimmett–Li [GL], Eq. (1.3), with p_u=p_c in the amenable case. Their main
strict-inequality theorem is not being relabelled as this identity.

**R: Critical box crossing.** For any fixed aspect ratio, critical NN
square-site crossing probabilities of axis-aligned rectangles are bounded
away from zero uniformly in scale. Zeng [Z], Theorem 1.1, states the specific
site version. Only its lower bound, at aspect ratios 4 and 14, is needed for
the new argument. The retrieved source is an arXiv preprint; journal status
and a quantitative value for its constant are not assumed.

**T: Finite matching duality.** For every configuration on an honest torus,

    r_NN(omega)+r_matching(omega^c)=2.                          (T)

For completeness, a topological route to T is as follows. Take a closed regular
neighbourhood U of the occupied NN graph. In a face, a white matching diagonal
can be replaced by a white boundary path unless its endpoints are the only
two white corners. The remaining diagonals do not cross, and facewise replacement
does not change the ambient homology image. The resulting white graph represents
the image of the complementary subsurface V=closure(T_Lambda\U). For
complementary subsurfaces, relative cohomology, excision and Poincare–Lefschetz
duality identify im H_1(V) with the intersection-orthogonal complement of
im H_1(U). The nondegenerate torus intersection form has dimension two, proving
T. This is the repository's digital-Alexander bridge, not a new result here.

We also use Harris positive association for increasing or decreasing events
under a product Bernoulli law, and continuity of probabilities of finite
cylinder events in p. The latter is simply finite polynomial continuity.

## 3. Period geometry without rotating the interaction

Choose a shortest nonzero u in Lambda and put S=|u|^2=ell^2. It is primitive
in Lambda: u=kz with z in Lambda and |k|>1 would contradict shortest length.
It need not be primitive in Z^2. Complete it to a basis (u,v), with

    det(u,v)=N>0,       |u dot v| <= ell^2/2.

Subtracting a nearest integer multiple of u from v achieves the second
condition. Since |v|>=ell, the transverse height satisfies

    h = N/ell >= sqrt(3) ell/2.                                (3.1)

Write n=(-u_y,u_x)/ell. The transverse coordinate on the continuous torus is
n dot x modulo h. On vertices it is det(u,x) modulo N, divided by ell.
Every physical matching edge has Euclidean length at most sqrt(2), so its
transverse displacement is at most sqrt(2). These coordinates do not change
the interaction graph.

Let g=gcd(|u_x|,|u_y|). The possible transverse coordinates of integer
translations have spacing

    delta = g/ell <= 1.                                        (3.2)

Indeed, Bezout gives z_0 in Z^2 with det(u,z_0)=g, and g divides N.
Thus the translation levels are delta times integers modulo h. This is the
point where ambient nonprimitivity must be retained, rather than assuming g=1.

A tube of transverse half-width W around R u / Z u embeds in T_Lambda when
2W<h. Two points differing by a nonzero multiple of v have transverse
separation at least h; the remaining identifications are exactly by u.

## 4. The root theorem: compare rare-sector rates

This section consolidates the core probability proof from #735. No new
transfer operator or census is required.

### Lemma 1 — a product lower bound for rank zero

For ell>=64 set rho=ell/64. Let a=a_rho(p) be the probability that an occupied
origin in the infinite graph has an occupied path to Euclidean distance rho.
Stop at first exit. Its support is contained in radius rho+sqrt(2), whose
diameter is less than ell. It therefore injects into every relevant quotient.
Let A_x be its translate to site x.

A nonzero ambient cycle has a lifted path escaping this ball. Consequently,
if no A_x occurs, r=0. Harris association applies to the overlapping decreasing
events A_x^c and gives

    P_0 >= Pr(intersection_x A_x^c) >= (1-a)^N.                 (4.1)

The local balls are not independent. The lower bound depends on positive
association, not on pretending their supports are disjoint.

### Lemma 2 — independent bands constrain rank two

In the transverse circle of circumference h, select

    k = floor(8N/ell^2)

successive half-open bands of physical width ell/8, shifting their boundaries
off vertices. Their vertex sets are disjoint; any leftover strip is unused.
In band j, define B_j to be an occupied path using only that band's sites,
from the lower sqrt(2) layer to the upper sqrt(2) layer.

Rank two implies a closed occupied walk with nonzero transverse winding.
Repeat its lift if necessary. For each band, take a first passage above its
upper boundary and the last preceding passage below its lower boundary.
Discard the endpoint edges crossing the boundaries. The intervening path is
in the band and starts/ends within sqrt(2) of the corresponding boundary.
This last-entry/first-exit construction permits arbitrary backtracking.
Hence {r=2} is contained in the intersection of all B_j.

Each B_j uses only its own site's variables, so these events are independent.
Physical edges between different bands are unused and do not alter this fact.

To count entrance vertices, centre a unit square at each lattice site. These
squares tile the flat torus with total area N. The squares of vertices in a
transverse layer of width sqrt(2) lie in a layer of width 2sqrt(2). The fibres
have length ell, so there are at most 2sqrt(2)ell entrance sites. We may use

    B = 4 ceil(ell).

A B_j crossing has transverse separation at least ell/8-2sqrt(2)>rho, so its
initial site witnesses a local arm. A union bound and independence give

    P_2 <= (B a)^k.                                            (4.2)

Equation (3.1) implies k>=4N/ell^2. Both (4.1) and (4.2) apply to NN and to
the matching graph, with their respective local-arm probabilities.

### Proof of Theorem A

Fix p<p_c(NN). By S, a<=C exp[-c ell/64]. For sufficiently large ell,

    log(Ba)<=-c ell/128,     a<=1/2,     2a<=c/(64ell).

Retaining the denominator lower bound, not replacing it by 1, yields

    log(P_2/P_0)
       <= k log(Ba)-N log(1-a)
       <= -cN/(32ell)+2Na
       <= -cN/(64ell).

This proves A1. At p>p_c(NN), D makes 1-p subcritical for the matching graph;
T exchanges its ranks 0 and 2 with the original ones. Apply the same estimate.
The two fixed parameters p_c-epsilon and p_c+epsilon then trap the unique
root uniformly over all Lambda with sufficiently large ell, proving A.

The conditional function H=P_2/(P_0+P_2) also has fixed interior quantiles
converging uniformly, since H<=P_2/P_0 below p_c and 1-H<=P_0/P_2 above it.
This conditional function is distinct from F.

A finite arithmetic consequence is useful but not an infinite-volume bound:
if an actual upper bound abar for a satisfies

    (4 ceil(sqrt(S)) abar)^4 < (1-abar)^S,

then P_2/P_0<1 on that finite torus. No near-critical arm bound or numerical
p_c enclosure is computed in this manuscript.

## 5. An oblique necklace of axis-aligned crossings

This is the new deterministic construction. It addresses the previously
missing necessity direction without assuming rotated RSW, conformal
invariance, a primitive vector in ambient Z^2, or a particular Smith type.

### Lemma 3 — explicit geometric construction

Let s>=8 be even and ell>=64s. Set

    n_0 = ceil(4ell/s),
    z_i = nearest_integer_coordinatewise(i u/n_0), 0<=i<=n_0.

Use rounding satisfying round(x+k)=round(x)+k for integer k, and extend
z_{i+n_0}=z_i+u. In particular z_0=0 and z_{n_0}=u. Successive centres satisfy

    |z_{i+1}-z_i|_infinity <= s/4+1 <= s/2.                    (5.1)

For z_i=(x_i,y_i), require four occupied NN crossings:

- horizontal crossings of [x_i-2s,x_i+2s] times [y_i+s,y_i+2s]
  and [x_i-2s,x_i+2s] times [y_i-2s,y_i-s];
- vertical crossings of [x_i-2s,x_i-s] times [y_i-2s,y_i+2s]
  and [x_i+s,x_i+2s] times [y_i-2s,y_i+2s].

The four crossings meet in the four corner squares. Their union contains an
occupied circuit C_i surrounding the inner square z_i+[-s,s]^2. This is the
usual planar annulus gluing: choose crossing subpaths joining the consecutive
corner intersections. The resulting closed walk travels successively through
top, right, bottom and left strips and has winding one about z_i; extracting
simple cycles leaves one enclosing the inner square's interior.

Connect C_i to C_{i+1} by requiring one more horizontal crossing of

    R_i = [min(x_i,x_{i+1})-3s, max(x_i,x_{i+1})+3s]
          times
          [max(y_i,y_{i+1})-s/2, min(y_i,y_{i+1})+s/2].         (5.2)

Its height is at least s/2 and width at most 13s/2, hence aspect ratio at most
13. A left-right crossing of R_i has endpoints outside both outer squares,
and at x=x_i and x=x_{i+1} passes strictly inside their inner squares.
Planarity of NN edges forces it to meet both C_i and C_{i+1} in occupied
vertices. This works even when the progression of centres is nearly vertical:
the connector is still horizontal, with a common central vertical interval.

There are exactly 5 n_0 crossing events. Every individual rectangle has
diameter less than 7s<ell and so injects into the quotient. A support point
is within 4s of its centre; the rounded centres are within distance 1 of the
line R u. All supports therefore lie in the transverse tube of half-width 6s
around that line. By (3.1), the tube embeds modulo u in T_Lambda.

### The periodic seam is closed, not assumed

Choose C_0,...,C_{n_0-1} in their planar lifts and define C_{n_0}=C_0+u.
The last rectangle in (5.2) connects C_{n_0-1} to this translated copy.
The connector paths and circuits connect a chosen vertex x on C_0 to x+u
on C_{n_0}. Projecting gives a closed occupied walk with ambient class u.
If the walk has repeated vertices, decompose it into cycles; their homology
classes sum to u, so at least one is nonzero.

Thus the intersection of the 5 n_0 crossing events implies r>0. It is a
sufficient event, not an equality with all winding configurations. Multiple
intersections, side contacts and path backtracking cannot destroy it.

### Lemma 4 — arbitrarily low exponential cost at a fixed subcritical p

There is c in (0,1), independent of s, orientation and u, bounding all these
individual crossing probabilities below at p_c. For the annuli use aspect 4
in R. For connectors compare in the infinite grid with a 7s-by-s/2 crossing
(aspect 14) and restrict to the shorter rectangle. The relevant event on the
torus has exactly the infinite-grid law because its own support injects.

Fix s. Only finitely many integer rectangle dimensions occur. Finite-product
continuity gives a single p_s in (0,p_c) at which each has probability at least
c/2. Integer translations do not change these probabilities. Consequently
Harris association on the overlapping crossing events gives

    Pr_{p_s}(necklace) >= (c/2)^(5 n_0)
                         >= exp[-25 log(2/c) ell/s].           (5.3)

For every eta>0, choose one even s large enough that 25 log(2/c)/s<=eta,
and then choose p_eta=p_s. For every shortest u with ell>=64s,

    Pr_{p_eta}(necklace around u) >= exp(-eta ell).             (5.4)

The quantifier order is eta -> s -> p_eta -> arbitrary large tori.
p_eta is a fixed parameter below p_c, not a size-dependent sequence.
No numerical RSW constant or estimate of p_c-p_eta is supplied.

## 6. Packing independent necklaces despite shear and nonprimitivity

### Lemma 5 — integer translations supply enough transverse space

Let delta=g/ell from (3.2), and take

    D = ceil((12s+2)/delta),          Delta = D delta.

Then 12s+2<=Delta<12s+3. Let z_0 have det(u,z_0)=g, and translate the whole
necklace by jD z_0 for j=0,...,b-1, where

    b=floor(h/Delta).

Their transverse centres are j Delta modulo h. Successive centres, including
the last-to-first circular gap, are at least Delta apart. Their tubes have
full width 12s, so the complete event supports have disjoint site sets.
No claim of independence is made for the crossings inside one necklace.

Since h>=sqrt(3)ell/2>=32sqrt(3)s and s>=8, h/Delta>2, and

    b >= h/(2Delta) >= h/(30s) = N/(30s ell).                  (6.1)

This explicit use of delta handles u=(w,0), u=(an,bn), and primitive ambient
vectors in the same construction. Taking a period basis with small display
entries is not required. A large longitudinal shift of a translated necklace
is immaterial to disjoint transverse support.

### Proposition 6 — a uniform upper bound for the rank-zero probability

For every eta>0 there are s and p_eta<p_c as in Lemma 4 such that every
integer-period torus with ell>=64s satisfies

    P_0^Lambda(p_eta)
       <= [1-exp(-eta ell)]^b
       <= exp[-N exp(-eta ell)/(30s ell)].                     (6.2)

Only one successful necklace is needed to prevent rank zero. Disjoint support
makes the b necklace events independent. Neither a success nor their union is
being asserted to imply rank two.

The appearance of both a geometric opportunity count and an exponential
cost is essential. An upper union bound failing to vanish would not yield
this lower-quantile obstruction.

## 7. Proof of the sharp full-law criterion

### Sufficiency

Suppose log(N_n)/ell_n -> 0. In particular ell_n -> infinity. At fixed p<p_c,
site sharpness and the embedded first-exit arm give

    Pr_p(r_n>0) <= C(p) N_n exp[-c(p) ell_n] -> 0.

Apply the same statement to the matching complement at fixed p>p_c and use T
to obtain Pr_p(r_n<2)->0. Hence F_n tends to the threshold step. This is the
previous sufficient argument of #613, reproduced to close the theorem.

### Necessity when ell is unbounded along a witnessing subsequence

If log(N_n)/ell_n does not tend to zero, choose a subsequence with
log(N_n)>=d ell_n for some fixed d>0. First suppose that along a further
subsequence ell_n tends to infinity. In Proposition 6 choose eta=d/2.
Then

    N_n exp(-eta ell_n)/(30s ell_n)
       >= exp(d ell_n/2)/(30s ell_n) -> infinity.

Therefore P_0(p_eta)->0 at the fixed subcritical p_eta, and

    F_n(p_eta) >= (1-P_0(p_eta))/2 -> at least 1/2.

For every fixed u<1/2, eventually Q_n(u)<=p_eta<p_c. Thus not all fixed
interior quantiles converge to p_c. This does not obstruct the median.

### Necessity when ell is bounded along a subsequence

Suppose instead ell_n<=L on a further subsequence. Choose a shortest u_n and
an axis-monotone NN path from 0 to u_n. Its length is
k_n=|u_{n,x}|+|u_{n,y}|<=sqrt(2)L. Every proper subpath has displacement of
Euclidean norm strictly less than ell_n, so it cannot be a nonzero period.
The path therefore projects to an essential cycle with k_n distinct sites.

For a site set S of k sites on a finite translation group of size N, each
translate meets at most k^2 distinct translates: an intersection requires
a translation difference in S-S. A greedy selection supplies at least
floor(N/k^2) disjoint translated supports. Fully occupying any such cycle
has probability p^k at a fixed p in (0,p_c). With k uniformly bounded and
N_n->infinity, independent repetition again gives P_0(p)->0 and the same
lower-quantile obstruction.

Every subsequence witnessing failure of the geometric condition has one of
these two further subsequences. This proves necessity without assuming ell_n
already tends to infinity.

### Equivalence of law and quantile formulations

Threshold-step convergence traps all Q_n(u), u in [epsilon,1-epsilon],
between p_c-delta and p_c+delta for sufficiently large n. It thus implies
uniform compact-quantile convergence, which implies pointwise convergence.
Conversely, if all fixed quantiles converge, then for fixed p<p_c and any
u>0, eventually p<Q_n(u), so limsup F_n(p)<=u; let u decrease to zero.
The argument above p_c is analogous. This also gives the usual weak
convergence equivalence to delta_{p_c}, since the state interval is compact.
No value or limit for F_n exactly at p_c is required.

## 8. Consequences that belong to this paper

### Root consistency and law concentration separate on genuinely tilted tori

Let

    u_n=(n,n),        v_n=(m_n,-m_n),        m_n=ceil(exp(n)).

These orthogonal period vectors have N_n=2n m_n and ell_n=sqrt(2)n, so
log(N_n)/ell_n -> 1/sqrt(2). Theorem A gives q_{Lambda_n}->p_c, while Theorem B
shows that the whole law does not concentrate; each fixed lower quantile is
bounded away from p_c along the witnessing construction. Both conclusions
concern the same microscopic site model, and u_n is nonprimitive in ambient
Z^2. No value of a displaced limiting quantile is asserted.

### Standard square and Gaussian-square families remain inside the good regime

For a square Gaussian period lattice generated by g and ig, ell=sqrt(N).
Thus log(N)/ell tends to zero, including nonprimitive Gaussian representatives.
Both root and full-law consistency hold. The obstruction is excessive volume
relative to the shortest period, not arithmetic nonprimitivity by itself.

### What remains outside the result

There is no quantitative near-critical rate, no value of a shifted quantile
on exponential-aspect sequences, and no L^-4 correction law. Critical RSW
contributes a qualitative, scale-uniform lower bound, followed by continuity
at a fixed block size. This argument does not identify a continuum operator,
original-U source map, generic-q tangent, or width-uniform spectral expansion.

Theorems A and B concern independent site occupation. Positive dependent
marks, arbitrary random-cluster measures and large externally imposed sources
need their own hypotheses; they are not inserted as corollaries here.

## 9. Relation to prior results and contribution boundary

The relevant comparison is between precise models and limit orders, not
whether a paper uses the name "matching" or "homological".

| Source / existing asset | Statement used or read | Relation to this manuscript |
|---|---|---|
| Mertens–Ziff [MZ], Section IV and matching identities | Finite matching function, square-sequence root/step convergence; empirical rapid root shift | Observable and square-sequence conclusions are prior art; not an all-period iff quoted from this source |
| Duncan–Kahle–Schweinhart [DKS], Section 1, Theorems 1–2 | Ambient-homology transitions on the specified growing cubical/permutohedral tori; i=1 cubical model is bond, 2D permutohedral site model is triangular | Prior art for ambient-image observables and homological transitions; not automatically the arbitrary-shape NN square-site law |
| Zeng [Z], Theorem 1.1 | Critical NN square-site box crossing | Imported lower bound; no new RSW theorem is claimed |
| Duminil-Copin–Tassion [DT], Theorem 1.1(3), Section 1.2 | Subcritical decay and explicit site adaptation | Imported sharpness, not relabelled bond criticality |
| Grimmett–Li [GL], introduction Eqs. (1.1)–(1.3) | Matching critical relation through uniqueness and amenability | Imported infinite-graph relation |
| Damron–Lam [DL], Section 2, Proposition 2.5 / Corollary 2.7 | Fixed-subcritical bond crossings of tall thin free rectangles, controlled by opportunity count times exp[-short side / correlation length] | Closest retrieved quantitative mechanism; ordinary crossing is not a seam-closed NN-site winding event on an arbitrary period quotient |
| #735 | Uniform root comparison on arbitrary integer-period tori | Consolidated, not counted as a new theorem in this delivery |
| #736 | Axial full-law iff | Extended by the oblique necklace and integer-translation packing |

[DL] explicitly relates its result to Grimmett's 1981 "Critical sponge
dimensions" and later subcritical-connectivity work. This is a meaningful
prior-art lead, not a claim that the broad opportunity-versus-cost mechanism
originates here. The original 1981 article was not independently read in this
delivery. No exhaustive citation graph or novelty certification was performed.
The direct primary texts read for [Z], [DT], [GL], [MZ] and [DKS] were their
stated theorem/definition sections; [DL] was read at its introduction and
Section 2 crossing statements. Full original proofs of all imported results
were not re-proved or independently refereed.

## 10. Finite controls and integration

The new script `scripts/oblique_winding_necklace.py` is standalone standard
Python. It uses exact integer centres and period cosets. Free-rectangle crossing
BFS is separate from its physical lifted-graph winding detector. Eight periods
include an axis shortest vector with ambient gcd 512, an oblique vector with
gcd 128, primitive oblique vectors, two genuinely reduced HNF examples and a
near-diagonal short-period choice. Every example is in ell>=64s with s=8.

Executed controls verify 12,990 specified crossing events, local rectangle
injectivity bounds, tube support, seam closure with primitive u-winding, and
disjoint *whole event supports* for selected first/second/last packed translates.
These are deterministic occupied witnesses, not Monte Carlo samples. A tiny
separate Fraction calculation checks Harris arithmetic but supplies no RSW
constant. Six local tests pass. No full repository CI was run.

Read this document as one probability manuscript, not as a new dispatch tree.
The working research boundary is now the validity and positioning of this
single combined statement. Keep the older narrow proofs and original raw
results in their branches; they remain useful if any proposed extension
needs revision. The repository's representation and original-U work is not
made dependent on this document and is not reopened by it.

## References

[Z] X. Zeng, *A Russo Seymour Welsh Theorem for critical site percolation on
Z^2*, arXiv:1309.2273v1, Theorem 1.1.
https://arxiv.org/html/1309.2273

[DT] H. Duminil-Copin and V. Tassion, *A new proof of the sharpness of the phase
transition for Bernoulli percolation and the Ising model*, arXiv:1502.03050v3,
Theorem 1.1 and Section 1.2.
https://arxiv.org/html/1502.03050v3

[GL] G. Grimmett and Z. Li, *Percolation critical probabilities of matching
lattice-pairs*, arXiv:2205.02734v3, introduction Eqs. (1.1)–(1.3), with the
companion proof identified there.
https://arxiv.org/html/2205.02734v3

[MZ] S. Mertens and R. M. Ziff, *Percolation in Finite Matching Lattices*,
arXiv:1603.07289v2. Equation numbering can differ between HTML and PDF.
https://arxiv.org/html/1603.07289v2

[DKS] P. Duncan, M. Kahle and B. Schweinhart, *Homological percolation on a
torus: plaquettes and permutohedra*, arXiv:2011.11903v4, Section 1.
https://arxiv.org/html/2011.11903v4

[DL] M. Damron and W.-K. Lam, *Asymptotics for first passage percolation on
logarithmic subgraphs of Z^2*, arXiv:2502.18235v2, Section 2, especially
Corollary 2.7.
https://arxiv.org/html/2502.18235v2

Repository inputs read: #735 at 9d29d014df28af7c635e6859d98a95ffe2b34d06;
#736 at 64d809b4404f80ff3f9adf9713337cc76008e92d; main AGENTS.md,
RESEARCH-FRONTIER.md and ROADMAP.md after #738. No unmerged runtime inputs are
needed for the new script or this proof.
