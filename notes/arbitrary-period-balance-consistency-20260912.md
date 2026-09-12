# Matching balance on arbitrary integer-period tori

2026-09-12. Completed proof continuation of #718/#732 and #276/#613.
The result concerns **balance roots and conditional rank odds**, not the full
birth-time mixture. It does not assert a new value, a near-critical rate, an
all-width transfer representation, or publication novelty.

## 1. Statement and external inputs

Let Lambda be any rank-two sublattice of Z^2. Put

    N = [Z^2:Lambda],     ell = min{|u|_2: 0 != u in Lambda}.

On the occupied square-site NN quotient, let r be the rational rank of its
ambient H1 image, P_j=Pr_p(r=j), and M=P_2-P_0. All sufficiently large ell
are honest square-cell tori. The matching graph has the eight steps
(+-1,0),(0,+-1),(+-1,+-1); parallel lifted edges are not collapsed.

**Theorem A.** If p_Lambda is the unique zero of M, then

    lim_(L->infinity) sup_{Lambda: ell(Lambda)>=L}
              |p_Lambda - p_c^site(Z^2)| = 0.                     (1)

There is no restriction on area N relative to ell, aspect, shear, orientation,
Smith class, or primitivity of a Gaussian representative. The local site
product law and actual period lattice are essential. No arbitrary graph family
or dependent random-cluster law is included.

More quantitatively, for each fixed p<p_c(NN), there are L(p), kappa(p)>0,
independent of Lambda, such that

    P_2^Lambda(p)/P_0^Lambda(p) <= exp[-kappa(p) N/ell], ell>=L(p). (2)

Above p_c the inverse ratio obeys the analogous estimate, with matching-side
constants. The conditional CDF H=P_2/(P_0+P_2) therefore has all its fixed
interior quantiles tending to p_c uniformly in Lambda. H is NOT F=(1+M)/2.

The imported probability inputs are the same as in #718:

* Subcritical **site** one-arm exponential decay on both infinite finite-range
  transitive graphs. Duminil-Copin--Tassion, arXiv:1502.03050v3, Thm 1.1(3)
  is printed for bonds; the explicit site-adaptation discussion in section 1.2
  is part of the provenance. We do not relabel the printed bond theorem as a
  printed site theorem, or substitute square-bond p_c=1/2.
* p_c^site(NN)+p_c^site(NN+NNN)=1. Grimmett--Li,
  arXiv:2205.02734v3, introduction Eq (1.3), together with amenability
  p_u=p_c, states the relation and identifies the companion proof.
* The repository's configurationwise honest-torus digital-Alexander identity
  r_NN(omega)+r_matching(omega^c)=2. It is an input, not inferred from this
  delivery's finite checks.

The new step beyond #718 is the oblique, integer-period slab construction.

## 2. Reduce the period lattice, not the physical interaction

Choose a shortest nonzero u in Lambda. It is primitive **in Lambda**: a proper
integer multiple would contradict shortness. Complete it to a basis (u,v),
orient det(u,v)=N>0, and replace v by v-ku so that

    |u dot v| <= ell^2/2,    |v|>=ell.

Let n=(-u_y,u_x)/ell and h=n dot v=N/ell. Pythagoras gives

    h >= sqrt(3) ell/2.                                         (3)

This is a change of period basis only. Physical edges remain the original NN
or NN+NNN steps, not rotated nearest-neighbour edges on a new grid.
Neither u nor v is required to be primitive in ambient Z^2. In particular an
axis period u=(w,0) is completely legitimate.

The map theta(x)=n dot x modulo h is well-defined on the continuous torus.
Its fibres have length ell. On lattice vertices it can be implemented exactly
as q(x)=det(u,x) modulo N; each local edge e changes its lifted q by det(u,e).
Every edge has physical length and projected displacement at most sqrt(2).
An ambient rank-two cycle space necessarily contains a cycle whose theta
winding is nonzero. Repeating a lift traverses every transverse band.

## 3. A finite-support local arm gives a lower bound for P_0

Set r=ell/64. On the infinite graph let a_r(p) be the probability that the
occupied origin has an occupied path reaching Euclidean distance at least r.
Stop the path on first exit. The event is measurable in radius r+sqrt(2).
For ell>=64, twice that support radius is smaller than ell, so it injects
into every quotient under consideration. Define the translated local event A_x
at each of the N torus vertices. Its probability is exactly a_r(p).

A nonzero ambient cycle has a lift escaping such a ball. Thus

    intersection_x A_x^c subset {r=0}.

These are overlapping, decreasing events. Harris association for independent
sites, iterated over the N indicators, gives

    P_0 >= (1-a_r)^N.                                           (4)

No independence of overlapping local balls is assumed. Association itself
follows by induction over Bernoulli coordinates: the conditional covariance
is nonnegative, and so is the covariance of the two monotone conditional means.

## 4. Vertex-disjoint oblique slabs give an upper bound for P_2

Choose physical slab width b=ell/8 and

    k=floor(h/b)=floor(8N/ell^2).

Translate their boundaries by a generic common offset to avoid vertices.
For an exact implementation, raw q-width is ell^2/8 and raw offset 1/17;
q-values are integers, so no boundary passes through a vertex. The k half-open
bands have disjoint vertex sets. Leave any residual strip unused.

In each band let B_j be existence of an occupied path contained in that band,
from its lower sqrt(2)-layer to its upper sqrt(2)-layer. Edges crossing a band
boundary are not part of that band's event. The last-entry/first-exit portion
of any nonzero theta-winding cycle provides a B_j crossing in every band.
This remains true for a cycle that backtracks. Therefore

    {r=2} subset intersection_j B_j.                            (5)

Each B_j is measurable using just that slab's site variables. These events ARE
independent. Diagonal edges or horizontal periodic winding do not create shared
site variables between bands.

**Uniform bound on possible entry vertices.** Centre one unit square at each
lattice vertex; these squares tile the torus and have total area N. Squares
centred in a band of physical width sqrt(2) are contained in its enlargement
of width 2sqrt(2). Since fibres have length ell, the number of entry vertices
is at most 2sqrt(2) ell, hence at most M=4 ceil(ell). This area argument handles
arbitrary orientation, nonprimitive u in Z^2, and twisted longitudinal seam.

A crossing's endpoint separation in the theta direction is at least
ell/8-2sqrt(2)>ell/64=r, for ell>=64. Its initial lift therefore witnesses A_x
at one of the entry vertices (possibly exiting the local ball sideways first).
By a union bound, Pr(B_j)<=M a_r. Consequently, when M a_r<1,

    P_2 <= (M a_r)^k.                                           (6)

Equation (3) implies 8N/ell^2>=4sqrt(3)>2, so

    k >= 4N/ell^2.                                             (7)

The new proof depends on the shortest period and transverse area, not on an
axis-aligned slicing of a chosen HNF display. A long HNF basis vector is not
itself evidence of a large systole.

## 5. Uniform comparison of rates

For fixed subcritical p, write a_r<=C exp(-c r)=C exp(-c ell/64).
For all sufficiently large ell, uniformly over orientation and N,

    log(M a_r)<=-c ell/128,   a_r<=1/2,
    2a_r<=c/(64ell).

Combining (4),(6),(7), using -log(1-a)<=2a,

    log(P_2/P_0)
       <= k log(M a_r) - N log(1-a_r)
       <= -c N/(32ell) + 2N a_r
       <= -c N/(64ell).

This proves (2). The negative logarithmic contribution in (4) is retained;
no exponentially small probability has been approximated by zero.

For p>p_c(NN), set p*=1-p<p_c(matching). The same geometry and argument apply
to the eight-step graph. Configurationwise duality exchanges rank 0 and rank 2
and yields the opposite odds estimate.

Finite r is monotone under adding occupied sites. Both increasing events
{r>0} and {r=2} are nonconstant. A positive pivotal configuration along an
empty-to-full chain has positive product probability at every interior p, so
M'=dPr(r>0)/dp+dPr(r=2)/dp>0. With endpoint values -1 and +1, the root is unique.
The two fixed p values p_c+-epsilon trap it for all lattices with sufficiently
large ell. This proves (1), without choosing a numerical p_c.

For H, H<=P_2/P_0 below p_c and 1-H<=P_0/P_2 above p_c. This traps its compact
interior quantiles as well. No claim about the behaviour exactly at p_c is used.

## 6. What does not follow

The theorem does NOT imply concentration of the birth-time mixture F. The
axis sequence w=j, m=ceil(exp(j^2)) is a special case of (1), but #716's
full/empty-row argument still gives L(T)=>0.5 delta_0+0.5 delta_1. Root consistency
and full-law concentration remain different statements in the SAME site model.

No L^-4 rate, critical exponent, sharpness constant near p_c, numerical interval,
or proof that a particular all-width pTL eigenvalue crossing equals a finite
balance root is supplied. An inner cylinder limit, if separately known to exist,
inherits consistency by this uniform theorem; existence of that limit is not
silently assumed for arbitrary widths.

## 7. Exact finite controls

The script uses integer Lagrange reduction, verifies shortness against ambient
lattice membership on every HNF with determinant <=50, and verifies height and
period identities. That is 2,080 period bases, not 2,080 percolation productions.

On HNFs (10,3,1), (13,5,1), (5,2,3), (4,1,4), every configuration is traversed
as a physical lifted graph for NN and matching. That gives 215,040 graph/config
checks and 107,520 complementary pair checks. Slab crossings are checked by a
separate boundary-restricted reachability search in the determinant coordinate.
Tiny controls use wider raw slabs than the asymptotic choice, retaining the
explicit separation and injectivity inequalities; they are not claimed to have
ell>=64. With r<1 their arm probability is exactly p[1-(1-p)^degree]. Rational
P_0 lower, P_2 upper and slab union bounds are checked at three p values.

Larger fixed masks include true oblique Gaussian ideals (63+16i) and (32+57i),
plus a twisted rectangle. Every chosen coarse slab is explicitly vertex-disjoint.
Two distant-subcritical sign examples use a certified simple-path union bound;
a Bernoulli-inequality check avoids exponentially large integer powers. They are
not critical-point enclosures. The proof above, not an extrapolation of controls,
is the reason the conclusion covers all integer period lattices.
