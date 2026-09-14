# Giant white components: bulk filling, a dual tail mass, and a Laplace response

2026-09-14. Mathematical continuation of #764/#739. This note develops the
previous gap-law result, rather than treating an eventual publication as a
prerequisite. Statements described as proved below are the arguments supplied
here; the finite calculations check their exact interfaces, not their
large-width limits.

## 1. The model, conventions and the results

On `(Z/wZ) x Z`, black sites have probability `p`, NN connectivity, and white
sites have probability `q=1-p`, matching (eight-neighbour) connectivity. Fix
`0<p<pc_NN` and let `w` increase. A component is essential when it has a
nonzero horizontal winding cycle. Every complete component is finite at
fixed width: monochromatic rows of either colour occur in both directions.
Use **one complete component once**, anchored at its lowest row. The common
black/white essential-component density per vertical row is `nu=nu_w(p)`.

For a white Palm component write

    L = max y - min y + 1,
    K = number of occupied WHITE sites,
    B = number of DISTINCT external BLACK matching-neighbour sites.

`B` counts sites, not contact edges. A black site adjacent to two distinct
white essential components contributes once to EACH component's B. Let

    theta(q) = Pr_plane(0 is in an infinite white matching component),
    beta(q)  = Pr_plane(0 is black and adjacent to an infinite white component).

Single-site insertion gives `beta(q)=(p/q)theta(q)`; see Section 5.

The new conclusions are:

1. A giant white component is asymptotically filled with deterministic bulk
   densities, though its total length remains random:

       (nu L, nu K/w, nu B/w) -> (E, theta E, beta E), E~Exp(1).       (1.1)

   All fixed mixed positive moments converge. More strongly, its longitudinal
   occupation and boundary measures converge to constant-density intervals.

2. The ultimate fixed-width WHITE span-tail mass is controlled by the same
   BLACK barrier density:

       gamma_white(w,q)/nu_black(w,p) -> 1.                         (1.2)

   This requires a gluing argument; it is not inferred from convergence of
   moments of nu L.

3. For the physical complete-cluster score

       S = K/q - B/p,       c = theta/(p q^2),

   there is a normal-exponential mixture law:

       (nu L, sqrt(nu/w) S) -> (E, sqrt(c E) Z),                    (1.3)

   where Z is standard normal, independent of E. The second marginal is
   Laplace with variance c, not normal. This follows from physical parameter
   tilting and the bulk law, without assuming a separate bulk CLT.

4. The longitudinal connectivity structure factor has a Lorentzian limit:

       (nu/w) sum_(x,y) exp(i k nu y)
           Pr((0,0) and (x,y) belong to the same white component)
           -> 2 theta^2/(1+k^2).                                  (1.4)

   This is a **connectivity** correlator, not the covariance of independent
   site labels. No continuum field identification is made.

5. On a torus with `m nu->t`, the volumes of macroscopic white components
   are theta times the Poisson circular-gap fractions. This gives explicit
   fragmentation and largest-component predictions, including the white
   cross when there are no black barriers.

All statements fix a strictly subcritical black p. No assertion is made
uniformly up to criticality. The arguments are locally uniform on compact
subintervals of `(0,pc_NN)`, which is needed for the small physical tilts in
Section 6. `q=3/4` is the intentional supercritical white control, NOT the
previous subcritical-white `q=1/8` experiment.

## 2. Reused inputs, explicitly separated from this continuation

The preceding `black-white-gap-law.md`, Sections 2--8, supplies:

- the alternating complete-component chain and common density nu;
- row projection counts `N_W(y)<=1+N_B(y)` and
  `E N_B(0)=epsilon_w=nu E_B L`;
- black root-volume exponential tails uniform in cylinder width, via the
  first-query discovery-tree lift and the plane site tail [AV];
- `a^w<=nu<=Cw exp(-cw)` and black component-Palm moments `E_B K^j=O(w^j)`;
- the black-anchor Poisson limit on the nu scale, the white gap law
  `nu L_W->Exp(1)`, and uniform exponential tails on that scale;
- the circular Poisson-gap law when `m nu->t`.

These are dependencies, not new independent theorems established by rerunning
old tables. The relevant old note is included unchanged under source_inputs.
The first-query lift preserves a discovery TREE, vertex count and vertical
reach; it does not lift winding cycles as cycles.

The estimates in that note are uniform on a fixed compact p interval: use
the black root-volume bound at its largest p, a uniform isolated-ring lower
bound on nu, and the upper bound `Cw exp(-cw)`. With localization H=w^2, its
Poisson error bounds are uniform because polynomial powers of H times nu
vanish, while `exp(-cH)` beats every fixed power of `1/nu`. Its independent
block argument then gives uniform exponential tails for the rescaled gaps.
The geometric gap/span sandwich transfers those tails to white spans.

The only additional external tools used here are matching boundary
connectivity [MZ], site Harris and the site volume tail [AV], and bounded-variable exponential concentration. For an independent X in [0,1],
the second derivative of log E exp(tX) is a tilted variance, at most 1/4.
Integrating twice proves E exp[t(X-E X)]<=exp(t^2/8); independence and
Chernoff give the inequality used below. Finite-range rows are handled by
residue classes. No additional asymptotic input is hidden in this step.

## 3. A local cage lemma: the bulk white density approaches theta

Let f_r(v) be the indicator that v is white and has a white matching path,
using only the square `v+[-r,r]^2`, to its vertex boundary. For `w>2r+2`,
this square injects in the cylinder, so

    E f_r = theta_r(q),                 theta_r(q) down to theta(q).

Let eta_w(v) indicate that v belongs to an essential white component. Then
`eta_w<=f_r`: a component entirely inside this injecting square cannot wind.

### Lemma 3.1 (nonessential white components are caged by a large black component)

There are fixed geometric constants A0,A1 such that a finite nonessential
white component C containing v and reaching distance r has an adjacent
black component D with

    |D| >= (r-A0)/A1,
    some z in D has dist_cylinder(z,v) <= A1(|D|+1).                 (3.1)

Also `|C|<=A1(|D|+1)^2` after increasing A1. The same conclusion holds for a
finite white component in the plane.

**Geometric proof.** Use a fixed subdivision of the square faces for the
4/8 complementary neighbourhood construction. At a diagonal white contact,
join the white corners; black arcs do not cross this connection. The boundary
of a neighbourhood is then a finite collection of disjoint polygonal curves.
Its black-side incident vertices along each curve are NN-connected: moving
along the curve either stays at one black vertex or crosses an NN black
edge. This is the matching-boundary observation of [MZ, Section II], with a
fixed local subdivision rather than an arbitrary chosen white path.

Since C has no essential cycle, its neighbourhood lies in a disc in the
annulus. The outer boundary is contractible and encloses C; extra boundary
curves surround holes and are not used. Lift this outer curve to a closed
plane curve. All black-side vertices project into one complete black
component D. A site of D has a bounded number of incident face sectors,
and each sector contributes a bounded number of unit-length segments of
this PARTICULAR boundary. Thus its length is at most `A1 |D|`, independently
of the number of white vertices and independently of w. One may use a loose
constant such as 128 after fixing the quarter-cell construction; the argument
uses only the existence of a universal constant.

A curve enclosing both v and a point reached at distance r has diameter at
least r-O(1), so its length forces the first bound. Since it encloses v, a
black incident vertex is within its length plus a local constant of v.
The enclosed area is at most a constant times its squared length, giving
the bound on |C|. The outer curve is contractible even when the FULL black
component D also has an essential cycle somewhere else. This distinction
is why we bound by the size of D, not by a presupposed simple black loop.

### Consequences of the cage

The uniform black root-volume tail and a dyadic union bound over possible
vertices z within `A1(2k+1)` of v give

    Pr(f_r=1, eta_w=0) <= C exp(-c r),      1<=r<w/8,               (3.2)

with constants depending on a compact subcritical black p interval, but
not on w. The number of possible z at a dyadic scale k is O(k^2), and
`Pr(|D(z)|>=k)<=C exp(-ck)`, so the union is summable. Apply the same argument
in the plane to finite white components to get

    0<=theta_r-theta<=C exp(-cr).

Since `theta_w:=E eta_w=(nu/w)E_W K`, (3.2) gives

    |theta_w-theta| <= C exp(-c w).                                (3.3)

The same proof implies a stretched-exponential **volume** tail for a
nonessential white root component, bounded by `C exp(-c sqrt(n))`, uniformly
in w. We do NOT claim an exponential volume tail for supercritical white
finite components. In particular all their root moments are uniformly
bounded. It also proves continuity of theta on the present compact q
intervals, by uniform approximation by the finite polynomials theta_r.

## 4. Filling a random, exponentially long interval without dividing by a rare event

Let `I(C)=[min y(C),max y(C)]`. For a white essential component define

    F_r(C)=sum_(y in I(C)) sum_x f_r(x,y).

Every site in C is included, so `F_r(C)-K(C)>=0`. Campbell counting and
`N_W<=1+N_B` give the exact expectation estimate

    (nu/w) E_W[F_r(C)-K(C)]
      = E[N_W(0) * (1/w)sum_x f_r(x,0)] - theta_w
      <= theta_r-theta_w+epsilon_w.                              (4.1)

This removes the difficult rare-conditioning factor from the **defect**.
Bounding an unconditional error and then blindly dividing it by nu would
not have done so.

The row field `X_y=(1/w)sum_x f_r(x,y)` lies in [0,1] and is dependent only
within 2r rows. Its mean is theta_r. For each deterministic n, split the
sum into 2r+1 independent residue classes and apply Hoeffding. Taking a union
over `0<=n<=M/nu` gives, for fixed M,delta>0,

    Pr(max_(n<=M/nu) |sum_(y=0)^(n-1)(X_y-theta_r)| > delta/nu)
      <= C(M/nu)(2r+1) exp[-c delta^2/(M nu(2r+1))].               (4.2)

Choosing `r=floor(w/8)-2` is allowed for large w. Even AFTER division by the
white anchor probability nu, the right side tends to zero. Thus the estimate
holds under white component Palm, despite the dependence of the component
length on the same labels. Use the uniform exponential tail of nu L to
remove the restriction L<=M/nu. Abel summation extends the partial-sum
estimate to bounded Lipschitz weights along the interval.

Combining (4.1), (3.3), epsilon_w->0 and (4.2) proves

    nu(K/w-theta L) -> 0 in L^1.                                 (4.3)

Since K<=wL and nu L has uniform exponential tails, this also gives every
fixed L^j convergence after using a higher uniform moment. More generally,
after anchoring the component at its lowest row,

    M_w := (nu/w) sum_(x,y in C) delta_(nu y)
       -> theta 1_[0,E](s) ds,      jointly with nu L->E.          (4.4)

The convergence is of finite measures, with their total-mass moments.
This is a bulk law for randomly delimited components, NOT a claim of
independent occupancy conditional on the component.

## 5. The boundary fills too, and the leading thermal score cancels

A complete component has physical activity `q^K p^B`. At fixed w the finite
transfer representation (or direct convergent component sum) gives

    d_q log nu = E_W S,
    d_q^2 log nu = Var_W(S) - E_W V,
    S=K/q-B/p,                V=K/q^2+B/p^2.                       (5.1)

By exact colour-density duality, `nu_white(q)=nu_black(p)`.
Black component moments and `B_black<=4K_black` give, locally uniformly,

    |d_p log nu_black|=O(w),      |d_p^2 log nu_black|=O(w^2).       (5.2)

Set `b_w=(nu/w)E_W B`. Equations (5.1)--(5.2) give

    theta_w/q - b_w/p = -(1/w)d_p nu_black,
    b_w = (p/q)theta_w + (p/w)d_p nu_black
        -> beta=(p/q)theta.                                      (5.3)

The coefficient beta has an independent plane meaning. Freeze the origin
black and ask whether a white neighbour belongs to an infinite white
component in the graph with the origin removed. Turning the origin white
makes it infinite exactly on that event: joining finitely many finite
components cannot make an infinite one. Independence of the origin label
therefore gives `theta=q Pr(A)` and `beta=p Pr(A)`.

For completeness, define `g_r(v)=1{v black}1{A_r(v)}`, where A_r asks for
such a neighbour path to the r-square boundary without using v. Then
`E g_r=(p/q)theta_r`. Every boundary site of an essential white C satisfies
g_r=1, and lies in the extended interval `I^+(C)=[min y-1,max y+1]`.
The expanded interval count exceeds N_W by at most an anchor and an upper
endpoint incidence; their expected total is 2nu. Consequently

    (nu/w)E_W[sum_(v in I^+(C))g_r(v)-B(C)]
       <= (p/q)theta_r - b_w + epsilon_w+2nu.                     (5.4)

Repeat the finite-range row argument. This proves

    nu(B/w-beta L)->0 in every fixed L^j,                        (5.5)

and the analogous boundary-measure limit. Equations (4.3),(5.5) prove (1.1).
The boundary and occupation leading vectors are collinear. In particular,
the apparently huge terms K/q and B/p cancel at order w/nu.

This does NOT make the remaining score zero: its mean is O(w), and its
fluctuations have the smaller but still diverging scale sqrt(w/nu).

## 6. A physical-tilt proof of the Gaussian--Laplace mixture

Write `s_w=sqrt(nu/w)` and `c=theta/(p q^2)>0`.
From (5.1)--(5.5),

    E(s_w S)->0,
    Var(s_w S) = (nu/w)E V + O(nu w) -> c.                        (6.1)

A variance limit alone is not a distribution theorem. To identify the law,
use an actual change of site probability, `q_h=q+h`, `p_h=p-h`.
For every complete cluster,

    log[(q_h/q)^K(p_h/p)^B]
      = h S - h^2 V/2 + O(|h|^3(K+B)),                           (6.2)

uniformly for h in a small fixed neighbourhood of zero. Take `h=t s_w`,
with real t near zero. Exact component reweighting yields

    E_q exp[-u nu L+t s_w S]
      = (nu(q_h)/nu(q))
        E_(q_h) exp[-u nu L+(t^2/2)s_w^2 V
                     + O(|t|^3 s_w^3(K+B))].                    (6.3)

The ratio of densities tends to one: by (5.2), its log is O(w|h|), and
`w s_w=sqrt(w nu)->0`. Under q_h, which approaches q, the locally uniform
bulk theorem gives

    (nu L,s_w^2 V)->(E,cE).

The cubic remainder is bounded by `C |t|^3 s_w nu L`, since B<=8K and
K<=wL. It vanishes, including under the exponential weight in (6.3).
For |t| and |u| sufficiently small, uniform exponential span tails provide
domination: the positive exponent is at most `C t^2 nu L+o(nu L)`.
Thus the joint transform converges on a neighbourhood of the origin to

    E exp[-u E+t Y] = 1/(1+u-c t^2/2).                            (6.4)

This identifies `(E,Y)=(E,sqrt(cE)Z)` with E exponential and Z an independent
standard normal. Exponential integrability near the origin also gives all
fixed joint moments. In particular,

    density_Y(y) = exp(-sqrt(2/c)|y|)/sqrt(2c),
    E Y^2=c,       E Y^4=6c^2.                                  (6.5)

This is an UNCONDITIONAL Laplace law for the score. The joint limiting
conditional law given E=e is N(0,ce). No local conditioning theorem for
every exact discrete L is asserted.

The key distinction is physical: the Gaussian fluctuations are mixed over
the exponentially distributed amount of bulk. This score projection can
be identified without first finding the entire covariance matrix of the
occupation/boundary residuals.

### The still-open, richer bulk fluctuation question

Put `rho_K,w=E K/(w E L)` and `rho_B,w=E B/(w E L)`. A natural next
conjecture is a joint stable Gaussian limit for

    sqrt(nu/w) (K-w rho_K,w L, B-w rho_B,w L),

with conditional covariance E Sigma(q). The score result constrains

    (1/q,-1/p) Sigma(q) (1/q,-1/p)^T = theta/(p q^2).               (6.6)

It does not determine Sigma. The finite-width centring is important:
using theta in its place requires an error smaller than sqrt(nu/w), which
is not supplied merely by theta_w->theta. The new task asks for this richer
limit, not a repeat of the first-order volume law.

## 7. The ultimate white pole is the black barrier density

Let `c_w(n)` be the probability of a white matching path crossing the cylinder
slab of n rows, from the first to the last row. Let `r_w(H)` be the probability
of a white horizontal essential cycle in H rows.

By the cylinder 4/8 crossing duality,

    c_w(n)=Pr(no black NN essential cycle in those n rows).        (7.1)

The prior localized black-anchor Poisson theorem, with H0=w^2, implies for
each fixed T>0

    c_w(floor(T/nu))->exp(-T).                                   (7.2)

To check the boundary issue: any slab cycle belongs either to a complete
black component of span <=H0 with anchor in an H0 enlargement, or to a black
component of size >H0 meeting the slab. The latter probability is at most
`Cwn exp(-cH0)`. Conversely, a short component anchored H0 away from the upper
edge has its cycle inside. Anchor enlargements change scaled length by o(1).

The crossing probabilities are submultiplicative on disjoint row blocks.
Therefore their height decay rate gamma exists and, for B>H>=1,

    -log c_w(B)/B <= gamma.                                     (7.3)

For the other direction, overlap successive B-row windows in H rows and
require a white horizontal cycle in each overlap. Both vertical white paths
must meet that cycle, so they join. In the matching graph a geometric
crossing of two white diagonals also joins their endpoints; it is not a
crossing of disjoint components. All the events are increasing in WHITE
sites. Harris on the overlapping supports gives

    c_w(k(B-H)+H) >= c_w(B)^k r_w(H)^(k-1),
    gamma <= [-log c_w(B)-log r_w(H)]/(B-H).                       (7.4)

This is not an independence assertion for overlapping windows.
The complementary event to r_w(H) is a black NN vertical crossing.
The uniform black reach tail hence gives

    1-r_w(H)<=Cw exp[-c(H-1)].                                   (7.5)

Use `B=floor(T/nu)`, `H=w^2`. Equations (7.2)--(7.5), `Hnu->0` and
`r_w(H)->1` squeeze gamma/nu to 1.

Finally gamma is the white complete-Palm span-tail rate. The anchored event
`L>=n` supplies a crossing of the first n rows. Conversely plant a full white
row below a crossing, and a full black row below that. The crossing joins a
white essential component with prescribed lowest row, at probability cost
`q^w p^w`, independent of its height. These two comparisons identify the
height exponents. Thus (1.2) is proved.

This removes the previous rare-mixture obstruction by a model-specific
uniform gluing bound, not by arguing that the first few moments determine
the tail.

## 8. Connectivity, form factors, and susceptibility

Let the anchored white component have row occupations k_y. Its exact
component form factor is

    A_C(k) = sum_y k_y exp(i k nu y).

The volume-measure convergence (4.4) and second moments give

    (nu/w)^2 E_W |A_C(k)|^2
      -> theta^2 E |integral_0^E exp(i k s) ds|^2
       = 2 theta^2/(1+k^2).                                    (8.1)

Mass transport identifies the left side with the scaled essential white
connectivity sum in (1.4). The nonessential white contribution is negligible
by the uniform root-volume moments after Lemma 3.1. At k=0 this yields

    chi_white ~ 2 w theta^2/nu.                                (8.2)

Equivalently, the rescaled, transverse-averaged connectivity measures tend
to `theta^2 exp(-|s|) ds`. This is a weak/integrated statement. A separate
pointwise microscopic-endpoint theorem is not silently included.

The exact rational numerical control uses `z=(1+s nu)^(-1)` and

    T_w(s) = nu^2/(2w^2 theta_w^2)
               E_W sum_(a,b) k_a k_b z^|a-b| -> 1/(1+s).         (8.3)

No transcendental evaluation is required for this finite computation.
For a direct activity transfer with matrices R_j weighting current-row
occupation k^j and exit vectors b_j, first solve

    y=(I-R_0)^(-1)b_0,
    f=(I-zR_0)^(-1)(b_1+R_1 y),
    h=(I-R_0)^(-1)(b_2+R_2 y+2z R_1 f).

Then `alpha h` is the unnormalised pair activity. The independent physical
shape enumeration checks this recursion at several finite heights; the
resolvent itself supplies all heights.

## 9. Finite-torus volume fragmentation

Reuse the earlier circular-gap result at `m nu->t>0`. Conditionally on
J=j>=1 black barriers, the j gap fractions have Dirichlet(1,...,1) law,
after a uniformly chosen barrier fixes the circular ordering. J converges
to Poisson(t). At J=0 there is one white cross, not zero white clusters.

The same local-cage argument and row concentration hold on the torus for
r<w/8. The entire volume is O(w/nu); the exponentially accurate local row
concentration dominates the number of possible macroscopic interval
endpoints. The black boundary zones have vanishing area fraction. Thus,
conditionally on each fixed j, the macroscopic WHITE component volume
fractions of all N=wm sites converge to

    theta (D_1,...,D_j),       (D_1,...,D_j)~Dirichlet(1,...,1),

or to the single fraction theta when j=0. Nonessential white components have
no macroscopic member, by their cage tail and a union bound. Their combined
fraction is q-theta; it is not generally zero.

Let P2vol denote the sum of squared macroscopic volume fractions, divided
by theta^2. The limiting expectation is

    E P2vol = 2(1-exp(-t))/t - exp(-t).                           (9.1)

Let M be the largest white component fraction divided by theta. Then

    Pr(M=1) = (1+t) exp(-t).                                    (9.2)

Both zero and one black barrier produce one macroscopic white component;
its topology differs in the two cases. For j>=2,

    Pr(M<=x | J=j)
      = sum_(l=0)^j (-1)^l binom(j,l)(1-lx)_+^(j-1),             (9.3)
    E[M | J=j] = H_j/j.

At J=0 or 1, M=1. The mixture is fully specified. These are derived limit
predictions, not newly simulated torus volumes. They illustrate that a
locally supercritical white model can have a random finite number of
macroscopic components on an exponentially elongated sequence.

## 10. Finite results and execution scope

The physical control uses black p=1/4, white q=3/4. Exact all-height results:

| w | theta_w | T_w(0) | T_w(1) | T_w(2) | T_w(4) |
|---|---:|---:|---:|---:|---:|
| 2 | .748890532544 | .944533333333 | .510668185243 | .356668711606 | .229827452010 |
| 3 | .749847650131 | .980867697451 | .501548687782 | .338740246690 | .207375743066 |
| 4 | .749969685902 | .993248030686 | .500162917714 | .334764577282 | .202069892807 |
| limit | theta | 1 | 1/2 | 1/3 | 1/5 |

The score variances scaled by nu/w at these widths are
`5.06833994408, 5.30253981269, 5.33091713810`.
The corresponding finite-density proxies `theta_w/(p q^2)` are
`5.32544378698, 5.33224995649, 5.33311776642`.
The table does not determine theta to the printed precision.
The mean scores are nonzero; they are retained in the JSON.

An independent no-black-ring transfer gives certified positive Perron
intervals. The derived gamma_white/nu_black values at widths 2,3,4,6 are
`1.13642432987, 1.03689306046, 1.01175415938, 1.00131735895`.
Every stored Perron enclosure is obtained from exact rational Collatz
quotients; displayed logarithms are high-precision approximations to exact
-log interval endpoints. The finite gluing inequalities are separately
checked by rational powers against slab probabilities.

Executed checks: 73,984 annulus configurations, 147,968 complementary
crossing/ring dichotomies; 250 complete nonessential cage controls;
2,189 connected winding shapes and 36 independent mass/pair identities;
exact dual first/second derivative identities; and 2,035 cyclic point
patterns for the discrete spacing participation identity. No new Monte
Carlo, external machine, or full repository test suite is run.

The engine is an unchanged input, Git blob
`52f3611990ce2b1331d9e5296e0262f5e402e0d7`. Numerical values of old moments
are controls, not independent new data. New outputs include pair transforms,
boundary/score covariance, finite gluing certificates and Poisson-volume
predictions. The broader two-dimensional bulk CLT remains a conjecture.

## References and reading scope

[MZ] Mertens--Ziff, *Percolation in Finite Matching Lattices*,
arXiv:1603.07289v2, Section II, surrounding-component connectivity and
single/spiral/cross classification. Relevant HTML body read this round.
https://arxiv.org/html/1603.07289v2

[AV] Antunovic--Veselic, *Sharpness of the phase transition and exponential
decay of the subcritical cluster size for percolation on quasi-transitive
graphs*, arXiv:0707.1089v3, Theorems 2--3 and Section 3. Site applicability
is stated explicitly. Relevant definitions and theorem statements read.
https://arxiv.org/html/0707.1089v3

Internal input: `black-white-gap-law.md` (2026-09-14), Sections 2--8.
The present all-width arguments explicitly depend on that note's topology
and rare-gap estimates. The local-cage, bulk filling, physical-tilt mixture,
finite gluing rate, and form-factor arguments are supplied here.

## Concurrent team connection

During this analysis PR #771 arrived; its read head was
`53b4ec6111f5f97e11cbcdf4606aeaed17d0a9ae`. Its
`structural-consequences-20260914.md`, Section 6, independently records the
complementary score identity; its
`supercritical-white-slab-bulk-20260914.md`, Section 1, gives the same exact
one-site identity beta=(p/q)theta and Sections 3--5 formulate the bulk and
vector-CLT questions. These coincident identities are NOT counted as a second
new mathematical discovery here. The present extension is the cage/Campbell
proof of the bulk law, physical-tilt score distribution, gluing tail theorem
and connectivity/fragmentation consequences. The richer vector question is
now coordinated at #772; no shared branch was overwritten.
