# A dilute winding crossover, before the fixed-p sewing theorem

2026-09-13. Continuation of the same geometric-balance paper (#739).
This note does not redo #740's literature investigation or #741's six stationary
computations. It gives a different, directly proved limit of the **actual NN
site model**, and a separate finite-state renewal lemma useful for interpreting
the pending sewing calculation. No fixed-p site Ornstein--Zernike prefactor,
critical exponent, new threshold value or novelty certification is asserted.

## 1. The actual-site result

On the infinite cylinder C_w x Z, independently occupy each vertex with
probability p. Count each complete component having nonzero horizontal
homology once, and let nu_w(p) be its mean count per vertical row. It is the
same density as in `winding-intensity-and-prefactor.md`, not a no-winding
survival exponent or a number of paths. In the range used below, 3p<1 already
ensures finite clusters by elementary path counting. Anchoring components at
their lowest row is one possible definition of the stationary density.

**Theorem 1 (uniform dilute crossover).** For w>=6 and 0<p<=1/8,

    max{0, 1 - 128 w p^2 - 64 w (3p)^(w-1)}
       <= nu_w(p) / [p^w I_0(2wp)]
       <= exp(6 w p^2).                                      (1)

Here

    I_0(2x) = sum_(r>=0) x^(2r)/(r!)^2
            = (1/2pi) integral_(-pi)^pi exp(2x cos t) dt.     (2)

The constants in (1) are deliberately conservative. In particular, if
w->infinity, p=p_w>0 and w*p_w^2->0, then

    nu_w(p_w) = p_w^w I_0(2w p_w) (1+o(1)).                  (3)

This includes three regimes, rather than presupposing a fixed Gaussian
amplitude:

- wp->0: nu_w/p^w -> 1, almost straight winding;
- wp->lambda in (0,infinity): nu_w/p^w -> I_0(2lambda);
- wp->infinity while wp^2->0:

    nu_w(p) = p^w exp(2wp)/sqrt(4pi wp)
              * [1+O(wp^2 + 1/(wp) + w(3p)^(w-1))].         (4)

Equation (4) describes a genuine NN-site overlap regime, for example
p=w^(-3/4). It does **not** say that the fixed-p mass equals -log p-2p:
terms of size wp^2 are precisely what stop being negligible at fixed p.
No interchange of those limits is made.

Only the standard Bessel identity and its large-real-argument expansion are
imported [B1,B2]. The component-versus-cycle bounds proving (1) are below;
there is no percolation OZ input and no surrogate renewal law in Theorem 1.

## 2. An upper bound from unrooted positive-winding walks

Every finite essential NN component in an annulus contains a simple essential
cycle. Such a cycle has winding +1 or -1, not a higher primitive winding;
orient it positively. A simple cycle with n vertices is occupied with
probability p^n. Summing over simple cycles therefore upper-bounds the
component count. Replacing simple cycles by all walks with lift displacement
(w,0), still assigning weight p^n, is a further positive upper bound.
Repeated walks are **not** assigned their actual occupancy probability here.
They are merely extra positive terms in a generating series bounding the
simple-cycle sum.

The density of rooted walks is w times the count rooted at one prescribed
site per row. Dividing by n removes the n choices of root for every simple
cycle. Counting only displacement +w removes the orientation factor two.
Thus, absolutely convergently for p<1/4,

    nu_w(p) <= U_w(p)
      := w [z^w y^0] sum_(n>=1) p^n/n
                     (z+z^(-1)+y+y^(-1))^n.                 (5)

The longitudinal density statement can equivalently be obtained by counting
in a finite-height cylinder, dividing by height, then sending the height to
infinity. Subcritical finite moments dispose of the end effects. This is
mass transport for marked cycle vertices, not a percolation independence
assumption.

Put b(t)=1-2p cos t and

    t_p(t) = 2p / [b(t)+sqrt(b(t)^2-4p^2)].                  (6)

Factoring

    1-2p cos t-p(z+z^(-1))
       = (p/t_p)(1-t_p z)(1-t_p z^(-1))

gives the exact walk-series evaluation

    U_w(p)=(1/2pi) integral_(-pi)^pi t_p(t)^w dt.             (7)

For 0<p<=1/8, the smaller quadratic root obeys 0<t_p<=2p. If
v=2p cos t+p t_p, then t_p/p=1/(1-v), |v|<=9p/4<=9/32.
Using

    |-log(1-v)-v| <= |v|^2/[2(1-|v|)],

we obtain, uniformly in t,

    |log(t_p/p)-2p cos t| <= 6p^2.                           (8)

Equations (2), (7), (8) prove the upper bound in (1). The planar walk series
is not confused with the actual site intensity: it is only an upper bound.

## 3. A lower bound which still counts components once

### 3.1 Separated directed winding cycles

Choose 2r columns on the circle of length w, any two at cyclic distance at
least three, and assign r up steps and r down steps. At each chosen column
make that single vertical step, and otherwise advance one step in the
positive horizontal direction. The path closes after w horizontal steps
and has w+2r distinct vertices. Its orientation is positive and it has no
occupied chord. Use the height of its incoming edge at column zero to
address vertical translates; this addresses each such cycle exactly once
per unit height, with no extra division by w or w+2r.

The number of admissible column sets of size k is

    S(w,0)=1,
    S(w,k)=w/(w-2k) * binom(w-2k,k),       w>=3k, k>=1,
    S(w,k)=0,                             w<3k.              (9)

This is the usual cyclic-gap count: distinguish a chosen column, subtract
two mandatory vacant gaps per chosen column, and then remove the k choices
of distinguished column. Consequently the number of cycles with 2r vertical
steps per unit height is

    N_r(w)=S(w,2r) * binom(2r,r).                            (10)

A form more convenient for a uniform bound follows without any factorial
approximation. Draw k labelled independent uniform columns. A given pair
is at cyclic distance 0,1 or 2 with probability 5/w (w>=6). The union bound
shows

    S(w,k) k! / w^k >= 1 - 5k(k-1)/(2w).

For k=2r this yields, for all r including those for which N_r=0,

    N_r(w) >= w^(2r)/(r!)^2 * (1-10r^2/w).                  (11)

A negative right side causes no difficulty.

### 3.2 Decorations are allowed; extra essential cycles are not

Requiring the entire outside boundary of each cycle to be closed would
insert an artificial factor approximately exp(-2wp), losing the crossover.
Instead, condition only on its n=w+2r sites being occupied, allow decorations,
and subtract an upper bound on the probability that another essential cycle
lies in the same complete component.

Call C good when it is the unique essential simple cycle of its component.
Then good cycles in a configuration are counted at most once per component.
Sufficient conditions for goodness are:

(i) no off-C component has two or more attachment edges to C; and
(ii) no off-C component attached to C itself has essential winding.

Condition (i) includes different attachment edges landing at the SAME cycle
vertex. Omitting that possibility would miss a second cycle attached at one
vertex. It does not forbid a tree or contractible loop attached by a single
edge.

Here is a union bound conditional on C being occupied; outside sites remain
independent Bernoulli(p).

**One-site return.** An off-C vertex adjacent to two C vertices is a missing
corner at a turn. A straight triple cannot cause one; and a return of the
height within two columns is excluded by the separation condition. There
are 4r turns, so at most 4r such sites. We use the looser bound 16r p.
The same classification follows by inspecting the isolated up-step and
down-step three-column patterns, not from numerical extrapolation.

**Longer exterior return.** Start at a C-to-outside edge (at most 4n choices)
and follow a self-avoiding off-C path with k>=2 vertices that has a further
attachment to C. Ignoring the final attachment restriction only increases
the count. There are at most 4n*3^(k-1) candidate paths, hence total
conditional probability at most

    sum_(k>=2) 4n*3^(k-1)*p^k = 12n p^2/(1-3p).

This also covers two attachments at the same C vertex when their exterior
endpoints differ. Every multiply attached off-C component has a shortest
path between two of its attachment vertices and is detected this way.

**An attached off-C winding component.** Such a component contains an
essential simple cycle of at least w vertices. From an off-C neighbour of
C, a shortest path to that cycle followed almost once around it gives a
self-avoiding off-C path with at least w vertices. Summing these paths and
using the 4n possible initial attachments is bounded by

    16n (3p)^(w-1).

This covers an essential cycle hanging from C through a stem, even though
there is only one attachment edge. These are distinct reasons for failure;
none is assumed independent of the others.

We have proved the finite lower bound

    Pr(C good | C occupied) >= 1-e_r,
    e_r=16r p + 12(w+2r)p^2/(1-3p)
                    +16(w+2r)(3p)^(w-1).                   (12)

Replacing 1-e_r by its positive part is valid. Summing good-cycle intensities
therefore gives the computable bound

    nu_w(p)/p^w >=
      sum_(0<=r<=floor(w/6)) N_r(w) p^(2r) (1-e_r)_+.        (13)

The sum in (13) is a rigorous finite lower bound, not an estimate from
independent cycle occurrences. The uniqueness condition is precisely what
makes the sum safe in the presence of overlaps.

### 3.3 Summing without assuming a bounded number of vertical steps

Write lambda=wp, A_r=lambda^(2r)/(r!)^2. Equation (2) gives

    sum A_r=I_0(2lambda),
    sum r^2 A_r=lambda^2 I_0(2lambda),
    sum r A_r<=lambda I_0(2lambda).                          (14)

The last inequality is Cauchy--Schwarz. Since N_r p^(2r)<=A_r, (11) loses at
most 10lambda^2/w times I_0. For the actual separated cycles, r<=w/6 and
w+2r<=4w/3. Equations (12)--(14) then show that the decoration subtraction
costs at most

    [16wp^2 + 16wp^2/(1-3p)
                         +(64/3)w(3p)^(w-1)] I_0(2lambda).

For p<=1/8 the sum of these losses is less than the conservative error
128wp^2+64w(3p)^(w-1) in (1). This proves the lower bound for all stated
w,p, not merely for fixed lambda. In particular the error tends to zero when
wp^2->0, even when wp grows without bound. This completes Theorem 1.

## 4. An exact finite-width series and its interpretation

For every fixed integer w>=3,

    nu_w(p)=p^w [1+w(w-3)p^2+O_w(p^3)],       p downarrow 0.  (15)

The coefficient of p^(w+1) is zero, rather than a negative perimeter term.
One proof counts the difference between complete winding components and
fully occupied rows. Fully occupied rows have mean density exactly p^w.
Components containing only one such row make zero contribution to the
difference regardless of tree decorations. Merging two full rows needs at
least 2w occupied vertices. A winding component without a full row needs at
least w+2 vertices: an essential simple NN cycle has w net horizontal steps
and must make at least two vertical steps if it leaves a row.

At size w+2 the cycle has exactly two vertical steps, on two adjacent rows,
and no negative horizontal step. Its two horizontal arc lengths are r and
w-r. Avoiding a full row requires 2<=r<=w-2, giving w(w-3) sets per unit
height. Each is an induced cycle. For w=3 the family is empty and the stated
coefficient is zero. For w=2, 2w=w+2 and full-row mergers intervene: its
separate retained rational control gives nu_2=p^2-p^4+... . It is not covered
by (15).

For fixed w the expansion is legitimate from the absolutely convergent
small-p cluster sum: the number of connected k-vertex sets rooted at a
specified site is at most exponential in k, and their external boundaries
have size at most 4k. No infinite-volume analytic continuation at pc is used.

The older exact rational functions at w=3,4 reproduce (15). An independent
physical BFS over two-row fixed-occupation sets at w=3,...,8 verifies the
minimal no-full-row counts 0,4,10,18,28,40. These checks are finite controls
of the classification, not the reason (15) holds at every fixed width.

For the mass-cancelling contrast in #741,

    R_w(p)=nu_w(p) nu_(3w)(p)/nu_(2w)(p)^2,
    beta_eff=log R_w/log(4/3),

(15) gives, at each fixed w>=3,

    log R_w(p)=2w^2p^2+O_w(p^3),
    beta_eff ->0 as p downarrow 0.                           (16)

Thus even the actual site model does not have beta_eff=1/2 uniformly in p
at a fixed finite set of widths.

A sharper result follows from Theorem 1. If w->infinity and wp_w->lambda,
using the SAME p_w at w,2w,3w,

    beta_eff(w,p_w) -> B(lambda)
      := log[I_0(2lambda) I_0(6lambda)/I_0(4lambda)^2]
                                                    /log(4/3). (17)

B(lambda)=2lambda^2/log(4/3)+O(lambda^4) near zero. At infinity,

    B(lambda)=1/2+1/[48lambda log(4/3)]+O(lambda^(-2)).        (18)

Consequently B can overshoot 1/2; it is not a monotone interpolation from
zero to one half. Direct high-precision evaluations include

| lambda | B(lambda) |
|---:|---:|
| 0.1 | 0.0653971727463187978 |
| 0.25 | 0.307660359689515652 |
| 0.5 | 0.602174810123850554 |
| 1 | 0.633453012816985342 |
| 4 | 0.520636654466422068 |

These are limiting crossover values in (17), NOT computed stationary
intensities at widths 4,8,12. In particular NN p=1/4 in #741 is outside the
explicit small-p estimate (1); no numerical prediction for that task is
manufactured from this table. The task should still return its original six
values and achieved errors.

### 4.1 Why the matching graph already differs in its leading small-p term

For each fixed width, the previously classified minimal matching cycles give

    nu_w^8(p)=c_w p^w+O_w(p^(w+1)),
    c_w=[y^0](1+y+y^(-1))^w.

At minimal size w there is one vertex in each column, each forward step has
vertical displacement -1,0,+1, and the displacements sum to zero. Addressing
the cycle by its column-zero height counts it once per vertical row. This
reuses the existing central-trinomial minimal-pattern classification; it is
not the NN dilute theorem (1) applied to diagonal edges.

At the exact fixed-width limit p downarrow 0,

    R_4^8 -> (19*73789)/1107^2 = 1401991/1225449,
    beta_eff^8(4,p) -> 0.4678291580386769565... .

The NN limit at the same three widths is instead zero. Neither is the
finite-p value assigned in #741, and no missing nu_12 was computed here.
The p and width limits must be distinguished on both graphs.

## 5. Internal-state renewals: the variance entering a sewn loop

This section is a SEPARATE explicit renewal theorem, not a claimed finite-state
representation of actual percolation irreducible clusters. It strengthens the
scalar coefficient calculation already supplied to #740.

Let A(z,y) be a finite matrix with nonnegative finite-support coefficients
a_ij(x,y), positive integer forward lengths x, and integer transverse increments
y. Suppose its spectral radius reaches one at (R,1), R>1, with a simple
Perron eigenvalue. Write right and left eigenvectors r,l, normalized l^T r=1.
The tilted edge kernel

    q_ij(x,y)=a_ij(x,y) R^x r_j/r_i

is stochastic after summing x,j. Let P be its state transition matrix and
pi_i=l_i r_i its invariant law. Assume P primitive; zero stationary transverse
drift; reflection symmetry, possibly involving a permutation of states; and
that det(I-A(R exp(it), exp(i theta))) has no zero on the unit two-torus
except (t,theta)=(0,0). These are explicit lattice-span/spectral assumptions.
A nonzero local second moment alone is not a diffusion hypothesis.

Define Q_k(i,j)=sum_(x,y) q_ij(x,y)y^k and mu=E_pi X. Solve

    (I-P)h=Q_1 1,       pi h=0.

Then the asymptotic transverse variance per renewal and per forward length are

    sigma_eff^2=pi Q_2 1+2pi Q_1 h,
    D=sigma_eff^2/mu.                                        (19)

The second term retains serial correlation carried by the internal state.
Equivalently sigma_eff^2 is the stationary mean of
(y+h_j-h_i)^2. This follows by decomposing the transverse additive functional
into a martingale plus the telescoping h boundary term. In particular it is
nonnegative; it vanishes precisely when every allowed increment is the
coboundary y=h_i-h_j. Require D>0 below.

For the explicitly defined loop object

    L_w = w[z^w y^0]{-log det(I-A(z,y))},

one obtains

    L_w=R^(-w)/sqrt(2pi D w) * (1+O(1/w)).                    (20)

**Proof.** Apply z d/dz to the logarithm. The singularity of
tr[(I-A)^(-1) z A_z] at its unique simple Perron root R(theta) has principal
part z/(R(theta)-z); the left/right eigenvector derivative cancels. The pole
coefficient is exactly one, not the number of internal states. Differentiating
the Perron equation gives

    log R(theta)=log R+(D/2)theta^2+O(theta^4).

Here the second derivative of the Perron log eigenvalue is the asymptotic
variance (19), not just E Y^2. The spectral unit-torus condition controls all
other phases exponentially. Fourier inversion and Laplace integration yield
(20), including the stated O(1/w) for finite support and reflection symmetry.

An open resolvent has a different normalization. For fixed vectors c,b,

    [z^w y^0] c^T(I-A)^(-1)b
       ~ (c^T r)(l^T b)/mu
                         * R^(-w)/sqrt(2pi D w),             (21)

when its displayed amplitude is nonzero. This is why a two-point endpoint
amplitude cannot simply be declared the cyclic component amplitude. An actual
site sewing theorem must first identify the correct object, boundary weights
and cut multiplicity. A nontrivial cyclic weight can change even the leading
amplitude. None is removed by naming the Perron root.

### 5.1 An exact two-state correlated control

Take P=[[3/4,1/4],[1/4,3/4]]. Independently choose X=1 or 2 with probability
one half each. Conditional on arrival in state + or -, choose Y=0 with
probability one half, otherwise Y=+1 or -1 respectively. Tilt back with R=2.
Thus

    A(z,y)=(z/4+z^2/8) P diag((1+y)/2,(1+y^(-1))/2).

The stationary one-step variance is 1/2. The lag-k covariance is
(1/4)(1/2)^k, so sigma_eff^2=1, mu=3/2 and D=2/3. Ignoring state correlations
would instead give D=1/3 and the wrong leading amplitude by a factor sqrt(2).

The determinant simplifies exactly:

    det(I-A(2z,y))
       =1+(2+y+y^(-1))(-6z-5z^2+2z^3+z^4)/32.              (22)

Coefficient recursion from (22) gives the whole exact loop sequence. An
independent direct matrix trace sum agrees through width ten. Near zero,

    log R(theta)=log2+theta^2/3+7theta^4/54+O(theta^6),

so in this specified model

    L_w=2^(-w)/sqrt((4pi/3)w) * [1-7/(8w)+O(w^(-2))].       (23)

Its beta_eff at width 4 has no obligation to be near one half. The exact
control approaches from below, whereas the actual dilute-site crossover
(17) can lie above one half. These are reasons to distinguish finite-width
corrections, not excuses to disregard a contradictory returned computation.

For clarity, memory can also kill diffusion altogether: with the same P,
X=1 and Y=g(i)-g(j) for two different g values, the single-step variance is
positive but every closed state cycle has transverse sum zero. Then D=0 and
L_w=R^(-w)tr(P^w), with no w^(-1/2) factor. That control violates the explicit
nondegeneracy assumptions of (20), not the theorem.

## 6. What changes in #740 and #741, without expanding either task

The fixed-p site sewing remains open in this delivery. Theorem 1 does not
supply the prefactor at NN p=1/4 or matching p=1/8; it identifies a directly
proved nonuniform dilute regime of the same NN-site model. Equation (20)
does not identify the actual site renewal state with a finite matrix.

The additions useful to the pending work are concrete:

1. A site prefactor claimed uniformly toward p=0 must reproduce (1)--(4),
   or explicitly exclude that regime. A fixed-p theorem need not be uniform.
2. A multistate sewing uses the long-run transverse variance (19), the actual
   lattice span and a specified closure, not an iid step variance and a
   borrowed endpoint amplitude.
3. Finite R_4 is always meaningful, but one triple cannot separate beta from
   corrections. If log nu=log A-kappa*w-beta log w+c1/w+O(w^(-2)), then

       beta_eff=beta+c1/[3w log(4/3)]+O(w^(-2)).             (24)

   Numerical errors in log nu at w,2w,3w bounded by e1,e2,e3 give beta error
   at most (e1+2e2+e3)/log(4/3). This is numerical error, not the correction
   in (24). With each error 1e-8 the former is about 1.3904e-7.

Do not enlarge the external width plan on the basis of this note. First obtain
its already requested values and the site-sewing result/obstruction.

## 7. Executed scope

The companion script independently traverses physical NN cylinders for:

- 43,743 fixed-occupation two-row sets at widths 3,...,8, checking the
  no-full-row minimum and first counts in (15);
- 2,443 separated cycles and 27 cyclic-gap families, checking inducedness,
  nonzero winding and short exterior contacts;
- an exact matrix-logdet coefficient recurrence against independent matrix
  trace expansion through ten forward lengths;
- retained NN rational density inputs at widths 2,3,4, and exact Taylor
  coefficients. Those rational inputs are reused from the preceding delivery,
  not presented as new stationary calculations.

The displayed integrals, Bessel ratios and Gaussian amplitude controls use
mpmath, with a higher-precision rerun. They are not interval quadrature.
The finite-sum lower bound (13) itself is computed as a Fraction; comparisons
to its Bessel normalization are displayed numerically. The proof, not the
finite checks, supplies Theorem 1. Thirteen local tests pass; full repository CI
has not run. No Monte Carlo, no new stationary width-eight/twelve run and no
new external issue are required by this continuation.

Commands:

    python -m unittest discover -s tests -p 'test_winding_dilute_crossover.py' -v
    python scripts/winding_dilute_crossover.py --output /tmp/dilute-new.json

The report refuses an existing output path. The script has no dependency on
an unmerged state table; mpmath is needed for numerical report generation.

## Sources and originality boundary

[B1] NIST Digital Library of Mathematical Functions, 10.32.1:
https://dlmf.nist.gov/10.32.E1 . The I_0 integral, not a percolation theorem.
[B2] NIST DLMF 10.40.1:
https://dlmf.nist.gov/10.40.E1 . Large positive argument expansion of I_nu.
Both primary reference pages were read this delivery.

[CIV] Campanino--Ioffe--Velenik, *Fluctuation theory of connectivities for
subcritical random cluster models*, Ann. Probab. 36 (2008), 1287--1321;
https://www.unige.ch/math/folks/velenik/papers/abs_CIV08.html . Author abstract
read this delivery, not a re-audit of the theorem body. It establishes the
context of irreducible chains/effective walks, not the SITE cylinder sewing
requested in #740. Matrix-additive variance and Gaussian coefficient methods
are standard; the explicit calculations here are not a claim to invent them.
No systematic prior-art search for Theorem 1 has been completed in this round.
