# Count the component before assigning its prefactor

2026-09-13. Continuation of the **same** geometric-balance paper (#739).

The earlier fixed-d Poisson/Gumbel note leaves the cylinder intensity
`nu_w^G(p)` as the microscopic object controlling the finite centre. This
note makes that object explicitly computable, identifies an exact
complementary-count constraint on the crossover, and separates a proved
renewal-loop calculation from the still missing SITE-cluster sewing theorem.
The current delivery does **not** certify a site prefactor exponent or a
near-critical Gumbel limit.

## 1. Results and one correction

**Actual site model, completed here.** A one-frontier transfer retains the
connectivity, horizontal lift gains, and existing winding flag of each active
component. A reward is issued only when a winding component is permanently
retired. For widths 2,3,4 this has respectively 6,14,38 reachable states for
each of G4=NN and G8=NN+NNN. Complete successor tables close exactly. Common
homogeneous row-weight/reward lumpings have 3,4,7 states. Their stationary
rewards give rational functions of p for the actual infinite-cylinder
component intensities. These are sufficient representations, not claimed
minimal positive or linear realizations.

**Prior topology, with useful consequences.** The classification already
printed in Mertens--Ziff [MZ, section II] gives, for an honest finite torus,

    W4(omega)-W8(omega^c)=r4(omega)-1,                           (1.1)

where W counts components with nonzero ambient homology. It implies
`nu4_w(p)=nu8_w(1-p)`, equality of complementary count pressures, and a
rigorous obstruction to independent two-colour Poisson counts at one common
parameter. Neither the classification nor the idea of matching topology is
claimed new.

**Explicit renewal object, completed here.** For a positive-length,
aperiodic, transversely symmetric finite-range renewal kernel, the
once-per-closed-loop coefficient with the correct horizontal translation
factor has prefactor `1/sqrt(2*pi*D*w)`. This is a direct coefficient/Laplace
calculation. It is NOT yet an identification of percolation components with
those renewal loops. That identification is the precise external task #740.

**Correction to the preceding conditional displacement formula.** Mere
C1 regularity at a does not justify an `o(1/w)` remainder after linearizing
at a displacement of size `log(w)/w`. The exact mass-inverse formula in
section 7 is valid under the previous prefactor hypothesis and a local
inverse-Lipschitz bound. Its linear version additionally needs a Taylor
remainder of order `o(1/w)` there, for instance C^{1,alpha} regularity. This
does not affect a finite-median-centred window of size `1/w`, which only
uses differentiability. A concrete C1, locally semiconcave countercontrol is
included.

## 2. A component-retirement transfer for the actual site graph

### 2.1 Boundary data

Process spatial rows from bottom to top on `C_w x Z`. Only the current row
is exposed; there is no stored first-row torus seam. Record:

* which current sites are occupied and their connected-component partition;
* in a component with no horizontal winding, the integer horizontal cut
  gain from its first current vertex to each other current vertex;
* a Boolean flag on a component once a nonzero horizontal cycle has appeared.

Within a winding component, relative gains may be replaced by zero: every
future component meeting it is already winding. This replacement is NOT
permitted in a nonwinding component. It would destroy the distinction, for
example, between row histories `[13,5,13,0]` and `[7,5,13,0]` at width four.
The first has zero horizontal-winding components, the second has one.

For the new row, add its occupied NN ring edges. Across the old/new interface
add vertical edges (G4), or vertical and both diagonal edges (G8). A cut
gain for an edge from column i to column i+dx is `floor((i+dx)/w)`. Lifted
parallel edges at width two are retained, not collapsed.

When a component has no representative on the new row, it can never meet
any later row: both graphs have vertical step at most one. Give reward 1
if that component winds and reward 0 otherwise. Then forget the old row.
Appending a deterministic empty row at the end flushes all active components.

### 2.2 Why it retains the right count for every continuation

Replace each processed nonwinding component by a gain-labelled tree between
its exposed vertices. Any future closed path can substitute an old path for
a tree path or vice versa; the difference is a previously existing closed
path. In a nonwinding component that difference has zero horizontal gain.
A winding component stays winding under arbitrary continuation, so its flag
is enough when it merges. Components never split when processed vertices are
forgotten, because the forgotten connecting paths are retained by the
partition. A completely unexposed component cannot reappear.

Induction on appended rows therefore proves: the cumulative reward plus the
number of active winding components is the number of complete winding
components of the entire finite free-height cylinder. The empty-row flush
returns the exact complete count. This is a statement about component counts,
not just the rank of their union.

The encoding works whenever its reachable-state exploration terminates. The
returned complete successor table proves termination/closure at each tested
width. It is not an all-width asymptotic state-count formula. For widths
2/3/4, both graphs give the stated 6/14/38 states, all gains in {-1,0,1}, and
1488 total row-transition entries over the six tables.

### 2.3 Row probabilities and the source-marked kernel

For row mask b, let `k(b)` be its occupation count and
`a_b(p)=p^{k(b)}(1-p)^{w-k(b)}`. If `tau(s,b)` is the new state and `R(s,b)`
the retirement reward, define

    K_z(s,t)=sum_{b:tau(s,b)=t} a_b(p) exp[z R(s,b)].            (2.1)

At z=0 the matrix is stochastic. Every state jumps to the empty state with
probability at least `(1-p)^w`; every listed state is reachable from empty.
Hence for 0<p<1 it is irreducible and aperiodic. Its stationary law pi is
unique. If `g(s)=sum_b a_b R(s,b)`, then

    nu_w^G(p)=pi*g.                                            (2.2)

The equality with the earlier bottom-anchor definition follows by counting
the same full components at their bottom versus their retirement row. In a
long window the discrepancy involves only components meeting its end rows,
of which there are at most O(w). For fixed w the empty-row regeneration
makes the mean cycle length and rewards integrable. Equivalently, use the
finite-window count identity and the stationary Markov reward theorem.

This does not require simulating an exponentially long torus. It does
require preserving gains and component identity, and its state space grows
with w.

Let `lambda(z)` be the Perron root of K_z and `psi(z)=log lambda(z)`. Then

    psi'(0)=nu,
    lim_{m->infinity} Var(W_m)/m=psi''(0).                       (2.3)

Changing the initial frontier or flushing the last components has only a
bounded boundary effect on the log generating function. Perron simplicity
near z=0 gives analyticity and the derivatives in (2.3) for these finite
chains. No claim is made that psi''(0)=nu at a fixed small width.

The implementation computes the variance with exact rational linear algebra.
Writing K_j=partial_z^j K_z|0 and 1 for the constant vector, solve

    (I-K_0)h=K_1*1-nu*1,      pi*h=0.

Then

    psi''(0)=pi*K_2*1+2*pi*K_1*h-nu^2.                         (2.4)

### 2.4 What the smaller matrices preserve

Partition refinement compares, for each old state, the exact integer counts
of row masks for every triple

    (new-row occupied count, reward, next block).

Equality is a polynomial identity in p, not a numerical test at a few p's.
The resulting 3/4/7-state kernels preserve the joint count/reward law for
any sequence of *row-wise homogeneous* probabilities.

They need NOT give the same answer configuration by configuration for an
unchanged sequence of fixed row masks. One can relabel/mix same-weight masks
inside a lump. Consequently the verification compares the unreduced transfer
with each physical configuration and the reduced transfer with the complete
occupation/count polynomial. Treating a stochastic lump as a deterministic
pathwise quotient would be a mistake.

## 3. Exact intensities, not a fitted prefactor

The NN functions for widths two and three are

    nu_2(p)= p^2(1-p)^2 (p^2+p+1)/(p^2-p+1),                  (3.1)

    nu_3(p)= p^3(1-p)^3 (p^6+p^3+2p^2+2p+1)
               /(p^6-3p^5+3p^4+p^3-p^2-p+1).                  (3.2)

The width-four rational function has a degree-27 numerator and degree-19
denominator. Its exact factored expression and all integer coefficients,
together with both matching-graph functions and the generating state tables,
are stored in `results/geometric-consistency/cylinder-winding-intensity.json`.
They are generated by the actual transfer, not supplied as fitted data.

At p=1/2:

| width | NN nu_w | count variance per vertical row |
|---|---:|---:|
| 2 | 7/48 | 343/6912 |
| 3 | 169/1984 | 4769011/122023936 |
| 4 | 323849/5576960 | 186754153229427053/6098108298338304000 |

For example the width-two asymptotic count Fano ratio is `49/144`, not one.
The finite-width count process is not exactly Poisson. A small-width control
cannot be called evidence for an exact Poisson process at that width.

Similarly, the no-horizontal-winding probability for a free strip at width
two is `(1-p^2)^m`, because a full row is necessary and sufficient for
horizontal winding there. Its decay rate is `-log(1-p^2)`, not (3.1).
At p=1/2 these are approximately 0.287682 and 0.145833 respectively. The two
objects can have the same leading rare-event rate as w grows without being
identical at fixed w. The component density is the one needed by the
preceding Poisson-window analysis.

All-p complement identities are checked symbolically. At p=1/4,1/2,3/4,
exact Fraction stationary solves independently agree with the symbolic
functions and with complementary-graph variance rates.

## 4. Complementary counts constrain the window merger

### 4.1 The finite topology already in the literature

[MZ] explicitly states that single-wrapping black and white matching clusters
occur in equal numbers, including spirals, and that a cross-wrapping cluster
is unique and occurs only when the other colour has no wrapping cluster.
On any honest torus the same elementary subsurface argument gives:

    (W4,W8) is (1,0), (0,1), or (k,k) with k>=1.               (4.1)

For completeness, in the rank-one case take a regular neighbourhood of the
occupied graph. Its essential boundary curves are parallel primitive circles.
An essential connected neighbourhood has exactly two essential boundaries;
any further holes bound discs in the torus. Cutting along all essential
boundaries yields an alternating cyclic order of black and white essential
regions. Each has two such boundaries, so their numbers agree. Contractible
regions do not affect this count. The embedded white reduction uses the same
facewise diagonal replacement as the digital-Alexander argument. Rank-two
and rank-zero cases follow from the complement-rank identity and disjointness
of independent essential curves. This proves (4.1) and (1.1).

This is a reformulation/use of the published classification, not a claim to
have discovered a new homology observable.

### 4.2 Exact intensity and pressure duality

Divide the expectation of (1.1) by m and take m->infinity at fixed width.
Closing an open cylinder by its vertical torus seam adds at most 3w edges;
each added edge changes the essential-component count by at most one. Thus
the open-cylinder and torus mean count densities have the same limit. Hence

    nu_w^4(p)=nu_w^8(1-p),             0<p<1.                   (4.2)

For real s,t the deterministic bound |W4-W8|<=1 gives

    exp(-|t|) E exp[(s+t)W4]
        <= E exp[sW4+tW8]
        <= exp(|t|) E exp[(s+t)W4].                             (4.3)

A count pressure per row exists at each fixed w: joining two free cylinders
changes the sum of their counts by at most 3w, so log moment-generating
functions are almost additive with a size-independent error. Dividing by
length and applying subadditivity to the upper/lower shifted sequences
proves existence. Boundary closure does not change the limit. Equation (4.3)
then gives

    psi_joint,w(p;s,t)=psi_4,w(p;s+t)
                      =psi_8,w(1-p;s+t).                       (4.4)

Whenever differentiating this pressure is justified, the two count variance
rates and their covariance rate coincide. The finite-state chains above
justify it for widths 2--4 and verify the rates explicitly. The extensive
count source is a different object from the bounded rank source X.

### 4.3 A hard boundary on independent-colour Poisson continuation

At the SAME occupation parameter p on the SAME labelled torus,
`(W4,W8)=(0,0)` is impossible. Two independent Poisson variables with finite
means lambda4,lambda8 put mass `exp[-lambda4-lambda8]` there. Therefore

    d_TV(Law(W4,W8),Poi(lambda4) x Poi(lambda8))
                  >= exp[-lambda4-lambda8].                    (4.5)

The TV convention is sup over events. This nonvanishing obstruction applies
at every finite size and any limit with both means bounded. It does not
contradict the earlier two-window theorem: that theorem measures BLACK at
one parameter and WHITE at a different, separated parameter. It does prove
that their independence cannot be extended to a common critical window merely
by sending d to zero inside the fixed-d theorem.

Even when each marginal admits some approximation, their joint law must
respect (4.1). In particular, shared-label coupling is not an optional
normalization detail in a crossover calculation.

## 5. A closed-renewal prefactor calculation, with its model boundary

This section identifies exactly what a successful microscopic renewal
mapping would buy. It is a theorem about the following explicit renewal
object, NOT an unproved substitution for the site component activity.

Let `a(x,y)>=0` have finite support with integer x>=1 and integer y. Assume
reflection symmetry in y, total mass below one, and positive weights at
(1,0),(1,1),(1,-1). Define

    A(z,y)=sum_{x,j} a(x,j) z^x y^j.

Let R>1 be the unique root A(R,1)=1, set kappa=log R, and normalize
`q(x,j)=a(x,j)R^x`. With expectation under q, put

    mu=E X>0,      sigma^2=E Y^2>0,      D=sigma^2/mu.

Define the closed-loop intensity with one unit of transverse length as

    L_w = w [z^w y^0] {-log(1-A(z,y))}
        = w sum_{n>=1} (1/n)
              sum_{sum x_i=w, sum y_i=0} product_i a(x_i,y_i).  (5.1)

The factor w is horizontal translation; the 1/n removes the marked renewal
cut with the usual weighted cyclic convention. The expression itself fixes
the object even for periodic words. A percolation sewing argument must show
that its own multiplicities match this convention or explicitly correct it.

**Proposition.** Under these hypotheses,

    L_w = exp(-kappa*w)/sqrt(2*pi*D*w) * (1+O(1/w)).             (5.2)

**Proof.** Take the y^0 coefficient by Fourier inversion. For each theta,

    w[z^w]{-log(1-A(z,e^{i theta}))}
        = [z^w] z A_z(z,e^{i theta})/(1-A(z,e^{i theta})).       (5.3)

Near theta=0 the implicit-function theorem gives a simple root R(theta)
near R. Symmetry makes it even and real for real theta sufficiently small.
Expanding A at the root gives

    log R(theta)=kappa+(D/2)theta^2+O(theta^4).                 (5.4)

Indeed, differentiating with respect to log z yields mu, while the second
Fourier derivative is -sigma². The pole in (5.3) has principal part
`z/(R(theta)-z)`, so its coefficient is exactly R(theta)^(-w), with leading
amplitude one. Other roots are uniformly farther away for small theta.

For theta away from zero the triangle inequality is strict at |z|=R: equality
would require identical phases at every support point. The three positive
(1,0),(1,±1) weights force both the phase of z and theta to be zero. Compactness
then supplies a uniform larger coefficient-contour radius away from that
neighbourhood. Those Fourier contributions are exponentially smaller.
Laplace integration of (5.4) yields

    (1/2pi) integral exp[-w(kappa+D theta²/2+O(theta⁴))] dtheta
       = R^(-w)/sqrt(2*pi*D*w) * (1+O(1/w)).

This proves (5.2). More general step supports require their actual lattice
span factors; those are not hidden in the constant. Finite-support is a
sufficient hypothesis here, not the expected final site-renewal class.

For a finite-state Markov renewal kernel the analogous singular object is
`-log det(I-A)`. A simple Perron crossing again isolates one logarithmic
singularity; however an actual site mapping may involve an infinite internal
state and nontrivial boundary weights. No such mapping is supplied here.

### 5.1 A completely exact control

For x=1 and Y uniform on {-1,0,1}, with a common killing weight t in (0,1),

    L_w=t^w c_w/3^w,
    c_w=sum_{k=0}^{floor(w/2)} binom(w,k) binom(w-k,k).

Here D=2/3 and (5.2) has amplitude sqrt(3)/(2sqrt(pi)). The code computes
all coefficients with integers and checks the mass-cancelling contrast below.
For w=8,16,32,64,128, the effective powers approach 1/2 from below; the last
is approximately 0.4983006988 and the normalized amplitude is 0.9985352834.
These are a renewal-model check, NOT site-percolation data.

Omitting the w factor in (5.1) would change the power from 1/2 to 3/2 without
changing kappa. A two-point decay rate alone cannot determine this counting
normalization. Degenerate transverse variance would also invalidate the
Gaussian prefactor rather than supply the same theorem with D=0.

## 6. The computation that distinguishes a power without fitting the mass

Suppose, at a fixed subcritical p,

    nu_w=A w^(-beta) exp(-kappa*w)(1+o(1)).                     (6.1)

Then

    R_w=nu_w*nu_{3w}/nu_{2w}^2 -> (4/3)^beta,                  (6.2)
    beta_eff(w)=log(R_w)/log(4/3) -> beta.

Both log A and the exponential mass cancel algebraically. No pc estimate or
amplitude fitting is used. The finite expression is always defined; its
interpretation as beta uses (6.1). Corrections and possible width arithmetic
oscillations can dominate a small-width contrast, so an agreement at three
widths is not a theorem about w->infinity.

The concrete external task #741 is exactly two graph/probability inputs,
G4 at p=1/4 and G8 at p=1/8, on w=4,8,12. Both inputs are safely subcritical
by elementary nonbacktracking path bounds. These are not a common-mass pair.
They avoid giant torus simulations and probe the proposed sewing power.

The reference Python builder was capacity-probed, not physically enumerated,
at widths 5--8. Full state counts were 102,282,786,2214; the all-mask transition
counts were 3264,18048,100608,566784. NN width eight took about 6.7 seconds
on this host and gave 90 stochastic reward blocks. No polynomial algorithm
or cheap width-24 calculation is inferred. At width 12 sparse/site-wise
factorization may be necessary.

### 6.1 A residual can be turned into an actual rare-density error bound

The row chain has the common empty-row reset `delta=(1-p)^w`. It is therefore
an L1 contraction by at most 1-delta on zero-mass signed measures. For any
normalized nonnegative candidate pi_hat, let `r=pi_hat*K-pi_hat`. Summing
its propagated residuals gives

    ||pi_hat-pi||_1 <= ||r||_1/delta.

Consequently

    |pi_hat*g-nu| <= ||g||_infinity ||r||_1/delta.               (6.3)

The implementation supplies this certificate in exact fractions; larger
floating runs must use outward error bounds including matrix-vector and
reward-evaluation rounding. The requested target is absolute error <=1e-8
on log nu where practical. Merely printing a tiny residual is not a bound on
relative error in a rare event. This is a usable numerical estimate, not a
request for an additional audit framework.

## 7. Correct centering at order 1/w

Let (6.1) hold locally uniformly near a with A continuous and positive,
constant beta, kappa(a)=d and kappa locally inverse-Lipschitz. Assume

    log m=d*w+gamma log w+c0+o(1).

For a fixed intensity level lambda>0 (lambda=-log(1-u) for a first-birth
u-quantile), the exact-mass centering is

    p_w = kappa^{-1}( d +
           [(gamma-beta)log w+c0+log A(a)-log lambda]/w )
             + o(1/w).                                        (7.1)

Here p_w denotes an intensity solution; the earlier uniform Poisson/window
argument transfers it to the corresponding true quantile at a regular
mass point. The proof first brackets p_w within O(log w/w) of a using the
positive lower slope, replaces log A(p_w) by log A(a)+o(1), and uses the
inverse-Lipschitz bound to turn the remaining o(1) logarithmic error into
an o(1/w) p error. This does not linearize kappa at the larger displacement.

If kappa is C^{1,alpha} near a, alpha>0, and v=-kappa'(a)>0, (7.1) simplifies to

    p_w=a+[(beta-gamma)log w-c0-log A(a)+log lambda]/(v*w)
                  +o(1/w).                                    (7.2)

Indeed w*(log w/w)^{1+alpha}->0. The same linear formula is valid with
mere differentiability when the logarithmic coefficient beta-gamma vanishes.
Without either that cancellation or an appropriate Taylor-remainder bound,
(7.2) has not been justified to o(1/w). This corrects that precision claim
in the preceding conditional-prefactor discussion, not its finite-median
Gumbel theorem.

**Countercontrol.** On h>0 put

    kappa(a+h)=d-h-h/log(e/h),

and on h<=0 put kappa(a+h)=d-h. On a small neighbourhood it is C1,
strictly decreasing and concave, with derivative -1 at a. Let
`nu_w(p)=w^{-1} exp[-w*kappa(p)]` and `log m=dw`. At intensity one the true
positive displacement h solves

    h+h/log(e/h)=log w/w.

Then h/(log w/w)->1, but

    w*(h-log w/w)=-w*h/log(e/h) -> -1.

Thus the linear formula's error is not o(1/w), even under the local
semiconcavity obtained in the prior argument. This is a mathematical
regularity control, not a claim that the actual site mass has this defect.

## 8. Near-critical crossover: what is now fixed and what remains open

Equation (4.5) is an unconditional constraint on simultaneous complementary
counts. It rules out copying the separated-window independent Poisson law
into a common critical window. It does not by itself find the crossover law.

The earlier Poisson proof used a compact subcritical p interval, hence fixed
localization constants C,c. If p=p_w approaches criticality, its error estimate
would require new uniform constants. Keeping them explicit, sufficient
conditions of the same form include a cutoff H_w with 4H_w<m, H_w>=w, and

    C_w m(w+1) exp(-c_w H_w) -> 0,
    m w^8 (H_w+w+1)^3 exp[-2(w-1)kappa(p_w)] -> 0.              (8.1)

These are conservative sufficient conditions inherited from the local
indicator proof, not a necessary physical crossover criterion. Neither
c_w nor their divergence scale is supplied by a fixed-p compactness argument.

There is one useful, explicitly CONDITIONAL scaling calculation. Suppose a
near-critical component-density theorem establishes

    w nu_w(p)=Psi(z)(1+o(1)),      z=w*kappa(p),
    Psi(z)~c_* sqrt(z) exp(-z),    z->infinity,                 (8.2)

with the uniformity needed for the same p_w and w limit. Then, for aspect
ratio R=m/w->infinity and intensity lambda,

    z=log R + (1/2)log log R + log(c_*/lambda)+o(1).            (8.3)

This follows by taking logarithms of R*Psi(z)=lambda. It shows why a
near-critical analysis may depend on log(m/w), not simply log m. At a bounded
aspect ratio the large-z approximation itself fails. In fixed-p notation,
(8.2) would require A(p)~c_*sqrt(kappa(p)); a fixed subcritical amplitude cannot
be held constant all the way to pc.

Neither (8.2) nor a site critical exponent is proved here. The literature
input most directly suggested by this question is [DM26], but it treats
BOND random-cluster configurations and its Theorem 1.1 gives two-sided
comparability of a TWO-POINT function, not the exact component amplitude.
The site model and the sewing/normalization remain genuine separate steps.
Do not label the heuristic `beta=1/2` or equation (8.3) as a demonstrated
near-critical square-site law.

## 9. Execution and source record

Commands from the repository root:

    python -m unittest discover -s tests -p 'test_cylinder_winding_intensity.py' -v
    python scripts/cylinder_winding_intensity.py --output /tmp/intensity-new.json

The report refuses to overwrite an existing path. Python stdlib suffices
for the transfer, Fraction stationary solves, controls and tests; report
regeneration also uses SymPy for the all-p rational functions. No existing
unmerged table or script is a dependency of this new calculation.

Executed controls: 139,776 open-cylinder graph/configuration comparisons
on both adjacencies (2x4,3x4,4x4); 66,064 complementary torus pairs
(2x2,3x3,4x4); 18 Fraction/symbolic intensity checks and their complementary
variance checks; all-p symbolic density duality at the three widths. The
physical oracle is a separate raw-(dx,dy)-potential BFS, not the transfer's
one-dimensional gain DSU. Reduced and full transfer occupation/count
polynomials agree, without a false configurationwise lumping assertion.
Twelve local mathematical tests passed. Full Matching-One CI was not run.

[MZ] S. Mertens and R. M. Ziff, *Percolation in Finite Matching Lattices*,
arXiv:1603.07289v2, section II: paragraphs on single/spiral counts and unique
cross-wrapping, before equations (6)--(11). Primary HTML read this round:
https://arxiv.org/html/1603.07289v2 . Its known classification supports section 4;
no prior-art absence is claimed for the count identity or ordinary transfer
methods.

[DM26] L. D'Alimonte and I. Manolescu, *Near-critical Ornstein--Zernike theory
for the planar random-cluster model*, arXiv:2510.13648v3, 23 June 2026.
Primary PDF printed pp3--4 (model definition and Theorem 1.1) read AND rendered:
https://arxiv.org/pdf/2510.13648 . The configuration variables are edges and
the symbol in the theorem is asymp. Current abstract/version record checked:
https://arxiv.org/abs/2510.13648 . Only these targeted sections were read,
not the entire 44-page proof. Some HTML versions displayed mismatched dates
and incomplete formulas; the cited scope is tied to the v3 PDF.

[CI/CIV leads] Campanino--Ioffe (2002), DOI 10.1214/AOP/1023481005; and
Campanino--Ioffe--Velenik, arXiv:math/0610100. Author/institution metadata or
abstracts checked in this round, not a full site-sewing theorem reading.
They are explicit starting points for #740, not imported site conclusions.

The only external jobs opened are #740 (the actual site sewing/prefactor
mapping) and #741 (the six fixed-p intensity values and one cancelling
contrast per graph). Both return to #739. No new GPU, exponential torus
simulation, or parallel mechanism programme has been launched.
