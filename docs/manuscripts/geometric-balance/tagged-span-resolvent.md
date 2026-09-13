# A complete-component span law without an age cutoff

2026-09-13. Continuation of #739, after reading the span-spectrum delivery at
`7226a2c6d099535f34486eccb1bd996f7affda13`. This is a new exact finite-state
construction and its executed finite controls. It is not a proof of the
fixed-p Brownian sewing conjecture, a publication-priority claim, or an
independent certification of the previous large solver.

The main change is to stop keeping the ages of *all* active components.
One candidate component is followed, with forbidden ancestors ensuring a
unique lowest-row anchor. A second, independent transfer sums the exact
complete-component activities with two-row boundary memory. Both describe
all heights using one fixed operator; neither has a depth parameter.

## 1. The object and the unique anchor

Use independent site percolation on `C_w x Z`, with NN or matching NN+NNN
adjacency and lifted horizontal edge displacements retained. All statements
below use `0<p<1`. Empty rows occur independently with positive probability,
so at fixed finite w every component is finite, even when p is above the
corresponding *planar* threshold.

For a complete winding component C put

    L(C) = max_y(C) - min_y(C) + 1.

Anchor C at its lowest row, with the smallest column label in that row as
a tie break. Let d_h be the mean number of anchored components per row with
L=h, and nu=sum_h d_h. The component-Palm law is d_h/nu; it is not a
vertex-Palm, seam-mark-Palm, or sample-of-configurations law.

We use a lossless frontier connectivity/winding representation. The included
implementation closes it by exhaustive successor exploration, and refuses
an incomplete state graph at its explicit cap. It closes at widths 2--8 in
this delivery. The abstract result below applies to any fixed width at
which such a finite lossless representation is supplied. It does not infer
an all-width complexity bound from these counts.

### 1.1 A two-row source, not a stationary past

Draw rows -1 and 0 with their usual independent site weights. Every occupied
component of row -1 is forbidden. In row 0, any component already connected
to row -1 is forbidden. For every remaining row-zero component j, create
one candidate: colour j red, colour all other row-zero components whose
smallest column is smaller than j's forbidden, and leave the rest neutral.

The source is a counting measure over these candidates. Its total weight
can exceed one; it is not asserted to be a probability distribution.

Why is a stationary distribution of the entire past unnecessary? A selected
component with minimum row zero cannot meet any occupied site in row -1.
Any path into an earlier row must visit row -1 first. All such paths are
forbidden. Connectivity among forbidden sites below the cut therefore
cannot affect a successful candidate: erase those earlier rows and retain
only the actual independent occupation of row -1 and its forbidden marks.
We do **not** force the whole row -1 empty.

### 1.2 Three colours suffice through every future merger

Append independent rows. Connectivity and winding are updated by lifted
union-find. On a merger, red dominates neutral, and forbidden dominates
neutral. A red/forbidden merger rejects the candidate permanently.
If red leaves the frontier, accept it only when its component winds;
otherwise reject. Retirement of neutral or forbidden components is ignored.
There is no time or age stored in the state.

The forbidden labels must persist through neutral descendants. Checking
only direct contact with the original anchor, or choosing every row-zero
candidate independently without the tie-break veto, would overcount.
For example, row zero `[0,2]` can consist of two components which merge much
later into a winding component. Only the candidate with smaller column may
survive. Conversely a smaller candidate which retires before it joins red
is harmless, since it can never join again.

### 1.3 Pathwise accounting theorem

For any occupied realization, each complete winding component with minimum
row zero produces exactly one accepted candidate, and no other candidate is
accepted.

Proof. The selected component has no vertex in row -1. Among its connected
pieces in row zero, choose the one with smallest column. It never meets a
forbidden ancestor; all other pieces in this final component are neutral
until they join it. It accepts upon final retirement. Any other candidate
which joins this component eventually meets the smaller forbidden piece
and rejects. A candidate that meets row -1 cannot have minimum row zero.
At retirement the component is complete: no future path can reach a
component that has no frontier vertex. These are deterministic statements
about the actual nearest-row interaction, including matching diagonals.

The argument handles components which wind before the anchor has joined
all its later branches. Winding flags are inherited at mergers, not computed
only at the last row. Summing unused occupation variables leaves exactly
p^|C| (1-p)^|boundary C| for each accepted complete component.

## 2. One finite matrix gives the entire height law

Let alpha(p) be the two-row source row vector; its entries are polynomials
of degree at most 2w in Bernstein form. Let R(p) be the transition matrix
among continuing tagged states and b(p) the next-row successful-retirement
probability. Failed candidates have no successor. Every entry of R and b
is a sum of ordinary w-site Bernoulli row weights.

An empty next row always terminates the candidate, successfully or not.
Consequently

    R >= 0,  R 1 <= (1-delta)1,    delta=(1-p)^w>0.

The same property holds after the all-p exit-law-preserving lumping used
in the implementation. This lumping preserves distributions, not the state
path of a fixed mask word. We make no minimal-realization claim.

Write Z=(I-R)^(-1). A component with span h continues for h-1 steps after
its row-zero source and retires on the next step. The pathwise theorem gives

    d_h = alpha R^(h-1) b,                                      (1)
    D(z) = sum_(h>=1) d_h z^h = z alpha (I-zR)^(-1) b,           (2)
    nu = alpha Z b.                                            (3)

These are identities for all heights, not fits to the first few bins. In
particular the complete-component span distribution is rational at each
closed finite width. General rational absorption-time formulae are classical
phase-type theory [PH]; the contribution here is the candidate/forbidden-
ancestor realization of the *correct complete site-component Palm law*.

### 2.1 All moments and all omitted-tail moments

Unnormalised raw moments of order zero, one and two are

    m0 = alpha Z b,
    m1 = alpha Z^2 b,
    m2 = alpha (2Z^3-Z^2)b.                                    (4)

For k>=1, the factorial moment is

    sum_h (h)_k d_h = k! alpha R^(k-1) Z^(k+1) b.               (5)

Thus E L=m1/m0 and Var L=m2/m0-(m1/m0)^2. There is no cutoff bias.
For a reporting cutoff H, rather than dismissing a small tail-bin mass,
compute its contributions exactly:

    sum_(h>H) d_h   = alpha R^H Z b,
    sum_(h>H) h d_h = alpha R^H (H Z+Z^2)b,
    sum_(h>H) h^2 d_h
       = alpha R^H [H^2 Z+(2H-1)Z^2+2Z^3]b.                   (6)

A small tail probability alone never bounds its first or second moment.
Here the missing moments are returned by the same resolvent.

### 2.2 The height-window activity is now directly computable

For the earlier complete-component window sum Xi_(w,H), including its
external boundary in the infinite cylinder,

    Xi_H = sum_(h<=H)(H-h+1)d_h
         = alpha[(H+1)Z-Z^2+R^(H+1)Z^2]b.                    (7)

This reproduces Xi_H-Xi_(H-1)=sum_(h<=H)d_h. At fixed w,

    Xi_H = nu[H+1-E L] + alpha R^(H+1)Z^2 b.                  (8)

The intercept of the large-height activity therefore contains the full
mean span, including branches. Equations (2) and (7) refer to propagation
in the **vertical height**, not a derived horizontal OZ renewal. They do
not prove a w^(-1/2) prefactor as w increases.

### 2.3 Exact Palm sampling, with no burn-in or rejection

Let h=Zb. Remove states with h_i=0. Then

    R*_ij=R_ij h_j/h_i,   b*_i=b_i/h_i,
    alpha*_i=alpha_i h_i/nu

define a stochastic absorbing chain with initial mass one. Its accepted
absorption time has distribution d_h/nu. For a physical sample retain the
individual source-row masks and individual transition masks, weighting each
by the corresponding h value; do not sample only a representative lumped
mask word.

The product of these probabilities telescopes to

    Q(source rows, anchor, future rows)
      = product_of_original_Bernoulli_row_weights / nu

for every successful trajectory. The unique-anchor theorem then supplies
the complete-component Palm distribution, not a seam- or size-biased one.
The included exact-Fraction sampler records the original rows; independent
lifted-graph BFS recovers the actual component and checks its full span.
A row safety cap raises rather than returning a censored sample.

This is an ordinary Doob conditioning of a finite chain after its model
mapping has been proved. It does not assert that building the finite chain
at arbitrarily large w is cheap.

## 3. Independent sewing through direct component activities

There is a second representation which does not introduce a random
surrounding environment or forbidden colours. A finite selected component
is processed row by row. The state stores its frontier connectivity and
winding, plus the preceding row's selected-site mask. Empty interior rows
are impossible for a connected nearest-row component. If any selected
component retires while another selected frontier component remains, reject:
it can never reconnect. At final retirement require exactly one connected
component and nonzero winding.

For row masks A,B,C (previous, current, next), write S(B) for horizontal
neighbours of B, and E4(B)=B, E8(B)=B union left(B) union right(B).
The distinct boundary sites in the current row are exactly

    boundary(A,B,C) = [S(B) union E_G(A) union E_G(C)] minus B.   (9)

Finalise the weight u^|B| v^|boundary(A,B,C)| when the next mask is known.
The source additionally contributes v^|E_G(first)| for the external bottom
row. The terminal transition includes v^|E_G(last)| for the external top
row. Thus a valid path has weight u^|C| v^|boundary C|, with each vacant
boundary site counted once. Other occupied components in the random
realization need not be enumerated at all.

Let Q_w(u,v) be the resulting nonnegative matrix, c_w its source and e_w
its successful exit vector. Then the formal series identity is

    Psi_w(z,u,v)
       = sum_(C anchored,minrow=0,winding,connected)
             z^L(C) u^|C| v^|boundary C|
       = z c_w(u,v) (I-zQ_w(u,v))^(-1) e_w(u,v).              (10)

At u=p,v=1-p, Psi_w=D(z) from (2). The two representations have different
state spaces and transition rules. Q is not generally substochastic, so
the tagged-chain reset certificate must NOT be blindly applied to Q.
Its physical series converges: fixed-width empty-row bounds give an
exponential tail for L, and |C|<=w L, |boundary C|<=8|C|. This also ensures
local convergence in the occupation and boundary fugacities near a fixed
physical point. A finite reachable/co-reachable nonnegative realization
therefore has spectral radius below one there.

The coefficients of (10) give joint complete-component observables, not
only the marginal span. For example, derivatives with respect to log u
and log v give E|C| and E|boundary C|, and their log-Hessian is the
corresponding covariance matrix. Along the physical p curve,

    d(log nu)/dp = E|C|/p - E|boundary C|/(1-p).                (11)

This is an exact fixed-width score identity. Passing its derivatives to a
w->infinity limit is a separate regularity question.

The implementation gives exact all-height moments of (L,K), K=|C|, using
linear solves for the derivatives of (10). For NN p=1/4, width four:

    E L = 3.1371029435733466...,
    E K = 7.238700819343702...,
    Corr(L,K)^2 = 0.8674117242489037....

All are exact rationals in the result file. A finite positive correlation
is neither an extensive-span law nor a refutation of a proposed asymptotic
shape/occupation decoupling.

## 4. Completed finite results and numerical certificates

Tagged reachable states, followed by the all-p exit-law lumps:

| w | tagged states | NN lumps | matching lumps |
|---|---:|---:|---:|
|2|5|3|2|
|3|13|5|3|
|4|43|10|7|
|5|131|17|15|
|6|411|36|33|
|7|1275|71|68|
|8|3963|161|152|

Each recorded state has every row successor resolved. Width-eight builds
and two-parameter certificates took seconds to tens of seconds here; this
is not a like-for-like benchmark against the team's other hardware. The
structural saving is the absence of all component ages and of D_MAX.

For NN p=1/4 and matching p=1/8 the all-height results are:

| graph | w | E L | Var(L)/(E L)^2 |
|---|---:|---:|---:|
|NN|4|3.1371029435733466|0.2114153618009875|
|NN|6|3.9380183817628125|0.1570569824332094|
|NN|8|4.56334396378038|0.1271635000680440|
|matching|4|3.2422208646738944|0.1778141825223629|
|matching|6|4.060109171060249|0.1398612019890094|
|matching|8|4.711805413158092|0.1160049348309105|

These reproduce the finite-width trend in the returned note, using a new
operator without height truncation. They still do not determine a limiting
shape. In particular none of these widths is claimed asymptotic.

At width two the full all-p generating functions are compact. Write q=1-p:

    D_NN(z) = p^2 q^4 z [1+p(1+p)z]^2
              /[(1-pqz)(1-pz-p^3qz^2)],
    D_matching(z) = p^2 q^4 z [1+q(2-p)z]
              /[(1-p(2-p)z)(1-pqz)].                         (12)

At p=1/2 both total densities are 7/48, but their span means are 76/21 and
100/21. Equal complementary densities do not imply equal shape laws.
The p=1/2 matching calculation is an exact fixed-cylinder control, NOT a
planar-subcritical matching example.

### 4.1 Certification at widths five through eight

Write A=I-R. Since ||A^(-1)||_infinity<=1/delta, approximate solutions
x1~A^(-1)b, x2~A^(-1)x1, x3~A^(-1)x2 have errors bounded recursively by

    e1=||b-Ax1||_inf/delta,
    e2=(||x1-Ax2||_inf+e1)/delta,
    e3=(||x2-Ax3||_inf+e2)/delta.                            (13)

Multiply by ||alpha||_1 to enclose m0,m1,m2, using errors e1,e2,2e3+e2.
Positive-interval division then encloses means and variances. The implementation
uses a floating solve, then one correction driven by the **exact rational**
residual, storing the sum of two dyadic rationals. The final residual is
recomputed exactly for that stored vector. It is not a certificate for a
different, pre-rounding vector. Width-eight CV^2 interval widths in the two
main controls are about 2.3e-24 and 1.6e-25. This is arithmetic control on
the new finite operator, not an asymptotic error estimate.

Widths two through four and the small direct-activity solves use Fraction
Gaussian elimination throughout. Mean-scaled Laplace values and Brownian
comparison integrals in the report are separately labelled floating/high-
precision diagnostics; they are not covered by (13).

## 5. Reading the returned span calculation accurately

The mean and CV trends in `span-spectrum-diagnostic-20260913.md` are useful.
Two interpretations and implementation details must not be carried forward:

1. A complete component's span is not the sum of irreducible-piece heights.
   It is a maximum minus a minimum after placing the pieces at their actual
   transverse positions. See the separate frontier note for the deterministic
   stability bound and the precise conditional theorem.
2. A decreasing finite-width CV^2 above pi/3-1 cannot identify its limiting
   value, whether zero or pi/3-1. Widths 2--4 at p=1/2 do not establish a
   non-diffusive asymptotic regime.
3. The `light=True` branch of `span_spectrum_solve.py` at the read commit
   passes its integer-weight matrix (converted to float32) into `stationary`
   without division by p_den^w. The normal branch does divide. Hence light
   is not using a row-stochastic kernel as written. This is source inspection,
   not evidence that a particular reported run used that branch.
4. Its certificate is computed for a rounded vector a/2^k, while the displayed
   probabilities are from pi. Add the explicit evaluation discrepancy or
   report the certified candidate. A single-bin reward bound is not by itself
   a total-moment bound. The new implementation avoids these ambiguities.
5. The headline `results/geometric-consistency/span-spectrum-20260913.json`
   was not present at the read commit; its nested directory contained only
   validation tables. The new report is fully present and independently
   regenerated. This does not erase the old data or assert that they never
   existed on a team machine.

No remote job was stopped and no old result file is overwritten by this work.

## 6. Reproduction and evidence boundary

Files in this additive delivery:

- this note and `span-resolvent-frontier.md`;
- `scripts/tagged_winding_span.py` and `scripts/tagged_span_controls.py`;
- `tests/test_tagged_winding_span.py`;
- `results/geometric-consistency/tagged-span-resolvent.json`.

From a checkout with these files:

    python -m unittest discover -s tests -p 'test_tagged_winding_span.py' -v
    OPENBLAS_NUM_THREADS=1 python scripts/tagged_span_controls.py \
        --output /tmp/tagged-span-new.json

Python 3.10+; core exact construction/sampling uses the standard library.
The full report/tests additionally use NumPy, SciPy, SymPy and mpmath, listed
here rather than silently adding a project-wide dependency. Results refuse
to overwrite an existing output. No binary tables or old source snapshots
are required.

Executed: 18,754 nonempty shape masks with independent graph-potential BFS;
11,904 complete source-row words checking the actual accepted anchors;
66 exact physical-activity parameter controls; 18 full two-fugacity
coefficient comparisons; all tagged closures through width eight; direct
activity/marked-resolvent agreement through width five; eight exact Palm
trajectory controls; 17 local tests. The small complete-word and shape
counts are controls, not Monte Carlo probability estimates. No full
Matching-One CI, fresh external peer acceptance or priority certification
is claimed.

## Sources

[PH] R. S. Maier, *The algebraic construction of phase-type distributions*,
Communications in Statistics—Stochastic Models 7 (1991), 573--602,
doi:10.1080/15326349108807207. General rational absorption-time representations
are prior art; this note proves the model-specific anchor mapping.

[CIV] M. Campanino, D. Ioffe, Y. Velenik, *Fluctuation theory of connectivities
for subcritical random cluster models*, Ann. Probab. 36 (2008), 1287--1321,
arXiv:math/0610100v2. Section 1.3.3, equation (1.10), Theorem C were read in
full and rendered; the model in section 1.1 is bond random cluster. They do
not automatically settle the periodic site-component sewing used here.
