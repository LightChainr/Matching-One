# Predictive state depends on the generator and forecast horizon

**Research assessment and exact controls, 2026-08-30.** Advances #400 and #403;
supplies an exact priority-integral control and conditional arithmetic for #405.

## Decision

Do not infer additional continuum fields merely because an additive spatial
Hankel block or a temporal covariance block acquires rank. First specify the
operation being represented and compare against its simplest physically
appropriate control. For the birth process, go beyond the already available
one-step H2 ceiling by measuring cooperative present-state triggers. For direct
birth corrections, count the newly distinguished carrier types instead of
fitting another unconstrained exponent.

These are three targeted discriminators within existing research tracks, not
three new explanations of the same residual. No existing frozen score is changed.

## What was reread

The initial main snapshot was `6f711e524d9b33693703109521e592ca893e9cb0`.
The navigation documents do not by themselves contain all of the active-branch
science. This assessment also used the following explicit result/protocol records:

| Track | Read record | Consequence for this contribution |
|---|---|---|
| P250 endpoint relations | [rank-five extension](https://github.com/LightChainr/Matching-One/commit/2646e8f), [rank-eight R2 bridge](https://github.com/LightChainr/Matching-One/commit/93eaab1), [five fixed maps](https://github.com/LightChainr/Matching-One/commit/df57b69) | Their exclusions stand; calibrate the translation-generator class before enlarging it. |
| P334 birth age | [controlled age result](https://github.com/LightChainr/Matching-One/commit/742a8b0), [two-time kernel](https://github.com/LightChainr/Matching-One/commit/5a7f2d9) | Neither observed age dependence nor covariance rank identifies a full-state memory field. |
| P334 current configuration | [current-k0 pilot freeze](https://github.com/LightChainr/Matching-One/blob/analysis/p334-current-k0-geometry-pilot-20260830/notes/p334-current-k0-geometry-pilot-freeze.md), read blob `814f90f6b7e0d874c7c6e3718166807ef6803575` | H2 is already in the pilot. This work adds a two-step obstruction, not that one-step identity again. |
| P337 direct births | [four-generation result](https://github.com/LightChainr/Matching-One/commit/b887ef3), [carrier theorem](https://github.com/LightChainr/Matching-One/commit/1c93bed), [typed gluing](https://github.com/LightChainr/Matching-One/commit/cfb3ead) | Split theta and figure-eight channels; ordinary arm events alone are not enough. |
| P333 scalar marks | [minimal multimark obstruction](https://github.com/LightChainr/Matching-One/commit/e7e6c80) | Do not suggest another scalar fugacity/block-count mark as though the current failure were untested. |
| Jordan identifiability | [PR #385](https://github.com/LightChainr/Matching-One/pull/385) | The previous diagonalizable-to-Jordan closure suggestion is already represented; it is not a new result here. |

No active-branch production data were downloaded or rescored by this contribution.
The result JSON contains independent exact controls only.

## 1. One scaling function, unbounded translation Hankel rank

For `a,beta>0`, define

\[
g_n=(n+a)^{-\beta},\qquad H_d[i,j]=g_{i+j},\quad 0\le i,j\le d.
\]

The elementary integral representation is

\[
g_n=\frac1{\Gamma(\beta)}\int_0^1t^{n+a-1}(-\log t)^{\beta-1}\,dt.
\]

For a nonzero polynomial `P(t)=sum_i v_i t^i`, the associated quadratic form is
an integral of `|P(t)|^2` against a strictly positive density. Hence `H_d` is
positive definite and has rank `d+1` for every d. This proves the obstruction
for the control; it does not impose positivity on the measured charged data.

The rational special case `g_n=1/(n+1)` gives Hilbert matrices. Exact elimination
and the independent Cauchy determinant product agree through matrix size ten.
Nevertheless the sequence satisfies the first-order variable-coefficient rule

\[
(n+2)g_{n+1}-(n+1)g_n=0.
\]

A rank-one constant multiplier fitted at n=0,1 first misses at n=2 by exactly
`1/12`. The failure is structural, not finite precision.

For the *same* function `g(r)=1/r`, geometric sampling `r_n=2^n` gives rank one.
The log control `g(r)=(1+log_2 r)/r` gives rank two and obeys

\[
g_{n+2}-g_{n+1}+\tfrac14g_n=0.
\]

In continuous coordinates the distinction is even sharper: `r^-beta` obeys
`(r d/dr+beta)g=0`, and its single-log extension obeys the square of that Euler
operator. Low differential/dilation complexity does not imply low rank under
constant additive shifts. The sampled-value and finite-Hankel results in [R1]
provide the appropriate algebraic setting; the newer completion equivalence
[R2] does not establish that a measured critical kernel lies in a finite
Artinian model class.

**Testable conjecture, not a finding about P250:** part of its rank growth may
be translation-generated descendant/kinematic content. A properly justified
conformal or dilation kernel may predict held-out mesoscopic separations with
few amplitudes even when constant-shift rank keeps increasing.

This requires a field/instrument dictionary, exact physical offsets, charge and
hand conventions, and an actual scaling window. No target-dependent interpolation
onto a geometric grid is allowed. On a finite torus the translation representation
is finite but can grow with N. A torus kernel need not be one plane power or a
finite collection of BPZ solutions. The latter require the relevant degenerate
field and sector assumptions. Neither this example nor endpoint data determines
an ordered cover commutator.

## 2. Even one birth clock has full temporal covariance rank

Let `U` be uniform on (0,1) and `X(u)=1[U<=u]`. Then

\[
\operatorname{Cov}(X(u),X(v))=\min(u,v)-uv.
\]

At any `0<u_1<...<u_m<1`, the covariance matrix is positive definite. Indeed,
a linear combination of the indicators has zero variance only if its jumps
at every distinct `u_i` vanish. Its determinant is

\[
u_1(u_2-u_1)\cdots(u_m-u_{m-1})(1-u_m)>0.
\]

The checked uniform grid gives determinant `(m+1)^(-(m+1))` and rank m through
m=10. Thus a one-clock process can have arbitrarily large temporal covariance
rank. This process has an elementary time-inhomogeneous two-state Markov
realization; its covariance is not a finite-dimensional state count.

This does not reproduce P334's measured eigenvalue fractions, explain them, or
invalidate its rejection of a *separable* covariance ansatz. It specifies a
necessary null benchmark before interpreting those fractions as a count of
physical or predictive states.

## 3. An exact obstruction beyond the H2 one-step ceiling

Use the N10 square-NN torus with period columns `(3,1),(-1,3)`. Label site j by
`(0,j)` modulo the period lattice. Physical +x and +y step the labels by +3 and
+1 modulo ten. The two current occupied sets are

\[
A=\{0,1,2,3,4\},\qquad B=\{0,1,2,3,5\}.
\]

Both have five occupied sites, ambient rank one, primitive period-basis line
`(0,1)`, and one vacant single-site rank-two trigger. Their one-step exit
probability is exactly `1/5`. Their two-step survival probabilities differ:

| Current set | Single trigger | Minimal safe-site trigger pairs | Two-step survival |
|---|---|---|---:|
| A | 7 | {5,8}, {6,9} | 2/5 |
| B | 8 | {4,7}, {6,7}, {6,9} | 3/10 |

An independent enumeration of all 120 future permutations per configuration
gives exit counts `(24,48,48)` for A and `(24,60,36)` for B at steps 1,2,3.
Consequently `(k,rank,ell,H2)` is insufficient for the full future birth law.
No claim of smallest possible counterexample over all geometries is made.

For a general current rank-one configuration A with q vacant sites, let b1
count singleton triggers and b2 count unordered pairs that are individually
safe but jointly reach rank two. Monotonicity gives exactly

\[
s_1(A)=1-b_1/q,\qquad
s_2(A)=\frac{\binom{q-b_1}{2}-b_2}{\binom q2}\quad(q\ge2).
\]

Here binomial(n,2)=0 for n<2. The extra coordinate is a cooperative present-state
topological trigger count, not an arbitrary learned temporal mode.

For any valid future horizon,

\[
s_m(A)=\frac{\#\{U\subset A^c:|U|=m,\ r(A\cup U)=1\}}{\binom qm}.
\]

Equivalently, kill the uniform insertion kernel when rank reaches two:

\[
(Q_kf)(A)=\frac1{N-k}\sum_{v\notin A,\ r(A+v)=1}f(A+v),\qquad
s_m=Q_k\cdots Q_{k+m-1}1.
\]

The entire survival vector determines the remaining unmarked rank path, which
has only one exit. It does not determine spatial marks or the full configuration.
The minimal triggering subsets form a reliability hypergraph; higher-layer
counts without their overlaps generally do not determine its reliability law.

Conditional on the current set, the past and future orders of a uniform
permutation are independent. Therefore residual age association under coarse
conditioning does not prove intrinsic full-state history dependence. The
full-state near-critical Markov construction in [R3,R4] provides a controlled
continuum comparison, not a theorem that a few square-site scalars close.

**Next discriminator:** finish the existing k0 pilot unchanged. In a separately
specified continuation, test whether b1+b2 or a small survival vector predicts
longer held-out horizons beyond the cheap geometry variables. Exact one- and
two-step identities are calibration ceilings, not discovery scores. Current-state
snapshots are required; old clock-only archives do not recover b2. Conditional
future replicas share the outer configuration and must retain cluster covariance.

## 4. Separate a typed correction from a free fitted exponent

For a typed direct `0->2` birth event `E_c(v,A)`, independent uniform priorities
imply the exact path-probability identity

\[
D_c(N)=\sum_v\int_0^1\Pr_p^{V\setminus\{v\}}(E_c(v,A))\,dp.
\]

There is no extra p or (1-p) factor: the central site's priority has uniform
density. A directed predecessor edge with k occupied sites contributes
`k!(N-k-1)!/N!`. On N10 the untyped direct-edge census has 40 edges at k=5
and 40 at k=6. Its weighted sum is `5/63`; a separate Boolean-lattice dynamic
program counts 288,000 direct-birth permutations out of 3,628,800 and agrees.
This is not a new typed carrier implementation; use the existing carrier oracle
for that decomposition.

For the latest theta/figure-eight split, impose the assumptions separately:
a window of width `L^-3/4`, a nonzero scale-stable global type fraction, an
integrable near-critical profile with negligible off-window tail, and transport
of the controlled triangular arm exponents to square-site percolation. Then
[R5,R6] motivate the conditional arithmetic

\[
D_j=L^{2-3/4-\alpha_j+o(1)},\qquad \alpha_j=(j^2-1)/12.
\]

With a genuine asymptotic profile and amplitude, this becomes

\[
D_{\theta}\sim A_\theta N^{-5/6},\qquad
D_{8}\sim A_8N^{-2},\qquad D_8/D_\theta\sim(A_8/A_\theta)N^{-7/6}.
\]

The topology-specific A8 is directly measurable and nonnegative. Subtract the
observed figure-eight component rather than choose its amplitude from a fit to
D_total. If it is too small to explain the total drift, reject that correction
story and study theta's own corrections/type frequency. A regular-variation
assumption is additionally needed to turn exponent arithmetic into fixed
doubling-ratio limits.

The exact gluing work has not yet proved a scale-uniform positive frequency of
the requisite global landing type. Ordinary six arms do not suffice. Also,
`alpha_8=21/4` is not an identification with spin-four thermal Q4: the unweighted
arm exponent and the spin/representation of that descendant are distinct data.

## Reproduction and verification boundary

```bash
python scripts/predictive_state_counterexamples.py \
  --output /tmp/predictive-state-controls.json
cmp /tmp/predictive-state-controls.json \
  results/exact-predictive-state-controls/20260830.json
python -m unittest discover -s tests \
  -p 'test_predictive_state_counterexamples.py' -v
python -m compileall -q scripts/predictive_state_counterexamples.py \
  tests/test_predictive_state_counterexamples.py
```

The focused suite has **19 passing tests**. Lifted BFS and independently coded
weighted union-find agree on all **1,024** N10 configurations. All **310**
rank-one states satisfy the one/two-step controls; **1,650** valid state/horizon
pairs agree between killed-kernel recursion and direct subset enumeration.
The two displayed witnesses additionally use explicit future permutations.

Runtime used: Python 3.13.5, standard library only. Python 3.9 grammar parsing
also passed; no claim of running a Python 3.9 interpreter is made. These tests
were run in an isolated contribution directory, not as the whole-repository
suite. The manifest records source/output hashes. No RNG or Monte Carlo samples
were used. No production protocol, historical result, navigation file, or
current claim ledger is modified.

## arXiv reading map and relevance

The following 13 primary records were checked during the research sweep.
Reading depth is explicit: selected relevant sections for the technical anchors,
abstract/metadata screening for the broader horizon; this is not a claim to have
verified every proof or exhaustively searched arXiv.

| ID | Work and checked version | Reading/use and boundary |
|---|---|---|
| R1 | B. Mourrain, [1609.05720v3](https://arxiv.org/abs/1609.05720v3), *Polynomial-exponential decomposition from moments*; revised 20 October 2017 | HTML Theorem 3.1 and sampled-value discussion in Section 5.3: finite constant-shift Hankel rank restricts the function class. An HTML rendering date is not a new paper revision. |
| R2 | A. Bernardi, J. Jelisiejew, O. Reig Fité, [2606.30600v1](https://arxiv.org/abs/2606.30600v1), *Hankel and Multiplication Tensor Completions for Cactus Rank*; 29 June 2026 | HTML completion formulation: useful algebraic tooling after the finite Artinian model is justified, not a physical justification. |
| R3 | C. Garban, G. Pete, O. Schramm, [1305.5526v4](https://arxiv.org/abs/1305.5526v4), *The scaling limits of near-critical and dynamical percolation* | HTML Section 8 and Theorems 11.3/11.4: full-state pivotal evolution and Markov property on quad-crossing space; triangular-model scope. |
| R4 | C. Garban, G. Pete, O. Schramm, [1008.1378v5](https://arxiv.org/abs/1008.1378v5), *Pivotal, cluster and interface measures for critical planar percolation* | Abstract/metadata: pivotal measure as the geometric continuum input; not scalar-state sufficiency. |
| R5 | S. Smirnov, W. Werner, [math/0109120v2](https://arxiv.org/abs/math/0109120v2), *Critical exponents for two-dimensional percolation* | HTML Theorem 4: plane polychromatic arm exponent. Does not prove square-site universality or constant global type fractions. |
| R6 | P. Nolin, [0711.4948](https://arxiv.org/abs/0711.4948), *Near-critical percolation in two dimensions* | HTML near-critical arm framework: scope/colour/scale conditions matter for the priority integral. |
| R7 | F. Camia, Y. Feng, [2508.16047v2](https://arxiv.org/abs/2508.16047v2), *The percolation energy field and its logarithmic partner*; revised 1 June 2026 | Abstract/revision check: externally defined triangular logarithmic control; not an automatic dictionary for a square-site birth source. |
| R8 | P. Roux, S. Ribault, J. L. Jacobsen, [2604.24491v1](https://arxiv.org/abs/2604.24491v1), *Torus one-point functions in critical loop models*; 27 April 2026 | Abstract: torus/sphere bootstrap correspondence and infinite conformal-block combinations. Few named primaries do not mean finite functional rank. |
| R9 | M. Ang et al., [2604.05503](https://arxiv.org/abs/2604.05503), *Exact solution of three-point functions in critical loop models*; 7 April 2026 | Abstract: proposed structure constants with bootstrap, lattice and continuum checks; use after charged-insertion normalization is specified. |
| R10 | F. Camia, V. F. Foit, R. Nivesvivat, [2605.04395](https://arxiv.org/abs/2605.04395), *Anchored random clusters and SLE excursions*; 6 May 2026 | Abstract: anchored/pivotal spatial shapes are relevant to existing birth-site work, but a torus-cut dictionary remains necessary. |
| R11 | C. Alves et al., [2606.11503](https://arxiv.org/abs/2606.11503), *Percolation on hierarchical lattices*; 9 June 2026 | Abstract: useful controlled recursive-geometry comparison; no direct square-site map is supplied. |
| R12 | S. Diskin et al., [2603.03257](https://arxiv.org/abs/2603.03257), *Supercritical sharpness of percolation*; 3 March 2026 | Abstract: broad finite-cluster tail control on transitive graphs; not a ready quantitative rate for the Matching-One root. |
| R13 | E. Borel et al., [2608.17073](https://arxiv.org/abs/2608.17073), *Near-critical percolation with sparse reinforcements*; 17 August 2026 | Abstract: a recent dependent-enhancement setting. Its ensemble differs from the iid target; it cannot simply be imported as a new iid threshold bound. |

The useful novelty is the collision of these mathematical controls with the
current live results, not relabelling already-known literature as a new field
identification. The proposed order is: **type the generator; calibrate rank;
measure present cooperative geometry; resolve carrier contributions; then use
modulus, charged and OPE data to name any surviving continuum sector.**
