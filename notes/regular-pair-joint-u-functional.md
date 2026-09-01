# Canonical regular-pair J2 is one exact source response

**Result.** For the fixed canonical Kreg and a vacant-site perturbation
`T_0 -> T_0+epsilon Kreg/N`, the next original-U readout is exactly

\[
 J_2(N)=\left.\partial_{\log Q}\partial_\epsilon^2U
       \right|_{Q=1,\epsilon=0}
       =\mathcal L_{U,N}[S_2],\qquad
 S_2(A)=\frac{2}{N^2}\sum_{x<y\,;\,x,y\notin A}g_{xy}(A).       \tag{1}
\]

Here g_xy is the **joint physical** two-kernel Q derivative, and L_U,N
is the existing complete one-source moving-root/slope functional.
For every translation-invariant O, including 1,q,E and their K-weighted
thermal moments, the exact origin reduction is

\[
 \langle O S_2\rangle
   =\left\langle O\,\frac{\mathbf1_{0\notin A}}N
       \sum_{y\ne0\,;\,y\notin A}g_{0y}(A)\right\rangle.        \tag{2}
\]

There is no remaining factor two in (2). Its integer source numerator is
`sum_y g16(0,y)` with scale **16N**, not 8N, 16N² or a conditional
vacancy normalization. Adjacent empty vertices must share their actual
common edge-port identifier. This note derives these identities without
counting configurations or scoring J2.

## 1. The fixed kernel and both vacancy conditions

Use the completed canonical `Kreg=average_C4 i(I-P1)i^dagger`, whose
singlet completion coefficient is identically one. The
[spatial kernel](https://github.com/LightChainr/Matching-One/blob/a237968f1d7a82d26b46e83c58179dbba7f1a908/notes/regular-pair-spatial-kernel.md)
orders eight ports as

```text
(xN,xE,xS,xW,yN,yE,yS,yW).
```

For two distinct vacant sites, let pi_xy(A) be the exact component
partition of those edge-nodes in the unmodified occupation's exterior
hypergraph, and b its number of blocks. Define

\[
 \beta_{xy}(Q,A)=Q^{-b}\sum_{\mathrm{colours}}
       K_{{\rm reg},x}K_{{\rm reg},y}D_{\pi_{xy}(A)},\qquad
 g_{xy}(A)=\partial_Q\beta_{xy}(1,A).                         \tag{3}
\]

Other exterior components cancel from this relative coefficient. If
either site is occupied its insertion coefficient is zero. The sparse
lookup stores `g16(pi)=16 g_pi`, with canonical first-appearance labels
in the stated common port orientation. A missing **valid canonical**
key is exactly zero. The kernel is signed and is symmetric under
exchanging the two sites, including the exchange of their four-port
groups. q and E always use the original occupation/rank, never virtual
connections supplied by a kernel diagram.

The two-kernel definition applies to all physical eight-port partitions,
including adjacent sites. The existing noncontact spatial production is
a narrower sampling contract, not a restriction on the kernel itself.

## 2. Why one Q derivative removes all first-insertion product terms

At fixed occupation A, divide the perturbed colour weight by its
unperturbed weight. Each original vacant vertex is affine in epsilon,
so the finite product has expansion

\[
 F_A(Q,\epsilon)
 =1+\frac{\epsilon}{N}\sum_{x\notin A}\beta_x(Q,A)
    +\frac{\epsilon^2}{N^2}
       \sum_{x<y\,;\,x,y\notin A}\beta_{xy}(Q,A)
    +O(\epsilon^3).                                        \tag{4}
\]

There is no x=y term: a second derivative of that one site's affine
tensor is zero. Using an exponential vertex interaction instead would
be a different model with extra local second-order terms.

Every nonempty Kreg insertion has zero Q1 contraction. Therefore
`beta_x(1,A)=beta_xy(1,A)=0`, for every A and every thermal parameter.
Writing delta=log Q, the full unnormalized occupation weight has the
first-delta expansion

\[
 w_A(Q,p,\epsilon)=w_A^B(p)\left[
  1+\delta\left(\ell_0(A,p)+\epsilon S_1(A)
                   +\frac{\epsilon^2}{2}S_2(A)
                   +O(\epsilon^3)\right)+O(\delta^2)\right]. \tag{5}
\]

Here ell_0 is the background Q-path derivative, while
`S1=N^(-1) sum_x partial_Q beta_x(1,A)`. Differentiating (4) twice
produces precisely the factor 2 in (1). The background occupation,
rank-factor and thermal-activity Q derivatives multiply zero nonempty
insertion coefficients; they cannot add an epsilon² term to (5).
At Q1, partial_Q and partial_logQ agree at this derivative order.

There is also no S1² term at first order in delta. The product of two
one-insertion coefficients is O(delta²), and differentiation once in Q
kills it. Equivalently, the first-Q derivative of log F_A has exactly
the same epsilon²/2 coefficient S2 as (5).

After summing occupations and normalizing, for any unchanged O,

\[
 \left.\partial_{\log Q}\partial_\epsilon^2\langle O\rangle
       \right|_{1,0}
   =\langle O S_2\rangle_B-\langle O\rangle_B\langle S_2\rangle_B
   =\operatorname{Cov}_B(O,S_2).                            \tag{6}
\]

All thermal derivatives of (6) obey the same identity: this is a finite
regular sum, and the zero-insertion statements hold for the entire
thermal curve. Thus normalization still supplies the ordinary source
centering in (6); what disappears are the *additional* background-Q and
products-of-first-insertion terms. In particular neither Cov(a_x,a_y)
nor Cov(t_x,t_y) represents S2. Multiplying already Q-activated one-mark
scores puts their product at the wrong Q order.

## 3. Complete original-U functional and the surviving root terms

For each geometry separately let
`m_g=<q>`, `e_g=<E>`, and define at its original common pooled root p0

\[
 M=\frac{m_a+m_b}{2},\quad Y=\frac{e_a-e_b}{\Delta_4},\quad
 M(p_0)=0,\quad D=M_p(p_0)\ne0,\quad R=Y_p/D,\quad
 A_N=N^{13/8}/2.                                           \tag{7}
\]

For the source S2 use the separately normalized covariances

```text
j_q,g = <q S2>_g - m_g <S2>_g,
j_E,g = <E S2>_g - e_g <S2>_g,
jM = (j_q,a+j_q,b)/2,
jY = (j_E,a-j_E,b)/Delta4.
```

The exact mixed root tangent and U response are

\[
 \partial_{\log Q}\partial_\epsilon^2p_\star\big|_{1,0}
       =-\frac{jM}{D},                                    \tag{8}
\]

\[
 \boxed{\frac{J_2}{A_N}
   =\frac{jY_p}{D}-\frac{Y_{pp}jM}{D^2}
       -R\frac{jM_p}{D}+R\frac{M_{pp}jM}{D^2}.}             \tag{9}
\]

These are the old direct, root-motion, source-slope and root-slope
terms. All are evaluated at the original Q1 pooled root, with the
new S2 source covariances. Although the baseline pure-epsilon derivatives
vanish at Q1, the mixed root movement in (8) generally does not vanish.
Neither that term nor the two slope responses may be dropped.

One way to prove that no other mixed terms occur is to regard the
complete root/slope statistic as a differentiable functional of the
finite occupation law. At Q1 the baseline law is independent of epsilon
for every p. Its first Q derivative is linear in the bracket of (5).
Taking epsilon² then selects only S2 under that same first-order
functional. First-Q background terms and first mixed epsilon terms
cannot form a nonlinear correction at first order in Q-1.

The equivalent intrinsic expression is

\[
 \frac{J_2}{A_N}
   =\frac1{M_p}\partial_p
       \left[jY-\frac{Y_p}{M_p}jM\right]_{p=p_0}.            \tag{10}
\]

Using the common thermal coordinate h=p/(1-p) replaces every p derivative
in (7)--(10) by h derivatives without changing J2. For the original N25
pair use `Delta4=1152/625` and its saved root, slope and q/E jets. S2 has
no rank1-only support theorem, so the direct q/E numerators must remain.

## 4. Translation reduction and the factor two

Extend g_xy by zero when either endpoint is occupied. Symmetry gives

\[
 S_2(A)=\frac1{N^2}\sum_x\sum_{y\ne x}g_{xy}(A).             \tag{11}
\]

On each Gaussian quotient, translations act transitively on the N sites
and preserve the physical edge directions, original rank, K and the
unperturbed Bernoulli law. For translation-invariant O every inner
expectation in (11) is identical after translating x to the origin.
Consequently

\[
 \langle O S_2\rangle
  =\frac1N\sum_{y\ne0}\langle O g_{0y}\rangle,               \tag{12}
\]

which is (2). This is an equality of first-source joint moments, not
a configurationwise identity between S2 and the origin-reduced source.
It holds separately at each K, so it survives all thermal derivatives.
It does not make different displacements independent samples or justify
replacing a higher source moment by a power of the reduced source.

For a configuration with the origin fixed vacant, set

```text
G16(A) = sum over y != 0, y vacant of g16(pi_0y(A)).
T0(A)  = G16(A)/(16*N).
```

Set T0=0 for origin-occupied configurations. Then the required source
moments are exactly `<T0>`, `<q T0>`, `<E T0>` with the **full** original
probability denominator. No extra factor 2 is applied to G16: the two
orders of each unordered pair already changed (1) into (11).

An implementation may traverse only the 2^(N-1) configurations with the
origin vacant, since all other numerator terms are zero. At N25 this is
2^24 configurations, but the baseline remains the full 2^25 occupation
population. Explicitly, accumulate integers

\[
 D_O(K)=\sum_{A:\,0\notin A,\,K(A)=K}O(A)G16(A),
       \qquad O=1,q,E.                                    \tag{13}
\]

For h=p/(1-p), the actual unconditioned source moments are

\[
 \boxed{\langle O S_2\rangle_p
 =\frac{\sum_{K=0}^{N-1}D_O(K)h^K}
             {16N(1+h)^N}.}                              \tag{14}
\]

Equivalently the fixed-K source moment is
`D_O(K)/(16N binom(N,K))`, not the denominator binom(N-1,K).
If a producer first forms a conditional origin-vacant mean, it must
restore its probability `1-p`, or `(N-K)/N` in a fixed-K slice.
At p=1/2 the full configuration denominator is 2^N rather than 2^(N-1).
For general p, the complete binomial polynomial in (14) is the relevant
denominator; the raw traversal count is not a probability normalization.

## 5. Adjacent vacant sites share an isolated physical edge

When y is the east neighbour of vacant x, xE and yW are the **same
physical edge-node**, even though they occupy two positions in the
eight-port list. Neither endpoint has an occupied site tensor, so this
edge is an isolated hypergraph component. Its partition block contains
both port occurrences and counts as one component shared by the marks.
No extra Q factor, local bond occupation or rank update is introduced.

For any incident port of a vacant mark:

- If its other endpoint is occupied, use that endpoint's occupied NN
  component identifier. All incident edges tied through that component
  share its colour.
- If its other endpoint is vacant, use a canonical identifier of the
  **physical undirected edge**, disjoint from the occupied-component ID
  namespace. Repeated appearances of the same edge get the same ID.

Only after assembling these IDs should the eight labels be canonicalized
in the fixed N,E,S,W order and packed for the lookup. In particular
inventing a different singleton for every port occurrence incorrectly
splits xE from yW. The shared-component-zero rule remains valid, but the
shared isolated edge counts toward it; one further shared occupied
component can bring an adjacent pair to the nonzero s=2 class. Adjacency
must not be discarded merely because the noncontact spatial experiment
did not sample it.

The existing spatial sampler at the pinned source uses `n+i` for an
isolated port and explicitly assumes all eight physical edges are
distinct. That is correct for its frozen separations r8/r16. It must
not be copied unchanged into an all-displacement J2 producer. On the
honest N25 quotients an undirected endpoint-pair edge ID suffices; a
more general multigraph would require the actual distinct edge identity.

## 6. The exact decision this permits

A model with first-Q effective log weight linear in epsilon has S2=0,
and therefore predicts J2=0. A nonzero J2 rejects that global additive
closure for the specified canonical interaction and original U. A zero
J2 need not imply S2=0: centering, directional projection, the root and
the slope can annihilate a nonzero source. In particular a common
`c0+c1 K` component has zero response under (9).

The completed spatial observer is `<g_xy>` at prescribed noncontact
separations. It has neither the all-displacement sum nor q/E/K source
crossmoments in (13), and does not establish J2. Positive conditional
Gram examples or positive selected-separation means also fix no sign
of (9). The present reduction adds no descriptor, completion parameter
or substitute covariance; it specifies the one missing joint source
and its exact original-observer transmission.

Before collecting the new joint moments, the same fixed comparison also
separates the four nearest-neighbour displacements from all other nonzero
displacements. Equations (11)--(14) apply to each translation-invariant
displacement class, and the linear functional (9) gives
`J2=J2_NN+J2_nonNN`. The identical 16N units apply to both parts. A model
in which only nearest-neighbour contact interactions transmit to original
U predicts `J2_nonNN=0`; the already positive spatial C does not decide
that null. This is one predetermined physical split in the same traversal,
not a selected support radius, a new source fit or an independent data
block. Non-NN on N25 does not mean macroscopic separation.

## Source pins and delivery

- Execution `a237968f1d7a82d26b46e83c58179dbba7f1a908`, `branch_only`:
  `notes/regular-pair-spatial-kernel.md`,
  `notes/regular-pair-spatial-observer.md`,
  `notes/regular-pair-activation-original-u.md` and
  `scripts/p337_regular_pair_spatial_sampler.cpp`.
- The pinned kernel lookup is `analysis/regular_pair_spatial_kernel.tsv`,
  SHA256 `36ae069d370b1d7a4398861c928afb41aa76885c8895c696b1bc0c97e9c314fd`;
  its unit is g16=16g, and valid omitted entries are zero.
- The canonical finite-network zero and old U definition are reused,
  not rescored. This delivery contains derivation only: no new occupation
  traversal, simulation, root search, numerical score, test or server job.
