# Square 4/8 birth transfer and the population-IQR clock

2026-09-29. Continuation of the
[two-time density derivation](two-time-pivotal-density-bound-20260929.md),
for #778 / #773. Author proof using named standard probability inputs;
not independent certification, an exact exponent calculation, or a novelty claim.

**Outcome.** The preceding argument transfers to square-site percolation on
its own pivotal clock. It does not have to wait for a proof of the exponent
5/4. The relevant NN/matching arm conventions are compatible, and the
required qualitative power bounds and near-critical comparisons belong to
the general-lattice theory. Every subsequential joint birth limit is locally
absolutely continuous and has two distinct finite births on a fixed
nondegenerate torus.

A second argument connects this clock to the **population** IQR of the
archived insertion births: the width is bounded above and below by constant
multiples of `N*a_eta`. Hence the IQR normalization cannot by itself manufacture
or remove a limiting zero-gap atom. No limiting width constant, copula
uniqueness, finite-sample guarantee or cross-lattice equality is asserted.

## 1. Exact graph and arm conventions

Black occupied sites use the square nearest-neighbour graph `G`. White sites
use the matching graph `G*`, obtained by adding both diagonals in every
square face. Only the original square vertices are random Bernoulli sites.
An auxiliary face centre used to draw white paths is **not** an independent
site, a second source or an extra occupation variable.

Let `A_(2j)(r,R)` mean alternating black NN and white matching arms in a
square annulus. Arms are truncated to proper inner-to-outer crossings, with
no other boundary visits; fixed lattice collars are harmless. Their cyclic
order is alternating. Set

```text
pi_4(r,R) = P_pc[A_4(r,R)],
pi_6(r,R) = P_pc[A_6(r,R)].
```

The critical parameter is square **site** `pc`, not the square-bond value
1/2. There is no need to insert a numerical estimate of pc in this proof.

### Compatibility with the planarised-white separator

The stronger-looking face-capacity convention in #823 is compatible with
these alternating events. Fix the selected black proper crossing arcs and
cut the annulus along them. Each white arm starts in a different intervening
strip, by alternation. A white matching edge cannot cross a black NN arc:
its diagonal lies inside a square face and its endpoints are white.
Therefore each selected white path stays in its own strip.

An open square-face interior belongs to just one such strip, since black
NN arcs run on its boundary. Two white arms in different strips cannot use
the same face interior. Inside each strip, erase loops and replace a white
diagonal by the two edges through a face centre. The resulting paths cannot
share face centres or white vertices with paths in another strip. Conversely,
a face-centre path projects to a white matching path. Thus planarisation
does not discard the proper alternating crossings used here. Fixed collar
changes handle contacts with the annulus boundaries; this is not an equality
of the earlier degenerate L3/L4 thin-annulus census conventions.

This argument uses the alternating **whole collection**. It does not say that
arbitrary vertex-disjoint white paths, without intervening black arms, are
geometrically disjoint: crossing white diagonals alone would be a counterexample.

The four-arm event rooted at a pivotal site's neighbours and the annular
event with fixed microscopic radius have comparable probabilities. Arm
separation and a fixed finite-energy extension join their endpoints in a
bounded collar. The collar cost is uniform because p stays in a compact
subinterval of `(0,1)`. No extra mesh-dependent power is introduced.

## 2. Probability inputs actually used

The primary sources read for this transfer are:

1. Harry Kesten, *Scaling relations for 2D-percolation* (1987),
   [DOI](https://doi.org/10.1007/BF01205674),
   [full text](https://math.bme.hu/~balint/oktatas/perkolacio/percolation_papers/kesten_2.pdf).
   Page 112 explicitly includes square site and matching graphs. Lemmas 4--6
   supply separation/gluing, Lemma 7 (2.59) supplies the two-scale gain,
   Lemma 8 compares pivotal probabilities below
   characteristic length, and (4.5), (4.28), (4.33), Theorem 4 give the
   characteristic-length relation and the two-sided comparison.
2. Pierre Nolin, *Near-critical percolation in two dimensions*,
   [arXiv:0711.4948](https://arxiv.org/abs/0711.4948).
   The downloaded text's Proposition 16, Theorem 26, Proposition 32 and
   Lemma 37 give the arm/length formulations; section 8.1 explicitly discusses
   their square/matching extension and the strict four-/six-arm power bounds.
   This transfer uses that general-lattice scope, not the triangular exact
   values. Journal numbering differs.
3. Van den Berg and Nolin, *On the four-arm exponent for 2D percolation at
   criticality*, [arXiv:2008.01606](https://arxiv.org/abs/2008.01606),
   [author full text](https://ir.cwi.nl/pub/31186/31186.pdf).
   Section 2 fixes the square NN/matching arm convention. The end of section 5
   proves `pi_4(3n)<=C*n^-1*sqrt(pi_2(n))`; RSW decay of `pi_2` yields the
   strict power above one used for the label/insertion comparison below.

Here is the precise package consumed, in lattice units. Constants need not
be known numerically. There are positive `kappa_4,kappa_6,beta` such that

```text
c (r/R)^(2-kappa_4) <= pi_4(r,R) <= C (r/R)^beta,
pi_6(r,R) <= C (r/R)^(2+kappa_6),                        (1)
pi_4(r,R) asymp pi_4(r,s)*pi_4(s,R),                    (2)
```

uniformly over nondegenerate annuli, with fixed-radius changes absorbed in
constants. The upper power beta may be decreased so that `beta<2`.
The statements are uniform inequalities, not assumptions that limiting
critical exponents exist.

For clarity, the lower bound in (1) is an **annular**, not just a rooted,
input. It can be recovered from Kesten's (2.59) (printed page 137, checked
against the rendered formula). For dyadic r=2^j<R=2^k, divide that estimate
by `|p-pc|` and let p tend to pc. The finite-box probabilities are
polynomials in p, and the critical square-crossing derivative obeys the
standard Russo/separation comparison `sigma'_pc(R) asymp R^2*pi_4(1,R)`.
Thus, for a fixed zeta>0,

```text
r^2*pi_4(1,r) <= C*(r/R)^zeta*R^2*pi_4(1,R).
```

Quasi-multiplicativity gives
`pi_4(r,R)>=c*(r/R)^(2-zeta)`. Dyadic rounding and fixed collars yield
the stated uniform form (decrease zeta if needed). This explicitly explains
why a mere claim about a limiting rooted exponent is not being substituted
for the two-scale estimate consumed by the pair sum.

In addition, below characteristic length the
corresponding near-critical arm probabilities are comparable to critical
ones, and

```text
|p-pc| L(p)^2 pi_4(1,L(p)) asymp 1.                    (3)
```

Subcritical rectangle crossings decay exponentially in size divided by
`L(p)`. The matching graph has its own critical density `pc*=1-pc`; (3) and
the crossing comparisons apply to both members. These are literature inputs,
not theorems rederived by the finite square census.

**What is not imported:** triangular colour-switching between arbitrary
polychromatic words, 5/4, 35/12, 4/3, SLE6 convergence, or square conformal
invariance. We use only the indicated alternating words and matching-pair
estimates. Colour exchange of the entire graph pair is different from a
claim that arbitrary colour words on the same non-self-matching graph agree.

## 3. Square joint-birth corollary

Let the mesh eta square tori approximate a fixed flat torus of positive
systole and finite area, so `N asymp eta^-2`. In rescaled arm notation set

```text
a_eta = eta^2/pi_4(eta,1),
lambda_i = (T_i-pc)/a_eta,                i=1,2.
```

Every rank-changing site forces four arms by the punctured-disk replacement
argument in the preceding note. A direct `0->2` forces six (or eight) arms by
the offset/spine argument in #773 and the square separator in
[#823](https://github.com/LightChainr/Matching-One/blob/3d0176a9f394f7d7228766709731493fd15095c5/analysis/queue-20260919/815/NOTE.md).
Section 1 identifies the arm events used for their probability bounds.

For each bounded time window choose the local annulus radius below the
minimum physical characteristic length in that window. Relation (3) and
the lower four-arm power guarantee that this radius is a positive constant
independent of eta. This avoids applying near-critical comparison above its
stated range. Equations (1)--(2) make the two-time pair sum integrable:

```text
sum_{dyadic eta<=rho<=R} rho^2/pi_4(rho,R) <= C_R.
```

For the same-site atom, the two qualitative powers give the concrete
vanishing estimate

```text
pi_6(eta,R)/pi_4(eta,1) <= C_R eta^(kappa_4+kappa_6).    (4)
```

This is not the triangular 5/3 rate and does not determine the optimal square
rate. The earlier density proof now gives, for every finite Lambda,

```text
f_eta(s,t) <= C_Lambda,          -Lambda<=s<t<=Lambda,
P(gap<=delta, both births in the window)
   <= C_Lambda delta + C_Lambda eta^(kappa_4+kappa_6).    (5)
```

For tightness, apply (3) and the arm powers as in the earlier note to make
the physical characteristic length tend to zero when `|lambda|` tends to
infinity. Cover the torus by finitely many embedded charts. An essential
black cycle at the negative endpoint forces a fixed-size rectangle crossing,
whose probability tends to zero by subcritical exponential decay.

At the positive endpoint use the **white matching** configuration, not white
NN: the finite digital Alexander identity gives

```text
r_black,NN(p) + r_white,matching(p) = 2.
```

White occupation density is `1-p=pc*-Lambda*a_eta`. Swapping colours and
the graph pair leaves the unrestricted alternating four-arm event the same,
up to the already declared microscopic conventions. Its pivotal clock is
therefore comparable to `a_eta`. Matching subcritical decay controls
`P(lambda_2>Lambda)`. No false square-site symmetry `p <-> 1-p` on the
**same** NN graph is invoked.

**Corollary (author proof from the above standard inputs).** The square
birth pairs are tight on their pivotal clock. Every subsequential limit has
a locally bounded two-dimensional density and distinct finite coordinates;
the corresponding rank paths have two separate unit jumps in Skorokhod J1.
This completes the square input transfer claimed as conditional in the
earlier note. A unique law, exact clock power or cross-lattice law equality
does not follow.

## 4. The archive clock is comparable: a separate probability argument

The stored simulation records permutation birth counts `J1,J2`. Label time
uses independent sorted uniform labels; the uniform ordering permutation is
independent of those order statistics. In this coupling,

```text
T_i = U_(J_i),
E[T_i | J_i=j] = j/(N+1),
Var(T_i | J_i=j) = j*(N+1-j)/[(N+1)^2*(N+2)]
                <= 1/[4*(N+2)].                       (6)
```

The bound holds even though `J_i` depends on the entire permutation.
Consequently, with `X_i=(J_i/(N+1)-pc)/a_eta`,

```text
E[(X_i-lambda_i)^2] <= 1/[4*(N+2)*a_eta^2].             (7)
```

The square bound from van den Berg--Nolin gives
`pi_4(eta,1)<=C eta^(1+epsilon)` for some epsilon>0. Thus
`sqrt(N)*a_eta asymp eta/pi_4(eta,1) -> infinity`, and (7) tends to zero.
For triangular site the standard four-arm bound also gives this property.
Hence both individual birth times, not only their difference, have the same
subsequential laws in insertion and label clocks. This step needs a power
above one; the earlier pair-density integrability needed a power below two.

Let H_eta be the equally weighted **population** mixture of `X1,X2`, and let
`w_eta=Q_.75(H_eta)-Q_.25(H_eta)`, with lower quantiles. We claim

```text
0 < c <= w_eta <= C < infinity                         (8)
```

for all sufficiently small eta.

**Upper bound.** Tightness of the pair and (7) put more than 3/4 of the
mixture inside a fixed bounded interval; choose the window with each tail
less than 1/4. Both quartiles are then inside that interval.

**Lower bound.** Fix a window whose pair escape probability is less than a
small tau. Equation (5)'s joint density implies, uniformly over intervals I
of length h, for the label mixture,

```text
P(mixture label coordinate in I)
   <= tau + 2Lambda*C_Lambda*h + C_Lambda*e_eta.         (9)
```

Transfer to the insertion mixture by enlarging I by b at each end and adding
`P(|X_i-lambda_i|>b)`, which vanishes for each fixed b by (7). Choose tau,
b and then a fixed positive h so that the resulting upper bound is strictly
less than 1/2. Any interval between the .25 and .75 lower quantiles contains
at least half the distribution, including in the discrete case. Its length
cannot be less than h. This proves (8); no unproved convergence or strict
positivity of a limiting density at a quartile was assumed.

In the simulation's count units, write W_eta for the population IQR of the
equal mixture `(J1,J2)`. Positive affine scaling of lower quantiles gives

```text
W_eta = (N+1)*a_eta*w_eta asymp N*a_eta.                 (10)
```

Therefore the population-normalized gaps

```text
G_J = (J2-J1)/W_eta,
G_T = (N+1)*(T2-T1)/W_eta
```

differ by `o_P(1)` using (7)--(8). Their denominators relative to the pivotal
clock stay in a compact positive interval. In particular, since soft mass is
monotone in its resolution,

```text
lim_{delta->0} limsup_{eta->0} E[(1-G_T/delta)_+] = 0,   (11)
```

and the same holds for G_J. To justify the latter directly, split on
`|G_T-G_J|<=delta` and use the small-ball criterion at resolution `2delta`;
then take eta->0 before delta->0. Do not interchange these limits.

This supplies a population-level justification for the archive's exponent-free
normalization. It does **not** prove that `w_eta` converges, calculate it,
replace the empirical W by a known deterministic value, or make square and
triangular primitive tori have the same modulus.

## 5. What can and cannot now be said about the existing data

- The code's `label_IQR` uses exactly `(N+1)*(T2-T1)/W`, with W estimated from
  the pooled insertion births. Its target now has a clock-comparison argument.
- The population proof does not validate a fixed finite sample's empirical
  quantiles or imply a confidence guarantee for 14 batches. Those remain
  sampling questions; the existing delete-one calculation re-estimates W.
- No result file or production sample is changed or regenerated for this note.
  No exponent is fitted and no larger-size job is launched.
- The next missing process object is the actual joint limit / copula and its
  torus dependence, not another scalar direct-jump exponent. The same-site
  arm suppression, two-time density, tightness and width comparison address
  different links and should stay separately visible.

The literature input is old mathematics, not a claim that this project has
newly proved square critical exponents. The contribution here is the explicit
transfer to the ambient-rank birth process and its stored insertion readout.
