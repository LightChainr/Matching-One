# Two-time pivotal density: an arm-based route to distinct limiting births

2026-09-29. Consumer: #778 / #773. Author derivation, not an independently
reviewed theorem or a claim of literature novelty. No new simulation is used.

**Main result.** On a fixed nondegenerate torus, two different rank-changing
sites have a locally uniformly bounded **joint near-critical time density**,
provided macro-pivotality forces four alternating arms and the usual
near-critical arm estimates hold. The proof uses independent labels in
disjoint spatial annuli, not independence of the two births. Together with
six-arm suppression of a direct double birth and tightness of the birth pair,
this excludes coalescence of two distinct insertions in a scaling limit.

The triangular-site corollary below supplies the probability inputs from
standard near-critical theory. The square-site statement remains explicitly
input-conditional here: the precise 4/8 arm contract must be carried through,
but **an exact square four-arm exponent is not required**. The result does not
identify the limiting copula, prove square/triangular universality, or identify
the original-U operator.

This advances the earlier [soft-gap criterion](birth-gap-noncoalescence-20260929.md):
the two-time upper bound is now derived, rather than merely requested.

## 1. Model, clock and explicit hypotheses

Let a fixed flat torus have positive systole and finite area. Approximate it
by a periodic planar site lattice of mesh `eta`, with `O(eta^-2)` sites and
the usual `O(rho^2/eta^2)` local packing bound. Work in embedded disks smaller
than the injectivity radius. The torus shape is fixed, not an arbitrarily
thinning family.

Give each site an independent uniform label `U_x`. The occupied graph at p
contains exactly sites with `U_x<=p` and all their lattice edges. Let

```text
r_eta(p) = rank im[H1(occupied graph) -> H1(torus)],
T_i      = inf{p:r_eta(p)>=i},                 i=1,2,
a_eta    = eta^2 / alpha_4(eta,1),
lambda_i = (T_i-p_c)/a_eta.
```

Here `alpha_j` is the critical probability of a specified alternating arm
event, with a fixed microscopic inner convention. Changing that convention
or the fixed reference radius 1 changes comparison constants. Filling
contractible occupied faces does not change this ambient image rank.

For each finite `Lambda`, assume:

**A1 — deterministic local importance.** A site that changes ambient rank
has four alternating arms in its punctured neighbourhood to a fixed radius
`R=R_Lambda>0`. A direct rank `0->2` update has at least six alternating arms
to that radius. The radius can be decreased; it does not depend on eta.

**A2 — compact-window arm stability.** For
`p=p_c+u*a_eta`, `|u|<=Lambda`, four- and six-arm probabilities in embedded
annuli with outer radius at most R are bounded above by a constant times
their critical probabilities, uniformly over the inner radius. Conditioning
the centre site's label does not alter an arm event outside the centre.

**A3 — critical quasi-multiplicativity and integrability.** Up to fixed
changes of radii, `alpha_4(r,t)` is comparable to
`alpha_4(r,s)*alpha_4(s,t)`. For some `kappa>0`,

```text
alpha_4(r,R) >= c_R (r/R)^(2-kappa),     eta<=r<=R.       (1)
```

It is enough to assume directly that the dyadic sum
`sum_{rho=R*2^-j, rho>=eta} rho^2/alpha_4(rho,R)` stays bounded.
No equality with a power law is needed.

**A4 — direct-jump suppression.**

```text
e_eta(R) := alpha_6(eta,R)/alpha_4(eta,1) -> 0.           (2)
```

**A5 — birth-pair tightness.** The pairs `(lambda_1,lambda_2)` are tight on
the finite real plane. Section 5 supplies a sufficient finite-size argument;
this is not inferred from a finite Monte Carlo window.

## 2. The two-time bound

### Proposition: uniformly bounded off-diagonal density

Under A1--A3, the law of `(lambda_1,lambda_2)` restricted to different birth
sites has a density `f_eta(s,t)` on `s<t`. For every finite Lambda,

```text
0 <= f_eta(s,t) <= C_Lambda,
              -Lambda<=s<t<=Lambda,                    (3)
```

uniformly for sufficiently small eta. A1--A4 also give

```text
P(lambda_1=lambda_2 in [-Lambda,Lambda])
    <= C_Lambda e_eta(R).                               (4)
```

Consequently, for `0<delta<=1`,

```text
P(0<=lambda_2-lambda_1<=delta,
  lambda_1,lambda_2 in [-Lambda,Lambda])
    <= C_Lambda delta + C_Lambda e_eta(R).               (5)
```

The constants depend on the window and geometry. This is not a global
Lambda-independent density or tail bound.

### Proof: fix the two labels before discarding the birth conditions

For distinct sites x,y, condition on
`U_x=p_c+s*a_eta`, `U_y=p_c+t*a_eta`, with `s<t`. Their joint label density
in `(s,t)` is exactly `a_eta^2`. All other labels remain independent uniforms.
The ordered birth-site events partition the different-site event, so

```text
f_eta(s,t) = a_eta^2 * sum_{x!=y}
  P(x causes birth 1 and y causes birth 2 | U_x=p(s), U_y=p(t)).   (6)
```

Conditional probabilities at prescribed label values can be defined by
integrating over the remaining finite set of labels; (6) holds almost
everywhere and gives a density representative with the asserted bound.

Let `rho=dist(x,y)`, initially `eta << rho << R`. The first pivotality
requires four arms around x at time s; the second requires four arms around
y at time t. Retain only these necessary subevents:

1. four arms around x from microscopic scale to `rho/8`, at s;
2. four arms around y from microscopic scale to `rho/8`, at t;
3. four arms around x from `2rho` to R, at s.

Their supports are spatially disjoint and avoid **both conditioned sites**.
Each is a function of the labels in its own support. They are therefore
independent even though the first and third use one threshold and the second
uses another. Arm stability and fixed-radius comparisons bound (6)'s
conditional probability by

```text
C_Lambda alpha_4(eta,rho)^2 alpha_4(rho,R)
    <= C_Lambda alpha_4(eta,R)^2 / alpha_4(rho,R).        (7)
```

This is not a factorisation of two overlapping global pivotal events, and
not a Markov approximation of the rank history. No new mixed-time arm theorem
is smuggled into A2: only single-time arm events on disjoint sets are used.

For `rho` within a fixed factor of eta, omit the two inner events and use the
outer event starting at a sufficiently large fixed multiple of eta. Its
bound agrees with (7) up to constants. For `rho>=R/8`, retain disjoint
four-arm disks of a smaller fixed radius around x and y. Their product gives
the same required uniform far-pair contribution.

For a close dyadic distance bin, the number of ordered pairs is at most
`C eta^-4 rho^2`. Substituting (7) in (6) and using
`alpha_4(eta,R) <= C_R alpha_4(eta,1)` gives

```text
f_eta(s,t)
  <= C_Lambda [1 + sum_{dyadic eta<=rho<=R}
                          rho^2/alpha_4(rho,R)]
  <= C_Lambda [1 + C_R sum_j (R*2^-j)^kappa]
  <= C'_Lambda.                                        (8)
```

The cancellation of `eta^-4` against `a_eta^2 alpha_4(eta,1)^2=eta^4`
is the useful clock normalization. The spatial singularity is integrable
because the four-arm power is strictly below two.

A same-site double birth has only one label density. Sum its six-arm bound
over sites and integrate over the time window:

```text
O(eta^-2) * O(Lambda*a_eta) * C_Lambda alpha_6(eta,R)
    <= C'_Lambda e_eta(R).                             (9)
```

Finally integrate (3) over the triangular band `0<t-s<=delta` inside the
window, whose area is at most `2Lambda*delta`, and add (9). This proves
(3)--(5).

## 3. Limiting consequences and the remaining unidentified object

Under A1--A5, every weak subsequential birth-pair limit has a locally bounded
two-dimensional density and is supported on `lambda_1<lambda_2`. To see the
density assertion, restrict to a compact square: the off-diagonal measures
are dominated by `C_Lambda` times area, while the diagonal measure tends to
zero by (4). Domination passes to the weak limit, using a slightly larger
square at its boundary. Tightness prevents mass escaping to infinity.

Equivalently, if `G_eta=lambda_2-lambda_1` and
`tau_eta(Lambda)=P((lambda_1,lambda_2) outside [-Lambda,Lambda]^2)`, then

```text
E[(1-G_eta/delta)_+]
   <= tau_eta(Lambda) + C_Lambda delta + C_Lambda e_eta(R).   (10)
```

Take `limsup_eta->0`, then `delta->0`, then `Lambda->infinity`. The result is
zero, supplying exactly the earlier noncoalescence criterion. One must not
exchange those limits or erase the tail term to claim a global O(delta) law.

There is also a path consequence: along a convergent birth-pair subsequence,
the two-jump rank paths converge in Skorokhod J1 on compact intervals whose
endpoints are not limiting birth times. Distinct ordered jumps can be matched
by a small time change. The local density makes every fixed endpoint a
non-birth almost surely. Thus the resulting subsequential paths have two
separate unit jumps, not a merged jump of size two.

**Not supplied:** uniqueness of that subsequential law; its copula or density
formula; a Markov property; a rank-one direction-mark limit; an explicit
constant for the measured soft mass. Distinct births do not determine their
dependence. Identification of the torus observable in the continuum process
is a separate, more informative next theorem.

## 4. Why the local arm hypotheses fit triangular ambient rank

Use regular triangular-site percolation, whose white matching connectivity
is also triangular NN connectivity. All local disks below are embedded in
the torus and have lattice polygon boundaries, up to harmless fixed changes
in radius.

**Four arms for any rank change.** Delete the centre x and examine black
connections inside a disk D around it. If all black branches incident to x
which reach the boundary of D are connected to each other within `D\{x}`,
every traversal through x in a global black cycle can be replaced by a
black traversal in that one local component. Branches not reaching the
boundary only add cycles supported in D. The replacement differs from the
original cycle by cycles in a disk, which have zero ambient homology.
Hence opening x cannot change ambient rank in this case.

A rank change therefore requires at least two different components in the
punctured disk, each incident to x and reaching its outer boundary. Select
two disjoint black representatives. Planar separation in the two intervening
strips supplies a white path from the microscopic hole to the outer boundary
in each strip. They are disjoint and alternate with the black paths: four
arms. This argument concerns local components, which can be connected outside
D; it does not wrongly require two different global clusters.

**Six arms for a direct double birth.** Before such a birth the ambient
black rank is zero. Each lifted black component and its nonzero deck
translates are disjoint. The rank increment from adjoining x is generated
by differences of deck offsets among its contacts in each base component.
Two independent differences require either three different lifted copies
of one component, or two pairs from different components. The attachment
paths joining the corresponding translated neighbourhoods reach a fixed
fraction of the systole; their initial branches give respectively three or
four different lifted black components to that scale. Alternating white
separators give six or eight arms, and eight implies six at the macroscopic
annular level. The same offset argument works for six triangular neighbours;
it does not use the square degree-four restriction.

These are deterministic disk/lift arguments, not an extrapolation from the
square L3--L5 census. The existing square-NN/white-matching version is recorded
in #773's [attachment/arm note](https://github.com/LightChainr/Matching-One/blob/9ef22d1cc531c09dedf9b1195ab6aad8ef79668f/notes/rank-jump-two-arm-extraction-20260914.md)
and #823's [separator proof](https://github.com/LightChainr/Matching-One/blob/3d0176a9f394f7d7228766709731493fd15095c5/analysis/queue-20260919/815/NOTE.md).
The triangular argument uses ordinary planar white paths; it does not identify
crossing square diagonals as disjoint white arms.

## 5. Probability inputs and triangular corollary

Two primary papers were read with the arXiv retrieval skill:

- Pierre Nolin, [Near-critical percolation in two dimensions](https://arxiv.org/abs/0711.4948),
  DOI [10.1214/EJP.v13-565](https://doi.org/10.1214/EJP.v13-565).
  In the downloaded arXiv text: Theorem 26 (near-critical arm comparison),
  Theorem 20 (polychromatic exponents), Proposition 32 / (7.13) (finite-size
  relation), Lemma 37 (uniform subcritical exponential decay), and section 8.1
  (general-lattice scope). Numbering differs in the journal version.
- Garban, Pete and Schramm,
  [The scaling limits of near-critical and dynamical percolation](https://arxiv.org/abs/1305.5526),
  [published JEMS version](https://ems.press/journals/jems/articles/15407),
  DOI 10.4171/JEMS/786. Definition 1.2 uses an exponential-clock convention;
  section 10.2 gives correlation-length control, and section 11.2 discusses
  near-critical arms. Its quad-crossing path convergence is **not** quoted as
  an already-proved torus-rank path continuity theorem.

### Local estimates

For a compact time window, choose R smaller than the torus injectivity radius
and the minimum physical characteristic length at its two endpoints. R can
depend on Lambda. Arm comparison then gives A2 without extrapolating a
below-correlation-length theorem to arbitrarily large radii. Standard critical
quasi-multiplicativity and the four-arm bound give A3. For the triangular
model, the four-/six-arm exponents yield

```text
e_eta(R) = eta^(35/12-5/4+o(1)) = eta^(5/3+o(1)).        (11)
```

Only its vanishing is consumed. This is a **direct-atom error bound**, not a
newly fitted gap exponent. GPS's weaker six-arm estimate also suffices.
For their exponential clock, passing to the affine p-clock used here has
derivative bounded above and below on each fixed window for small eta
(at p_c=1/2 its limiting time factor is 2). No clock equality is assumed.

### Tightness from finite-size decay, not from mean gap

Here is the torus reduction needed for A5. Write `n=eta^-1`, and let `m` be
the subcritical characteristic length in lattice units at
`p=p_c-Lambda*a_eta`. The finite-size relation and quasi-multiplicativity give,
when `m<=n`,

```text
1 asymp Lambda * (m/n)^2 / alpha_4(m,n).                 (12)
```

Usual critical annular bounds
`c s^(2-kappa) <= alpha_4(sn,n) <= C s^beta`, with
`0<beta<2`, imply

```text
c Lambda^(-1/kappa) <= m/n <= C Lambda^(-1/(2-beta)).    (13)
```

Constants can be enlarged for bounded Lambda. For large Lambda the same
finite-size relation forces `m<n`; otherwise `m^2 alpha_4(1,m)` is, by the
four-arm lower bound and quasi-multiplicativity, at least a fixed multiple
of `n^2 alpha_4(1,n)`, contradicting (12)'s defining relation at large Lambda.
Thus the physical characteristic length stays positive on each fixed window
and tends to zero as Lambda grows. Exact 4/3 scaling is unnecessary here.

Cover the fixed torus by finitely many embedded charts of a fixed radius
smaller than its systole. An essential black cycle must cross one of a finite
collection of annuli in these charts. Such a crossing forces a crossing of
one of finitely many rectangles with fixed positive dimensions. Uniform
subcritical exponential decay therefore makes its probability tend to zero
as `m/n->0`. There is no union over `eta^-2` individual starting sites.
It follows that

```text
lim_{Lambda->infinity} limsup_{eta->0} P(lambda_1<-Lambda)=0.
```

For triangular site percolation, complementary ambient ranks sum to two and
the white graph has the same lattice. Applying the same argument to white
at `1-p` gives the upper tail `P(lambda_2>Lambda)`. The other two tails are
smaller by ordering. This proves A5.

**Triangular corollary (author proof using the stated standard inputs).**
On any fixed nondegenerate triangular torus approximation, the near-critical
birth pairs are tight; every subsequential limit has a locally bounded
two-dimensional density and almost surely two distinct finite births.
The rank paths consequently have only unit jumps in those subsequential
limits. Unique continuum law is not asserted.

### Square-site consequence is conditional, but exact exponents are not the gap

**Follow-up:** [Square birth transfer and the population-IQR clock](square-birth-transfer-and-iqr-clock-20260929.md)
supplies the convention/standard-input transfer and the nondegenerate population
width comparison left open in this section and item 3 below. The original
conditional formulation is retained here to separate the two proof steps.

For square black NN / white matching connectivity, the same theorem consumes
the actual alternating 4/8 arms. Nolin's general-lattice discussion indicates
the relevant RSW-based ingredients, including four-arm power below two and
six-arm power above two. It is therefore too strong to say that this strategy
must wait for a proved square exponent of 5/4 or square conformal invariance.

This note does not silently identify that general discussion with every
endpoint and planarised-white convention in #823. A square upgrade should
state the arm comparisons and quasi-multiplicativity under that convention,
and use the square/matching complementary finite-size tails. Once those
inputs are supplied, (6)--(10) require no new two-time argument. It would
prove noncoalescence on the model's own pivotal clock, **not** equality with
the triangular clock, amplitudes or limiting copula.

## 6. What changes in the research allocation

1. The triangular positive-control target is now a proof-level density result,
   not just a request for more samples or a numerical extrapolation.
2. The next joint-process work is the square 4/8 input transfer and continuum
   torus-observable identification / uniqueness. This is narrower than
   reproving the entire near-critical construction.
3. The measured width-normalized gap uses the mixture IQR of `(J1,J2)`.
   Equating it with the canonical pivotal clock requires a nondegenerate
   width comparison; this note does not insert an unmeasured conversion
   constant. The existing pilot is compatible evidence, not proof of (3).
4. Do not infer a Markov process, independent births, a universal density,
   or an original-U mechanism from this anti-concentration result.

No Monte Carlo block, exact torus census, old test suite or cloud session was
run for this derivation. External literature was retrieved to check the
actual probability inputs. The arXiv search API did not return promptly;
known-paper retrieval succeeded, so no completeness-of-search claim is made.
