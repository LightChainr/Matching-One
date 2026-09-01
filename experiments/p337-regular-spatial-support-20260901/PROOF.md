# Canonical regular pair kernel: a two-component spatial support lemma

**Conclusion.** At a fixed exterior occupation, if the two marked port groups share at most one exterior equality component, their normalized double-insertion coefficient is exactly the product of their normalized single-insertion coefficients. Consequently its first `d_logQ` derivative at Q=1 is zero. For two nonadjacent single-site marks with disjoint edge-port sets, every shared exterior component contains an occupied connection, so a nonzero activated pair coefficient requires **at least two distinct occupied components reaching both neighbourhoods**. This is a necessary connectivity condition, not a sufficient one or an arm-exponent assertion.

Read-only source: commit `2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb`, especially `notes/local-pair-two-insertion-algebra.md`, `notes/regular-pair-two-site-q-susceptibility.md`, and `notes/regular-pair-activation-original-u.md`. This proof uses no new occupation population or Monte Carlo calculation. The accompanying fixed Bell8 algebra and its independent verification are documented in RESULT.md.

## 1. Definition and exact factorization

Fix a finite graph and the original occupation A with both marked sites vacant. Remove their two prospective insertion tensors. All remaining local tensors are ordinary equality tensors at occupied sites and constant vacant tensors. Their connected equality components carry independent uniform colours in `{1,...,Q}` at integer Q>=4. Spectator components can be summed and cancelled against the unmarked colour weight.

Let `C_x,C_y` be the sets of exterior components incident to the two port groups. A component can touch several ports within either group. Set `s=|C_x intersection C_y|`. Let `k_x,k_y` be the canonical C4-averaged Kreg tensors evaluated on the component colours of their ports, including every repeated argument prescribed by that fixed connectivity.

The normalized coefficients are

```
b_x(Q,A) = E_col[k_x],
b_y(Q,A) = E_col[k_y],
c_xy(Q,A) = E_col[k_x k_y].
```

Distinct exterior components are allowed to have the same colour. The expectation is over all independent colour assignments, not injective assignments or assignments of different colours to different components.

If s=0, the colour variables in the two kernels are disjoint, giving `c_xy=b_x b_y` immediately.

If s=1, call the shared component colour a. Summing the colours private to each mark defines

```
f_x(a;Q) = Q^(-n_x) sum_(private x colours) k_x,
f_y(a;Q) = Q^(-n_y) sum_(private y colours) k_y.
```

Given a, the two private colour collections are independent, so

```
c_xy = Q^(-1) sum_a f_x(a;Q) f_y(a;Q).
```

The tensor Kreg is invariant under simultaneously permuting all colour labels. Applying a permutation carrying a to any other label is a bijection of every private-colour summation. Thus `f_x(a;Q)` is independent of a, and the same holds for y. Their common values are exactly b_x and b_y. Therefore

```
                       c_xy(Q,A) = b_x(Q,A) b_y(Q,A).       (1)
```

Multiple ports attached to the shared component do not alter this argument: they are repeated occurrences of the same a, and global colour permutations preserve all repetitions. Arbitrary connections among private ports likewise only identify summation variables within the corresponding factor.

## 2. Continuation to Q=1 does not use a singular conditioning operation

The equality-diagram expansion of Kreg has coefficients regular at Q=1. Its explicit unaveraged form is

```
Areg(a,b;c,d) = (1-delta_ab)(1-delta_cd)/2 *
 [delta_ac delta_bd + delta_ad delta_bc
  -(delta_ac+delta_ad+delta_bc+delta_bd)/(Q-2)
  +4/(Q(Q-2))].
```

Finite colour sums of these diagrams, normalized by Q powers, are rational functions regular near one. Equation (1) holds for every integer Q>=4 and therefore is an identity of their rational continuations. There is no literal conditioning on a nonexistent noninteger colour set. In particular, no factor `1/(Q-1)` has been introduced: the one shared-colour normalization is `1/Q`, regular at one.

On the single-colour set every delta is one, so either single insertion vanishes. Equivalently, the completed Bell4 table gives `b_x(1,A)=b_y(1,A)=0` for every exterior partition. Regularity and (1) imply

```
b_x(Q,A) = O(Q-1),       b_y(Q,A) = O(Q-1),
c_xy(Q,A) = O((Q-1)^2),
[d_logQ c_xy(Q,A)]_(Q=1) = 0.                             (2)
```

This argument uses the fixed canonical completion and needs no new contraction values. More generally, the proof applies to any colour-permutation invariant, entrywise-regular insertion whose every closed single coefficient vanishes at Q=1. Original-occupation scalar factors analytic at one, such as `Q^(-r(A)/2)`, cancel in the conditional colour ratio or multiply zeros and do not defeat (2). A rank assigned separately to virtual diagram joins would be a different definition and is not covered.

## 3. Conditional logarithm and occupation-summed spatial coefficient

With independent local parameters lambda_x,lambda_y, the relative fixed-occupation colour weight is

```
F_A = 1 + lambda_x b_x + lambda_y b_y + lambda_x lambda_y c_xy.
```

For s<=1, (1) makes `F_A=(1+lambda_x b_x)(1+lambda_y b_y)` exactly. Thus the mixed conditional logarithmic coefficient `c_xy-b_x b_y` is identically zero as a rational function of Q, not only at Q=1. In particular its first Q activation is zero.

Now let mu_Q be the normalized, analytic, finite-volume **original-occupation** law. This can depend on Q; no independence of the two neighbourhoods under mu_Q is assumed. The mixed coefficient of the logarithm after summing occupations is

```
G_xy(Q) = E_(mu_Q)[c_xy]
          -E_(mu_Q)[b_x] E_(mu_Q)[b_y].
```

For every occupation, all nonempty Kreg coefficients vanish at Q=1, including `c_xy(1,A)=0` even when s>=2. Consequently the derivative of the base law multiplies zero, and the derivative of the product of single means also vanishes. Hence

```
[d_logQ G_xy]_(Q=1)
 = E_(mu_1)[ a_xy(A) ],
a_xy(A) = [d_logQ c_xy(Q,A)]_(Q=1)
        = 0 whenever s(A)<=1.                             (3)
```

Equation (3) is a coefficientwise spatial support statement. Ordinary correlations between the single-mark occupations cannot bypass it at first order in Q-1; their product enters at least second order. For an unchanged original observable O, the corresponding normalized mixed response has source `Cov_(mu_1)(O,a_xy)`. Thermal differentiation and transmission to the original common-root U still require their specified covariance, root-motion and denominator terms; (3) does not evaluate those quantities.

## 4. Essential qualification on what counts as a bridge

The exact count s concerns **all exterior equality components of the port edge-nodes**, including a singleton edge-node not joined by an occupied site. If the marks are nonadjacent and their edge-port sets are disjoint, a component reaching both necessarily includes occupied equality vertices; shared components are then exactly the occupied exterior connections relevant to the proposed geometric formulation.

If two marks directly share an edge-node, that singleton is already one shared colour component even when no occupied path joins them. It must be counted. The theorem must not be restated for adjacent/overlapping port groups using only the number of occupied components. Likewise, if additional exterior tensors couple colours beyond ordinary equality components, the conditional-independence hypothesis must first be established rather than inferred from occupied connectivity alone.

At least two shared components is necessary, not sufficient. Their detailed port incidences and the kernel contraction can still cancel the activation. The existing four-path witness has four distinct shared components and the already proved derivative 13/8, fully consistent with this lower support threshold; it has not been recomputed here.

Finally, distinct components reaching both neighbourhoods are a finite-graph connectivity event. Calling this a specific planar alternating-arm event, assigning a square-lattice arm exponent, establishing decay with separation, or identifying a CFT field would require further geometric and probabilistic arguments not supplied by the lemma. Singular projectors/normalizers at Q=1, fixed-colour symmetry-breaking boundaries, uncontinued permutation seams, and taking infinite-volume limits before Q=1 are outside its hypotheses.

## 5. Annealed bound using the independently checked canonical maximum

Include the vacancy indicators in the definition of `a_xy`: it is zero unless both marks are vacant. The separately checked canonical Bell8 table establishes

```
max_(exterior port partitions) |a_xy| = 43/16.
```

The colour-coarsening computation and the independent equality-diagram calculation agree on all 4140 entries; see [the comparison receipt](review/COMPARISON.json). With that checked input, the support lemma yields the pointwise inequality

```
|a_xy(A)| <= (43/16) 1_B(A),
B = {x,y vacant and at least two shared exterior components}.
```

For disjoint edge-port sets at nonadjacent sites the shared components in B are occupied exterior components. The Q=1 original law is a positive Bernoulli measure, so the triangle inequality gives exactly

```
|C_xy| = |E_(mu_1)[a_xy]| <= E_(mu_1)[|a_xy|]
       <= (43/16) Pr_(mu_1)(B).                            (4)
```

This is an unconditional probability including the two vacancy events, not a conditional probability given their vacancy. No positivity of a_xy, independence of the two neighbourhoods, or positivity of the analytically continued Q!=1 interaction is used. The coefficient C_xy here is precisely the first Q-activated mixed logarithmic partition coefficient in (3).

Equation (4) does not immediately provide the same bound on the original-U transmission: the latter uses covariances with the source and their thermal derivatives, then divides by the common-root slope. Those operations are specified next.

## 6. The original-U mixed response is exactly the existing W[a_xy]

Use the same two geometries, the same original rank observables `q(A),E(A)=q(A)^2`, the same Delta and the same fixed size factor A_N. Match the two marked pairs between geometries as stipulated by the source definition. The source `a_xy` denotes that pair of geometry-specific functions; it is not silently averaged over sites or over pairs.

Let `eta=log Q`. For the two-mark family, the coefficientwise Q=1 zero theorem gives, for every original occupation and thermal activity h,

```
F_A(1;lambda_x,lambda_y)=1,
partial_eta log F_A|_(Q=1)
 = lambda_x a_x(A) + lambda_y a_y(A)
   + lambda_x lambda_y a_xy(A).                           (5)
```

In particular the Q=1 normalized occupation law, all original q/E means, and the selected common root are independent of both lambdas. Differentiating the complete original-U functional once in eta is therefore its first, linear variation at this same base law. The base Q dependence adds a lambda-independent tangent and disappears under `partial_lambda_x partial_lambda_y`. The product of the single-insertion tangents cannot contribute at this order: the undifferentiated single coefficients vanish at Q=1, and the entire source-induced change of the base law starts at first order in eta. Hence

```
partial_lambda_x partial_lambda_y partial_logQ U
 |_(Q=1,lambda_x=lambda_y=0)
 = W[a_xy].                                               (6)
```

Here W is exactly the existing linear original-U response, with no new observer or numerical score. In the existing h chart define at the unchanged Q=1 common root h0

```
M = (<q>_1+<q>_2)/2 = 0,
Y = (<E>_1-<E>_2)/Delta,
D = M_h,                 R = Y_h/D,                 U = A_N R,
jM = [Cov_1(q,a_xy^(1))+Cov_2(q,a_xy^(2))]/2,
jY = [Cov_1(E,a_xy^(1))-Cov_2(E,a_xy^(2))]/Delta.
```

The superscripts on `a_xy^(1)` and `a_xy^(2)` label geometries. With thermal derivatives taken before evaluation at h0,

```
W[a_xy]/A_N
 = jY_h/D                    # normalized direct thermal response
   -Y_hh*jM/D^2              # common-root motion
   -R*jM_h/D                 # source change of the denominator
   +R*M_hh*jM/D^2,           # root change of the denominator

partial_lambda_x partial_lambda_y partial_logQ h0
 = -jM/D.                                                   (7)
```

All four terms remain necessary. Covariance, rather than a raw `<O a_xy>` moment, accounts for each geometry's own normalization. No term involving `a_x*a_y` needs to be added to (7). Independent lambdas carry no factor of two; if a separately stipulated common coupling is substituted later, its ordinary chain-rule combinatorial factor must be retained.

The analytic domain is a finite graph, an interior thermal root, finite equality-diagram coefficients regular near Q=1, nonzero partition normalizers, and a simple chosen common root `D!=0`. The implicit-function theorem then supplies a joint analytic root branch near `(Q,lambda_x,lambda_y)=(1,0,0)`, so the derivatives commute. At Q=1 the Bernoulli law is positive; positivity of the continued interaction away from one is not required. These are local finite-volume statements: no uniform-in-volume neighbourhood, infinite-volume interchange, global uniqueness for all couplings, or missing numerical U value is asserted.

## 7. Which statements survive a different regular completion

The exact factorization `c_xy=b_x b_y` for at most one shared component uses colour-permutation invariance and the exterior equality-component structure. It applies to every such insertion in a coherent finite equality-diagram continuation, not only the canonical Kreg. The resulting **raw first-Q selection rule** additionally uses `b_x(1)=b_y(1)=0`.

Every entrywise-regular completion retaining the original unequal-pair support has that Q=1 zero. In particular all regular `K2+c(Q)K0` completions with `c(1)=1` retain the rule: writing `c=1+(Q-1)a(Q)` adds the regular tensor `a(Q)*average[(1-delta_ab)(1-delta_cd)]/Q`, which still vanishes on the one-colour set. Analytic coefficient functions cause no continuation problem: expand the finite diagram contractions first, where the factorization is a formal diagram identity, then multiply by those analytic coefficients.

Regularity and S_Q invariance **alone**, with no Q=1-zero/support requirement, are insufficient for the raw-coefficient statement. For example `Ktilde=Kreg+1` is regular and invariant but has `btilde(1)=1`. At two colour-disjoint all-free marks,

```
ctilde(Q) = [1+(Q-1)/Q^3]^2,
ctilde'(1) = 2.
```

It still has an exactly factorized conditional logarithm, but it changes the Q=1 interaction and is not an endpoint-invisible completion of the stipulated family. This is why the zero hypothesis must not be dropped when stating (2), (3), or (6).

The support selection rule is therefore completion-independent **within the endpoint-invisible regular invariant class**. The constant 43/16, and any ensuing numerical transmission, belong to the fixed canonical completion only; a different admissible counterterm can change its activation coefficients and their maximum.

## 8. Why the finite colour-coarsening calculation is exact

Let pi be the exterior partition of the eight ports and b its block count.
The independent colours on its b components may coincide. Their complete
colour-equality partition rho therefore runs over all coarsenings of pi.
If rho has k blocks, its multiplicity is the falling factorial `(Q)_k`.
For the fixed canonical rational kernel this gives exactly

```
c_xy(Q;pi) = Q^(-b) sum_(rho >= pi) (Q)_k
                         Kreg(Q;rho_x) Kreg(Q;rho_y).
```

At k=1 both kernels vanish identically in Q: all four arguments have
one colour and retain the unequal-pair prefactor. For k>=2, `(Q)_k`
has a simple zero at Q=1 and derivative
`(-1)^(k-2)*(k-2)!`. Thus neither the kernel derivatives nor the
normalizing factor contributes to the first derivative, and

```
a_xy(pi) = sum_(rho >= pi, k>=2) (-1)^(k-2)*(k-2)!
                                 Kreg(1;rho_x) Kreg(1;rho_y).
```

Each `4*Kreg(1;rho_x)` is an integer. The main computation therefore
accumulates `16*a_xy` as an integer, without finite-Q interpolation,
roundoff, or a limit estimated from nearby Q values. Enumerating every
restricted growth string of length eight gives 4140 pi; enumerating the
partitions of each pi's block set gives 167894 total pairs `(pi,rho)`.

The independent algorithm instead expands the rational kernel into
delta diagrams, joins them to pi, and differentiates both diagram
coefficients and their powers of Q. It does not use coarsenings or
falling factorials. Its exact agreement and symmetry controls are
recorded in the companion result. Exhaustiveness is over these eight
ports only, with the originally fixed kernel; it introduces no family
of new sources or physical configurations.
