# Fixed local pair insertion: what its size response can decide

**Derived predictions, not a new score.** Execution's `branch_only`
[`923f66b9` result](https://github.com/LightChainr/Matching-One/blob/923f66b979a6b6132875f783106c041ed3c0c1a9/notes/local-four-port-transmission-result.md)
already gives the N25 site-average response
`V_av(25)=+0.0018155512845251097`. Its local tensor is configurationwise
different from the full seam trace. This note preserves that tensor and
derives one concrete size comparison; it performs no enumeration,
Monte Carlo, root search, source fit or tests.

The immediate quantitative question is whether the **extensive response**
`W_N=N V_av(N)` grows like `sqrt(N)` or stays of order one on a fixed
homothetic geometry family. Under the explicitly stated single-field
loading assumptions these are the different predictions of dimensions
17/4 and 21/4. A constant W alone does not identify the latter field:
changing an existing anisotropy amplitude can give the same prediction.

## 1. Exact units and an intrinsic form of the old U

Keep the frozen Q1 mark

```text
t_x = I_(NS|EW) + (I_(NE|SW)+I_(NW|ES))/2,
S_av = -sum_x t_x/N,       S_ext = -sum_x t_x.
```

The site-average coefficient epsilon changes each vertex by epsilon/N;
an extensive coefficient g changes each vertex by g. At first order,
epsilon=N g. Linearity of the entire response, including its moving root
and both slope terms, gives the exact finite-N relation

\[
 W_N:=\partial_g U_N=N\,\partial_\epsilon U_N=N V_{av}(N).       \tag{1}
\]

Every one of the four reported response terms and the root tangent has
this same factor N. Translation invariance permits one fixed origin for
these first-source moments; it does not replace a two-insertion moment.
Equation (1) requires no continuum limit or scaling hypothesis.

For an equal-area direction pair, write

\[
 M=\tfrac12(\langle q\rangle_a+\langle q\rangle_b),\qquad
 Y=\frac{\langle E\rangle_a-\langle E\rangle_b}{\Delta_4},\qquad
 A_N=\frac{N^{13/8}}2,\quad
 U_N=A_N\frac{Y_p}{M_p}\bigg|_{M=0},                         \tag{2}
\]

where `q=r-1`, `E=q²`, and `Delta4=cos(4 theta_a)-cos(4 theta_b)`.
Use the same thermal coordinate in both geometries; p may be replaced
by h=p/(1-p) without changing this ratio or its complete source response.

Whenever `D=M_p` is nonzero, let `z=M(p,g)` and
`X_N(z,g)=Y(p(z,g),g)`. Thus

\[
 U_N=A_N\partial_z X_N(0,g).
\]

For the extensive source, put `jM=partial_g M`, `jY=partial_g Y`,
with separately normalized covariance responses in each geometry. Its
intrinsic, fixed-z deformation is exactly

\[
 \partial_g X_N=jY-\frac{Y_p}{M_p}jM,\qquad
 \boxed{\frac{W_N}{A_N}
  =\frac1{M_p}\partial_p\left[jY-\frac{Y_p}{M_p}jM\right]_{M=0}.} \tag{3}
\]

Expanding (3) gives the four terms in the completed interface. This
identity makes clear what the response measures: the slope of the
directional E-versus-pooled-q curve, after removing the common thermal
motion. It is not the root shift alone or an uncentered local density.

### Exact thermal redundancy is already separated

If both geometries' entire q/E curves change only through a common
thermal reparameterization `p -> p+g c_N(p)`, then
`jM=c_N M_p`, `jY=c_N Y_p`. The bracket in (3) vanishes identically,
so `W_N=V_av(N)=0` at every finite N. Geometry-dependent normalization
constants have zero covariance and do not alter this statement.
In particular, a common K source gives this exact null.

The positive N25 result already excludes that exact finite-family
redundancy for the fixed local source. It does **not** exclude a
microscopic perturbation whose dominant long-distance component is
thermal but whose subleading components change anisotropy. Subtracting
any common K column would leave W_N unchanged; it cannot create an
additional mechanism discriminator.

## 2. Derivation of the fixed power comparison

Set `L=sqrt(N)` in the original lattice units. The repository's dimension
convention for an integrated local field is `g L^(2-x)`; its thermal
exponent is `y_t=3/4`. These are the conventions used in
[`operator-mixing-identifiability-boundary.md`](operator-mixing-identifiability-boundary.md),
not an exponent inferred from the N25 response.

Assume the following specific mechanism on a fixed-shape, fixed-angle
homothetic family:

1. The fixed tensor has an N-independent, nonzero coupling to one leading
   field of dimension x in the directional channel under consideration.
2. After the common thermal motion has been removed, its intrinsic
   response has `partial_g X_N(z)=a L^(2-x) Phi(z)+o(L^(2-x))`
   near z=0, with a controlled differentiated remainder.
3. `Phi'(0)` is nonzero. This is an observer/parity loading assumption,
   not a consequence of the field existing in the spectrum.
4. A same-order logarithmic partner, another leading field, or a changing
   microscopic normalization is not part of this single-power claim.

Since the p derivative in the numerator and M_p in the denominator both
bring the same thermal factor `L^(3/4)`, it cancels. Equations (1)--(3)
then give

\[
 \boxed{W_N\sim \frac{a\Phi'(0)}2 N^{21/8-x/2},\qquad
 V_{av}(N)\sim\frac{a\Phi'(0)}2 N^{13/8-x/2}.}                \tag{4}
\]

The `13/8` prefactor is part of the observable definition, not proof that
the inserted field has dimension 21/4. In particular, the normalization
adds one extra power of N when moving from the reported site-average
response to the physically fixed per-vertex perturbation.

| Candidate, conditional on the nonzero U loading above | W_N | V_av(N) | W_(4N)/W_N | V_av(4N)/V_av(N) |
|---|---|---|---|---|
| Four-leg spin +/-4, x=17/4 | N^(1/2) | N^(-1/2) | 2 | 1/2 |
| Thermal level-4 spin +/-4, x=21/4 | N^0 | N^(-1) | 1 | 1/4 |
| A spin-4 field with x=4, **only if this U slope is allowed** | N^(5/8) | N^(-3/8) | 2^(5/4) | 2^(-3/4) |

The exact dimensions 17/4 and 21/4 are documented in
[`q1-spin4-competitor-preflight.md`](q1-spin4-competitor-preflight.md)
and [`thermal-level4-spin4-candidate.md`](thermal-level4-spin4-candidate.md).
The familiar identity-family x=4 assignment is matching-even. Under that
parity hypothesis its leading E scaling function is even in the thermal
coordinate, so `Phi'(0)=0`: the last row is then **absent**, not a
prediction of N^(5/8) growth in U. More generally, a forbidden or vanishing
slope invalidates a row's use as a field exclusion; it does not change x.
The new tensor's C4 averaging has not supplied a pure matching-parity or
continuum-spin projection.

### One concrete dimensionless comparison

Preserve the original geometry family

```text
axis: (5k,0),       tilted: (4k,3k),       N=25 k²,
Delta4=1152/625,    continuum modulus tau=i.
```

For integer k and 2k, use exactly the same local tensor and all four
response terms. If the reference response is nonzero, score

\[
 \mathcal R(k)=\frac{W_{100k^2}}{W_{25k^2}}
             =4\frac{V_{av}(100k^2)}{V_{av}(25k^2)}.          \tag{5}
\]

The two principal predictions are `R -> 2` for an allowed x=17/4 loading
and `R -> 1` for an allowed x=21/4 loading. No free radial exponent or
post-result rescaling is needed. More generally an area multiplier r
gives the dimensionless fixed-model residual

\[
 \mathcal D_x(N,r)
 =r^{x/2-21/8}\frac{W_{rN}}{W_N}-1\longrightarrow0.           \tag{6}
\]

N25 is a completed finite anchor, not an established scaling window.
A ratio at one microscopic dilation cannot prove an asymptotic exponent;
a sign change or failed ratio rejects the declared single-leading-field
description of that window. It does not exclude all occurrences of the
field in larger systems or other observers. An assumed logarithmic
factor would instead multiply a ratio by
`[a+b log(sqrt(rN))]/[a+b log(sqrt(N))]`; fitting it after any failed
two-size comparison would remove the intended distinction.

## 3. What same-area axis/tilted quotients do and do not isolate

Both tori in (5) have the same area and continuum modulus, so a leading
rotation-scalar response cancels from Y. Their finite translation groups
are nevertheless different: the axis quotient has Smith factors
`(5k,5k)`, while the tilted quotient has `(k,25k)`. This is a paired
orientation comparison, not an exact graph isomorphism or a licence to
reuse one canonical cyclic permutation on both graphs. These Smith
classes stay fixed along k-dilation; short wrapping/contact corrections
at k=1 do not establish the common continuum limit's asymptotic accuracy.

The source is a reflection-even C4 invariant. Such an invariant permits
continuum spins 0, +/-4, +/-8, ... . With only these two orientations,
every allowed harmonic contributes through the scalar number

\[
 H_s=\frac{\cos(s\theta_a)-\cos(s\theta_b)}
           {\cos(4\theta_a)-\cos(4\theta_b)}.
\]

H4 is one, but H8 and higher factors need not vanish. Keeping the same
two angles at more sizes adds radial information, not angular rank.
Therefore (5) can compare specified dimension/loading mechanisms; it
cannot certify pure spin 4, distinguish fields with the same dimension,
or establish a local four-leg primary from the pair-representation name.

A common Gaussian multiplication by `1+i` rotates both directions by
pi/4 while doubling area. For a pure H4 term the raw directional
difference changes sign, **and its exact Delta4 denominator also changes
sign**. The normalized V and W follow (4) without an additional minus
sign. Importing the raw matching Gaussian-doubling sign into this
normalized response would therefore be incorrect.

### Why a constant W does not eliminate scalar mediation

A C4-scalar microscopic perturbation may change an existing anisotropy
coupling: `u_j(g)=u_j(0)+g u_j'(0)+...`. If the baseline directional
contribution has dimension x_j, that change gives
`W_N proportional to N^(21/8-x_j/2)`. For x_j=21/4 it is constant,
just like a direct field of that dimension. This is a genuine change of
the intrinsic E-versus-q curve and is not the exact thermal null.

If a source is instead a pure scalar RG eigenfield of dimension x_s
with no such direct anisotropic coefficient, its leading directional
effect can be a product with a pre-existing spin sector x_j. The
conditional product scaling is
`W_N proportional to N^(29/8-(x_s+x_j)/2)`, obtained by multiplying
`L^(2-x_s) L^(2-x_j)` before the old A_N normalization. Neither its
coefficient nor its nonzero U slope has been established here. Thus
the proposed ratio can disfavor a **specified direct x=17/4-dominant
mechanism**; an order-one result cannot by itself choose direct thermal
Q4 over modulation of existing anisotropy.

## 4. The two-insertion interface has a specific obstruction

A genuine separated two-insertion contraction of this same tensor could
add information that the integrated response lacks. It is the connected
coefficient `Z_xy/Z - Z_x Z_y/Z²` with Kbar at two distinct vacant
vertices, contracted jointly in the physical colour/port network before
Q1 continuation. The one-insertion proof alone does not establish its
removable Q1 limit or identify it with `Cov(t_x,t_y)`. In particular,
the two vertices can share exterior colour components, so multiplying
their separately closed marks is not a justified tensor contraction.

If that joint local contraction exists and loads a nonzero ordinary
two-point field, the plane regime `1 << r << sqrt(N)` has radial power
`r^(-2x)`. The scalar thermal field x=5/4 gives `r^(-5/2)`; the two
spin-4 candidates give `r^(-17/2)` and `r^(-21/2)`. A real +/-s pair has
the same-spin angular factor `cos(2s phi)` under the usual diagonal
two-point selection, so a real spin-4 primary gives an eighth angular
harmonic, not automatically a fourth one. Axis versus diagonal
separations alias this factor. At separations of order the torus length,
the plane formula must be replaced by the fixed torus scaling function.
Zero norms/logarithmic partners can change the usable correlator and add
logs; a vanishing ordinary two-point function alone does not eliminate
a field in this c=0 problem.

The subsequent [exact two-hole closure](local-pair-two-insertion-obstruction.md)
now shows that this unrenormalized tensor has a Q1 simple pole, with
Gram residue1/2, in a realizable four-path exterior. The connected
conditional response retains it. Thus the preceding plane two-point
powers cannot be assigned directly to an already regular Q1 Kbar field;
an explicit finite/confluent completion is needed for that interpretation.
This does not invalidate the finite occupation tangent or its response.

Consequently the next defined linear discriminator is
(5), with unchanged first-source moments and tensor. A local-field
identification stronger than its conditional radial conclusion would
require this **specific joint two-insertion closure**, not a covariance
of separately closed marks, another N25 mixture, or a generic finite
interface. The linked two-insertion result excludes the unrenormalized
regular conditional-family option; it does not select its replacement.

## Source pins and lifecycle

- Execution snapshot `7681eedd938019d977ede41a7d74ee1b88ffbc50`,
  `branch_only`: `notes/local-four-port-original-u-interface.md` and
  `analysis/p337_local_pair_insertion_contract.json`; fixed source,
  average/extensive units and complete original-U functional.
- `923f66b979a6b6132875f783106c041ed3c0c1a9`, `branch_only`:
  `notes/local-four-port-transmission-result.md`; completed positive
  N25 result, not an unrun request.
- Kernel `3c3fe12f4e3b6212797980b2981d0bd5506d2c07`:
  `notes/closed-source-local-four-port-pair-kernel.md`; C4 recoupling,
  one-insertion-only removable continuation and lack of pure-cut identity.
- Topology `ab402605f54d5aa0f4e2209dbad715cb9fb159da`:
  `notes/local-four-port-pair-insertion.md`; physical mark, old-rank
  convention and strict separation from the global seam.
- `6bd46ad30bc8f583c3ca1f1c8a1b95e7d90571bc`:
  `notes/norm4-common-source-response-determinants.md`; actual U,
  source units and common thermal null.
- `a38f76b903644434fe438211264c7c8a23aa141c`:
  `notes/operator-mixing-identifiability-boundary.md`; integrated-field
  scaling, y_t and observer/parity loading boundary.
- `d4d101311f1b850656849c8ba9a4f30338e2fe00`:
  `notes/q1-spin4-competitor-preflight.md`; x=17/4 and distinct overlap gates.
- `ba1ca3097175dc441df61f57c98ad6e34d8e2984`:
  `notes/thermal-level4-spin4-candidate.md`; x=21/4 and the thermal exponent.

Equations (1)--(3) are exact finite-family identities. Equations (4)--(6)
are conditional scaling predictions, not measured exponents. No new
source, exact population, random block, score, test or cloud job was run.
