# N25 admits an explicit stable colour completion

**Completed evaluation:** [5c1f9d3b](../results/n25-stable-colour-q1/REPORT.md)
now gives `B=−0.001904836180602413` and
`B_logQ=+0.005036496028411871`, both with exact enclosures excluding zero.
The same saved counts and Q1 root suffice; no new enumeration or root
search was performed. The proof and the pre-score definition follow.

**Result.** For the two scored N25 tori `(5,0)` and `(4,3)`, every
rank-one occupation configuration with `c` essential components and
first-period winding coordinate `u` satisfies `c |u| <= 2`. Consequently
its full S4 `[22]` seam contraction already lies in the stable range of
the generic `[Q-2,2]` character calculation. The completion below agrees
with the completed Q4 insertion and has an explicit Q1 jet, recoverable
from the existing per-component seam counts.

This is a geometry-specific positive result. The
[finite-colour counterexamples](colour-specialization-gap.md), committed
as `f81f70c8866e9165c2c9a3c0802b270e66133c24` and calculated in
`ef7c3f58b9bbe2e21b1b4a9f5d8f57688ef121c2`, remain valid for larger
tori. Their failure of automatic specialization does not occur in this
N25 packet. Conversely, agreement here does not equate a full isotypic
trace with a particular local four-leg field.

The [completed Q4 score](https://github.com/LightChainr/Matching-One/blob/54352b2eefa651ca482ca84837053c792e82c71e/results/p337-s4-trace-transmission/score/score.json)
is an input, not recalculated here. This note is a proof and a definition
of the next derivative object; it reports no new score, enumeration or
test.

## 1. The closure and its stable range

Use the honest square-cell quotient and the first-period colour seam of
the [homology-decorated closure](https://github.com/LightChainr/Matching-One/blob/4a4390f2aeff3e79855fb11ef1185ba52c80a43a/notes/closed-source-finite-torus-pair-closure.md).
The period columns are `(a,b),(-b,a)`. A rank-one configuration has `c>=1`
essential occupied NN components with a common primitive homology
direction `(u,v)`, up to sign, and `c0` zero-image hypergraph components.
The seam class function, after removing common factors, is

```text
F_(u,c)(pi) = Fix(pi^u)^c.
```

For `u!=0`, write `X_j` for the number of j-cycles of pi. Then

```text
Fix(pi^u) = sum_(j divides |u|) j X_j,
chi_[Q-2,2] = (X_1)_2/2 + X_2 - X_1,
d_[Q-2,2](Q) = Q(Q-3)/2.
```

At integer Q>=4 the full central coefficient is `d(Q) m_Q(u,c)`, where
`m_Q=<chi,F>` is the normalized S_Q class sum. It multiplies the physical
factor `v_act^K Q^(c0-1/2)`; `v_act` denotes hyperedge activity here, not
the second winding coordinate.

For a uniform permutation, the exact factorial-moment rule is

```text
E prod_j (X_j)_(n_j) = prod_j j^(-n_j), if sum_j j n_j <= Q;
                      0, otherwise.
```

Assign weighted degree j to X_j. The product `chi F` has weighted
degree at most `|u| c+2`. Its factorial moments therefore equal their
stable values whenever

```text
Q >= |u| c+2, with Q>=4.
```

This supplies an explicit finite-polynomial stabilization argument,
not an inference that infinitely many integer values determine an
arbitrary analytic function. If `u=0`, `F=Q^c` is constant in pi and its
nontrivial central projection is identically zero; no bound on c is
needed for that case.

## 2. Packing lemma: each black essential component costs a white separator

Let `w=P(u,v)` be the lifted primitive displacement of the common
essential direction. For c distinct essential black NN components,

\[
\boxed{N\ \ge\ c\bigl(\|w\|_1+\|w\|_\infty\bigr).}
\tag{1}
\]

Here the white graph is the vacant-site matching graph with NN and
diagonal adjacency. The following argument establishes the component
count in (1), not only the existence of one complementary wrapping
component.

### Black core curves and the intervening annuli

Choose one embedded essential simple NN cycle in each essential black
component. Such a cycle exists by decomposing a nonzero cycle into
simple cycles. An essential simple curve on a torus is primitive.
Disjoint essential simple curves have zero intersection number, so the
c chosen curves have the same primitive direction, with orientation
chosen consistently. Call them `gamma_1,...,gamma_c` in their cyclic
transverse order.

Cutting the torus along these curves leaves c annuli. When c=1 this
means one annulus with the two boundary copies of the same gamma; it
does not mean that the complement is a disk. For c>=2 the boundaries
of each intervening annulus belong to consecutive distinct black
components.

There is no black NN crossing between the two boundaries of any of
these cut annuli. For c>=2 such a crossing would connect two components
declared distinct. For c=1, a crossing joining the two boundary copies,
closed by an arc on gamma after regluing, would have transverse
intersection number one with gamma. The black ambient image would then
have rank two, contradicting the rank-one hypothesis. Boundary contacts
are interpreted in a thin regular neighbourhood, so paths running along
one boundary are not transverse crossings.

### White separators in each annulus

Use the [digital Alexander lattice bridge](digital-alexander-duality-proof.md):
a regular neighbourhood U of the black NN graph has a complementary
subsurface whose embedded white 1-skeleton is obtained from the matching
graph. Retain the diagonal in an opposite-white-pair face; replace
redundant matching diagonals by white face-boundary paths. This embedded
reduction preserves ambient homology and has only white vertices.
It applies componentwise in the cut annuli.

In one annulus take the black neighbourhood connected to its first
boundary. It cannot reach the second boundary, by the preceding
no-crossing argument. The boundary of that neighbourhood facing the
second boundary contains an essential separating circle. Equivalently,
this is annular black-crossing/white-separator duality: contractible
boundary circles alone cannot separate the two boundary components of
an annulus. Push the separating circle into the white complementary
subsurface and use its embedded white 1-skeleton to represent its
nonzero homology class. Decomposing that white closed walk into simple
cycles yields an essential simple white matching cycle in the same
annulus.

Its primitive class is `(u,v)` up to sign. An essential curve in that
annulus is parallel to its boundary cores; it cannot acquire an
independent transverse direction. The c annular interiors are disjoint,
and a white vertex cannot lie on a black core used as a cut. The c
separator cycles obtained this way therefore use disjoint sets of white
vertices, even if extra contractible white components or black islands
are present. This argument also gives the required single separator
when c=1.

### Counting vertices

Each of the c vertex-disjoint black simple cycles has a lift with
displacement w and NN steps, hence at least `||w||_1` edges and vertices.
Each of the c white simple separator cycles has the same primitive
displacement and steps whose coordinate increments have absolute value
at most one, hence at least `||w||_infty` vertices. Black and white
vertex sets are disjoint. Their sum is bounded by N, proving (1).

The proof uses an honest embedded square-cell torus and complementary
4/8 adjacency. It does not silently cover short-period local aliases,
arbitrary abstract graphs or a different vacant-site adjacency.

## 3. The two N25 bounds

For the axis period `(5,0)`,

```text
w=(5u,5v),
||w||_1+||w||_infty
  =5(|u|+|v|+max(|u|,|v|)) >= 10|u|.
```

Equation (1) gives `25>=10c|u|`. Since `c|u|` is an integer,

```text
c|u| <= 2.
```

For the tilted period `(4,3)`, write

```text
w=(x,y)=(4u-3v,3u+4v),     25u=4x+3y.
```

For all real x,y,

```text
|4x+3y| <= 4|x|+3|y|
         <= (7/3)(|x|+|y|+max(|x|,|y|)).
```

The last inequality follows separately from `|x|>=|y|` and its reverse;
the respective right-minus-left differences are
`(2/3)(|x|-|y|)` and `(5/3)(|y|-|x|)`. Therefore (1) implies

```text
25 c|u| <= (7/3)c(||w||_1+||w||_infty) <= 175/3,
c|u| <= 7/3,
```

and again `c|u|<=2` by integrality. In both geometries every nonconstant
seam class function is thus stable already at Q4. In particular, neither
`(|u|,c)=(1,3)` nor `(3,1)` from the larger-torus counterexamples can
occur in this packet.

## 4. The full stable insertion and its saved-count labels

For `u!=0` the only possibilities are `(1,1)`, `(1,2)` and `(2,1)`.
The full central multiplicities are

| `(|u|,c)` | Seam class function | Stable multiplicity |
|---|---|---:|
| `(1,1)` | `X_1` | 0 |
| `(1,2)` | `X_1^2` | 1 |
| `(2,1)` | `X_1+2X_2` | 1 |

These follow either from the factorial moments in section1 or the
decomposition of the point and two-point colour carriers. Rank0 is
constant in the seam. Rank2 has one full-rank essential component and
the point-colour character, so both also have zero `[Q-2,2]` contraction.

Define geometry-dependent configuration indicators

```text
I12 = 1[rank=1, |u|=1, c=2],
I21 = 1[rank=1, |u|=2, c=1].
```

After dividing by the untwisted colour factor Q^c, the complete stable
central insertion relative to the original occupation weight is

\[
\boxed{\beta(A;Q)
 =I_{12}(A)\frac{Q-3}{2Q}
  +I_{21}(A)\frac{Q-3}{2}.}
\tag{2}
\]

This defines the positive-real-Q continuation of the full stable
isotypic trace for these two geometries, including its multiplicity
contraction. It does not merely attach d(Q) to an unspecified finite
colour label. At Q4 its two nonzero values are 1/8 and1/2, exactly the
canonical S4 insertion already scored.

The saved seam statistics are sufficient. Let b2 indicate that some
component has a nonzero first-coordinate homology image modulo2, and
let n3 count components whose first-coordinate image is nonzero
modulo3. The N25 bound proves the exact identifications

```text
I12 = 1[rank=1, b2=1, n3=2],
I21 = 1[rank=1, b2=0, n3=1].
```

The rank filter is retained explicitly. No inference about unseen
higher winding or higher component count is made from these compressed
labels: the packing theorem has excluded those possibilities first.
The existing counts also retain K and g, so they can be regrouped and
reweighted without recovering individual occupation configurations.

Equation (2) gives the explicit operator jets

\[
\boxed{\beta(A;1)=-I_{12}-I_{21},\qquad
\partial_{\log Q}\beta(A;Q)|_{1}
 =\tfrac32I_{12}+\tfrac12I_{21}.}
\tag{3}
\]

These are signed colour-sector attributions, not Q1 probabilities or
literal one-colour irreducible dimensions. The original real-Q
occupation family remains positive. Since beta is bounded near Q1,
the bookkeeping insertion below also has positive weights for
sufficiently small absolute epsilon. It is a specified trace-coefficient
perturbation, not an identification of a local field or a new fitted
physical source.

## 5. The derivative object to compute from the existing counts

Use `Q=m^2`, `h=y/m`, and the same g stored in the closed-source packet.
Up to an observer-independent common factor, the original occupation
weight is

```text
w_0(A;h,Q)=h^K Q^(-g/2).
```

Define one two-parameter family, with beta fixed by (2):

\[
\boxed{w_{\epsilon,Q}(A;h)
 =h^{K(A)}Q^{-g(A)/2}[1+\epsilon\beta(A;Q)].}
\tag{4}
\]

Normalize (4) separately in each geometry. Write its expectations as
`q_j(h,epsilon,Q)` and `E_j(h,epsilon,Q)`, j=axis,tilted, and set

```text
M=(q_axis+q_tilted)/2,
Delta4=cos(4 theta_axis)-cos(4 theta_tilted)=1152/625,
Y=(E_axis-E_tilted)/Delta4,
A_N=25^(13/8)/2.
```

Let `h0(epsilon,Q)` be the continuation of the original regular pooled
root `M=0`. The entire original observable, including its slope, is

\[
\mathcal U(\epsilon,Q)
 =A_N\left.\frac{\partial_hY}{\partial_hM}
       \right|_{h=h_0(\epsilon,Q)}.
\tag{5}
\]

The two requested outputs are specifically

\[
\boxed{B=\partial_\epsilon\mathcal U(0,1),\qquad
 B_{\log Q}=\left.\partial_{\log Q}\partial_\epsilon
                         \mathcal U(\epsilon,Q)\right|_{\epsilon=0,Q=1}.}
\tag{6}

The second object differentiates the trace response, including beta's
explicit Q dependence. It is not `partial_logQ U(0,Q)`, not the total
closed-source tangent, and not a local four-leg matrix element. The
finite sums and the regular-root implicit-function theorem make these
local derivatives well-defined; the two differentiation orders agree.

Because beta is supported on rank1, its q and E unnormalized
numerators vanish for every h,Q. If `f_j=<beta>_(j,epsilon=0)`, then

```text
q_j(epsilon)=q_j(0)/(1+epsilon f_j),
E_j(epsilon)=E_j(0)/(1+epsilon f_j).
```

This gives the same normalization-only transmission mechanism as the
completed Q4 score, but at a different specified colour point and with
the mixed derivative in (6). Its root and thermal-denominator terms
must not be dropped. The moving pooled root need not be an individual
zero of either q_j.

For computing the log-Q derivative with the saved g,K counts, the
baseline score at fixed h is

```text
S_h = partial_logQ log w_0 |_h = -g/2.
```

In particular, at fixed h and Q1,

```text
partial_logQ f_j
 = E_j[ (3/2)I12+(1/2)I21 ]
   + Cov_j( -I12-I21, -g/2 ).
```

This is a fractional-trace jet at fixed thermal coordinate; equation(6)
still includes the Q motion of the baseline q/E curves, the pooled root
and all thermal derivatives. It is not obtained by retaining this one
covariance term alone.

At fixed original odds y, `h=y Q^(-1/2)`, so the baseline score is

```text
S_y = -(g+K)/2 = S_h-K/2.
```

The difference is the common K tilt corresponding to that thermal
reparameterization. At every fixed epsilon,Q, both `Y_y` and `M_y`
acquire the same factor `Q^(-1/2)`, and the corresponding roots obey
`y0=Q^(1/2) h0`. Thus (5), and hence both complete derivatives in (6),
are unchanged by using y in place of h. Dropping the moving-root or
slope response would destroy this coordinate-invariance check.

The completed calculation uses the existing exact counts once with
these fixed derivative objects. It took6.759455seconds locally, with no
Q4 rescore, root search or new seam counts. The two channel contributions
to B are `−0.001945570733316785` from I12 and
`+0.00004073455271437206` from I21. Hence a two-distinct-essential-cluster
restriction omits a real, oppositely signed part of this full trace.
These exact terms are not fitted fractions or independent evidence.

The explicit beta(Q) derivative contributes `+0.0028979888236179917`
to B_logQ; the remaining measure/root/slope contribution is
`+0.002138507204793879`. The complete derivative includes both.
All exact coefficient jets and rational intervals are saved in
[the result packet](../results/n25-stable-colour-q1/latest.json).

An interpretation as a particular local continuum operator still requires
separate port, spectator and multiplicity intertwiners; the present
completion does not assume that further identity. In particular B is
already nonzero at Q1: this is not a baseline-zero activation established
solely by differentiating a projector. Nor do the different Q1 and Q4
signs alone establish a crossing on one globally regular root branch.
