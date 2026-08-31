# Crossing recoupling activates every two-standard-leg colour block

**Result.** In the fixed `NE | SW` cut, the C4-averaged local tensor
`Kbar` has a nonzero singlet block, a negative standard block, and
nonzero `[Q-2,2]` and `[Q-2,1,1]` blocks. It still annihilates an
individually constant colour leg. These statements are compatible:
two standard legs can couple to a singlet. The bare all-equal overlap
zero is a cancellation between its singlet and standard contributions,
not a colour selection rule removing every thermal coupling.

The same block calculation gives a pole in the simplest two-insertion
closure: `Tr(Kbar^2)` has residue1/2 at Q1. Single-insertion removable
closure therefore does not establish an arbitrary-strength local Q1
family. It does not invalidate the already computed first response.

The starting kernel is the
[rational four-port projector](https://github.com/LightChainr/Matching-One/blob/9dc3c4269c1c44a0a0a82b15f8cd13f922a8b2d4/notes/closed-source-local-four-port-pair-kernel.md).
The [local original-U result](https://github.com/LightChainr/Matching-One/blob/923f66b979a6b6132875f783106c041ed3c0c1a9/notes/local-four-port-transmission-result.md)
is an input: `V_local=+0.0018155512845251097`. Neither that result nor
the distinct stable seam response is recalculated here. This note uses
symbolic finite-colour algebra; it contains no new occupation enumeration
or numerical mechanism score.

## 1. The fixed-cut space and the surviving one-leg zero

Initially take integer Q>=4. Let V be the Q-dimensional point-colour
permutation space, with normalized constant vector s, and put

```text
H = s^perp,       dim H=Q-1,
E = I-|s><s|.
```

The ordered-pair cut space is `V tensor V`. Denote the unrotated
embedded pair projector by P and its quarter-turn reshuffling by R:

```text
P_(ab;cd)=Pi2(a,b;c,d),
R_(ab;cd)=Pi2(b,c;d,a),
Kbar=(P+R)/2.
```

The rotation is a reshuffling of four indices, not conjugation by an
operator on this fixed pair cut. In particular it need not preserve
which irrep carries the original projector.

The incidence vectors `v_b=sum_(a!=b)|{a,b}>` are annihilated by the
unordered pair projector P2. Consequently summing Pi2 over **any one**
of its four colour indices gives zero, with the other three fixed.
Reshuffling and averaging preserve these four identities. Therefore

\[
\boxed{\overline K=(E\otimes E)\overline K(E\otimes E).}
\tag{1}
\]

The entire `1 tensor 1`, `1 tensor H` and `H tensor 1` subspace is
killed. This is stronger than its fully contracted all-free zero, and
shows that the bare constant matrix is still annihilated on both sides.
It does not remove the singlet contained in `H tensor H`.

Global colour relabelling commutes with Kbar. On its support the
multiplicity-free decomposition is

```text
H tensor H = [Q] + [Q-1,1] + [Q-2,2] + [Q-2,1,1].
```

The first three summands lie in the symmetric pair space; the last
lies in the antisymmetric pair space. Kbar commutes with simultaneous
input/output pair swap. No matrix element between inequivalent S_Q
irreps is being asserted: crossing recoupling produces a sum of fixed-cut
blocks, not a failure of colour invariance.

## 2. Complete fixed-cut spectral resolution

Write

```text
D=(Q-1)(Q-2),             k=Q(Q-3)/(4D).
```

Let `P_lambda^HH` be the orthogonal projector onto the indicated
summand of `H tensor H`. Then

\[
\boxed{\overline K=
 \frac{Q(Q-3)}{4(Q-1)}P_{[Q]}^{HH}
 -kP_{[Q-1,1]}^{HH}
 +\frac{3Q^2-9Q+8}{4D}P_{[Q-2,2]}^{HH}
 +kP_{[Q-2,1,1]}^{HH}.}
\tag{2}
\]

In the full ordered-pair space this means:

| Irrep | Total multiplicity in `V tensor V` | Kbar eigenvalues on its multiplicity space |
|---|---:|---|
| `[Q]` | 2 | `0`, `Q(Q-3)/(4(Q-1))` |
| `[Q-1,1]` | 3 | `0`, `0`, `-k` |
| `[Q-2,2]` | 1 | `(3Q^2-9Q+8)/(4D)` |
| `[Q-2,1,1]` | 1 | `+k` |

All four coefficients in (2) are nonzero at integer Q>=4. Hence
`rank Kbar=(Q-1)^2`, whereas `rank P=Q(Q-3)/2`. The averaged tensor is
neither a pure fixed-cut pair projector nor an orthogonal projector.
Its standard block is strictly negative even though its added vacant
vertex tensor can retain nonnegative entries at sufficiently small
perturbation strength.

### Singlet and standard multiplicity matrices

The two normalized singlet vectors of the full pair space can be chosen
as

```text
d = sum_a |aa>/sqrt(Q),
o = sum_(a!=b) |ab>/sqrt(Q(Q-1)).
```

In this basis the singlet multiplicity matrix is

\[
\overline K_{\rm singlet}
 =\frac{Q-3}{4}
 \begin{pmatrix}
 1&-1/\sqrt{Q-1}\\
 -1/\sqrt{Q-1}&1/(Q-1)
 \end{pmatrix}.
\tag{3}
\]

The zero eigenvector is the individually constant state
`s tensor s=(d+sqrt(Q-1)o)/sqrt(Q)`. The nonzero singlet is

```text
t=(sqrt(Q-1)d-o)/sqrt(Q)
 = sum_(a,b)(delta_ab-1/Q)|ab>/sqrt(Q-1).
```

It is a correlated two-standard-leg state, not an individually
constant colour leg.

For a real vector z with `sum z_a=0`, `sum z_a^2=1`, choose the two
symmetric standard copies

```text
d_z=sum_a z_a |aa>,
o_z=sum_(a<b)(z_a+z_b)(|ab>+|ba>)/sqrt(2(Q-2)).
```

Their matrix is

\[
\overline K_{\rm standard,sym}
 =\frac{Q-3}{4(Q-1)}
 \begin{pmatrix}
 -1&\sqrt{2/(Q-2)}\\
 \sqrt{2/(Q-2)}&-2/(Q-2)
 \end{pmatrix}.
\tag{4}
\]

It has eigenvalues0 and `-k`. The third standard copy, in
`1 tensor H-H tensor 1`, is killed by (1). Thus the diagonal/unequal-pair
basis exhibits genuine multiplicity mixing, while the invariant
decomposition `1+H` makes the active blocks transparent.

### The two higher pair blocks

For completeness, use the normalized symmetric unequal-pair basis.
The R entries depend only on the intersection size of the two pairs:

| Pair intersection | R entry in this symmetric pair subspace |
|---:|---:|
| 2 | `(Q-3)/(2(Q-1))` |
| 1 | `-(Q-3)/(2D)` |
| 0 | `2/D` |

On the pair-incidence kernel `[Q-2,2]`, the intersection-one adjacency
matrix has eigenvalue-2 and the disjoint-pair adjacency matrix has
eigenvalue1. Therefore

```text
R_[Q-2,2]=(Q^2-3Q+4)/(2D),
Kbar_[Q-2,2]=(1+R_[Q-2,2])/2.
```

In the antisymmetric unequal-pair basis, let B be the oriented
vertex-edge incidence matrix of the complete colour graph. Put
`alpha=(Q-3)/(2(Q-1))`, `beta=(Q-3)/(2D)`. The R matrix is

```text
R_antisym=(alpha+2 beta)I-beta B^T B.
```

The incidence image has eigenvalue Q for `B^T B` and is annihilated.
Its kernel is `wedge^2 H=[Q-2,1,1]`; R has eigenvalue `2k` there.
Together with (3)-(4), this proves the complete resolution (2).

## 3. Bare thermal overlap zero is a cancellation, not a scaling selection rule

In the fixed cut the two original vertex matrices are

```text
V0_(ab;cd)=1,
V1=sum_a |aa><aa|.
```

Equation (1) implies `Kbar V0=V0 Kbar=0` exactly. For V1, however,
(3)-(4) give

```text
<d,Kbar d>=(Q-3)/4,
<d_z,Kbar d_z>=-(Q-3)/(4(Q-1)).
```

Consequently its Frobenius overlap vanishes by cancellation:

\[
\operatorname{Tr}(\overline K V_1)
 =\frac{Q-3}{4}
 +(Q-1)\left[-\frac{Q-3}{4(Q-1)}\right]=0.
\tag{5}
\]

Neither of the two terms is zero at Q>=4. In particular Kbar has a
nonzero singlet block and does not annihilate V1 as an operator.

A small exact counterexample shows why the zero is not preserved in
an arbitrary colour-invariant propagation context. Set

```text
T_eta=I+eta |d><d|,          eta>0.
```

This is S_Q-invariant, positive definite and entrywise nonnegative.
Its extra term is an equality tensor `delta_ab delta_cd/Q`. Yet

\[
\boxed{\operatorname{Tr}(\overline K T_\eta V_1)
 =\eta\frac{Q-3}{4}\ne0.}
\tag{6}
\]

Thus even positivity and colour symmetry do not turn (5) into a
context-independent thermal orthogonality theorem. Equation (6) is
a transfer-context counterexample to that algebraic implication; it
does not assert that T_eta is the actual long-distance transfer matrix
of the original lattice family.

What (5) does prove is the previously established microscopic statement:
at fixed integer Q, Kbar is outside `span{V0,V1}` in the ordinary
positive Frobenius metric. This gives an additional bare interaction
coordinate, not an RG eigenoperator. A microscopic coordinate outside
the bare activity line can still overlap a thermal scaling field after
propagation and renormalization. Conversely the nonzero singlet block
does not prove that the actual thermal scaling coefficient is nonzero;
that coefficient remains a dynamical question.

The exact one-leg constant zero in (1) and the original regular endpoint
zero remain intact. They do not kill the correlated singlet t. Likewise
C4 invariance does not select one continuum spin-four state: it permits
the different spatial harmonics invariant under a quarter-turn. A
fixed-cut colour block and a continuum scaling-field assignment must
not be substituted for each other.

## 4. The two-insertion Gram closure and the Q1 limitation

For integer Q>=4, P is an orthogonal projector of trace
`d2=Q(Q-3)/2`. Reshuffling preserves its Frobenius norm, so

```text
Tr(P^2)=Tr(R^2)=d2,
Tr(PR)=d2 R_[Q-2,2].
```

The cross contraction and the averaged square are therefore

\[
\boxed{\langle P,R\rangle
 =\frac{Q(Q-3)(Q^2-3Q+4)}{4(Q-1)(Q-2)},}
\tag{7}
\]

\[
\boxed{G(Q)=\operatorname{Tr}(\overline K^2)
 =\frac{Q(Q-3)(3Q^2-9Q+8)}{8(Q-1)(Q-2)}.}
\tag{8}
\]

Equivalently `G=d2 Kbar_[Q-2,2]`. Writing `s=Q-1`, the Laurent
expansions begin

```text
<P,R> = 1/s+1+O(s),
G(Q)  = 1/(2s)+O(s).
```

In particular the Q1 Gram residue is1/2. This is the rational
continuation of a fully closed two-insertion tensor network, not
substitution of one colour into a singular projector entry. The
single-insertion Bell4 contraction cancellation does not extend
automatically to this two-insertion closure.

All statements about literal colour irreps and positive Frobenius norms
above initially concern integer Q>=4. The separate block projectors,
basis normalizations and eigenvalues cannot be specialized individually
at Q1 as ordinary Hilbert-space objects. Their rational closed
contractions are the appropriate objects to continue; (8) shows one
such contraction with a genuine pole.

The physical force of this obstruction must be matched to its scope.
A realizable exterior that wires two marked vacant vertices into this
Gram closure supplies a conditional two-site obstruction. With
independent site activities and marked couplings it also isolates a
specific multivariate coefficient. A pole in one such coefficient does
not, without a further argument, prove that the homogeneous,
unconditioned partition has a pole: summing configurations or summing
different marked pairs can produce cancellations. Nor does it undo the
already finite and nonzero first derivative of original U.

## 5. What mechanism space is actually reduced

- **Excluded identification:** C4 Kbar as a pure `[Q-2,2]` fixed-cut
  propagation projector. Its support is all `H tensor H`, with the four
  coefficients in (2).
- **Excluded inference:** bare constant/all-equal Frobenius orthogonality
  implies no thermal coupling in every physical propagation context.
  The singlet/standard cancellation and (6) identify the failure of
  that inference explicitly.
- **Not established by the linear result:** a pole-free finite-strength
  Q1 local interaction with arbitrary repeated insertions. Equation (8)
  requires that nonlinear completion to be specified and examined in
  physical exterior contexts.
- **Still established:** the prescribed single local tensor insertion
  has a finite closed Q1 definition and the completed nonzero original-U
  response. It is distinct from the global seam contraction, and no
  new N25 source fit is needed to preserve that result.

The useful next prediction must therefore concern the fixed local
interaction's long-distance response or a specified nonlinear
completion. Assigning its response to a single colour block, calling
its bare orthogonality a scaling selection rule, or automatically
exponentiating the formal Q1 vertex would each change the scientific
claim beyond what has been proved.
