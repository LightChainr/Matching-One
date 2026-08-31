# The closed source has a hypergraph RC family and an exact five-partition q/E readout

**Result.** The local colour gas is exactly a four-body hypergraph
random-cluster model. Summing a specified family of flat colour twists
restores its missing ambient-rank weight. A further full-rank twist then
reconstructs the original `q=r-1` and `E=q²`, with no fitted amplitude or
replacement observer. At `m=exp(t)=2`, **five nonnegative local partition
functions suffice**. The proposed coefficients below are correct.

The new step is the **local hypergraph realization of the already declared
closed-source Q-family**, and its fixed-m partition-function readout.
The annihilator count itself is already proved in `4b179559`:
[finite-abelian twist tomography](https://github.com/LightChainr/Matching-One/blob/4b17955946207f39762f836156aa063f64fbf67d/notes/finite-abelian-twist-tomography.md).
That note reconstructs a fixed underlying rank law from order-2 and
order-3 constraint averages. Here every T, R and I uses **one fixed m,
one local interaction and one activity**; neither the source nor its
colour number is scanned. The positive occupation Q-lift and its weighted
homology root are also already specified by `e4965d47`:
`notes/closed-source-q-lift-and-thermal-quotient.md`.

Thus the gap closed here is constructive: how that particular Q-family
and the original q/E are realized by local hyperedge/twist partitions.
It does not identify standard Potts universality. This note uses algebra
only; no enumeration or new samples enter.

Source convention: `85fd492312b597b3fa102ea913e4bcc7aeae2acf`,
[local colour gas](closed-source-local-colour-gas.md), together with the
[two-current representation](closed-source-two-current-representation.md).
The integral input is §4.1 of
`c1a72e5f0568b1e423ece5d870635030d83f2c7d`:
[digital Alexander saturation](https://github.com/LightChainr/Matching-One/blob/c1a72e5f0568b1e423ece5d870635030d83f2c7d/notes/digital-alexander-unrestricted-degenerate-quotients.md).

## 1. Fixed action, graph and normalization

Use the honest square torus of the source identity, with N vertices,
2N NN edge orbits and N unit faces. For an occupied subset A, write K,
B, C_B and r for its vertex count, occupied NN edge count, occupied
component count and ambient homology image rank. In particular, B is
not the vacant-edge count. The latter is

\[
B_{\rm vac}=2N-4K+B.
\]

The fixed source and the unnormalized occupation-odds partition are

\[
S_*=2B+2C_B-5K-r+2N+1,\qquad
Z_*(y,m)=\sum_A y^K m^{S_*(A)}.
\tag{1}
\]

Set `Q=m²`, `a=y/m⁵`, and initially take integer `m>=2` and `y>0`.
The local gas has a vacant state 0 and Q active colours. Its edge
matrix is `W00=W0c=Wc0=1`, `Wcd=Q δcd` for active c,d. Consequently

\[
Z_{\rm col}=\sum_A a^K Q^{B+C_B},\qquad
Z_*=m^{2N+1}\sum_A a^KQ^{B+C_B}m^{-r}.
\tag{2}
\]

The common Bernoulli factor `(1-p)^N`, with `y=p/(1-p)`, can be restored
at the end. It cancels from all normalized observers in this note.

## 2. Gram factorization gives a genuine hypergraph RC expansion

Index a Q-dimensional edge-colour space by active colours c. The vectors

\[
v_0=Q^{-1/2}(1,\ldots,1),\qquad v_c=Q^{1/2}e_c
\]

have Gram matrix W. Factor each NN edge weight using these vectors, then
sum the state at each original vertex. A vacant vertex contributes
`Q^-2`; an active vertex contributes `a Q²` if its four incident edge
colours agree and zero otherwise. Thus, with `v=aQ⁴=y m³`,

\[
\boxed{Z_{\rm col}
 =Q^{-2N}\sum_{\{c_e\}}
       \prod_{x\in V}\left[1+v\,
         \mathbf1\{c_{e_1(x)}=\cdots=c_{e_4(x)}\}\right].}
\tag{3}
\]

The hypergraph has **2N edge-nodes** and one four-node hyperedge for
each original site. Activating the hyperedges in A imposes those equality
constraints. Every occupied NN component makes one equality component;
an edge with two vacant endpoints is one remaining isolated edge-node.
Therefore its hypergraph cluster count is exactly

\[
c_H(A)=C_B(A)+B_{\rm vac}(A)=C_B+2N-4K+B.
\tag{4}
\]

Expanding (3) proves

\[
Z_{\rm col}=Q^{-2N}\mathcal Z_H(Q,v),\qquad
\mathcal Z_H(Q,v)=\sum_{A\subseteq V}v^{|A|}Q^{c_H(A)}.
\tag{5}
\]

This is a hyperedge random-cluster partition with a single activation
activity v, not an ordinary independent-bond FK expansion. Equations
(3)–(5) explain the factors `Q^-2N` and `v=y m³`; omitting either changes
the occupation law.

The equality-constraint graph also retains precisely the occupied NN
ambient image. To see this, introduce site-centres and edge-midpoints.
Every selected site contributes its four half-edges. A midpoint shared
by two occupied sites subdivides an occupied NN edge; a midpoint with
one occupied endpoint is a dangling leaf. Removing leaves and undoing
subdivision preserves all winding cycles. The isolated vacant-vacant
midpoints add none. This identifies the topology needed in the twists,
not only the cluster count in (4).

## 3. Flat translations select annihilators, not arbitrary rank weights

Write the Q colours as `(σ,τ) in (Z/mZ)²`. Fix an integer period basis
`P=(P1,P2)`. Winding is always expressed in **deck coordinates in P**,
not Cartesian displacement reduced modulo m.

For `α=(α1,α2) in (Z/mZ)²`, impose the flat colour translations

\[
c(x+P_i)=(\sigma(x),\tau(x)+\alpha_i).
\tag{6}
\]

Call the partition (3) with these seam identifications `Z_α`. All local
weights are unchanged and nonnegative. For a fixed A, transporting its
common colour around winding h changes τ by `α·h`. A nonzero translation
has no fixed colour. Hence

\[
Z_\alpha=\sum_A a^KQ^{B+C_B}
  \mathbf1\{\alpha\cdot H_A=0\pmod m\},
\tag{7}
\]

where `H_A` is the integral image of the whole occupied NN graph. A
compatible equality component has exactly Q choices, including those
that wind; an incompatible component has none.

The cited saturation theorem makes `H_A` zero, a primitive line or all
`Z²`. For the whole graph, disjoint essential rank-one components have
the same primitive slope; a rank-two component already supplies `Z²`.
Thus saturation also holds for their union. A unimodular integer basis
change now gives the exact count

\[
\#\{\alpha:\alpha\cdot H_A=0\pmod m\}=m^{2-r(A)}.
\tag{8}
\]

Equation (8) is the previously established finite-abelian annihilator
lemma, used here inside the new local partition realization. It works
for **every integer m>=2**, including composite m:
the coefficient group is cyclic, not an assumed finite field. Without
saturation, (8) would depend on nonunit Smith factors and rank alone
would not suffice.

Let

\[
T=Z_{(0,0)},\quad R=\sum_{\alpha\ne(0,0)}Z_\alpha,\quad D=T+R.
\]

Equations (2), (7) and (8) give the exact normalization

\[
\boxed{Z_*(y,m)=m^{2N-1}D(y,m).}
\tag{9}
\]

In particular, D represents the **projected closed-source law**, while
T alone represents the unprojected local colour law. Their normalized
rank observables must not be interchanged.

## 4. One full-rank twist completes the original q/E reconstruction

Define I using the same colour gas with the two commuting translations

\[
c(x+P_1)=(\sigma(x)+1,\tau(x)),\qquad
c(x+P_2)=(\sigma(x),\tau(x)+1).
\tag{10}
\]

A winding h now translates the colour by `(h1,h2) mod m`. Compatibility
requires `H_A subset m Z²`. A nonzero saturated subgroup cannot satisfy
this, so I selects exactly r=0. No Fourier sign weight is being inserted;
I is another nonnegative local partition function.

For transparent inversion define the untwisted occupation sector sums

\[
L_r=\sum_{A:r(A)=r}a^K Q^{B+C_B}.
\]

Then the complete three-sector linear system is

\[
T=L_0+L_1+L_2,\quad
R=(m^2-1)L_0+(m-1)L_1,\quad I=L_0,
\quad D=m^2L_0+mL_1+L_2.
\tag{11}
\]

The closed-source rank probabilities are

\[
P_0=\frac{m^2I}{D},\quad
P_1=\frac{m[R/(m-1)-(m+1)I]}{D},\quad
P_2=\frac{T-R/(m-1)+mI}{D}.
\tag{12}
\]

Since q=r−1 and E=q², the proposed formulas follow without amendment:

\[
\boxed{\langle q\rangle_*
 =\frac{T-R/(m-1)-m(m-1)I}{T+R},\qquad
\langle E\rangle_*
 =\frac{T-R/(m-1)+m(m+1)I}{T+R}.}
\tag{13}
\]

For m=2 the five partitions are the four `Z_α` and I, and

\[
\boxed{\langle q\rangle_*=(T-R-2I)/(T+R),\qquad
\langle E\rangle_*=(T-R+6I)/(T+R).}
\tag{14}
\]

All five partitions are strictly positive for finite y>0, because the
empty occupation subset is compatible with every twist. Numerator
coefficients may have either sign; their ratios are the original
expectations, not five positive numerator observables. For example,
`P0=4I/D`, `P1=2(R-3I)/D`, `P2=(T-R+2I)/D` are nonnegative by (11).
Five suffice; no claim of minimality among every possible representation
is needed.

The identities hold for all y in each geometry. They therefore determine
the original p-jets, pooled matching root and slope-normalized U wherever
that root and normalization are defined. They do not replace these by a
twisted-sector crossing, or prove the root assumptions afresh.
For a two-geometry pooled root, form each normalized ratio in (13) before
averaging; adding unnormalized numerators is not equivalent. In one
geometry the numerator is precisely `L2-Q L0`, agreeing with the already
declared weighted homology section in `e4965d47`.

## 5. Hypergraph coordinates for the already specified real-Q completion

The polynomial (5) is well-defined for all real Q>0 and v>=0 without
colours. Its **topologically completed** occupation partition is

\[
\Psi(Q,v)=\sum_A v^{K(A)}Q^{c_H(A)-r(A)/2}.
\tag{15}
\]

Using (4), the exact continuation of this particular closed source is

\[
\boxed{Z_*(y,m)=Q^{1/2-N}\Psi(Q,yQ^{3/2}),\qquad Q=m^2>0.}
\tag{16}
\]

This is the hypergraph coordinate form of the previously specified
occupation lift `a^K exp(eta B) Q^(C_B-r/2)` on `eta=log Q`,
`a=y Q^(-5/2)`; it is not a second candidate lift. At fixed y,
the normalized response of a parameter-independent observer O obeys

\[
\frac{d}{dt}\langle O\rangle_*
 =\left(2Q\partial_Q+3v\partial_v\right)
     \langle O\rangle_\Psi
 =\operatorname{Cov}_*(O,S_*),
\qquad Q=e^{2t},\ v=ye^{3t}.
\tag{17}
\]

Indeed the logarithmic subset weight derivative in Ψ is
`2c_H-r+3K=S*+2N-1`; the difference is constant. This expresses the fixed
lift and its thermal-activity drift in the new hyperedge activity
coordinate; a bare Q derivative of the unprojected gas would be a
different source response.

There are three distinct assertions here:

- Hypergraph RC (5) has a positive occupation expansion for arbitrary
  real Q>0. Its colour realization requires integer Q.
- The particular flat-translation construction and (13) as **literal
  colour-twist sums** require Q=m² with integer m>=2. At m=2 it is a
  statement at t=log2, not by itself a derivative at t=0.
- Equation (15) specifies a positive real-Q topological continuation.
  Its rank-sector sums also define q/E for every Q>0. These are not
  finite sets of colour twists at noninteger sqrt(Q).

For example, sector polynomials `H_r(Q,v)=sum_(A:r) v^K Q^c_H` can
continue the combinations in (11) algebraically after multiplying by
`Q^-2N`. That is an explicit occupation-sector prescription, not a
noninteger colour-counting proof. At Q=1 all rank fugacities become1
and (15) reduces to `(1+v)^N`; the three separate H_r still retain the
nontrivial q/E law. Formula (13) is not evaluated by substituting m=1
into its singular `1/(m-1)` factors.

## Scientific consequence

The former local-colour representation carried the correct bulk law
only after an external `m^-r` correction. Equations (3), (9) and (13)
now realize that correction and both original topology observers through
named local positive partitions. Equations (16)–(17) identify the resulting
hypergraph representation with the already supplied real-Q occupation
continuation and source direction. The new result is this constructive
local realization and fixed-m reconstruction, not the old annihilator
count, a second Q-lift discovery, or a new adjustable source.

The hyperedges are four-body stars and the topological completion is
essential. None of these equalities establishes a standard Potts critical
line, a continuum operator identity, universality, or uniqueness among
all other possible generic-Q lifts. Saturation, honest-cell counting and
the original q/E observer convention are part of the theorem's scope.
