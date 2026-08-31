# The closed source as two positive currents, one with zero ambient winding

**Result.** At `exp(t)=m`, every prime-power m admits an exact two-current
representation of the existing closed source. One current has arbitrary
winding; the other has zero **total** ambient winding. The common site
activity is `y/m^3`. For the honest square tori of this action, the repository's
integer saturation theorem makes the representation valid in every
characteristic, including2. This is not merely a characteristic-zero
nullity argument. A face-height form makes the second current local.

No new source, simulation, enumeration, test or response calculation is used.
The source identity is taken from `bc17b81d:notes/decimation-closed-cluster-gas-action.md`.
The needed integer theorem is §4.1 of
`c1a72e5f0568b1e423ece5d870635030d83f2c7d:notes/digital-alexander-unrestricted-degenerate-quotients.md`.

## 1. Two integer matrices for precisely the existing action

For an occupied set A, let G_A be its occupied **NN graph**, with K vertices,
e edges and c connected components. Retain the torus edge orbits and their
gains. Orient each edge once. Let B_A be its integer incidence matrix and
G_A^wind its two-row matrix of integer deck gains in a period basis P.
Specifically, for an edge u→v of lifted unit displacement d_e choose vertex
representatives x_v and set

\[
P g_e=x_u+d_e-x_v,\qquad g_e\in\mathbb Z^2.
\]

These are **deck coordinates**; the Cartesian unit-displacement rows are
not their substitute modulo a prime. Changing vertex representatives adds
a coboundary to the gain rows and does not change their action on cycles.
Put M_A=[B_A;G_A^wind]. Then, over Q,

\[
\beta_1=\operatorname{nullity}B_A=e-K+c,\qquad
\beta_{\rm null}=\operatorname{nullity}M_A=\beta_1-r,
\]

where r is the ambient image rank. Thus the given identity is exactly

\[
S_*(A)=\beta_1(A)+\beta_{\rm null}(A)-3K+2N+1.
\tag{1}
\]

The zero-ambient kernel allows opposite windings on different components
to cancel. It does not require every loop, or every component current, to
have zero winding separately. These stronger restrictions would change (1).

## 2. Why reduction modulo any prime preserves these nullities

Incidence matrices have no rank-changing primes: a spanning-forest minor
has determinant±1. More explicitly, an integer unimodular edge change of
basis splits edges into forest coordinates and a basis C_A of integral
graph cycles. Integer row operations transform M_A to a block with

\[
I_{K-c}\quad\hbox{and}\quad H_A=G_A^{\rm wind}C_A.
\]

The image of H_A is the integral ambient homology image. It is saturated
in Z² by the cited carrier theorem. For completeness, on the present honest
tori a regular neighborhood of each NN component is an embedded subsurface.
A genus-one component has an intersection-one pair and hence full image Z²;
a genus-zero component has zero image or one primitive boundary direction.
Disjoint essential components cannot have different slopes, since their
nonzero intersection would contradict disjointness. Hence the image of the
**whole** occupied graph is also zero, a primitive line, or Z².

It follows that H_A has only unit nonzero Smith factors. So does M_A after
the forest split. Therefore, for every prime ℓ and every subset A,

\[
\operatorname{rank}_{\mathbb F_\ell}B_A=K-c,\qquad
\operatorname{rank}_{\mathbb F_\ell}M_A=K-c+r.
\tag{2}
\]

There is no residual bad-prime set for this specified graph model. This
also covers m=2 without an odd-characteristic qualification. The unrestricted
repository theorem retains lifted loops and parallel edges; a simplified
quotient which discards them is not licensed by that theorem. The action
identity itself remains under its original honest-cell counting convention.

The distinction is substantive. An abstract one-loop gain graph with gain
(2,0) has a zero-dimensional ambient-zero kernel over Q and a one-dimensional
one over F2. Likewise, on the axis period P=5I, reducing Cartesian winding
displacements modulo5 would erase every period. These are failures of the
assumptions/convention, not exceptions to (2). Generic integer gain matrices
would instead require avoiding primes dividing their nonunit Smith factors.

## 3. Exact positive two-current partition function

For `m=ℓ^a` let

\[
\mathcal C_A=\ker(B_A:\mathbb F_m^{e}\to\mathbb F_m^{K}),\qquad
\mathcal C_A^0=\ker(M_A:\mathbb F_m^{e}\to\mathbb F_m^{K+2}).
\]

Equation (2) gives `|C_A|=m^beta1` and `|C_A^0|=m^beta_null`.
Consequently the unnormalized occupation-odds partition function obeys

\[
\boxed{\sum_A y^{|A|}m^{S_*(A)}
 =m^{2N+1}\sum_A (y/m^3)^{|A|}
       \sum_{x\in\mathcal C_A}\sum_{z\in\mathcal C_A^0}1.}
\tag{3}
\]

Conditional on A the two currents are independent uniform choices from
their respective spaces. Summing over the shared occupation set couples
them; unconditional independence is not asserted. Every local weight is
nonnegative for y>0. For the Bernoulli convention `y=p/(1-p)`, multiply both
sides of (3) by `(1-p)^N` before normalizing.

This is not two unrestricted circulations: their product would overcount
by `m^r`. The zero-ambient projection is exactly the rank correction in S*.
At fixed K, adding one independent graph cycle multiplies the cycle weight
by `m^(2−Δr)`: m² if ambient rank stays fixed, m if it increases by one.
The former includes differences of parallel essential cycles, not only
individually contractible loops.

The unit Smith factors give slightly more than the field statement:
for **any integer m≥2**, the same cardinalities hold for currents in the
cyclic coefficient group Z/mZ. This is a group-current representation,
not a nonexistent finite field of composite non-prime-power order.
At m=1 use the trivial group. For noninteger exp(t), the already-defined
positive cycle fugacity remains meaningful, but no finite-field counting
interpretation is claimed.

## 4. Replace the zero-winding constraint by local face heights

Use the full square torus cellulation, not only fully occupied faces.
Let h assign an F_m value to every face, with one reference face fixed to0.
The boundary current z=∂₂h has zero divergence and zero ambient winding.
Conversely every such current on the torus is a face boundary, uniquely
after this constant-height gauge fixing. Requiring support on occupied
edges means simply

\[
h_{f_L(e)}=h_{f_R(e)}\quad\hbox{for each absent occupied-NN edge }e.
\]

Thus (3) has the entirely local form

\[
m^{2N+1}\!\sum_{n_v\in\{0,1\}}(y/m^3)^{\sum_v n_v}
 \sum_{x_e\in\mathbb F_m}\sum_{h_f\in\mathbb F_m:\,h_{f_0}=0}
 \prod_v\mathbf1\{(Bx)_v=0\}
 \prod_{e=uv:\,n_un_v=0}
     \mathbf1\{x_e=0\}\mathbf1\{h_{f_L(e)}=h_{f_R(e)}\}.
\tag{4}
\]

One flow is divergence-free with unrestricted topology; the other is the
boundary of a face-height configuration. All edge support and vertex
constraints are local. Without fixing h at one face, divide the height
sum by m, changing the displayed external prefactor to m^(2N).
The same construction works with Z/mZ coefficients by the cellular
homology of the torus. A zero-homology current need not bound a region
inside the occupied graph; heights live on the full ambient faces.

## Interpretation and boundary

Equations (3)–(4) supply a named positive **two-current/one-height cycle gas**
for the unique previously selected closed source. They explain both its
cycle reward and its ambient-rank correction without adding tunable source
coefficients. They do not claim a new generic-Q lift, an unproved rank
stability, a continuum field, or a prediction for global U. No new numerical
evidence or sampling block is introduced.
