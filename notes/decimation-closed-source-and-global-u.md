# The decimation source closes after its forced local completion

## A completed global-U mechanism result

**Subsequent completed decisions:** the positive count version of this source
is `Sstar=Ŝ=C+F4+Bvac`. Execution's [full-source result `ec01768f`](https://github.com/LightChainr/Matching-One/blob/ec01768f520e85f1acfd9d3fde9bcf855477254e/results/p337-closed-source-n25/REPORT.md)
gives V_Sstar=+.126165363414 at N25. [The exact cycle/rank subtraction](../results/decimation-cycle-rank/REPORT.md)
leaves 2V_beta_null=+.072917828300 after removing the fixed q term.
The [one-hole result `f5c4a74a`](https://github.com/LightChainr/Matching-One/blob/f5c4a74a20bad8589c39e1034cfb209462110dbe/results/p337-endpoint-defect/score/REPORT.md)
now gives Xi=−10.755718407564 and R=+27.766563581230 at N50, excluding
thermal-only and source-independent gain extensions. The endpoint identities
proved below remain exact. The original one-hole question at the end is retained
as its derivation; its first calculation is no longer an open assignment.

The [exhaustive N25 calculation](../results/decimation-plaquette-u/score/REPORT.md)
has already measured the missing term:

\[
V^{F4}_{25}=+0.19441468646090693,\qquad
2^{13/8}V^{F4}_{25}=+0.5996568681566026.
\]

The exact rational enclosure of the reduced response excludes zero. Thus
bare cluster fugacity fails the specified endpoint **original-U** transport,
not just configuration-level source closure. This is a completed finite-pair
mechanism discriminator, with no Monte Carlo uncertainty or fitted amplitude.
The source is bulk `exp(tS)`, and U includes its moving matching root and
thermal-slope denominator. The frozen contract is
[`2bfe9b904bec0a21b244942ffb7c553813e4e740:analysis/decimation_plaquette_u_contract.json`](https://github.com/LightChainr/Matching-One/blob/2bfe9b904bec0a21b244942ffb7c553813e4e740/analysis/decimation_plaquette_u_contract.json).

There is nevertheless a precise closed completion. In the fixed basis
`(Ctot,F4,T_NN,K,constant)`, define

\[
S_N^*=C_N+F_N+T_N-4K_N,\qquad \widehat S_N=S_N^*+2N.
\]

At the checkerboard saturation/complement/rotation endpoint,
`S_parent*=S_child*−N_parent`, and **Ŝ_parent=Ŝ_child exactly**.
Two legal nested endpoint reductions send bare parent C to grandchild Ŝ.
These are finite endpoint identities, not an interior RG or critical-field
claim. The machine-readable counterpart is
[decimation_source_dictionary.json](../analysis/decimation_source_dictionary.json).

## Incidence counting fixes the entire five-coordinate dictionary

Let the parent area be `N=2M`, its even sublattice A fully occupied, and O
the occupied subset of its odd sublattice B. Require both parent periods
to preserve parity: `Λ⊂(1+i)Z[i]`. Map B by `f(w)=(w−1)/(1+i)` and let
the child occupied set be the complement `U=(Z[i]/Λ')\f(O)`, where
`Λ'=Λ/(1+i)` has M vertices.

On either lattice, C counts occupied-NN plus vacant-matching connected
components, T counts occupied NN edge orbits once each, K counts occupied
vertices, and F counts unit-face orbits with all four corners occupied.
T and F retain the lifted incidence multiplicities. A child loop contributes
twice to endpoint incidence, even though its edge orbit is counted once.

Each parent unit face has two A corners and two B diagonal corners. The B
pair maps to a child NN edge. Parent faces based at even sites give one
child edge direction; faces based at odd sites give the other. The map is
bijective on face/edge orbits: there are 2M of each. With child indicators η,

\[
F_p=\sum_{\langle ij\rangle_c}(1-\eta_i)(1-\eta_j)
    =2M-4K_c+T_c.
\]

Every occupied B site contributes four occupied parent edges, so
`T_p=4(M−K_c)`, while `K_p=2M−K_c`. The already established component
bijection gives `C_p=C_c+F_c`: each extra parent component is an isolated
filled A site and corresponds to one fully occupied child face.
That component/homology argument is pinned at
[`56838d5f068f6f0ba7795926dc9343229bdd28ce:notes/square-checkerboard-endpoint-homology.md`](https://github.com/LightChainr/Matching-One/blob/56838d5f068f6f0ba7795926dc9343229bdd28ce/notes/square-checkerboard-endpoint-homology.md).

Thus, with rows interpreted as images of the parent observables,

\[
\begin{pmatrix}C_p\\F_p\\T_p\\K_p\\1\end{pmatrix}
=
\begin{pmatrix}
1&1&0&0&0\\
0&0&1&-4&2M\\
0&0&0&-4&4M\\
0&0&0&-1&2M\\
0&0&0&0&1
\end{pmatrix}
\begin{pmatrix}C_c\\F_c\\T_c\\K_c\\1\end{pmatrix}.
\]

For a parent weight `exp(cC+fF+τT+kK+a)`, the child coupling vector is

\[
(c,f,\tau,k)\longmapsto(c,c,f,-4f-4\tau-k),\qquad
a' = a+M(2f+4\tau+2k).
\]

This is an exact identity of weights, including the additive constant;
no source coefficient is estimated from data.

## The fixed completion and its two-level meaning

Substitution yields

\[
S_{2M,p}^*=S_{M,c}^*-2M,\qquad
\widehat S_{2M,p}=\widehat S_{M,c}.
\]

The coupling vector `(1,1,1,−4)` is fixed; the multiplicative partition
factor `exp(−2Mt)` cancels from normalized expectations. Up to scale and
constants, this is the fixed source in this declared four-source span.
It is not uniqueness among all possible microscopic sources.

After quotienting constants and a common K source, the endpoint operator is

\[
\mathcal D:C\mapsto C+F,\quad F\mapsto T,\quad T\mapsto0,
\qquad
\mathcal D^2(aC+bF+cT)=a(C+F+T).
\]

The two transient directions disappear after two applications. This source
algebra is not a time-evolution operator or an LCFT Jordan classification.
Restoring the thermal and constant terms, a genuine two-level configuration
map with grandchild size L obeys

\[
C_{4L,p}=C_{L,g}+F_{L,g}+T_{L,g}-4K_{L,g}+2L
         =\widehat S_{L,g}.
\]

Two reductions require **both** parent and child periods to preserve parity:
`Λ⊂(1+i)^2 Z[i]=2Z[i]`. Area divisibility alone is insufficient. Further
levels require `Λ⊂(1+i)^r Z[i]` and the corresponding nested saturation
constraints. They are not successive applications of a single interior s.

The actual N25 pair has generators `(5,0)` and `(4,3)` and does not admit
checkerboard coloring. Consequently the measured N50→N25 reduction cannot
be iterated starting at N25. A legal two-level extension instead has direct
Gaussian generators

```text
N100: (0,10), (-6,8)
  -> N50: (5,5), (1,7)
  -> N25: (5,0), (4,3).
```

This is an exact geometry/source statement; no N50 or N100 enumeration is
claimed or required for the dictionary.

## The same dictionary transports the original root-normalized U

Use `p_A=s+(1−s)p, p_B=p`, with thermal differentiation holding s and all
source couplings fixed. At s=1 only B changes with p. Let
`Q=mean_g E[q_g]`, `Y=P4 E[E_g]`, `q=rank−1`, `E=q²`, and

\[
U_N=\frac{N^{13/8}}2\frac{Y_p}{Q_p}\bigg|_{Q=0}.
\]

The child Bernoulli probability is `1−p`. Complement reverses the E thermal
derivative; rotation by −π/4 reverses the fixed Δcos4 projector. These signs
cancel. With couplings g mapped to g′ by the matrix above, normalized finite
profiles therefore give

\[
U_{2M}^{end}(g)=2^{13/8}U_M(g').
\]

This statement uses corresponding regular root branches. Positivity of a
finite exponential weight does not by itself assert global root uniqueness
for every coupled model; the homogeneous root has a regular continuation
for sufficiently small couplings. The topology/thermal convention is pinned at
[`207436518db46dd13ef0ec91168cb1c99d52eaea:notes/p337-checkerboard-decimation-global-u.md`](https://github.com/LightChainr/Matching-One/blob/207436518db46dd13ef0ec91168cb1c99d52eaea/notes/p337-checkerboard-decimation-global-u.md).

A constant source cancels on normalization. A common K tilt changes only
the Bernoulli odds; its root-comoving U derivative is zero, because the
thermal Jacobian multiplies numerator and denominator equally. The K
coefficient must be common across the paired geometries. Thus at zero
couplings, in common bulk source units,

\[
\begin{aligned}
V_C^{end}&=2^{13/8}(V_C+V_F),\\
V_F^{end}&=2^{13/8}V_T,\\
V_T^{end}&=V_K^{end}=0,\\
V_{S^*}^{end}&=2^{13/8}V_{S^*}.
\end{aligned}
\]

Each V includes its own source-induced root and slope motion. In particular
the measured positive F4 term is exactly the missing bare-cluster endpoint
correction. The fixed family also has the finite-coupling identity
`U_(2M)^end,S*(t)=2^(13/8) U_M^S*(t)` on corresponding regular roots;
Ŝ and S* define identical normalized measures. At two nested endpoints,
bare C transports to this completed family with the factor `4^(13/8)`.
For this last statement hold both saturation constraints fixed and vary
only the remaining Bernoulli probability. It is not the thermal derivative
of an unconstrained one-level family with all its B sites still variable.

## What the result does and does not settle

The finite N25 counterexample uses the square quotients `Z5×Z5` and `Z25`,
which have different Smith classes. Its signed direction contrast and exact
`Δcos4=1152/625` are fixed by the contract; it is not an asymptotic H4
measurement or a copy of the N65/N85 production family. Its interval is
rational numerical enclosure conditional on exhaustive graph counts, not
a sampling confidence interval. The endpoint source correction and its
closed completion are now results, not requests for another first test.

For aliased quotients retain lattice edge/face incidences and winding gains;
collapsing parallel edges or discarding nonzero-gain loops can invalidate
the topology or the local count dictionary. The actual N25 pair has honest
unit cells. These statements do not supply an interior decimation semigroup,
a critical fixed point, a continuum field identity or a revival of the
completed P154/P334 source decisions.

## The completed one-hole question: does a saturation defect leave the fixed family?

The endpoint-fixed source suggests a specific interior question, without
choosing new features or launching production. Put `ε=1−s`, and use the
same fixed Ŝ coupling t. For each saturated A site a let `X_+` be the
endpoint configuration and `X_-a` the same configuration with a removed;
write `Δ_a Ŝ=Ŝ(X_-a)−Ŝ(X_+)`. At fixed p the exact one-hole insertion is

\[
J_{O,\epsilon}(p,t)
=(1-p)\sum_{a\in A}
\mathbb E_{sat,t}\!\left[
e^{t\Delta_a\widehat S}
\{O(X_{-a})-\mathbb E_{sat,t}[O(X_+)]\}\right].
\]

It follows by differentiating the normalized finite measure: the loss of
the all-A configuration cancels between numerator and partition function.
The four-terminal connectivity around the removed star is retained, so
this is a concrete defect operator rather than an assumed changed child
occupation probability.

A stronger **interior thermal-only extension** would require these q/E
insertions, in both geometries, to equal one common `b(p,t) ∂p E[O]` near
the mapped root. It predicts the specific mixed response

\[
\Xi_N=\left.\partial_t\partial_\epsilon
U_N(s=1-\epsilon,t)\right|_{\epsilon=t=0}=0.
\]

Whether the exact one-hole operator instead gives a nonzero Ξ after all
root/slope terms was the specified question. It is now answered by the pinned
result above; endpoint closure alone did not decide it. Nonzero Ξ locates
source-dependent interior transmission in this finite chart, not a full
intermediate curve or a continuum operator.

## A fixed next discriminator: weighted rank jumps versus baseline reweighting

The [exact defect operator, `bc17b81d`](https://github.com/LightChainr/Matching-One/blob/bc17b81d502fb1ca3323f5c20f63c544bb31602d/notes/checkerboard-single-defect-source.md)
gives, for each geometry and one translated A site,

`j_epsilon,O=M(1−p)[E(w*Delta O)+Cov(w,O_plus)]`,

with `w=exp(t*Delta S)`, `Delta S=S_minus−S_plus=3−2k_null−ell`.
Only an alternating child face can lose ambient rank, and ell≤1. Both
rank-changing and rank-preserving holes may contribute to the covariance
term; it is not a population share assigned only to the latter.

The **weighted-rank-jump-only response model** drops that covariance term.
It predicts zero mixed original-U contribution from baseline reweighting.
Since the covariance term is identically zero for every p at t=0, its
t derivative is simply

`h_O=M(1−p) Cov_0(Delta S,O_plus)`.

On the same baseline root put hQ=mean(hq), hY=P4(hE), D=Q_p and r=Y_p/D.
The source derivative of the projection operator/root multiplies an
identically zero function and drops out. The required mixed response is

`Xi_reweight=A/D * partial_p(hY−r*hQ)`

`=A[hY_p/D−Y_pp*hQ/D^2−Y_p*hQ_p/D^2+Y_p*Q_pp*hQ/D^3]`.

The four terms retain the mixed root displacement `−hQ/D`, the slope
response and the derivative of the dose M(1−p). Only two additional
per-K_B cross sums are needed beyond the existing exact packets:

```text
sum(S_minus*q_plus), sum(S_minus*E_plus).
Cov(Delta S,O_plus)
 = E(S_minus*O_plus)−E(S_plus*O_plus)
   −[E(S_minus)−E(S_plus)] E(O_plus).
```

A nonzero rational Xi_reweight enclosure would exclude this weighted-jump-only
model without fitting a share or changing the source. Zero would not prove
full profile closure: projection cancellations remain possible. This is a
fixed two-term identity with one missing cross-moment packet, not a request
for a generic descriptor, source or production search. Keep the original
Sstar and saturation chart; removing a K term without transforming that
two-parameter chart can change its mixed derivative.
