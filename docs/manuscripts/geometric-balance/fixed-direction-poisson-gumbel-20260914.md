# Fixed primitive directions: no-prefactor Poisson--Gumbel windows

2026-09-14.  Directional extension of `poisson-birth-windows.md` for a **fixed primitive ambient integer direction**.  The extension is exact at the lattice level through an `SL_2(Z)` coordinate change; it does not rotate the physical interaction or import a directional OZ prefactor.

The conclusion is intentionally narrower than `varying-direction-exponential-centres-20260914.md`: centres are now known for varying directions, but the full component-Poisson/Gumbel machinery is proved here only when the integer direction is fixed, so that the transformed finite-range interaction has fixed range.

## 1. Straightening a fixed integer direction without rotating the model

Fix a primitive ambient integer vector

\[
u=(a,b)\in\mathbb Z^2,
\qquad \gcd(a,b)=1.                                           \tag{1.1}
\]

Choose `v in Z^2` such that

\[
\det(u,v)=1,                                                   \tag{1.2}
\]

and let

\[
A=[u\ v]\in SL_2(\mathbb Z).                                 \tag{1.3}
\]

Use integer coordinates

\[
x=Ay.                                                          \tag{1.4}
\]

Then the original torus with periods

\[
nu,\qquad mv                                                    \tag{1.5}
\]

is graph-isomorphic to a rectangular `n x m` torus in `y` coordinates.

The price is not a rotation of NN edges.  The original NN step set becomes the fixed finite-range set

\[
\mathcal S_{4,u}=\{\pm A^{-1}e_1,\pm A^{-1}e_2\},             \tag{1.6}
\]

and the matching graph becomes the corresponding fixed transform of NN+NNN,

\[
\mathcal S_{8,u}=A^{-1}\mathcal S_8.                          \tag{1.7}
\]

All coordinates remain integer.  Since `u,v` are fixed, the horizontal and vertical jump ranges

\[
R_x(u),R_y(u)<\infty                                           \tag{1.8}
\]

are fixed constants independent of `n,m`.

Thus a fixed oblique direction is exactly an axial problem for one fixed anisotropic finite-range SITE graph.

## 2. Directional mass and centre

Let

\[
\xi_{G,u}(p)=\tau_{G,p}(u)                                    \tag{2.1}
\]

be the inverse-correlation cost per period copy.  On tori with

\[
n\to\infty,
\qquad
\frac{\log m}{n}\to d>0,                                    \tag{2.2}
\]

the directional centre theorem gives unique

\[
\xi_{4,u}(a_u(d))=d,                                         \tag{2.3}
\]

\[
\xi_{8,u}(c_u(d))=d,
\qquad
b_u(d)=1-c_u(d).                                               \tag{2.4}
\]

The two births satisfy

\[
T_1\to a_u(d),
\qquad
T_2\to b_u(d)                                                  \tag{2.5}
\]

in probability.

## 3. Complete winding components on the transformed cylinder

Consider the infinite cylinder

\[
C_n\times\mathbb Z                                             \tag{3.1}
\]

for the transformed graph `G_u`.

A whole empty **slab of thickness `2R_y+1`** disconnects the graph vertically.  Such slabs occur independently at well-separated heights with positive probability.  Hence every occupied component is vertically finite almost surely at each fixed `n` and `p<1`.

For a component with nonzero horizontal homology, anchor it at its lowest `y_2` coordinate, with a deterministic horizontal tie-break.  Define

\[
\nu^{G,u}_n(p)
=E[\text{number of complete winding-component anchors at level }0].\tag{3.2}
\]

This is the fixed-direction analogue of the axial component density.

To define a local truncated anchor with height `H`, use guard **slabs** of thickness at least `R_y` below and above the candidate window.  If the candidate component misses both guard slabs, finite range guarantees that it is the full cylinder component.

Set again

\[
H=n^2.                                                         \tag{3.3}
\]

The local anchor indicator then depends on `O(nH)` SITE variables, with constants depending only on `u` and the chosen graph.

## 4. Uniform localization is unchanged at exponential scale

Fix a compact subcritical parameter interval `I` for the transformed graph.  The graph is exactly isomorphic to the original finite-range graph, so planar one-arm decay is available uniformly on `I`.

Partition a long vertical traversal into slabs of height `c_u n` separated enough to have disjoint SITE supports.  A traversal of one slab forces a planar arm of radius `c'_u n` from one of `O_u(n)` entry sites.  Therefore

\[
P(\text{cross one prescribed }c_un\text{ slab})
\le e^{-c_I n}.                                               \tag{4.1}
\]

Crossing vertical distance `H` requires a linear number `H/n` of disjoint such slab events, hence

\[
\boxed{
P(\text{component height}>H)
\le C_I e^{-c_I H}}                                           \tag{4.2}
\]

for `H>=n`, after changing constants.  With `H=n^2`, localization error is superexponentially small on the `n` scale.

The same block construction gives the cylinder cluster-volume tail used later:

\[
\boxed{
P_p^{C_n\times\mathbb Z}(|C_v|\ge k)
\le C_Ie^{-c_Ik}}                                             \tag{4.3}
\]

uniformly for large `n`, `p in I`.

The only change from the axis proof is a fixed enlargement of blocks/guards by the transformed interaction range.

## 5. Component density has the directional mass exponent

The local upper bound uses the fixed-direction torus/cylinder first-exit estimate: a winding component in an `H`-window contains a nonzero horizontal homology witness and therefore costs

\[
\exp[-(\xi_{G,u}(p)-o(1))n]                                  \tag{5.1}
\]

up to polynomial factors.

For the lower bound, use the fixed finite directional seed from the centre theorem and repeat it around the circumference, closing the seam.  Unless its full component has height greater than `H`, one of at most `H+O_u(1)` anchor levels contains the resulting component.  The long-component error is controlled by (4.2).

Therefore, uniformly for moving `p_n->p` inside a compact subcritical interval,

\[
\boxed{
-\frac1n\log\nu^{G,u}_n(p_n)\to\xi_{G,u}(p).}                \tag{5.2}
\]

No directional OZ amplitude is used.

## 6. Chen--Stein Poisson approximation survives finite-range memory

The truncated anchor window has vertical size `H+O_u(1)` and horizontal circumference `n`.  Two anchors are independent when these windows are vertically separated by more than `2H+O_u(1)`.

The dependency neighbourhood size is therefore

\[
D_n=O_u(nH)=O_u(n^3).                                         \tag{6.1}
\]

For overlapping anchor windows, two distinct complete winding components yield two disjoint occupied winding witnesses.  Enclose each by the increasing event that the enlarged local band contains a horizontal winding.  Site BK gives the same square bound as in the axial proof.

Because every enclosing winding probability is

\[
B_n(p)\le \operatorname{poly}_u(n,H)
 e^{-(\xi_{G,u}(p)-o(1))n},                                  \tag{6.2}
\]

the Arratia--Goldstein--Gordon `b1,b2` terms satisfy

\[
b_1+b_2
\le m\operatorname{poly}_u(n)B_n(p)^2.                       \tag{6.3}
\]

Near the centre, choose a compact interval on which

\[
2\inf_I\xi_{G,u}>d.                                           \tag{6.4}
\]

Then with `log m/n->d`, (6.3) tends to zero exponentially.  Localization is smaller still.

Hence the true complete winding-component count satisfies

\[
\boxed{
Z_{n,m}^{G,u}(p)
\overset{TV}=\operatorname{Poi}(m\nu_n^{G,u}(p))+o(1)}         \tag{6.5}
\]

uniformly on that compact interval.

The process version gives a homogeneous Poisson process of anchor positions in the transformed longitudinal coordinate whenever `m nu_n -> lambda`.

## 7. Same-label lower/upper windows still decouple

Use one common uniform label field.  Near the lower NN birth, black means `U<=p_1`; near the upper birth, white matching means `U>p_2`, with `p_1<p_2`.

The transformed black and white graphs are still the exact matching pair of the transformed cellulation.  The combined dependency graph is local with the same fixed range.

For an overlapping black/white pair, use enclosing winding events.  They are monotone in opposite directions in the common labels, so the same product-measure correlation argument as in the axial proof bounds the joint occurrence by the product of the two marginal enclosing probabilities.

Therefore the two separated component processes converge jointly to independent Poisson processes when their means remain finite.

This is independence of **separated parameter windows**, not same-parameter black/white independence.

## 8. Semiconvex log intensity and directional mass slope

Let

\[
F_n(p)=\frac1n\log\nu^{G,u}_{n,n^2}(p).                       \tag{8.1}
\]

For a complete component `C`, its activity is still

\[
p^{|C|}(1-p)^{|\partial C|}.                                \tag{8.2}
\]

The derivative score is

\[
\frac{|C|}{p}-\frac{|\partial C|}{1-p}.                      \tag{8.3}
\]

The transformed graph has fixed degree, so

\[
|\partial C|\le\Delta_u|C|.                                  \tag{8.4}
\]

The uniform cylinder cluster-volume tail (4.3), together with an explicit fixed ring of activity at least `e^{-K_In}`, gives under the component-Palm law

\[
E_*|C|=O_I(n),
\qquad
E_*|\partial C|=O_I(n).                                      \tag{8.5}
\]

Exactly as in the axial activity calculation, the nonnegative score variance may be discarded to obtain

\[
\boxed{F_n''(p)\ge-C_I}                                      \tag{8.6}
\]

on compact subcritical intervals.

Since by (5.2)

\[
F_n(p)\to-\xi_{G,u}(p),                                      \tag{8.7}
\]

the limit is locally semiconvex, hence `xi_{G,u}` locally semiconcave.  At every differentiability point `a`, for every sequence `p_n->a`,

\[
\boxed{
F_n'(p_n)\to-\xi_{G,u}'(a).}                                 \tag{8.8}
\]

Thus the directional mass slope required for affine parameter scaling is obtained from the component intensity itself; no prefactor expansion is required.

## 9. Median-centred Gumbel theorem in a fixed direction

Call `d` regular for the fixed direction `u` if both

\[
\xi_{4,u}'(a_u(d))
\quad\text{and}\quad
\xi_{8,u}'(c_u(d))                                             \tag{9.1}
\]

exist.  Define positive slopes

\[
v_{4,u}=-\xi_{4,u}'(a_u(d)),
\qquad
v_{8,u}=-\xi_{8,u}'(c_u(d)).                                 \tag{9.2}
\]

Let `a_{n,m}` and `b_{n,m}` be the **true finite medians** of the first and second NN rank births.

Then the same intensity-clock argument as in `poisson-birth-windows.md` gives the joint limit

\[
\boxed{
X_n=v_{4,u}n(T_1-a_{n,m}),
\qquad
Y_n=v_{8,u}n(T_2-b_{n,m}),}                                   \tag{9.3}
\]

with

\[
\boxed{
P(X\le x)=1-2^{-e^x},
\qquad
P(Y\le y)=2^{-e^{-y}},}                                      \tag{9.4}
\]

and `X,Y` independent in the limit.

The same uniform-integrability proof gives the mean, variance, covariance and IQR consequences with `w` replaced by `n` and `v_G` by `v_{G,u}`.

In physical Euclidean width

\[
\ell=n|u|,                                                     \tag{9.5}
\]

we have

\[
v_{G,u}n
=\ell\,[-\partial_p\tau_{G,p}(e)]                            \tag{9.6}
\]

at differentiability points, so the theorem has the expected coordinate-free scaling.

## 10. Exceptional set

Local semiconcavity implies that each one-dimensional function `p -> xi_{G,u}(p)` is differentiable except at at most countably many `p` values.

Since the mass is a strict bijection, the corresponding exceptional `d` set is at most countable.  No SITE `p`-analyticity is asserted.

At an exceptional `d`, the natural intensity clock remains valid; the unique affine `1/n` scaling in `p` is not promoted without a derivative.

## 11. Why varying directions are not automatically covered

For `u_n` changing with `n`, an `SL_2(Z)` straightening produces transformed step sets whose ranges may grow with `n`.  The guard-slab thickness, local dependency range, block-tail constants and component-Palm boundary constant are then no longer uniform for free.

The centre theorem survives because the first-exit Wulff argument is coordinate-free.  The Poisson/Gumbel proof needs stronger local uniformity and is therefore kept separate.

A future varying-direction window theorem should either:

1. prove uniform finite-range controls in a geometrically bounded straightening scheme; or
2. work directly in the physical coordinates with a transverse anchor coordinate and uniform local-component estimates.

Neither step is silently assumed here.

## 12. Claim boundary

This note is an author-level extension of the already author-level no-prefactor Poisson/Gumbel proof.  The scientific novelty of the extension is the exact integer-coordinate reduction to a fixed anisotropic finite-range SITE graph and the observation that every ingredient in the axial proof is stable under that fixed finite-range change.

No OZ prefactor, critical near-window uniformity, or square-site `p`-analyticity is used.