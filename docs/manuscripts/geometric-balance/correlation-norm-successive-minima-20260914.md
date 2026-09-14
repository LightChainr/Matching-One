# Correlation-norm successive minima and homology-rank upper bounds

2026-09-14.  General-period consequence of the first-exit torus theta bound.  The correct fixed-subcritical geometry is expressed in the planar correlation norm, not Euclidean period length.

## 1. Correlation-norm lattice geometry

Fix a translation-invariant finite-range site graph `G` at a subcritical parameter `p`.  Let

\[
\tau_p(x)
\]

be its inverse-correlation norm and

\[
B_p=\{x:\tau_p(x)\le1\}.                                     \tag{1.1}
\]

Let `Lambda` be a rank-two period lattice.

Define the first correlation-norm minimum

\[
\boxed{
\rho_1(\Lambda;p)
=\min_{\lambda\in\Lambda\setminus0}\tau_p(\lambda).}          \tag{1.2}
\]

Choose any projective minimizing line `ell_*` represented by a minimizer `lambda_*`.

Define its cheapest nonparallel period cost

\[
\boxed{
\rho_\perp(\ell_*;\Lambda,p)
=\min_{\lambda\in\Lambda:\lambda\notin\operatorname{span}(\ell_*)}
\tau_p(\lambda).}                                             \tag{1.3}
\]

The invariant second successive minimum is

\[
\boxed{
\rho_2(\Lambda;p)
=\min\{R:\operatorname{span}_{\mathbb R}
(\Lambda\cap R B_p)=\mathbb R^2\}.}                           \tag{1.4}
\]

If the shortest projective line is unique, `rho_2=rho_perp(ell_*)`.  If several nonparallel shortest vectors tie, then `rho_2=rho_1`.

## 2. Norm-ball packing controls the period theta sum

Because every nonzero lattice vector has `tau`-norm at least `rho_1`, the translated open norm balls

\[
\lambda+(\rho_1/2)B_p,
\qquad\lambda\in\Lambda,                                     \tag{2.1}
\]

are disjoint.

Let

\[
N_\Lambda(R)=|\{\lambda\in\Lambda:\tau_p(\lambda)\le R\}|.  \tag{2.2}
\]

Every small ball centered at a point counted by `N_Lambda(R)` lies inside

\[
(R+\rho_1/2)B_p.                                               \tag{2.3}
\]

Comparing Euclidean areas and using homothetic scaling of norm balls gives

\[
\boxed{
N_\Lambda(R)
\le\left(1+\frac{2R}{\rho_1}\right)^2.}                      \tag{2.4}
\]

No determinant estimate or Euclidean angle is needed.

For `a>0`, shelling at multiples of `rho_1` yields

\[
\sum_{\lambda\in\Lambda\setminus0}e^{-a\tau_p(\lambda)}
\le
C\sum_{k\ge1}(k+1)^2e^{-ak\rho_1}.                           \tag{2.5}
\]

Hence

\[
\boxed{
\sum_{\lambda\ne0}e^{-a\tau_p(\lambda)}
\le
C_a(\rho_1)e^{-a\rho_1},}                                    \tag{2.6}
\]

where, for example,

\[
C_a(\rho_1)
\le C(1-e^{-a\rho_1})^{-3}.                                  \tag{2.7}
\]

In particular, as `rho_1->infinity`, the theta sum has the same exponential rate as its cheapest period.

## 3. Positive-rank upper bound

`first-exit-torus-winding-upper-20260914.md` proves that for every `epsilon in (0,1)`, on every sufficiently large honest torus,

\[
P_p(r>0)
\le
N C_{p,\epsilon}
\sum_{\lambda\ne0}
 e^{-(1-\epsilon)\tau_p(\lambda)}.                            \tag{3.1}
\]

Combining with (2.6),

\[
\boxed{
P_p(r>0)
\le
N\,\widetilde C_{p,\epsilon}(\rho_1)
 e^{-(1-\epsilon)\rho_1},}                                   \tag{3.2}
\]

with `log Ctilde=o(rho_1)` as `rho_1->infinity`.

Therefore, if along a sequence

\[
\rho_{1,n}\to\infty,
\qquad
\limsup\frac{\log N_n}{\rho_{1,n}}<1,                        \tag{3.3}
\]

then

\[
\boxed{P_p(r>0)\to0.}                                        \tag{3.4}
\]

More quantitatively, if

\[
\frac{\log N_n}{\rho_{1,n}}\to\alpha<1,                      \tag{3.5}
\]

then

\[
\boxed{
\limsup\frac1{\rho_{1,n}}\log P_p(r>0)
\le-(1-\alpha).}                                              \tag{3.6}
\]

The `epsilon` in the finite certificate is sent to zero after the sequence limit.

## 4. Rank-two requires a nonparallel homology class

Fix a projective line `ell_*`.  If the ambient homology image has rank two, then it contains a nonzero class not lying in `ell_*`.

Choose a minimum-length nonzero-homology closed walk whose class is outside `ell_*`.  The same simple-cycle/first-exit argument as for positive rank gives

\[
P_p(r=2)
\le
N C_{p,\epsilon}
\sum_{\lambda\in\Lambda:\lambda\notin\ell_*}
 e^{-(1-\epsilon)\tau_p(\lambda)}.                            \tag{4.1}
\]

The unrestricted packing estimate still bounds the number of nonparallel vectors in each shell.  Therefore

\[
\boxed{
P_p(r=2)
\le
N\,\widetilde C_{p,\epsilon}(\rho_1,\rho_\perp)
 e^{-(1-\epsilon)\rho_\perp},}                               \tag{4.2}
\]

where the prefactor is subexponential in `rho_perp` whenever `rho_1` is bounded below proportionally to `rho_perp`.

If the shortest line is unique and the lattice shape is nondegenerate in the correlation norm, this yields the clean second-minimum criterion

\[
\boxed{
\log N< (1-o(1))\rho_2
\quad\Longrightarrow\quad P_p(r=2)\to0.}                     \tag{4.3}
\]

Even without proportional successive minima, (4.1) is a rigorous restricted theta bound and is the safer statement.

## 5. Unique rank-one slope from a successive-minimum gap

Suppose one projective line `ell_*` has cost `rho_1` while

\[
\rho_\perp-\log N\to+\infty.                                 \tag{5.1}
\]

Then (4.2) gives

\[
P_p(r=2)\to0.                                                  \tag{5.2}
\]

Moreover, a rank-one configuration with slope different from `ell_*` also contains a nonparallel class relative to `ell_*`, so

\[
\boxed{
P_p(r=1,L\ne\ell_*)\to0.}                                    \tag{5.3}
\]

Thus whenever positive rank has nontrivial probability and the second successive direction remains exponentially suppressed, the rank-one projective slope is asymptotically deterministic.

`exponential-homology-class-selection-20260914.md` is the extreme case where `rho_perp` is of order the exponentially large transverse height.

## 6. A correlation-norm homological free-energy upper criterion

The combination appearing in (3.2) is

\[
\boxed{\mathcal F_1=\rho_1-\log N.}                           \tag{6.1}
\]

At fixed subcritical `p`,

- `F1->+infinity` forces rank zero;
- `F1` of order one is the scale where the union bound no longer decides the event;
- `F1->-infinity` by itself does not prove positive rank, because the actual number of approximately independent placements may be smaller than `N`.

So `rho_1-log N` is a rigorous **upper-side homological free energy**.  A matching lower theorem must replace `N` by a geometrically justified opportunity count from disjoint translates.

For an exponentially elongated shortest-period torus that opportunity count is `h/poly(ell)`, whose logarithm agrees with `log N` at the `ell` scale.  This is why the upper and lower exponents close there.

## 7. Relation to the Euclidean full-law criterion

On a fixed compact subcritical parameter interval, finite-range norm equivalence gives

\[
c_p|x|\le\tau_p(x)\le C_p|x|.                               \tag{7.1}
\]

Hence

\[
c_p\ell(\Lambda)
\le\rho_1(\Lambda;p)
\le C_p\ell(\Lambda).                                       \tag{7.2}
\]

The geometric manuscript's condition

\[
\log N/\ell\to0                                              \tag{7.3}
\]

therefore implies

\[
\log N/\rho_1(p)\to0                                        \tag{7.4}
\]

at every fixed subcritical parameter.  The correlation-norm formulation is sharper when period directions/anisotropy matter, while the Euclidean condition remains the clean parameter-uniform geometric statement near criticality.

## 8. Claim boundary

The lattice-point packing and theta estimates are deterministic convex geometry.  The probabilistic inputs are exactly the first-exit torus upper bound and existence of the subcritical correlation norm.  The rank-two restricted-theta bound uses only the deterministic fact that rank two contains a class outside any fixed projective line; it does not require two disjoint essential components.