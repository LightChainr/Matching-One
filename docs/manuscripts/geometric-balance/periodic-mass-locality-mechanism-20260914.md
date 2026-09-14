# Why periodic cylinder mass corrections should be exponential: the zero mode is exact

2026-09-14.  Conditional structural lemma for #760.  It identifies the **only** place a cylinder-versus-plane mass correction can enter in a Markov-additive/OZ description.

The note does not by itself prove that the actual SITE tagged-span mass has the required renewal/decorative representation.  It reduces the desired exponential locality theorem to a concrete exponentially rare wrap-defect estimate.

## 1. Translation-invariant Markov-additive kernel on the plane

Let

\[
A(z,y)=\sum_{x\ge1}\sum_{j\in\mathbb Z}A_{x,j}z^xy^j          \tag{1.1}
\]

be a positive matrix kernel for irreducible long-connection pieces.  The integer `j` is the transverse displacement and `x` the longitudinal advance.  Allow countably many `j`, with an exponential moment

\[
\sum_{x,j}\|A_{x,j}\|R^x e^{c|j|}<\infty                    \tag{1.2}
\]

near the physical root.

Suppose the plane inverse mass is determined by the zero-transverse-character Perron equation

\[
\rho(A(R,1))=1,
\qquad \kappa=\log R.                                         \tag{1.3}
\]

This is the standard form of a translation-invariant Markov-renewal/OZ skeleton.

## 2. Pure periodization does **not** move the zero mode

Put the transverse coordinate modulo `w`.  The periodized displacement matrices are

\[
\bar A_{x,r}^{(w)}
=\sum_{k\in\mathbb Z}A_{x,r+kw},
\qquad r\in\mathbb Z/w\mathbb Z.                              \tag{2.1}
\]

At transverse Fourier character `theta=2pi l/w`, the periodized transform is

\[
\bar A^{(w)}(z,e^{i\theta})
=\sum_{x,r}\bar A_{x,r}^{(w)}z^x e^{i\theta r}.                \tag{2.2}
\]

For the zero character,

\[
\boxed{
\bar A^{(w)}(z,1)
=\sum_{x,r,k}A_{x,r+kw}z^x
=A(z,1).}                                                      \tag{2.3}
\]

Therefore the Perron root equation of the zero mode is **identical for every circumference**:

\[
\boxed{\bar R_w=R,\qquad \bar\kappa_w=\kappa.}                \tag{2.4}
\]

This is exact.  Folding the transverse displacement modulo `w` is not itself a source of finite-width mass drift.

The same conclusion holds with finite internal memory: only translation invariance in the transverse displacement and the same local piece weights are used.

## 3. Consequence: all mass drift is a wrap-defect / decoration effect

Suppose the actual cylinder kernel `A^{cyl}_w` is not the pure periodization because an irreducible piece or its attached decoration can interact with a periodic copy of itself, changing its validity/weight.  Write schematically

\[
A^{cyl}_w=\bar A^{(w)}+E_w.                                   \tag{3.1}
\]

Any local object whose lifted transverse diameter is strictly smaller than `w/2-O(1)` embeds in the cylinder exactly as it does in the plane.  Therefore `E_w` is supported only on pieces/decorations that reach transverse scale `Omega(w)`.

If under the critical Perron tilt the full irreducible object has an exponential transverse-diameter tail,

\[
P_*(\operatorname{diam}_\perp\ge r)\le C e^{-cr},             \tag{3.2}
\]

with the corresponding weighted kernel estimate, then

\[
\boxed{\|E_w\|\le C'e^{-c'w}}                                \tag{3.3}
\]

in any operator norm compatible with the Perron perturbation argument.

Thus the strong #760 conjecture is reduced to a geometric locality statement about the **full weighted irreducible object**, not to another comparison of fitted finite-width eigenvalues.

## 4. Simple Perron perturbation gives exponential mass locality

Assume `A(R,1)` has a simple isolated Perron eigenvalue one and a nonzero longitudinal derivative

\[
\mu=\partial_{\log z}\log\rho(A(z,1))|_{z=R}>0.                \tag{4.1}
\]

If (3.3) holds uniformly on a neighbourhood of `R`, analytic/Kato perturbation gives

\[
\rho(A_w^{cyl}(R,1))=1+O(e^{-c'w}).                            \tag{4.2}
\]

The implicit-function theorem then shifts the root by the same order:

\[
R_w-R=O(e^{-c'w}),                                             \tag{4.3}
\]

and hence

\[
\boxed{
\gamma_w-\kappa=O(e^{-c'w}).}                                \tag{4.4}
\]

In particular

\[
w(\gamma_w-\kappa)\to0,                                     \tag{4.5}
\]

which is exactly the condition needed to substitute a computable cylinder mass for the plane mass in a leading amplitude diagnostic.

The sign `gamma_w>=kappa`, when separately established by the covering comparison, is compatible with (4.4) and sharpens it to

\[
0\le\gamma_w-\kappa\le Ce^{-cw}.                              \tag{4.6}
\]

## 5. Why a `1/w^2` Brownian confinement correction is the wrong default here

A Brownian path confined between **hard transverse boundaries** has a Dirichlet ground-state cost of order `1/w^2`.  That is a different geometry.

The cylinder is periodic.  For a translation-invariant effective path, the constant transverse Fourier mode survives exactly and has zero transverse Laplacian eigenvalue.  Equation (2.3) is the discrete Markov-additive version of this fact.

Therefore a polynomial `1/w^2` correction should not be attributed to mere diffusive transverse wandering on a periodic cylinder.  It would signal either:

- an additional constraint that effectively imposes a boundary/killing condition;
- a source/readout that removes the zero transverse mode;
- a nonlocal component condition not captured by pure periodization;
- or failure of the assumed massive/local renewal description.

This distinction is useful when interpreting finite-width spectral data.

## 6. Interface to the actual complete winding component

The actual SITE component problem adds two nontrivial layers absent from the abstract kernel:

1. complete-component external-boundary weights;
2. the global condition that the selected object is one connected winding component.

`sewing-with-memory.md` shows that the local weight itself has finite three-column memory, and `matrix-sewing-unit-residue-20260914.md` shows finite memory does not alter the pure cyclic zero-mode residue.  What is still needed for #760 is a representation in which the **full closed component/long arm** has exponentially localized decorations so that changing the plane to a periodic cylinder only changes pieces that see a periodic image.

A suitable proof can follow the architecture of open-connection OZ locality:

- define regeneration pieces in the lift;
- prove exponential tails for piece diameter and decorations under the appropriate tilted/component law;
- show pieces of diameter `<w/2-C` have identical plane/cylinder weights;
- sum the wrap-defect tail to obtain (3.3);
- apply the Perron perturbation above.

This is substantially more specific than the statement `gamma_w -> kappa`.

## 7. Numerical consequence

Finite-width data should not be used to **prove** (4.4), but they can falsify particular remainder scales once `kappa` is independently enclosed.

If the theory is correct, the quantity

\[
e^{cw}(\gamma_w-\kappa)
\]

should remain bounded for some `c>0`, while `w(\gamma_w-\kappa)` must tend to zero.  A stable nonzero `w(\gamma_w-\kappa)` would falsify the required locality for amplitude substitution even if plain convergence remains true.

The appropriate computation is therefore an independent `kappa` interval plus the existing `gamma_w` certificate, exactly as #760 already requests; no extra density fit is needed.

## 8. Claim boundary

Equations (2.3)--(2.4) are exact for the stated translation-invariant Markov-additive kernel.  Equations (3.3)--(4.4) are a conditional perturbation theorem under exponential localization of the full cylinder/plane kernel difference.  Establishing that kernel representation and localization for the actual SITE complete-component tail is the remaining model-specific theorem.
