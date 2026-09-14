# Poisson tessellation consequence for alternating black/white essential components

2026-09-14.  This note extracts a concrete consequence of the existing component-process Chen--Stein estimate together with the componentwise black/white alternation recorded in `structural-consequences-20260914.md`.  It is meant to close the bookkeeping gap in #764 without invoking an abstract Palm-convergence theorem as a black box.

Status: the probability input is still the author-level `poisson-birth-windows.md` proof on the #739 branch.  The argument below is an explicit consequence of that input, not a new independent validation of it.

## 1. Set the torus height so the expected count diverges slowly

Fix `p < p_c(NN)` and write `nu_w=nu_w^4(p)` for the once-per-complete-black-winding-component row intensity.  The existing rate theorem gives

\[
-\frac1w\log\nu_w\to\kappa_4(p)>0.
\]

Choose any deterministic sequence

\[
\lambda_w\to\infty,
\qquad
\log\lambda_w=o(w),
\]

and take

\[
m_w=\left\lfloor\lambda_w/\nu_w\right\rfloor.
\]

Then `log m_w / w -> kappa_4(p)`.  In the fixed-`p` specialization of the Chen--Stein bound, the pair term has the form

\[
m_w\,\operatorname{poly}(w)e^{-2\kappa_4(p)w+o(w)}
=\lambda_w\operatorname{poly}(w)e^{-\kappa_4(p)w+o(w)}\to0,
\]

and the height-localization term is smaller.  Thus the entire black anchor point process on the cyclic height coordinate is close in total variation to a homogeneous Poisson process of mean `lambda_w`, even though the mean itself diverges.

Rescale the cyclic vertical coordinate by `nu_w`; the circle then has length `lambda_w+o(1)` and the comparison Poisson process has unit intensity.

## 2. Uniform-anchor gap without abstract Palm convergence

Given a finite point configuration on the circle, apply the following Markov kernel:

1. if there is no point, return a cemetery symbol;
2. otherwise choose one of the points uniformly;
3. return the clockwise distance from that point to the next point, in the `nu_w`-rescaled coordinate.

Total variation distance contracts under a Markov kernel.  Therefore the output law for the percolation anchor process is at most the existing point-process TV error away from the corresponding output law for the Poisson process.

For a homogeneous Poisson process on a circle of length `lambda`, conditional on `N=n>=2`, the `n` cyclic spacings divided by `lambda` have the Dirichlet `(1,...,1)` law.  A uniformly selected spacing is therefore

\[
\lambda\,B_{1,n-1},
\]

where `B_{1,n-1}` is Beta `(1,n-1)`.  With `N~Poi(lambda)` and `lambda->infinity`,

\[
\lambda B_{1,N-1}\Rightarrow Exp(1),
\]

while `P(N<2)->0`.  Consequently, if `G_w` is the gap following a uniformly selected black winding-component anchor,

\[
\boxed{\nu_w G_w\Rightarrow Exp(1).}
\]

This proves the gap law by a finite-circle kernel argument.  No renewal independence of the percolation components is assumed.

The exact stationary mass-transport identity remains stronger at the first-moment level:

\[
E^{Palm}G_w=1/\nu_w
\]

for every fixed `w`.  Hence the mean does not need to be recovered by uniform integrability from the weak limit.

## 3. White component span is the same gap at leading scale

Let `B_i,B_{i+1}` be consecutive black essential components and `W_i` the unique intervening white matching essential component.  The componentwise complement lemma gives a fixed lattice constant `C_0` such that

\[
|L(W_i)-G_i|\le L(B_i)+L(B_{i+1})+C_0.
\]

The existing fixed-subcritical component-Palm volume bound gives `E L(B_i)=O(w)`.  Since `nu_w` is exponentially small,

\[
\nu_w L(B_i)\to0
\]

in `L^1`, hence also in probability.  Slutsky therefore yields

\[
\boxed{\nu_w L(W_i)\Rightarrow Exp(1).}
\]

Taking expectations in the deterministic comparison and using the exact mean gap identity gives independently

\[
\boxed{\nu_w E L(W_i)\to1.}
\]

The mean statement does not rely on the Poisson approximation at all.

## 4. Uniform-location protocol: two independent exponential sides

There is a second useful Markov kernel.  Pick an independent uniform vertical location on the rescaled circle and record the distances `D_-` and `D_+` to the nearest black anchors below and above it.

For a unit-rate Poisson process on a circle whose length tends to infinity,

\[
(D_-,D_+)\Rightarrow(E_1,E_2),
\]

where `E_1,E_2` are independent `Exp(1)`.  Total variation contraction transfers this to the black barrier process.  Therefore

\[
\boxed{\nu_w(D_-+D_+)\Rightarrow Gamma(2,1).}
\]

Moreover

\[
\boxed{U=\frac{D_-}{D_-+D_+}\Rightarrow Uniform(0,1),}
\]

and for independent exponentials the ratio is independent of the total.  Thus the stationary-location version of the white gap has a joint parameter-free prediction:

\[
(\nu_w L_{\rm containing},U)
\Rightarrow(Gamma(2,1),Uniform(0,1))
\]

with asymptotic independence, up to the same vanishing black-component thickness error.

This `U` is a **white-gap location coordinate**.  It should not be conflated with the #762 branch/core coordinate until a specific map between the two sampling protocols is proved.

## 5. Higher spacings

The same uniform-anchor kernel can return the sum of the next `k` cyclic gaps.  For fixed `k`, the Poisson limit gives

\[
\boxed{\nu_w(G_i+\cdots+G_{i+k-1})\Rightarrow Gamma(k,1).}
\]

Thus the entire fixed-order spacing hierarchy follows from the already available point-process approximation; it is not a collection of new phenomenological fits.

## 6. Scaling-limit interpretation

At fixed black-subcritical `p`, the two-colour essential geometry has a natural coarse limit under vertical scaling by `nu_w`:

- black essential components have vanishing scaled thickness and become the points of a unit-rate Poisson process;
- the unique white matching essential component between consecutive black barriers occupies, up to a vanishing boundary error, the complementary Poisson interval.

This is a **Poisson interval tessellation** of the vertical line, not two independent black/white Poisson clouds.  It packages the Exp component-Palm span, Gamma location-biased span, uniform relative position, and higher-gap Gamma laws in one object.

## 7. Remaining publication-grade checks

Only two local presentation items remain beyond the existing author-level probability input:

1. state the regular-neighbourhood/componentwise white-complement lemma with one explicit universal lattice constant `C_0`;
2. quote the precise process-TV version of the Chen--Stein estimate and note explicitly that the bound remains `o(1)` for `lambda_w->infinity` with `log lambda_w=o(w)`.

No new Monte Carlo, width scan, or independent-colour approximation is needed for this consequence.
