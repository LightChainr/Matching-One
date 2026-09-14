# Pell approximants to the spin-four nodal direction: an alternating w^-6 charge-root test

2026-09-14.  Sharp arithmetic/CFT falsification target for the sector-odd spin-four anisotropy hypothesis.

The square-lattice hypothesis predicts an orientation-dependent semi-infinite charge-root shift

\[
p_{\ell}^{ch}(\theta)-p_c
\sim -A_4\cos(4\theta)\,\ell^{-4},                            \tag{1}
\]

where `ell` is the physical Euclidean circumference and `theta` is the angle of the short period relative to a lattice axis.  Reflection fixes a cosine rather than a general phase.

The spin-four nodes are at

\[
\theta=\pi/8,3\pi/8,\ldots,                                  \tag{2}
\]

but `tan(pi/8)=sqrt(2)-1` is irrational, so no nonzero integer period lies exactly on the node.  This makes Pell approximants a feature, not a nuisance: their angular error is precisely strong enough to convert the `ell^-4` law into a clean alternating `ell^-6` prediction.

## 1. Integer spin-four harmonic

For a primitive integer period

\[
u=(a,b),
\qquad
\ell^2=a^2+b^2,                                               \tag{1.1}
\]

the real spin-four harmonic is exactly

\[
\boxed{
\cos4\theta
=\frac{a^4-6a^2b^2+b^4}{(a^2+b^2)^2}.}                 \tag{1.2}
\]

Factor the numerator:

\[
\boxed{
a^4-6a^2b^2+b^4
=(a^2-2ab-b^2)(a^2+2ab-b^2).}                                \tag{1.3}
\]

The nodal equation is therefore

\[
a^2-2ab-b^2=0,                                               \tag{1.4}
\]

whose positive slope is `b/a=sqrt(2)-1`.

## 2. Pell sequence

Choose coprime positive integer solutions of

\[
\boxed{a_n^2-2a_nb_n-b_n^2=(-1)^n}                           \tag{2.1}
\]

(up to an indexing shift).  Equivalently, with

\[
x_n=a_n-b_n,                                                 \tag{2.2}
\]

we have the standard Pell equations

\[
x_n^2-2b_n^2=\pm1.                                          \tag{2.3}
\]

One convenient sequence is

\[
(a,b)=(2,1),(5,2),(12,5),(29,12),(70,29),(169,70),\ldots      \tag{2.4}
\]

with alternating sign in (2.1), and

\[
\frac{b_n}{a_n}\to\sqrt2-1.                                 \tag{2.5}
\]

## 3. Exact angular leakage scale

On the Pell sequence, (1.3) and (2.1) give

\[
\ell_n^2\cos4\theta_n
=(-1)^n
\frac{a_n^2+2a_nb_n-b_n^2}{a_n^2+b_n^2}.                    \tag{3.1}
\]

Taking the nodal ratio limit,

\[
\boxed{
\ell_n^2\cos4\theta_n\longrightarrow(-1)^n\sqrt2.}          \tag{3.2}
\]

Numerically the right side is already visible:

```text
(a,b)       ell^2 cos(4 theta)
(2,1)       -1.4000000000
(5,2)       +1.4137931034
(12,5)      -1.4142011834
(29,12)     +1.4142131980
(70,29)     -1.4142135516
(169,70)    +1.4142135621
```

Thus the best Diophantine approach to the spin-four node naturally supplies a `1/ell^2` residual harmonic with an alternating universal geometric coefficient.

## 4. Predicted pseudo-critical root

Insert (3.2) into the spin-four root law (1):

\[
p_n^{ch}-p_c
\sim
-A_4\cos4\theta_n\,\ell_n^{-4}.                              \tag{4.1}
\]

Then

\[
\boxed{
p_n^{ch}-p_c
\sim (-1)^{n+1}A_4\sqrt2\,\ell_n^{-6}.}                      \tag{4.2}
\]

So the Pell sequence predicts THREE simultaneous signatures:

1. apparent correction exponent `6` rather than `4`;
2. sign alternation with Pell parity;
3. after multiplying by `(-1)^{n+1} ell^6`, the amplitude should approach `sqrt(2) A_4`, where `A_4` is the axial spin-four amplitude in physical-circumference normalization.

A generic scalar correction cannot reproduce the parity alternation tied to the angle error.

## 5. Diagonal direction gives an even simpler sign test

Before attempting the nodal Pell sequence, use the exact diagonal direction

\[
\theta=\pi/4.                                                 \tag{5.1}
\]

Then

\[
\cos4\theta=-1.                                              \tag{5.2}
\]

Thus the spin-four hypothesis predicts a leading **sign reversal** relative to the axial cylinder:

\[
\boxed{
\operatorname{sign}(p_{diag}^{ch}-p_c)
=-\operatorname{sign}(p_{axis}^{ch}-p_c).}                   \tag{5.3}
\]

Since the observed axial sequence has `p_axis^ch < p_c`, a diagonal short period should have

\[
\boxed{p_{diag}^{ch}>p_c}                                    \tag{5.4}
\]

for large circumference.

If the integer diagonal period is `n(1,1)`, its physical circumference is

\[
\ell=\sqrt2 n.                                                \tag{5.5}
\]

Hence equal physical spin-four amplitude predicts the raw-integer scaling

\[
p_{diag}^{ch}-p_c
\sim +\frac{A_4}{4n^4},                                      \tag{5.6}
\]

versus `-A_4/n^4` on the axis.

This is the cheapest direct falsification target.

## 6. Why the Pell node is more discriminating than an ordinary w^-6 fit

The square axial root already has subleading powers `w^-6,w^-8,...`.  An ordinary fit cannot tell whether the `w^-6` coefficient is

- a higher descendant of the same spin-four family;
- a different scalar sector-odd field;
- an analytic lattice correction.

Along the nodal Pell sequence, every correction carrying the SAME spin-four angular factor receives the extra

\[
\cos4\theta_n=O(\ell_n^{-2}).                                \tag{6.1}
\]

Thus an axial spin-four descendant that would produce `ell^-6` becomes `ell^-8` on the Pell node.  The leading `ell^-6` term (4.2) instead comes from angular leakage of the **leading** spin-four field itself.

If an additional orientation-independent sector-odd contribution also exists at exponent six, the scaled Pell data should split as

\[
\ell^6(p_n^{ch}-p_c)
=B_6+(-1)^{n+1}\sqrt2 A_4+o(1).                              \tag{6.2}
\]

Then averaging even/odd subsequences isolates `B_6`, while their difference isolates `A_4`.  This makes the Pell design a direct detector of an otherwise hidden scalar/angle-independent correction.

## 7. Period geometry for a clean cylinder

For each Pell short vector `u_n=(a_n,b_n)`, one clean finite-torus realization is

\[
\Lambda_{n,m}
=\langle
(a_n,b_n),
 m(-b_n,a_n)
\rangle.                                                      \tag{7.1}
\]

The long period is Euclidean-orthogonal to the short one, and

\[
N=m(a_n^2+b_n^2)=m\ell_n^2.                                  \tag{7.2}
\]

Taking `m` sufficiently large produces the semi-infinite-cylinder charge root without introducing a continuum shear.  The physical NN interaction is NOT rotated; only the quotient period is oblique, exactly as required by the directional programme.

For a direct transfer implementation it may be computationally preferable to use a Bezout longitudinal step plus a shear bookkeeping variable.  The final root should be independent of that transfer basis after the physical circumference/orientation are declared correctly.

## 8. Interface to #765 and the projective-slope controls

The directional mass work on this branch already proves that fixed or converging integer directions are legitimate probability objects and that nonparallel period classes separate in exponential elongation.

The present test is finer: it concerns the **irrelevant finite-width correction to the critical charge-sector crossing**, not the leading directional mass.  The spin-four phase is the same embedded harmonic

\[
\frac{(a+ib)^4}{(a^2+b^2)^2}                                 \tag{8.1}
\]

used in the projective-slope positive control, but the two observables should not be identified.  They merely share the same geometric spin transformation.

## 9. Claim boundary

The Pell algebra (1.2)--(3.2) is exact.  The sign reversal and alternating `ell^-6` root law are consequences of the sector-odd spin-four hypothesis, not established percolation theorems.  Their main value is falsifiability: an oblique semi-infinite charge transfer can reject the spin-four mechanism without relying on another lattice or on a multi-parameter CFT fit.
