# Sewing amplitude diagnostic: band multiplicity versus genuine insertion

2026-09-14.  Consequence of `matrix-sewing-unit-residue-20260914.md`.  The goal is to make the remaining SITE prefactor ambiguity falsifiable without confusing three distinct mechanisms: diffusion scale, multiplicity of soft bands, and a genuine component/mark insertion.

## 1. One simple band

For the pure cyclic object

\[
L_w=w[z^wy^0]\{-\log\det(I-A(z,y))\},                         \tag{1.1}
\]

one isolated simple Perron band with root `R=e^kappa` and transverse diffusion `D` gives

\[
L_w\sim\frac{e^{-\kappa w}}{\sqrt{2\pi D w}}.                 \tag{1.2}
\]

There is no extra endpoint/Perron-overlap residue.  Such factors enter only the analytic part of the determinant factorization and disappear from the coefficient of the logarithmic singularity.

Thus, after `D` is defined by the same twisted band, the dimensionless residue is exactly one.

## 2. Several dominant simple bands

Suppose instead that the determinant has finitely many dominant simple bands `lambda_j` with the same real exponential rate `R`, no Jordan singularity, and aperiodic real saddles at `theta=0`, with

\[
\log R_j(\theta)=\log R+\frac12D_j\theta^2+O(\theta^4).         \tag{2.1}
\]

The logarithm splits additively:

\[
-\log\det(I-A)
=\sum_j-\log(1-\lambda_j)+\text{analytic}.                    \tag{2.2}
\]

Each simple band contributes unit logarithmic residue, so

\[
\boxed{
L_w\sim
\frac{e^{-\kappa w}}{\sqrt{2\pi w}}
\sum_j D_j^{-1/2}.}                                           \tag{2.3}
\]

If all dominant bands have the same `D`, the apparent residue relative to one-band normalization is the **integer band multiplicity**.

For periodic/complex saddles, phases and residue-class oscillations must be kept explicitly; averaging them into an arbitrary positive amplitude would lose information.  This is another reason to diagnose lattice periodicity before fitting one scalar `zeta`.

## 3. A genuine insertion changes the amplitude continuously

Now consider a closed coefficient with an insertion,

\[
J_w=[z^wy^0]\operatorname{tr}
\left[B(z,y)(I-A(z,y))^{-1}\right],                            \tag{3.1}
\]

or the derivative of the log determinant with respect to a physical/mark source.  Near a simple Perron pole,

\[
(I-A)^{-1}
\sim\frac{r\,l^T}{1-\lambda},                                 \tag{3.2}
\]

so the singular amplitude contains

\[
l^T B r                                                        \tag{3.3}
\]

and the derivative of the pole location.  These quantities can vary smoothly and nontrivially with `p`.  Unlike the pure logarithmic coefficient, they are not forced to one.

Therefore a nontrivial prefactor function does **not** refute a diffusive transverse mode.  It diagnoses that the actual SITE observable corresponds to an inserted/marked cyclic object rather than the pure unrooted determinant.

## 4. Application to complete-component sewing

The current SITE problem has three logically distinct possibilities.

### A. Pure cyclic determinant

A canonical regeneration state makes each complete component exactly one unmarked closed cycle of a Markov-additive kernel.  Then

\[
\nu_w\sim\frac{e^{-\kappa w}}{\sqrt{2\pi D w}}
\]

(up to an integer/periodic band multiplicity).  `zeta=1` for a single aperiodic band.

### B. Marked/insertion object

A convenient cut, seam, regeneration mark, or external-boundary bookkeeping leaves a source `B` on the cyclic kernel.  Then the `w^{-1/2}` power and the same `D` can survive while the amplitude is a genuine smooth function.  The exact Palm mark-unbiasing formulas in `sewing-with-memory.md` belong here.

### C. No isolated finite/quasi-compact band

The effective state remains infinite with another soft mode, no spectral gap, or a continuum of near-leading states.  Then even the `w^{-1/2}` power or Brownian range law may fail.

This classification is sharper than asking only whether a fitted `beta_eff` is near `1/2`.

## 5. Falsification table

If future rigorous/numerical certificates establish:

- `beta != 1/2`: rule out A and the simple form of B; look for C or a different transverse scaling;
- `beta=1/2`, Brownian/curvature `D` consistent, but a smooth non-unit `zeta(p)`: rule out pure A, favour B;
- `beta=1/2`, `zeta` equal to an integer after common `D` normalization: inspect dominant-band multiplicity/periodicity before introducing a new physical amplitude;
- `beta=1/2`, `zeta=1`, and mark corrections vanish under exact unrooting: A becomes a plausible exact sewing target.

Three-width interpolation cannot distinguish these alternatives.  The right evidence is a certified mass, a separately defined diffusion/curvature scale, and a controlled component/mark normalization.

## 6. Relation to the unit-residue conjecture

The previous branch conjecture `U` asked whether the actual SITE sewing factor might equal one.  The matrix lemma and this diagnostic refine it:

> unit residue is automatic **after** the component has been identified with a pure cyclic log determinant having one simple soft band.

The difficult theorem is therefore the component-to-cyclic-object identification, not a second amplitude calculation once that identification is complete.
