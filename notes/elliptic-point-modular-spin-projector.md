# Elliptic-point modular spin projector

Status: **C0 theory/design**.

Related: #103, #138, #57.

## Core claim

The hexagonal/equianharmonic torus should be treated as an exact continuum spin selector, not merely as a place where `E4(tau)` happens to vanish.

For a scalar/topological torus observable `Q`, suppose a first-order finite-size correction comes from a continuum field of scaling dimension `x` and conformal spin `s`:

\[
\delta Q_s(L,\tau)=g_s L^{2-x}F_s(\tau,\bar\tau)+\cdots.
\]

At a torus automorphism with rotation angle `phi`, the observable is invariant but the spin-s insertion acquires `exp(i s phi)`. Therefore

\[
F_s(\tau_*)=e^{is\phi}F_s(\tau_*),
\]

so

\[
\boxed{F_s(\tau_*)=0\quad\text{unless}\quad e^{is\phi}=1.}
\]

This does not assume any particular Eisenstein-series formula.

## Square versus hexagonal elliptic points

At `tau=i`, the order-4 automorphism has `phi=pi/2`, so a scalar response keeps only

\[
s=0\pmod4.
\]

At

\[
\tau_\hex=-\frac12+i\frac{\sqrt3}{2},
\]

the order-6 automorphism has `phi=pi/3`, so a scalar response keeps only

\[
s=0\pmod6.
\]

Square-lattice anisotropy supplies `s=0 mod 4`. At the hexagonal point the intersection is therefore

\[
\boxed{s=0\pmod{12}.}
\]

Hence, to first order in the corresponding irrelevant coupling,

```text
H4   -> forced zero
H8   -> forced zero
H12  -> allowed
H16  -> forced zero
H20  -> forced zero
H24  -> allowed
...
```

This gives an operator-spin discriminator independent of the norm-5 Gaussian test. In particular, H4 and H12 behave qualitatively differently at the elliptic point.

## Generic zero order

The stabilizer also rotates the tangent space of moduli. For H4 the constant term is forbidden, while one power of the shape deformation can carry the compensating character. Thus a simple zero is symmetry-allowed and is the generic expectation:

\[
F_4(\tau_\hex+\delta\tau)
=K\,\delta\tau+\bar K\,\delta\bar\tau+O(|\delta\tau|^2).
\]

If an additional selection rule removes this linear term, suppression will be even stronger. The frozen first model should be a simple zero rather than a fitted arbitrary zero order.

The same geometric spin selection applies to an LCFT Jordan partner with the same spin. Logarithmic mixing modifies radial scale dependence but does not rescue a forbidden spin character at the elliptic point.

## Pell approach from both sides

Choose square-lattice periods

\[
v_1=(2q,0),\qquad v_2=(-q,p),
\]

with

\[
p^2-3q^2=\eta,\qquad \eta\in\{1,-2\}.
\]

Then

\[
\tau_{p,q}=-\frac12+i\frac{p}{2q},
\qquad
N=|\det(v_1,v_2)|=2pq.
\]

Define

\[
\delta_{p,q}=\operatorname{Im}\tau_{p,q}-\frac{\sqrt3}{2}.
\]

Using the Pell equation,

\[
\delta_{p,q}
=\frac{\eta}{2q(p+\sqrt3 q)},
\]

and therefore

\[
\boxed{N\delta_{p,q}=\frac{\eta p}{p+\sqrt3 q}\to\frac{\eta}{2}.}
\]

So the two exact subsequences approach the modular zero from opposite sides with shape error `O(N^-1)`:

```text
eta=+1: N*delta -> +1/2
eta=-2: N*delta -> -1
```

Useful sizes are

```text
eta=+1:
(p,q)=(7,4)       N=56
       (26,15)    N=780
       (97,56)    N=10864
       (362,209)  N=151316

eta=-2:
(p,q)=(5,3)       N=30
       (19,11)    N=418
       (71,41)    N=5822
       (265,153)  N=81090
```

Each subsequence advances by multiplication of `p+q sqrt(3)` by `2+sqrt(3)`.

## Parameter-free normalization

If a spin-4 contribution has the fixed-shape form

\[
X_N\sim A N^{-\alpha}F_4(\tau_N),
\]

then along the Pell sequence a simple modular zero gives

\[
X_N\sim A K\,\delta_N N^{-\alpha}
\sim \widetilde A N^{-(\alpha+1)}.
\]

The shape-normalized amplitude

\[
\boxed{\mathcal A_\hex(N)=\frac{N^\alpha X_N}{\delta_N}}
\]

should converge to the same constant from both sides.

Equivalently,

\[
N^{\alpha+1}X_N\to C\frac{\eta}{2},
\]

so the asymptotic side amplitudes obey

\[
\boxed{A_{+1}/A_{-2}\to-1/2.}
\]

This ratio cancels the unknown microscopic spin-4 coupling and the unknown derivative of the continuum shape function.

## Current candidate exponents

### Matching-even identity-like H4

For `x=4`, `alpha=1`.

Generic fixed-shape law:

```text
N^-1 = L^-2
```

Pell-to-hex simple-zero law:

```text
N^-2 = L^-4
```

### Matching-odd thermal H4

For `x=21/4`, `alpha=13/8`.

Generic central residual:

```text
N^-13/8 = L^-13/4
```

Pell-to-hex simple-zero law:

```text
N^-21/8 = L^-21/4
```

Since the thermal slope scales as `N^(3/8)`, the induced root shift becomes

\[
\boxed{\Delta p^*_{\hex}\sim N^{-3}=L^{-6}.}
\]

Thus a magic-torus sequence can turn the present `L^-4` root mechanism into a parameter-free `L^-6` superconvergence test.

### H12 adversary

A genuine spin-12 contribution is allowed by the hexagonal stabilizer and does not gain the extra `N^-1`. Therefore

```text
H4-like sector  -> exponent gains +1 in N
H12-like sector -> no elliptic-point acceleration
```

This is orthogonal to #57, which separates the harmonics through Gaussian-multiplier phase/magnitude.

## Clean first control: exact square bond

Use square-bond percolation at exact `p_c=1/2` and a modular-invariant topology observable (`cross`, `either/nontrivial`, or another homology combination) whose continuum critical value `R_infinity(tau)` can be computed from the Pinson/Ziff torus formulas.

For each Pell torus measure

\[
E_N(\tau_N)=R_N(1/2;\tau_N)-R_\infty(\tau_N).
\]

This avoids fitting a shape-dependent continuum intercept.

Frozen first tests:

1. generic square-lattice anisotropy has an `N^-1` component away from the elliptic point;
2. Pell-to-hex H4 response is proportional to `delta_N N^-1`;
3. `N E_N/delta_N` approaches one common amplitude from both sides;
4. `N^2 E_N` has the asymptotic side-amplitude ratio `-1/2`.

A pass would justify promoting the same geometry to square-site production.

## Square-site target

After the exact-threshold control:

1. test the central matching residual against
   `M_N ~ delta_N N^-13/8`;
2. compare with an unsuppressed H12-like `N^-13/8` adversary;
3. score the pseudo-critical root against `N^-3` versus the ordinary `N^-2` law;
4. report `N^(13/8) M_N/delta_N` and `N^2 DeltaRoot/delta_N` as the two-sided collapse variables.

If the central residual and root both accelerate, that is strong evidence that the existing root-moving field is genuinely in the H4 spin class. If H4 is suppressed but an `N^-13/8` residual remains, that is a direct reason to search the H12/multiple-of-12 sector rather than merely adding free radial corrections.

## Research significance

This changes torus-modulus spectroscopy from a broad shape scan into a symmetry-selected experiment. The hexagonal point can simultaneously provide:

- an independent H4/H12 discriminator;
- a parameter-free exponent acceleration;
- a two-sided sign test;
- a potentially superconvergent matching-root estimator;
- a control of whether the current finite-size mechanism is truly a first-order spin-4 coupling.

## Literature anchors

- Pinson (1994): critical torus homology probabilities.
- Feng, Deng, Blote (2008): orientation-sensitive square-lattice correction sector.
- Gaberdiel and Lang (2009): modular covariance and modular differential equations for torus one-point amplitudes.
- Roux, Ribault and Jacobsen (2026), arXiv:2604.24491: modular-covariant torus one-point functions in critical loop models, including logarithmic blocks.

The earlier `E4(tau_hex)=0` observation remains a useful example of the same stabilizer selection rule, but is no longer the foundation of the proposal.
