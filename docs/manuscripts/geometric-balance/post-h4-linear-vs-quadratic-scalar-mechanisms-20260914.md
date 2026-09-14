# Post-H4 residual: linear scalar insertion versus quadratic H4×H−4 mixing

Date: 2026-09-14

Status: mechanism split / analysis proposal.  It responds to the latest #808/#47 update and to the warning that a one-point angular zero does not eliminate second-order scalar mixing.  No field identity is asserted.

## 1. Why a nonzero H4-null residual has at least two qualitatively different origins

Suppose the leading square correction is carried by a real H4/spin-four coupling `g4` and an exact angular projector removes the **linear** H4 response.

A residual angular scalar can arise from at least two distinct mechanisms:

### Linear scalar channel

A separate scalar scaling field/source `g0` contributes at first order:

```text
delta O_linear = g0 L^-q0 + ... .
```

The current `x=33/4` / `V_<1,4>` candidate belongs to this class if its coupling is nonzero.

### Quadratic spin-four mixing

Two spin-four insertions can fuse to total spin zero:

```text
(+4)+(-4)=0.
```

A second-order connected response can therefore produce an angular scalar even when the linear H4 one-point term is projected to zero:

```text
delta O_quad
 ~ g4 g_-4 * integral-connected-two-point + contact/thermal counterterms.
```

An H4-null **linear** projector does not remove this scalar second-order contribution.

Thus

```text
post-H4 scalar != automatically new scalar primary.
```

## 2. Radial exponents can distinguish the two stories only after contact terms are typed

If the leading H4 root correction scales as `L^-4`, a naive product of two root amplitudes would suggest `L^-8`.  That is generally **not** the correct second-order field-theory exponent because integrated insertions bring volume powers, OPE singularities, thermal retuning and contact counterterms.

The right object is the connected second derivative of the physical topological free-energy/root functional with respect to a genuine H4 microscopic coupling.

Schematically, for a source `g` with first-order H4 character,

```text
Theta_gg
 = connected Green--Kubo / integrated two-point
   + explicit contact term
   + source-normalization term.
```

The scalar component of `Theta_gg` is the quadratic adversary to a new linear scalar field.

Therefore a measured `L^-7` residual cannot be accepted/rejected as quadratic H4 mixing by dimensional multiplication alone.

## 3. A clean lattice experiment: sign reversal of the H4 source

Let a microscopic anisotropy coupling `g` reverse the sign of the leading H4 insertion while preserving all scalar bare couplings.  Then expand the root or charge free energy:

```text
R(g)=R0 + a1 g + a2 g^2 + a3 g^3 + ... .
```

Form

```text
R_odd(g)  =[R(g)-R(-g)]/2 = a1 g+a3 g^3+...,
R_even(g) =[R(g)+R(-g)]/2-R(0)=a2 g^2+a4 g^4+....
```

The linear H4 channel lives in `R_odd`; a quadratic H4×H−4 scalar contribution lives in `R_even`.

A pre-existing linear scalar field independent of `g` instead contributes to `R0` and does not grow as `g^2` under this source.

This supplies a mechanism discriminator that does not require a gigantic same-circle angular tomography.

## 4. Use the normalized-source formalism

The source `g` must be a normalized physical perturbation.  For an unnormalized transfer source, include the row pressure/normalizer before differentiating.

At a root, the first derivative is the residualized covariance response.  The second derivative must include:

```text
integrated score covariance,
explicit second derivative/contact of log weight,
thermal/root-motion terms.
```

Do not square one-point amplitudes or multiply two first derivatives as a surrogate.

This is exactly the failure mode that the #802 normalization audit warns against.

## 5. Relation to angular projection

There are now two independent decompositions:

```text
geometry : H0, H4, H8, ...
source sign : even/odd under g -> -g.
```

A genuine new linear scalar correction is

```text
angular H0,
source-linear in its own scalar coupling.
```

Quadratic H4 mixing is

```text
angular H0,
even/quadratic in the H4 source.
```

Therefore angular scalarity alone cannot distinguish them.

## 6. A practical bounded test before N377/N1105-scale production

Use an existing small-width safe operator for which a physical orientation/anisotropy source can be inserted without changing the state space.  At two or three small amplitudes `±g` compute:

```text
Theta(pc,g),
Theta_p(pc,g),
p_root(g),
```

with the same normalization conventions and tight solver.

Check:

```text
odd part / g  -> existing H4 response;
even part/g^2 -> finite nonzero limit or zero within error.
```

If the quadratic scalar term is already sizeable and has a radial trend compatible with the post-H4 residual, the simple “new V_<1,4> linear field” narrative must remain mixed.

If the even second-order coefficient is tiny/forbidden while #808 finds a robust H0 `ell^-7`, the linear scalar interpretation becomes substantially stronger.

## 7. Continuum/OPE question after the lattice test

Only after observing a nonzero quadratic scalar response should one ask which OPE channels of the H4 correction feed it.  At `c=0` one must keep logarithmic collisions/contact terms; equality of total dimensions does not by itself determine the fusion coefficient.

Conversely a zero finite response for one microscopic H4 source does not prove all quadratic H4 mixing vanishes: it can be source/map specific.

## 8. Claim boundary

Exact algebra:

- source-sign even/odd decomposition;
- a linear H4 projector does not eliminate second-order H4×H−4 scalar response.

Programme:

- measure the connected normalized second response of an actual H4 source;
- use it as an adversary to the post-H4 linear-scalar/V14 interpretation.

This should be resolved before a pure `ell^-7` angular residual is promoted to a unique scalar operator.