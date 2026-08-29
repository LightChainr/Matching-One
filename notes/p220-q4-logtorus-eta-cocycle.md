# A gauge-free eta-cocycle fingerprint for the Q4 logarithmic torus block

Issue #220 asks for a modulus fingerprint stronger than an affine `log L`.  The full top-field intercept is basis-dependent, but one part of it can already be isolated exactly from the generic Potts energy block.

## 1. Differentiate the block, not an ansatz

The generic Potts energy one-point block is

```text
F_h(tau)=eta(tau)^(2h),
```

so the full non-chiral block is `|eta|^(4h)`.  A homogeneous level-4 torus Ward operation supplies a weight-4 factor.  Along any generic normalization path write

```text
H_Q4(u,tau)=C_Q4(u) g2(tau) |eta(tau)|^(4h(u)).
```

Therefore

```text
R_u(tau) := partial_u H_Q4/H_Q4
          = partial_u log C_Q4 + 4 h'(u) log|eta(tau)|.
```

Since the scaling dimension is `x=2h`, the modulus difference obeys

```text
[R_u(tau1)-R_u(tau0)]/x'(u)
 = 2 log|eta(tau1)/eta(tau0)|.                 (1)
```

Equation (1) is independent of the velocity, the Q4 normalization path, and the top-field basis shift

```text
H_top -> H_top + alpha H_bottom.
```

The shift adds the same constant to `H_top/H_bottom` at every modulus and disappears in the difference.  This is the useful observable; the absolute top intercept is not.

## 2. Audit of the fixed repository Q4 path

For the level-2 degenerate family,

```text
(L_-2 - 3/[2(2h+1)] L_-1^2)|h>=0.
```

Torus translation then gives

```text
<L_-2^2>/<L_-4> = 3/(2h+1),
<L_-3 L_-1>/<L_-4> = -2.
```

Keeping the repository coefficients `Q4=40L_-2^2-60L_-3L_-1-9L_-4` fixed,

```text
C_Q4(h)=h[111+120/(2h+1)]/20.
```

At `h=5/8`,

```text
C_Q4=493/96,
d_h C_Q4=3637/540,
d_h log C_Q4=29096/22185.
```

The first value independently recovers the existing Ward identity.  The derivatives depend on how Q4 is continued away from `Q=1`, which is exactly why they are removed from (1).

## 3. Three CM moduli give an exact rational discriminator

At the square CM point, theta/eta identities give

```text
eta(2i)/eta(i)=2^(-3/8).
```

One derivation uses

```text
theta2(tau)=2 eta(2tau)^2/eta(tau),
lambda(i)=theta2(i)^4/theta3(i)^4=1/2,
eta(i)^3=theta2(i)theta3(i)theta4(i)/2,
theta2(i)=theta4(i).
```

For the modularly related shear

```text
tau_s=(1+i)/2=i/(i+1),
```

Dedekind covariance gives

```text
|eta(tau_s)|/eta(i)=|1+i|^(1/2)=2^(1/4).
```

Define the gauge-free normalized difference

```text
Xi(tau1,tau0)
 = [H_top/H_bottom(tau1)-H_top/H_bottom(tau0)]/x_collision_velocity.
```

For the energy-block derivative mechanism,

```text
Xi(2i,i)             = -3 log(2)/4,
Xi((1+i)/2,i)        =  log(2)/2,
Xi(2i,i)/Xi((1+i)/2,i) = -3/2.                (2)
```

The rational `-3/2` is the sharpest target: it cancels the unknown top amplitude, bottom admixture, collision velocity, and common lattice normalization.  It uses one genuine rectangular deformation and one exact modular/shear cocycle.

## 4. The hexagonal node is not lifted by this derivative

Because the complete differentiated energy block retains the weight-4 factor,

```text
H_bottom(rho)=0,
partial_u H_bottom(rho)=0,
rho=exp(i pi/3).
```

Thus the simplest collision derivative does **not** generate a finite top intercept at the `E4` zero.  A reproducible nonzero top signal at `rho` would demand an additional logarithmic torus block, a different combinatorial-map sector, or a lattice field outside this energy-block derivative.  It cannot be repaired by changing the top/bottom basis.

## 5. Minimal score

Use the same typed bottom and top readouts at `i`, `2i`, and `(1+i)/2`.  Reconstruct `H_top/H_bottom` inside each common replicate, then score the single vector relation

```text
2 [R(2i)-R(i)] + 3 [R((1+i)/2)-R(i)] = 0.
```

The individual differences should also have opposite signs.  The hexagonal geometry is a separate joint-zero/null-sector discriminator, not part of the ratio because `H_bottom(rho)=0`.

## Boundary

This identity is exact for the derivative of the generic Potts energy eta block after a homogeneous level-4 Ward operation.  The bridge to a lattice observable remains conditional.  An arbitrary Q4 Jordan top may contain extra torus blocks; equation (2) is designed to detect precisely that failure rather than absorb it into a free intercept.
