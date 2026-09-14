# Generic-Q normalization target for the mixed thermal/vacuum spin-four response

2026-09-14.

Status: exact Virasoro normalization algebra plus a concrete generic-Q programme. This note does not identify the measured safe-root H8 response with a particular continuum tensor; it specifies what must be computed before such an identification is meaningful at c=0.

## 1. Two parent spin-four sectors are now separately visible on the lattice

The finite evidence points to two different square-spin-four corrections:

```text
T4: matching odd, thermal-family candidate, omega=13/4,
I4: matching even/common, omega=2.
```

Their radial powers and matching parities were already visible in the historical threshold-rank projectors. New safe-root Gaussian controls additionally resolve a q=2 spin-eight response consistent with a second-order `T4 x I4` channel.

The temptation is to multiply two ordinary c=0 descendants. That is not a well-defined normalization prescription because both relevant c=0 families are degenerate/logarithmic.

## 2. Vacuum level-four spin-four quasiprimary has norm proportional to c

At generic central charge, use the vacuum quasiprimary

```text
Lambda4 = (L_-2^2 - 3/5 L_-4)|0>.
```

The level-four vacuum Gram entries are

```text
<L_-2^2|L_-2^2> = c(c+8)/2,
<L_-2^2|L_-4>   = 3c,
<L_-4|L_-4>     = 5c.
```

Therefore

```text
<Lambda4|Lambda4>
 = c(c+8)/2 - (6/5)(3c) + (9/25)(5c)
 = c(5c+22)/10.
```

In particular,

```text
<Lambda4|Lambda4> ~ (11/5)c
```

as `c->0`. The ordinary vacuum spin-four state is zero-norm at percolation.

## 3. Explicit Q->1 slope of the vacuum norm

Use the standard critical FK-Potts parametrization

```text
Q = 4 cos^2[pi/(m+1)],
c = 1 - 6/[m(m+1)].
```

Percolation is `m=2`, `Q=1`, `c=0`. Direct differentiation gives

```text
dc/dQ |_(Q=1) = 5 sqrt(3)/(4 pi).
```

Hence

```text
<Lambda4|Lambda4>
 = [11 sqrt(3)/(4 pi)] (Q-1)
   + O((Q-1)^2).
```

The coefficient is approximately `1.5161544624` in the displayed state normalization.

Thus the even-spin4 branch has an explicit simple zero in its ordinary vacuum Gram norm.

## 4. Thermal Q4 has finite relative norm but inherits the bottom-field normalization problem

The repository exact checker at `c=0,h=5/8`, with the thermal bottom state normalized to one, gives

```text
Q4 = 40 L_-2^2 - 60 L_-3 L_-1 - 9 L_-4,
<Q4|Q4> / <h|h> = 4930.
```

So the level-four descendant itself is non-null in the ordinary Kac quotient. However the physical percolation energy/Kac field is a zero-norm bottom field of a logarithmic multiplet at c=0. Therefore any physical normalization zero of the bottom field is inherited by Q4 with a finite relative factor.

This separates two questions that were previously conflated:

```text
ordinary descendant exists and is non-null modulo Kac nulls: YES;
physical c=0 field has a non-singular standalone two-point normalization: NO.
```

## 5. The mixed H8 response needs a joint Q->1 limit

At generic Q define physical fields/couplings schematically as

```text
u_T(Q) * T4(Q),
u_I(Q) * I4(Q).
```

The second-order spin-eight response of a declared observable `O` contains

```text
u_T nu_I * R8(Q),

R8(Q)
 ~ integral integral
    < O T4^(+4) I4^(+4) >_conn
 + contact/counterterm contributions.
```

At Q=1 the individual Gram norms vanish or collide with logarithmic partners. The finite lattice coefficient therefore depends on the combined limit of

```text
lattice coupling normalization,
field two-point normalization,
mixed OPE coefficient,
contact subtraction,
possible Jordan collision coefficients.
```

There is no justification for obtaining this limit by multiplying two separately normalized ordinary c=0 one-point functions.

## 6. What generic Q should compute

The minimal useful calculation is a table for Q near one containing

```text
c(Q),
h_energy(Q),
N_T(Q)=two-point normalization of the thermal Q4 branch,
N_I(Q)=c(5c+22)/10 for the vacuum Q4 branch,
C_8(Q)=declared spin8 mixed matrix element/OPE coefficient,
C_0(Q)=declared spin0 mixed matrix element/OPE coefficient,
sector/projector factors for the safe magnetic/rank observable.
```

Then determine whether

```text
C_8 / sqrt(N_T N_I),
C_0 / sqrt(N_T N_I)
```

remain finite, diverge, or vanish as Q->1, and how the microscopic lattice couplings compensate those behaviours.

A `1/(Q-1)` or derivative collision would naturally generate logarithmic/contact terms. A regular limit would support an ordinary mixed-response interpretation.

## 7. Direct interface to the measured safe-root coefficients

The new deterministic root controls give, in one fixed lattice normalization,

```text
c46 ~= -0.116  [same-H4 q=2 dressing],
c88 ~= -0.138  [spin8 q=2 response].
```

Equal-circumference magnetic-gap controls separately resolve the common omega=2 parent into a dominant scalar coefficient and a smaller but nonzero spin-four coefficient.

This means a generic-Q calculation no longer needs to predict an arbitrary finite-size curve. It can target two specific dimensionless response ratios after the lattice/CFT normalization dictionary is declared:

```text
R46 = mixed(T4,S0) / [linear T4 * linear S0],
R88 = mixed(T4,I4) / [linear T4 * linear I4].
```

These ratios cancel the microscopic coupling amplitudes in the ideal scaling-field factorization and are therefore substantially more identifying than another exponent fit.

At present they should **not** be estimated by naively dividing coefficients from different observables or moduli; the same-observable Feynman--Hellmann normalization is required first.

## 8. Consequence for the E8 modular guess

Because the vacuum spin-four branch becomes zero-norm/logarithmic at c=0, the simple ordinary-module statement

```text
spin8 mixed torus shape proportional to E8=E4^2
```

is only a positive control. The actual percolation limit may acquire logarithmic, derivative, contact, or map-sector terms.

The right strategy is:

1. derive the generic-Q mixed torus response;
2. impose modular covariance before taking Q->1;
3. take the singular limit with the physical normalization fixed;
4. only then freeze a hexagonal/Pell shape target.

This prevents a second round of post-reveal modular-ray fitting.
