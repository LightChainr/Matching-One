# The noncommon correction as a difference-tangent amplitude hierarchy

Date: 2026-09-14

Status: research reformulation based on the exact square-site Euler source and angular/radial separation.  It weakens unsupported continuum parity language into concrete amplitude-derivative statements that can in principle be checked by transfer/Ward methods.

## 1. Two-colour source coordinates

The exact square-site Euler source has cluster-fugacity component

```text
q_b = log Q_black,
q_w = log Q_white,
```

with the Matching-One source tangent

```text
q_diff=(q_b-q_w)/2,
partial_h q_diff=1,
```

plus fixed local site/edge/plaquette counterterms.

Introduce also

```text
q_common=(q_b+q_w)/2.
```

The physical `h` source is a declared tangent in this enlarged two-colour parameter space, not a postulated scalar parity on continuum fields.

## 2. Expand finite-size blocks by angular irrep and scaling dimension

For a fixed geometry/modulus, write schematically the relevant finite-size free-energy/rank correction as

```text
sum_(x,a) A_a^(x)(q_common,q_diff,...) L^(2-x) H_a(theta),
```

where

```text
H_a(theta)=cos(a theta),
a=0,4,8,...
```

for the square/reflection-even sectors under discussion.

The root response to the Euler source depends on the **difference-tangent amplitude**

```text
D_a^(x)
 := [partial_(q_diff) A_a^(x)]_(0)
    + declared local-counterterm contribution.          (2.1)
```

The statement that a block is “common” is simply

```text
D_a^(x)=0.
```

No OPE-level matching parity is needed to state or test this.

## 3. The large common `x≈4`, H4 correction becomes a precise null condition

Individual/symmetric magnetic-sector data show a substantial ordinary square anisotropy compatible with

```text
x≈4,
a=4.
```

If its difference-tangent coefficient were nonzero, then its root shift would scale with exponent

```text
x-x_t = 4-5/4 = 11/4.
```

No such leading `L^-11/4` Matching-One root correction is observed; the actual leading shift is much smaller, `L^-4`.

Thus the first structural statement to prove is not

```text
“the x=4 field is matching even”
```

but the narrower condition

```text
boxed:
D_4^(4)=0.                                             (3.1)
```

This can arise from exact source symmetry, a Ward identity, cancellation with the local Euler counterterm, or a combination.  The mechanism is an open theorem/interface question.

## 4. The observed leading block says `D_4^(21/4) != 0`

The deterministic oblique-cylinder result identifies a leading angular H4 difference response with root exponent four.

Using

```text
x-root exponent = x-x_t,
```

root exponent four corresponds to

```text
x=21/4.
```

The observable statement is therefore

```text
boxed:
D_4^(21/4) != 0                                      (4.1)
```

for the actual square lattice/source, up to the usual asymptotic interpretation of the finite data.

The competing `eight-arm scalar` and `thermal spin-four descendant` stories are then different explanations of which continuum block produces this nonzero derivative.  The lattice fact is the nonzero H4 amplitude derivative, not either field name.

## 5. Post-H4 scalar becomes the next amplitude derivative

If #808/future leakage-controlled data establish an angular H0 root correction near exponent seven, its linear-scalar interpretation is

```text
D_0^(33/4) != 0,                                      (5.1)
```

because

```text
33/4-5/4=7.
```

But the quadratic H4 mechanism gives a different object: it is second order in the `x=21/4` coupling and produces H0/H8 at root exponent `29/4`, not a new first derivative `D_0^(33/4)`.

This cleanly separates the two theories:

```text
linear scalar:
  first difference tangent of a new x=33/4 H0 block;

quadratic composite:
  second response built from the existing x=21/4 H4 block.
```

## 6. A hierarchy of theorem targets

Instead of assigning continuum parities globally, aim to establish the small list

```text
D_4^(4)      = 0,
D_4^(21/4)   != 0,
D_0^(x<33/4) = 0 or bounded/identified,
D_0^(33/4)   ?
```

with angular and source definitions fixed.

A lower-dimensional scalar/interchiral competitor matters precisely if it gives a nonzero `D_0^(x)` in the actual Euler-source tangent.  This is the correct lattice-side kill test for #61/#768.

## 7. Ward identity target for the first null

The current Ward programme proposes a spin-four finite-size structure involving a quasiprimary `U4+Ubar4` and an `E4(tau)` thermal/coordinate term.

In the present language the useful theorem is to differentiate the **actual rank/Euler projection** along the `q_diff` source and show that the `x=4` H4 combination has zero net coefficient after its contact/seam/local-counterterm pieces are included.

That is a single mixed Ward identity:

```text
partial_(q_diff) [x=4 H4 amplitude of the rank projection] = 0.
```

It is weaker than deriving an OPE automorphism and exactly as strong as needed to explain why the `11/4` root term is absent.

## 8. Empirical response-block version

If the theorem is inaccessible, estimate a small response matrix whose columns are residualized microscopic sources and whose rows are angular/radial blocks.

For each candidate block `(x,a)` fit

```text
R_(x,a),source(L)
```

and test whether the difference-tangent column vanishes or persists.

Degenerate continuum blocks should be treated matrix-valuedly; only the amplitude derivative of the measured block is required.

## 9. Claim boundary

Exact/source-defined:

- Euler difference tangent and local counterterms;
- notion `D_a^(x)` once a finite-size block decomposition is declared.

Empirical/scaling:

- current evidence for `D_4^(21/4) != 0`;
- absence of a visible `11/4` root term as evidence for `D_4^(4)=0`.

Open:

- proof mechanism for the `x=4` null;
- whether the post-H4 scalar is a new linear `D_0^(33/4)` block, a quadratic H4 composite, or a logarithmic/mixed combination.

This reformulation keeps the selection problem attached to the actual source and observable instead of assigning parity to fields before the tangent map is known.