# Angular irrep first, radial spectrum second: a corrected research map

Date: 2026-09-14

Status: synthesis after the reverse audit.  This is a research-ordering note, not a theorem or publication plan.

## 1. Why the order matters

The recent programme repeatedly used one observed finite-size power to answer three different questions at once:

```text
which square-lattice angular representation is present?
which continuum scaling dimension controls the radial decay?
which matching/pair-exchange source direction does the correction occupy?
```

These are independent axes.  The current evidence and errata show that mixing them is a major source of over-interpretation.

The corrected order is:

```text
1. angular irrep / geometry,
2. normalized microscopic source direction,
3. radial exponent / logarithmic structure,
4. continuum module / field identification,
5. only then OPE/RG parity claims.
```

## 2. Angular layer: treat `H_(4n)=cos(4n theta)` as the primary basis

For the square model at fixed physical circumference, write

```text
p_ch(ell,theta)
 = P0(ell)
 + P4(ell) H4(theta)
 + P8(ell) H8(theta)
 + P12(ell) H12(theta)
 + ... .
```

This is just a D4/reflection-compatible Fourier decomposition.

The strongest current deterministic result is that the leading orientation-dependent coefficient is overwhelmingly H4-like.  Independent oblique automata now validate the underlying safe-root calculations.

The next task is not to fit another power to raw roots.  It is to isolate `P0,P8,...` after H4 is removed.

### Current gates

- axis/(3,4) equal-ell pair: clean H4 projector at `ell=10`, but with nonzero higher-harmonic leakage;
- #808 N377 pair: exact H4 notch and near H8 notch, much stronger angular-scalar gate;
- N1105 four-angle tomography: attractive algebraically but current automaton hits a prohibitive state-space wall.

Thus #808 is the present high-information experiment; N1105 is an algorithm-development target, not default production.

## 3. Radial layer: each angular coefficient has its own expansion

After angular separation, fit/derive each coefficient independently:

```text
P4(ell) = a4 ell^-4 [1+c4,2 ell^-2+...],
P0(ell)-pc = a0 ell^-q0 [1+...],
P8(ell) = a8 ell^-q8 [1+...],
...
```

The old “4,6,8,... correction ladder” should not be treated as one operator series before this decomposition.

In particular:

- `ell^-6` inside `P4` can be a radial dressing of the leading H4 coefficient;
- an `ell^-7` signal in an H4-null/angular-scalar channel is a different object;
- equality of exponents across angular sectors would be a resonance/mixing clue, not proof of identical mechanism.

## 4. Current `ell^-7` candidate should be typed conservatively

The axis/(3,4) H4-null residual at `ell=5,10` is compatible with a reference-free `ell^-7` extrapolation.  #808 is designed to test this with a geometry that strongly suppresses H8.

A positive result should first be called

```text
post-H4 angular-scalar ell^-7 correction.
```

Only after the angular leakage budget and radial law are stable should it be compared to

```text
x = x_t + 7 = 33/4
```

candidates such as `V_<1,4>` / its Q=1 logarithmic collision block.

The historical #47 warning remains important: an `ell^-7` root correction is not the ordinary next thermal spin-four quasiprimary in the simple tower.

## 5. Source layer is a separate coordinate

A local perturbation should first be represented by its normalized likelihood score and quotiented by constant/thermal nuisance directions.

Use the exact Bernoulli-chaos grading as a microscopic source basis, not as a continuum parity theorem.

The relevant data object is a response matrix

```text
angular channel × residualized source direction × width.
```

A continuum interpretation should explain this matrix, not merely one column/root sequence.

## 6. A three-index correction notation

To avoid future aliasing, label a finite-size correction schematically by

```text
C[a, q, s]
```

where

```text
a : angular harmonic index (0,4,8,...),
q : measured radial root exponent or correction exponent,
s : microscopic residualized source class / intervention label.
```

Only after theory identifies a stable continuum block attach representation/module names.

Example current hypotheses:

```text
C[4,4,uniform-root] : strong deterministic support;
C[0,7,H4-null]      : candidate under #808;
C[8,? ,...]         : nuisance/adversary until isolated;
```

This notation is deliberately phenomenological.  It prevents “spin 4” from being inferred from exponent 4 or “matching odd” from one S/D observable.

## 7. What would count as a real big step now

A high-value result is one of:

1. #808 demonstrates a robust angular-scalar residual whose higher-harmonic leakage is too small to explain it;
2. a residualized degree-3 microscopic source reveals a new response direction not aliased with thermal/pair sources;
3. a transfer/RG construction maps the finite pair-exchange source basis into a continuum tangent block;
4. a direct tagged-operator calculation assigns a distinct topological action to the #800 slow doublet;
5. an algorithm compresses large-memory oblique safe states enough to make 3+ same-circle angles feasible.

A lower-value continuation is another root width or another derived moment that does not change one of these distinctions.

## 8. Claim boundary

Exact/factual inputs:

- D4 angular harmonic algebra;
- independent oblique transfer verification;
- #802 source-normalization/alias erratum;
- N325/N1105 cost wall;
- exact H4-null leakage formulas.

Research organisation:

- angular-first/source-second/radial-third ordering;
- conservative naming of `ell^-7` residual;
- three-index phenomenological correction labels.

The purpose is to make new computations answer one typed question at a time.