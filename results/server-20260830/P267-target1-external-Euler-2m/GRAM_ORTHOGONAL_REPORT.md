# P267 Target 1: source-orthogonal continuation

## Decision

The bulk Euler coupling survives source-plane orthogonalization essentially
unchanged. The frozen same-next-site Gram projection removes only `0.121%` of
the far-D complex-vector magnitude at N325 and `0.086%` at N425. Therefore the
previous D/S transfer-phase lock is not evidence that JD is merely JS in a
different normalization.

This advances the mechanism one step: the lattice rank-birth source has a
component independent of JS in the recorded Gram metric, and `O_far` couples
strongly to it. It still does not identify that source component with the CFT
Q4-epsilon field.

## Frozen rule

The rule was committed before this score as `cd0b3ca`; the covariance-name
implementation correction `c0cca42` changed no statistic:

```text
beta_(N,orientation) = Re <J_D,J_S> / <|J_S|^2>,
J_D_perp = J_D - beta J_S,
Cov(O_far,J_D_perp) = Cov(O_far,J_D) - beta Cov(O_far,J_S).
```

Beta uses only the source Gram at the intrinsic root. It never sees `O_far`,
the N425/N325 transfer or a continuum model. Every delete-one replicate
recomputes the root and both orientation betas after omitting its batch.

## Source coordinates

| size | beta first | beta second | Gram imaginary | normal-equation residual |
|---|---:|---:|---:|---:|
| N325 | `-0.0364052 +/- 0.0002857` | `-0.0366208 +/- 0.0002708` | 0 | at most `4.28e-50` |
| N425 | `-0.0312734 +/- 0.0002714` | `-0.0314657 +/- 0.0002775` | 0 | 0 |

The real-beta rule is not an approximation: the stored imaginary Gram is
exactly zero at scorer precision.

## Orthogonal coupling

| size | complex `P4 Cov(O_far,J_D_perp)` | 2D chi-square | retained far-D magnitude |
|---|---:|---:|---:|
| N325 | `-21.88322(472)+21.91452(500)i` | 641994 | 0.998792 |
| N425 | `-27.16586(613)-18.10241(552)i` | 589232 | 0.999144 |

The transfer is

```text
C_perp(425)/C_perp(325) = 0.206200 + 1.033723 i,
amplitude = 1.054088 +/- 0.001956,
phase = 1.373908 +/- 0.002803 rad.
```

For comparison, unprojected far-D had amplitude transfer `1.053716` and phase
`1.373908`. The orthogonal score is therefore not a small subtraction hiding
inside large correlated errors; it is the dominant production vector.

## Sufficiency boundary

The stored path aggregates are sufficient for beta, the connected coupling,
their full delete-one covariance and the two-size transfer. They do not store
`|J_D|^2`, so the norm `|J_D_perp|^2` and a fractional source-energy removal
cannot be reconstructed. No energy-fraction claim is made. This missing norm
does not enter the coupling score above.

## Mechanism consequence and next experiment

The combined result is now:

1. q-only contact algebra is escaped;
2. a fixed root-R2 Euler neighborhood does not explain the signal;
3. projection onto JS in the recorded source Gram does not explain the signal;
4. nevertheless JD and JS share the same N325-to-N425 complex transport phase.

Thus the surviving ambiguity is no longer bulk versus contact or JD versus a
direct JS projection. It is whether a genuinely independent JD source is
transported by the same geometric C4 phase as JS, or whether it realizes the
thermal Q4 bridge.

Do not extend N along the same one-observer row. The next acquisition should
create a two-observer coupling determinant: retain `O_far` and add one
macroscopically separated, axis+diagonal typed local-H4/arm observer satisfying
the `83e98fc` two-orbit gate. Score both observers against `J_D_perp` and JS in
one common flow. Rank two would separate source directions; rank one would
identify a shared projective transport lane. A charged seam remains a
separate sector test. qJ remains control-only.
