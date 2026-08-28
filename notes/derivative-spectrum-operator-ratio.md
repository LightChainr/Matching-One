# Derivative-spectrum operator-ratio test

## Purpose

Issue #48 found the expected matching-parity pattern retrospectively, but the pure
`P4[S'] ~ N^-5/4` law drifted upward.  Rather than add an unconstrained log or free
exponent, use the same operator competitors already frozen for the central
matching-odd channel in #43.

Let a matching-odd spin-4 irrelevant field have scaling dimension `x`.  At the
finite-size critical center,

```text
P4[D]  ~ N^-beta0,   beta0 = (x-2)/2,
P4[S'] ~ N^-beta1,   beta1 = (x-2-y_t)/2,
```

with percolation `y_t=3/4`, hence

```text
beta1 = beta0 - 3/8.
```

This relation is more restrictive than fitting the two channels independently.

## Frozen operator competitors

The three central-D candidates already declared before N=185/265 target data give:

| candidate | x | beta_D | beta_Sprime |
|---|---:|---:|---:|
| thermal-family level-4 H4 | 21/4 | 13/8 | 5/4 |
| V_<1,3> H4 parity-failure competitor | 14/3 | 4/3 | 23/24 |
| W(2,2) H4 logarithmic/non-singlet leakage competitor | 17/4 | 9/8 | 3/4 |

The old P33 retrospective drift of `P4[S']` is qualitatively compatible with a
slower-decaying competitor, but those data must not be used as confirmatory evidence.

## Parameter-free N=185 -> 265 ratios

Let `q=265/185`.  After the exact spin-4 angular projection, amplitudes cancel:

```text
r_D  = P4[D](265)  / P4[D](185)  = q^-beta_D,
r_Sp = P4[S'](265) / P4[S'](185) = q^-beta_Sprime.
```

Frozen numerical values:

```text
x=21/4: r_D=0.5576728652710071, r_Sp=0.6381272901890036
x=14/3: r_D=0.6193000845309001, r_Sp=0.7086453535146922
x=17/4: r_D=0.6674466940342344, r_Sp=0.7637379846385327
```

These ratios require no source-amplitude refit.

## Operator-identity closure independent of x

If the same field dominates both channels, then

```text
r_Sp / r_D = q^(3/8)
           = 1.1442681362646159.
```

This is independent of the candidate `x`.  It tests the stronger statement

```text
central D and first-thermal-derivative S' are two Taylor coefficients of the
same matching-odd scaling field.
```

Failure of this closure would mean that the S' drift cannot be explained merely by
choosing a different radial dimension for the same operator; it would require
multi-field contamination, logarithmic/Jordan mixing, center-definition effects, or
failure of the proposed matching-parity assignment.

## Scoring protocol

Use fresh N=185/265 threshold-rank statistics from #43/#48.

1. Preserve the original #43 central-D scoring order; do not alter it.
2. Form `P4[D]` and `P4[S']` from the same batch histograms and full covariance.
3. Score the three fixed `(r_D,r_Sp)` pairs above with no target amplitude fit.
4. Independently score the x-free closure `r_Sp/r_D = q^(3/8)` using batchwise/jackknife propagation; do not form a naive ratio of noisy point estimates if the denominator is weak.
5. Only after those scores consider `N^-beta(A+B log N)`, mixed fields, or free exponents.

## Interpretation

- x=21/4 wins both central D and S': strengthens the thermal-family level-4 H4 assignment.
- x=14/3 wins: the V_<1,3> parity-failure route becomes serious and the matching/OPE parity map must be revised.
- x=17/4 wins: investigate special Q=1 logarithmic/non-singlet leakage.
- no single x wins but the x-free closure passes: likely the same field with non-asymptotic radial corrections.
- x-free closure fails: at least two fields or a nontrivial logarithmic/Jordan mixture is required.
